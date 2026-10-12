// SPDX-License-Identifier: LicenseRef-CRISPRuno — port of CRISPRunoCore.add_umi_from_umi_file / analyze_UMIs_initial / compute_gini (CRISPRuno 0.1.22, (c) MGH)
'use strict';
(function (root) {
const C = root.CuCore || require('./core.js');
const { strip, read4, str, F, DefaultInt } = C;

/* add_umi_from_umi_file: id line + ':' + UMI sequence. rd: {r1, r2|null, umi} readers; w1/w2 sinks */
async function addUmi({ r1, r2, umi, w1, w2, names, progress }) {
	let n = 0;
	for (;;) {
		const [id1, s1, p1, q1] = await read4(r1);
		if (!q1) break;
		if (!p1.startsWith('+')) throw new Error(`Fastq ${names.r1} cannot be parsed (${strip(id1)}${s1}${p1}${q1}) `);
		n += 1;
		const [, us, , uq] = await read4(umi);
		const u = strip(us);
		if (!uq) throw new Error(`Insufficient number of UMI entries in the UMI file ${names.umi}`);
		w1.write(strip(id1) + ':' + u + '\n' + s1 + p1 + q1);
		if (r2) { const [id2, s2, p2, q2] = await read4(r2); w2.write(strip(id2) + ':' + u + '\n' + s2 + p2 + q2); }
		if (progress && n % 50000 === 0) progress(n);
	}
	if (n === 0) throw new Error(`UMI command failed. Got no reads from ${names.r1} and ${names.r2}`);
	return n;
}

/* compute_gini(np.array(counts)): sum_{i<j} |xi-xj| / (n^2 * mean) */
function gini(counts) {
	const n = counts.length;
	if (n === 0) return null;
	const a = Float64Array.from(counts).sort();
	let total = 0n, prefix = 0n, S = 0;
	for (let j = 0; j < n; ++j) { const x = BigInt(a[j]); total += x * BigInt(j) - prefix; prefix += x; S += a[j]; }
	const mean = S / n;
	return Number(total) / ((n * n) * mean);
}

function umiRegex(umiRegexStr) {
	let s = umiRegexStr;
	for (const [k, v] of [['I', '([ATCG])'], ['N', '([ATCG])'], ['R', '([AG])'], ['Y', '([CT])'], ['S', '([GC])'], ['W', '([AT])'], ['K', '([GT])'], ['M', '([AC])'], ['B', '([CGT])'], ['D', '([AGT])'], ['H', '([ACT])'], ['V', '([ACG])']]) s = s.split(k).join(v);
	return s;
}

/* analyze_UMIs_initial. open() -> fresh readers {r1, r2|null} (the R1 file is read twice). */
async function analyzeUmis({ open, umi_regex, dedup_input_on_UMI, min_umi_seen_to_keep_read, write_UMI_counts, w1, w2, wCounts, names, outNames, progress, log }) {
	const head = '# experiment_had_UMIs (bool): whether the experiment had UMIs or not\n';
	const cols = ['experiment_had_UMIs', 'fastq_analyze_UMI_r1', 'fastq_analyze_UMI_r2', 'tot_read_count', 'count_with_regex', 'post_analyze_UMI_count', 'post_analyze_UMI_read_count', 'UMI_input_GINI', 'UMI_count'];
	// check the first 10 reads for UMIs
	let rd = await open();
	let hasUmi = true;
	for (let k = 1; ; ++k) {
		const [id, s, p, q] = await read4(rd.r1);
		if (k > 10) break;
		const qq = strip(q);
		if (!qq) break;
		if (!p.startsWith('+')) throw new Error(`Fastq ${names.r1} cannot be parsed (${strip(id)}${strip(s)}${p}${qq}) `);
		const u = strip(id).split(':').pop().toUpperCase();
		const bad = [...new Set(u)].filter((c) => !'ATCGN'.includes(c));
		if (u.length < 3 || bad.length) hasUmi = false;
	}
	if (!hasUmi) {
		log('UMIs were not detected in this sample');
		return { had: false, info: head + cols.join('\t') + '\n' + ['False', ...new Array(8).fill('None')].join('\t') + '\n' };
	}
	let re = null, reStr = '.*';
	if (umi_regex !== null) { reStr = umiRegex(umi_regex); re = new RegExp('^(?:' + reStr + ')$'); }

	let tot = 0, withRegex = 0, collision = 0, collisionReads = 0, tooFew = 0, tooFewReads = 0, postCount = 0, postReadCount = 0;
	const keyMost = new Map(), keyCounts = new DefaultInt(), keyBest = new Map(), keyBestQ = new Map(), seqCounts = new DefaultInt();
	rd = await open();
	for (;;) {
		const [a, b, c, d] = await read4(rd.r1);
		const id1 = strip(a), s1 = strip(b), p1 = c, q1 = strip(d);
		if (!q1) break;
		if (!p1.startsWith('+')) throw new Error(`Fastq ${names.r1} cannot be parsed (${id1}${s1}${p1}${q1}) `);
		tot += 1;
		let id2 = '', s2 = '', p2 = '', q2 = '';
		if (rd.r2) { const r = await read4(rd.r2); id2 = strip(r[0]); s2 = strip(r[1]); p2 = r[2]; q2 = strip(r[3]); }
		const u = id1.split(':').pop().toUpperCase();
		if (re) { if (re.test(u)) withRegex += 1; else continue; } // Python re.match('^...$'): '$' also matches before a final '\n' (UMI has none)
		const key = u + ' # ' + s1 + s2;
		if (dedup_input_on_UMI) {
			let qs = 0; const qq = q1 + q2; for (let i = 0; i < qq.length; ++i) qs += qq.charCodeAt(i) & 0xff;
			const rec = [id1 + '\n' + s1 + '\n' + p1 + q1, id2 + '\n' + s2 + '\n' + p2 + q2];
			if (!keyCounts.has(key)) { keyCounts.inc(key); keyBestQ.set(key, qs); keyBest.set(key, rec); }
			else { keyCounts.inc(key); if (keyBestQ.get(key) < qs) { keyBestQ.set(key, qs); keyBest.set(key, rec); } }
			if (!seqCounts.has(u)) { seqCounts.set(u, 1); keyMost.set(u, key); }
			else { seqCounts.inc(u); if (keyCounts.get(key) > keyCounts.get(keyMost.get(u))) keyMost.set(u, key); }
		} else {
			postCount += 1; postReadCount += 1; seqCounts.inc(u);
			w1.write(id1 + '\n' + s1 + '\n' + p1 + q1 + '\n');
			if (rd.r2) w2.write(id2 + '\n' + s2 + '\n' + p2 + q2 + '\n');
		}
		if (progress && tot % 50000 === 0) progress(tot);
	}
	if (tot === 0) throw new Error(`UMI dedup failed. Got no reads from ${names.r1} and ${names.r2}`);
	const umiList = [...seqCounts.keys()].sort((x, y) => seqCounts.get(x) - seqCounts.get(y)); // stable, like sorted(key=count)
	const umiCount = umiList.length;
	log(`Read ${tot} reads`); log(`Processed ${umiCount} UMIs`);
	if (umiCount === 1 && umiList[0] === '') throw new Error("Error: only the empty barcode '' was found.");
	const g = umiCount > 0 ? gini(umiList.map((k) => seqCounts.get(k))) : null;
	if (dedup_input_on_UMI) {
		for (const u of umiList) {
			const mk = keyMost.get(u), mc = keyCounts.get(mk);
			if (mc < min_umi_seen_to_keep_read) { tooFew += 1; tooFewReads += seqCounts.get(u); }
			else {
				if (mc !== seqCounts.get(u)) { collision += 1; collisionReads += seqCounts.get(u) - mc; }
				const [x1, x2] = keyBest.get(mk);
				w1.write(x1 + '\n'); if (rd.r2) w2.write(x2 + '\n');
				postCount += 1; postReadCount += keyCounts.get(mk);
			}
		}
	}
	if (write_UMI_counts && dedup_input_on_UMI && wCounts) {
		// as in CRISPRuno: iterates the UMI/read-sequence keys and looks them up as UMIs (defaultdicts -> 0)
		const info = umi_regex !== null ? `Only valid UMIs matching the regex ${umi_regex} (${reStr}) are reported here` : '';
		wCounts.write('#UMI: the UMI sequence. ' + info + '\n#seen_count: the number of times the UMI was seen in the input (the same UMI may be paired with reads with different sequences)\n#printed_count: the number of times the UMI/read pair printed was seen in the input\nUMI\tseen_count\tprinted_count\n');
		const keys = [...keyCounts.keys()].sort(pyStrCmp);
		for (const k of keys) {
			const mk = keyMost.has(k) ? keyMost.get(k) : (keyMost.set(k, 0), 0);
			const seen = seqCounts.get(k), printed = keyCounts.get(mk);
			wCounts.write(`${k}\t${seen}\t${printed}\n`);
		}
		log('Wrote UMI counts to ' + outNames.counts);
	}
	if (umi_regex === null) withRegex = tot;
	log(`Wrote ${postCount} reads with UMIs`);
	const info = head +
		'# fastq_analyze_UMI_r1 (str): file containing analyzed reads - filtered for correct UMIs if umi_regex is provided, and deduplicated if dedup_input_on_UMI is provided\n' +
		'# fastq_analyze_UMI_r2 (str): file containing analyzed r2 reads or None if input is single-ended\n' +
		'# tot_read_count (int): number of reads in input\n' +
		`# count_with_regex (int): number of reads matching the umi_regex '${str(umi_regex)}'\n` +
		'# post_analyze_UMI_count (int): number of reads post UMI analysis (after asserting UMI regex and/or initial deduplication on UMI)\n' +
		'# post_analyze_UMI_read_count (int): number of reads that contributed to the deduplicated reads (post_analyze_UMI_read_count reads were seen that were deduplicated to post_analyze_UMI_count reads. If post_analyze_UMI_count == post_analyze_UMI_read_count then every UMI-read pair was seen once.)\n' +
		'# UMI_input_GINI (float): GINI index for distribution of reads among UMIs (passing the umi_regex filter if umi_regex is set). (0 = perfect equality and 1 = perfect inequality in reads/UMI)\n' +
		'# UMI_count (int): Number of unique UMI sequences observed\n' +
		cols.join('\t') + '\n' + [true, outNames.r1, rd.r2 ? outNames.r2 : null, tot, withRegex, postCount, postReadCount, g === null ? null : F(g), umiCount].map(str).join('\t') + '\n';
	return { had: true, tot, withRegex, postCount, postReadCount, gini: g, umiCount, collision, collisionReads, tooFew, tooFewReads, info };
}
function pyStrCmp(a, b) { // Python str ordering (code points)
	const n = Math.min(a.length, b.length);
	for (let i = 0; i < n; ++i) { const x = a.codePointAt(i), y = b.codePointAt(i); if (x !== y) return x - y; }
	return a.length - b.length;
}

const api = { addUmi, analyzeUmis, gini, umiRegex, pyStrCmp };
if (typeof module !== 'undefined' && module.exports) module.exports = api;
root.CuUmi = api;
})(typeof self !== 'undefined' ? self : globalThis);
