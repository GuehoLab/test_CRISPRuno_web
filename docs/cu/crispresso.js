// SPDX-License-Identifier: LicenseRef-CRISPRuno — port of CRISPRunoCore.prep_crispresso2 / run_and_aggregate_crispresso (CRISPRuno 0.1.22)
// and of the quantification core of CRISPResso2 2.3.4 for the call CRISPRuno makes (both (c) MGH; academic non-commercial license)
'use strict';
/* CRISPRuno runs, per fragment junction with enough reads:
 *   CRISPResso -o DATA -n NAME --default_min_aln_score MIN -a ORIGIN+GENOMIC -g ORIGIN[-15:]+GENOMIC[:15] -wc -15 -w W -r1 READS.fq
 *              -gn 'Fragment junction' --fastq_output --quantification_window_coordinates A-B
 * and keeps from it only: the alignment counts (N_TOT_READS, counts_total/modified/unmodified/discarded and the
 * insertion/deletion/substitution combinations) and, per read, CLASS=Reference_MODIFIED or not (fastq_output).
 * crispressoQuant() reproduces exactly that part of CRISPResso2 (seeds, global_align with the cut-site gap incentive,
 * find_indels_substitutions in the quantification window, read caching, the reverse-complement merge, and the
 * conditions under which a CRISPResso run fails: then CRISPRuno writes NA). No CRISPResso HTML reports/plots.
 */
(function (root) {
const C = root.CuCore || require('./core.js');
const A = root.CuAlign || require('./cralign.js');
const { str, faidx, revcomp } = C;

class CrispressoError extends Error {}
const rcShared = (s) => { // CRISPRessoShared.reverse_complement: upper() first, KeyError for non ACGTN_- characters
	const T = { A: 'T', C: 'G', G: 'C', T: 'A', N: 'N', _: '_', '-': '-' }; const u = s.toUpperCase(); let o = '';
	for (let i = u.length - 1; i >= 0; --i) { const c = T[u[i]]; if (c === undefined) throw new CrispressoError(`Unexpected error: KeyError '${u[i]}'`); o += c; }
	return o;
};
const reEsc = (s) => s.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');

/* CRISPResso2 for one amplicon/guide; reads: array of [id, seq]; returns stats or throws CrispressoError */
function crispressoQuant({ amplicon, guide, wc, w, coords, minAlnScore, reads, matrix }) {
	const SEED_COUNT = 5, SEED_LEN = 10, SEED_MIN = 2, GAP_OPEN = -20, GAP_EXT = -2, GAP_INC = 1, EXCL = 15, PLOT_WIN = 20;
	// input checks
	if (!reads.length) throw new CrispressoError('The input contains 0 reads.');
	for (let i = 0; i < Math.min(25, reads.length); ++i) if (!reads[i][1].trim().length) throw new CrispressoError('File is not in fastq format! (empty line)');
	if (C.findWrongNt(guide).length) throw new CrispressoError('The sgRNA sequence contains bad characters');
	const ampUp = amplicon.toUpperCase();
	if (!ampUp.length) throw new CrispressoError('Amplicon 1 has length 0.');
	if (C.findWrongNt(ampUp).length) throw new CrispressoError('Reference amplicon sequence 0 (Reference) contains invalid characters');
	const ref = ampUp.trim().toUpperCase(), L = ref.length;
	// get_amplicon_info_for_guides
	const exclude = new Set();
	for (let i = 0; i < Math.min(EXCL, L); ++i) exclude.add(i);
	for (let i = Math.max(0, L - EXCL); i < L; ++i) exclude.add(i);
	const cutPoints = [];
	const offFw = wc + guide.length - 1, offRc = -wc - 1;
	const findAll = (pat) => { const out = []; const re = new RegExp('(?=(' + reEsc(pat) + '))', 'gi'); let m; while ((m = re.exec(ref))) { out.push(m.index); re.lastIndex = m.index + 1; } return out; };
	let guideInclude = [];
	if (guide !== '') {
		for (const s of findAll(guide)) { const c = s + offFw; cutPoints.push(c); if (w > 0) for (let x = Math.max(0, c - w + 1); x < Math.min(L - 1, c + w + 1); ++x) guideInclude.push(x); }
		for (const s of findAll(rcShared(guide))) { const c = s + offRc; cutPoints.push(c); if (w > 0) for (let x = Math.max(0, c - w + 1); x < Math.min(L - 1, c + w + 1); ++x) guideInclude.push(x); }
		if (!cutPoints.length) throw new CrispressoError(`The guide sequence 0 (${guide}) provided is not present in the amplicon sequences!`);
	}
	let include, given;
	if (coords !== null && coords !== '0') {
		given = [];
		for (const c of String(coords).split('_')) {
			const m = /^(\d+)-(\d+)$/.exec(c);
			if (!m) throw new CrispressoError(`Cannot parse analysis window coordinate '${c}'`);
			const a = parseInt(m[1], 10), b = parseInt(m[2], 10) + 1;
			if (b > L) throw new CrispressoError(`End coordinate ${b} for '${c}' is longer than the sequence length (${L})`);
			for (let x = a; x < b; ++x) given.push(x);
		}
		include = given;
	} else if (cutPoints.length && guideInclude.length) { include = guideInclude; given = guideInclude; }
	else { include = Array.from({ length: L }, (_, i) => i); given = []; }
	const incSet = new Set(include.filter((x) => !exclude.has(x)));
	for (const x of given) if (!incSet.has(x)) throw new CrispressoError('The quantification window has been partially exluded by the --exclude_bp_from_left or --exclude_bp_from_right parameters.');
	if (incSet.size === 0) throw new CrispressoError(include.length ? 'The quantification window around the sgRNA is excluded.' : 'The entire sequence has been excluded.');
	const gapInc = new Int32Array(L + 1);
	for (const c of cutPoints) { // numpy item assignment: negative index wraps, out of range raises
		let i = c + 1;
		if (i >= L + 1 || i < -(L + 1)) throw new CrispressoError(`Unexpected error: index ${i} is out of bounds for axis 0 with size ${L + 1}`);
		if (i < 0) i += L + 1;
		gapInc[i] = GAP_INC;
	}
	// seeds
	const refRc = rcShared(ref), seeds = [], rcSeeds = [];
	for (let st = EXCL; st < L - EXCL - SEED_LEN; st += SEED_COUNT) {
		let tries = 0, s0 = st, pot = ref.slice(s0, s0 + SEED_LEN);
		while (refRc.includes(pot) || seeds.includes(pot)) {
			tries += 1; if (tries > 100) break;
			if (st > L - SEED_LEN) s0 = 0;
			s0 += 1; pot = ref.slice(s0, s0 + SEED_LEN);
		}
		const prc = rcShared(pot);
		if (ref.includes(prc)) continue;
		if (!refRc.includes(pot)) { seeds.push(pot); rcSeeds.push(prc); }
	}
	// reads: cache by sequence
	const cache = new Map();
	for (const [, s] of reads) cache.set(s, (cache.get(s) || 0) + 1);
	const variants = new Map(), notAln = new Set();
	for (const [seq, count] of cache) {
		let fw = 0, rv = 0;
		for (let i = 0; i < SEED_COUNT && i < seeds.length; ++i) { if (seq.includes(seeds[i])) fw += 1; if (seq.includes(rcSeeds[i])) rv += 1; }
		let s1, s2, score;
		if (fw > SEED_MIN && rv === 0) [s1, s2, score] = A.globalAlign(seq, ref, matrix, gapInc, GAP_OPEN, GAP_EXT);
		else if (fw === 0 && rv > SEED_MIN) [s1, s2, score] = A.globalAlign(rcShared(seq), ref, matrix, gapInc, GAP_OPEN, GAP_EXT);
		else {
			const f = A.globalAlign(seq, ref, matrix, gapInc, GAP_OPEN, GAP_EXT), r = A.globalAlign(rcShared(seq), ref, matrix, gapInc, GAP_OPEN, GAP_EXT);
			[s1, s2, score] = r[2] > f[2] ? r : f;
		}
		if (score > minAlnScore && score > 0) {
			const p = A.findIndelsSubstitutions(s1, s2, incSet);
			variants.set(seq, { count, p, modified: p.deletion_n > 0 || p.insertion_n > 0 || p.substitution_n > 0 });
		} else notAln.add(seq);
	}
	// statistics (with CRISPResso's merge of a sequence and its reverse complement)
	const st = { N_TOT_READS: reads.length, total: 0, modified: 0, unmodified: 0, discarded: 0, insertion: 0, deletion: 0, substitution: 0,
		only_insertion: 0, only_deletion: 0, only_substitution: 0, insertion_and_deletion: 0, insertion_and_substitution: 0, deletion_and_substitution: 0, insertion_and_deletion_and_substitution: 0 };
	let nTotal = 0;
	for (const [seq, v] of variants) {
		let n = v.count;
		if (n === 0) continue;
		const rc = rcShared(seq), o = variants.get(rc);
		if (o && o.count > 0) { n += o.count; o.count = 0; v.count = n; }
		nTotal += n;
		st.total += n;
		if (v.modified) st.modified += n; else st.unmodified += n;
		const ins = v.p.insertion_n > 0, del = v.p.deletion_n > 0, sub = v.p.substitution_n > 0;
		if (ins) st.insertion += n; if (del) st.deletion += n; if (sub) st.substitution += n;
		if (del) { if (ins) { if (sub) st.insertion_and_deletion_and_substitution += n; else st.insertion_and_deletion += n; } else if (sub) st.deletion_and_substitution += n; else st.only_deletion += n; }
		else if (ins) { if (sub) st.insertion_and_substitution += n; else st.only_insertion += n; }
		else if (sub) st.only_substitution += n;
	}
	if (nTotal === 0) throw new CrispressoError('No alignments were found');
	// fastq_output: per read class
	const perRead = [];
	let modReads = 0;
	for (const [id, seq] of reads) {
		const v = notAln.has(seq) ? null : variants.get(seq);
		const m = v ? v.modified : false;
		if (m) modReads += 1;
		perRead.push([id, m ? 'M' : 'U']);
	}
	return { st, perRead, modReads };
}

/* prep_crispresso2: reads the R1 FASTQ (CRISPRuno's curr_r1_file) and groups the reads of each junction */
async function prepCrispresso({ reader, readIds, cutCounts, cls, cutAnnotations, avReadLength, originSeq, genome, S, root, fixSpaceIds, log }) {
	const dataDir = root + '_data';
	const printed = new Map(), names = new Map();
	let total = 0, printedN = 0, spaceIds = 0;
	const count = (k) => cutCounts.peek(k);
	for (;;) {
		const [i1, s1, p1, q1] = await C.read4(reader);
		const id = C.strip(i1), seq = C.strip(s1), qual = C.strip(q1);
		if (!qual) break;
		total += 1;
		let rid = id.split('\t')[0].slice(1);
		if (rid.includes(' ')) { spaceIds += 1; if (fixSpaceIds) rid = rid.split(' ').join('_'); }
		const cut = readIds.get(rid);
		if (cut === undefined) continue;
		let n = count(cut);
		const e = cut.split(':'), pos = e[0] + ':' + e[1];
		n = count(pos + ':left') + count(pos + ':right');
		if (n < S.crispresso_min_count || count(cut) === 0) continue;
		if (!names.has(cut)) names.set(cut, (cls.get(cut) + ':' + cut).split(' ').join('_').split(':').join('_'));
		printedN += 1;
		if (!printed.has(cut)) printed.set(cut, []);
		printed.get(cut).push([C.own(fixSpaceIds ? rid : id.slice(1)), C.own(seq.toUpperCase()), C.own(id), C.own(p1), C.own(qual)]);
	}
	if (spaceIds && !fixSpaceIds) log(`[crispresso] ${spaceIds.toLocaleString()} read names contain a space: CRISPRuno 0.1.22 replaces spaces with '_' in the aligned read names but not when it looks the reads up for CRISPResso, so these reads are not analyzed (same here; the web option "Fix read IDs with spaces" changes this)`, 'warn');
	const half = avReadLength / 1.5;
	const infos = [];
	for (const [cut, rs] of printed) {
		const [chr, ps, dir] = cut.split(':'); const pos = parseInt(ps, 10);
		let amp;
		if (dir === 'right') amp = faidx(genome, chr, pos, Math.trunc(pos + half));
		else if (dir === 'left') amp = revcomp(faidx(genome, chr, Math.trunc(pos - half), pos - 1));
		else throw new Error('Got unexpected cut direction: ' + dir + ' from cut ' + cut);
		const ampSeq = originSeq + amp, guide = originSeq.slice(-15) + amp.slice(0, 15);
		const qs = originSeq.length - S.crispresso_quant_window_size;
		const coords = `${qs}-${qs + 2 * S.crispresso_quant_window_size - 1}`;
		const name = names.get(cut), readsFile = dataDir + '/' + name + '.fq';
		const cmd = `${S.crispresso_command} -o ${dataDir} -n ${name} --default_min_aln_score ${S.crispresso_min_aln_score} -a ${ampSeq} -g ${guide} -wc -15 -w ${S.crispresso_quant_window_size} -r1 ${readsFile} -gn 'Fragment junction' --fastq_output --quantification_window_coordinates ${coords} ${S.n_processes > 4 ? '--n_processes 4' : ''} > ${readsFile}.log 2>&1`;
		let run = 'True';
		if (!S.run_crispresso_on_novel_sites && cls.get(cut).includes('Novel')) run = 'False: novel location';
		infos.push({ cut, name, type: cls.get(cut), cut_loc: chr + ':' + pos, amp_seq: ampSeq, guide, coords, reads: rs, printed_read_count: rs.length, command: cmd, run, readsFile });
	}
	return { infos, total, printed: printedN };
}

/* run_and_aggregate_crispresso */
async function runAndAggregate({ infos, matrix, S, finalLines, writer, root, log, progress }) {
	const modData = new Map();
	const head = 'name\tcut_annotation\tcut_type\tcut_location\treference\treads_printed\tn_total\treads_aligned\treads_unmod\treads_mod\treads_discarded\t' +
		'reads_insertion\treads_deletion\treads_substitution\treads_only_insertion\treads_only_deletion\treads_only_substitution\treads_insertion_and_deletion\t' +
		'reads_insertion_and_substitution\treads_deletion_and_substitution\treads_insertion_and_deletion_and_substitution\tamplicon_sequence\n';
	const rows = [], cmds = [], runs = [];
	let k = 0;
	for (const x of infos) {
		cmds.push(x.command);
		let v = new Array(15).fill('NA'), refName = 'NA';
		if (x.run === 'True') {
			try {
				const r = crispressoQuant({ amplicon: x.amp_seq, guide: x.guide, wc: -15, w: S.crispresso_quant_window_size, coords: x.coords, minAlnScore: S.crispresso_min_aln_score, reads: x.reads.map((q) => [q[0], q[1]]), matrix });
				refName = 'Reference';
				const s = r.st;
				v = [s.N_TOT_READS, s.total, s.unmodified, s.modified, s.discarded, s.insertion, s.deletion, s.substitution, s.only_insertion, s.only_deletion, s.only_substitution,
					s.insertion_and_deletion, s.insertion_and_substitution, s.deletion_and_substitution, s.insertion_and_deletion_and_substitution];
				for (const [id, m] of r.perRead) modData.set(id, m);
				if (s.modified !== r.modReads) throw new Error("Couldn't get correct number of modified reads from the fastq_output file of " + x.name);
				runs.push({ name: x.name, cut: x.cut, stats: s });
				log(`[crispresso] ${x.name}: ${s.N_TOT_READS} reads, ${s.total} aligned, ${s.modified} modified`);
			} catch (e) {
				if (!(e instanceof CrispressoError)) throw e;
				log(`[crispresso] ${x.name}: CRISPResso would fail (${e.message}) — NA values, as in CRISPRuno`, 'warn');
				v = new Array(15).fill('NA');
			}
		} else v = new Array(15).fill('NA');
		rows.push([x.name, x.cut, x.type, x.cut_loc, refName, x.run, x.printed_read_count, ...v, x.amp_seq].map(str).join('\t') + '\n');
		if (progress) progress(++k / infos.length);
	}
	await writer(root + '.summary.txt', head + rows.join(''));
	await writer(root + '.commands.txt', cmds.map((c) => c + '\n').join(''));
	// annotated final assignments
	const counts = new C.DefaultInt();
	let out = [], first = true;
	for (const line of finalLines()) {
		if (first) {
			const h = line.split('\t');
			if (h[19] !== 'classification') throw new Error("Couldn't parse final assignment file");
			out.push(line + '\tcrispresso_status\n'); first = false; continue;
		}
		const e = line.split('\t');
		let status = '(not analyzed)';
		if (modData.has(e[0])) status = modData.get(e[0]) === 'M' ? 'short indels' : 'no indels';
		if (e[19] !== 'Unidentified') counts.inc(e[19] + ' ' + status);
		out.push(line + '\t' + status + '\n');
		if (out.length >= 50000) { await writer.append(root + '.annotated_final_assignments', out.join('')); out = []; }
	}
	await writer.append(root + '.annotated_final_assignments', out.join(''), true);
	const labels = [...counts.keys()].sort((a, b) => (a < b ? -1 : a > b ? 1 : 0)), vals = labels.map((l) => counts.get(l));
	await writer(root + '.classificationSummary.txt', labels.join('\t') + '\n' + vals.join('\t') + '\n');
	const sum = vals.reduce((a, b) => a + b, 0);
	let plot = null;
	if (sum > 0 && !S.suppress_plots) plot = new C.PlotObject(root + '.classificationSummary', 'Read classification including CRISPResso analysis',
		'CRISPResso classification<br>' + labels.map((l, i) => `${l} N=${vals[i]}`).join('<br>') + '<br>Reads labeled "(not analyzed)" had too few reads for CRISPResso analysis or indels (based on genome alignment) that were too long. ' +
		'You can adjust the minimum number of reads to run CRISPResso with the parameter "crispresso_min_count". You can adjust the maximum indel length for analysis by CRISPResso with the parameter "crispresso_max_indel_size".',
	[['CRISPResso classification summary', root + '.classificationSummary.txt']], 40, { type: 'pie', items: labels.map((l, i) => [l, vals[i]]) });
	await writer(root + '.info', (plot ? JSON.stringify({ name: plot.name, title: plot.title, label: plot.label, datas: plot.datas, order: plot.order, plot: plot.plot }) : 'None') + '\n\nRun Name\tRun sub-html\n');
	return { runs, plot };
}

const api = { crispressoQuant, prepCrispresso, runAndAggregate, CrispressoError };
if (typeof module !== 'undefined' && module.exports) module.exports = api;
root.CuCrispresso = api;
})(typeof self !== 'undefined' ? self : globalThis);
