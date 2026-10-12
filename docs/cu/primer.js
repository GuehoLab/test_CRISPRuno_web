// SPDX-License-Identifier: LicenseRef-CRISPRuno — port of CRISPRunoCore.filter_on_primer / trim_primers_single / trim_primers_pair (CRISPRuno 0.1.22, (c) MGH)
'use strict';
/* One pass over the reads: primer/origin trimming with SSW (exactly as CRISPRuno) followed by the cutadapt
 * transposase-adapter step on read 1 (cu/cutadapt.js; `cutadapt -a ADAPTER --minimum-length N [--pair-filter=first]`).
 * Outputs <root>.trimPrimers.[r1.|r2.]trimmed.fq (+ .has_primer.fq if keep_intermediate), .counts.txt, .info and a
 * cutadapt-style report <root>.trimPrimers.cutadapt.transposase_adapter.log.
 */
(function (root) {
const C = root.CuCore || require('./core.js');
const SSW = root.CuSsw || require('./ssw.js');
const CA = root.CuCutadapt || require('./cutadapt.js');
const { strip, read4, str, revcomp } = C;

const fmtN = (n) => n.toLocaleString('en-US');
const pct1 = (a, b) => (b ? (Math.round((1000 * a) / b) / 10).toFixed(1) : '0.0');

async function filterOnPrimer({ ssw, r1, r2, origin_seq, min_primer_aln_score, allow_indels_in_origin_aln, min_primer_length, min_read_length, transposase_adapter_seq,
	out, names, keep_intermediate, progress }) {
	const paired = !!r2;
	const a = ssw.primer(origin_seq), arc = paired ? ssw.primer(revcomp(origin_seq)) : null;
	const adapter = transposase_adapter_seq !== null && transposase_adapter_seq !== undefined ? new CA.BackAdapter(transposase_adapter_seq) : null;
	let tot = 0, postTrim = 0, tooShortOrigin = 0, untrimmed = 0;
	let caWithAdapter = 0, caTooShort = 0, caWritten = 0, bpIn = 0, bpOut = 0;
	try {
		for (;;) {
			const [i1, s1r, p1, q1r] = await read4(r1);
			const id1 = strip(i1), s1 = strip(s1r), q1 = strip(q1r);
			if (!q1) break;
			if (!p1.startsWith('+')) throw new Error(`Fastq ${names.r1} cannot be parsed (${id1}${s1}${p1}${q1}) `);
			tot += 1;
			let id2, s2, p2, q2;
			if (paired) { const r = await read4(r2); id2 = strip(r[0]); s2 = strip(r[1]); p2 = r[2]; q2 = strip(r[3]); }
			const [primerSeq, trimmed, trimmedPos] = SSW.trimLeft(s1, a, min_primer_aln_score, allow_indels_in_origin_aln);
			if (primerSeq.length > min_primer_length) {
				postTrim += 1;
				const nq1 = q1.slice(primerSeq.length);
				const nid = (id1 + '_primer_len=' + (trimmedPos + 1)).split(' ').join('_');
				let t2 = null, nq2 = null;
				if (paired) { [t2] = SSW.trimRight(s2, arc, min_primer_aln_score, allow_indels_in_origin_aln); nq2 = q2.slice(0, t2.length); }
				if (keep_intermediate && out.hp1) { out.hp1.write(nid + '\n' + trimmed + '\n' + p1 + nq1 + '\n'); if (paired) out.hp2.write(nid + '\n' + t2 + '\n' + p2 + nq2 + '\n'); }
				// cutadapt on read 1
				bpIn += trimmed.length;
				let keep = trimmed.length;
				if (adapter) { const m = adapter.matchTo(trimmed); if (m) { caWithAdapter += 1; keep = m.rstart; } }
				if (adapter && keep < min_read_length) { caTooShort += 1; }
				else {
					caWritten += 1; bpOut += keep;
					out.t1.write(nid + '\n' + trimmed.slice(0, keep) + '\n' + p1 + nq1.slice(0, keep) + '\n');
					if (paired) out.t2.write(nid + '\n' + t2 + '\n' + p2 + nq2 + '\n');
				}
			} else if (primerSeq.length > 0) tooShortOrigin += 1;
			else untrimmed += 1;
			if (progress && tot % 20000 === 0) progress(tot);
		}
	} finally { a.destroy(); if (arc) arc.destroy(); }

	let post = postTrim, adapterTooShort = 'NA', report = null;
	if (postTrim > 0 && adapter) {
		post = caWritten; adapterTooShort = caTooShort;
		const R = paired ? 'Pairs' : 'Reads', w = (s, n, tot) => (s + ':').padEnd(35) + fmtN(n).padStart(12) + (tot === undefined ? '' : ` (${pct1(n, tot)}%)`) + '\n';
		report = `This is cutadapt 5.2 (web port: 3' adapter search re-implemented in JavaScript, same alignment algorithm)\n` +
			`Command line parameters: -a ${transposase_adapter_seq} -o ${names.t1}${paired ? ' -p ' + names.t2 : ''} --minimum-length ${min_read_length}${paired ? ' --pair-filter=first' : ''} ${names.hp1}${paired ? ' ' + names.hp2 : ''}\n\n` +
			'=== Summary ===\n\n' + w(paired ? 'Total read pairs processed' : 'Total reads processed', postTrim) +
			(paired ? w('  Read 1 with adapter', caWithAdapter, postTrim) : w('Reads with adapters', caWithAdapter, postTrim)) +
			'\n== Read fate breakdown ==\n' + w(`${R} that were too short`, caTooShort, postTrim) + w(`${R} written (passing filters)`, caWritten, postTrim) +
			`\nTotal basepairs processed (read 1): ${fmtN(bpIn)} bp\nTotal written (filtered, read 1):   ${fmtN(bpOut)} bp\n`;
	}
	const labels = ['Filtered: no primer', 'Filtered: origin match too short', 'Filtered: trimmed read too short', 'Retained'];
	const values = [untrimmed, tooShortOrigin, adapterTooShort, post];
	const counts = labels.join('\t') + '\n' + values.map(str).join('\t') + '\n';
	return { tot, post, untrimmed, tooShortOrigin, adapterTooShort, report, counts, labels, values, usedAdapter: postTrim > 0 && !!adapter };
}

const api = { filterOnPrimer };
if (typeof module !== 'undefined' && module.exports) module.exports = api;
root.CuPrimer = api;
})(typeof self !== 'undefined' ? self : globalThis);
