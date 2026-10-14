// SPDX-License-Identifier: LicenseRef-CRISPRuno — CRISPRunoCore.processCRISPRuno split into steps (CRISPRuno 0.1.22, (c) MGH; academic non-commercial license)
'use strict';
/* Steps (each reads what the previous ones wrote in the run folder, like CRISPRuno's own cache files):
 *   prep       prep_input: primer/guide positions (minibwa alignment of the sequences, see locateFastq), origin sequence,
 *              read length/count -> <root>.primerInfo.info
 *   umi        add_umi_from_umi_file (if fastq_umi) + analyze_UMIs_initial -> <root>.addUMI.*, <root>.analyzeInputUMI.*
 *   primer     filter_on_primer (SSW + cutadapt) -> <root>.trimPrimers.*
 *   (align)    minibwa on the trimmed reads (align-worker.js), then sortAlignment -> <root>.genomeAlignment.sam
 *   assign     make_final_read_assignments -> <root>.final.*
 *   crispresso prep_crispresso2 + CRISPResso2 quantification + run_and_aggregate_crispresso -> <root>.CRISPResso.*
 *   summary    make_final_summary -> <root>.summary.txt (the HTML report is made by cu/report.js)
 * io: see cu-worker.js (OPFS) and test/node-io.js (files).
 * Cas-OFFinder is not part of this port.
 */
(function (root) {
const C = root.CuCore || require('./core.js');
const UMI = root.CuUmi || require('./umi.js');
const PR = root.CuPrimer || require('./primer.js');
const SAM = root.CuSam || require('./sam.js');
const AS = root.CuAssign || require('./assign.js');
const CR = root.CuCrispresso || require('./crispresso.js');
const AL = root.CuAlign || require('./cralign.js');
const { str, F, strip } = C;

const stateName = (S) => S.root + '.web_state.json';
async function loadState(io, S) { const t = await io.readText(stateName(S)); return t ? JSON.parse(t) : {}; }
async function saveState(io, S, st) { const w = await io.writer(stateName(S)); w.write(JSON.stringify(st, null, 1)); w.close(); }
async function writeText(io, name, text) { const w = await io.writer(name); w.write(text); w.close(); }
const plotJson = (p) => (p ? JSON.stringify({ name: p.name, title: p.title, label: p.label, datas: p.datas, order: p.order, plot: p.plot }) : 'None');

/* ---------- prep_input ---------- */
function locateFastq(S) { // the sequences CRISPRuno aligns with bowtie2 --end-to-end -c, as reads for minibwa
	const recs = [['primer', S.primer_seq], ...S.guide_sequences.map((g, i) => ['guide' + i, g])];
	return recs.map(([n, s]) => `@${n}\n${s}\n+\n${'I'.repeat(s.length)}\n`).join('');
}
function parseLocateSam(text) {
	const hits = new Map();
	for (const l of text.split('\n')) {
		if (!l || l[0] === '@') continue;
		const f = l.split('\t'), fl = +f[1];
		if (fl & 0x900) continue; // primary records only
		if (fl & 4) { hits.set(f[0], null); continue; }
		const lead = /^(\d+)S/.exec(f[5]);
		hits.set(f[0], { chr: f[2], pos: parseInt(f[3], 10) - (lead ? parseInt(lead[1], 10) : 0), rc: (fl & 16) !== 0, mapq: parseInt(f[4], 10), cigar: f[5], line: l });
	}
	return hits;
}
async function stepPrep(io, S, { locateSam }) {
	const log = io.log, hits = parseLocateSam(locateSam);
	const cutSites = [], ann = new Map();
	let genomic = false, primerChr = 'exogenous', primerLoc = 'None', primerRc = false;
	log('Finding genomic coordinates of primer ' + S.primer_seq);
	const ph = hits.get('primer');
	if (ph) {
		primerChr = ph.chr; primerLoc = ph.pos; primerRc = ph.rc;
		if (primerRc) primerLoc += S.primer_seq.length;
		if (ph.mapq === 0) log('Primer alignment resulted in an ambiguous alignment (MAPQ 0)', 'warn');
		log(`Setting genomic coordinates for primer aligned to ${primerChr}:${primerLoc} (minibwa ${ph.cigar}, MAPQ ${ph.mapq})`);
		genomic = true;
	} else log('Primer is not aligned to genome. Assuming exogenous primer sequence');
	S.guide_sequences.forEach((g, i) => {
		log('Finding genomic coordinates of guide ' + g);
		const h = hits.get('guide' + i);
		if (!h) throw new Error('Could not find unique genomic coordinates for guide ' + g);
		let cut, dir;
		if (!h.rc) { cut = h.pos + g.length + S.cleavage_offset; dir = 'FW'; } else { cut = h.pos - S.cleavage_offset; dir = 'RC'; }
		if (h.mapq === 0) log('Guide alignment resulted in an ambiguous alignment (MAPQ 0)', 'warn');
		log(`Setting cut site for guide ${g} to ${h.chr}:${cut} (minibwa ${h.cigar}, MAPQ ${h.mapq})`);
		const key = h.chr + ':' + cut;
		cutSites.push(key); ann.set(key, ['Programmed', 'Not-origin', dir, g]);
	});
	let closest = -1, closestDist = Infinity;
	if (genomic) for (const site of cutSites) {
		const [chr, p] = site.split(':'); const pos = parseInt(p, 10);
		if (chr === primerChr) {
			if (closest === -1) { closestDist = Math.abs(primerLoc - pos); closest = pos; }
			else if (Math.abs(primerLoc - pos) < closestDist) { closestDist = Math.abs(primerLoc - pos); closest = pos; }
		}
	}
	let origin = S.primer_seq;
	if (genomic && closest === -1) { genomic = false; log(`Primer was found to be genomic (${primerChr}:${primerLoc}) but no cut sites were found on chromosome ${primerChr}. Assuming primer is the origin sequence`, 'warn'); }
	else if (genomic) {
		const genome = await io.getGenome();
		let os = primerLoc, oe = closest;
		const k = primerChr + ':' + closest;
		ann.get(k)[1] = 'Origin:left';
		if (primerRc) { os = closest; oe = primerLoc; ann.get(k)[1] = 'Origin:right'; }
		if (os > oe) throw new Error('Primer does not point toward the closest cut site. Please check the orientation of the primer sequence parameter');
		origin = C.faidx(genome, primerChr, os, oe - 1);
		if (origin === '') throw new Error('Cannot retrieve sequence information from fasta file');
		if (primerRc) origin = C.revcomp(origin);
		log('Got origin sequence: ' + origin);
	}
	if (origin.length < S.min_primer_length) throw new Error(`The min_primer_length parameter must be less than the length of the primer/origin sequence (Origin: ${origin} (length ${origin.length}), min_primer_length: ${S.min_primer_length})`);
	// get_av_read_len / get_num_reads_fastq (awk over all reads; wc -l)
	const rd = await io.userReader(S.fastq_r1);
	let lines = 0, idx = 0, sum = 0, n = 0; // awk 'NR%4 == 2 {s+=length($0);n++}' and wc -l
	for (;;) {
		if (!rd.has(1)) await rd.ensure(1);
		const l = rd.readline(); if (!l) break;
		const nl = l.endsWith('\n');
		if (idx % 4 === 1) { sum += nl ? l.length - 1 : l.length; n += 1; }
		if (nl) lines += 1;
		idx += 1;
		if (idx % 400000 === 0) io.progress({ reads: idx / 4, label: 'reads' });
	}
	if (n === 0) throw new Error('Cannot compute the average read length of ' + S.fastq_r1 + ' (no reads)');
	const avLen = Math.trunc(sum / n), nReads = Math.trunc((lines + 1) / 4);
	const info = ['origin_seq', 'primer_chr', 'primer_loc', 'primer_seq', 'primer_is_genomic', 'cut_site_count', 'av_read_length', 'num_reads_input'].join('\t') + '\n' +
		[origin, primerChr, primerLoc, S.primer_seq, genomic, cutSites.length, avLen, nReads].map(str).join('\t') + '\n' + cutSites.map((c) => c + '\t' + ann.get(c).join(',') + '\n').join('');
	await writeText(io, S.root + '.primerInfo.info', info);
	log(`${nReads} reads in input; average read length ${avLen}`);
	return { origin, cutSites: cutSites.length, nReads };
}
async function readPrimerInfo(io, S) {
	const t = await io.readText(S.root + '.primerInfo.info');
	if (!t) throw new Error(S.root + '.primerInfo.info not found — run the prep step first');
	const L = t.split('\n'), v = L[1].split('\t');
	const ann = new Map(), sites = [];
	for (const l of L.slice(2)) { if (!l) continue; const e = l.split('\t'); sites.push(e[0]); ann.set(e[0], e[1].split(',')); }
	return { origin_seq: v[0], primer_chr: v[1] === 'None' ? null : v[1], primer_loc: v[2], primer_is_genomic: v[4] === 'True', av_read_length: parseFloat(v[6]), num_reads_input: parseInt(v[7], 10), cut_sites: sites, cut_annotations: ann };
}

/* ---------- UMIs ---------- */
async function stepUmi(io, S) {
	const log = io.log, st = await loadState(io, S);
	let cur1 = { name: S.primer_in_r2 ? S.fastq_r2 : S.fastq_r1, user: true }, cur2 = (S.primer_in_r2 ? S.fastq_r1 : S.fastq_r2) ? { name: S.primer_in_r2 ? S.fastq_r1 : S.fastq_r2, user: true } : null;
	const open = async (x) => (x ? (x.user ? io.userReader(x.name) : io.reader(x.name)) : null);
	if (S.fastq_umi) {
		log('Adding UMIs to input fastqs');
		const n1 = S.root + '.addUMI.r1.fq', n2 = cur2 ? S.root + '.addUMI.r2.fq' : null;
		const w1 = await io.writer(n1), w2 = n2 ? await io.writer(n2) : null;
		const n = await UMI.addUmi({ r1: await open(cur1), r2: await open(cur2), umi: await io.userReader(S.fastq_umi), w1, w2, names: { r1: cur1.name, r2: cur2 ? cur2.name : 'None', umi: S.fastq_umi }, progress: (k) => io.progress({ reads: k, label: 'add UMI' }) });
		w1.close(); if (w2) w2.close();
		await writeText(io, S.root + '.addUMI.log', 'fastq_r1_with_UMI\tfastq_r2_with_UMI\ttot_read_count\n' + [n1, n2, n].map(str).join('\t') + '\n');
		log(`Added UMIs for ${n} reads`);
		cur1 = { name: n1, user: false }; cur2 = n2 ? { name: n2, user: false } : null;
	}
	const a1 = S.root + '.analyzeInputUMI.r1.fq', a2 = cur2 ? S.root + '.analyzeInputUMI.r2.fq' : null;
	const w1 = await io.writer(a1), w2 = a2 ? await io.writer(a2) : null;
	const wc = S.write_UMI_counts && S.dedup_input_on_UMI ? await io.writer(S.root + '.analyzeInputUMI.umiCounts.txt') : null;
	const r = await UMI.analyzeUmis({ open: async () => ({ r1: await open(cur1), r2: await open(cur2) }), umi_regex: S.umi_regex, dedup_input_on_UMI: S.dedup_input_on_UMI,
		min_umi_seen_to_keep_read: S.min_umi_seen_to_keep_read, write_UMI_counts: S.write_UMI_counts, w1, w2, wCounts: wc, names: { r1: cur1.name, r2: cur2 ? cur2.name : 'None' },
		outNames: { r1: a1, r2: a2, counts: S.root + '.analyzeInputUMI.umiCounts.txt' }, progress: (k) => io.progress({ reads: k, label: 'UMI' }), log });
	w1.close(); if (w2) w2.close(); if (wc) wc.close();
	await writeText(io, S.root + '.analyzeInputUMI.info', r.info);
	if (!r.had) { await io.remove(a1); if (a2) await io.remove(a2); }
	const nIn = (await readPrimerInfo(io, S)).num_reads_input;
	st.experiment_had_UMIs = r.had;
	st.post_UMI_regex_count = r.had ? r.withRegex : nIn;
	st.post_initial_dedup_count = r.had ? r.postCount : nIn;
	if (r.had) { cur1 = { name: a1, user: false }; cur2 = a2 ? { name: a2, user: false } : null; }
	st.curr_r1 = cur1; st.curr_r2 = cur2;
	st.intermediate = [S.fastq_umi ? S.root + '.addUMI.r1.fq' : null, S.fastq_umi && cur2 ? S.root + '.addUMI.r2.fq' : null, r.had ? a1 : null, r.had ? a2 : null].filter(Boolean);
	await saveState(io, S, st);
	return { had: r.had, reads: r.had ? r.postCount : nIn };
}

/* ---------- primer filtering ---------- */
async function stepPrimer(io, S, { ssw }) {
	const log = io.log, st = await loadState(io, S);
	if (!st.curr_r1) throw new Error('run the UMI step first');
	const P = await readPrimerInfo(io, S);
	const open = async (x) => (x ? (x.user ? io.userReader(x.name) : io.reader(x.name)) : null);
	const paired = !!st.curr_r2, R = S.root + '.trimPrimers';
	const names = paired ? { t1: R + '.r1.trimmed.fq', t2: R + '.r2.trimmed.fq', hp1: R + '.r1.has_primer.fq', hp2: R + '.r2.has_primer.fq' } : { t1: R + '.trimmed.fq', hp1: R + '.has_primer.fq' };
	const out = { t1: await io.writer(names.t1), t2: paired ? await io.writer(names.t2) : null };
	if (S.keep_intermediate) { out.hp1 = await io.writer(names.hp1); out.hp2 = paired ? await io.writer(names.hp2) : null; }
	log('Filtering and trimming primers from reads');
	const r = await PR.filterOnPrimer({ ssw, r1: await open(st.curr_r1), r2: await open(st.curr_r2), origin_seq: P.origin_seq, min_primer_aln_score: S.min_primer_aln_score,
		allow_indels_in_origin_aln: S.allow_indels_in_origin_aln, min_primer_length: S.min_primer_length, min_read_length: S.min_read_length, transposase_adapter_seq: S.transposase_adapter_seq,
		out, names: { ...names, r1: st.curr_r1.name }, keep_intermediate: S.keep_intermediate, progress: (k) => io.progress({ reads: k, label: 'primer' }) });
	for (const w of Object.values(out)) if (w) w.close();
	if (r.report) { log('Trimming transposase sequence from reads'); await writeText(io, R + '.cutadapt.transposase_adapter.log', r.report); }
	await writeText(io, R + '.counts.txt', r.counts);
	let plot = null;
	if (!S.suppress_plots) {
		const label = "Reads were filtered to contain primer origin sequence <p class='text-break text-monospace'>" + P.origin_seq + '</p>' + 'Reads must contain at least ' + S.min_primer_length +
			' bases that match the primer/origin sequence<br>Reads containing the primer/origin sequence but shorter than ' + S.min_read_length + 'bp were filtered as too short.';
		plot = new C.PlotObject(R + '.counts', 'Read filtering based on primer/origin presence', label + '<br>' + r.labels.map((l, i) => `${l} N=${str(r.values[i])}`).join('<br>'),
			[['Primer/origin information', R + '.info'], ['Primer filtering statistics', R + '.counts.txt']], 3, { type: 'pie', items: r.labels.map((l, i) => [l, r.values[i]]), title: 'Read filtering based on primer/origin presence' });
	}
	const f1 = r.usedAdapter ? names.t1 : names.hp1, f2 = paired ? (r.usedAdapter ? names.t2 : names.hp2) : null;
	await writeText(io, R + '.info', 'filtered_on_primer_fastq_r1\tfiltered_on_primer_fastq_r2\tpost_trim_read_count\n' + [f1, f2, r.post].map(str).join('\t') + '\n' + plotJson(plot) + '\n');
	log(`Keeping ${r.post} reads with primer sequence`);
	st.post_filter_on_primer_read_count = r.post;
	st.aligned_inputs = paired ? [names.t1, names.t2] : [names.t1];
	await saveState(io, S, st);
	return { post: r.post, tot: r.tot, inputs: st.aligned_inputs };
}

/* ---------- alignment post-processing ---------- */
async function stepSortAlignment(io, S, { headerText, sources, alignerLog }) {
	const log = io.log, st = await loadState(io, S);
	const name = S.root + '.genomeAlignment.sam';
	const out = await io.writer(name);
	out.write(headerText.endsWith('\n') || !headerText ? headerText : headerText + '\n');
	const r = await SAM.sortSam({ headerText, sources, keepSecondary: S.analyze_multimap_assignments, out, tmp: io.tmp, progress: (f) => io.progress({ frac: f, label: 'sort' }) });
	out.close();
	await writeText(io, S.root + '.genomeAlignment.alignerLog', alignerLog || '');
	await writeText(io, S.root + '.genomeAlignment.info', 'mapped_bam_file\n' + name + '\n');
	log(`[sort] ${r.kept.toLocaleString()} records (${r.supplementary.toLocaleString()} supplementary and ${r.secondary.toLocaleString()} secondary records dropped)`);
	if (!S.keep_intermediate) {
		for (const f of [...(st.intermediate || []), ...(st.aligned_inputs || [])]) {
			if (st.curr_r1 && f === st.curr_r1.name) continue; // the CRISPResso input
			await io.remove(f);
		}
	}
	return r;
}

/* ---------- final read assignments ---------- */
async function stepAssign(io, S) {
	const st = await loadState(io, S), P = await readPrimerInfo(io, S);
	const src = await io.source(S.root + '.genomeAlignment.sam');
	if (!src) throw new Error(S.root + '.genomeAlignment.sam not found — run the align step first');
	const IO = root.GsIO || require('../lib/textio.js');
	const sam = function* () { for (const l of IO.lines(src)) if (l.length && l.charCodeAt(0) !== 64) yield l; };
	let regionLines = null;
	if (S.cut_region_annotation_file) {
		const rd = await io.userReader(S.cut_region_annotation_file); regionLines = [];
		for (;;) { if (!rd.has(1)) await rd.ensure(1); const l = rd.readline(); if (!l) break; regionLines.push(l); }
	}
	const genome = S.suppress_homology_detection ? null : await io.getGenome();
	const r = await AS.makeFinalReadAssignments({ io: { ...io, sam }, S, genome, origin_seq: P.origin_seq, cut_sites: P.cut_sites, cut_annotations: P.cut_annotations,
		experiment_had_UMIs: !!st.experiment_had_UMIs, log: io.log, progress: (k) => io.progress({ reads: k, label: 'assign' }), root: S.root + '.final', regionLines, keepIntermediate: S.keep_intermediate });
	// (CRISPRuno deletes the alignment here unless keep_intermediate; the web version keeps it in the run folder so the
	// assign step can be re-run with other filter settings without aligning again. 'Delete run' frees the space.)
	return { final: r.finalTotal, classes: r.clsCounts };
}

/* ---------- CRISPResso ---------- */
async function stepCrispresso(io, S, { ednafull, fixSpaceIds }) {
	const log = io.log, st = await loadState(io, S), P = await readPrimerInfo(io, S);
	const finalName = S.root + '.final.final_assignments.txt';
	const fsrc = await io.source(finalName);
	if (!fsrc) throw new Error(finalName + ' not found — run the assign step first');
	const IO = root.GsIO || require('../lib/textio.js');
	const lines = function* () { for (const l of IO.lines(fsrc)) if (l.length) yield l; };
	// as CRISPRuno rebuilds them from a finished final_assignments file
	const ids = new C.StrMap(1), cuts = [], cutIdx = new Map(), cutCounts = new C.DefaultInt(); // read name -> index of its final cut key
	let first = true;
	for (const l of lines()) {
		if (first) { first = false; continue; }
		const e = l.split('\t');
		cutCounts.inc(e[16]);
		if (Math.abs(parseInt(e[17], 10)) <= S.crispresso_max_indel_size) { let c = cutIdx.get(e[16]); if (c === undefined) { c = cuts.length; cuts.push(C.own(e[16])); cutIdx.set(cuts[c], c); } let i = ids.find(e[0]); if (i < 0) i = ids.insert(e[0]); ids.setV(i, 0, c); }
	}
	const readIds = { get: (k) => { const i = ids.find(k); return i < 0 ? undefined : cuts[ids.v(i)]; } };
	const info = (await io.readText(S.root + '.final.info')).split('\n');
	const cls = new Map(Object.entries(JSON.parse(info[6])));
	log('Preparing reads for analysis by CRISPResso2');
	const genome = await io.getGenome();
	const rd = st.curr_r1.user ? await io.userReader(st.curr_r1.name) : await io.reader(st.curr_r1.name);
	const prep = await CR.prepCrispresso({ reader: rd, readIds, cutCounts, cls, cutAnnotations: P.cut_annotations, avReadLength: P.av_read_length, originSeq: P.origin_seq, genome, S,
		root: S.root + '.CRISPResso_r1', fixSpaceIds, log });
	log(`Running and analyzing ${prep.infos.length} alignments using CRISPResso2 (web port of the quantification)`);
	const R = S.root + '.CRISPResso';
	const writer = async (name, text) => writeText(io, name, text);
	const appenders = new Map();
	writer.append = async (name, text, close) => {
		let w = appenders.get(name); if (!w) { w = await io.writer(name); appenders.set(name, w); }
		w.write(text); if (close) { w.close(); appenders.delete(name); }
	};
	const res = await CR.runAndAggregate({ infos: prep.infos, matrix: AL.readMatrix(ednafull), S, finalLines: lines, writer, root: R, log, progress: (f) => io.progress({ frac: f, label: 'CRISPResso' }) });
	if (S.keep_intermediate) for (const x of prep.infos) await writeText(io, x.readsFile.split('/').join('.'), x.reads.map((q) => `${q[2]}\n${q[1]}\n${q[3]}${q[4]}\n`).join(''));
	if (!S.keep_intermediate && st.curr_r1 && !st.curr_r1.user) await io.remove(st.curr_r1.name);
	return { sites: prep.infos.length, analyzed: res.runs.length };
}

/* ---------- final summary ---------- */
async function stepSummary(io, S) {
	const st = await loadState(io, S), P = await readPrimerInfo(io, S);
	const info = (await io.readText(S.root + '.final.info') || '').split('\n');
	if (info.length < 7) throw new Error('run the assign step first');
	const v = info[1].split('\t').map((x) => parseInt(x, 10));
	const discarded = [['Unaligned', v[4]], ['Duplicate read', v[5]], ['Not supported by R2', v[6]], ['Bad alignment', v[7]]];
	const zip = (a, b) => a.split('\t').map((x, i) => [x, parseInt(b.split('\t')[i], 10)]);
	const clsCounts = info[2] ? zip(info[2], info[3]) : [], clsIndel = info[4] ? zip(info[4], info[5]) : [];
	const nIn = P.num_reads_input, umiR = st.post_UMI_regex_count, dd = st.post_initial_dedup_count, prim = st.post_filter_on_primer_read_count, fin = v[8];
	const pct = (a, b) => (b > 0 ? C.pctStr(100 * a / b) : 'None');
	let s = 'Total input reads: ' + nIn + '\n';
	const head = ['total_input_reads'], vals = [nIn];
	s += `\tSurvived assertion of UMI regex: ${umiR}/${nIn} (${pct(umiR, nIn)}%)\n`; head.push('post_UMI_regex_reads'); vals.push(umiR);
	s += `\tSurvived initial UMI deduplication: ${dd}/${umiR} (${pct(dd, umiR)}%)\n`; head.push('post_initial_dedup_reads'); vals.push(dd);
	s += `\tSurvived filtering for primer/origin presence: ${prim}/${dd} (${pct(prim, dd)}%)\n`; head.push('post_primer_filter_reads'); vals.push(prim);
	const dTot = discarded.reduce((a, x) => a + x[1], 0), surv = prim - dTot;
	s += `\tSurvived quality filtering: ${surv}/${prim} (${pct(surv, prim)}%)\n`; head.push('post_quality_filter_reads'); vals.push(surv);
	s += `\tFailed quality filtering: ${dTot}/${prim} (${pct(dTot, prim)}%)\n`; head.push('discarded_quality_filter_reads'); vals.push(dTot);
	for (const [c, n] of discarded) { s += `\t\t${c}: ${n}/${dTot} (${pct(n, dTot)}%)\n`; head.push(c); vals.push(n); }
	s += `\tReads used for final analysis: ${fin}/${prim} (${pct(fin, prim)}%)\n`; head.push('analyzed_read_count'); vals.push(fin);
	for (const [c, n] of clsCounts) { s += `\t\t${c}: ${n}/${fin} (${pct(n, fin)}%)\n`; head.push(c); vals.push(n); }
	s += `\tReads used for final analysis (including breakdown of indels): ${fin}/${prim} (${pct(fin, prim)}%)\n`;
	for (const [c, n] of clsIndel) { s += `\t\t${c}: ${n}/${fin} (${pct(n, fin)}%)\n`; head.push(c); vals.push(n); }
	await writeText(io, S.root + '.summary.txt', head.join('\t') + '\n' + vals.join('\t') + '\n' + s + '\n');
	const labels = ['Did not pass UMI regex filter', 'Deduplicated by initial UMI deduplication', 'Did not contain origin/primer seqeunce', 'Failed quality filtering', 'Used for final analyses'];
	const values = [nIn - umiR, umiR - dd, dd - prim, prim - fin, fin];
	const plot = S.suppress_plots ? null : new C.PlotObject(S.root + '.summary', 'Read summary', 'Assignment summary for N=' + nIn + ' input reads:<br>' + labels.map((l, i) => `${l} N=${values[i]}`).join('<br>'),
		[['Read assignment summary', S.root + '.summary.txt']], 1, { type: 'pie', items: labels.map((l, i) => [l, values[i]]), title: 'Read Assignment Summary' });
	await writeText(io, S.root + '.summary.plot.json', plotJson(plot) + '\n');
	return { final: fin };
}

const api = { locateFastq, parseLocateSam, stepPrep, stepUmi, stepPrimer, stepSortAlignment, stepAssign, stepCrispresso, stepSummary, readPrimerInfo, loadState };
if (typeof module !== 'undefined' && module.exports) module.exports = api;
root.CuPipeline = api;
})(typeof self !== 'undefined' ? self : globalThis);
