// SPDX-License-Identifier: LicenseRef-CRISPRuno — port of CRISPRunoCore.make_final_read_assignments, collapse_cut_points,
// getLeftRightMismatchesMDZ, get_guide_match_from_aln, compute_gini (CRISPRuno 0.1.22, (c) MGH; academic non-commercial license)
'use strict';
/* Input: the coordinate-sorted SAM (as `samtools view` prints it from the sorted BAM).
 * Output: the same text files as CRISPRuno (<root>.final.*); plots are described by PlotObjects and drawn as SVG by cu/plots.js.
 * Differences on purpose: SAM lines are split without their trailing newline (in CRISPRuno an MD:Z tag at the very end of a
 * line would keep the '\n' and count 0 right-side matches; bowtie2 never ends a line with MD:Z, minibwa can).
 */
(function (root) {
const C = root.CuCore || require('./core.js');
const A = root.CuAlign || require('./cralign.js');
const U = root.CuUmi || require('./umi.js');
const { str, F, DefaultInt, revcomp, faidx, PlotObject, pyJson } = C;
const pyCmp = U.pyStrCmp;

function leftRightMatchesMDZ(md) {
	let i = 0; while (i < md.length && md.charCodeAt(i) >= 48 && md.charCodeAt(i) <= 57) i++;
	const left = i ? parseInt(md.slice(0, i), 10) : 0;
	let j = md.length; while (j > 0 && md.charCodeAt(j - 1) >= 48 && md.charCodeAt(j - 1) <= 57) j--;
	const right = j < md.length ? parseInt(md.slice(j), 10) : 0;
	return [left, right];
}
function guideMatchFromAln(g, r, cutInd) { // -> [match, mismatch, gaps, ref_ind_cut_site, potential_guide_in_ref]
	if (g.length !== r.length) throw new Error('Sequence lengths must be the same');
	let first = -1, last = -1, gi = 0, ri = 0, cutSite = -1;
	for (let i = 0; i < g.length; ++i) {
		if (g[i] !== '-') { gi += 1; last = i; if (first === -1) first = i; }
		if (r[i] !== '-') ri += 1;
		if (gi === cutInd) cutSite = ri;
	}
	let m = 0, mm = 0, gaps = 0;
	for (let i = first; i <= last; ++i) {
		if (i < 0) continue;
		if (g[i] === '-' || r[i] === '-') gaps += 1; else if (g[i] === r[i]) m += 1; else mm += 1;
	}
	const pot = r.slice(first < 0 ? Math.max(0, r.length + first) : first, last + 1).split('-').join('');
	return [m, mm, gaps, cutSite, pot];
}

/* collapse_cut_points (including its quirks: the list is modified while iterating, so some points
 * that were assigned to a known site are still clustered as novel points) */
function collapseCutPoints(novelDist, knownDist, originDist, cutPointsByChr, cutAnnotations, homology) {
	let originChr = null, originCutPos = null, originDirection = null;
	const known = new Map(), finalByChr = new Map(), lookup = new Map();
	for (const [site, anno] of cutAnnotations) {
		const [chr, p] = site.split(':'); const pos = parseInt(p, 10);
		if (anno[1] === 'Origin:right') { originChr = chr; originCutPos = pos; originDirection = 'right'; }
		else if (anno[1] === 'Origin:left') { originChr = chr; originCutPos = pos; originDirection = 'left'; }
		if (!known.has(chr)) known.set(chr, []);
		known.get(chr).push(pos);
		if (!finalByChr.has(chr)) finalByChr.set(chr, new DefaultInt());
		finalByChr.get(chr).set(pos, 0);
	}
	const closestOKH = new Map();
	let oMin = -1, oMax = -1;
	if (originCutPos !== null) { oMin = originCutPos - originDist; oMax = originCutPos + originDist; }
	if (!cutPointsByChr.has(originChr)) cutPointsByChr.set(originChr, new DefaultInt()); // defaultdict access
	for (const cp of cutPointsByChr.get(originChr).keys()) if (cp >= oMin && cp <= oMax) closestOKH.set(originChr + ':' + cp, originChr + ':' + originCutPos);
	const kh = new Map();
	for (const [chr, arr] of known) kh.set(chr, arr.slice());
	for (const key of homology.keys()) {
		const [chr, p] = key.split(':'); const pos = parseInt(p, 10);
		if (!kh.has(chr)) kh.set(chr, []);
		if (!kh.get(chr).includes(pos)) kh.get(chr).push(pos);
	}
	const num = (a, b) => a - b;
	for (const [chr, pts] of cutPointsByChr) {
		const these = [...pts.keys()].sort(num);
		const khs = (kh.get(chr) || (kh.set(chr, []), kh.get(chr))).slice().sort(num);
		if (khs.length === 0) { for (const cp of these) closestOKH.set(chr + ':' + cp, null); continue; }
		let ci = 0;
		for (const cp of these) {
			const dCur = Math.abs(cp - khs[ci]);
			if (ci + 1 >= khs.length) closestOKH.set(chr + ':' + cp, dCur <= knownDist ? chr + ':' + khs[ci] : null);
			else {
				const dNext = Math.abs(cp - khs[ci + 1]);
				if (dCur <= dNext) closestOKH.set(chr + ':' + cp, dCur <= knownDist ? chr + ':' + khs[ci] : null);
				else { ci += 1; closestOKH.set(chr + ':' + cp, dNext <= knownDist ? chr + ':' + khs[ci] : null); }
			}
		}
	}
	let found = 0;
	const closest = new Map();
	for (const [chr, pts] of cutPointsByChr) {
		found += pts.size;
		const these = [...pts.keys()].sort(num);
		for (let it = 0; it < these.length;) { // Python: for x in list: ... list.remove(x)  (skips the element after a removed one)
			const cp = these[it]; it += 1;
			const k = chr + ':' + cp;
			if (closestOKH.get(k) !== null) { these.splice(these.indexOf(cp), 1); closest.set(k, closestOKH.get(k)); } // 'it' is not moved back: the next element is skipped, as in Python
		}
		if (these.length > 0) {
			let cur = [these[0]];
			for (let i = 1; i <= these.length; ++i) {
				let sum = 0, tot = 0;
				for (const p of cur) { sum += p * pts.get(p); tot += pts.get(p); }
				const mean = tot > 0 ? sum / tot : 0;
				if (i === these.length || Math.abs(these[i] - mean) > novelDist) {
					const pos = Math.trunc(mean);
					for (const p of cur) closest.set(chr + ':' + p, chr + ':' + pos);
					cur = [];
				}
				if (i < these.length) cur.push(these[i]);
			}
		}
	}
	let nAll = 0, nHom = 0, nOrig = 0, nColl = 0;
	for (const [chr, pts] of cutPointsByChr) {
		for (const cp of pts.keys()) {
			nAll += 1;
			const key = chr + ':' + cp;
			let asg;
			if (closestOKH.get(key) !== null && closestOKH.get(key) !== undefined) { asg = closestOKH.get(key); nHom += 1; }
			else if (chr === originChr && cp >= oMin && cp <= oMax) { asg = originChr + ':' + originCutPos; nOrig += 1; }
			else { asg = closest.get(key); nColl += 1; }
			const apos = parseInt(asg.split(':')[1], 10);
			lookup.set(key, asg);
			if (!finalByChr.has(chr)) finalByChr.set(chr, new DefaultInt());
			finalByChr.get(chr).inc(apos, pts.get(cp));
		}
	}
	return { originChr, originCutPos, originDirection, lookup, finalByChr, found, log: [`Processed ${nAll} observed cut points`, `${nHom} cut points were collapsed to known/given/homologous cut sites`, `${nOrig} cut points were collapsed to the origin site`, `${nColl} cut points were collapsed to nearby cut sites`] };
}

/* main. io: { sam(): iterator of SAM record lines (sorted, no header), writer(name) -> {write, close}, tmpWriter(name) -> {write, close() -> lines iterator factory},
 *             remove(name) }   S: settings; genome: {has, length, slice}; regionLines: lines of cut_region_annotation_file or null */
async function makeFinalReadAssignments(ctx) {
	const { io, S, genome, origin_seq, cut_sites, cut_annotations, experiment_had_UMIs, log, progress } = ctx;
	const R = ctx.root; // '<root>.final'
	const P = {
		r1_r2_support_max_distance: S.r1_r2_support_max_distance, novel: S.novel_cut_merge_distance, known: S.known_cut_merge_distance, origin: S.origin_cut_merge_distance,
		short_indel_length_cutoff: 50, gaps: S.guide_homology_max_gaps, mismatches: S.guide_homology_max_mismatches, armMin: S.arm_min_matched_start_bases, armClip: S.arm_max_clipped_bases,
		crispresso_max_indel_size: S.crispresso_max_indel_size, suppressDedup: S.suppress_dedup_on_aln_pos_and_UMI_filter, dedupFinal: S.dedup_by_final_cut_assignment_and_UMI,
		suppressR2: S.suppress_r2_support_filter, minQ: S.min_alignment_quality_score, minQ2: S.min_r2_alignment_quality_score, suppressPoor: S.suppress_poor_alignment_filter,
		multimap: S.analyze_multimap_assignments, writeDiscarded: S.write_discarded_read_info, suppressHom: S.suppress_homology_detection, suppressPlots: S.suppress_plots,
	};
	const sup = 'Supported by R2', notOri = 'Not supported by R2 - different R1/R2 orientations', notDist = 'Not supported by R2 - alignment separation more than maximum',
		notChr = 'Not supported by R2 - different chrs', notAln = 'Not supported by R2 - unaligned', notLow = 'Not supported by R2 - R2 alignment score too low';
	const supLabels = [sup, notOri, notDist, notChr, notAln, notLow];

	const seen = new Map(); let seenCounts = new Map();
	let totalReads = 0, totalR1 = 0, dUnaligned = 0, dDup = 0, dPoor = 0, dR2 = 0;
	const cutPointsByChr = new Map();
	const cp = (chr) => { if (!cutPointsByChr.has(chr)) cutPointsByChr.set(chr, new DefaultInt()); return cutPointsByChr.get(chr); };
	const tlens = new DefaultInt(), chrCounts = new DefaultInt(), supDist = new DefaultInt(), notSupDist = new DefaultInt();
	const discardedName = R + '.discarded_read_info.txt';
	const fd = P.writeDiscarded ? await io.writer(discardedName) : null;

	// multimapping: choose one alignment per read (closest to an Origin/Programmed cut on equal AS)
	const mmAssign = new Map(), mmAS1 = new Map(), mmAS2 = new Map(), mmR1Count = new DefaultInt();
	if (P.multimap) {
		log('Analyzing multimapping assignments');
		const prio = [];
		for (const site of cut_sites) { const a = cut_annotations.get(site); if (a[0].includes('Programmed') || a[1].includes('Origin')) prio.push(site.split(':')); }
		for (const line of io.sam()) {
			const f = line.split('\t'), fl = +f[1], chr = f[2], pos = parseInt(f[3], 10), sec = fl & 0x100, isR2 = fl & 0x80;
			if (!isR2) mmR1Count.inc(f[0]);
			let as = 0; for (const x of f) if (x.startsWith('AS:i:')) { as = parseInt(x.slice(5), 10); break; }
			let best = null, bestD = Infinity;
			for (const p of prio) if (p[0] === chr) { const d = Math.abs(parseInt(p[1], 10) - pos); if (d < bestD) { bestD = d; best = p; } }
			const M = isR2 ? mmAS2 : mmAS1, cur = M.get(f[0]);
			if (!cur) M.set(f[0], [as, bestD, line]);
			else if (as > cur[0] || (as === cur[0] && bestD < cur[1])) M.set(f[0], [as, bestD, line]);
			if (!mmAssign.has(f[0])) mmAssign.set(f[0], [chr, pos, bestD, sec]);
			else {
				const o = mmAssign.get(f[0]);
				if (best !== null) { if (o[2] === -1) mmAssign.set(f[0], [chr, pos, bestD, sec]); else if (bestD < o[2]) mmAssign.set(f[0], [chr, pos, bestD, sec]); }
				else if (o[2] === -1 && !sec) mmAssign.set(f[0], [chr, pos, -1, sec]);
			}
		}
	}
	let mmSecUsed = 0, mmSkipped = 0;
	// R2 info
	const r2Info = new Map();
	for (const line of io.sam()) {
		const f = line.split('\t'), fl = +f[1];
		if (!(fl & 0x80)) continue;
		if (P.multimap && mmAssign.has(f[0])) {
			if (mmAS2.get(f[0])[2] !== line) { mmSkipped += 1; if (!(fl & 0x100)) mmSecUsed += 1; continue; }
		}
		const q = parseInt(f[4], 10);
		let nm = null; for (const x of f) if (x.startsWith('NM:i')) { nm = parseInt(x.slice(5), 10); break; }
		r2Info.set(f[0].split(' ')[0], [fl & 0x4, q, nm]);
		totalReads += 1;
	}
	const supCounts = new DefaultInt();
	const tmpName = R + '.final_assignments.tmp';
	const tmp = await io.writer(tmpName);
	tmp.write(['read_id', 'aln_pos', 'cut_pos', 'cut_direction', 'del_primer', 'ins_target', 'insert_size', 'r1_r2_support_status', 'r1_r2_support_dist', 'r1_r2_orientation_str',
		'r1_alignment_quality', 'r1_alignment_mismatches', 'r2_alignment_quality', 'r2_alignment_mismatches', 'umi_pos_dedup_key'].join('\t') + '\n');
	const CIG = /(\d+)([MIDNSHPX=])/g;
	let nLines = 0;
	for (const line of io.sam()) {
		const f = line.split('\t'), fl = +f[1];
		if (fl & 0x80) continue;
		if (progress && ++nLines % 50000 === 0) progress(nLines);
		if (P.multimap && mmAssign.has(f[0])) {
			if (mmAS1.get(f[0])[2] !== line) { mmSkipped += 1; if (!(fl & 0x100)) mmSecUsed += 1; continue; }
		}
		totalReads += 1;
		const multi = fl & 0x1, isR1 = fl & 0x40;
		if (multi && !isR1) continue;
		totalR1 += 1;
		const mmSel = P.multimap && mmR1Count.peek(f[0]) > 1;
		if ((fl & 0x100) && !mmSel) { dUnaligned += 1; if (fd) fd.write(f[0] + '\tUnaligned\t' + f[1] + '\n'); continue; }
		const q = parseInt(f[4], 10);
		if (q < P.minQ && !mmSel) { dPoor += 1; if (fd) fd.write(f[0] + '\tLow quality alignment\t' + f[4] + '\n'); continue; }
		if (fl & 0x4) { dUnaligned += 1; if (fd) fd.write(f[0] + '\tUnmapped\t' + f[1] + '\n'); continue; }
		const parts = f[0].split('_primer_len=');
		if (parts.length !== 2) throw new Error('Error parsing primer length from read name: ' + f[0]);
		const lineInfo = parts[0], primerLen = parseInt(parts[1], 10);
		if (!/^\s*[+-]?\d+\s*$/.test(parts[1])) throw new Error(`ValueError: invalid literal for int() with base 10: '${parts[1]}'`);
		let left = 0, right = 0;
		for (const x of f) if (x.startsWith('MD:Z')) { [left, right] = leftRightMatchesMDZ(x.slice(5)); break; }
		let nm = null; for (const x of f) if (x.startsWith('NM:i')) { nm = parseInt(x.slice(5), 10); break; }
		const cig = f[5];
		let startClipped = 0, endClipped = 0, m;
		if ((m = /^(?:\d+H)?(\d+)S/.exec(cig))) startClipped = parseInt(m[1], 10);
		if ((m = /(\d+)S(?:\d+H)?$/.exec(cig))) endClipped = parseInt(m[1], 10);
		const isRc = fl & 0x10;
		if (!P.suppressPoor && ((startClipped > P.armClip && endClipped > P.armClip) || left < P.armMin || right < P.armMin)) {
			dPoor += 1;
			if (fd) fd.write(f[0] + '\tPoor alignment\tclipped: [' + startClipped + ',' + endClipped + '] matches: [' + left + ',' + right + ']\n');
			continue;
		}
		const chr = f[2];
		let start = parseInt(f[3], 10), span = 0;
		CIG.lastIndex = 0; while ((m = CIG.exec(cig))) if ('MDNX='.includes(m[2])) span += parseInt(m[1], 10);
		let end = start + span - 1, dir = 'right';
		if (isRc) { [start, end] = [end, start]; [startClipped, endClipped] = [endClipped, startClipped]; dir = 'left'; }
		const alnPos = `${chr}:${start}-${end}`, cutPos = `${chr}:${start}`;
		const delPrimer = origin_seq.length - primerLen, insTarget = startClipped;
		let status = null, dist = null, ori = null;
		const [, r2Q, r2NM] = r2Info.get(f[0].split(' ')[0]) || [null, null, null];
		if (isR1) {
			const o1 = isRc ? '-' : '+', c2 = f[6], o2 = (fl & 0x20) ? '-' : '+';
			dist = Math.abs(parseInt(f[8], 10));
			const mateUnmapped = fl & 0x8;
			if (r2Q !== null && r2Q < P.minQ2) status = notLow;
			else if (chr === '*' || c2 === '*' || mateUnmapped) status = notAln;
			else if (c2 === '=' || chr === c2) status = dist <= P.r1_r2_support_max_distance ? (o1 !== o2 ? sup : notOri) : notDist;
			else status = notChr;
			ori = `${o1}/${o2}`;
			supCounts.inc(status);
			if (status === sup) supDist.inc(dist); else if (status === notDist) notSupDist.inc(dist);
			if (!P.suppressR2 && status !== sup) { dR2 += 1; if (fd) fd.write(f[0] + '\tNot supported by R2\t' + status + '\n'); continue; }
		}
		const barcode = lineInfo.split(':').pop();
		const tlen = Math.abs(parseInt(f[8], 10));
		let key = `b:${barcode} 1:${f[2]}:${f[3]} 2:${f[6]}:${f[7]} t:${tlen} dp:${delPrimer} it:${insTarget}`;
		if (experiment_had_UMIs) {
			if (seen.has(key)) {
				seenCounts.set(key, seenCounts.get(key) + 1);
				if (!(P.suppressDedup || P.dedupFinal)) { dDup += 1; if (fd) fd.write(f[0] + '\tDuplicate\tDuplicate of ' + seen.get(key) + '(' + key + ')\n'); continue; }
			} else { seen.set(key, f[0]); seenCounts.set(key, 1); }
		} else key = '(None)';
		cp(chr).inc(start);
		let insertSize = 'None';
		if (isR1) { insertSize = Math.abs(parseInt(f[8], 10)); tlens.inc(insertSize); }
		chrCounts.inc(chr);
		tmp.write([lineInfo, alnPos, cutPos, dir, delPrimer, insTarget, insertSize, status, dist, ori, q, nm, r2Q, r2NM, key].map(str).join('\t') + '\n');
	}
	tmp.close();
	let tmpLines = () => io.linesOf(tmp);
	log('Finished first pass through alignments');

	const W = async (name, text) => { const w = await io.writer(name); w.write(text); w.close(); };
	const plots = [];
	// chr alignments
	const chrKeys = [...chrCounts.keys()].sort(pyCmp);
	await W(R + '.chr_alignments.txt', 'chr\tnumReads\n' + chrKeys.map((k) => k + '\t' + chrCounts.get(k) + '\n').join(''));
	if (!P.suppressPlots) plots.push(new PlotObject(R + '.chr_alignments', 'Genome alignment summary', chrKeys.length === 0 ? '(No reads aligned to the genome)' :
		'Bar plot showing raw read alignment location<br>Note: this plot includes reads that may be descarded at a later point including duplicates or poor alignments.',
	[['Genome alignment summary', R + '.chr_alignments.txt']], 20, { type: 'bar', file: R + '.chr_alignments.txt', title: 'Location of read alignments', ylabel: 'Number of Reads', categorical: true }));
	let tlenPlot = null;
	if (tlens.size > 0) {
		const ks = [...tlens.keys()].sort((a, b) => a - b);
		await W(R + '.insertSizes.txt', 'insertSize\tnumReads\n' + ks.map((k) => k + '\t' + tlens.get(k) + '\n').join(''));
		if (!P.suppressPlots) tlenPlot = new PlotObject(R + '.insertSizes', 'Genome Alignment Insert Size Summary', 'Bar plot showing insert size of reads aligned to genome',
			[['Genome alignment insert size summary', R + '.insertSizes.txt']], 15, { type: 'bar', file: R + '.insertSizes.txt', title: 'Insert size', ylabel: 'Number of Reads' });
	}
	if (tlenPlot) plots.push(tlenPlot);

	// homology at cut points
	const homology = new Map();
	if (!P.suppressHom) {
		log('Analyzing sequence homology at cut points');
		const w = await io.writer(R + '.cut_homology.txt');
		w.write(['guide_seq', 'cut_chr', 'cut_point', 'is_valid_homology_site', 'potential_guide', 'best_score', 'match_direction', 'n_matches', 'n_mismatches', 'n_gaps', 'aln_guide', 'aln_ref'].join('\t') + '\n');
		const mat = A.makeMatrix(2, -1), pad = 30;
		for (const guide of S.guide_sequences) {
			const gU = guide.toUpperCase(), gRc = C.revcomp(gU);
			for (const [chr, pts] of cutPointsByChr) {
				for (const p of [...pts.keys()].sort((a, b) => a - b)) {
					const seq = faidx(genome, chr, p - pad, p + pad).toUpperCase();
					const inc = new Int32Array(seq.length + 1);
					const [f1, f2, fs] = A.globalAlign(gU, seq, mat, inc, -5, -2);
					const [r1, r2, rs] = A.globalAlign(gRc, seq, mat, inc, -5, -2);
					let isRc = false, best = fs, b1 = f1, b2 = f2, nM, nMM, nG, cind, pot;
					if (fs > rs) [nM, nMM, nG, cind, pot] = guideMatchFromAln(f1, f2, guide.length + S.cleavage_offset - 1);
					else { isRc = true; best = rs; b1 = r1; b2 = r2; [nM, nMM, nG, cind, pot] = guideMatchFromAln(r1, r2, -1 * S.cleavage_offset - 1); pot = revcomp(pot); }
					const valid = nG <= P.gaps && nMM <= P.mismatches;
					if (valid) {
						const k = chr + ':' + (p - pad + cind + 1);
						if (!homology.has(k) || best > homology.get(k)[0]) homology.set(k, [best, pot, nM, nMM, nG, isRc]);
					}
					w.write([guide, chr, p, valid, pot, F(best), isRc ? 'Reverse' : 'Forward', nM, nMM, nG, b1, b2].map(str).join('\t') + '\n');
				}
			}
		}
		w.close();
	}

	const col = collapseCutPoints(P.novel, P.known, P.origin, cutPointsByChr, cut_annotations, homology);
	log('Merging cut points'); for (const l of col.log) log(l);
	const { originChr, originCutPos, originDirection, lookup, finalByChr } = col;

	const regions = new Map();
	if (ctx.regionLines) {
		let nl = 0, na = 0;
		for (const l of ctx.regionLines) {
			nl += 1;
			const e = C.strip(l).split('\t');
			if (e.length > 3) {
				const a = parseInt(e[1], 10), b = parseInt(e[2], 10);
				if (/^\s*[+-]?\d+\s*$/.test(e[1]) && /^\s*[+-]?\d+\s*$/.test(e[2])) { if (!regions.has(e[0])) regions.set(e[0], []); regions.get(e[0]).push([a, b, e[3]]); na += 1; }
			}
		}
		log(`Read ${nl} from cut region annotation file. Added ${na}/${nl} annotations.`);
	}
	const annoOf = (chrom, pos, finalKey) => {
		let reg = 'Novel';
		for (const [a, b, n] of regions.get(chrom) || []) if (pos >= a && pos <= b) { reg = n; break; }
		let anno = [reg, 'Not-origin', 'NA'];
		if (homology.has(finalKey)) anno = ['Homologous', 'Not-origin', 'NA'];
		if (cut_annotations.has(finalKey)) anno = cut_annotations.get(finalKey);
		return anno;
	};
	let finalCutTotal = 0;
	const cls = new Map();
	const chroms = [...finalByChr.keys()].sort(pyCmp);
	for (const chrom of chroms) {
		const fp = finalByChr.get(chrom);
		finalCutTotal += fp.size;
		for (const pos of [...fp.keys()].sort((a, b) => a - b)) {
			const anno = annoOf(chrom, pos, chrom + ':' + pos);
			const L = `${chrom}:${pos}:left`, Rt = `${chrom}:${pos}:right`;
			if (originDirection === null) { cls.set(L, anno[0]); cls.set(Rt, anno[0]); }
			else if (chrom !== originChr) { cls.set(L, anno[0] + ' translocation'); cls.set(Rt, anno[0] + ' translocation'); }
			else if (originDirection === 'left') {
				if (pos === originCutPos) { cls.set(Rt, 'Linear'); cls.set(L, 'Chimera'); }
				else { cls.set(Rt, anno[0] + ' large deletion'); cls.set(L, anno[0] + ' large inversion'); }
			} else if (pos === originCutPos) { cls.set(L, 'Linear'); cls.set(Rt, 'Chimera'); }
			else { cls.set(L, anno[0] + ' large deletion'); cls.set(Rt, anno[0] + ' large inversion'); }
		}
	}
	log(`${finalCutTotal} final cut points after collapsing`);
	for (const a of S.cut_classification_annotations) { const e = a.split(':'); cls.set(e.slice(0, 3).join(':'), e[3]); }
	await W(R + '.cut_classification.txt', [...cls.keys()].sort(pyCmp).map((k) => k + '\t' + cls.get(k) + '\n').join(''));
	{
		let t = 'chr\tpos\tcount\tannotation\n';
		for (const chrom of chroms) { const fp = finalByChr.get(chrom); for (const pos of [...fp.keys()].sort((a, b) => a - b)) t += `${chrom}\t${pos}\t${fp.get(pos)}\t${annoOf(chrom, pos, chrom + ':' + pos).join(',')}\n`; }
		await W(R + '.final_cut_points.txt', t);
	}
	{
		const keys = [...lookup.keys()].sort((a, b) => { const [ca, pa] = a.split(':'), [cb, pb] = b.split(':'); return pyCmp(ca, cb) || (parseInt(pa, 10) - parseInt(pb, 10)); });
		let t = 'discovered cut point\tfinal mapped cut point\tcount\tannotation\n';
		for (const k of keys) { const [chrom, ps] = k.split(':'); const pos = parseInt(ps, 10); t += `${k}\t${lookup.get(k)}\t${cutPointsByChr.get(chrom).get(pos)}\t${annoOf(chrom, pos, lookup.get(k)).join(',')}\n`; }
		await W(R + '.final_cut_points_lookup.txt', t);
	}

	// dedup by final cut assignment
	let totalForDedup = totalR1;
	if (experiment_had_UMIs && P.dedupFinal) {
		log('Deduplicating based on final cut assignments');
		const seen2 = new Map(); seenCounts = new Map();
		let readN = 0, printed = 0; dDup = 0;
		const dt = await io.writer(R + '.final_assignments.dedup.tmp');
		let first = true;
		for (const l of tmpLines()) {
			if (first) { dt.write(l + '\n'); first = false; continue; }
			readN += 1;
			const e = l.split('\t');
			const fcp = lookup.get(e[2]);
			const k = `${e[0].split(':').pop()} ${fcp} ${e[3]}`;
			if (seen2.has(k)) { seenCounts.set(k, seenCounts.get(k) + 1); dDup += 1; if (fd) fd.write(e[0] + '\tDuplicate\tDuplicate of ' + seen2.get(k) + '(' + k + ')\n'); continue; }
			seen2.set(k, e[0]); seenCounts.set(k, 1);
			e[14] = k; dt.write(e.join('\t') + '\n'); printed += 1;
		}
		dt.close(); tmpLines = () => io.linesOf(dt);
		totalForDedup = readN;
		log(`Printed ${printed}/${readN} reads after deduplication based on final cut point assignment`);
	}
	if (fd) fd.close();

	// final assignments
	const finalName = R + '.final_assignments.txt';
	let finalTotal = 0;
	const clsCounts = new DefaultInt(), clsIndelCounts = new DefaultInt(), readIds = new Map(), cutCounts = new DefaultInt();
	const oppositeDir = originDirection === 'left' ? 'right' : 'left';
	const indelLens = new DefaultInt(), invLens = new DefaultInt(), delLens = new DefaultInt();
	const LONG = 100000, oLen = origin_seq.length;
	const h100 = new Float64Array(LONG + 2), h100p = new Float64Array(oLen + 2), hL = new Float64Array(LONG + 2), hLp = new Float64Array(oLen + 2); // histograms of covered lengths
	let t100 = 0, tL = 0;
	const dedupSum = new DefaultInt(), dedupCnt = new DefaultInt();
	const fw = await io.writer(finalName);
	let first = true;
	for (const l of tmpLines()) {
		if (first) {
			const h = C.strip(l).split('\t');
			fw.write([...h.slice(0, 15), 'reads_with_same_umi_pos', 'final_cut_pos', 'final_cut_indel', 'indel_classification', 'classification'].join('\t') + '\n');
			first = false; continue;
		}
		finalTotal += 1;
		const e = l.split('\t');
		const fcp = lookup.get(e[2]), dir = e[3], fcd = fcp + ':' + dir;
		const [fchr, fps] = fcp.split(':'); const fstart = parseInt(fps, 10);
		const alnStart = parseInt(e[1].split(':')[1].split('-')[0], 10);
		const alnIndel = dir === 'right' ? -1 * (fstart - alnStart) : fstart - alnStart - 1;
		const delPrimer = parseInt(e[4], 10), insTarget = parseInt(e[5], 10);
		const indel = insTarget - delPrimer - alnIndel;
		cutCounts.inc(fcd);
		if (Math.abs(indel) <= P.crispresso_max_indel_size) readIds.set(e[0], fcd);
		if (!cls.has(fcd)) throw new Error(`KeyError: '${fcd}'`);
		const fc = cls.get(fcd);
		clsCounts.inc(fc);
		let indelStr = 'no indels';
		if (Math.abs(indel) > 0) indelStr = Math.abs(indel) > P.short_indel_length_cutoff ? 'long indels' : 'short indels';
		clsIndelCounts.inc(fc + ' ' + indelStr);
		let nSeen = 1;
		if (experiment_had_UMIs) { nSeen = seenCounts.get(e[14]); if (nSeen === undefined) throw new Error(`KeyError: '${e[14]}'`); dedupSum.inc(fc + ' ' + indelStr, nSeen); }
		else dedupSum.inc(fc + ' ' + indelStr, 1);
		dedupCnt.inc(fc + ' ' + indelStr);
		fw.write([...e.slice(0, 15), nSeen, fcd, indel, indelStr, fc].join('\t') + '\n');
		if (fchr === originChr) {
			if (fstart === originCutPos && dir === oppositeDir) indelLens.inc(indel);
			const d = alnStart - originCutPos;
			if (dir === originDirection) invLens.inc(d); else delLens.inc(d);
			let dl = null;
			if (originDirection === 'left' && dir === 'right' && alnStart >= originCutPos) dl = alnStart - originCutPos;
			else if (originDirection === 'right' && dir === 'left' && alnStart <= originCutPos) dl = (alnStart - originCutPos) * -1;
			if (dl !== null) {
				const dp = Math.min(delPrimer, oLen + 1);
				if (dl <= 100) { t100 += 1; h100[Math.min(dl, LONG + 1)] += 1; h100p[dp] += 1; }
				if (dl <= LONG) { tL += 1; hL[Math.min(dl, LONG + 1)] += 1; hLp[dp] += 1; }
			}
		}
	}
	fw.close();
	if (!ctx.keepIntermediate) { await io.remove(tmpName); await io.remove(R + '.final_assignments.dedup.tmp'); }
	log(`Processed ${finalTotal} reads`);

	// classification labels
	let clsLabels = originChr !== null ? ['Linear'] : [];
	for (const k of [...clsCounts.keys()].sort(pyCmp)) if (!clsLabels.includes(k) && clsCounts.get(k) > 0) clsLabels.push(k);
	const indelCats = [' no indels', ' short indels', ' long indels'];
	let clsIndelLabels = [];
	for (const l of clsLabels) for (const c of indelCats) clsIndelLabels.push(l + c);

	// deduplication summary
	const dupCounts = new DefaultInt();
	for (const v of seenCounts.values()) dupCounts.inc(v);
	const dupKeys = [...dupCounts.keys()].sort((a, b) => a - b);
	let g = null, gStr = '#Final UMI GINI:\tNA\n';
	if (seenCounts.size > 0) { g = U.gini([...seenCounts.values()]); gStr = '#Final UMI GINI:\t' + C.pyFloat(g) + '\n'; }
	if (!experiment_had_UMIs) gStr = '#Final UMI GINI:\tNA\t(UMIs not used in this experiment)';
	const dedupRoot = R + (P.dedupFinal ? '.deduplication_by_UMI_and_final_cut_assignment' : '.deduplication_by_UMI_and_aln_pos');
	await W(dedupRoot + '.txt', gStr + 'Reads per UMI\tNumber of UMIs\n' + dupKeys.map((k) => k + '\t' + dupCounts.get(k) + '\n').join(''));
	await W(dedupRoot + '.reads_per_umi_by_classification.txt', 'Classification\tNumber of reads\tNumber of UMIs\tAverage reads/umi\n' + clsIndelLabels.map((c) => {
		const n = dedupCnt.get(c), s = dedupSum.get(c);
		return c + '\t' + s + '\t' + n + '\t' + (n > 0 ? C.pyFloat(s / n) : '0') + '\n';
	}).join(''));
	let dedupPlot = null;
	if (totalReads > 0 && !P.suppressPlots) {
		const gNote = 'Final UMI GINI score: ' + (g === null ? 'None' : C.pyFloat(g)) + ' (0 = perfect equality and 1 = perfect inequality in reads/UMI).';
		let label = 'Number of reads that were duplicates based on alignment position and UMI. ' + gNote + ' Note that deduplication of reads by alignment position and UMI is performed by default. To disable deduplication, set the parameter --suppress_dedup_on_aln_pos_and_UMI_filter.';
		if (P.suppressDedup) label = 'Number of reads that were duplicates based on alignment position and UMI. ' + gNote + ' Deduplication based on alignment and UMI position is not performed because the flag --suppress_dedup_on_aln_pos_and_UMI_filter is set.';
		if (P.dedupFinal) label = 'Number of reads that were duplicates based on UMI and final cut assignment. ' + gNote + " Note that deduplication was performed by removing reads with the same UMI and final cut assignment (because the '--dedup_by_final_cut_assignment_and_UMI' flag was set). This may result in deduplication of reads that have the same UMI but different alignment positions.";
		if (!experiment_had_UMIs) label = '(No UMIs used in this experiment)';
		dedupPlot = new PlotObject(dedupRoot, 'Duplicate read counts', label, [['UMI counts', dedupRoot + '.txt'], ['Read assignments', finalName]], 16,
			{ type: 'dedup', pie: [['Not duplicate', totalForDedup - dDup], ['Duplicate', dDup]], bars: dupKeys.map((k) => [k, dupCounts.get(k)]) });
	}
	if (tlenPlot) plots.push(dedupPlot); // (sic) CRISPRuno adds the dedup plot only when the insert size plot exists

	// translocation lists
	const txKeys = [...cutCounts.keys()].sort(pyCmp);
	{
		let t = 'cut\tcut_annotation\tcount_left\tcount_right\tcount_total\n';
		const written = new Set();
		for (const k of txKeys) {
			const cut = k.split(':').slice(0, 2).join(':');
			if (written.has(cut)) continue;
			const ann = cut_annotations.has(cut) ? cut_annotations.get(cut)[0] : 'Novel';
			const l = cutCounts.get(cut + ':left'), r = cutCounts.get(cut + ':right');
			t += `${cut}\t${ann}\t${l}\t${r}\t${l + r}\n`; written.add(cut);
		}
		await W(R + '.translocation_list.txt', t);
	}
	const fragName = R + '.fragment_translocation_list.txt';
	await W(fragName, 'fragment\tcut_annotation\tfragment_annotation\tcount\n' + txKeys.map((k) => {
		const cut = k.split(':').slice(0, 2).join(':');
		return `${k}\t${cut_annotations.has(cut) ? cut_annotations.get(cut)[0] : 'Novel'}\t${cls.has(k) ? cls.get(k) : 'Unknown'}\t${cutCounts.get(k)}\n`;
	}).join(''));
	const sortedTx = [...cutCounts.keys()].sort((a, b) => cutCounts.get(b) - cutCounts.get(a));
	const top = sortedTx.slice(0, 20);
	let otherN = 0, otherReads = 0;
	for (let i = 20; i < sortedTx.length; ++i) if (cutCounts.get(sortedTx[i]) > 0) { otherN += 1; otherReads += cutCounts.get(sortedTx[i]); }

	const clsValues = clsLabels.map((x) => clsCounts.get(x));
	const clsStr = clsLabels.join('\t') + '\n' + clsValues.join('\t') + '\n';
	await W(R + '.classifications.txt', clsStr);
	clsIndelLabels = []; const clsIndelValues = [];
	const inner = [], outer = []; let innerOther = 0; const outerOther = [0, 0, 0];
	const sumInner = clsValues.reduce((a, b) => a + b, 0), indelColors = ['silver', 'pink', 'crimson'];
	for (const l of clsLabels) {
		const these = indelCats.map((c) => { clsIndelLabels.push(l + c); const v = clsIndelCounts.get(l + c); clsIndelValues.push(v); return v; });
		const ts = these.reduce((a, b) => a + b, 0);
		if (sumInner > 0) {
			if (ts / sumInner > 0.05) { inner.push([l, clsCounts.get(l)]); these.forEach((v, i) => { if (v > 0) outer.push([v, indelColors[i]]); }); }
			else { innerOther += ts; these.forEach((v, i) => { outerOther[i] += v; }); }
		}
	}
	if (innerOther > 0) { inner.push(['Other', innerOther]); outerOther.forEach((v, i) => { if (v > 0) outer.push([v, indelColors[i]]); }); }
	const clsIndelStr = clsIndelLabels.join('\t') + '\n' + clsIndelValues.join('\t') + '\n';
	await W(R + '.classifications_with_indels.txt', clsIndelStr);
	const truncated = (n) => (n < cutCounts.size ? '. Note that only the top ' + n + '/' + cutCounts.size + ' events are shown. All events and counts can be found in the Fragment translocation list.' : '');
	if (top.length > 0 && !P.suppressPlots) {
		const bars = top.map((k) => {
			const cut = k.split(':').slice(0, 2).join(':');
			let a = cut_annotations.has(cut) ? cut_annotations.get(cut)[0] : 'Novel';
			if (a.includes('Cas-OFFinder')) a = 'Cas-OFFinder';
			return { left: cls.get(k), right: k, count: cutCounts.get(k), cut, anno: a };
		});
		if (otherN > 0) bars.push({ left: 'Other', right: '(' + otherN + ' locations)', count: otherReads, other: true });
		plots.push(new PlotObject(R + '.fragment_translocation_list', 'Translocation Summary', 'Final translocation counts' + truncated(sortedTx.length), [['Fragment translocation list', fragName]], 38, { type: 'txorder', bars }));
		const heights = new Map(); let maxCount = 0;
		for (const k of top) {
			const e = k.split(':'), cut = e[0] + ':' + e[1];
			if (!heights.has(e[0])) heights.set(e[0], new Map());
			heights.get(e[0]).set(parseInt(e[1], 10), { left: cutCounts.get(cut + ':left'), right: cutCounts.get(cut + ':right') });
			if (cutCounts.get(k) > maxCount) maxCount = cutCounts.get(k);
		}
		plots.push(new PlotObject(R + '.fragment_translocation_list.counts', 'Translocation Counts', 'Final translocation positions and counts' + truncated(top.length), [['Fragment translocation list', fragName]], 40,
			{ type: 'txcounts', chroms: [...heights].map(([c, m]) => [c, [...m].sort((a, b) => a[0] - b[0])]), maxCount, originChr, originCutPos, originDirection }));
		const tx = new DefaultInt(); for (const [k, v] of cutCounts) { const e = k.split(':'); tx.inc(e[0] + ':' + e[1], v); }
		plots.push(new PlotObject(R + '.fragment_translocation_list.circos', 'Translocation Circos Plot', 'Final translocation positions and counts', [['Fragment translocation list', fragName]], 41,
			{ type: 'circos', sites: [...tx].filter(([, v]) => v > 0).sort((a, b) => b[1] - a[1]).slice(0, 100), originChr, originCutPos }));
		plots.push(new PlotObject(R + '.classifications', 'Read Classification', 'Read classification<br>' + clsLabels.map((l, i) => `${l} N=${clsValues[i]}`).join('<br>'),
			[['Alignment classifications', R + '.classifications.txt'], ['Read assignments', finalName]], 36, { type: 'pie', items: clsLabels.map((l, i) => [l, clsValues[i]]), title: 'Read classifications', all: true }));
		plots.push(new PlotObject(R + '.classifications_with_indels', 'Read Classification including indels',
			'Read classification including annotations for presence of indels<br>(Short indels: <=50bp, Long indels: >50bp)<br>' + clsIndelLabels.map((l, i) => `${l} N=${clsIndelValues[i]}`).join('<br>'),
			[['Alignment classifications with indels', R + '.classifications_with_indels.txt'], ['Alignment classifications', R + '.classifications.txt'], ['Read assignments', finalName]], 2,
			{ type: 'donut', inner, outer, title: 'Read classifications' }));
	}

	// origin histograms
	const ik = [...indelLens.keys()].sort((a, b) => a - b);
	await W(R + '.origin_indel_histogram.txt', 'Indel_length\tnumReads\n' + ik.map((k) => k + '\t' + indelLens.get(k) + '\n').join(''));
	if (indelLens.size > 0 && !P.suppressPlots) plots.push(new PlotObject(R + '.origin_indel_histogram', 'Indel length at origin',
		'Bar plot showing lengths of deletions (negative) or insertions (positive) for linear reads supporting the origin cut site at ' + originChr + ': ' + originCutPos,
		[['Indel length at origin histogram', R + '.origin_indel_histogram.txt']], 25, { type: 'bar', file: R + '.origin_indel_histogram.txt', title: 'Indel length at origin', ylabel: 'Number of Reads' }));
	const invK = [...invLens.keys()].sort((a, b) => a - b), delK = [...delLens.keys()].sort((a, b) => a - b);
	const idRoot = R + '.origin_inversion_deletion_histogram';
	await W(idRoot + '.inversions.txt', 'Inversion_distance\tnumReads\n' + invK.map((k) => k + '\t' + invLens.get(k) + '\n').join(''));
	await W(idRoot + '.deletions.txt', 'Deletion_distance\tnumReads\n' + delK.map((k) => k + '\t' + delLens.get(k) + '\n').join(''));
	if (!P.suppressPlots && (delLens.size > 0 || invLens.size > 0)) { // CRISPRuno makes this plot but never adds it to the report (variable name typo); shown here
		plots.push(new PlotObject(idRoot, 'Inversion/Deletion distance relative to origin cut',
			(invK.length ? 'Bar plot showing distance from origin cut of reads supporting inversions (up to 100kb). The origin cut is at ' + originChr + ':' + originCutPos + ' and the origin extends ' + originDirection + ' from the cut to the primer. Reads supporting inversions extend ' + originDirection + ' after genomic alignment. Distances greater than 100000bp are not shown.' : '(No reads supporting inversions found)') + '<br/>' +
			'Bar plot showing distance from origin cut of reads supporting deletions (up to 100kb). The origin cut is at ' + originChr + ':' + originCutPos + ' and the origin extends ' + originDirection + ' from the cut to the primer. Reads supporting deletions extend ' + oppositeDir + ' after genomic alignment. Distances greater than 100000bp are not shown.',
			[['Inversion distance relative to origin cut (all values)', idRoot + '.inversions.txt'], ['Deletion distance relative to origin cut (all values)', idRoot + '.deletions.txt']], 26,
			{ type: 'invdel', files: [idRoot + '.inversions.txt', idRoot + '.deletions.txt'] }));
	}
	// read depth around the origin: depth[i] = total - #reads with deletion length >= i (numpy float arrays in CRISPRuno)
	const depthText = (hist, histP, total, n) => {
		const cover = (h, len) => { const c = new Float64Array(len); let acc = 0; for (let i = h.length - 1; i >= 0; --i) { acc += h[i]; if (i < len) c[i] = acc; } return c; };
		const post = cover(hist, n + 1), pri = cover(histP, oLen + 1);
		let t = 'bpFromOriginCut\treadDepth\n';
		for (let x = oLen; x >= 1; --x) t += -x + '\t' + C.pyFloat(total - pri[x]) + '\n';
		for (let x = 1; x <= n; ++x) t += x + '\t' + C.pyFloat(total - post[x]) + '\n';
		return t;
	};
	const dRoot = R + '.origin_depth_counts';
	await W(dRoot + '_100bp.txt', depthText(h100, h100p, t100, 100));
	await W(dRoot + '_100kb.txt', depthText(hL, hLp, tL, LONG));
	if (!P.suppressPlots) {
		const primerInfo = originChr !== null ? 'Origin cut is at ' + originChr + ':' + originCutPos + ' and the origin extends ' + originDirection + ' from the cut to the primer. The vertical red dotted line shows the origin location. ' : 'The vertical red dotted line shows the end of the primer sequence. ';
		plots.push(new PlotObject(dRoot, 'Read depth around origin', 'Read depth surrounding origin. The top plot shows reads with deletions smaller than 100bp (N=' + t100 + '). The bottom plot shows reads with deletions smaller than 100kb (N=' + tL + '). ' + primerInfo +
			'The vertical black dotted line shows the minimum requred primer length (' + S.min_primer_length + 'bp) set by the parameter --min_primer_length.',
		[['Depth around origin (100bp)', dRoot + '_100bp.txt'], ['Depth around origin (100kb)', dRoot + '_100kb.txt']], 28, { type: 'depth', files: [dRoot + '_100bp.txt', dRoot + '_100kb.txt'], minPrimerLength: S.min_primer_length }));
	}
	// R1/R2 support
	for (const k of supCounts.keys()) if (!supLabels.includes(k)) throw new Error('R1/R2 support label ' + k + ' not found');
	const supVals = supLabels.map((x) => supCounts.peek(x));
	await W(R + '.r1_r2_support.txt', supLabels.join('\t') + '\n' + supVals.join('\t') + '\n');
	if (totalReads > 0 && !P.suppressPlots) {
		const ds = P.suppressR2 ? '(Reads without R1/R2 support are currently included in the final analyses because the --suppress_r2_support_filter flag is set.)' : '(Reads without R1/R2 support are currently excluded from the final analyses. Set the --suppress_r2_support_filter flag to include these reads.)';
		plots.push(new PlotObject(R + '.r1_r2_support', 'R1/R2 Support Classification', 'R1/R2 support classification<br>' + supLabels.map((l, i) => `${l} N=${supVals[i]}`).join('<br>') + '<br>' + ds,
			[['R1/R2 support classifications', R + '.r1_r2_support.txt']], 29, { type: 'pie', items: supLabels.map((l, i) => [l, supVals[i]]), title: 'R1/R2 support status' }));
	}
	const dsk = [...supDist.keys()].sort((a, b) => a - b), dnk = [...notSupDist.keys()].sort((a, b) => a - b);
	await W(R + '.r1_r2_distances_supported.txt', 'R1_R2_distance\tnumReads\n' + dsk.map((k) => k + '\t' + supDist.get(k) + '\n').join(''));
	await W(R + '.r1_r2_distances_not_supported.txt', 'R1_R2_distance\tnumReads\n' + dnk.map((k) => k + '\t' + notSupDist.get(k) + '\n').join(''));
	if (totalReads > 0 && !P.suppressPlots) plots.push(new PlotObject(R + '.r1_r2_distances', 'Distance between R1/R2 Reads',
		dsk.length === 0 && dnk.length === 0 ? '(No R1/R2 reads)' : '(top) Bar plot showing distance between R1 and R2 reads for reads with R2 supporting the R1 alignment.<br>(bottom) Bar plot showing distance between R1 and R2 reads for reads with R2 NOT supporting the R1 alignment due to incorrect read orientation or exceeding maximum alignment distance (max distance for support is ' + P.r1_r2_support_max_distance + 'bp)',
		[['R1/R2 supported distance histogram', R + '.r1_r2_distances_supported.txt'], ['R1/R2 not supported distance histogram', R + '.r1_r2_distances_not_supported.txt']], 30,
		{ type: 'bar2', files: [R + '.r1_r2_distances_supported.txt', R + '.r1_r2_distances_not_supported.txt'], titles: ['Distance between R1/R2 for read pairs in which R2 supports the R1 alignment', 'Distance between R1/R2 for read pairs in which R2 does not support the R1 alignment (but R1 and R2 are on the same chromosome)'] }));
	// discarded
	const dLabels = ['Unaligned', 'Duplicate', 'Poor alignment', 'Not supported by R2'], dVals = [dUnaligned, dDup, dPoor, dR2];
	await W(R + '.discarded_reads.txt', dLabels.join('\t') + '\n' + dVals.join('\t') + '\n');
	let discPlot = null;
	if (!P.suppressPlots) {
		const nDisc = dVals.reduce((a, b) => a + b, 0);
		let dupS = P.suppressDedup ? '(Duplicate reads based on UMI/alignment position are currently included in the final analyses because the --suppress_dedup_on_aln_pos_and_UMI_filter flag is set.)' : '(Reads that may be duplicates based on UMI/alignment position are currently excluded from the final analyses. Set the --suppress_dedup_on_aln_pos_and_UMI_filter flag to include these reads.)';
		if (P.dedupFinal) dupS = "(Reads that may be duplicates based on UMI/cut assignment are currently excluded from the final analyses because the '--dedup_by_final_cut_assignment_and_UMI' flag is set. This may result in deduplication of reads that have the same UMI but different alignment positions.)";
		const poorS = P.suppressPoor ? '(Reads with poor alignment are currently included in the final analyses because the --suppress_poor_alignment_filter flag is set.)' : '(Reads with poor alignment are currently excluded from the final analyses. Set the --suppress_poor_alignment_filter flag to include these reads.)';
		const r2S = P.suppressR2 ? '(Reads without R1/R2 support are currently included in the final analyses because the --suppress_r2_support_filter flag is set.)' : '(Reads without R1/R2 support are currently excluded from the final analyses. Set the --suppress_r2_support_filter flag to include these reads.)';
		discPlot = new PlotObject(R + '.discarded_reads', 'Quality filtering counts', nDisc + '/' + totalR1 + ' reads were discarded.<br>' + finalTotal + ' R1 reads were used for the final analysis. <br>' + dLabels.map((l, i) => `${l} N=${dVals[i]}`).join('<br>') + '<br>' +
			dupS + '<br>' + poorS + '<br>' + r2S + '<br>Poor alignments are alignments of reads such that the front AND back of the read have greater than ' + P.armClip + 'bp clipped, or the start or the end of the read have fewer than ' + P.armMin + 'bp matching the reference. These cutoffs can be adjusted using the --arm_max_clipped_bases and --arm_min_matched_start_bases parameters.',
		[['Discarded reads summary', R + '.discarded_reads.txt']], 3, { type: 'pie', items: dLabels.map((l, i) => [l, dVals[i]]) });
		if (P.writeDiscarded) discPlot.datas.push(['Discarded reads', discardedName]);
		plots.push(discPlot);
	}
	const infoName = R + '.info';
	const plotJ = (p) => (p ? p.toJson() : 'None');
	await W(infoName, ['total_reads_processed', 'total_multialignments_skipped', 'total_secondary_multialignments_used', 'total_r1_processed', 'discarded_reads_unaligned', 'discarded_reads_duplicates',
		'discarded_reads_not_supported_by_R2', 'discarded_reads_poor_alignment', 'final_read_count', 'discovered_cut_point_count', 'final_cut_point_count'].join('\t') + '\n' +
		[totalReads, mmSkipped, mmSecUsed, totalR1, dUnaligned, dDup, dR2, dPoor, finalTotal, col.found, finalCutTotal].join('\t') + '\n' +
		clsStr + clsIndelStr + pyJson(cls) + '\n' + JSON.stringify(plots.filter(Boolean).map((p) => ({ name: p.name, title: p.title, label: p.label, datas: p.datas, order: p.order, plot: p.plot }))) + '\n');
	return {
		finalName, readIds, cutCounts, cls, finalTotal,
		discarded: [['Unaligned', dUnaligned], ['Duplicate read', dDup], ['Not supported by R2', dR2], ['Bad alignment', dPoor]],
		clsCounts: clsLabels.map((l, i) => [l, clsValues[i]]), clsIndelCounts: clsIndelLabels.map((l, i) => [l, clsIndelValues[i]]),
		plots: plots.filter(Boolean), originChr, originCutPos, originDirection,
	};
}

const api = { makeFinalReadAssignments, collapseCutPoints, leftRightMatchesMDZ, guideMatchFromAln };
if (typeof module !== 'undefined' && module.exports) module.exports = api;
root.CuAssign = api;
})(typeof self !== 'undefined' ? self : globalThis);
