// SPDX-License-Identifier: LicenseRef-CRISPResso2 — port of CRISPResso2 2.3.4 CRISPResso2Align.pyx / CRISPRessoCOREResources.pyx
// (CRISPResso2 is (c) MGH, academic non-commercial license; see README)
'use strict';
/* global_align(seqj, seqi, matrix, gapIncentive, gapOpen, gapExtend) -> [alignedJ, alignedI, score]
 *   Needleman-Wunsch with affine gaps exactly as the Cython code (same tie-breaking, the "last row / last column
 *   without gap-open penalty" rule, traceback that steps to index 0), score = round(100*matches/aln_len, 3).
 * readMatrix(text) (EDNAFULL), makeMatrix(match, mismatch, nMismatch, nMatch)
 * findIndelsSubstitutions(readAln, refAln, includeIdxs) -> payload (CRISPRessoCOREResources.find_indels_substitutions)
 * pyRound(x, nd) == Python round(float, nd)
 */
(function (root) {
const MARRAY = 1, IARRAY = 2, JARRAY = 3;

function pyRound(x, nd) { // CPython float round(): exact decimal rounding (half-even on exact ties), then nearest double
	if (!Number.isFinite(x) || x === 0) return x;
	const buf = new DataView(new ArrayBuffer(8)); buf.setFloat64(0, x);
	const hi = buf.getUint32(0), lo = buf.getUint32(4);
	const neg = hi >>> 31, ex = (hi >>> 20) & 0x7ff;
	let mant = (BigInt(hi & 0xfffff) << 32n) | BigInt(lo), e;
	if (ex === 0) e = -1074; else { mant |= 1n << 52n; e = ex - 1075; }
	const p10 = 10n ** BigInt(nd);
	let num = mant * p10, den = 1n; // |x| * 10^nd = num / den
	if (e >= 0) num <<= BigInt(e); else den <<= BigInt(-e);
	let q = num / den; const r = num - q * den;
	if (2n * r > den || (2n * r === den && (q & 1n))) q += 1n;
	const v = Number(q) / Number(p10);
	return neg ? -v : v;
}

function makeMatrix(match = 5, mismatch = -4, nMismatch = -2, nMatch = -1) {
	const size = 'T'.charCodeAt(0) + 1, a = new Int32Array(128 * 128); // indices beyond size read 0 here (C reads past the array)
	const nuc = ['A', 'T', 'C', 'G'].map((c) => c.charCodeAt(0)), N = 78;
	for (const x of nuc) for (const y of nuc) a[x * 128 + y] = x === y ? match : mismatch;
	for (const x of nuc) { a[x * 128 + N] = nMismatch; a[N * 128 + x] = nMismatch; }
	a[N * 128 + N] = nMatch;
	a.size = size;
	return a;
}
function readMatrix(text) { // NCBI format, mat[ord(row), ord(col)]
	const a = new Int32Array(128 * 128);
	const lines = text.split('\n');
	let li = 0, headers = null;
	while (headers === null) { const l = lines[li++].trim(); if (l[0] === '#') continue; headers = l.split(' ').filter(Boolean).map((x) => x.charCodeAt(0)); }
	let ai = 0;
	for (; li < lines.length; ++li) {
		const l = lines[li]; if (!l) continue;
		const vals = l.split(' ').slice(1).filter(Boolean).map((x) => parseInt(x, 10));
		for (let k = 0; k < Math.min(vals.length, headers.length); ++k) a[headers[ai] * 128 + headers[k]] = vals[k];
		ai += 1;
	}
	return a;
}

// work arrays reused between calls
let cap = 0, mS, iS, jS, mP, iP, jP;
function ensure(n) {
	if (n <= cap) return;
	cap = Math.max(n, cap * 2);
	mS = new Int32Array(cap); iS = new Int32Array(cap); jS = new Int32Array(cap);
	mP = new Uint8Array(cap); iP = new Uint8Array(cap); jP = new Uint8Array(cap);
}

function globalAlign(seqj, seqi, matrix, gapIncentive, gapOpen = -1, gapExtend = -1) {
	const maxJ = seqj.length, maxI = seqi.length;
	if (gapIncentive.length !== maxI + 1) throw new Error('Mismatch in gap_incentive length');
	const W = maxJ + 1;
	ensure((maxI + 1) * W);
	const sj = new Uint8Array(maxJ + 1), si = new Uint8Array(maxI + 1);
	for (let j = 0; j < maxJ; ++j) sj[j] = seqj.charCodeAt(j) & 0x7f;
	for (let i = 0; i < maxI; ++i) si[i] = seqi.charCodeAt(i) & 0x7f;
	const G = gapIncentive, M = matrix;
	const minScore = gapOpen * maxJ * maxI;
	// init (np.empty cells that are never written stay at whatever we hold; the recurrences below never read them)
	for (let j = 1; j <= maxJ; ++j) { mS[j] = minScore; mP[j] = IARRAY; }
	for (let i = 1; i <= maxI; ++i) { mS[i * W] = minScore; mP[i * W] = JARRAY; }
	mS[0] = 0; mP[0] = 0;
	for (let j = 1; j <= maxJ; ++j) { iS[j] = gapExtend * j + G[0]; iP[j] = IARRAY; }
	for (let i = 0; i <= maxI; ++i) iS[i * W] = minScore;
	for (let i = 1; i <= maxI; ++i) { jS[i * W] = gapExtend * i + G[0]; jP[i * W] = JARRAY; }
	for (let j = 0; j <= maxJ; ++j) jS[j] = minScore;
	for (let i = 0; i <= maxI; ++i) iP[i * W] = 0; // never written in the Cython code (np.empty): 0 makes a (never expected) read fail loudly
	for (let j = 0; j <= maxJ; ++j) jP[j] = 0;

	const cell = (i, j, open) => {
		const ci = si[i - 1], cj = sj[j - 1], c = i * W + j;
		const iFromM = open + mS[c - 1] + G[i], iExt = gapExtend + iS[c - 1] + G[i];
		if (iFromM > iExt) { iS[c] = iFromM; iP[c] = MARRAY; } else { iS[c] = iExt; iP[c] = IARRAY; }
		const jFromM = open + mS[c - W] + G[i - 1], jExt = gapExtend + jS[c - W];
		if (jFromM > jExt) { jS[c] = jFromM; jP[c] = MARRAY; } else { jS[c] = jExt; jP[c] = JARRAY; }
		const s = M[ci * 128 + cj], d = c - W - 1;
		const mVal = mS[d] + s, iVal = iS[d] + s, jVal = jS[d] + s;
		if (mVal > jVal) { if (mVal > iVal) { mS[c] = mVal; mP[c] = MARRAY; } else { mS[c] = iVal; mP[c] = IARRAY; } }
		else if (jVal > iVal) { mS[c] = jVal; mP[c] = JARRAY; } else { mS[c] = iVal; mP[c] = IARRAY; }
	};
	for (let i = 1; i < maxI; ++i) for (let j = 1; j < maxJ; ++j) cell(i, j, gapOpen);
	if (maxJ >= 1) for (let i = 1; i < maxI; ++i) cell(i, maxJ, gapExtend); // last column: no gap open penalty
	if (maxI >= 1) for (let j = 1; j <= maxJ; ++j) cell(maxI, j, gapExtend); // last row

	// traceback
	let i = maxI, j = maxJ;
	let ci = si[i > 0 ? i - 1 : 0], cj = sj[j > 0 ? j - 1 : 0];
	const c0 = i * W + j;
	let cur;
	if (mS[c0] > jS[c0]) cur = mS[c0] > iS[c0] ? MARRAY : IARRAY;
	else cur = jS[c0] > iS[c0] ? JARRAY : IARRAY;
	const aj = [], ai = [];
	let matchCount = 0;
	while (i > 0 || j > 0) {
		const c = i * W + j;
		if (cur === MARRAY) {
			cur = mP[c]; aj.push(cj); ai.push(ci);
			if (cj === ci) matchCount += 1;
			if (i > 1) { i -= 1; ci = si[i - 1]; } else { i = 0; ci = si[0]; }
			if (j > 1) { j -= 1; cj = sj[j - 1]; } else { j = 0; cj = sj[0]; }
		} else if (cur === JARRAY) {
			cur = jP[c]; aj.push(45); ai.push(ci);
			if (i > 1) { i -= 1; ci = si[i - 1]; } else { i = 0; ci = si[0]; }
		} else if (cur === IARRAY) {
			cur = iP[c]; aj.push(cj); ai.push(45);
			if (j > 1) { j -= 1; cj = sj[j - 1]; } else { j = 0; cj = sj[0]; }
		} else throw new Error('wtf4!:pointer: ' + i);
	}
	const n = aj.length;
	let sJ = '', sI = '';
	for (let k = n - 1; k >= 0; --k) { sJ += String.fromCharCode(aj[k]); sI += String.fromCharCode(ai[k]); }
	return [sJ, sI, pyRound(100 * matchCount / n, 3)];
}

/* CRISPRessoCOREResources.find_indels_substitutions (non-legacy) */
function findIndelsSubstitutions(read, ref, includeIdxs) {
	const inc = includeIdxs instanceof Set ? includeIdxs : new Set(includeIdxs);
	const refPositions = [], allSubPos = [], subPos = [], allSubVal = [], subVal = [];
	const allDelPos = [], allDelCoords = [], delPos = [], delCoords = [], delSizes = [];
	const allInsPos = [], allInsLeft = [], insPos = [], insCoords = [], insSizes = [];
	let startDel = -1, startIns = -1, idx = 0, curInsSize = 0;
	const L = ref.length;
	const anyIn = (a, b) => { for (let x = a; x < b; ++x) if (inc.has(x)) return true; return false; };
	for (let ic = 0; ic < L; ++ic) {
		const c = ref[ic], r = read[ic];
		if (c !== '-') {
			refPositions.push(idx);
			if (c !== r && r !== '-' && r !== 'N') {
				allSubPos.push(idx); allSubVal.push(r);
				if (inc.has(idx)) { subPos.push(idx); subVal.push(r); }
			}
			if (startIns !== -1) {
				allInsLeft.push(startIns); allInsPos.push(startIns, idx);
				if (inc.has(startIns) && inc.has(idx)) { insCoords.push([startIns, idx]); insPos.push(startIns, idx); insSizes.push(curInsSize); }
				startIns = -1;
			}
			curInsSize = 0;
			idx += 1;
		} else {
			refPositions.push(idx === 0 ? -1 : -idx);
			if (idx > 0 && startIns === -1) startIns = idx - 1;
			curInsSize += 1;
		}
		if (r === '-' && startDel === -1) startDel = ic - 1 >= 0 ? refPositions[ic] : 0;
		else if (r !== '-' && startDel !== -1) {
			const endDel = refPositions[ic];
			for (let x = startDel; x < endDel; ++x) allDelPos.push(x);
			allDelCoords.push([startDel, endDel]);
			if (anyIn(startDel, endDel)) { for (let x = startDel; x < endDel; ++x) delPos.push(x); delCoords.push([startDel, endDel]); delSizes.push(endDel - startDel); }
			startDel = -1;
		}
	}
	if (startDel !== -1) {
		const endDel = refPositions[L - 1];
		for (let x = startDel; x < endDel + 1; ++x) allDelPos.push(x);
		allDelCoords.push([startDel, endDel + 1]);
		if (anyIn(startDel, endDel + 1)) { for (let x = startDel; x < endDel + 1; ++x) delPos.push(x); delCoords.push([startDel, endDel + 1]); delSizes.push(endDel + 1 - startDel); }
	}
	const sum = (a) => a.reduce((s, v) => s + v, 0);
	return {
		all_insertion_positions: allInsPos, all_insertion_left_positions: allInsLeft, insertion_positions: insPos, insertion_coordinates: insCoords,
		insertion_sizes: insSizes, insertion_n: sum(insSizes),
		all_deletion_positions: allDelPos, all_deletion_coordinates: allDelCoords, deletion_positions: delPos, deletion_coordinates: delCoords,
		deletion_sizes: delSizes, deletion_n: sum(delSizes),
		all_substitution_positions: allSubPos, substitution_positions: subPos, all_substitution_values: allSubVal, substitution_values: subVal,
		substitution_n: subPos.length, ref_positions: refPositions,
	};
}

const api = { globalAlign, makeMatrix, readMatrix, findIndelsSubstitutions, pyRound };
if (typeof module !== 'undefined' && module.exports) module.exports = api;
root.CuAlign = api;
})(typeof self !== 'undefined' ? self : globalThis);
