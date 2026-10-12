// SPDX-License-Identifier: MIT — port of the cutadapt 5.2 3' adapter search (Marcel Martin et al., MIT)
'use strict';
/* What CRISPRuno runs:  cutadapt -a ADAPTER -o out1 [-p out2 --pair-filter=first] --minimum-length N in1 [in2]
 * = one regular 3' adapter (BackAdapter: max error rate 0.1, min overlap 3, indels allowed, no read wildcards,
 *   adapter wildcards only if the adapter has non-ACGT IUPAC codes), trimmed from read 1 only; read 2 is untouched;
 *   pairs/reads are discarded when read 1 is shorter than N afterwards.
 * Ported: kmer_heuristic.create_positions_and_kmers + KmerFinder.kmers_present (prefilter) and _align.Aligner.locate
 * (flags BACK = QUERY_START|QUERY_STOP|REFERENCE_END) with the same tie-breaking and Ukkonen cut-off.
 */
(function (root) {
const MATCH = 1, MISMATCH = -1, INS = -2, DEL = -2;

/* ---------- k-mer prefilter (cutadapt/kmer_heuristic.py) ---------- */
function kmerChunks(seq, chunks) { // set of near-equal chunks (insertion order kept)
	const size = Math.floor(seq.length / chunks), rem = seq.length % chunks;
	const sizes = []; for (let i = 0; i < rem; ++i) sizes.push(size + 1); for (let i = rem; i < chunks; ++i) sizes.push(size);
	const out = new Set(); let off = 0;
	for (const s of sizes) { out.add(seq.slice(off, off + s)); off += s; }
	return out;
}
function backOverlapSearchsets(adapter, minOverlap, errorRate) {
	const L = adapter.length, errorLengths = [], sets = [];
	let maxError = 0;
	for (let i = 0; i <= L; ++i) if (Math.trunc(i * errorRate) > maxError) { errorLengths.push([maxError, i - 1]); maxError += 1; }
	errorLengths.push([maxError, L]);
	let minimumLength = minOverlap;
	for (const [maxErrors, length] of errorLengths) {
		if (minimumLength > length) continue;
		if (maxErrors === 0) {
			const minKmer = 5;
			if (minimumLength < minKmer) {
				for (let i = minimumLength; i < minKmer; ++i) sets.push([-i, null, new Set([adapter.slice(0, i)])]);
				minimumLength = minKmer;
			}
		}
		sets.push([-length, null, kmerChunks(adapter.slice(0, minimumLength), maxErrors + 1)]);
		minimumLength = length + 1;
	}
	return sets;
}
function minimizeKmerSearchList(list) {
	const byKmer = new Map();
	for (const [kmer, start, stop] of list) { if (!byKmer.has(kmer)) byKmer.set(kmer, []); byKmer.get(kmer).push([start, stop]); }
	const out = [];
	for (const [kmer, positions] of byKmer) {
		if (positions.length === 1) { out.push([kmer, positions[0][0], positions[0][1]]); continue; }
		if (positions.some(([a, b]) => a === 0 && b === null)) { out.push([kmer, 0, null]); continue; }
		const front = positions.filter(([a]) => a === 0), back = positions.filter(([, b]) => b === null);
		if (positions.some(([a, b]) => a !== 0 && b !== null)) throw new Error('NotImplementedError: searches starting in the middle');
		if (front.length) out.push([kmer, 0, Math.max(...front.map(([, b]) => b))]);
		if (back.length) out.push([kmer, Math.min(...back.map(([a]) => a)), null]);
	}
	return out;
}
function createPositionsAndKmers(adapter, minOverlap, errorRate, backAdapter, frontAdapter, internal = true) {
	const maxErrors = Math.trunc(adapter.length * errorRate);
	const sets = [];
	if (backAdapter) sets.push(...backOverlapSearchsets(adapter, minOverlap, errorRate));
	if (frontAdapter) {
		const rev = adapter.split('').reverse().join('');
		for (const [start, , ks] of backOverlapSearchsets(rev, minOverlap, errorRate)) sets.push([0, -start, new Set([...ks].map((k) => k.split('').reverse().join('')))]);
	}
	if (internal) sets.push([0, null, kmerChunks(adapter, maxErrors + 1)]);
	const list = [];
	for (const [start, stop, ks] of sets) for (const k of ks) list.push([k, start, stop]);
	const res = new Map();
	for (const [k, start, stop] of minimizeKmerSearchList(list)) { const key = start + ',' + stop; if (!res.has(key)) res.set(key, [start, stop, []]); res.get(key)[2].push(k); }
	return [...res.values()];
}
function kmersPresent(searches, seqUpper) { // KmerFinder.kmers_present (no wildcards: case-insensitive equality)
	const n = seqUpper.length;
	for (const [s0, s1, kmers] of searches) {
		let start = s0, stop = s1 === null ? 0 : s1;
		if (start < 0) { start = n + start; if (start < 0) start = 0; } else if (start > n) continue;
		if (stop < 0) { stop = n + stop; if (stop <= 0) continue; } else if (stop === 0) stop = n;
		if (stop - start <= 0) continue;
		const w = start === 0 && stop === n ? seqUpper : seqUpper.slice(start, stop);
		for (const k of kmers) if (w.includes(k)) return true;
	}
	return false;
}

/* ---------- Aligner.locate (cutadapt/_align.pyx), reference = adapter, query = read ---------- */
class Aligner {
	constructor(reference, maxErrorRate, flags, minOverlap, indelCost = 1) {
		this.ref = reference; this.m = reference.length; this.maxErrorRate = maxErrorRate;
		this.startInRef = !!(flags & 1); this.startInQuery = !!(flags & 2); this.stopInRef = !!(flags & 4); this.stopInQuery = !!(flags & 8);
		this.minOverlap = minOverlap; this.insCost = indelCost; this.delCost = indelCost;
		this.s1 = new Uint8Array(this.m); for (let i = 0; i < this.m; ++i) this.s1[i] = reference.charCodeAt(i);
		this.cost = new Int32Array(this.m + 1); this.score = new Int32Array(this.m + 1); this.origin = new Int32Array(this.m + 1);
	}
	locate(queryUpper) {
		const m = this.m, n = queryUpper.length, s1 = this.s1, maxErr = this.maxErrorRate;
		const C = this.cost, S = this.score, O = this.origin;
		const k = Math.trunc(maxErr * m);
		let maxN = n, minN = 0;
		if (!this.startInQuery) maxN = Math.min(n, m + k);
		if (!this.stopInQuery) minN = Math.max(0, n - m - k);
		const dc = this.delCost, ic = this.insCost;
		if (!this.startInRef && !this.startInQuery) for (let i = 0; i <= m; ++i) { S[i] = i * DEL; C[i] = Math.max(i, minN) * dc; O[i] = 0; }
		else if (this.startInRef && !this.startInQuery) for (let i = 0; i <= m; ++i) { S[i] = 0; C[i] = minN * dc; O[i] = Math.min(0, minN - i); }
		else if (!this.startInRef && this.startInQuery) for (let i = 0; i <= m; ++i) { S[i] = i * DEL; C[i] = i * dc; O[i] = Math.max(0, minN - i); }
		else for (let i = 0; i <= m; ++i) { S[i] = 0; C[i] = Math.min(i, minN) * dc; O[i] = minN - i; }
		let bRefStop = m, bQueryStop = n, bCost = m + n + 1, bOrigin = 0, bScore = 0;
		let last = Math.min(m, k + 1);
		if (this.startInRef) last = m;
		let lastFilled = 0;
		const originInc = this.startInQuery ? 1 : 0, insCostInc = this.startInQuery ? 0 : ic, insScoreInc = this.startInQuery ? 0 : INS;
		const half = Math.trunc(m / 2);
		let cost = 0, origin = 0, score = 0; // function-level like the cdef ints in _align.pyx ('origin' is reused below)
		for (let j = minN + 1; j <= maxN; ++j) {
			let dCost = C[0], dScore = S[0], dOrigin = O[0];
			O[0] += originInc; C[0] += insCostInc; S[0] += insScoreInc;
			const cj = queryUpper.charCodeAt(j - 1);
			for (let i = 1; i <= last; ++i) {
				if (s1[i - 1] === cj) { cost = dCost; origin = dOrigin; score = dScore + MATCH; }
				else {
					const cDiag = dCost + 1, cIns = C[i] + ic, cDel = C[i - 1] + dc;
					if (cDiag <= cDel && cDiag <= cIns) { cost = cDiag; origin = dOrigin; score = dScore + MISMATCH; }
					else if (cDel <= cIns) { cost = cDel; origin = O[i - 1]; score = S[i - 1] + DEL; }
					else { cost = cIns; origin = O[i]; score = S[i] + INS; }
				}
				dCost = C[i]; dScore = S[i]; dOrigin = O[i];
				C[i] = cost; O[i] = origin; S[i] = score;
			}
			lastFilled = last;
			while (last >= 0 && C[last] > k) last -= 1;
			if (last < m) last += 1;
			else if (this.stopInQuery) {
				cost = C[m]; score = S[m]; origin = O[m];
				const length = m + Math.min(origin, 0);
				const ok = length >= this.minOverlap && cost <= length * maxErr;
				const bestLength = m + Math.min(bOrigin, 0);
				if (ok && (bCost === m + n + 1 || (origin <= bOrigin + half && score > bScore) || (length > bestLength && score > bScore))) {
					bScore = score; bCost = cost; bOrigin = origin; bRefStop = m; bQueryStop = j;
					if (cost === 0 && origin >= 0) break;
				}
			}
		}
		if (maxN === n) {
			const firstI = this.stopInRef ? 0 : m;
			for (let i = lastFilled; i >= firstI; --i) {
				const length = i + Math.min(O[i], 0);
				cost = C[i]; score = S[i];
				const ok = length >= this.minOverlap && cost <= length * maxErr;
				const bestLength = bRefStop + Math.min(bOrigin, 0);
				// as in _align.pyx: 'origin' is not updated in this loop (value left over from the column loop)
				if (ok && (bCost === m + n + 1 || (origin <= bOrigin + half && score > bScore) || (length > bestLength && score > bScore))) {
					bScore = score; bCost = cost; bOrigin = O[i]; bRefStop = i; bQueryStop = n;
				}
			}
		}
		if (bCost === m + n + 1) return null;
		let refStart, queryStart;
		if (bOrigin >= 0) { refStart = 0; queryStart = bOrigin; } else { refStart = -bOrigin; queryStart = 0; }
		return [refStart, bRefStop, queryStart, bQueryStop, bScore, bCost];
	}
}

class BackAdapter {
	constructor(sequence, { maxErrors = 0.1, minOverlap = 3 } = {}) {
		this.sequence = sequence.toUpperCase().replace(/U/g, 'T').replace(/I/g, 'N');
		if (!this.sequence) throw new Error('Adapter sequence is empty');
		if (!/^[ACGT]+$/.test(this.sequence)) throw new Error('adapters with IUPAC wildcards are not supported by this port (use A, C, G, T only)');
		this.maxErrorRate = maxErrors; this.minOverlap = Math.min(minOverlap, this.sequence.length);
		this.aligner = new Aligner(this.sequence, this.maxErrorRate, 2 | 8 | 4, this.minOverlap, 1);
		this.kmers = createPositionsAndKmers(this.sequence, this.minOverlap, this.maxErrorRate, true, false);
	}
	matchTo(seq) { // -> {rstart, rstop, astart, astop, score, errors} or null
		const up = seq.toUpperCase();
		if (!kmersPresent(this.kmers, up)) return null;
		const a = this.aligner.locate(up);
		if (!a) return null;
		return { astart: a[0], astop: a[1], rstart: a[2], rstop: a[3], score: a[4], errors: a[5] };
	}
	trim(seq) { const m = this.matchTo(seq); return m ? m.rstart : seq.length; } // length of read 1 kept (RemoveAfterMatch.trimmed)
}

const api = { BackAdapter, Aligner, createPositionsAndKmers, kmersPresent };
if (typeof module !== 'undefined' && module.exports) module.exports = api;
root.CuCutadapt = api;
})(typeof self !== 'undefined' ? self : globalThis);
