// SPDX-License-Identifier: MIT — SAM post-processing for the CRISPRuno web port
'use strict';
/* Equivalent of the end of CRISPRuno's align_reads():  aligner | samtools view [-F 256] -Shu - | samtools sort
 *   - records with flag 0x800 (supplementary) are dropped: bowtie2 never writes them, minibwa does
 *   - secondary records (0x100) are dropped unless keepSecondary (analyze_multimap_assignments, bowtie2 -k N)
 *   - RNEXT equal to RNAME is printed as '=' (as samtools view does)
 *   - samtools sort order: reference index (@SQ order, unmapped last), POS, reverse strand; ties keep input order
 * Large inputs: records are first appended to one temporary file per reference ("bucket"), each bucket is then
 * loaded and sorted in memory (memory ~ size of the largest bucket).
 */
(function (root) {
const IO = root.GsIO || require('../lib/textio.js');
const enc = new TextEncoder();

function headerRefs(headerText) {
	const tid = new Map();
	for (const l of headerText.split('\n')) if (l.startsWith('@SQ\t')) { const m = /\tSN:([^\t]+)/.exec(l); if (m && !tid.has(m[1])) tid.set(m[1], tid.size); }
	return tid;
}

/* opts: { headerText, sources: [byte sources of SAM text (chunks)], keepSecondary, tmp: { create(name) -> {write(Uint8Array), close() -> source}, remove(name) },
 *         out: { raw(Uint8Array) }, progress(frac) }  -> stats */
async function sortSam(opts) {
	const tid = headerRefs(opts.headerText), nRef = tid.size;
	const st = { records: 0, kept: 0, supplementary: 0, secondary: 0, unmappedR1: 0 };
	const buckets = new Map(); // tid -> { w, parts, n }
	const bucketOf = async (t) => {
		let b = buckets.get(t);
		if (!b) { b = { w: await opts.tmp.create('bucket.' + (t < 0 ? 'unmapped' : t)), parts: [], n: 0 }; buckets.set(t, b); }
		return b;
	};
	const flushB = (b) => { if (b.parts.length) { b.w.write(enc.encode(b.parts.join(''))); b.parts = []; b.n = 0; } };
	const totalBytes = opts.sources.reduce((a, s) => a + s.size, 0);
	let done = 0;
	for (const src of opts.sources) {
		for (const l of IO.lines(src, { onProgress: (p) => opts.progress && opts.progress(0.5 * (done + p) / Math.max(1, totalBytes)) })) {
			if (!l.length || l.charCodeAt(0) === 64) continue;
			st.records += 1;
			let t1 = l.indexOf('\t'), t2 = l.indexOf('\t', t1 + 1);
			const flag = +l.slice(t1 + 1, t2);
			if (flag & 0x800) { st.supplementary += 1; continue; }
			if ((flag & 0x100) && !opts.keepSecondary) { st.secondary += 1; continue; }
			const f = l.split('\t', 8);
			let line = l;
			if (f[6] === f[2] && f[2] !== '*') { const p = l.split('\t'); p[6] = '='; line = p.join('\t'); }
			const t = f[2] === '*' ? -1 : tid.get(f[2]);
			if (t === undefined) throw new Error(`SAM: reference ${f[2]} not in the header`);
			const b = await bucketOf(t);
			b.parts.push(line, '\n'); b.n += line.length + 1;
			if (b.n > (1 << 20)) flushB(b);
			st.kept += 1;
		}
		done += src.size;
	}
	const order = [...buckets.keys()].filter((t) => t >= 0).sort((a, b) => a - b);
	if (buckets.has(-1)) order.push(-1);
	let wrote = 0;
	const ob = new Uint8Array(4 << 20); let on = 0; // output batching
	const emit = (u) => {
		if (on + u.length > ob.length) { opts.out.raw(ob.slice(0, on)); on = 0; }
		if (u.length > ob.length) { opts.out.raw(u.slice()); return; }
		ob.set(u, on); on += u.length;
	};
	for (const t of order) {
		const b = buckets.get(t); flushB(b);
		const src = b.w.close();
		const bytes = src.read(0, src.size).slice();
		// index the lines
		const starts = [], ends = [], keys = [];
		let s = 0;
		while (s < bytes.length) {
			let e = bytes.indexOf(10, s); if (e < 0) e = bytes.length;
			// POS = 4th field, FLAG = 2nd
			let k = s, tab = 0, flag = 0, pos = 0;
			while (k < e && tab < 4) {
				const c = bytes[k++];
				if (c === 9) { tab++; continue; }
				if (tab === 1) flag = flag * 10 + (c - 48);
				else if (tab === 3) pos = pos * 10 + (c - 48);
			}
			starts.push(s); ends.push(e + 1); keys.push(pos * 2 + ((flag & 16) ? 1 : 0));
			s = e + 1;
		}
		const idx = Array.from(keys.keys());
		if (t >= 0) idx.sort((x, y) => keys[x] - keys[y] || x - y);
		const nl = enc.encode('\n');
		for (const i of idx) { emit(bytes.subarray(starts[i], Math.min(ends[i], bytes.length))); if (ends[i] > bytes.length) emit(nl); }
		wrote += idx.length;
		if (opts.tmp.remove) await opts.tmp.remove('bucket.' + (t < 0 ? 'unmapped' : t));
		if (opts.progress) opts.progress(0.5 + 0.5 * wrote / Math.max(1, st.kept));
	}
	if (on) opts.out.raw(ob.slice(0, on));
	return st;
}

const api = { sortSam, headerRefs };
if (typeof module !== 'undefined' && module.exports) module.exports = api;
root.CuSam = api;
})(typeof self !== 'undefined' ? self : globalThis);
