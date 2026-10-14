// SPDX-License-Identifier: MIT — SAM post-processing for the CRISPRuno web port
'use strict';
/* Equivalent of the end of CRISPRuno's align_reads():  aligner | samtools view [-F 256] -Shu - | samtools sort
 *   - records with flag 0x800 (supplementary) are dropped: bowtie2 never writes them, minibwa does
 *   - secondary records (0x100) are dropped unless keepSecondary (analyze_multimap_assignments, bowtie2 -k N)
 *   - RNEXT equal to RNAME is printed as '=' (as samtools view does)
 *   - samtools sort order: reference index (@SQ order, unmapped last), POS, reverse strand; ties keep input order
 * Large inputs: records are first appended to one temporary file per reference ("bucket"). A bucket up to runBytes
 * (256 MB) is sorted in memory; a larger one (e.g. amplicon data where almost every read is on one chromosome) is
 * sorted in runs of runBytes that are then merged, so memory stays ~runBytes whatever the input size.
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
	const out = batcher((u) => opts.out.raw(u));
	const nl = enc.encode('\n');
	const RUN = opts.runBytes || (256 << 20); // a bucket larger than this is sorted in runs of this size and merged
	for (const t of order) {
		const b = buckets.get(t); flushB(b);
		const src = b.w.close(), bname = 'bucket.' + (t < 0 ? 'unmapped' : t);
		if (src.size <= RUN) {
			wrote += sortSegment(src.read(0, src.size), t >= 0, out.emit, nl);
		} else {
			// external merge sort: sorted runs (consecutive pieces of the bucket, so ties keep input order) ...
			const runs = [];
			for (let pos = 0, r = 0; pos < src.size; ++r) {
				let seg = src.read(pos, RUN);
				if (pos + seg.length < src.size) { const e = seg.lastIndexOf(10); if (e < 0) throw new Error('SAM line longer than the sort buffer'); seg = seg.subarray(0, e + 1); }
				pos += seg.length;
				const name = `${bname}.run${r}`, w = await opts.tmp.create(name), rb = batcher((u) => w.write(u));
				sortSegment(seg, t >= 0, rb.emit, nl); rb.flush();
				runs.push({ name, src: w.close() });
				seg = null;
			}
			if (opts.tmp.remove) await opts.tmp.remove(bname);
			// ... merged by (key, run index)
			const rd = runs.map((x) => lineReader(x.src));
			const cur = rd.map((r) => r.next()), key = cur.map((l) => (l ? keyOf(l, 0, l.length) : 0));
			for (;;) {
				let m = -1;
				for (let i = 0; i < cur.length; ++i) if (cur[i] && (m < 0 || (t >= 0 && key[i] < key[m]))) m = i;
				if (m < 0) break;
				out.emit(cur[m]); ++wrote;
				cur[m] = rd[m].next(); if (cur[m]) key[m] = keyOf(cur[m], 0, cur[m].length);
				if ((wrote & 0xffff) === 0 && opts.progress) opts.progress(0.5 + 0.5 * wrote / Math.max(1, st.kept));
			}
			for (const x of runs) if (opts.tmp.remove) await opts.tmp.remove(x.name);
		}
		if (src.size <= RUN && opts.tmp.remove) await opts.tmp.remove(bname);
		if (opts.progress) opts.progress(0.5 + 0.5 * wrote / Math.max(1, st.kept));
	}
	out.flush();
	return st;
}

/* sort key of a SAM line bytes[s, e): POS * 2 + reverse strand (samtools sort; the reference is fixed per bucket) */
function keyOf(bytes, s, e) {
	let k = s, tab = 0, flag = 0, pos = 0;
	while (k < e && tab < 4) {
		const c = bytes[k++];
		if (c === 9) { tab++; continue; }
		if (tab === 1) flag = flag * 10 + (c - 48);
		else if (tab === 3) pos = pos * 10 + (c - 48);
	}
	return pos * 2 + ((flag & 16) ? 1 : 0);
}
/* emit the lines of bytes (newline-terminated SAM text) in stable key order (input order if !sorted); returns the line count */
function sortSegment(bytes, sorted, emit, nl) {
	const starts = [], ends = [], keys = [];
	let s = 0;
	while (s < bytes.length) {
		let e = bytes.indexOf(10, s); if (e < 0) e = bytes.length;
		starts.push(s); ends.push(e + 1); keys.push(keyOf(bytes, s, e));
		s = e + 1;
	}
	const idx = Array.from(keys.keys());
	if (sorted) idx.sort((x, y) => keys[x] - keys[y] || x - y);
	for (const i of idx) { emit(bytes.subarray(starts[i], Math.min(ends[i], bytes.length))); if (ends[i] > bytes.length) emit(nl); }
	return idx.length;
}
/* 4 MB write batching */
function batcher(write) {
	const ob = new Uint8Array(4 << 20); let on = 0;
	const flush = () => { if (on) { write(ob.slice(0, on)); on = 0; } };
	const emit = (u) => {
		if (on + u.length > ob.length) flush();
		if (u.length > ob.length) { write(u.slice()); return; }
		ob.set(u, on); on += u.length;
	};
	return { emit, flush };
}
/* sequential reader of newline-terminated lines (each line includes its '\n'; valid until the next call) */
function lineReader(src, block = 4 << 20) {
	let pos = 0, buf = new Uint8Array(0), off = 0;
	return { next() {
		let e = buf.indexOf(10, off);
		while (e < 0) {
			if (pos >= src.size) return null; // runs always end with '\n'
			const rest = buf.subarray(off), nb = src.read(pos, block); pos += nb.length;
			const m = new Uint8Array(rest.length + nb.length); m.set(rest); m.set(nb, rest.length);
			buf = m; off = 0; e = buf.indexOf(10, off);
		}
		const line = buf.subarray(off, e + 1); off = e + 1; return line;
	} };
}


const api = { sortSam, headerRefs };
if (typeof module !== 'undefined' && module.exports) module.exports = api;
root.CuSam = api;
})(typeof self !== 'undefined' ? self : globalThis);
