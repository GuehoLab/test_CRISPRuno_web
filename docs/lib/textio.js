// SPDX-License-Identifier: AGPL-3.0 — part of the GUIDE-seq port (derived from tsailabSJ/guideseq V2, AGPL-3.0)
'use strict';
/* Byte sources, line reading, FASTA random access and buffered text output for the GUIDE-seq steps.
 * A "source" is { size, read(pos, len) -> Uint8Array } (synchronous: File via FileReaderSync in a worker,
 * an OPFS SyncAccessHandle, or a Node.js file in the tests).
 */
(function (root) {
const dec = new TextDecoder();
const enc = new TextEncoder();

function blobSource(blob) {
	const fr = new FileReaderSync();
	return { size: blob.size, read: (pos, len) => new Uint8Array(fr.readAsArrayBuffer(blob.slice(pos, Math.min(blob.size, pos + len)))) };
}
function handleSource(h) { // FileSystemSyncAccessHandle
	const size = h.getSize();
	return { size, read: (pos, len) => { const b = new Uint8Array(Math.max(0, Math.min(len, size - pos))); const n = h.read(b, { at: pos }); return n === b.length ? b : b.subarray(0, n); } };
}
// block cache for many small random reads (FASTA windows, .l2b intervals): FileReaderSync costs ~1 ms per call
function cachedSource(src, block = 1 << 20, maxBlocks = 64) {
	const cache = new Map();
	const get = (k) => {
		let b = cache.get(k);
		if (b) { cache.delete(k); cache.set(k, b); return b; }
		b = src.read(k * block, block).slice();
		cache.set(k, b);
		if (cache.size > maxBlocks) cache.delete(cache.keys().next().value);
		return b;
	};
	return { size: src.size, read: (pos, len) => {
		len = Math.max(0, Math.min(len, src.size - pos));
		const k0 = Math.floor(pos / block), k1 = Math.floor((pos + len - 1) / block);
		if (len === 0) return new Uint8Array(0);
		if (k0 === k1) { const b = get(k0), o = pos - k0 * block; return b.subarray(o, o + len); }
		const out = new Uint8Array(len); let w = 0;
		for (let k = k0; k <= k1; ++k) { const b = get(k), o = k === k0 ? pos - k0 * block : 0, n = Math.min(b.length - o, len - w); out.set(b.subarray(o, o + n), w); w += n; }
		return out;
	} };
}
function bytesSource(u8) { return { size: u8.length, read: (pos, len) => u8.subarray(pos, Math.min(u8.length, pos + len)) }; }

// iterate lines (without the trailing '\n'; a final line without '\n' is returned too). onProgress(bytesDone)
function* lines(src, opts = {}) {
	const CH = opts.chunk || 8 << 20;
	let pos = 0, carry = '';
	while (pos < src.size) {
		const b = src.read(pos, CH); if (!b.length) break;
		pos += b.length;
		let last = b.lastIndexOf(10);
		let text;
		if (last < 0) { carry += dec.decode(b, { stream: true }); continue; }
		text = carry + dec.decode(b.subarray(0, last + 1));
		carry = last + 1 < b.length ? dec.decode(b.subarray(last + 1), { stream: true }) : '';
		let s = 0;
		for (;;) { const e = text.indexOf('\n', s); if (e < 0) break; yield text.slice(s, e); s = e + 1; }
		if (opts.onProgress) opts.onProgress(pos);
	}
	carry += dec.decode();
	if (carry.length) yield carry;
}
// like lines() but also gives byte offsets: yields [line, startOffset, byteLengthIncludingNewline]
function* linesWithOffsets(src, opts = {}) {
	const CH = opts.chunk || 8 << 20;
	let pos = 0, base = 0, rest = new Uint8Array(0);
	while (pos < src.size || rest.length) {
		let b = pos < src.size ? src.read(pos, CH) : new Uint8Array(0);
		pos += b.length;
		if (rest.length) { const t = new Uint8Array(rest.length + b.length); t.set(rest); t.set(b, rest.length); b = t; }
		const eof = pos >= src.size;
		let s = 0;
		for (;;) {
			let e = b.indexOf(10, s);
			if (e < 0) { if (eof && s < b.length) e = b.length; else break; }
			yield [dec.decode(b.subarray(s, e)), base + s, Math.min(e + 1, b.length) - s];
			s = e + 1;
			if (s >= b.length) break;
		}
		base += Math.min(s, b.length);
		rest = s < b.length ? b.slice(s) : new Uint8Array(0);
		if (eof) break;
		if (opts.onProgress) opts.onProgress(pos);
	}
}

// buffered text writer over sink(Uint8Array)
class TextSink {
	constructor(sink, cap = 4 << 20) { this.sink = sink; this.parts = []; this.n = 0; this.cap = cap; this.bytes = 0; }
	write(s) { this.parts.push(s); this.n += s.length; if (this.n >= this.cap) this.flush(); }
	flush() { if (!this.parts.length) return; const b = enc.encode(this.parts.join('')); this.parts = []; this.n = 0; this.bytes += b.length; this.sink(b); }
}

/* FASTA with a samtools-style .fai (pyfaidx compatible slicing). */
class Fasta {
	constructor(src, faiText) {
		this.src = src; this.idx = new Map();
		if (faiText) this.parseFai(faiText); else this.buildFai();
	}
	parseFai(t) {
		for (const l of t.split('\n')) {
			if (!l.trim()) continue;
			const f = l.split('\t');
			this.idx.set(f[0], { len: +f[1], off: +f[2], lb: +f[3], lB: +f[4] });
		}
	}
	buildFai(onProgress) { // one sequential pass over the bytes; same rules as samtools faidx (name = header up to the first white space)
		const CH = 16 << 20, size = this.src.size;
		let cur = null, pos = 0, lineStart = 0, partial = null;
		const fin = () => { if (cur) this.idx.set(cur.name, cur); };
		const onLine = (startOff, bytes, nl) => { // bytes: Uint8Array of the line without '\n'; nl: 1 if terminated
			if (bytes.length && bytes[0] === 62) { // '>'
				fin();
				const name = new TextDecoder().decode(bytes).slice(1).split(/[ \t\r]/)[0];
				cur = { name, len: 0, off: startOff + bytes.length + nl, lb: 0, lB: 0, lastShort: false };
				return;
			}
			if (!cur) return;
			const L = bytes.length && bytes[bytes.length - 1] === 13 ? bytes.length - 1 : bytes.length;
			if (cur.lastShort && L) throw new Error(`FASTA ${cur.name}: different line length in sequence (re-wrap the FASTA, e.g. samtools faidx would fail as well)`);
			if (cur.lb === 0) { cur.lb = L; cur.lB = bytes.length + nl; }
			else if (L !== cur.lb) cur.lastShort = true;
			if (L > cur.lb) throw new Error(`FASTA ${cur.name}: inconsistent line length`);
			cur.len += L;
		};
		while (pos < size) {
			const b = this.src.read(pos, CH); if (!b.length) break;
			let s = 0;
			for (;;) {
				const e = b.indexOf(10, s);
				if (e < 0) break;
				if (partial) { const t = new Uint8Array(partial.length + e - s); t.set(partial); t.set(b.subarray(s, e), partial.length); onLine(lineStart, t, 1); partial = null; }
				else onLine(pos + s, b.subarray(s, e), 1);
				s = e + 1; lineStart = pos + s;
			}
			if (s < b.length) { const rest = b.subarray(s); if (partial) { const t = new Uint8Array(partial.length + rest.length); t.set(partial); t.set(rest, partial.length); partial = t; } else { partial = rest.slice(); lineStart = pos + s; } }
			pos += b.length;
			if (onProgress) onProgress(pos);
		}
		if (partial && partial.length) onLine(lineStart, partial, 0);
		fin();
	}
	has(name) { return this.idx.has(name); }
	length(name) { const r = this.idx.get(name); return r ? r.len : 0; }
	faiText() { return [...this.idx].map(([n, r]) => `${n}\t${r.len}\t${r.off}\t${r.lb}\t${r.lB}\n`).join(''); }
	// genome[name][start:end] with Python slice semantics (negative indices count from the end, clamped)
	slice(name, start, end) {
		const r = this.idx.get(name);
		if (!r) throw new Error(`KeyError: '${name}' not in the reference FASTA`);
		const L = r.len;
		const norm = (v) => (v < 0 ? Math.max(0, v + L) : Math.min(v, L));
		const a = norm(start), b = norm(end);
		if (a >= b) return '';
		const o1 = r.off + Math.floor(a / r.lb) * r.lB + (a % r.lb);
		const o2 = r.off + Math.floor((b - 1) / r.lb) * r.lB + ((b - 1) % r.lb) + 1;
		const bytes = this.src.read(o1, o2 - o1);
		let s = '';
		for (let i = 0; i < bytes.length; ++i) { const c = bytes[i]; if (c !== 10 && c !== 13) s += String.fromCharCode(c); }
		return s;
	}
}

const api = { blobSource, handleSource, bytesSource, cachedSource, lines, linesWithOffsets, TextSink, Fasta, enc, dec };
if (typeof module !== 'undefined' && module.exports) module.exports = api;
root.GsIO = api;
})(typeof self !== 'undefined' ? self : globalThis);

/* Python-style readline() over a byte source, optionally gzip (multi-member, via fflate), with universal newlines
 * ('\r\n' and '\r' become '\n', like open(..., 'r') / gzip.open(..., 'rt')). */
(function (root) {
const IO = root.GsIO;
class LineReader {
	constructor(src, opts = {}) {
		this.src = src; this.pos = 0; this.CH = opts.chunk || 4 << 20;
		this.gz = opts.gz === undefined ? null : opts.gz;
		this.lines = []; this.li = 0; this.carry = ''; this.eof = false; this.out = [];
		this.dec = new TextDecoder();
		if (this.gz === null) { const h = src.read(0, 2); this.gz = h.length === 2 && h[0] === 0x1f && h[1] === 0x8b; }
		if (this.gz) {
			const ff = root.fflate || (typeof require !== 'undefined' ? require('./fflate.min.js') : null);
			this.inflater = new ff.Gunzip((chunk) => this.out.push(chunk));
		}
	}
	get bytesRead() { return this.pos; }
	fill() { // decode more text into this.lines; false at EOF
		while (!this.eof) {
			let data;
			if (this.pos >= this.src.size) {
				this.eof = true;
				if (this.gz) { this.inflater.push(new Uint8Array(0), true); data = this.drain(); } else data = new Uint8Array(0);
			} else {
				const b = this.src.read(this.pos, this.CH); this.pos += b.length;
				if (this.gz) { this.inflater.push(b, false); data = this.drain(); } else data = b;
			}
			let text = this.carry + this.dec.decode(data, { stream: !this.eof });
			this.carry = '';
			if (!this.eof && text.endsWith('\r')) { this.carry = '\r'; text = text.slice(0, -1); }
			if (text.indexOf('\r') >= 0) text = text.replace(/\r\n?/g, '\n');
			const parts = text.split('\n');
			const last = parts.pop();
			if (this.li < this.lines.length) { /* not reached: fill() only when empty */ }
			this.lines = parts.map((p) => p + '\n'); this.li = 0;
			if (this.eof) { if (last.length) this.lines.push(last); }
			else if (last.length) this.carry = last + this.carry; // partial line: in front of the next chunk
			if (this.lines.length) return true;
		}
		return false;
	}
	drain() {
		if (this.out.length === 1) { const o = this.out[0]; this.out = []; return o; }
		let n = 0; for (const c of this.out) n += c.length;
		const o = new Uint8Array(n); let k = 0; for (const c of this.out) { o.set(c, k); k += c.length; }
		this.out = []; return o;
	}
	readline() {
		if (this.li >= this.lines.length && !this.fill()) return '';
		return this.lines[this.li++];
	}
	has(n) { return true; } // synchronous reader: readline() fills on demand
	async ensure() {}
}
IO.LineReader = LineReader;

/* Asynchronous variant for a Blob/File (user files): gzip through the browser's native DecompressionStream
 * (~4x faster than fflate); if that fails (DecompressionStream stops at the end of the first member of a
 * multi-member gzip file), it continues with fflate, skipping the bytes already delivered.
 * Usage: `if (!r.has(4)) await r.ensure(4);` then up to 4 synchronous readline() calls. */
class BlobLineReader {
	constructor(blob, opts = {}) {
		this.blob = blob; this.CH = opts.chunk || 4 << 20; this.pos = 0;
		this.lines = []; this.li = 0; this.carry = ''; this.eof = false; this.dec = new TextDecoder();
		this.mode = null; this.out = 0; this.skip = 0; this.pending = [];
	}
	get bytesRead() { return this.pos; }
	has(n) { return this.lines.length - this.li >= n || this.eof; }
	async ensure(n) { while (!this.has(n)) await this.fill(); }
	readline() { return this.li < this.lines.length ? this.lines[this.li++] : ''; }
	async next() { // next decompressed chunk or null at EOF
		if (this.mode === null) {
			const h = new Uint8Array(await this.blob.slice(0, 2).arrayBuffer());
			const gz = h.length === 2 && h[0] === 0x1f && h[1] === 0x8b;
			if (!gz) this.mode = 'plain';
			else if (typeof DecompressionStream !== 'undefined') { this.mode = 'native'; this.rd = this.blob.stream().pipeThrough(new DecompressionStream('gzip')).getReader(); this.nIn = 0; }
			else this.startFflate(0);
		}
		if (this.mode === 'plain') {
			if (this.pos >= this.blob.size) return null;
			const b = new Uint8Array(await this.blob.slice(this.pos, this.pos + this.CH).arrayBuffer()); this.pos += b.length; return b;
		}
		if (this.mode === 'native') {
			try {
				const { done, value } = await this.rd.read();
				if (done) return null;
				this.out += value.length; this.pos = Math.min(this.blob.size, Math.round(this.out / 4)); // progress estimate
				return value;
			} catch (e) { this.startFflate(this.out); }
		}
		for (;;) { // fflate
			if (this.pending.length) return this.pending.shift();
			if (this.fpos >= this.blob.size) { if (this.fdone) return null; this.fdone = true; this.inflater.push(new Uint8Array(0), true); continue; }
			const b = new Uint8Array(await this.blob.slice(this.fpos, this.fpos + this.CH).arrayBuffer()); this.fpos += b.length; this.pos = this.fpos;
			this.inflater.push(b, false);
		}
	}
	startFflate(skip) {
		const ff = (typeof self !== 'undefined' && self.fflate) || (typeof require !== 'undefined' ? require('./fflate.min.js') : null);
		this.mode = 'fflate'; this.fpos = 0; this.fdone = false; this.skip = skip;
		this.inflater = new ff.Gunzip((c) => {
			if (this.skip >= c.length) { this.skip -= c.length; return; }
			if (this.skip) { c = c.subarray(this.skip); this.skip = 0; }
			this.pending.push(c);
		});
	}
	async fill() {
		const data = await this.next();
		if (data === null) { this.eof = true; const t = this.carry + this.dec.decode(); this.carry = ''; const tt = t.replace(/\r\n?/g, '\n'); this.push(tt, true); return; }
		let text = this.carry + this.dec.decode(data, { stream: true });
		this.carry = '';
		if (text.endsWith('\r')) { this.carry = '\r'; text = text.slice(0, -1); }
		if (text.indexOf('\r') >= 0) text = text.replace(/\r\n?/g, '\n');
		this.push(text, false);
	}
	push(text, final) {
		const parts = text.split('\n');
		const last = parts.pop();
		const rest = this.lines.slice(this.li);
		this.lines = rest.concat(parts.map((p) => p + '\n')); this.li = 0;
		if (final) { if (last.length) this.lines.push(last); }
		else if (last.length) this.carry = last + this.carry;
	}
}
IO.BlobLineReader = BlobLineReader;
})(typeof self !== 'undefined' ? self : globalThis);

/* Reference sequence from a minibwa .l2b (2-bit packed + soft-mask and ambiguity intervals), same slice() as Fasta.
 * Differences from the FASTA: every ambiguous base reads as N (IUPAC codes other than N are lost), and a soft-masked
 * or N run that reaches the very end of a contig is not recorded in .l2b (minibwa l2b_add_seq), so such bases come
 * back upper case / as the random base minibwa stored. */
(function (root) {
const IO = root.GsIO;
class L2bGenome {
	constructor(src) {
		this.src = src;
		const h = new DataView(src.read(0, 64).slice().buffer);
		if (String.fromCharCode(h.getUint8(0), h.getUint8(1), h.getUint8(2)) !== 'L2B') throw new Error('not a .l2b file');
		const u = (o) => Number(h.getBigUint64(o, true));
		const nCtg = u(8), nAmbi = u(24), nMask = u(32), lenName = u(40), nPac = u(56);
		let o = 64;
		const lens = new DataView(src.read(o, nCtg * 8).slice().buffer); o += nCtg * 8;
		this.ambiOff = o; this.nAmbi = nAmbi; o += nAmbi * 16;
		this.maskOff = o; this.nMask = nMask; o += nMask * 16;
		this.pacOff = o; o += nPac * 8;
		const names = new TextDecoder().decode(src.read(o, lenName)).split('\0');
		this.idx = new Map(); let off = 0;
		for (let i = 0; i < nCtg; ++i) { const L = Number(lens.getBigUint64(i * 8, true)); this.idx.set(names[i], { len: L, off }); off += L; }
	}
	has(n) { return this.idx.has(n); }
	length(n) { const r = this.idx.get(n); return r ? r.len : 0; }
	contigs() { return [...this.idx].map(([n, r]) => [n, r.len]); }
	ivAt(base, i) { const d = new DataView(this.src.read(base + i * 16, 16).slice().buffer); return [Number(d.getBigUint64(0, true)), Number(d.getBigUint64(8, true))]; }
	overlapping(base, n, g0, g1) { // intervals (sorted, disjoint) overlapping [g0, g1)
		let lo = 0, hi = n;
		while (lo < hi) { const m = (lo + hi) >> 1; if (this.ivAt(base, m)[1] > g0) hi = m; else lo = m + 1; }
		const out = [];
		for (let i = lo; i < n; ++i) { const iv = this.ivAt(base, i); if (iv[0] >= g1) break; out.push(iv); }
		return out;
	}
	slice(name, start, end) {
		const r = this.idx.get(name);
		if (!r) throw new Error(`KeyError: '${name}' not in the reference index`);
		const L = r.len, norm = (v) => (v < 0 ? Math.max(0, v + L) : Math.min(v, L));
		const a = norm(start), b = norm(end);
		if (a >= b) return '';
		const g0 = r.off + a, g1 = r.off + b;
		const b0 = Math.floor(g0 / 4), bytes = this.src.read(this.pacOff + b0, Math.floor((g1 - 1) / 4) - b0 + 1);
		const s = new Array(g1 - g0);
		// offsets relative to byte b0 (genome coordinates can exceed 2^31, where >> and & would wrap)
		for (let d = g0 - b0 * 4, i = 0; i < g1 - g0; ++i, ++d) s[i] = 'ACGT'[(bytes[d >> 2] >> ((d & 3) * 2)) & 3];
		for (const [x, y] of this.overlapping(this.ambiOff, this.nAmbi, g0, g1)) for (let g = Math.max(x, g0); g < Math.min(y, g1); ++g) s[g - g0] = 'N';
		for (const [x, y] of this.overlapping(this.maskOff, this.nMask, g0, g1)) for (let g = Math.max(x, g0); g < Math.min(y, g1); ++g) s[g - g0] = s[g - g0].toLowerCase();
		return s.join('');
	}
}
IO.L2bGenome = L2bGenome;
})(typeof self !== 'undefined' ? self : globalThis);
