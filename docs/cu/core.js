// SPDX-License-Identifier: LicenseRef-CRISPRuno — helpers for the CRISPRuno web port (CRISPRuno (c) MGH; academic non-commercial license, see README)
'use strict';
/* Python-compatible helpers shared by the CRISPRuno step modules. */
(function (root) {
const PY = root.GsPy || require('../lib/pyfmt.js');
const { pyFloat } = PY;

const NT = { A: 'T', C: 'G', G: 'C', T: 'A', N: 'N', a: 't', c: 'g', g: 'c', t: 'a', n: 'n', _: '_', '-': '-' };
function revcomp(seq) { // CRISPRunoCore.reverse_complement (KeyError for other characters)
	let o = '';
	for (let i = seq.length - 1; i >= 0; --i) { const c = NT[seq[i]]; if (c === undefined) throw new Error(`KeyError: '${seq[i]}' (reverse_complement of ${seq})`); o += c; }
	return o;
}
const NT2 = { A: 'T', C: 'G', G: 'C', T: 'A', N: 'N', _: '_', '-': '-' };
function revcompShared(seq) { // CRISPRessoShared.reverse_complement (upper case only)
	let o = '';
	for (let i = seq.length - 1; i >= 0; --i) { const c = NT2[seq[i]]; if (c === undefined) throw new Error(`KeyError: '${seq[i]}'`); o += c; }
	return o;
}
// Python str.strip() (ASCII whitespace + the unicode spaces JS also knows)
const strip = (s) => s.replace(/^\s+|\s+$/g, '');
const rstripNl = (s) => (s.endsWith('\n') ? s.slice(0, -1) : s);
function findWrongNt(seq) { const s = new Set(seq.toUpperCase()); for (const c of 'ATCGN') s.delete(c); return [...s]; }

/* str() of the values CRISPRuno writes: ints as JS numbers, floats wrapped with F(), None, bools */
const F = PY.F, PyF = PY.PyF;
function str(v) {
	if (v === null || v === undefined) return 'None';
	if (v instanceof PyF) return pyFloat(v.v);
	if (typeof v === 'boolean') return v ? 'True' : 'False';
	return String(v);
}
// Python round(x, 2) as printed by '%s' (float)
function pctStr(x) { return pyFloat((root.CuAlign || require('./cralign.js')).pyRound(x, 2)); }

/* defaultdict(int) with Python's insertion-on-read semantics */
/* A substring (split/slice) of a long string can keep the whole parent alive in V8 (sliced strings), e.g. a read name
 * from IO.lines() retains the 8 MB text chunk it was cut from. Strings kept for the whole step (Map keys/values) are
 * copied with own(); intern() also shares one copy among equal values (chromosome names, cut keys). */
const own = (s) => (typeof s === 'string' && s.length > 12 ? (' ' + s).slice(1) : s);
function makeIntern() { const m = new Map(); return (s) => { let v = m.get(s); if (v === undefined) { v = own(s); m.set(v, v); } return v; }; }
/* Compact string-keyed hash map for millions of keys (read names, UMI/position keys): the key bytes live in 64 MB
 * arena blocks and the entries in typed arrays (~keyLength + 30 bytes per entry, vs ~300+ for a Map of strings).
 * Each entry has nVals int32 values and optionally one auxiliary string. Iteration is in insertion order. */
class StrMap {
	constructor(nVals = 1, withAux = false) {
		this.nVals = nVals; this.withAux = withAux; this.size = 0; this.cap = 0;
		this.blocks = []; this.bpos = 0; this.BLK = 64 << 20;
		this.koff = null; this.klen = null; this.vals = null; this.aoff = null; this.alen = null;
		this.mask = 1023; this.table = new Int32Array(1024).fill(-1);
		this.buf = new Uint8Array(1 << 16); this.enc = new TextEncoder(); this.dec = new TextDecoder();
		this._grow(1024);
	}
	_grow(n) {
		const g = (A, old, k = 1) => { const a = new A(n * k); if (old) a.set(old); return a; };
		this.koff = g(Float64Array, this.koff); this.klen = g(Int32Array, this.klen); this.vals = g(Int32Array, this.vals, this.nVals);
		if (this.withAux) { this.aoff = g(Float64Array, this.aoff); this.alen = g(Int32Array, this.alen); }
		this.cap = n;
	}
	_bytes(s) { // UTF-8 bytes of s in this.buf (returns length)
		const n = s.length;
		if (n * 3 > this.buf.length) this.buf = new Uint8Array(n * 3 + 16);
		const b = this.buf;
		for (let i = 0; i < n; ++i) { const c = s.charCodeAt(i); if (c > 127) return this.enc.encodeInto(s, b).written; b[i] = c; }
		return n;
	}
	_hash(n) { let h = 0x811c9dc5; const b = this.buf; for (let i = 0; i < n; ++i) h = Math.imul(h ^ b[i], 16777619); return h >>> 0; }
	_store(n) { // copy this.buf[0, n) into the arena -> global offset
		if (!this.blocks.length || this.bpos + n > this.BLK) { this.blocks.push(new Uint8Array(Math.max(this.BLK, n))); this.bpos = 0; }
		const bi = this.blocks.length - 1; this.blocks[bi].set(this.buf.subarray(0, n), this.bpos);
		const off = bi * this.BLK + this.bpos; this.bpos += n; return off;
	}
	_eq(i, n) { const off = this.koff[i]; if (this.klen[i] !== n) return false; const blk = this.blocks[Math.floor(off / this.BLK)], o = off % this.BLK, b = this.buf; for (let k = 0; k < n; ++k) if (blk[o + k] !== b[k]) return false; return true; }
	_str(off, n) { const blk = this.blocks[Math.floor(off / this.BLK)], o = off % this.BLK; return this.dec.decode(blk.subarray(o, o + n)); }
	find(s) { const n = this._bytes(s); let h = this._hash(n) & this.mask; for (;;) { const i = this.table[h]; if (i < 0) return -1; if (this._eq(i, n)) return i; h = (h + 1) & this.mask; } }
	has(s) { return this.find(s) >= 0; }
	insert(s, aux) { // s must not be present; returns the entry index
		if ((this.size + 1) * 2 > this.mask + 1) this._rehash((this.mask + 1) * 2);
		if (this.size === this.cap) this._grow(this.cap * 2);
		const i = this.size++;
		const n = this._bytes(s); let h = this._hash(n) & this.mask;
		while (this.table[h] >= 0) h = (h + 1) & this.mask;
		this.table[h] = i; this.koff[i] = this._store(n); this.klen[i] = n;
		if (this.withAux) { const m = this._bytes(aux === undefined ? '' : String(aux)); this.aoff[i] = this._store(m); this.alen[i] = m; }
		return i;
	}
	_rehash(sz) {
		const t = new Int32Array(sz).fill(-1), mask = sz - 1;
		for (let i = 0; i < this.size; ++i) {
			const off = this.koff[i], n = this.klen[i], blk = this.blocks[Math.floor(off / this.BLK)], o = off % this.BLK;
			let h = 0x811c9dc5; for (let k = 0; k < n; ++k) h = Math.imul(h ^ blk[o + k], 16777619);
			h = (h >>> 0) & mask; while (t[h] >= 0) h = (h + 1) & mask; t[h] = i;
		}
		this.table = t; this.mask = mask;
	}
	v(i, j = 0) { return this.vals[i * this.nVals + j]; }
	setV(i, j, x) { this.vals[i * this.nVals + j] = x; }
	key(i) { return this._str(this.koff[i], this.klen[i]); }
	aux(i) { return this._str(this.aoff[i], this.alen[i]); }
	*values(j = 0) { for (let i = 0; i < this.size; ++i) yield this.vals[i * this.nVals + j]; }
}
class DefaultInt extends Map {
	get(k) { if (!super.has(k)) { super.set(own(k), 0); } return super.get(k); }
	peek(k) { return super.has(k) ? super.get(k) : 0; } // .get() without inserting (dict.get(k, 0))
	inc(k, n = 1) { super.set(k, this.get(k) + n); }
}

/* json.dumps(obj, separators=(',', ':')) with ensure_ascii */
function pyJson(v) {
	if (v === null || v === undefined) return 'null';
	if (typeof v === 'boolean') return v ? 'true' : 'false';
	if (typeof v === 'number') return Number.isInteger(v) ? String(v) : pyFloat(v);
	if (v instanceof PyF) return pyFloat(v.v);
	if (typeof v === 'string') return '"' + v.replace(/[\\"\u0000-\u001f\u007f-￿]/g, (c) => {
		const m = { '\\': '\\\\', '"': '\\"', '\n': '\\n', '\r': '\\r', '\t': '\\t', '\b': '\\b', '\f': '\\f' }[c];
		return m || '\\u' + c.charCodeAt(0).toString(16).padStart(4, '0');
	}) + '"';
	if (Array.isArray(v)) return '[' + v.map(pyJson).join(',') + ']';
	if (v instanceof Map) return '{' + [...v].map(([k, x]) => pyJson(String(k)) + ':' + pyJson(x)).join(',') + '}';
	return '{' + Object.keys(v).map((k) => pyJson(k) + ':' + pyJson(v[k])).join(',') + '}';
}

/* PlotObject: in this port the plots are SVG files (<name>.svg) drawn by cu/plots.js */
class PlotObject {
	constructor(name, title, label, datas, order = 50, plot = null) { this.name = name; this.title = title; this.label = label; this.datas = datas; this.order = order; this.plot = plot; }
	toJson() { return pyJson({ plot_name: this.name, plot_title: this.title, plot_label: this.label, plot_datas: this.datas, plot_order: this.order }); }
	static fromObj(o) { return new PlotObject(o.name, o.title, o.label, o.datas, o.order, o.plot); }
}

/* FASTQ record reading with Python's readline() semantics (lines keep their '\n'; '' at EOF).
 * reader: { readline(), has(n), ensure(n) } (GsIO.LineReader or BlobLineReader) */
async function read4(r) { if (!r.has(4)) await r.ensure(4); return [r.readline(), r.readline(), r.readline(), r.readline()]; }

/* samtools faidx -n 10000 genome chr:start-end | tail -n 1 (1-based inclusive; for regions <= 10000 bp) */
function faidx(genome, chr, start, end) {
	if (!genome.has(chr)) throw new Error(`samtools faidx: could not retrieve ${chr}:${start}-${end}`);
	const a = Math.max(1, start), b = Math.min(genome.length(chr), end);
	if (b < a) return '';
	const s = genome.slice(chr, a - 1, b);
	if (s.length <= 10000) return s;
	const last = s.length % 10000 || 10000; // output wrapped at 10000 bp; tail -n 1 keeps the last line only
	return s.slice(s.length - last);
}

const base = (p) => String(p).split(/[\\/]/).pop();
const api = { own, makeIntern, StrMap, revcomp, revcompShared, strip, rstripNl, findWrongNt, str, F, PyF, pctStr, DefaultInt, pyJson, PlotObject, read4, faidx, base, pyFloat };
if (typeof module !== 'undefined' && module.exports) module.exports = api;
root.CuCore = api;
})(typeof self !== 'undefined' ? self : globalThis);
