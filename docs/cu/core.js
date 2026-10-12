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
class DefaultInt extends Map {
	get(k) { if (!super.has(k)) super.set(k, 0); return super.get(k); }
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
const api = { revcomp, revcompShared, strip, rstripNl, findWrongNt, str, F, PyF, pctStr, DefaultInt, pyJson, PlotObject, read4, faidx, base, pyFloat };
if (typeof module !== 'undefined' && module.exports) module.exports = api;
root.CuCore = api;
})(typeof self !== 'undefined' ? self : globalThis);
