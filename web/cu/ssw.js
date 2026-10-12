// SPDX-License-Identifier: MIT — wrapper around ssw.wasm (Striped Smith-Waterman, Mengyao Zhao et al., MIT; ssw.c as shipped in CRISPRuno lib/ssw)
'use strict';
/* Same calls as CRISPRuno's SSWPrimerAlign (ctypes on libssw.so):
 *   5-letter alphabet A C G T N (lower case too; other characters -> N), match 2, mismatch -5, N scores 0,
 *   ssw_init(primer, score_size 2), ssw_align(profile, read, gap_open 5, gap_extend 5, flag 1, 0, 0, maskLen 0).
 * In SSW terms the primer/origin is the "read" (query profile) and the sequencing read is the "reference":
 *   nRefBeg/nRefEnd = positions in the sequencing read, nQryBeg/nQryEnd = positions in the primer/origin.
 * Usage: const S = await CuSsw.load(bytes); const a = S.primer(origin); a.align(read) -> {nScore, nRefBeg, nRefEnd, nQryBeg, nQryEnd}
 */
(function (root) {
const CODE = new Uint8Array(256).fill(4);
'ACGT'.split('').forEach((c, i) => { CODE[c.charCodeAt(0)] = i; CODE[c.toLowerCase().charCodeAt(0)] = i; });
CODE[78] = CODE[110] = 4;

function instantiate(module) {
	let mem = null;
	const wasi = {
		fd_write: (fd, iovs, n, nw) => { if (mem) new DataView(mem.buffer).setUint32(nw, 0, true); return 0; }, // stderr messages of ssw.c are dropped
		fd_close: () => 0, fd_seek: () => 70, random_get: (p, n) => { const b = new Uint8Array(mem.buffer, p, n); for (let i = 0; i < n; ++i) b[i] = (Math.random() * 256) | 0; return 0; },
	};
	const inst = new WebAssembly.Instance(module, { wasi_snapshot_preview1: wasi });
	mem = inst.exports.memory;
	if (inst.exports._initialize) inst.exports._initialize();
	return inst;
}

class SswLib {
	constructor(module) { this.inst = instantiate(module); this.x = this.inst.exports; this.refPtr = 0; this.refCap = 0; }
	u8() { return new Uint8Array(this.x.memory.buffer); }
	primer(seq, opts = {}) { return new SswPrimerAlign(this, seq, opts); }
	ensureRef(n) {
		if (n <= this.refCap) return this.refPtr;
		if (this.refPtr) this.x.free(this.refPtr);
		this.refCap = Math.max(1024, n * 2); this.refPtr = this.x.malloc(this.refCap);
		return this.refPtr;
	}
}
class SswPrimerAlign {
	constructor(lib, primerSeq, { match = 2, mismatch = 5, gapOpen = 5, gapExtend = 5, nFlag = 1 } = {}) {
		this.lib = lib; this.primer_seq = primerSeq; this.gapOpen = gapOpen; this.gapExtend = gapExtend; this.nFlag = nFlag;
		const x = lib.x, n = 5;
		this.matPtr = x.malloc(n * n);
		let u8 = lib.u8();
		for (let i = 0; i < n * n; ++i) u8[this.matPtr + i] = 0;
		for (let i = 0; i < n - 1; ++i) for (let j = 0; j < n - 1; ++j) u8[this.matPtr + i * n + j] = (i === j ? match : -mismatch) & 0xff;
		this.qPtr = x.malloc(Math.max(1, primerSeq.length));
		u8 = lib.u8();
		for (let i = 0; i < primerSeq.length; ++i) u8[this.qPtr + i] = CODE[primerSeq.charCodeAt(i) & 0xff];
		this.prof = x.ssw_init(this.qPtr, primerSeq.length, this.matPtr, n, 2);
	}
	align(read) {
		const lib = this.lib, x = lib.x, L = read.length;
		const p = lib.ensureRef(Math.max(1, L));
		let u8 = lib.u8();
		for (let i = 0; i < L; ++i) u8[p + i] = CODE[read.charCodeAt(i) & 0xff];
		const r = x.ssw_align(this.prof, p, L, this.gapOpen, this.gapExtend, this.nFlag, 0, 0, 0);
		if (!r) throw new Error('ssw_align failed');
		const dv = new DataView(x.memory.buffer);
		const res = { nScore: dv.getUint16(r, true), nScore2: dv.getUint16(r + 2, true), nRefBeg: dv.getInt32(r + 4, true), nRefEnd: dv.getInt32(r + 8, true),
			nQryBeg: dv.getInt32(r + 12, true), nQryEnd: dv.getInt32(r + 16, true) };
		x.align_destroy(r);
		return res;
	}
	destroy() { const x = this.lib.x; x.init_destroy(this.prof); x.free(this.qPtr); x.free(this.matPtr); }
}

async function load(src) { // src: ArrayBuffer/Uint8Array, a URL, or a WebAssembly.Module
	let mod = src;
	if (!(src instanceof WebAssembly.Module)) {
		let bytes = src;
		if (typeof src === 'string') bytes = await (await fetch(src)).arrayBuffer();
		mod = await WebAssembly.compile(bytes);
	}
	return new SswLib(mod);
}
function loadSync(bytes) { return new SswLib(new WebAssembly.Module(bytes)); }

/* CRISPRunoCore.trim_left_primer_from_read / trim_right_primer_from_read */
function trimLeft(read, a, minScore = 40, allowIndels = false) {
	const r = a.align(read);
	if (!allowIndels) {
		if (r.nRefBeg !== 0) return ['', read, 0];
		if (r.nRefEnd - r.nRefBeg + 1 !== r.nQryEnd - r.nQryBeg + 1) return ['', read, 0];
	}
	if (r.nScore < minScore) return ['', read, 0];
	return [read.slice(0, r.nRefEnd + 1), read.slice(r.nRefEnd + 1), r.nQryEnd];
}
function trimRight(read, a, minScore = 40, allowIndels = false) {
	const r = a.align(read);
	if (!allowIndels) {
		if (r.nRefEnd - r.nRefBeg + 1 !== r.nQryEnd - r.nQryBeg + 1) return [read, ''];
	}
	if (r.nScore < minScore) return [read, ''];
	return [read.slice(0, r.nRefBeg), read.slice(r.nRefBeg)];
}

const api = { load, loadSync, trimLeft, trimRight };
if (typeof module !== 'undefined' && module.exports) module.exports = api;
root.CuSsw = api;
})(typeof self !== 'undefined' ? self : globalThis);
