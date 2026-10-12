// SPDX-License-Identifier: AGPL-3.0 — part of the GUIDE-seq port (derived from tsailabSJ/guideseq V2, AGPL-3.0)
'use strict';
/* Python / numpy number formatting and arithmetic used by guideseq outputs:
 *   pyFloat(x)   == repr(float) / str(float) in Python 3 (shortest round-trip, '.0', e-notation rules)
 *   pyStr(v)     == str(v) for int / float / str / None / bool values as guideseq passes them to print()
 *   npStd(ints)  == numpy.std(list_of_ints) (float64, pairwise summation like numpy's add.reduce)
 */
(function (root) {
function pyFloat(x) {
	if (Number.isNaN(x)) return 'nan';
	if (x === Infinity) return 'inf';
	if (x === -Infinity) return '-inf';
	if (x === 0) return Object.is(x, -0) ? '-0.0' : '0.0';
	const neg = x < 0, s = Math.abs(x).toExponential(); // shortest round-trip digits
	const k = s.indexOf('e'), mant = s.slice(0, k).replace('.', ''), exp = parseInt(s.slice(k + 1), 10);
	let out;
	if (exp >= -4 && exp < 16) {
		if (exp >= 0) {
			const ip = mant.length > exp + 1 ? mant.slice(0, exp + 1) : mant + '0'.repeat(exp + 1 - mant.length);
			const fp = mant.length > exp + 1 ? mant.slice(exp + 1) : '0';
			out = ip + '.' + fp;
		} else out = '0.' + '0'.repeat(-exp - 1) + mant;
	} else {
		const e = (exp < 0 ? '-' : '+') + String(Math.abs(exp)).padStart(2, '0');
		out = (mant.length > 1 ? mant[0] + '.' + mant.slice(1) : mant) + 'e' + e;
	}
	return neg ? '-' + out : out;
}
// value tagging: guideseq mixes ints and floats; JS numbers carry no type, so floats are wrapped
class PyF { constructor(v) { this.v = v; } toString() { return pyFloat(this.v); } }
const F = (v) => new PyF(v);
function pyStr(v) {
	if (v === null || v === undefined) return 'None';
	if (v instanceof PyF) return pyFloat(v.v);
	if (typeof v === 'boolean') return v ? 'True' : 'False';
	return String(v);
}

// numpy DOUBLE_pairwise_sum
function pairwise(a, off, n) {
	if (n < 8) {
		let r = -0.0; // numpy starts from -0.0 (since 1.24; 0. before) -- same result for non-empty sums of finite values except sign of zero
		for (let i = 0; i < n; ++i) r += a[off + i];
		return r;
	} else if (n <= 128) {
		const r = [a[off], a[off + 1], a[off + 2], a[off + 3], a[off + 4], a[off + 5], a[off + 6], a[off + 7]];
		let i;
		for (i = 8; i < n - (n % 8); i += 8) for (let j = 0; j < 8; ++j) r[j] += a[off + i + j];
		let res = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]));
		for (; i < n; ++i) res += a[off + i];
		return res;
	}
	let n2 = Math.floor(n / 2);
	n2 -= n2 % 8;
	return pairwise(a, off, n2) + pairwise(a, off + n2, n - n2);
}
function npSum(a) { // np.add.reduce of a contiguous float64 array (one pairwise pass; verified against numpy for n <= 20000)
	return 0 + pairwise(a, 0, a.length);
}
function npStd(list) {
	const n = list.length, a = Float64Array.from(list);
	const mean = npSum(a) / n;
	const x = new Float64Array(n);
	for (let i = 0; i < n; ++i) { const d = a[i] - mean; x[i] = d * d; }
	return Math.sqrt(npSum(x) / n);
}

const api = { pyFloat, pyStr, PyF, F, npStd, npSum };
if (typeof module !== 'undefined' && module.exports) module.exports = api;
root.GsPy = api;
})(typeof self !== 'undefined' ? self : globalThis);
