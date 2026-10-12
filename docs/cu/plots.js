// SPDX-License-Identifier: MIT — SVG versions of the CRISPRuno report plots (matplotlib figures in the original)
'use strict';
/* Each function returns a standalone <svg> string. The data comes from the text files CRISPRuno writes
 * (the plots show the same numbers; layout and styling are new). */
(function (root) {
const PAL = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'];
const SET1 = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00']; // Programmed, Off-target, Known, Cas-OFFinder, Novel
const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
const trim0 = (x) => (x.includes('.') ? x.replace(/0+$/, '').replace(/\.$/, '') : x);
const fmt = (v) => (Math.abs(v) >= 1e6 ? trim0((v / 1e6).toPrecision(3)) + 'M' : Math.abs(v) >= 1e4 ? trim0((v / 1e3).toPrecision(3)) + 'k' : String(Math.round(v * 100) / 100));
const svgOpen = (w, h) => `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${w} ${h}" width="${w}" height="${h}" font-family="Helvetica,Arial,sans-serif" font-size="12" style="max-width:100%;height:auto;background:#fff">`;
const title = (w, t) => (t ? `<text x="${w / 2}" y="22" text-anchor="middle" font-size="${Math.max(9, Math.min(15, (w - 20) / (0.56 * String(t).length))).toFixed(1)}" font-weight="600">${esc(t)}</text>` : '');
function niceTicks(lo, hi, n = 5) {
	if (!(hi > lo)) { hi = lo + 1; }
	const span = hi - lo, step0 = span / n, mag = 10 ** Math.floor(Math.log10(step0)), r = step0 / mag;
	const step = (r < 1.5 ? 1 : r < 3 ? 2 : r < 7 ? 5 : 10) * mag;
	const out = []; for (let v = Math.ceil(lo / step) * step; v <= hi + step * 1e-9; v += step) out.push(Math.round(v / step) * step);
	return out;
}
function logTicks(lo, hi) { const out = []; for (let e = Math.floor(Math.log10(lo)); e <= Math.ceil(Math.log10(hi)); ++e) out.push(10 ** e); return out.filter((v) => v >= lo * 0.999 && v <= hi * 1.001); }

/* pie with labels "name\n(n)" and percentages, like ax.pie(..., autopct='%1.2f%%'); values <= 0 dropped */
function pie(items, { title: t = '', w = 900, h = 460, colors = PAL } = {}) {
	const it = items.filter(([, v]) => Number.isInteger(v) && v > 0);
	const tot = it.reduce((a, [, v]) => a + v, 0);
	if (!tot) return svgOpen(w, h) + title(w, t) + `<text x="${w / 2}" y="${h / 2}" text-anchor="middle" fill="#888">(no reads)</text></svg>`;
	// label sides do not depend on the position of the pie: lay out the labels first, then fit the pie between them
	const labels = []; let a = Math.PI / 2; // 12 o'clock, clockwise
	for (const [lab, v] of it) { const da = (2 * Math.PI * v) / tot; labels.push({ lab, v, da, a0: a, mid: a - da / 2, right: Math.cos(a - da / 2) >= 0 }); a -= da; }
	const tw = (l) => 6.2 * (String(l.lab).length + String(l.v).length + 11) + 34; // ~11px Helvetica
	const sideW = (r) => Math.max(0, ...labels.filter((l) => l.right === r).map(tw));
	const lw = sideW(false), rw = sideW(true);
	let R = Math.min(h / 2 - 60, (w - lw - rw - 20) / 2);
	if (R < 110) { R = Math.min(110, h / 2 - 60); w = Math.ceil(lw + rw + 2 * R + 20); }
	const cx = 10 + lw + R, cy = h / 2 + 14;
	let s = svgOpen(w, h) + title(w, t);
	const arc = (b, e) => { const large = e - b > Math.PI ? 1 : 0; const p = (x) => [cx + R * Math.cos(x), cy - R * Math.sin(x)]; const [x0, y0] = p(b), [x1, y1] = p(e); return `M${cx},${cy} L${x0.toFixed(2)},${y0.toFixed(2)} A${R},${R} 0 ${large} 0 ${x1.toFixed(2)},${y1.toFixed(2)} Z`; };
	labels.forEach((l, i) => {
		const col = colors[i % colors.length], pct = (100 * l.v / tot);
		s += labels.length === 1 ? `<circle cx="${cx}" cy="${cy}" r="${R}" fill="${col}" stroke="#fff"/>` : `<path d="${arc(l.a0 - l.da, l.a0)}" fill="${col}" stroke="#fff" stroke-width="1.5"><title>${esc(l.lab)}: ${l.v}</title></path>`;
		l.y = cy - (R + 16) * Math.sin(l.mid);
		const f = l.da * R; // arc length available for the percentage
		if (f > 60) { const px = cx + R * 0.62 * Math.cos(l.mid), py = cy - R * 0.62 * Math.sin(l.mid); s += `<text x="${px.toFixed(1)}" y="${py.toFixed(1)}" text-anchor="middle" font-size="11" fill="#fff" font-weight="600">${pct.toFixed(2)}%</text>`; }
		else if (f > 22) { const px = cx + R * 0.8 * Math.cos(l.mid), py = cy - R * 0.8 * Math.sin(l.mid); s += `<text x="${px.toFixed(1)}" y="${py.toFixed(1)}" text-anchor="middle" font-size="9" fill="#fff">${pct.toFixed(1)}%</text>`; }
	});
	let bottom = 0;
	for (const side of [true, false]) { // keep labels on each side at least 20 px apart
		const L = labels.filter((l) => l.right === side).sort((p, q) => p.y - q.y);
		for (let k = 1; k < L.length; ++k) if (L[k].y - L[k - 1].y < 20) L[k].y = L[k - 1].y + 20;
		for (const l of L) {
			const ex = cx + (R + 4) * Math.cos(l.mid), ey = cy - (R + 4) * Math.sin(l.mid);
			const lx = cx + (side ? 1 : -1) * (R + 26);
			s += `<polyline points="${ex.toFixed(1)},${ey.toFixed(1)} ${(lx - (side ? 6 : -6)).toFixed(1)},${(l.y - 4).toFixed(1)} ${lx.toFixed(1)},${(l.y - 4).toFixed(1)}" fill="none" stroke="#999" stroke-width="0.7"/>`;
			s += `<text x="${(lx + (side ? 3 : -3)).toFixed(1)}" y="${l.y.toFixed(1)}" text-anchor="${side ? 'start' : 'end'}" font-size="11">${esc(l.lab)} <tspan fill="#555">(${l.v}, ${(100 * l.v / tot).toFixed(2)}%)</tspan></text>`;
			bottom = Math.max(bottom, l.y + 8);
		}
	}
	if (bottom > h) s = s.replace(/viewBox="0 0 ([\d.]+) [\d.]+" width="([\d.]+)" height="[\d.]+"/, `viewBox="0 0 $1 ${Math.ceil(bottom)}" width="$2" height="${Math.ceil(bottom)}"`);
	return s + '</svg>';
}
/* classification donut: inner pie by class, outer ring by indel category (silver/pink/crimson) */
function donut(inner, outer, { title: t = 'Read classifications', w = 680, h = 480 } = {}) {
	let s = svgOpen(w, h) + title(w, t);
	const cx = w / 2 - 40, cy = h / 2 + 10, R = Math.min(w - 80, h) / 2 - 70, Ri = R * 0.9;
	const COL = { silver: '#c0c0c0', pink: '#ffc0cb', crimson: '#dc143c' };
	const ring = (items, r0, r1, colorOf, labels) => {
		const tot = items.reduce((a, x) => a + x.v, 0); if (!tot) return '';
		let o = '', a = Math.PI / 2;
		items.forEach((x, i) => {
			const da = 2 * Math.PI * x.v / tot, b = a - da, mid = a - da / 2;
			const P = (r, ang) => `${(cx + r * Math.cos(ang)).toFixed(2)},${(cy - r * Math.sin(ang)).toFixed(2)}`;
			const large = da > Math.PI ? 1 : 0;
			o += items.length === 1 ? `<circle cx="${cx}" cy="${cy}" r="${(r0 + r1) / 2}" fill="none" stroke="${colorOf(x, i)}" stroke-width="${r1 - r0}"/>`
				: `<path d="M${P(r1, a)} A${r1},${r1} 0 ${large} 1 ${P(r1, b)} L${P(r0, b)} A${r0},${r0} 0 ${large} 0 ${P(r0, a)} Z" fill="${colorOf(x, i)}" stroke="#fff" stroke-width="2"><title>${esc(x.label || '')} ${x.v}</title></path>`;
			if (labels && da > 0.12) { const px = cx + (r0 + r1) / 2 * Math.cos(mid), py = cy - (r0 + r1) / 2 * Math.sin(mid); o += `<text x="${px.toFixed(1)}" y="${py.toFixed(1)}" text-anchor="middle" font-size="11">${esc(x.label)}<tspan x="${px.toFixed(1)}" dy="13">(${x.v}) ${(100 * x.v / tot).toFixed(2)}%</tspan></text>`; }
			if (!labels && da > 0.08) { const px = cx + (r1 + 14) * Math.cos(mid), py = cy - (r1 + 14) * Math.sin(mid); o += `<text x="${px.toFixed(1)}" y="${py.toFixed(1)}" text-anchor="middle" font-size="10" fill="#444">${(100 * x.v / tot).toFixed(2)}%</text>`; }
			a = b;
		});
		return o;
	};
	s += ring(inner.map(([l, v]) => ({ label: l, v })), 0, Ri, (x, i) => PAL[i % PAL.length], true);
	s += ring(outer.map(([v, c]) => ({ v })), Ri + 2, R + 2, (x, i) => COL[outer[i][1]] || outer[i][1], false);
	const L = [['No indels', 'silver'], ['Short indels <=50bp', 'pink'], ['Long indels >50bp', 'crimson']];
	L.forEach(([l, c], i) => { s += `<rect x="${w - 170}" y="${50 + i * 20}" width="14" height="14" fill="${COL[c]}"/><text x="${w - 150}" y="${62 + i * 20}">${esc(l)}</text>`; });
	return s + '</svg>';
}
/* bar chart. xs numeric (or categorical labels), ys counts */
function bars(xs, ys, { title: t = '', ylabel = '', xlabel = '', w = 820, h = 360, categorical = false, logx = false, fill = false, vlines = [], ylog = false } = {}) {
	const ml = 70, mr = 20, mt = 36, mb = categorical ? 90 : 46, pw = w - ml - mr, ph = h - mt - mb;
	let s = svgOpen(w, h) + title(w, t);
	if (!xs.length) return s + `<text x="${w / 2}" y="${h / 2}" text-anchor="middle" fill="#888">(no data)</text></svg>`;
	let ymax = Math.max(...ys, 0) * 1.05 || 1, ymin = 0;
	const yt = niceTicks(0, ymax);
	ymax = Math.max(ymax, yt[yt.length - 1]);
	const Y = (v) => mt + ph - (ph * (v - ymin)) / (ymax - ymin);
	let X, bw;
	if (categorical) { bw = pw / xs.length; X = (i) => ml + bw * (i + 0.5); }
	else if (logx) { const lo = Math.log10(Math.max(1, Math.min(...xs))), hi = Math.log10(Math.max(...xs) + 1); X = (v) => ml + pw * (Math.log10(Math.max(v, 1)) - lo) / Math.max(1e-9, hi - lo); bw = 1; }
	else { let lo = Math.min(...xs), hi = Math.max(...xs); if (hi === lo) { lo -= 1; hi += 1; } const pad = (hi - lo) * 0.02 + 0.5; lo -= pad; hi += pad; X = (v) => ml + pw * (v - lo) / (hi - lo); bw = Math.max(1, pw / (hi - lo) * 0.8); X.lo = lo; X.hi = hi; }
	for (const v of yt) s += `<line x1="${ml}" x2="${ml + pw}" y1="${Y(v).toFixed(1)}" y2="${Y(v).toFixed(1)}" stroke="#eee"/><text x="${ml - 6}" y="${(Y(v) + 4).toFixed(1)}" text-anchor="end" font-size="11">${fmt(v)}</text>`;
	if (fill) {
		let d = `M${X(xs[0]).toFixed(2)},${Y(0).toFixed(2)}`;
		xs.forEach((x, i) => { d += ` L${X(x).toFixed(2)},${Y(ys[i]).toFixed(2)}`; });
		d += ` L${X(xs[xs.length - 1]).toFixed(2)},${Y(0).toFixed(2)} Z`;
		s += `<path d="${d}" fill="${PAL[0]}" fill-opacity="0.85" stroke="${PAL[0]}"/>`;
	} else xs.forEach((x, i) => {
		const cxp = categorical ? X(i) : X(x), y = Y(ys[i]);
		s += `<rect x="${(cxp - bw * (categorical ? 0.4 : 0.5)).toFixed(2)}" y="${y.toFixed(2)}" width="${(categorical ? bw * 0.8 : bw).toFixed(2)}" height="${Math.max(0, Y(0) - y).toFixed(2)}" fill="${PAL[0]}"><title>${esc(x)}: ${ys[i]}</title></rect>`;
	});
	for (const [xv, col, lab] of vlines) { const xx = X(xv); s += `<line x1="${xx}" x2="${xx}" y1="${mt}" y2="${mt + ph}" stroke="${col}" stroke-dasharray="3,3"/>` + (lab ? `<text x="${xx + 4}" y="${mt + 12}" fill="${col}" font-size="10">${esc(lab)}</text>` : ''); }
	s += `<line x1="${ml}" x2="${ml + pw}" y1="${mt + ph}" y2="${mt + ph}" stroke="#333"/><line x1="${ml}" x2="${ml}" y1="${mt}" y2="${mt + ph}" stroke="#333"/>`;
	if (categorical) xs.forEach((x, i) => { s += `<text transform="translate(${X(i).toFixed(1)},${mt + ph + 10}) rotate(45)" font-size="11">${esc(x)}</text>`; });
	else if (logx) { const lo = Math.max(1, Math.min(...xs)), hi = Math.max(...xs); for (const v of logTicks(lo, hi)) s += `<line x1="${X(v)}" x2="${X(v)}" y1="${mt + ph}" y2="${mt + ph + 4}" stroke="#333"/><text x="${X(v)}" y="${mt + ph + 16}" text-anchor="middle" font-size="11">${fmt(v)}</text>`; }
	else for (const v of niceTicks(X.lo, X.hi, 8)) s += `<line x1="${X(v)}" x2="${X(v)}" y1="${mt + ph}" y2="${mt + ph + 4}" stroke="#333"/><text x="${X(v).toFixed(1)}" y="${mt + ph + 16}" text-anchor="middle" font-size="11">${fmt(v)}</text>`;
	if (ylabel) s += `<text transform="translate(16,${mt + ph / 2}) rotate(-90)" text-anchor="middle">${esc(ylabel)}</text>`;
	if (xlabel) s += `<text x="${ml + pw / 2}" y="${h - 6}" text-anchor="middle">${esc(xlabel)}</text>`;
	return s + '</svg>';
}
/* translocation summary: left label | colored box | right label, and log-scale count bars */
function txOrder(items) {
	const n = items.length, rowH = 26, w = 980, h = 70 + n * rowH + 40;
	let s = svgOpen(w, h) + title(w, 'Translocation summary (top fragments)');
	const cats = ['Programmed', 'Off-target', 'Known', 'Cas-OFFinder', 'Novel'];
	const positions = [...new Set(items.filter((x) => !x.other).map((x) => x.cut))].sort();
	const vir = (t) => { const c = [[68, 1, 84], [59, 82, 139], [33, 145, 140], [94, 201, 98], [253, 231, 37]]; const f = t * (c.length - 1), i = Math.min(c.length - 2, Math.floor(f)), r = f - i; return `rgb(${c[i].map((v, k) => Math.round(v + (c[i + 1][k] - v) * r)).join(',')})`; };
	const maxC = Math.max(5, ...items.map((x) => x.count));
	const bx0 = 620, bx1 = w - 60, X = (v) => bx0 + (bx1 - bx0) * Math.log10(Math.max(1, v)) / Math.log10(maxC);
	items.forEach((x, i) => {
		const y = 50 + i * rowH;
		const fillc = x.other ? '#cccccc' : vir(positions.length > 1 ? positions.indexOf(x.cut) / positions.length : 0);
		const outl = x.other ? 'none' : SET1[Math.max(0, cats.indexOf(x.anno))] || SET1[4];
		s += `<text x="236" y="${y + 16}" text-anchor="end">${esc(x.left)}</text><rect x="244" y="${y + 2}" width="90" height="${rowH - 6}" fill="${fillc}" stroke="${outl}" stroke-width="2.5"/>`;
		s += `<text x="342" y="${y + 16}" class="mono">${esc(x.right)}</text>`;
		s += `<rect x="${bx0}" y="${y + 7}" width="${Math.max(1, X(Math.max(1, x.count)) - bx0).toFixed(1)}" height="${rowH - 14}" fill="${PAL[0]}"/><text x="${(X(Math.max(1, x.count)) + 4).toFixed(1)}" y="${y + 16}" font-size="11">${x.count}</text>`;
	});
	const yb = 50 + n * rowH + 6;
	s += `<line x1="${bx0}" x2="${bx1}" y1="${yb}" y2="${yb}" stroke="#333"/>`;
	for (const v of logTicks(1, maxC)) s += `<line x1="${X(v)}" x2="${X(v)}" y1="${yb}" y2="${yb + 4}" stroke="#333"/><text x="${X(v)}" y="${yb + 16}" text-anchor="middle" font-size="11">${fmt(v)}</text>`;
	s += `<text x="${(bx0 + bx1) / 2}" y="${yb + 32}" text-anchor="middle">Number of reads</text>`;
	cats.forEach((c, i) => { s += `<rect x="${20 + i * 110}" y="${yb + 20}" width="14" height="14" fill="none" stroke="${SET1[i]}" stroke-width="2.5"/><text x="${38 + i * 110}" y="${yb + 31}" font-size="11">${esc(c)}</text>`; });
	return s + '</svg>';
}
/* reads at left/right of each top cut site, one panel per chromosome (log y) */
function txCounts(chroms, maxCount, originChr, originCutPos) {
	const panelH = 170, w = 980, h = 40 + chroms.length * panelH;
	let s = svgOpen(w, h) + title(w, 'Number of reads aligned at left and right sides of cut sites');
	const maxCuts = Math.max(1, ...chroms.map(([, l]) => l.length));
	chroms.forEach(([chr, locs], ci) => {
		const top = 40 + ci * panelH, ml = 80, pw = w - ml - 30, ph = panelH - 70, slot = pw / maxCuts;
		const Y = (v) => top + ph - ph * Math.log10(Math.max(1, v)) / Math.log10(Math.max(2, maxCount * 2));
		s += `<text transform="translate(18,${top + ph / 2}) rotate(-90)" text-anchor="middle" font-size="12">${esc(chr)}</text><line x1="${ml}" x2="${ml}" y1="${top}" y2="${top + ph}" stroke="#333"/><line x1="${ml}" x2="${ml + pw}" y1="${top + ph}" y2="${top + ph}" stroke="#333"/>`;
		for (const v of logTicks(1, Math.max(2, maxCount * 2))) s += `<text x="${ml - 5}" y="${Y(v) + 4}" text-anchor="end" font-size="10">${fmt(v)}</text><line x1="${ml}" x2="${ml + pw}" y1="${Y(v)}" y2="${Y(v)}" stroke="#f0f0f0"/>`;
		locs.forEach(([pos, lr], i) => {
			const x0 = ml + slot * i + slot * 0.1, bw = Math.min(30, slot * 0.38);
			for (const [k, col, dx] of [['left', PAL[0], 0], ['right', PAL[1], bw]]) { const v = lr[k]; if (v > 0) s += `<rect x="${(x0 + dx).toFixed(1)}" y="${Y(v).toFixed(1)}" width="${bw.toFixed(1)}" height="${(top + ph - Y(v)).toFixed(1)}" fill="${col}"><title>${chr}:${pos} ${k}: ${v}</title></rect>`; }
			s += `<text transform="translate(${(x0 + bw).toFixed(1)},${top + ph + 8}) rotate(60)" font-size="9">${pos}${chr === originChr && pos === originCutPos ? ' (origin)' : ''}</text>`;
		});
	});
	s += `<rect x="${w - 150}" y="8" width="12" height="12" fill="${PAL[0]}"/><text x="${w - 134}" y="18">left</text><rect x="${w - 90}" y="8" width="12" height="12" fill="${PAL[1]}"/><text x="${w - 74}" y="18">right</text>`;
	return s + '</svg>';
}
/* circos-style overview: chromosomes on a circle, chords from the origin to each site (top 100), origin chromosome enlarged */
function circos(chromLens, sites, originChr, originCutPos) {
	const w = 720, h = 720, cx = w / 2, cy = h / 2 + 10, R = 280;
	let s = svgOpen(w, h) + title(w, 'Translocations from the origin');
	const tx = new Set(sites.map(([k]) => k.split(':')[0]));
	let chrs = chromLens.filter(([c]) => !c.includes('_') || tx.has(c));
	let oc = originChr, opos = originCutPos;
	const lens = new Map(chrs);
	if (oc === null) { const mn = Math.min(...chrs.map(([, l]) => l)); chrs = [...chrs, ['primer', Math.floor(mn / 2)]]; lens.set('primer', Math.floor(mn / 2)); oc = 'primer'; opos = 40; }
	const maxL = Math.max(...chrs.map(([, l]) => l)), mult = Math.max(1, Math.floor(maxL * 4 / (lens.get(oc) || 1)));
	const eff = chrs.map(([c, l]) => [c, c === oc ? l * mult : l]);
	const tot = eff.reduce((a, [, l]) => a + l, 0), gap = 0.012, avail = 2 * Math.PI - gap * eff.length;
	const start = new Map(); let a = -Math.PI / 2;
	eff.forEach(([c, l], i) => { start.set(c, [a, l]); const da = avail * l / tot; const a1 = a + da;
		const P = (r, ang) => `${(cx + r * Math.cos(ang)).toFixed(2)},${(cy + r * Math.sin(ang)).toFixed(2)}`;
		s += `<path d="M${P(R, a)} A${R},${R} 0 ${da > Math.PI ? 1 : 0} 1 ${P(R, a1)} L${P(R + 16, a1)} A${R + 16},${R + 16} 0 ${da > Math.PI ? 1 : 0} 0 ${P(R + 16, a)} Z" fill="${PAL[i % PAL.length]}" opacity="0.85"><title>${esc(c)}</title></path>`;
		if (da > 0.06) { const m = (a + a1) / 2; s += `<text x="${(cx + (R + 30) * Math.cos(m)).toFixed(1)}" y="${(cy + (R + 30) * Math.sin(m)).toFixed(1)}" text-anchor="middle" font-size="10">${esc(c.replace(/^chr/, ''))}</text>`; }
		a = a1 + gap; });
	const ang = (c, p) => { const st = start.get(c); if (!st) return null; const [a0, l] = st; return a0 + avail * Math.min(p * (c === oc ? mult : 1), l) / tot; };
	const maxC = Math.max(1, ...sites.map(([, v]) => v));
	const ao = ang(oc, opos);
	for (const [k, v] of [...sites].reverse()) {
		const [c, p] = k.split(':'); const at = ang(c, parseInt(p, 10)); if (at === null || ao === null) continue;
		const r = R - 4, t = Math.max(0.1, v / maxC), g = Math.round(220 - 200 * t);
		const x0 = cx + r * Math.cos(ao), y0 = cy + r * Math.sin(ao), x1 = cx + r * Math.cos(at), y1 = cy + r * Math.sin(at);
		s += `<path d="M${x0.toFixed(1)},${y0.toFixed(1)} Q${cx},${cy} ${x1.toFixed(1)},${y1.toFixed(1)}" fill="none" stroke="rgb(${g},${g},${g})" stroke-width="${(0.6 + 3 * t).toFixed(2)}" opacity="0.9"><title>${esc(k)}: ${v}</title></path>`;
		const rr = R - 18 - 40 * Math.log(v) / Math.log(maxC + 1);
		s += `<circle cx="${(cx + rr * Math.cos(at)).toFixed(1)}" cy="${(cy + rr * Math.sin(at)).toFixed(1)}" r="2.2" fill="#000"/>`;
	}
	return s + '</svg>';
}

const api = { pie, donut, bars, txOrder, txCounts, circos, esc };
if (typeof module !== 'undefined' && module.exports) module.exports = api;
root.CuPlots = api;
})(typeof self !== 'undefined' ? self : globalThis);
