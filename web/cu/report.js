// SPDX-License-Identifier: LicenseRef-CRISPRuno — HTML report like CRISPRunoCore.make_report (CRISPRuno 0.1.22, (c) MGH), with SVG plots
'use strict';
/* buildReport(io, S, { chromLens }) -> writes <root>.<plot>.svg for every plot and the self-contained report <root>.html
 * Plot objects come from the step outputs (.trimPrimers.info line 3, .final.info line 8, .CRISPResso.info line 1,
 * .summary.plot.json), as CRISPRuno's cached PlotObjects; they are ordered like the original report. */
(function (root) {
const PL = root.CuPlots || require('./plots.js');
const esc = PL.esc;

function table(text) { // 2-column numeric tables written by CRISPRuno (skip '#' lines and the header)
	const xs = [], ys = [];
	let head = true;
	for (const l of String(text || '').split('\n')) {
		if (!l || l[0] === '#') continue;
		if (head) { head = false; continue; }
		const [a, b] = l.split('\t'); xs.push(a); ys.push(parseFloat(b));
	}
	return [xs, ys];
}
function stack(svgs, w = 820) { // vertical concatenation of standalone SVGs
	let y = 0, inner = '';
	for (const s of svgs) {
		const m = /viewBox="0 0 ([\d.]+) ([\d.]+)"/.exec(s); const hh = m ? parseFloat(m[2]) : 300, ww = m ? parseFloat(m[1]) : w;
		inner += s.replace('<svg ', `<svg x="0" y="${y}" `).replace(/ style="[^"]*"/, '').replace(/ width="[\d.]+" height="[\d.]+"/, ` width="${ww}" height="${hh}"`);
		y += hh; w = Math.max(w, ww);
	}
	return `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${w} ${y}" width="${w}" height="${y}" style="max-width:100%;height:auto;background:#fff">${inner}</svg>`;
}
function side(a, b) {
	const dim = (s) => { const m = /viewBox="0 0 ([\d.]+) ([\d.]+)"/.exec(s); return [parseFloat(m[1]), parseFloat(m[2])]; };
	const [w1, h1] = dim(a), [w2, h2] = dim(b);
	const fix = (s, x, w, h) => s.replace('<svg ', `<svg x="${x}" y="0" `).replace(/ style="[^"]*"/, '').replace(/ width="[\d.]+" height="[\d.]+"/, ` width="${w}" height="${h}"`);
	return `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${w1 + w2} ${Math.max(h1, h2)}" width="${w1 + w2}" height="${Math.max(h1, h2)}" style="max-width:100%;height:auto;background:#fff">${fix(a, 0, w1, h1)}${fix(b, w1, w2, h2)}</svg>`;
}

async function render(p, io, ctx) {
	const P = p.plot || {}, rd = async (n) => (await io.readText(n)) || '';
	switch (P.type) {
		case 'pie': return PL.pie(P.items, { title: P.title || p.title });
		case 'donut': return PL.donut(P.inner, P.outer, { title: P.title });
		case 'bar': { const [xs, ys] = table(await rd(P.file)); return PL.bars(P.categorical ? xs : xs.map(Number), ys, { title: P.title, ylabel: P.ylabel, categorical: !!P.categorical }); }
		case 'bar2': {
			const out = [];
			for (let i = 0; i < 2; ++i) { const [xs, ys] = table(await rd(P.files[i])); out.push(PL.bars(xs.map(Number), ys, { title: P.titles[i], ylabel: 'Number of Reads', h: 320 })); }
			return stack(out);
		}
		case 'dedup': return side(PL.pie(P.pie, { title: 'Duplicate read counts', w: 560, h: 380 }), PL.bars(P.bars.map((x) => x[0]), P.bars.map((x) => x[1]), { title: 'Reads per UMI', ylabel: 'Number of UMIs', w: 520, h: 380 }));
		case 'txorder': return PL.txOrder(P.bars);
		case 'txcounts': return PL.txCounts(P.chroms, P.maxCount, P.originChr, P.originCutPos);
		case 'circos': return PL.circos(ctx.chromLens || [], P.sites, P.originChr, P.originCutPos);
		case 'invdel': {
			const [ix, iy] = table(await rd(P.files[0])), [dx, dy] = table(await rd(P.files[1]));
			const f = (xs, ys) => { const a = [], b = []; xs.forEach((x, i) => { if (Math.abs(+x) < 100000) { a.push(+x); b.push(ys[i]); } }); return [a, b]; };
			const [a1, b1] = f(ix, iy), [a2, b2] = f(dx, dy);
			return side(PL.bars(a1, b1, { title: 'Inversion distance from origin', ylabel: 'Number of Reads', w: 520, h: 420 }), PL.bars(a2, b2, { title: 'Deletion distance from origin', ylabel: 'Number of Reads', w: 520, h: 420 }));
		}
		case 'depth': {
			const [x1, y1] = table(await rd(P.files[0])), [x2, y2] = table(await rd(P.files[1]));
			const xs = x1.map(Number), primerLen = xs.filter((x) => x < 0).length;
			const vl = [[0, '#d00', 'origin location']];
			if (P.minPrimerLength < primerLen) vl.unshift([-primerLen + P.minPrimerLength + 0.5, '#000', 'min primer length']);
			const top = PL.bars(xs, y1, { title: 'Read depth of 100bp after origin', ylabel: 'Read depth (Number of reads)', xlabel: 'Position (bp)', vlines: vl });
			// bottom: positions 1..100000, log-spaced subset like CRISPRuno (250 points)
			const pos = x2.map(Number), keep = [], seen = new Set();
			const post = []; pos.forEach((x, i) => { if (x >= 1) post.push([x, y2[i]]); });
			if (post.length) {
				const endPow = Math.log10(post.length), step = endPow / 250;
				for (let k = 0; k < 250; ++k) { const ind = k > 0 ? Math.trunc(10 ** (k * step)) : 1; if (!seen.has(ind) && ind < post.length) { seen.add(ind); keep.push(ind); } }
				keep.sort((a, b) => a - b); keep.push(post.length - 1);
			}
			const bx = keep.map((i) => post[i][0]), by = keep.map((i) => post[i][1]);
			const bottom = PL.bars(bx, by, { title: 'Read depth of 100kb after origin', ylabel: 'Read depth (Number of reads)', xlabel: 'Position (Log scale bp)', logx: true, fill: true });
			return stack([top, bottom]);
		}
		default: return null;
	}
}

const ORDER = ['.trimPrimers.counts', '.final.chr_alignments', '.final.insertSizes', '.final.deduplication_by_UMI_and_aln_pos', '.final.deduplication_by_UMI_and_final_cut_assignment',
	'.final.classifications', '.final.classifications_with_indels', '.final.fragment_translocation_list', '.final.fragment_translocation_list.counts', '.final.fragment_translocation_list.circos',
	'.final.origin_indel_histogram', '.final.origin_inversion_deletion_histogram', '.final.origin_depth_counts', '.final.r1_r2_support', '.final.r1_r2_distances', '.final.discarded_reads',
	'.CRISPResso.classificationSummary', '.summary'];

async function collectPlots(io, S) {
	const out = [];
	const add = (o) => { if (o && o !== 'None') out.push(o); };
	const line = async (name, i) => { const t = await io.readText(name); if (!t) return null; const l = t.split('\n')[i]; return l && l !== 'None' ? l : null; };
	let l;
	if ((l = await line(S.root + '.trimPrimers.info', 2))) add(JSON.parse(l));
	if ((l = await line(S.root + '.final.info', 7))) JSON.parse(l).forEach(add);
	if ((l = await line(S.root + '.CRISPResso.info', 0))) add(JSON.parse(l));
	if ((l = await line(S.root + '.summary.plot.json', 0))) add(JSON.parse(l));
	const rank = (p) => { const k = ORDER.findIndex((s) => p.name === S.root + s); return k < 0 ? 99 : k; };
	return out.sort((a, b) => a.order - b.order || rank(a) - rank(b));
}

async function buildReport(io, S, ctx = {}) {
	const plots = await collectPlots(io, S);
	const cards = [];
	for (const [i, p] of plots.entries()) {
		let svg = null;
		try { svg = await render(p, io, ctx); } catch (e) { io.log(`[report] ${p.name}: ${e.message}`, 'warn'); }
		const svgName = p.name + '.svg';
		if (svg) { const w = await io.writer(svgName); w.write('<?xml version="1.0" encoding="utf-8"?>\n' + svg); w.close(); }
		const base = (n) => String(n).split('/').pop();
		cards.push(`<div class="card" id="plot${i}"><div class="card-header"><h5>${esc(p.title)}</h5></div><div class="card-body">` +
			(svg ? `<div class="fig">${svg}</div><p class="m-0"><small><a href="${esc(base(svgName))}">${esc(base(svgName))}</a></small></p>` : '') +
			`<label>${p.label}</label>` + p.datas.map(([lab, f]) => `<p class="m-0"><small>Data: <a href="${esc(base(f))}">${esc(lab)}</a></small></p>`).join('') + '</div></div>');
	}
	const nav = plots.map((p, i) => `<li><a href="#plot${i}">${esc(p.title)}</a></li>`).join('');
	const html = `<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>Report</title>
<style>
html,body{height:100%}body{margin:0;padding-bottom:40px;background:#f5f5f5;font-family:Lato,"Helvetica Neue",Arial,sans-serif;color:#212529;font-size:15px}
nav{width:200px;position:fixed;top:0;left:0;height:100%;overflow-y:auto;background:#343a40;padding:10px;box-sizing:border-box}
nav .brand{color:#fff;font-size:20px;display:block;margin:6px 0 12px;text-decoration:none}nav ul{list-style:none;margin:0;padding:0}
nav a{color:rgba(255,255,255,.75);text-decoration:none;display:block;padding:6px 0;font-size:14px}nav a:hover{color:#fff}
.container{margin-left:220px;max-width:1000px;padding:0 20px}.head{text-align:center;padding:30px 0 10px}.head h1{font-size:56px;font-weight:300;margin:10px 0}
.card{background:#fff;border:1px solid rgba(0,0,0,.125);border-radius:4px;margin-bottom:10px;text-align:center}.card-header{background:rgba(0,0,0,.03);border-bottom:1px solid rgba(0,0,0,.125);padding:10px}
.card-header h5{margin:0;font-size:18px;font-weight:500}.card-body{padding:16px}.fig>svg{max-width:100%;height:auto}.m-0{margin:0}label{display:block;margin:8px 0;font-size:14px}
.text-monospace{font-family:monospace}.text-break{word-break:break-all}a{color:#18bc9c}
@media (max-width:760px){nav{position:static;width:auto;height:auto}.container{margin-left:0}}
</style></head><body>
<nav><a class="brand" href="#">CRISPRuno</a><ul>${nav}</ul></nav>
<div class="container"><div class="head"><h1>CRISPRuno</h1><hr><h2>Report</h2><p style="color:#777;font-size:13px">CRISPRuno Web (port of CRISPRuno 0.1.22; alignment: minibwa) · ${esc(S.root)}</p></div>
${cards.join('\n')}
</div></body></html>
`;
	const w = await io.writer(S.root + '.html'); w.write(html); w.close();
	return { plots: plots.length };
}

const api = { buildReport, collectPlots, render };
if (typeof module !== 'undefined' && module.exports) module.exports = api;
root.CuReport = api;
})(typeof self !== 'undefined' ? self : globalThis);
