// SPDX-License-Identifier: LicenseRef-CRISPRuno — CRISPRuno web port: page logic
'use strict';
/* Orchestrates the steps (one worker per step, run one after another).
 * OPFS: minibwa/crispruno/<run>/ (all outputs, named <root>.* like CRISPRuno), minibwa/cu-jobs/<run>.<job>/ (alignment chunks, temporary)
 */
const STEPS = [
	['prep', 'Prep', 'prep_input — primer/guide location (bowtie2 -c → minibwa), origin sequence, read count and length'],
	['umi', 'UMI', 'add_umi_from_umi_file · analyze_UMIs_initial'],
	['primer', 'Primer/adapter trimming', 'filter_on_primer (SSW · cutadapt -a)'],
	['align', 'Align (minibwa)', 'align_reads (bowtie2 → minibwa · samtools sort)'],
	['assign', 'Read assignment', 'make_final_read_assignments · collapse_cut_points'],
	['crispresso', 'CRISPResso quantification', 'prep_crispresso2 · CRISPResso2 (quantification) · run_and_aggregate_crispresso'],
	['summary', 'Summary · report', 'make_final_summary · make_report (SVG figures)'],
];
let OPT = null; // option form (cu/options.js)
let busy = false, curWorker = null, curReject = null, curIsAlign = false, cancelled = false;

/******** settings & inputs ********/
function fileMap() { const m = {}; for (const f of $('fastqs').files) m[f.name] = f; return m; }
function runName() { return sanitizeName($('run-name').value.trim()); }
function parsed() {
	const run = runName();
	return CuSettings.parse($('settings').value, { root: run || undefined, files: fileMap(), cores: CORES, lenient: true });
}
function renderSettings() {
	const el = $('settings-info');
	try { localStorage.setItem('cu-settings', $('settings').value); } catch (_) {}
	try {
		const { settings: S, warnings } = parsed();
		const rows = [['root', S.root], ['fastq_r1', S.fastq_r1], ['fastq_r2', S.fastq_r2 || '—'], ['fastq_umi', S.fastq_umi || '—'], ['primer_seq', S.primer_seq],
			['guide_sequences', S.guide_sequences.join(', ') || '—'], ['min_primer_length', S.min_primer_length], ['primer_in_r2', S.primer_in_r2]];
		const changed = CuSettings.DEFAULTS.filter(([k, d]) => !['root', 'fastq_r1', 'fastq_r2', 'fastq_umi', 'primer_seq', 'guide_sequences', 'min_primer_length', 'primer_in_r2', 'genome', 'aligner_genome'].includes(k) && JSON.stringify(S[k]) !== JSON.stringify(k === 'n_processes' ? 1 : d));
		el.innerHTML = `<table><tbody>${rows.map(([k, v]) => `<tr><th class="mono">${esc(k)}</th><td class="mono">${esc(v)}</td></tr>`).join('')}</tbody></table>` +
			(changed.length ? `<p class="note">Settings that differ from the defaults: ${changed.map(([k]) => `<span class="mono">${esc(k)}=${esc(JSON.stringify(S[k]))}</span>`).join(', ')}</p>` : '') +
			warnings.map((w) => `<p class="warn">${esc(w)}</p>`).join('');
		if (!$('run-name').value.trim()) $('run-name').placeholder = S.root;
		if (OPT) OPT.sync(S, Object.keys(fileMap()));
	} catch (e) { el.innerHTML = `<span class="err">Settings error: ${esc(e.message)}</span>`; }
	renderInputs();
}
function refFiles() {
	const fs = [...$('fasta').files];
	return { fasta: fs.find((f) => !f.name.endsWith('.fai')) || null, fai: fs.find((f) => f.name.endsWith('.fai')) || null };
}
function renderInputs() {
	const F = fileMap(), R = refFiles();
	$('input-info').innerHTML = `<b>Selected files:</b> ${Object.keys(F).map(esc).join(', ') || 'none'}<br><b>Reference sequence:</b> ${R.fasta ? esc(R.fasta.name) + (R.fai ? ' + .fai' : ' (.fai is built when running)') : 'read from the index .l2b'}`;
}
async function selectedIndex() {
	if (document.querySelector('input[name=src]:checked').value === 'file') {
		const fs = [...$('idx-files').files], l2b = fs.find((f) => f.name.endsWith('.l2b')), mbw = fs.find((f) => f.name.endsWith('.mbw'));
		if (!l2b || !mbw) throw new Error('Select both index files, .l2b and .mbw.');
		return { name: sanitizeName(l2b.name.replace(/\.l2b$/, '')), source: 'file', files: { l2b, mbw }, info: await readIndexHeaders(mbw, l2b) };
	}
	const v = $('opfs-idx').value, x = (await opfsIndexes()).find((i) => i.name === v);
	if (!x) throw new Error('There is no index in browser storage. Build one on the Build index page or select the index files.');
	return { name: x.name, source: 'opfs', info: x };
}
async function refreshIndexes() {
	const list = await opfsIndexes(), sel = $('opfs-idx'), prev = sel.value;
	sel.innerHTML = list.map((i) => `<option value="${esc(i.name)}">${esc(i.name)} (${fmtBytes(i.genome)} bp, SA 1/${2 ** i.saBit})</option>`).join('') || '<option value="">(none)</option>';
	if (prev) sel.value = prev;
}

/******** OPFS ********/
async function cuTop() { return (await opfsTop()).getDirectoryHandle('crispruno', { create: true }); }
async function runDirH(run) { return (await cuTop()).getDirectoryHandle(run, { create: true }); }
async function jobsTop() { return (await opfsTop()).getDirectoryHandle('cu-jobs', { create: true }); }
async function fileIn(dir, name) { try { return await (await dir.getFileHandle(name)).getFile(); } catch (_) { return null; } }

/******** workers ********/
function runWorker(step, run, cfg) {
	return new Promise((resolve, reject) => {
		const w = curWorker = new Worker('cu-worker.js');
		curReject = reject; curIsAlign = false;
		w.onmessage = (e) => {
			const m = e.data;
			if (m.type === 'log') { log(m.text, m.cls); return; }
			if (m.type === 'progress') { if (m.frac !== undefined) $('prog').value = m.frac; $('run-msg').textContent = `${step}${m.label ? ' · ' + m.label : ''}${m.reads ? ' · ' + fmtInt(m.reads) + ' reads' : ''}`; return; }
			w.terminate(); curWorker = null;
			if (m.type === 'done') resolve(m.result); else { if (m.stack) console.error(m.stack); reject(new Error(m.message)); }
		};
		w.onerror = (ev) => { w.terminate(); curWorker = null; reject(new Error(ev.message || 'worker error')); };
		w.postMessage({ type: 'run', step, run, cfg });
	});
}
function runAlign({ job, inputs, reads, args, idx, batch, label }) {
	const E = engine(), W = alignWasm();
	let peak = 0;
	return new Promise((resolve, reject) => {
		const w = curWorker = new Worker('align-worker.js');
		curReject = reject; curIsAlign = true;
		w.onmessage = (e) => {
			const x = e.data;
			if (x.mem) { peak = Math.max(peak, x.mem); gauge(x.mem, peak, `align worker (${E.name})`); }
			if (x.type === 'log') { log(x.text.replace(/\n$/, '')); return; }
			if (x.type === 'loaded') { log(`[${label}] index loaded in ${fmtSec(x.sec)}`); return; }
			if (x.type === 'progress') { $('prog').value = x.inFrac; $('run-msg').textContent = `${label} · ${fmtInt(x.records)} records · ${fmtInt(x.records / Math.max(x.mapSec, 1e-3))}/s`; return; }
			if (x.type === 'jobdone') return;
			w.terminate(); curWorker = null; gauge(0, peak, null);
			if (x.type === 'done') resolve({ ...x, wasm: W.label }); else reject(new Error(x.message));
		};
		w.onerror = (ev) => { w.terminate(); curWorker = null; gauge(0, peak, null); reject(new Error(ev.message || 'worker error')); };
		w.postMessage({ type: 'start', jobsDir: 'cu-jobs', job, reads, inputs, index: idx.name, indexSource: idx.source, indexFiles: idx.files, args, batch,
			cacheMb: E.m64 ? 1024 : 256, wasmUrl: W.url, threads: W.threads });
	});
}
async function jobText(job) { // header + chunks of a small alignment job as one string (job folder removed)
	const jt = await jobsTop(), jd = await jt.getDirectoryHandle(job);
	const names = (await listDir(jd)).filter((f) => f.name.startsWith('chunk.')).map((f) => f.name).sort();
	let t = await (await (await jd.getFileHandle('header')).getFile()).text();
	for (const n of names) t += await (await (await jd.getFileHandle(n)).getFile()).text();
	await jt.removeEntry(job, { recursive: true });
	return t;
}
const splitOpts = (s) => s.trim().split(/\s+/).filter(Boolean);

async function stepPrep(run, S, idx, cfg) {
	const fq = [['primer', S.primer_seq], ...S.guide_sequences.map((g, i) => ['guide' + i, g])].map(([n, s]) => `@${n}\n${s}\n+\n${'I'.repeat(s.length)}\n`).join('');
	const job = `${run}.locate`;
	try { await (await jobsTop()).removeEntry(job, { recursive: true }); } catch (_) {}
	const args = [...splitOpts($('opts-loc').value), '-b', 'MD'];
	await runAlign({ job, inputs: { 'locate.fq': new File([fq], 'locate.fq') }, reads: ['locate.fq'], args, idx, batch: 1000, label: 'locate' });
	const sam = await jobText(job);
	for (const l of sam.split('\n')) if (l && l[0] !== '@') log('[locate] ' + l.split('\t').slice(0, 6).join('\t'));
	return runWorker('prep', run, { ...cfg, locateSam: sam });
}
async function stepAlign(run, S, idx, cfg) {
	const dir = await runDirH(run);
	const st = JSON.parse((await (await fileIn(dir, S.root + '.web_state.json'))?.text()) || '{}');
	if (!st.aligned_inputs) throw new Error('Run the primer step first (no trimmed FASTQ)');
	const inputs = {};
	for (const n of st.aligned_inputs) { const f = await fileIn(dir, n); if (!f) throw new Error(n + ' not found — re-run the primer step'); inputs[n] = f; }
	const job = `${run}.genome`;
	try { await (await jobsTop()).removeEntry(job, { recursive: true }); } catch (_) {}
	const args = [...splitOpts($('opts').value), '-b', 'MD'];
	if (S.analyze_multimap_assignments) args.push(`--outn=${Math.max(0, S.max_multimap_alignments - 1)}`);
	const t0 = performance.now();
	const r = await runAlign({ job, inputs, reads: st.aligned_inputs, args, idx, batch: Math.max(1000, parseInt($('batch').value, 10) || 20000), label: 'align' });
	const alignerLog = `\nAligning ${st.aligned_inputs.length === 2 ? 'paired' : 'single-end'} reads\n===\nCommand used:\n===\nminibwa map ${args.join(' ')} ${idx.name} ${st.aligned_inputs.join(' ')} (WebAssembly, ${r.wasm})\n===\nOutput:\n===\n${fmtInt(r.records)} SAM records in ${fmtSec((performance.now() - t0) / 1000)}\n`;
	log(`[align] ${fmtInt(r.records)} records, index load ${fmtSec(r.loadSec)}, mapping ${fmtSec(r.mapSec)}`, 'ok');
	return runWorker('sort', run, { ...cfg, jobsDir: 'cu-jobs', job, alignerLog });
}

async function runSteps() {
	if (busy) return;
	let S, warnings;
	try { ({ settings: S, warnings } = parsed()); } catch (e) { $('run-msg').textContent = 'Settings error: ' + e.message; return; }
	const run = runName() || sanitizeName(S.root);
	S.root = run;
	const want = STEPS.map(([k]) => k).filter((k) => $('chk-' + k).checked);
	if (!want.length) { $('run-msg').textContent = 'Select the steps to run.'; return; }
	let idx = null;
	if (want.includes('prep') || want.includes('align') || !refFiles().fasta) { try { idx = await selectedIndex(); } catch (e) { if (want.includes('prep') || want.includes('align')) { $('run-msg').textContent = e.message; return; } } }
	busy = true; cancelled = false; $('btn-run').disabled = true; $('btn-cancel').disabled = false; $('prog').classList.remove('hidden');
	try { if (navigator.storage.persist) await navigator.storage.persist(); } catch (_) {}
	const dir = await runDirH(run);
	await writeText(dir, S.root + '.settings.txt', $('settings').value);
	await writeText(dir, S.root + '.settingsUsed.txt', CuSettings.settingsUsed(S));
	const R = refFiles();
	const reference = R.fasta ? { fasta: R.fasta, fai: R.fai } : idx ? { l2b: idx.source === 'file' ? idx.files.l2b : idx.name + '.l2b' } : null;
	const cfg = { settings: S, files: fileMap(), reference, fixSpaceIds: $('fix-ids').checked };
	const T0 = performance.now();
	log(`=== run ${run}: ${want.join(' → ')} ===`, 'ok');
	for (const w of warnings) log('[settings] ' + w, 'warn');
	for (const k of want) setStatus(k, 'queued');
	let cur = null;
	try {
		for (const k of want) {
			if (cancelled) throw new Error('cancelled');
			cur = k; setStatus(k, 'running…'); $('prog').value = 0;
			const t0 = performance.now();
			if (k === 'prep') await stepPrep(run, S, idx, cfg);
			else if (k === 'align') await stepAlign(run, S, idx, cfg);
			else await runWorker(k, run, cfg);
			setStatus(k, `<span class="ok">done ${fmtSec((performance.now() - t0) / 1000)}</span>`);
		}
		$('run-msg').textContent = `Done: ${fmtSec((performance.now() - T0) / 1000)}`;
		$('prog').value = 1;
	} catch (e) {
		if (cur) setStatus(cur, `<span class="err">${e.message.includes('cancelled') ? 'cancelled' : 'failed'}</span>`);
		$('run-msg').textContent = e.message.includes('cancelled') ? 'Cancelled' : 'Failed — see the log';
		log('[run] ' + e.message, 'err');
	} finally {
		busy = false; $('btn-run').disabled = false; $('btn-cancel').disabled = true;
		await refreshRuns(run);
		checkEnv();
	}
}
function setStatus(k, html) { $('st-' + k).innerHTML = html; }
$('btn-cancel').addEventListener('click', () => {
	cancelled = true;
	if (curWorker && curIsAlign) curWorker.postMessage({ type: 'cancel' });
	else if (curWorker) { curWorker.terminate(); curWorker = null; if (curReject) curReject(new Error('cancelled')); }
	$('run-msg').textContent = 'Cancelling…';
});

/******** results ********/
async function walk(dir, prefix = '') {
	const out = [];
	for (const e of await listDir(dir)) {
		if (e.kind === 'directory') out.push(...await walk(e.h, prefix + e.name + '/'));
		else out.push({ path: prefix + e.name, h: e.h, size: e.size });
	}
	return out;
}
async function refreshRuns(select) {
	const runs = (await listDir(await cuTop())).filter((e) => e.kind === 'directory').map((e) => e.name);
	const sel = $('run-sel'), prev = select || sel.value;
	sel.innerHTML = runs.map((r) => `<option>${esc(r)}</option>`).join('') || '<option value="">(none)</option>';
	if (prev && runs.includes(prev)) sel.value = prev;
	await renderResults();
}
let reportUrl = null;
async function renderResults() {
	const run = $('run-sel').value;
	$('files').tBodies[0].innerHTML = ''; $('summary-box').textContent = ''; $('report-frame').classList.add('hidden');
	if (!run) return;
	const files = await walk(await runDirH(run));
	$('files').tBodies[0].innerHTML = files.map((f, i) => `<tr><td class="mono">${esc(f.path)}</td><td class="num">${fmtBytes(f.size)}</td><td class="num"><button data-i="${i}">Download</button></td></tr>`).join('');
	$('files').onclick = (e) => { const i = e.target.dataset.i; if (i !== undefined) download(files[+i].h, files[+i].path.split('/').pop()); };
	const sum = files.find((f) => f.path.endsWith('.summary.txt') && !f.path.includes('CRISPResso'));
	if (sum) { const t = await (await sum.h.getFile()).text(); $('summary-box').textContent = t.split('\n').slice(2).join('\n').trim(); }
	const rep = files.find((f) => f.path.endsWith('.html'));
	if (rep) {
		const html = await (await rep.h.getFile()).text();
		if (reportUrl) URL.revokeObjectURL(reportUrl);
		reportUrl = URL.createObjectURL(new Blob([html], { type: 'text/html' }));
		$('report-frame').srcdoc = html; $('report-frame').classList.remove('hidden');
	}
}
$('run-sel').addEventListener('change', renderResults);
$('btn-report').addEventListener('click', () => { if (reportUrl) window.open(reportUrl, '_blank'); });
$('btn-zip').addEventListener('click', async () => {
	const run = $('run-sel').value; if (!run) return;
	const files = (await walk(await runDirH(run))).filter((f) => !/\.(fq|sam)$/.test(f.path) && !f.path.startsWith('tmp/') && f.size < 512e6);
	const z = {};
	for (const f of files) z[run + '/' + f.path] = new Uint8Array(await (await f.h.getFile()).arrayBuffer());
	download(new Blob([fflate.zipSync(z, { level: 6 })], { type: 'application/zip' }), run + '_CRISPRuno.zip');
});
$('btn-del-run').addEventListener('click', async () => {
	const run = $('run-sel').value;
	if (!run || busy || !confirm(`Delete all files of run "${run}"?`)) return;
	await (await cuTop()).removeEntry(run, { recursive: true });
	await refreshRuns(); checkEnv();
});

/******** start ********/
(async () => {
	initPersistButton();
	$('steps').tBodies[0].innerHTML = STEPS.map(([k, name, orig]) => `<tr><td><input type="checkbox" id="chk-${k}" checked></td><td>${esc(name)}</td><td class="note mono" style="margin:0">${esc(orig)}</td><td id="st-${k}"></td></tr>`).join('');
	let t = null; try { t = localStorage.getItem('cu-settings'); } catch (_) {}
	$('settings').value = t || CuSettings.TEMPLATE;
	OPT = CuOptions.build($('options'), { getText: () => $('settings').value, setText: (t) => { $('settings').value = t; renderSettings(); }, defaults: CuSettings.DEFAULTS });
	$('settings').addEventListener('input', renderSettings);
	$('run-name').addEventListener('input', renderSettings);
	$('settings-file').addEventListener('change', async () => { const f = $('settings-file').files[0]; if (f) { $('settings').value = await f.text(); if (!$('run-name').value) $('run-name').value = f.name.replace(/\.[^.]*$/, '') + '.CRISPRuno'; renderSettings(); } });
	for (const id of ['fastqs', 'fasta']) $(id).addEventListener('change', renderSettings);
	document.querySelectorAll('input[name=src]').forEach((r) => r.addEventListener('change', () => {
		const file = document.querySelector('input[name=src]:checked').value === 'file';
		$('idx-files').classList.toggle('hidden', !file); $('opfs-idx').classList.toggle('hidden', file);
	}));
	$('btn-run').addEventListener('click', runSteps);
	$('threads').value = String(Math.min(CORES, 16));
	const onThreads = () => { $('threads-note').textContent = !(MT_OK && engine().m64) ? 'Multi-threading unavailable — aligning with a single thread' : `using ${alignWasm().threads} of ${CORES} CPU cores`; };
	$('threads').addEventListener('change', onThreads);
	initEngine(() => { $('threads').disabled = !(MT_OK && engine().m64); onThreads(); });
	try { const jt = await jobsTop(); for await (const [n] of jt.entries()) await jt.removeEntry(n, { recursive: true }); } catch (_) {}
	renderSettings();
	await checkEnv();
	await refreshIndexes();
	await refreshRuns();
})().catch((e) => log(String(e), 'err'));
