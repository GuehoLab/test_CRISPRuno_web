'use strict';
/* minibwa Index page: FASTA -> .l2b + .mbw in OPFS (one worker per stage, see worker.js), then save to disk. */
let busy = false;

function updateSaNote() {
	const E = engine(), sa = +$('sa-bit').value;
	const t2t = estimateMapMem(6234551004, sa, E.m64 ? 1024 : 256);
	$('sa-note').textContent = `WASM memory for aligning to a human genome (T2T) ≈ ${fmtBytes(t2t)} (${E.name}, ${E.m64 ? 1024 : 256} MB cache)` +
		(t2t > E.limit ? ` — exceeds the ${E.name} limit. Build with 1/128 or use wasm64.` : '') +
		' · SA sampling does not change alignment results, only speed and memory.';
	$('sa-note').className = 'note' + (t2t > E.limit ? ' warn' : '');
}

/******** stored indexes ********/
async function refreshIndexes(highlight) {
	const info = await opfsIndexes(), E = engine(), cache = E.m64 ? 1024 : 256;
	$('idx-table').tBodies[0].innerHTML = info.length ? info.map((x) => {
		const est = x.ok ? estimateMapMem(x.seqLen, x.saBit, cache) : 0, over = est > E.limit;
		const m = x.meta ? `${(x.meta.date || '').slice(0, 10)} · ${x.meta.engine || ''} · ${x.meta.totalSec ? fmtSec(x.meta.totalSec) : ''}` : 'imported index';
		return `<tr${x.name === highlight ? ' class="ok"' : ''}><td class="mono">${esc(x.name)}</td><td class="num">${(x.genome / 1e6).toFixed(1)} Mb</td>` +
			`<td class="num">${x.ok ? '1/' + 2 ** x.saBit : '?'}</td><td class="num">${fmtBytes(x.bytes)}</td>` +
			`<td class="num ${over ? 'err' : ''}">${x.ok ? fmtBytes(est) + (over ? ` (exceeds ${E.name} limit)` : '') : 'invalid format'}</td><td class="note">${esc(m)}</td>` +
			`<td class="num"><button class="primary" data-act="folder" data-name="${esc(x.name)}">Save to folder</button> ` +
			`<button data-act="l2b" data-name="${esc(x.name)}">.l2b</button> <button data-act="mbw" data-name="${esc(x.name)}">.mbw</button> ` +
			`<button data-act="rm" data-name="${esc(x.name)}">Delete</button></td></tr>`;
	}).join('') : '<tr><td colspan="7" class="note">None — build one above.</td></tr>';
	return info;
}

async function saveToFolder(x) {
	if (!window.showDirectoryPicker) { // Firefox/Safari: one download per file
		for (const ext of ['l2b', 'mbw']) await download(x.files[ext], `${x.name}.${ext}`);
		return;
	}
	let dir;
	try { dir = await window.showDirectoryPicker({ id: 'mbw-index', mode: 'readwrite' }); } catch (e) { if (e.name === 'AbortError') return; throw e; }
	const list = [[x.name + '.l2b', x.files.l2b], [x.name + '.mbw', x.files.mbw]];
	if (x.files.json) list.push([x.name + '.build.json', x.files.json]);
	const files = await Promise.all(list.map(async ([n, h]) => [n, await h.getFile()]));
	const total = files.reduce((a, [, f]) => a + f.size, 0);
	let done = 0;
	const t0 = performance.now();
	$('save-box').classList.remove('hidden');
	for (const [n, f] of files) {
		const w = await (await dir.getFileHandle(n, { create: true })).createWritable();
		await streamTo(f, w, (b) => {
			done += b;
			$('save-prog').value = done / total;
			$('save-msg').textContent = `Saving ${dir.name}/${n} · ${fmtBytes(done)} / ${fmtBytes(total)}`;
		});
	}
	$('save-msg').textContent = `Saved ${files.map(([n]) => n).join(', ')} to ${dir.name}/ (${fmtBytes(total)}, ${fmtSec((performance.now() - t0) / 1000)})`;
	log(`[save] ${x.name} → ${dir.name}/ (${fmtBytes(total)})`, 'ok');
}

$('idx-table').addEventListener('click', async (e) => {
	const { act, name } = e.target.dataset;
	if (!act || busy) return;
	const x = (await opfsIndexes()).find((i) => i.name === name);
	if (!x) return;
	try {
		if (act === 'folder') await saveToFolder(x);
		if (act === 'l2b' || act === 'mbw') await download(x.files[act], `${name}.${act}`);
		if (act === 'rm' && confirm(`Delete index ${name} from this browser's storage? (Downloaded files are not affected.)`)) { await removeIndex(name); await refreshIndexes(); checkEnv(); }
	} catch (err) { log(`[save] ${err.message}`, 'err'); $('save-msg').textContent = 'Failed: ' + err.message; }
});

/******** index build (one worker per stage, see worker.js) ********/
function runStage(label, msg) {
	return new Promise((resolve, reject) => {
		const w = new Worker('worker.js'), W = stageWasm();
		let peak = 0;
		const t0 = performance.now();
		w.onmessage = (e) => {
			const m = e.data;
			if (m.mem) { peak = Math.max(peak, m.mem); gauge(m.mem, peak, label); }
			if (m.type === 'log') { log(m.text.replace(/\n$/, '')); return; }
			w.terminate(); gauge(0, peak, null); // terminate = WASM memory released
			if (m.type === 'done' && m.rc === 0) resolve({ ...m, peak, wall: (performance.now() - t0) / 1000 });
			else reject(new Error(m.type === 'error' ? m.message : `exit code ${m.rc}`));
		};
		w.onerror = (ev) => { w.terminate(); reject(new Error(ev.message || 'worker error')); };
		log(`$ minibwa ${msg.args.join(' ')}  [${W.label}]`, 'warn');
		gauge(0, 0, label);
		w.postMessage({ wasmUrl: W.url, threads: W.threads, opfsDir: 'minibwa', ...msg });
	});
}
$('fa').addEventListener('change', () => {
	const f = $('fa').files[0];
	if (f && !$('idx-name').value) $('idx-name').value = sanitizeName(f.name.replace(/\.gz$/, '').replace(/\.(fa|fasta|fna)$/, ''));
});
$('sa-bit').addEventListener('change', updateSaNote);
$('btn-build').addEventListener('click', async () => {
	const fa = $('fa').files[0], name = sanitizeName($('idx-name').value);
	if (!fa || !name || busy) return;
	if ((await opfsIndexes()).some((x) => x.name === name) && !confirm(`Index ${name} already exists. Overwrite it?`)) return;
	const sa = $('sa-bit').value, estBp = fa.name.endsWith('.gz') ? fa.size * 3.3 : fa.size, E = engine(), W = stageWasm(), T = String(W.threads);
	let method = $('method').value;
	if (method === 'auto') method = estBp > 300e6 ? 'staged' : 'inmem';
	busy = true; $('btn-build').disabled = true; $('stages').tBodies[0].innerHTML = '';
	const P = '/opfs/' + name, inputs = { [fa.name]: fa }, stats = [], T0 = performance.now();
	const steps = method === 'staged' ? [
		['prep: FASTA → .l2b + .pac', { args: ['prep', '/input/' + fa.name, P], write: [name + '.l2b', name + '.pac'], inputs }],
		['genraw: BWT (bwtgen)', { args: ['genraw', '-b', '100m', '-t', T, P + '.pac', P + '.raw'], read: [name + '.pac'], write: [name + '.raw'] }, [name + '.pac']],
		['raw2mbw: BWT encoding + SA', { args: ['raw2mbw', '-u', sa, '-t', T, P + '.raw', P + '.mbw'], read: [name + '.raw'], write: [name + '.mbw'] }, [name + '.raw']],
	] : [['index (libsais)', { args: ['index', '-u', sa, '/input/' + fa.name, P], write: [name + '.l2b', name + '.mbw'], inputs }]];
	try {
		if (navigator.storage.persist) await navigator.storage.persist();
		try { await (await opfsTop()).removeEntry(name + '.build.json'); } catch (_) {}
		for (const [label, msg, cleanup] of steps) {
			const tr = document.createElement('tr');
			tr.innerHTML = `<td>${esc(label)}</td><td class="warn">running…</td><td class="num"></td><td class="num"></td>`;
			$('stages').tBodies[0].appendChild(tr);
			$('build-msg').textContent = label;
			try {
				const r = await runStage(label, msg);
				tr.cells[1].className = 'ok'; tr.cells[1].textContent = 'done';
				tr.cells[2].textContent = fmtSec(r.wall); tr.cells[3].textContent = 'peak ' + fmtBytes(r.peak);
				stats.push({ stage: label, sec: r.wall, peakBytes: r.peak });
			} catch (err) { tr.cells[1].className = 'err'; tr.cells[1].textContent = 'failed'; throw err; }
			const dir = await opfsTop();
			for (const f of cleanup || []) { try { await dir.removeEntry(f); } catch (_) {} }
		}
		const total = (performance.now() - T0) / 1000;
		await writeText(await opfsTop(), name + '.build.json', JSON.stringify({ name, fasta: fa.name, fastaBytes: fa.size, method, saBit: +sa, engine: E.name, threads: W.threads,
			totalSec: total, stages: stats, userAgent: navigator.userAgent, date: new Date().toISOString() }, null, 2));
		$('build-msg').textContent = `Done: ${fmtSec(total)} — save it to a folder under "Saved indexes" below.`;
		log(`[index] ${name} done in ${fmtSec(total)}`, 'ok');
	} catch (err) {
		$('build-msg').textContent = 'Failed — see the log';
		log('[index] ' + err.message, 'err');
	} finally { busy = false; $('btn-build').disabled = false; refreshIndexes(name); checkEnv(); }
});

/******** start ********/
(async () => {
	initPersistButton();
	$('threads').value = String(Math.min(CORES, 16));
	initEngine(() => {
		$('sa-bit').value = engine().m64 ? '4' : '7';
		$('threads').disabled = !(MT_OK && engine().m64);
		$('threads-note').textContent = MT_OK && engine().m64 ? `${CORES} CPU cores. Part of the BWT construction and the SA sampling run in parallel; the output files are the same for any number of threads.`
			: 'Multi-threading unavailable (see the environment table) — building with a single thread.';
		updateSaNote(); refreshIndexes();
	});
	await checkEnv();
})().catch((e) => log(String(e), 'err'));
