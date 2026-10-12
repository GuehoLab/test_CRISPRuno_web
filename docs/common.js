'use strict';
/* Shared by build.html (index) and align.html (alignment).
 * OPFS layout (directory "minibwa"):
 *   <index>.l2b / <index>.mbw / <index>.build.json   reference indexes
 *   jobs/<job>/header, chunk.000000 ...             per-batch alignment output (temporary)
 *   results/<job>.sam|paf (+ .json summary)         merged results
 */
const GIB = 1024 ** 3;
const $ = (id) => document.getElementById(id);
const fmtBytes = (b) => b >= 1e9 ? (b / 1e9).toFixed(2) + ' GB' : b >= 1e6 ? (b / 1e6).toFixed(1) + ' MB' : b >= 1e3 ? (b / 1e3).toFixed(1) + ' kB' : b + ' B';
const fmtSec = (s) => { if (s < 59.95) return `${s.toFixed(1)} s`; const t = Math.round(s); return t >= 3600 ? `${Math.floor(t / 3600)}h ${Math.floor((t % 3600) / 60)}m` : `${Math.floor(t / 60)}m ${t % 60}s`; };
const fmtInt = (x) => Math.round(x).toLocaleString();
const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
const sanitizeName = (s) => String(s || '').replace(/[^A-Za-z0-9._-]/g, '_');

/******** engine (wasm32 / wasm64) ********/
const MEM64 = (() => { try { return WebAssembly.validate(new Uint8Array([0, 0x61, 0x73, 0x6d, 1, 0, 0, 0, 5, 3, 1, 4, 1])); } catch (_) { return false; } })();
// wasm64 (Memory64, minibwa64.wasm) lifts the 4 GiB limit of wasm32 so that a human index can keep SA 1/16.
function engine() {
	const sel = $('engine'), v = sel ? sel.value : 'auto', m64 = MEM64 && (v === 'wasm64' || v === 'auto');
	return m64 ? { m64, url: 'minibwa64.wasm', limit: 16 * GIB, name: 'wasm64', limitText: '17.2 GB (16 GiB, wasm64)' }
		: { m64, url: 'minibwa.wasm', limit: 4 * GIB, name: 'wasm32', limitText: '4.29 GB (4 GiB, wasm32)' };
}
// Multi-threading needs shared memory (SharedArrayBuffer), which browsers allow only on cross-origin isolated
// pages (COOP/COEP headers: serve.py sends them; on GitHub Pages coi-serviceworker.js emulates them).
// Only the wasm64 build has a multi-threaded variant (minibwa64-mt.js); it is used for the index build.
const MT_OK = MEM64 && self.crossOriginIsolated === true && typeof SharedArrayBuffer !== 'undefined';
const CORES = Math.max(1, navigator.hardwareConcurrency || 1);
function threadsWanted() { // value of the #threads input (build page), 1 when multi-threading is unavailable
	const el = $('threads');
	if (!MT_OK || !engine().m64) return 1;
	return Math.max(1, Math.min(64, parseInt(el ? el.value : CORES, 10) || 1));
}
function stageWasm() { // wasm for index stages: multi-threaded wasm64 when possible
	const E = engine(), n = threadsWanted();
	return n > 1 ? { url: 'minibwa64-mt.js', threads: n, label: `${E.name}, ${n} threads` } : { url: E.url, threads: 1, label: `${E.name}, 1 thread` };
}
function alignWasm() { // wasm for alignment: multi-threaded wasm64 when possible (whole reference in memory)
	const E = engine(), n = threadsWanted();
	return n > 1 ? { url: 'minibwa64-mt.js', threads: n, label: `${E.name}, ${n} threads` } : { url: E.url, threads: 1, label: `${E.name}, 1 thread` };
}
const MT_THREAD_BYTES = 16e6; // per mapping thread in the multi-threaded build (2 MB stack + buffers; measured ~6 MB at batch 10000 x 150 bp)
function initEngine(onChange) { // the choice is shared by both pages (localStorage, same origin)
	const sel = $('engine');
	try { const v = localStorage.getItem('mbw-engine'); if (v && (v !== 'wasm64' || MEM64)) sel.value = v; } catch (_) {}
	if (!MEM64) sel.querySelector('option[value=wasm64]').disabled = true;
	sel.addEventListener('change', () => {
		try { localStorage.setItem('mbw-engine', sel.value); } catch (_) {}
		gauge(0, 0, null);
		if (onChange) onChange();
	});
	gauge(0, 0, null);
	if (onChange) onChange();
}

/******** log & gauge ********/
function log(text, cls) {
	const el = $('log'), span = document.createElement('span');
	if (cls) span.className = cls;
	span.textContent = text.endsWith('\n') ? text : text + '\n';
	el.appendChild(span);
	while (el.childNodes.length > 3000) el.removeChild(el.firstChild);
	el.scrollTop = el.scrollHeight;
}
function gauge(cur, peak, label) {
	const g = $('gauge'), E = engine();
	g.firstElementChild.style.width = Math.min(100, (cur / E.limit) * 100) + '%';
	g.lastElementChild.style.left = Math.min(100, (peak / E.limit) * 100) + '%';
	$('gauge-txt').textContent = label === null ? `실행 중인 worker 없음${peak ? ' (메모리 반환됨)' : ''} · 한도 ${E.limitText}`
		: `${label} · 현재 ${fmtBytes(cur)} · 최대 ${fmtBytes(peak)} · 한도 ${E.limitText}`;
}

/******** OPFS ********/
async function opfsTop() { return (await navigator.storage.getDirectory()).getDirectoryHandle('minibwa', { create: true }); }
async function sub(name) { return (await opfsTop()).getDirectoryHandle(name, { create: true }); }
async function listDir(dir) {
	const out = [];
	for await (const [name, h] of dir.entries()) out.push({ name, h, kind: h.kind, size: h.kind === 'file' ? (await h.getFile()).size : 0 });
	return out.sort((a, b) => a.name.localeCompare(b.name));
}
async function writeText(dir, name, text) { const w = await (await dir.getFileHandle(name, { create: true })).createWritable(); await w.write(text); await w.close(); }
async function streamTo(file, writable, onBytes) { // copy a File/Blob into a FileSystemWritableFileStream with progress
	await file.stream().pipeThrough(new TransformStream({ transform(chunk, ctl) { onBytes(chunk.byteLength); ctl.enqueue(chunk); } })).pipeTo(writable);
}
async function download(fh, name) { // one OPFS file -> user's disk (streamed in Chrome/Edge, normal download elsewhere)
	const file = fh instanceof Blob ? fh : await fh.getFile();
	if (window.showSaveFilePicker) {
		try {
			const h = await window.showSaveFilePicker({ suggestedName: name });
			await file.stream().pipeTo(await h.createWritable());
			log(`[save] ${name} → disk (${fmtBytes(file.size)})`, 'ok');
			return;
		} catch (e) { if (e.name === 'AbortError') return; }
	}
	const a = document.createElement('a');
	a.href = URL.createObjectURL(file); a.download = name; a.click();
	setTimeout(() => URL.revokeObjectURL(a.href), 600000);
}

/******** index files ********/
// .mbw header: "MBW\2", sa_bit (u32), primary, L2[1..4] (u64 each) -> seq_len = L2[4] at byte 40
// .l2b header: magic, dummy, n_ctg (u64 @8), tot_len (u64 @16)
async function readIndexHeaders(mbw, l2b) { // Blob/File objects
	const dv = new DataView(await mbw.slice(0, 48).arrayBuffer()), dl = new DataView(await l2b.slice(0, 24).arrayBuffer());
	const ok = dv.byteLength >= 48 && String.fromCharCode(dv.getUint8(0), dv.getUint8(1), dv.getUint8(2)) === 'MBW'
		&& dl.byteLength >= 24 && String.fromCharCode(dl.getUint8(0), dl.getUint8(1), dl.getUint8(2)) === 'L2B';
	const saBit = ok ? dv.getUint32(4, true) : -1, seqLen = ok ? Number(dv.getBigUint64(40, true)) : 0;
	const genome = ok ? Number(dl.getBigUint64(16, true)) : 0, nCtg = ok ? Number(dl.getBigUint64(8, true)) : 0;
	return { ok, saBit, seqLen, genome, nCtg, bytes: mbw.size + l2b.size };
}
function estimateMapMem(seqLen, saBit, cacheMb) {
	const data = (Math.floor((seqLen + 127) / 128) * 8 + 4) * 8;
	const sa = Math.floor((seqLen + 2 ** saBit) / 2 ** saBit) * 8;
	const pac = Math.ceil((seqLen - 1) / 2 / 4); // 2-bit reference of one strand; the page cache never grows beyond it
	return data + sa + 4 ** 10 * 32 + Math.min(cacheMb * 2 ** 20, pac) + 150e6; // + 10-mer cache, page cache, batch/working memory
}
async function opfsIndexes() { // complete indexes stored in OPFS
	const dir = await opfsTop(), files = await listDir(dir), byName = {};
	for (const f of files) if (f.kind === 'file') byName[f.name] = f;
	const out = [];
	for (const f of files) {
		if (!f.name.endsWith('.mbw') || !byName[f.name.slice(0, -4) + '.l2b']) continue;
		const name = f.name.slice(0, -4), l2b = byName[name + '.l2b'];
		const info = await readIndexHeaders(await f.h.getFile(), await l2b.h.getFile());
		let meta = null;
		if (byName[name + '.build.json']) { try { meta = JSON.parse(await (await byName[name + '.build.json'].h.getFile()).text()); } catch (_) {} }
		out.push({ name, ...info, meta, files: { l2b: l2b.h, mbw: f.h, json: byName[name + '.build.json'] ? byName[name + '.build.json'].h : null } });
	}
	return out;
}
async function removeIndex(name) {
	const dir = await opfsTop();
	for (const ext of ['.l2b', '.mbw', '.build.json']) { try { await dir.removeEntry(name + ext); } catch (_) {} }
}

/******** environment ********/
async function checkEnv() {
	const rows = [];
	const simd = WebAssembly.validate(new Uint8Array([0, 97, 115, 109, 1, 0, 0, 0, 1, 5, 1, 96, 0, 1, 123, 3, 2, 1, 0, 10, 10, 1, 8, 0, 65, 0, 253, 15, 253, 98, 11]));
	rows.push(['WebAssembly SIMD128', simd ? '지원' : '미지원 (실행 불가)', simd ? 'ok' : 'err']);
	rows.push(['WebAssembly Memory64 (wasm64)', MEM64 ? '지원 — 4 GiB 한도 없음 (최대 16 GiB)' : '미지원 — wasm32 (4 GiB)로 실행', MEM64 ? 'ok' : 'warn']);
	rows.push(['멀티스레드 (SharedArrayBuffer)', MT_OK ? `가능 — CPU 코어 ${CORES}개` : !MEM64 ? '미지원 (wasm64 필요) — 단일 스레드로 실행'
		: '불가 — 페이지가 cross-origin isolated가 아님 (COOP/COEP 헤더 필요). 단일 스레드로 실행', MT_OK ? 'ok' : 'warn']);
	const sah = await new Promise((res) => {
		try {
			const w = new Worker(URL.createObjectURL(new Blob(["postMessage(typeof FileSystemFileHandle !== 'undefined' && 'createSyncAccessHandle' in FileSystemFileHandle.prototype)"], { type: 'text/javascript' })));
			w.onmessage = (e) => { res(e.data); w.terminate(); };
			w.onerror = () => { res(false); w.terminate(); };
		} catch (_) { res(false); }
	});
	rows.push(['OPFS + SyncAccessHandle (Worker)', sah ? '지원' : '미지원 (실행 불가)', sah ? 'ok' : 'err']);
	if (navigator.deviceMemory) rows.push(['기기 메모리', `≥ ${navigator.deviceMemory} GB (최대 8로 표시됨) · 인간 게놈은 16 GB 이상 권장`, '']);
	if (navigator.storage && navigator.storage.estimate) {
		const est = await navigator.storage.estimate();
		rows.push(['저장소 사용량 / 할당량', `${fmtBytes(est.usage)} / ${fmtBytes(est.quota)}`, est.quota - est.usage < 10e9 ? 'warn' : '']);
	}
	if (navigator.storage && navigator.storage.persisted) rows.push(['persist', (await navigator.storage.persisted()) ? 'persisted' : '아님 (브라우저가 정리할 수 있음)', '']);
	$('env').tBodies[0].innerHTML = rows.map(([k, v, c]) => `<tr><th>${esc(k)}</th><td class="${c}">${esc(v)}</td></tr>`).join('');
}
function initPersistButton() {
	$('btn-persist').addEventListener('click', async () => {
		$('persist-msg').textContent = (navigator.storage.persist && await navigator.storage.persist()) ? 'persisted' : '브라우저가 거절함';
		checkEnv();
	});
}
