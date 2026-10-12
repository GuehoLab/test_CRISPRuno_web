'use strict';
/* minibwa Index page: FASTA -> .l2b + .mbw in OPFS (one worker per stage, see worker.js), then save to disk. */
let busy = false;

function updateSaNote() {
	const E = engine(), sa = +$('sa-bit').value;
	const t2t = estimateMapMem(6234551004, sa, E.m64 ? 1024 : 256);
	$('sa-note').textContent = `인간 게놈(T2T) 기준 정렬 시 WASM 메모리 ≈ ${fmtBytes(t2t)} (${E.name}, 캐시 ${E.m64 ? 1024 : 256} MB)` +
		(t2t > E.limit ? ` — ${E.name} 한도 초과. 1/128로 만들거나 wasm64를 쓰세요.` : '') +
		' · SA 샘플링은 정렬 결과에는 영향이 없고 속도와 메모리만 바꿉니다.';
	$('sa-note').className = 'note' + (t2t > E.limit ? ' warn' : '');
}

/******** stored indexes ********/
async function refreshIndexes(highlight) {
	const info = await opfsIndexes(), E = engine(), cache = E.m64 ? 1024 : 256;
	$('idx-table').tBodies[0].innerHTML = info.length ? info.map((x) => {
		const est = x.ok ? estimateMapMem(x.seqLen, x.saBit, cache) : 0, over = est > E.limit;
		const m = x.meta ? `${(x.meta.date || '').slice(0, 10)} · ${x.meta.engine || ''} · ${x.meta.totalSec ? fmtSec(x.meta.totalSec) : ''}` : '가져온 인덱스';
		return `<tr${x.name === highlight ? ' class="ok"' : ''}><td class="mono">${esc(x.name)}</td><td class="num">${(x.genome / 1e6).toFixed(1)} Mb</td>` +
			`<td class="num">${x.ok ? '1/' + 2 ** x.saBit : '?'}</td><td class="num">${fmtBytes(x.bytes)}</td>` +
			`<td class="num ${over ? 'err' : ''}">${x.ok ? fmtBytes(est) + (over ? ` (${E.name} 한도 초과)` : '') : '형식 오류'}</td><td class="note">${esc(m)}</td>` +
			`<td class="num"><button class="primary" data-act="folder" data-name="${esc(x.name)}">폴더에 저장</button> ` +
			`<button data-act="l2b" data-name="${esc(x.name)}">.l2b</button> <button data-act="mbw" data-name="${esc(x.name)}">.mbw</button> ` +
			`<button data-act="rm" data-name="${esc(x.name)}">삭제</button></td></tr>`;
	}).join('') : '<tr><td colspan="7" class="note">없음 — 위에서 만드세요.</td></tr>';
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
			$('save-msg').textContent = `${dir.name}/${n} 저장 중 · ${fmtBytes(done)} / ${fmtBytes(total)}`;
		});
	}
	$('save-msg').textContent = `저장 완료: ${dir.name}/ 에 ${files.map(([n]) => n).join(', ')} (${fmtBytes(total)}, ${fmtSec((performance.now() - t0) / 1000)})`;
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
		if (act === 'rm' && confirm(`이 브라우저 저장소에서 인덱스 ${name}을(를) 삭제할까요? (내려받은 파일은 그대로)`)) { await removeIndex(name); await refreshIndexes(); checkEnv(); }
	} catch (err) { log(`[save] ${err.message}`, 'err'); $('save-msg').textContent = '실패: ' + err.message; }
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
	if ((await opfsIndexes()).some((x) => x.name === name) && !confirm(`${name} 인덱스가 이미 있습니다. 덮어쓸까요?`)) return;
	const sa = $('sa-bit').value, estBp = fa.name.endsWith('.gz') ? fa.size * 3.3 : fa.size, E = engine(), W = stageWasm(), T = String(W.threads);
	let method = $('method').value;
	if (method === 'auto') method = estBp > 300e6 ? 'staged' : 'inmem';
	busy = true; $('btn-build').disabled = true; $('stages').tBodies[0].innerHTML = '';
	const P = '/opfs/' + name, inputs = { [fa.name]: fa }, stats = [], T0 = performance.now();
	const steps = method === 'staged' ? [
		['prep: FASTA → .l2b + .pac', { args: ['prep', '/input/' + fa.name, P], write: [name + '.l2b', name + '.pac'], inputs }],
		['genraw: BWT (bwtgen)', { args: ['genraw', '-b', '100m', '-t', T, P + '.pac', P + '.raw'], read: [name + '.pac'], write: [name + '.raw'] }, [name + '.pac']],
		['raw2mbw: BWT 인코딩 + SA', { args: ['raw2mbw', '-u', sa, '-t', T, P + '.raw', P + '.mbw'], read: [name + '.raw'], write: [name + '.mbw'] }, [name + '.raw']],
	] : [['index (libsais)', { args: ['index', '-u', sa, '/input/' + fa.name, P], write: [name + '.l2b', name + '.mbw'], inputs }]];
	try {
		if (navigator.storage.persist) await navigator.storage.persist();
		try { await (await opfsTop()).removeEntry(name + '.build.json'); } catch (_) {}
		for (const [label, msg, cleanup] of steps) {
			const tr = document.createElement('tr');
			tr.innerHTML = `<td>${esc(label)}</td><td class="warn">실행 중…</td><td class="num"></td><td class="num"></td>`;
			$('stages').tBodies[0].appendChild(tr);
			$('build-msg').textContent = label;
			try {
				const r = await runStage(label, msg);
				tr.cells[1].className = 'ok'; tr.cells[1].textContent = '완료';
				tr.cells[2].textContent = fmtSec(r.wall); tr.cells[3].textContent = '최대 ' + fmtBytes(r.peak);
				stats.push({ stage: label, sec: r.wall, peakBytes: r.peak });
			} catch (err) { tr.cells[1].className = 'err'; tr.cells[1].textContent = '실패'; throw err; }
			const dir = await opfsTop();
			for (const f of cleanup || []) { try { await dir.removeEntry(f); } catch (_) {} }
		}
		const total = (performance.now() - T0) / 1000;
		await writeText(await opfsTop(), name + '.build.json', JSON.stringify({ name, fasta: fa.name, fastaBytes: fa.size, method, saBit: +sa, engine: E.name, threads: W.threads,
			totalSec: total, stages: stats, userAgent: navigator.userAgent, date: new Date().toISOString() }, null, 2));
		$('build-msg').textContent = `완료: ${fmtSec(total)} — 아래 "저장된 인덱스"에서 폴더에 저장하세요.`;
		log(`[index] ${name} done in ${fmtSec(total)}`, 'ok');
	} catch (err) {
		$('build-msg').textContent = '실패 — 로그 확인';
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
		$('threads-note').textContent = MT_OK && engine().m64 ? `CPU 코어 ${CORES}개. BWT 생성의 일부와 SA 샘플링이 병렬로 돌고, 결과 파일은 스레드 수와 무관하게 같습니다.`
			: '멀티스레드 불가 (환경 표 참고) — 단일 스레드로 만듭니다.';
		updateSaNote(); refreshIndexes();
	});
	await checkEnv();
})().catch((e) => log(String(e), 'err'));
