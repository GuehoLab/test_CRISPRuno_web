/* One minibwa stage per worker.
 * The main thread creates a new Worker for every stage and terminates it when
 * the stage ends: WASM linear memory cannot shrink, so this is how memory is
 * returned between stages. Files are exchanged through OPFS (dir "minibwa").
 *
 * message in : { args: [...], read: [opfs names], write: [opfs names], inputs: {name: File}, wasmUrl, threads }
 *   wasmUrl 'minibwa.wasm' / 'minibwa64.wasm': single-threaded standalone module
 *   wasmUrl 'minibwa64-mt.js': multi-threaded wasm64 build (emscripten glue + minibwa64-mt.wasm). Needs a
 *     cross-origin isolated page; `threads` thread workers are started (nested in this worker) before the
 *     stage runs, and are terminated together with this worker.
 * message out: { type: 'log', stream, text, mem } | { type: 'done', rc, sec, mem, io } | { type: 'error', message }
 */
importScripts('wasi-shim.js');
const { Wasi, OpfsDir, BlobDir } = self.MinibwaWasi;

onmessage = async (e) => {
	const { args, read = [], write = [], inputs = {}, wasmUrl = 'minibwa.wasm', opfsDir = 'minibwa', threads = 1 } = e.data;
	const handles = {};
	let inst = null, result;
	const mem = () => (inst ? inst.exports.memory.buffer.byteLength : 0);
	try {
		const root = await navigator.storage.getDirectory();
		const dir = await root.getDirectoryHandle(opfsDir, { create: true });
		for (const name of read) handles[name] = await (await dir.getFileHandle(name)).createSyncAccessHandle();
		for (const name of write) handles[name] = await (await dir.getFileHandle(name, { create: true })).createSyncAccessHandle();
		const opts = {
			args: ['minibwa', ...args],
			preopens: { '/opfs': new OpfsDir(handles), '/input': new BlobDir(inputs) },
			stdout: (text) => postMessage({ type: 'log', stream: 'stdout', text, mem: mem() }),
			stderr: (text) => postMessage({ type: 'log', stream: 'stderr', text, mem: mem() }),
		};
		let wasi;
		if (wasmUrl.endsWith('.js')) { // multi-threaded wasm64
			importScripts(wasmUrl);
			wasi = new Wasi({ ...opts, mem64: true });
			const M = await self.createMinibwa({ mbwIO: wasi, mbwThreads: threads, mainScriptUrlOrBlob: wasmUrl });
			inst = { exports: { memory: M.getMemory(), mbw_alloc: M._mbw_alloc, mbw_free: M._mbw_free, mbw_run: M._mbw_run } };
		} else {
			const mod = await WebAssembly.compileStreaming(fetch(wasmUrl));
			wasi = new Wasi({ ...opts, mem64: Wasi.isMem64(mod) });
			inst = await WebAssembly.instantiate(mod, wasi.getImports());
		}
		const t0 = performance.now();
		const rc = wasi.start(inst);
		result = { type: 'done', rc, sec: (performance.now() - t0) / 1000, mem: mem(), io: wasi.stats };
	} catch (err) {
		result = { type: 'error', message: String((err && err.stack) || err), mem: mem() };
	}
	for (const h of Object.values(handles)) { try { h.flush(); h.close(); } catch (_) {} } // release OPFS locks before reporting
	postMessage(result);
};
