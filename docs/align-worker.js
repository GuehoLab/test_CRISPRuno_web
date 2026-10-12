/* Persistent align worker.
 *
 * One worker = one alignment job:
 *   1. load minibwa.wasm once and keep the reference index in WASM memory
 *      (BWT + sampled SA; the 2-bit reference stays on OPFS behind a page cache)
 *   2. map the input `batch` records at a time; every batch is written to its own
 *      OPFS chunk file, so per-batch memory stays small and results are on disk
 *   3. report 'done'; the main thread then terminates this worker (frees ~4 GB for a
 *      human index with SA 1/128, ~6.7 GB with 1/16 on wasm64) and merges the chunks
 *      into the final SAM/PAF.
 * wasmUrl selects the engine: minibwa.wasm (wasm32, <= 4 GiB) or minibwa64.wasm (Memory64, <= 16 GiB), or
 * minibwa64-mt.js (multi-threaded wasm64, `threads` mapping threads; needs a cross-origin isolated page).
 * Multi-threaded mode loads the whole 2-bit reference into memory (-P; T2T 0.78 GB, no page cache), so the
 * mapping threads never touch files; reading, pairing statistics and SAM output stay on this thread.
 *
 * in : { type:'start', job, index, indexSource:'opfs'|'file', indexFiles:{l2b: File, mbw: File} (file only),
 *        args:[map options], inputs:{name: File}, reads:[names], batch, cacheMb, wasmUrl,
 *        jobsDir:'jobs' (OPFS directory under minibwa/), more:[{job, reads:[names]}] (optional: further read sets
 *        mapped with the same loaded index, e.g. the samples of a GUIDE-seq run; their files are in `inputs` too) }
 *      indexSource 'file': the user's .mbw is read in place (FileReaderSync, once, sequentially). The .l2b is first
 *      streamed into a temporary OPFS file of the job (T2T: 1 GB, a few seconds) because the reference page cache
 *      does many small random reads: ~0.15 ms each through a SyncAccessHandle vs ~1-1.5 ms through FileReaderSync.
 *      { type:'cancel' }
 * out: { type:'log', text } | { type:'progress', ... } | { type:'loaded', ... } | { type:'jobdone', job, ... } (per read set)
 *      | { type:'done', ... } | { type:'error', message }
 */
importScripts('wasi-shim.js');
const { Wasi, OpfsDir, BlobDir } = self.MinibwaWasi;

let cancel = false;
const pad = (k) => String(k).padStart(6, '0');

async function run(m) {
	const root = await navigator.storage.getDirectory();
	const top = await root.getDirectoryHandle('minibwa', { create: true });
	const jobsTop = await top.getDirectoryHandle(m.jobsDir || 'jobs', { create: true });
	let jobs = await jobsTop.getDirectoryHandle(m.job, { create: true });
	const handles = {}, blobs = new BlobDir(m.inputs);
	const opfs = new OpfsDir(handles);
	const open = async (dir, name, create) => (handles[name] = await (await dir.getFileHandle(name, { create })).createSyncAccessHandle());
	const close = (name) => { try { handles[name].flush(); handles[name].close(); } catch (_) {} delete handles[name]; };
	let inst = null;
	const mem = () => (inst ? inst.exports.memory.buffer.byteLength : 0);
	const fromFile = m.indexSource === 'file', mt = (m.wasmUrl || '').endsWith('.js');
	const preopens = { '/opfs': opfs, '/input': blobs };
	if (fromFile && mt) { // both files read once, sequentially, from the user's files
		preopens['/index'] = new BlobDir({ [m.index + '.l2b']: m.indexFiles.l2b, [m.index + '.mbw']: m.indexFiles.mbw });
	} else if (fromFile) { // .mbw from the user's file, .l2b from the job's temporary OPFS copy
		const mbwDir = new BlobDir({ [m.index + '.mbw']: m.indexFiles.mbw });
		preopens['/index'] = { open: (name, oflags, writable) => (name.endsWith('.l2b') ? opfs : mbwDir).open(name, oflags, writable) };
	}
	const idxPath = (fromFile ? '/index/' : '/opfs/') + m.index;
	const log = { stdout: (text) => postMessage({ type: 'log', text, mem: mem() }), stderr: (text) => postMessage({ type: 'log', text, mem: mem() }) };
	let wasi, mod = null, M = null;
	if (mt) {
		importScripts(m.wasmUrl);
		wasi = new Wasi({ args: [], mem64: true, preopens, ...log });
	} else {
		mod = await WebAssembly.compileStreaming(fetch(m.wasmUrl || 'minibwa.wasm'));
		wasi = new Wasi({ args: [], mem64: Wasi.isMem64(mod), preopens, ...log });
	}
	try {
		const t0 = performance.now();
		if (!fromFile) {
			await open(top, m.index + '.l2b', false);
			await open(top, m.index + '.mbw', false);
		} else if (!mt) {
			const tc = performance.now(), fh = await jobs.getFileHandle('index.l2b', { create: true });
			await m.indexFiles.l2b.stream().pipeTo(await fh.createWritable());
			handles[m.index + '.l2b'] = await fh.createSyncAccessHandle();
			postMessage({ type: 'log', text: `[align] ${m.index}.l2b → OPFS temp file (${(m.indexFiles.l2b.size / 1e9).toFixed(2)} GB, ${((performance.now() - tc) / 1000).toFixed(1)} s; for random access to the reference)\n` });
		}
		if (mt) {
			M = await self.createMinibwa({ mbwIO: wasi, mbwThreads: m.threads, mainScriptUrlOrBlob: m.wasmUrl });
			inst = { exports: { memory: M.getMemory() } };
			for (const f of ['mbw_alloc', 'mbw_free', 'mbw_memory', 'mbw_align_open', 'mbw_align_reads', 'mbw_align_header', 'mbw_align_batch', 'mbw_align_close']) inst.exports[f] = M['_' + f];
		} else inst = await WebAssembly.instantiate(mod, wasi.getImports());
		wasi.initialize(inst);
		const argv = [...(mt ? ['-P'] : ['-C', String(m.cacheMb)]), 'map', ...(mt ? ['-t', String(m.threads)] : []), ...m.args, idxPath, ...m.reads.map((r) => '/input/' + r)];
		postMessage({ type: 'log', text: `$ minibwa ${argv.join(' ')}  [batch=${m.batch} records/file, ${mt ? `wasm64, ${m.threads} threads` : wasi.m64 ? 'wasm64' : 'wasm32'}]\n` });
		if (wasi.call(inst, 'mbw_align_open', argv) !== 0) throw new Error('failed to load the index / parse options (see log)');
		if (!fromFile) close(m.index + '.mbw'); // BWT+SA are in memory now; .l2b stays open for the page cache
		const tLoad = (performance.now() - t0) / 1000;
		postMessage({ type: 'loaded', sec: tLoad, mem: mem() });

		const sets = [{ job: m.job, reads: m.reads }, ...(m.more || [])];
		let k = 0, records = 0, peak = mem(), tMap = 0, totRecords = 0, totBatches = 0;
		for (let si = 0; si < sets.length; ++si) {
			const set = sets[si];
			if (si > 0) {
				jobs = await jobsTop.getDirectoryHandle(set.job, { create: true });
				if (wasi.call(inst, 'mbw_align_reads', set.reads.map((r) => '/input/' + r)) !== 0) throw new Error('mbw_align_reads failed');
				postMessage({ type: 'log', text: `[align] ${set.job}: ${set.reads.join(' + ')} (same index in memory)\n` });
			}
			await open(jobs, 'header', true);
			wasi.callStr(inst, 'mbw_align_header', '/opfs/header');
			close('header');
			const totalIn = set.reads.reduce((a, r) => a + m.inputs[r].size, 0);
			const tSet = performance.now(), tMap0 = tMap;
			k = 0; records = 0;
			for (;;) {
				if (cancel) throw new Error('cancelled');
				const name = 'chunk.' + pad(k);
				await open(jobs, name, true); // async: lets 'cancel' messages in between batches
				const tb = performance.now();
				const n = wasi.callStr(inst, 'mbw_align_batch', '/opfs/' + name, m.batch);
				tMap += (performance.now() - tb) / 1000;
				close(name);
				if (n < 0) throw new Error('mbw_align_batch failed: ' + n);
				if (n === 0) { await jobs.removeEntry(name); break; }
				records += n; ++k; peak = Math.max(peak, mem());
				const readIn = set.reads.reduce((a, r) => a + (blobs.progress[r] || 0), 0);
				postMessage({ type: 'progress', job: set.job, set: si, sets: sets.length, batches: k, records, mapSec: tMap - tMap0, wallSec: (performance.now() - tSet) / 1000,
					inFrac: totalIn ? readIn / totalIn : 0, mem: mem(), peak, io: wasi.stats });
			}
			totRecords += records; totBatches += k;
			postMessage({ type: 'jobdone', job: set.job, batches: k, records, mapSec: tMap - tMap0, wallSec: (performance.now() - tSet) / 1000, peak });
		}
		k = totBatches; records = totRecords;
		inst.exports.mbw_align_close();
		return { type: 'done', batches: k, records, loadSec: tLoad, mapSec: tMap, wallSec: (performance.now() - t0) / 1000, peak, io: wasi.stats };
	} finally {
		for (const name of Object.keys(handles)) close(name);
	}
}

onmessage = async (e) => {
	if (e.data.type === 'cancel') { cancel = true; return; }
	if (e.data.type !== 'start') return;
	let res;
	try { res = await run(e.data); } catch (err) { res = { type: 'error', message: String((err && err.stack) || err) }; }
	postMessage(res);
};
