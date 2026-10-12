// SPDX-License-Identifier: LicenseRef-CRISPRuno — CRISPRuno web port: step runner (OPFS)
/* Runs one pipeline step (cu/pipeline.js) on the files of a run folder in OPFS: minibwa/crispruno/<run>/.
 * in : { type:'run', step, run, cfg }   cfg: { settings, files:{name: File}, reference:{fasta, fai} | {l2b} , ... step specific }
 * out: { type:'log'|'progress'|'done'|'error', ... }
 */
importScripts('lib/fflate.min.js', 'lib/pyfmt.js', 'lib/textio.js', 'cu/core.js', 'cu/cralign.js', 'cu/ssw.js', 'cu/cutadapt.js', 'cu/umi.js', 'cu/primer.js', 'cu/sam.js',
	'cu/assign.js', 'cu/crispresso.js', 'cu/pipeline.js', 'cu/plots.js', 'cu/report.js');
const IO = self.GsIO, PL = self.CuPipeline;
const log = (text, cls) => postMessage({ type: 'log', text, cls });
let lastProg = 0;
const progress = (p) => { const t = performance.now(); if (t - lastProg < 200 && p.frac !== 1) return; lastProg = t; postMessage({ type: 'progress', ...p }); };
const tick = () => new Promise((r) => setTimeout(r, 0));
async function syncHandle(fh) { await tick(); return fh.createSyncAccessHandle(); } // yield first (Chrome could stall right after other handles were closed)

async function topDir() { return (await navigator.storage.getDirectory()).getDirectoryHandle('minibwa', { create: true }); }
async function runDir(run) { return (await (await (await topDir()).getDirectoryHandle('crispruno', { create: true })).getDirectoryHandle(run, { create: true })); }

function makeIo(dir, cfg) {
	const handles = new Map(); // name -> SyncAccessHandle (kept open until the end of the step)
	const extra = []; // other handles (tmp files, index .l2b)
	const closeName = (name) => { const h = handles.get(name); if (h) { try { h.flush(); h.close(); } catch (_) {} handles.delete(name); } };
	const getH = async (name, create) => {
		if (handles.has(name)) return handles.get(name);
		let fh; try { fh = await dir.getFileHandle(name, { create }); } catch (_) { return null; }
		const h = await syncHandle(fh); handles.set(name, h); return h;
	};
	let genome = null;
	const io = {
		async writer(name) {
			closeName(name);
			const h = await getH(name, true); h.truncate(0);
			let at = 0;
			const sink = new IO.TextSink((b) => { h.write(b, { at }); at += b.length; });
			sink.name = name; sink.h = h;
			sink.raw = (u8) => { sink.flush(); h.write(u8, { at }); at += u8.length; };
			sink.close = () => { sink.flush(); h.flush(); };
			return sink;
		},
		linesOf: (w) => IO.lines(IO.handleSource(w.h)),
		async source(name) { const h = await getH(name, false); return h ? IO.handleSource(h) : null; },
		async reader(name) { const h = await getH(name, false); if (!h) throw new Error(name + ' not found'); return new IO.LineReader(IO.handleSource(h)); },
		async userReader(name) { const f = cfg.files && cfg.files[name]; if (!f) throw new Error('input file not selected: ' + name); return new IO.BlobLineReader(f); },
		async readText(name) {
			if (handles.has(name)) { const s = IO.handleSource(handles.get(name)); return new TextDecoder().decode(s.read(0, s.size)); }
			try { return await (await (await dir.getFileHandle(name)).getFile()).text(); } catch (_) { return null; }
		},
		async exists(name) { try { await dir.getFileHandle(name); return true; } catch (_) { return false; } },
		async remove(name) { closeName(name); try { await dir.removeEntry(name); } catch (_) {} },
		async getGenome() {
			if (genome) return genome;
			const R = cfg.reference || {};
			if (R.fasta) {
				const fai = R.fai ? await R.fai.text() : null, t0 = performance.now();
				genome = new IO.Fasta(IO.cachedSource(IO.blobSource(R.fasta)), fai);
				if (!fai) log(`[reference] ${R.fasta.name}: .fai index built in ${((performance.now() - t0) / 1000).toFixed(1)} s (${genome.idx.size} sequences)`);
			} else if (R.l2b) {
				let src;
				if (R.l2b instanceof Blob) src = IO.blobSource(R.l2b);
				else { const h = await syncHandle(await (await topDir()).getFileHandle(R.l2b)); extra.push(h); src = IO.handleSource(h); }
				log('[reference] sequences from the index .l2b (soft-masking kept, IUPAC codes read as N); select the FASTA for exact samtools faidx output', 'warn');
				genome = new IO.L2bGenome(IO.cachedSource(src, 256 << 10, 256));
			} else throw new Error('reference sequence (FASTA or index .l2b) not selected');
			return genome;
		},
		tmp: {
			async create(name) {
				const td = await dir.getDirectoryHandle('tmp', { create: true });
				const h = await syncHandle(await td.getFileHandle(name, { create: true })); extra.push(h); h.truncate(0);
				let at = 0;
				return { write: (u8) => { h.write(u8, { at }); at += u8.length; }, close: () => { h.flush(); return IO.handleSource(h); } };
			},
			async remove(name) {
				for (let i = extra.length - 1; i >= 0; --i) { /* closed with the others at the end */ }
			},
		},
		log, progress,
		closeAll() {
			for (const name of [...handles.keys()]) closeName(name);
			for (const h of extra) { try { h.close(); } catch (_) {} }
			extra.length = 0;
		},
	};
	return io;
}

async function readChunks(io, jobsDirName, job) { // header + chunk.* of an align-worker job as byte sources
	const jobs = await (await topDir()).getDirectoryHandle(jobsDirName);
	const jd = await jobs.getDirectoryHandle(job);
	const names = [];
	for await (const [n] of jd.entries()) if (n.startsWith('chunk.')) names.push(n);
	names.sort();
	const headerText = await (await (await jd.getFileHandle('header')).getFile()).text();
	const sources = [];
	for (const n of names) { const h = await syncHandle(await jd.getFileHandle(n)); io._chunkHandles = (io._chunkHandles || []).concat(h); sources.push(IO.handleSource(h)); }
	return { headerText, sources };
}

const STEPS = {
	prep: (io, S, cfg) => PL.stepPrep(io, S, { locateSam: cfg.locateSam }),
	umi: (io, S) => PL.stepUmi(io, S),
	primer: async (io, S) => PL.stepPrimer(io, S, { ssw: await CuSsw.load('cu/ssw.wasm') }),
	sort: async (io, S, cfg) => {
		const { headerText, sources } = await readChunks(io, cfg.jobsDir, cfg.job);
		try { return await PL.stepSortAlignment(io, S, { headerText, sources, alignerLog: cfg.alignerLog }); }
		finally {
			for (const h of io._chunkHandles || []) { try { h.close(); } catch (_) {} }
			try { await (await (await topDir()).getDirectoryHandle(cfg.jobsDir)).removeEntry(cfg.job, { recursive: true }); } catch (_) {}
		}
	},
	assign: (io, S) => PL.stepAssign(io, S),
	crispresso: async (io, S, cfg) => PL.stepCrispresso(io, S, { ednafull: await (await fetch('cu/EDNAFULL')).text(), fixSpaceIds: !!cfg.fixSpaceIds }),
	summary: async (io, S) => {
		const r = await PL.stepSummary(io, S);
		let chromLens = [];
		try { const g = await io.getGenome(); chromLens = g.contigs ? g.contigs() : [...g.idx].map(([n, x]) => [n, x.len]); } catch (_) {}
		const rep = await CuReport.buildReport(io, S, { chromLens });
		log(`[report] ${S.root}.html (${rep.plots} plots)`, 'ok');
		return r;
	},
};

onmessage = async (e) => {
	const m = e.data;
	if (m.type !== 'run') return;
	let io = null;
	try {
		const dir = await runDir(m.run);
		io = makeIo(dir, m.cfg);
		const t0 = performance.now();
		const r = await STEPS[m.step](io, m.cfg.settings, m.cfg);
		io.closeAll();
		try { const td = await dir.getDirectoryHandle('tmp'); await dir.removeEntry('tmp', { recursive: true }); } catch (_) {}
		postMessage({ type: 'done', step: m.step, result: r, sec: (performance.now() - t0) / 1000 });
	} catch (err) {
		if (io) io.closeAll();
		postMessage({ type: 'error', step: m.step, message: String((err && err.message) || err), stack: String(err && err.stack) });
	}
};
