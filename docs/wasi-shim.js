/* Minimal WASI preview1 implementation for minibwa.wasm.
 *
 * Runs in a dedicated Web Worker (browser) or in Node.js (tests). Files are
 * served by "directory backends" mounted as WASI preopens:
 *   OpfsDir  - OPFS files through FileSystemSyncAccessHandle (read/write, synchronous)
 *   BlobDir  - user-selected File/Blob objects (read-only, FileReaderSync)
 *   NodeDir  - a real directory (Node.js only, for testing)
 * Reads go straight into WASM linear memory (no extra copy for OPFS).
 *
 * Two memory models:
 *   wasm32 (minibwa.wasm, wasi-libc): pointers are unsigned 32-bit numbers, so memory above 2 GB works
 *   wasm64 (minibwa64.wasm, emscripten -sMEMORY64 -sSTANDALONE_WASM; `mem64: true`): pointers and
 *          size_t are BigInt i64 at the JS boundary, iovecs are 16 bytes, and open() arrives as
 *          env.mbw_open(path, musl O_* flags) because emscripten's libc has no path_open.
 */
(function (root) {
'use strict';

const E = { SUCCESS: 0, ACCES: 2, BADF: 8, EXIST: 20, INVAL: 28, IO: 29, ISDIR: 31, NOENT: 44, NOSYS: 52, NOTSUP: 58, SPIPE: 70 };
const FT = { UNKNOWN: 0, CHAR: 2, DIR: 3, REG: 4 };
const O = { CREAT: 1, DIRECTORY: 2, EXCL: 4, TRUNC: 8 };
const FDFLAG_APPEND = 1;
const MO = { WRONLY: 1, RDWR: 2, CREAT: 0o100, EXCL: 0o200, TRUNC: 0o1000, APPEND: 0o2000 }; // musl open() flags (wasm64 build)

class WasiError extends Error { constructor(errno) { super('wasi errno ' + errno); this.errno = errno; } }
class WasiExit extends Error { constructor(code) { super('exit ' + code); this.code = code; } }

/******** backends ********/

class OpfsDir { // handles: { name: FileSystemSyncAccessHandle } opened by the caller (creation is async)
	constructor(handles) { this.handles = handles; }
	open(name, oflags) {
		const h = this.handles[name];
		if (!h) throw new WasiError((oflags & O.CREAT) ? E.ACCES : E.NOENT); // new files must be pre-opened by the worker
		if (oflags & O.TRUNC) h.truncate(0);
		return {
			read: (u8, pos) => h.read(u8, { at: pos }),
			write: (u8, pos) => h.write(u8, { at: pos }),
			size: () => h.getSize(),
			close: () => h.flush(),
		};
	}
}

class BlobDir { // files: { name: Blob }; progress[name] = furthest byte read
	constructor(files) { this.files = files; this.reader = new FileReaderSync(); this.progress = {}; }
	open(name, oflags) {
		const b = this.files[name];
		if (!b) throw new WasiError(E.NOENT);
		if (oflags & (O.CREAT | O.TRUNC)) throw new WasiError(E.ACCES);
		return {
			read: (u8, pos) => {
				this.progress[name] = Math.max(this.progress[name] || 0, Math.min(b.size, pos + u8.length));
				if (pos >= b.size) return 0;
				const ab = this.reader.readAsArrayBuffer(b.slice(pos, pos + u8.length));
				u8.set(new Uint8Array(ab));
				return ab.byteLength;
			},
			write: () => { throw new WasiError(E.ACCES); },
			size: () => b.size,
			close: () => {},
		};
	}
}

class NodeDir { // Node.js only
	constructor(dir) { this.dir = dir; this.fs = require('fs'); this.path = require('path'); }
	open(name, oflags, writable) {
		const fs = this.fs, p = this.path.join(this.dir, name);
		let flags;
		if (oflags & O.CREAT) flags = (oflags & O.TRUNC) ? 'w+' : (fs.existsSync(p) ? 'r+' : 'w+');
		else flags = writable ? 'r+' : 'r';
		let fd;
		try { fd = fs.openSync(p, flags); } catch (e) { throw new WasiError(e.code === 'ENOENT' ? E.NOENT : E.IO); }
		return {
			read: (u8, pos) => fs.readSync(fd, u8, 0, u8.length, pos),
			write: (u8, pos) => fs.writeSync(fd, u8, 0, u8.length, pos),
			size: () => fs.fstatSync(fd).size,
			close: () => fs.closeSync(fd),
		};
	}
}

/******** WASI ********/

class Wasi {
	constructor({ args = [], env = [], preopens = {}, stdout = null, stderr = null, mem64 = false }) {
		this.args = args; this.env = env; this.m64 = mem64;
		this.ps = mem64 ? 8 : 4; // sizeof(void*) == sizeof(size_t)
		this.stdout = stdout || ((s) => (typeof process !== 'undefined' ? process.stdout.write(s) : console.log(s)));
		this.stderr = stderr || ((s) => (typeof process !== 'undefined' ? process.stderr.write(s) : console.error(s)));
		this.fds = [{ type: 'stdin' }, { type: 'out', fn: (s) => this.stdout(s) }, { type: 'out', fn: (s) => this.stderr(s) }];
		for (const [path, dir] of Object.entries(preopens)) this.fds.push({ type: 'dir', path, dir });
		this.dec = new TextDecoder();
		this.memory = null;
		this.stats = { reads: 0, readBytes: 0, readMs: 0, writes: 0, writeBytes: 0, writeMs: 0 };
	}
	static isMem64(mod) { return WebAssembly.Module.imports(mod).some((i) => i.module === 'env' && i.name === 'mbw_open'); } // minibwa64.wasm
	static memory64Supported() { // validate a module with a single i64 memory
		try { return WebAssembly.validate(new Uint8Array([0, 0x61, 0x73, 0x6d, 1, 0, 0, 0, 5, 3, 1, 4, 1])); } catch (_) { return false; }
	}
	get dv() { return new DataView(this.memory.buffer); }
	P(x) { return typeof x === 'bigint' ? Number(x) : x < 0 ? x >>> 0 : x; } // pointer/size_t from WASM (i32 may be negative, i64 is BigInt) -> Number; never truncates numbers >= 2^32
	toPtr(n) { return this.m64 ? BigInt(n) : n; }             // Number -> pointer argument of an export
	setSize(ptr, v) { if (this.m64) this.dv.setBigUint64(this.P(ptr), BigInt(v), true); else this.dv.setUint32(this.P(ptr), v, true); }
	u8(ptr, len) { return new Uint8Array(this.memory.buffer, this.P(ptr), this.P(len)); }
	iovs(iovs, n) {
		const dv = this.dv, out = [], ps = this.ps;
		iovs = this.P(iovs); n = this.P(n);
		for (let i = 0; i < n; ++i) {
			const a = iovs + i * 2 * ps;
			out.push(this.m64 ? [Number(dv.getBigUint64(a, true)), Number(dv.getBigUint64(a + 8, true))] : [dv.getUint32(a, true), dv.getUint32(a + 4, true)]);
		}
		return out;
	}
	cstrAt(ptr) { // NUL-terminated string (paths: short); never view the whole >4 GB memory at once
		const p = this.P(ptr), u = new Uint8Array(this.memory.buffer, p, Math.min(65536, this.memory.buffer.byteLength - p));
		const e = u.indexOf(0);
		return this.dec.decode(u.slice(0, e < 0 ? u.length : e));
	}
	resolve(path) { // path -> [preopen fd entry, name inside it]
		let best = null;
		for (const f of this.fds) {
			if (!f || f.type !== 'dir') continue;
			const pre = f.path.replace(/\/+$/, '');
			if (pre === '.' && !path.startsWith('/')) { if (!best) best = [f, path.replace(/^\.\//, '')]; continue; }
			if (path.startsWith(pre + '/') && (!best || pre.length > best[0].path.length)) best = [f, path.slice(pre.length + 1)];
		}
		if (!best) throw new WasiError(E.NOENT);
		return best;
	}
	fd(fd) { const f = this.fds[fd]; if (!f) throw new WasiError(E.BADF); return f; }
	wrap(fn) { // convert WasiError exceptions to errno return values
		return (...a) => {
			try { return fn(...a) || 0; } catch (e) {
				if (e instanceof WasiError) return e.errno;
				throw e;
			}
		};
	}
	start(instance) { // command-style: run this.args once
		this.memory = instance.exports.memory;
		try {
			if (instance.exports._start) { instance.exports._start(); return 0; }
			this.initialize(instance);
			return this.call(instance, 'mbw_run', this.args.slice(1));
		} catch (e) {
			if (e instanceof WasiExit) return e.code;
			throw e;
		}
	}
	initialize(instance) { // reactor modules: once, before calling any export
		this.memory = instance.exports.memory;
		if (instance.exports._initialize) instance.exports._initialize();
	}
	cstr(instance, strs) { // NUL-separated strings in WASM memory (free with mbw_free)
		const enc = new TextEncoder(), parts = strs.map((x) => enc.encode(String(x) + '\0'));
		const n = parts.reduce((a, b) => a + b.length, 0), p = this.P(instance.exports.mbw_alloc(this.toPtr(n)));
		let o = 0;
		for (const b of parts) { this.u8(p + o, b.length).set(b); o += b.length; }
		return p;
	}
	call(instance, fn, strs) { // fn(argc, buf) with NUL-separated strs
		const p = this.cstr(instance, strs);
		try { return instance.exports[fn](strs.length, this.toPtr(p)); } finally { instance.exports.mbw_free(this.toPtr(p)); }
	}
	callStr(instance, fn, str, ...pre) { // fn(...pre, cstring)
		const p = this.cstr(instance, [str]);
		try { return instance.exports[fn](...pre, this.toPtr(p)); } finally { instance.exports.mbw_free(this.toPtr(p)); }
	}
	getImports() {
		const self = this, w = (f) => this.wrap(f);
		const enc = new TextEncoder();
		const strs = (list) => list.map((s) => enc.encode(s + '\0'));
		return { wasi_snapshot_preview1: {
			args_sizes_get: w((pc, ps) => { const a = strs(self.args); self.setSize(pc, a.length); self.setSize(ps, a.reduce((x, s) => x + s.length, 0)); }),
			args_get: w((pv, pb) => { pv = self.P(pv); pb = self.P(pb); for (const s of strs(self.args)) { self.setSize(pv, pb); self.u8(pb, s.length).set(s); pv += self.ps; pb += s.length; } }),
			environ_sizes_get: w((pc, ps) => { const a = strs(self.env); self.setSize(pc, a.length); self.setSize(ps, a.reduce((x, s) => x + s.length, 0)); }),
			environ_get: w((pv, pb) => { pv = self.P(pv); pb = self.P(pb); for (const s of strs(self.env)) { self.setSize(pv, pb); self.u8(pb, s.length).set(s); pv += self.ps; pb += s.length; } }),
			clock_time_get: w((id, prec, ptr) => {
				const ns = id === 0 ? BigInt(Date.now()) * 1000000n : BigInt(Math.round(performance.now() * 1e6));
				self.dv.setBigUint64(self.P(ptr), ns, true);
			}),
			clock_res_get: w((id, ptr) => { self.dv.setBigUint64(self.P(ptr), 1000n, true); }),
			random_get: w((ptr, len) => { const u = self.u8(ptr, len); for (let i = 0; i < u.length; ++i) u[i] = (Math.random() * 256) | 0; }),
			sched_yield: w(() => 0),
			proc_exit: (code) => { throw new WasiExit(code); },
			fd_prestat_get: w((fd, ptr) => {
				const f = self.fds[fd];
				if (!f || f.type !== 'dir') return E.BADF;
				self.dv.setUint8(self.P(ptr), 0);
				self.setSize(self.P(ptr) + self.ps, enc.encode(f.path).length);
			}),
			fd_prestat_dir_name: w((fd, ptr, len) => {
				const f = self.fds[fd];
				if (!f || f.type !== 'dir') return E.BADF;
				self.u8(ptr, len).set(enc.encode(f.path).subarray(0, len));
			}),
			fd_fdstat_get: w((fd, ptr) => {
				const f = self.fd(fd), dv = self.dv; ptr = self.P(ptr);
				dv.setUint8(ptr, f.type === 'dir' ? FT.DIR : f.type === 'file' ? FT.REG : FT.CHAR);
				dv.setUint16(ptr + 2, f.append ? FDFLAG_APPEND : 0, true);
				dv.setBigUint64(ptr + 8, 0xffffffffffffffffn, true);
				dv.setBigUint64(ptr + 16, 0xffffffffffffffffn, true);
			}),
			fd_fdstat_set_flags: w((fd, flags) => { const f = self.fd(fd); f.append = !!(flags & FDFLAG_APPEND); }),
			fd_filestat_get: w((fd, ptr) => {
				const f = self.fd(fd), dv = self.dv; ptr = self.P(ptr);
				for (let i = 0; i < 64; i += 8) dv.setBigUint64(ptr + i, 0n, true);
				dv.setUint8(ptr + 16, f.type === 'dir' ? FT.DIR : f.type === 'file' ? FT.REG : FT.CHAR);
				dv.setBigUint64(ptr + 24, 1n, true);
				if (f.type === 'file') dv.setBigUint64(ptr + 32, BigInt(f.file.size()), true);
			}),
			path_open: w((dirfd, dirflags, pptr, plen, oflags, rb, ri, fdflags, fdptr) => {
				const d = self.fd(dirfd);
				if (d.type !== 'dir') return E.BADF;
				let name = self.dec.decode(self.u8(pptr, plen).slice()).replace(/^\.\//, '');
				if (name === '' || name === '.' || (oflags & O.DIRECTORY)) return E.ISDIR;
				const writable = (BigInt(rb) & (1n << 6n)) !== 0n; // RIGHTS_FD_WRITE
				const file = d.dir.open(name, oflags, writable);
				self.fds.push({ type: 'file', name, file, pos: 0, append: !!(fdflags & FDFLAG_APPEND) });
				self.dv.setUint32(self.P(fdptr), self.fds.length - 1, true);
			}),
			fd_close: w((fd) => {
				const f = self.fd(fd);
				if (f.type === 'file') f.file.close();
				if (fd > 2 && f.type !== 'dir') self.fds[fd] = undefined;
			}),
			fd_seek: w((fd, off, whence, ptr) => {
				const f = self.fd(fd);
				if (f.type !== 'file') return E.SPIPE;
				const o = Number(off);
				const np = whence === 0 ? o : whence === 1 ? f.pos + o : f.file.size() + o;
				if (np < 0) return E.INVAL;
				f.pos = np;
				self.dv.setBigUint64(self.P(ptr), BigInt(np), true);
			}),
			fd_tell: w((fd, ptr) => {
				const f = self.fd(fd);
				if (f.type !== 'file') return E.SPIPE;
				self.dv.setBigUint64(self.P(ptr), BigInt(f.pos), true);
			}),
			fd_read: w((fd, iovs, n, ptr) => {
				const f = self.fd(fd);
				let tot = 0;
				if (f.type === 'file') {
					const t0 = performance.now();
					for (const [b, l] of self.iovs(iovs, n)) {
						if (l === 0) continue;
						const r = f.file.read(self.u8(b, l), f.pos);
						f.pos += r; tot += r;
						if (r < l) break;
					}
					self.stats.reads++; self.stats.readBytes += tot; self.stats.readMs += performance.now() - t0;
				} else if (f.type !== 'stdin') return E.BADF;
				self.setSize(ptr, tot);
			}),
			fd_pread: w((fd, iovs, n, off, ptr) => {
				const f = self.fd(fd);
				if (f.type !== 'file') return E.SPIPE;
				let tot = 0, pos = Number(off);
				const t0 = performance.now();
				for (const [b, l] of self.iovs(iovs, n)) {
					const r = f.file.read(self.u8(b, l), pos);
					pos += r; tot += r;
					if (r < l) break;
				}
				self.stats.reads++; self.stats.readBytes += tot; self.stats.readMs += performance.now() - t0;
				self.setSize(ptr, tot);
			}),
			fd_write: w((fd, iovs, n, ptr) => {
				const f = self.fd(fd);
				let tot = 0;
				if (f.type === 'out') {
					let s = '';
					for (const [b, l] of self.iovs(iovs, n)) { s += self.dec.decode(self.u8(b, l).slice()); tot += l; }
					f.fn(s);
				} else if (f.type === 'file') {
					const t0 = performance.now();
					if (f.append) f.pos = f.file.size();
					for (const [b, l] of self.iovs(iovs, n)) {
						if (l === 0) continue;
						const r = f.file.write(self.u8(b, l), f.pos);
						f.pos += r; tot += r;
					}
					self.stats.writes++; self.stats.writeBytes += tot; self.stats.writeMs += performance.now() - t0;
				} else return E.BADF;
				self.setSize(ptr, tot);
			}),
		}, env: { // wasm64 build only (see wasm/em64_sys.c)
			mbw_open: (pptr, flags) => { // returns a new fd or -errno
				try {
					const [d, name] = self.resolve(self.cstrAt(pptr));
					if (name === '' || name === '.') return -E.ISDIR;
					const of = ((flags & MO.CREAT) ? O.CREAT : 0) | ((flags & MO.EXCL) ? O.EXCL : 0) | ((flags & MO.TRUNC) ? O.TRUNC : 0);
					const file = d.dir.open(name, of, (flags & 3) !== 0);
					self.fds.push({ type: 'file', name, file, pos: 0, append: !!(flags & MO.APPEND) });
					return self.fds.length - 1;
				} catch (e) {
					if (e instanceof WasiError) return -e.errno;
					throw e;
				}
			},
			emscripten_notify_memory_growth: () => {},
		} };
	}
}

const api = { Wasi, WasiError, WasiExit, OpfsDir, BlobDir, NodeDir, E };
if (typeof module !== 'undefined' && module.exports) module.exports = api;
root.MinibwaWasi = api;
})(typeof self !== 'undefined' ? self : globalThis);
