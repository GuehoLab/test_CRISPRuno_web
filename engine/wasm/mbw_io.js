// [mt] I/O imports of the multi-threaded wasm64 build, routed to the WASI shim (web/wasi-shim.js)
// object that the main worker puts in Module.mbwIO. All file I/O happens on the main thread
// (kt_pipeline is serial, and the align worker loads the whole 2-bit reference with -P so that
// mapping threads never read files). If a worker thread ever writes to stdout/stderr (e.g. a
// warning), it has no shim object, so the text goes to the browser console instead.
addToLibrary({
  $mbwIo: () => {
    const w = Module['mbwIO'];
    if (!w) return null; // a pthread: no shim (and no files) here
    w.memory = wasmMemory;
    return w._imp || (w._imp = w.getImports());
  },
  $mbwThreadWrite: (fd, iov, n, pnum) => { // stdout/stderr from a worker thread -> console
    if (fd !== 1 && fd !== 2) return 8; // EBADF
    const dv = new DataView(wasmMemory.buffer), dec = new TextDecoder();
    let s = '', tot = 0;
    iov = Number(iov); n = Number(n);
    for (let i = 0; i < n; ++i) {
      const p = Number(dv.getBigUint64(iov + i * 16, true)), l = Number(dv.getBigUint64(iov + i * 16 + 8, true));
      s += dec.decode(new Uint8Array(wasmMemory.buffer, p, l).slice());
      tot += l;
    }
    console.error('[minibwa thread] ' + s);
    dv.setBigUint64(Number(pnum), BigInt(tot), true);
    return 0;
  },
  fd_write__deps: ['$mbwIo', '$mbwThreadWrite'],
  fd_write: (fd, iov, n, pnum) => { const io = mbwIo(); return io ? io.wasi_snapshot_preview1.fd_write(fd, iov, n, pnum) : mbwThreadWrite(fd, iov, n, pnum); },
  fd_read__deps: ['$mbwIo'], fd_read: (fd, iov, n, pnum) => { const io = mbwIo(); return io ? io.wasi_snapshot_preview1.fd_read(fd, iov, n, pnum) : 8; },
  fd_pread__deps: ['$mbwIo'], fd_pread: (fd, iov, n, off, pnum) => { const io = mbwIo(); return io ? io.wasi_snapshot_preview1.fd_pread(fd, iov, n, off, pnum) : 8; },
  fd_seek__deps: ['$mbwIo'], fd_seek: (fd, off, wh, pnew) => { const io = mbwIo(); return io ? io.wasi_snapshot_preview1.fd_seek(fd, off, wh, pnew) : 8; },
  fd_close__deps: ['$mbwIo'], fd_close: (fd) => { const io = mbwIo(); return io ? io.wasi_snapshot_preview1.fd_close(fd) : 8; },
  mbw_open__deps: ['$mbwIo'], mbw_open__sig: 'ipi', mbw_open: (p, f) => { const io = mbwIo(); return io ? io.env.mbw_open(p, f) : -8; },
});
