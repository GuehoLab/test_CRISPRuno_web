// [mt] --post-js: give the page access to the (shared) linear memory for the WASI shim
Module['getMemory'] = () => wasmMemory;
