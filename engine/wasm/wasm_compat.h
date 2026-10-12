/* Force-included (-include) when building for wasm32.
 * In wasm32 `long` is 32-bit, so fseek()/ftell() cannot address files >2GB.
 * Route them to the 64-bit off_t variants. */
#ifndef MB_WASM_COMPAT_H
#define MB_WASM_COMPAT_H
#include <stdio.h>
#include <sys/types.h>
#define fseek(fp, off, whence) fseeko((fp), (off_t)(off), (whence))
#define ftell(fp) ((long long)ftello(fp))
#ifdef __wasi__
#include <sys/mman.h>
#ifndef MAP_POPULATE
#define MAP_POPULATE 0
#endif
#endif
#endif
