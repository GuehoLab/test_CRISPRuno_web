#!/usr/bin/env bash
# Build the Striped Smith-Waterman library (ssw.c, as shipped with CRISPRuno lib/ssw; MIT) as a WASI reactor module
# -> ../../web/cu/ssw.wasm. SSE2 intrinsics map to WASM SIMD128 through the emscripten compat headers in compat/ (MIT).
# Toolchain: zig 0.13 (pip install ziglang==0.13.0) or CC="clang --target=wasm32-wasi" from wasi-sdk.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
CC="${CC:-python3 -m ziglang cc -target wasm32-wasi}"
$CC -O3 -msimd128 -D__SSE__ -D__SSE2__ -I"$HERE/compat" -mexec-model=reactor \
  -Wl,--export=ssw_init -Wl,--export=ssw_align -Wl,--export=align_destroy -Wl,--export=init_destroy \
  -Wl,--export=malloc -Wl,--export=free -Wl,--export-memory -Wl,--strip-all \
  -o "$HERE/../../web/cu/ssw.wasm" "$HERE/ssw.c"
ls -l "$HERE/../../web/cu/ssw.wasm"
