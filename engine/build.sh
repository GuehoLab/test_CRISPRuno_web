#!/usr/bin/env bash
# Build minibwa as a WASI WebAssembly module (single-threaded, SIMD128) -> web/minibwa.wasm
#
# Toolchain: zig (clang 18 + wasi-libc). The easiest way to get it on any OS:
#     pip install ziglang==0.13.0
# Alternatively set CC to a wasi-sdk clang, e.g. CC="/opt/wasi-sdk/bin/clang --target=wasm32-wasi".
#
# Usage:
#   ./build.sh                          # minibwa sources from ../git/minibwa
#   MB_SRC=/path/to/minibwa ./build.sh
#   ZLIB_SRC=/path/to/zlib-1.3.1 ./build.sh   # otherwise zlib is downloaded from zlib.net
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
MB_SRC="${MB_SRC:-$HERE/../git/minibwa}"
BUILD="$HERE/build"
OUT="$HERE/../web"
CC="${CC:-python3 -m ziglang cc -target wasm32-wasi}"
ZLIB_VER=1.3.1

rm -rf "$BUILD/src" && mkdir -p "$BUILD/src" "$BUILD/zlib"
cp "$MB_SRC"/*.c "$MB_SRC"/*.h "$BUILD/src/"
( cd "$BUILD/src" && patch -p1 -s < "$HERE/wasm/minibwa-wasm.patch" && patch -p1 -s < "$HERE/wasm/end-bonus.patch" )
cp "$HERE"/wasm/*.c "$HERE"/wasm/*.h "$BUILD/src/"

# --- zlib ---
if [[ -z "${ZLIB_SRC:-}" ]]; then
  ZLIB_SRC="$BUILD/zlib-$ZLIB_VER"
  [[ -d "$ZLIB_SRC" ]] || ( cd "$BUILD" && curl -sSL "https://zlib.net/fossils/zlib-$ZLIB_VER.tar.gz" | tar xz )
fi
if [[ ! -f "$BUILD/zlib/libz.a" ]]; then
  for f in adler32 crc32 deflate inflate inffast inftrees trees zutil gzlib gzread gzwrite gzclose compress uncompr; do
    $CC -O3 -DZ_HAVE_UNISTD_H -c "$ZLIB_SRC/$f.c" -o "$BUILD/zlib/$f.o"
  done
  ar rcs "$BUILD/zlib/libz.a" "$BUILD"/zlib/*.o 2>/dev/null || python3 -m ziglang ar rcs "$BUILD/zlib/libz.a" "$BUILD"/zlib/*.o
fi

# --- minibwa ---
# SSE intrinsics are mapped to WASM SIMD128 by the emscripten compat headers in wasm/compat (MIT).
SSE="-D__SSE__ -D__SSE2__ -D__SSE3__ -D__SSSE3__ -D__SSE4_1__ -I$HERE/wasm/compat"
CFLAGS="-std=gnu99 -O3 -Wno-absolute-value -Wno-macro-redefined -msimd128 $SSE -DUSE_GPL -DHAVE_KALLOC -DMB_BATCH_API -DMBW_REACTOR -D_WASI_EMULATED_MMAN -D_WASI_EMULATED_PROCESS_CLOCKS
        -include wasm_compat.h -I$ZLIB_SRC -Wno-unused-result -Wno-deprecated-declarations"
SRCS="kommon.c kalloc.c bwt.c l2bit.c options.c seed.c map-algo.c lchain.c align.c pe.c cs.c format.c
      ksw2_extz2_sse.c ksw2_extd2_sse.c ksw2_ll_sse.c libsais.c libsais64.c index.c bseq.c map-main.c
      fastmap.c QSufSort.c bwtgen.c kthread_serial.c wasm_main.c mb_malloc.c"
cd "$BUILD/src"
OBJS=""
for f in $SRCS; do
  EXTRA=""
  [[ "$f" == ksw2_ll_sse.c ]] && EXTRA="-include emmintrin.h"   # it only includes <xmmintrin.h>
  $CC $CFLAGS $EXTRA -c "$f" -o "${f%.c}.o"
  OBJS="$OBJS ${f%.c}.o"
done
# 8 MiB stack. zig 0.13 cannot parse --max-memory >= 1 GiB, so the memory limits
# (64 MiB initial, 4 GiB = 65536 pages maximum, the wasm32 limit) are patched in afterwards.
# WASI "reactor": no _start; JS calls _initialize once, then the mbw_* exports (see wasm/wasm_main.c)
$CC -O3 -s -mexec-model=reactor $OBJS "$BUILD/zlib/libz.a" -o "$BUILD/minibwa.wasm" \
  -lwasi-emulated-mman -lwasi-emulated-process-clocks -Wl,-z,stack-size=8388608
python3 "$HERE/wasm/set_memory.py" "$BUILD/minibwa.wasm" "$OUT/minibwa.wasm" --initial-pages 1024 --max-pages 65536
ls -la "$OUT/minibwa.wasm"
