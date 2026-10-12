#!/usr/bin/env bash
# Build minibwa as a wasm64 (Memory64) module -> web/minibwa64.wasm
#
# Same sources and patch as build.sh (wasm32), but compiled with emscripten in
# standalone mode with -sMEMORY64=1: pointers are 64-bit and linear memory can
# grow past 4 GiB (up to 16 GiB here, the current browser limit for memory64).
# The module only imports WASI preview1 functions plus env.mbw_open (see
# wasm/em64_sys.c), so it runs on the same JS shim (web/wasi-shim.js, mem64 mode).
#
# Toolchain: emscripten >= 3.1.60 (tested: 3.1.69 with LLVM 19, binaryen 120).
#   EMCC=emcc ./build64.sh          (default: emcc on PATH)
#   MB_SRC=/path/to/minibwa ZLIB_SRC=/path/to/zlib-1.3.1 ./build64.sh
# Browsers: Chrome/Edge >= 133, Firefox >= 134 (WebAssembly Memory64). Not Safari (2026-10).
#
# MT=1 ./build64.sh -> web/minibwa64-mt.js + web/minibwa64-mt.wasm: the same code with pthreads
# (shared memory; needs a cross-origin isolated page, i.e. COOP/COEP headers). It is not standalone:
# emscripten's JS glue creates the thread workers; file I/O still goes to web/wasi-shim.js through
# wasm/mbw_io.js. kt_for() is multi-threaded (index build: bwtgen + SA sampling), kt_pipeline() and
# all I/O stay on the main thread (wasm/kthread_mt.c).
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
MB_SRC="${MB_SRC:-$HERE/../git/minibwa}"
MT="${MT:-0}"
BUILD="$HERE/build64"; [[ "$MT" == 1 ]] && BUILD="$HERE/build64mt"
OUT="$HERE/../web"
EMCC="${EMCC:-emcc}"
EMAR="${EMAR:-${EMCC%cc*}ar}"
ZLIB_VER=1.3.1

rm -rf "$BUILD/src" && mkdir -p "$BUILD/src" "$BUILD/zlib"
cp "$MB_SRC"/*.c "$MB_SRC"/*.h "$BUILD/src/"
( cd "$BUILD/src" && patch -p1 -s < "$HERE/wasm/minibwa-wasm.patch" && patch -p1 -s < "$HERE/wasm/end-bonus.patch" )
cp "$HERE"/wasm/wasm_main.c "$HERE"/wasm/kthread_serial.c "$HERE"/wasm/kthread_mt.c "$HERE"/wasm/em64_sys.c "$BUILD/src/"

M64="-sMEMORY64=1 -Wno-experimental"
[[ "$MT" == 1 ]] && M64="$M64 -pthread"
# --- zlib ---
if [[ -z "${ZLIB_SRC:-}" ]]; then
  ZLIB_SRC="$BUILD/zlib-$ZLIB_VER"
  [[ -d "$ZLIB_SRC" ]] || ( cd "$BUILD" && curl -sSL "https://zlib.net/fossils/zlib-$ZLIB_VER.tar.gz" | tar xz )
fi
if [[ ! -f "$BUILD/zlib/libz.a" ]]; then
  for f in adler32 crc32 deflate inflate inffast inftrees trees zutil gzlib gzread gzwrite gzclose compress uncompr; do
    $EMCC $M64 -O3 -DZ_HAVE_UNISTD_H -c "$ZLIB_SRC/$f.c" -o "$BUILD/zlib/$f.o"
  done
  $EMAR rcs "$BUILD/zlib/libz.a" "$BUILD"/zlib/*.o
fi

# --- minibwa ---
# emscripten maps SSE..SSE4.1 intrinsics to WASM SIMD128 itself (-msse4.1).
CFLAGS="$M64 -std=gnu99 -O3 -msimd128 -msse4.1 -Wno-absolute-value -Wno-macro-redefined -DUSE_GPL -DHAVE_KALLOC -DMB_BATCH_API -DMBW_REACTOR
        -I$ZLIB_SRC -Wno-unused-result -Wno-deprecated-declarations"
SRCS="kommon.c kalloc.c bwt.c l2bit.c options.c seed.c map-algo.c lchain.c align.c pe.c cs.c format.c
      ksw2_extz2_sse.c ksw2_extd2_sse.c ksw2_ll_sse.c libsais.c libsais64.c index.c bseq.c map-main.c
      fastmap.c QSufSort.c bwtgen.c wasm_main.c em64_sys.c"
[[ "$MT" == 1 ]] && SRCS="$SRCS kthread_mt.c" || SRCS="$SRCS kthread_serial.c"
cd "$BUILD/src"
OBJS=""
for f in $SRCS; do
  EXTRA=""
  [[ "$f" == ksw2_ll_sse.c ]] && EXTRA="-include emmintrin.h"   # it only includes <xmmintrin.h>
  $EMCC $CFLAGS $EXTRA -c "$f" -o "${f%.c}.o"
  OBJS="$OBJS ${f%.c}.o"
done
if [[ "$MT" == 1 ]]; then
  # JS glue (MODULARIZE: createMinibwa) + wasm. The thread pool (Module.mbwThreads workers) is started
  # before the module resolves, because a pthread cannot start while the main thread is blocked.
  EXPORTS=_mbw_alloc,_mbw_free,_mbw_memory,_mbw_run,_mbw_align_open,_mbw_align_reads,_mbw_align_header,_mbw_align_batch,_mbw_align_close
  $EMCC $M64 -O3 --no-entry -sMODULARIZE=1 -sEXPORT_NAME=createMinibwa -sENVIRONMENT=worker -sFILESYSTEM=0 \
    -sALLOW_MEMORY_GROWTH=1 -sINITIAL_MEMORY=64MB -sMAXIMUM_MEMORY=16GB -sSTACK_SIZE=8MB -sDEFAULT_PTHREAD_STACK_SIZE=2MB -sMEMORY_GROWTH_GEOMETRIC_STEP=0 \
    -sPTHREAD_POOL_SIZE='Module["mbwThreads"]' -sEXPORTED_FUNCTIONS=$EXPORTS -Wno-pthreads-mem-growth \
    --js-library "$HERE/wasm/mbw_io.js" --post-js "$HERE/wasm/mbw_post.js" $OBJS "$BUILD/zlib/libz.a" -o "$OUT/minibwa64-mt.js"
  # emscripten 3.1.69 still creates i64 memories with the old descriptor key "index"; browsers use "address"
  sed -i 's/index:"i64"/address:"i64"/g' "$OUT/minibwa64-mt.js"
  ls -la "$OUT/minibwa64-mt.js" "$OUT/minibwa64-mt.wasm"
  exit 0
fi
# Standalone reactor: no JS glue; JS calls _initialize once, then the mbw_* exports.
# 64 MiB initial, grows on demand up to 16 GiB; 8 MiB stack.
$EMCC $M64 -O3 -sSTANDALONE_WASM=1 --no-entry -sALLOW_MEMORY_GROWTH=1 -sINITIAL_MEMORY=64MB -sMAXIMUM_MEMORY=16GB \
  -sSTACK_SIZE=8MB $OBJS "$BUILD/zlib/libz.a" -o "$OUT/minibwa64.wasm"
ls -la "$OUT/minibwa64.wasm"
