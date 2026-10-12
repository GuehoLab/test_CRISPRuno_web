/* malloc for the wasm32 build: dlmalloc on top of a sbrk that can grow the
 * heap by more than 2 GB at once (wasi-libc's sbrk takes a signed increment,
 * and its default allocator cannot allocate past ~2 GB). */
#include <stddef.h>
#include <stdint.h>
#include <errno.h>

extern unsigned char __heap_base;
static uintptr_t mb_brk = 0;

void *mb_morecore(ptrdiff_t inc_)
{
	size_t inc = (size_t)inc_; // dlmalloc never trims (MORECORE_CANNOT_TRIM); sizes may exceed 2 GB
	uint64_t new_end, cur_end;
	uintptr_t old;
	if (mb_brk == 0) mb_brk = (uintptr_t)&__heap_base;
	old = mb_brk;
	if (inc == 0) return (void*)old;
	new_end = (uint64_t)old + inc;
	if (new_end > 0xffffffffULL) { errno = ENOMEM; return (void*)-1; }
	cur_end = (uint64_t)__builtin_wasm_memory_size(0) << 16;
	if (new_end > cur_end) {
		size_t pages = (size_t)((new_end - cur_end + 65535) >> 16);
		if (__builtin_wasm_memory_grow(0, pages) == (size_t)-1) { errno = ENOMEM; return (void*)-1; }
	}
	mb_brk = (uintptr_t)new_end;
	return (void*)old;
}

#define HAVE_MMAP 0
#define HAVE_MREMAP 0
#define MORECORE mb_morecore
#define MORECORE_CONTIGUOUS 1
#define MORECORE_CANNOT_TRIM 1
#define USE_LOCKS 0
#define NO_MALLOC_STATS 1
#define MALLOC_ALIGNMENT 16
#define DEFAULT_GRANULARITY ((size_t)65536U)
#define malloc_getpagesize ((size_t)65536U)
#define LACKS_SCHED_H 1
#define LACKS_SYS_MMAN_H 1
#define MB_MAX_SYS_REQ ((size_t)0xfc000000U) /* ~3.94 GB */
#include "dlmalloc.c"

/* public entry points (the emscripten copy of dlmalloc only aliases them under __EMSCRIPTEN__) */
void *malloc(size_t n) { return dlmalloc(n); }
void free(void *p) { dlfree(p); }
void *calloc(size_t n, size_t s) { return dlcalloc(n, s); }
void *realloc(void *p, size_t n) { return dlrealloc(p, n); }
void *memalign(size_t a, size_t n) { return dlmemalign(a, n); }
int posix_memalign(void **pp, size_t a, size_t n) { return dlposix_memalign(pp, a, n); }
void *aligned_alloc(size_t a, size_t n) { return dlmemalign(a, n); }
size_t malloc_usable_size(void *p) { return dlmalloc_usable_size(p); }
