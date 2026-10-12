/* System glue for the wasm64 (Memory64) build made with emscripten -sSTANDALONE_WASM.
 *
 * In standalone mode emscripten's libc does all file I/O through WASI imports
 * (fd_read, fd_write, fd_seek, fd_pread, fd_close), but open() ends in a weak
 * __syscall_openat() that only knows /dev/std*. We override it with a call to
 * the JS host (env.mbw_open), which resolves the path against the same OPFS /
 * Blob "preopens" as the wasm32 build and returns a WASI fd (or -errno).
 */
#include <errno.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <sys/resource.h>

__attribute__((import_module("env"), import_name("mbw_open")))
int mbw_js_open(const char *path, int flags);

int __syscall_openat(int dirfd, intptr_t path, int flags, ...) // mode (optional vararg) is ignored
{
	const char *p = (const char*)path;
	(void)dirfd;
	if (strcmp(p, "/dev/stdin") == 0) return 0;
	if (strcmp(p, "/dev/stdout") == 0) return 1;
	if (strcmp(p, "/dev/stderr") == 0) return 2;
	return mbw_js_open(p, flags);
}

int __syscall_getrusage(int who, intptr_t usage) // only the CPU time is used (kom_cputime)
{ // no CPU clock here: report the wall time since the first call (for %CPU in the logs, ~100% per busy thread is not visible)
	static int64_t t0 = -1;
	struct rusage *r = (struct rusage*)usage;
	struct timespec ts;
	int64_t us;
	(void)who;
	memset(r, 0, sizeof(*r));
	clock_gettime(CLOCK_MONOTONIC, &ts);
	us = (int64_t)ts.tv_sec * 1000000 + ts.tv_nsec / 1000;
	if (t0 < 0) t0 = us;
	us -= t0;
	r->ru_utime.tv_sec = us / 1000000;
	r->ru_utime.tv_usec = us % 1000000;
	return 0;
}
