/* Entry point for the WebAssembly build of minibwa.
 *
 * Each subcommand is meant to run in its own Web Worker (one fresh WASM
 * instance per stage). WASM linear memory never shrinks, so terminating the
 * worker after each stage is the only way to give memory back. Files are
 * handed from one stage to the next through OPFS (mounted at /opfs).
 *
 * Index pipeline for large genomes (e.g. T2T-CHM13, 3.1 Gb):
 *   prep     FASTA -> .l2b + .pac (both strands)        ~1.5 GB peak
 *   genraw   .pac  -> .raw  (BWT-SW / bwtgen, GPL)       ~2.3 GB peak (100m block)
 *   raw2mbw  .raw  -> .mbw  (stream-encode + sample SA)  ~3.6 GB peak (-u 7)
 * Small genomes can use `index` (libsais, all in memory, ~10 bytes/bp).
 *
 * `map` loads .mbw into memory but keeps the 2-bit reference (.l2b pac[]) on
 * disk/OPFS and reads it on demand, saving ~N/4 bytes of heap.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "kommon.h"
#include "mbpriv.h"
#include "ketopt.h"

int main_index(int argc, char *argv[]);
int main_map(int argc, char *argv[]);
int main_mem(int argc, char *argv[]);
int main_genraw(int argc, char *argv[]);
int main_fastmap(int argc, char *argv[]);

#ifdef __wasi__
#include <sys/mman.h>
int madvise(void *addr, size_t len, int advice) { (void)addr, (void)len, (void)advice; return 0; } // used by --mmap; not provided by wasi-libc
#endif

static void report_heap(const char *what)
{
#ifdef __wasm__
	fprintf(stderr, "[wasm] %s: linear memory %.3f GB\n", what, (double)__builtin_wasm_memory_size(0) * 65536.0 / 1e9);
#else
	(void)what;
#endif
}

/******************************************
 * prep: FASTA -> .l2b + .pac (streaming) *
 ******************************************/

typedef struct {
	FILE *fp;
	uint8_t *buf;
	uint64_t n_buf, x; // x: number of bases written so far
	uint8_t cur;
} pac_writer_t;

static void pw_put(pac_writer_t *w, int c)
{
	w->cur |= (uint8_t)(c & 3) << ((~w->x & 3) * 2);
	if ((++w->x & 3) == 0) {
		w->buf[w->n_buf++] = w->cur;
		w->cur = 0;
		if (w->n_buf == 1<<20) {
			fwrite(w->buf, 1, w->n_buf, w->fp);
			w->n_buf = 0;
		}
	}
}

static int write_pac_stream(const char *fn, const l2b_t *l2b) // identical output to l2b_save_pac(fn, l2b, 1)
{
	pac_writer_t w;
	int64_t i;
	uint8_t ct;
	memset(&w, 0, sizeof(w));
	w.fp = fopen(fn, "wb");
	if (w.fp == 0) return -1;
	w.buf = (uint8_t*)malloc(1<<20);
	for (i = 0; i < (int64_t)l2b->tot_len; ++i) pw_put(&w, l2b_get0(l2b, i));
	for (i = l2b->tot_len - 1; i >= 0; --i) pw_put(&w, 3 - l2b_get0(l2b, i));
	if (w.x & 3) w.buf[w.n_buf++] = w.cur;
	fwrite(w.buf, 1, w.n_buf, w.fp);
	if (w.x % 4 == 0) {
		ct = 0;
		fwrite(&ct, 1, 1, w.fp);
	}
	ct = w.x % 4;
	fwrite(&ct, 1, 1, w.fp);
	fclose(w.fp);
	free(w.buf);
	return 0;
}

static int main_prep(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	uint64_t seed = 11;
	int c;
	l2b_t *l2b;
	char *fn;
	while ((c = ketopt(&o, argc, argv, 1, "s:", 0)) >= 0)
		if (c == 's') seed = atol(o.arg);
	if (argc - o.ind < 2) {
		fprintf(stderr, "Usage: minibwa prep [-s seed] <in.fa> <out.prefix>\n");
		return 1;
	}
	l2b = l2b_import(argv[o.ind], seed);
	if (l2b == 0) {
		fprintf(stderr, "[E::%s] failed to read the FASTA\n", __func__);
		return 1;
	}
	fprintf(stderr, "[M::%s] %lld contigs, %lld bp\n", __func__, (long long)l2b->n_ctg, (long long)l2b->tot_len);
	fn = (char*)calloc(strlen(argv[o.ind+1]) + 10, 1);
	strcat(strcpy(fn, argv[o.ind+1]), ".l2b");
	if (l2b_save(fn, l2b) != 0) return 1;
	strcat(strcpy(fn, argv[o.ind+1]), ".pac");
	if (write_pac_stream(fn, l2b) != 0) return 1;
	free(fn);
	report_heap("prep");
	l2b_destroy(l2b);
	return 0;
}

/****************************************************
 * raw2mbw: stream raw bwtgen BWT -> encoded .mbw   *
 ****************************************************/

#define BWT_CNT_SHIFT 56

static int main_raw2mbw(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	int c, sa_bit = 4;
	const uint64_t chunk_words = 1<<22; // 16 MB of raw input per read
	uint64_t L2[5], primary, len, i, k, c4[4], x[4], *last_c = 0, w_base = 0, w_n = 0;
	uint32_t *chunk;
	mb_bwt_t *bwt;
	FILE *fp;
	double t;

	while ((c = ketopt(&o, argc, argv, 1, "u:t:", 0)) >= 0) {
		if (c == 'u') sa_bit = atoi(o.arg);
		else if (c == 't') mb_n_threads = atoi(o.arg); // [par] SA sampling threads (multi-threaded build only)
	}
	if (argc - o.ind < 2) {
		fprintf(stderr, "Usage: minibwa raw2mbw [-u sa_bit] [-t threads] <in.raw> <out.mbw>\n");
		return 1;
	}
	fp = fopen(argv[o.ind], "rb");
	if (fp == 0) return 1;
	if (fread(&primary, 8, 1, fp) != 1 || fread(&L2[1], 8, 4, fp) != 4) return 1;
	L2[0] = 0, len = L2[4];

	// the same encoding as mb_bwt_init_from_raw() but reading the raw BWT chunk by chunk
	bwt = mb_bwt_init();
	bwt->primary = primary;
	bwt->seq_len = len;
	bwt->data_len = ((len + 127) / 128) * 4 + ((len + 127) / 128 + 1) * 4;
	bwt->data = kom_calloc(uint64_t, bwt->data_len);
	kom_assert(bwt->data, "failed to allocate the BWT");
	chunk = (uint32_t*)malloc(chunk_words * 4);
	memset(c4, 0, 32);
	t = kom_realtime();
	for (i = k = 0; i < len; ++i) {
		uint64_t wi = i >> 4;
		uint8_t j, a;
		if (wi >= w_base + w_n) { // load the next chunk
			w_base = wi;
			w_n = fread(chunk, 4, chunk_words, fp);
			kom_assert(w_n > 0, "truncated raw BWT");
		}
		a = chunk[wi - w_base] >> ((~i & 0xf) << 1) & 3;
		if ((i & 0x7f) == 0) {
			if (i > 0) {
				memcpy(&bwt->data[k], x, 32);
				k += 4;
			}
			last_c = &bwt->data[k];
			memcpy(&bwt->data[k], c4, 32);
			k += 4;
			memset(x, 0, 32);
		} else if ((i & 0x3f) == 0 && last_c) {
			for (j = 0; j < 4; ++j)
				last_c[j] |= (c4[j] - last_c[j]) << BWT_CNT_SHIFT;
		}
		++c4[a];
		x[(i&0x7f)>>5] |= (uint64_t)a << ((i&0x1f)<<1);
	}
	memcpy(&bwt->data[k], x, 32); k += 4;
	memcpy(&bwt->data[k], c4, 32); k += 4;
	kom_assert(k == bwt->data_len, "inconsistent BWT length");
	for (i = 0, bwt->L2[0] = 0; i < 4; ++i)
		bwt->L2[i+1] = bwt->L2[i] + c4[i];
	kom_assert(bwt->L2[4] == len, "inconsistent base counts");
	free(chunk);
	fclose(fp);
	fprintf(stderr, "[M::%s] encoded BWT of %lld symbols in %.1f sec\n", __func__, (long long)len, kom_realtime() - t);

	t = kom_realtime();
	mb_bwt_gen_sa(bwt, sa_bit);
	fprintf(stderr, "[M::%s] sampled SA (1/%d) in %.1f sec with %d thread(s)\n", __func__, 1<<sa_bit, kom_realtime() - t, mb_n_threads);
	report_heap("raw2mbw");
	if (mb_bwt_save(argv[o.ind+1], bwt) != 0) return 1;
	mb_bwt_destroy(bwt);
	return 0;
}

/********
 * main *
 ********/

int mbw_main(int argc, char *argv[])
{
	int ret = 0;
	kom_realtime();
	if (argc > 2 && strcmp(argv[1], "-B") == 0) { // testing: put a ballast below everything else so that buffers sit above 2GB
		size_t sz = (size_t)atol(argv[2]) << 20;
		void *p = malloc(sz);
		fprintf(stderr, "[wasm] ballast of %ld MB at %p\n", (long)(sz>>20), p);
		argc -= 2, argv += 2;
	}
	if (argc > 2 && strcmp(argv[1], "-C") == 0) { // size of the reference page cache (MB) used by map in lazy mode
		l2b_cache_bytes = (int64_t)atol(argv[2]) << 20;
		argc -= 2, argv += 2;
	}
	if (argc == 1) {
		fprintf(stderr, "Usage: minibwa [-C cacheMB] <prep|genraw|raw2mbw|index|map|mem|fastmap|version> ...\n");
		return 1;
	}
	if (strcmp(argv[1], "prep") == 0) ret = main_prep(argc-1, argv+1);
	else if (strcmp(argv[1], "genraw") == 0) ret = main_genraw(argc-1, argv+1), report_heap("genraw");
	else if (strcmp(argv[1], "raw2mbw") == 0) ret = main_raw2mbw(argc-1, argv+1);
	else if (strcmp(argv[1], "index") == 0) ret = main_index(argc-1, argv+1), report_heap("index");
	else if (strcmp(argv[1], "map") == 0) {
		l2b_lazy_pac = 1, ret = main_map(argc-1, argv+1), report_heap("map");
		fprintf(stderr, "[wasm] reference page cache (%lld MB): %lld hits, %lld misses\n", (long long)(l2b_cache_bytes>>20), (long long)l2b_cache_hit, (long long)l2b_cache_miss);
	}
	else if (strcmp(argv[1], "map-inmem") == 0) ret = main_map(argc-1, argv+1), report_heap("map");
	else if (strcmp(argv[1], "mem") == 0) l2b_lazy_pac = 1, ret = main_mem(argc-1, argv+1), report_heap("mem");
	else if (strcmp(argv[1], "fastmap") == 0) ret = main_fastmap(argc-1, argv+1);
	else if (strcmp(argv[1], "version") == 0) printf("%s-wasm\n", MB_VERSION);
	else {
		fprintf(stderr, "ERROR: unknown command '%s'\n", argv[1]);
		return 1;
	}
	if (kom_verbose >= 3 && ret == 0)
		fprintf(stderr, "[M::main] Real time: %.3f sec; CPU: %.3f sec\n", kom_realtime(), kom_cputime());
	fflush(stdout); fflush(stderr);
	return ret;
}

#ifndef MBW_REACTOR
int main(int argc, char *argv[]) { return mbw_main(argc, argv); }
#else
/****************************************************************************
 * Reactor exports (the module is instantiated once per Web Worker)          *
 *                                                                          *
 *   mbw_run(argc, buf)          run one CLI-style command (index stages)     *
 *   mbw_align_open(argc, buf)   parse `map` options, load the index and keep *
 *                               it in memory (BWT + SA; reference on OPFS)   *
 *   mbw_align_header(path)      write the SAM header                         *
 *   mbw_align_batch(n, path)    map the next n records per input file and    *
 *                               write them to path; returns #records, 0=EOF  *
 *   mbw_align_close()           free everything                              *
 * Strings are passed as NUL-separated arguments in memory from mbw_alloc(). *
 ****************************************************************************/
#define EXPORT(name) __attribute__((export_name(#name)))
extern int mb_batch_keep;
int mb_batch_write_header(const char *fn);
int64_t mb_batch_run(int64_t max_reads, const char *fn_out);
void mb_batch_close(void);
int mb_batch_set_reads(int32_t n, const char **fn);

static char **split_args(int argc, char *buf)
{
	char **argv = (char**)calloc(argc + 2, sizeof(char*));
	int i;
	argv[0] = "minibwa";
	for (i = 0; i < argc; ++i) {
		argv[i + 1] = buf;
		buf += strlen(buf) + 1;
	}
	return argv;
}

EXPORT(mbw_alloc) void *mbw_alloc(size_t n) { return malloc(n); }
EXPORT(mbw_free) void mbw_free(void *p) { free(p); }
EXPORT(mbw_memory) double mbw_memory(void) { return (double)__builtin_wasm_memory_size(0) * 65536.0; }

EXPORT(mbw_run) int mbw_run(int argc, char *buf)
{
	char **argv = split_args(argc, buf);
	int ret = mbw_main(argc + 1, argv);
	free(argv);
	return ret;
}

EXPORT(mbw_align_open) int mbw_align_open(int argc, char *buf) // buf: [-C cacheMB | -P] "map" followed by map options, index prefix and reads
{
	char **argv = split_args(argc, buf);
	int ret, i = 1;
	kom_realtime();
	l2b_lazy_pac = 1;
	if (argc > 2 && strcmp(argv[i], "-C") == 0) l2b_cache_bytes = (int64_t)atol(argv[i+1]) << 20, i += 2;
	if (argc > 1 && strcmp(argv[i], "-P") == 0) l2b_lazy_pac = 0, l2b_skip_mask = 1, ++i; // [mt] whole 2-bit reference in memory: worker threads never do I/O
	mb_batch_keep = 1;
	ret = main_map(argc + 1 - i, argv + i);
	mb_batch_keep = 0;
	free(argv);
	report_heap("index loaded");
	return ret;
}

EXPORT(mbw_align_reads) int mbw_align_reads(int argc, char *buf) // next read files (1 or 2) for the loaded index
{
	char **argv = split_args(argc, buf);
	int ret = mb_batch_set_reads(argc, (const char**)argv + 1);
	free(argv);
	return ret;
}

EXPORT(mbw_align_header) int mbw_align_header(const char *fn) { return mb_batch_write_header(fn); }

EXPORT(mbw_align_batch) double mbw_align_batch(double n, const char *fn) // double: exact for these sizes and avoids BigInt in JS
{
	int64_t r = mb_batch_run((int64_t)n, fn);
	fflush(stderr);
	return (double)r;
}

EXPORT(mbw_align_close) void mbw_align_close(void)
{
	fprintf(stderr, "[wasm] reference page cache (%lld MB): %lld hits, %lld misses\n", (long long)(l2b_cache_bytes>>20), (long long)l2b_cache_hit, (long long)l2b_cache_miss);
	report_heap("align");
	mb_batch_close();
	fflush(stderr);
}
#endif
