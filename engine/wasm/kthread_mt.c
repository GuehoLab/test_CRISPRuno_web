/* [mt] Threads for the multi-threaded wasm64 build (emscripten pthreads):
 * kt_for() is the real multi-threaded one from kthread.c, kt_pipeline() stays serial so that all
 * file I/O (done by the pipeline steps) happens on the main thread, where the OPFS handles live. */
#define kt_pipeline kt_pipeline_threaded_unused
#include "kthread.c"
#undef kt_pipeline

void kt_pipeline(int n_threads, void *(*func)(void*, int, void*), void *shared_data, int n_steps)
{
	(void)n_threads;
	for (;;) {
		int step;
		void *in = func(shared_data, 0, 0);
		if (in == 0) break;
		for (step = 1; step < n_steps && in; ++step)
			in = func(shared_data, step, in);
	}
}
