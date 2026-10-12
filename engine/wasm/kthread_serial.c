/* Serial drop-in replacement of kthread.c for single-threaded WebAssembly builds.
 * kt_for() runs the loop in the calling thread; kt_pipeline() runs the steps
 * one after another for each batch. Output is identical to the threaded version. */
#include "kthread.h"

void kt_for(int n_threads, void (*func)(void*,long,int), void *data, long n)
{
	long i;
	(void)n_threads;
	for (i = 0; i < n; ++i) func(data, i, 0);
}

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
