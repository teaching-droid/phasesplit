/* A minimal parallel-for. Enough to spread independent work across cores
 * without pulling in a threading library.
 */
#ifndef PS_THREAD_H
#define PS_THREAD_H

/* Number of hardware threads, or 1 if that cannot be determined. */
int ps_cpu_threads(void);

/* Calls fn(i, ctx) for i in [0, count), spread over `threads` workers.
 * Each index is handled exactly once, so fn must only touch data belonging
 * to its own index. Runs inline when threads <= 1 or count is small. */
void ps_parallel_for(int count, void (*fn)(int, void *), void *ctx, int threads);

#endif
