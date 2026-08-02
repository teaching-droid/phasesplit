#include "thread.h"

#include <stdlib.h>

#if defined(_WIN32)
  #define WIN32_LEAN_AND_MEAN
  #include <windows.h>
#else
  #include <pthread.h>
  #include <unistd.h>
#endif

int ps_cpu_threads(void) {
#if defined(_WIN32)
    SYSTEM_INFO si;
    GetSystemInfo(&si);
    return si.dwNumberOfProcessors > 0 ? (int)si.dwNumberOfProcessors : 1;
#else
    long n = sysconf(_SC_NPROCESSORS_ONLN);
    return n > 0 ? (int)n : 1;
#endif
}

typedef struct {
    void  (*fn)(int, void *);
    void   *ctx;
    int     count;
    int     next;        /* next index to hand out */
#if defined(_WIN32)
    CRITICAL_SECTION lock;
#else
    pthread_mutex_t  lock;
#endif
} ps_pool;

/* Work is taken one index at a time rather than split up front, so an uneven
 * piece cannot leave one worker running long after the others have finished. */
static int take_next(ps_pool *p) {
    int i;
#if defined(_WIN32)
    EnterCriticalSection(&p->lock);
    i = (p->next < p->count) ? p->next++ : -1;
    LeaveCriticalSection(&p->lock);
#else
    pthread_mutex_lock(&p->lock);
    i = (p->next < p->count) ? p->next++ : -1;
    pthread_mutex_unlock(&p->lock);
#endif
    return i;
}

static void run_worker(ps_pool *p) {
    for (;;) {
        int i = take_next(p);
        if (i < 0) return;
        p->fn(i, p->ctx);
    }
}

#if defined(_WIN32)
static DWORD WINAPI worker_entry(LPVOID arg) { run_worker((ps_pool *)arg); return 0; }
#else
static void *worker_entry(void *arg) { run_worker((ps_pool *)arg); return NULL; }
#endif

void ps_parallel_for(int count, void (*fn)(int, void *), void *ctx, int threads) {
    if (count <= 0 || !fn) return;
    if (threads > count) threads = count;

    if (threads <= 1) {
        for (int i = 0; i < count; i++) fn(i, ctx);
        return;
    }

    ps_pool p;
    p.fn = fn; p.ctx = ctx; p.count = count; p.next = 0;

#if defined(_WIN32)
    InitializeCriticalSection(&p.lock);
    HANDLE *th = (HANDLE *)malloc((size_t)(threads - 1) * sizeof(HANDLE));
    int started = 0;
    if (th) {
        for (int i = 0; i < threads - 1; i++) {
            th[i] = CreateThread(NULL, 0, worker_entry, &p, 0, NULL);
            if (th[i]) started++;
            else break;
        }
    }
    run_worker(&p);                       /* the calling thread works too */
    for (int i = 0; i < started; i++) {
        WaitForSingleObject(th[i], INFINITE);
        CloseHandle(th[i]);
    }
    free(th);
    DeleteCriticalSection(&p.lock);
#else
    pthread_mutex_init(&p.lock, NULL);
    pthread_t *th = (pthread_t *)malloc((size_t)(threads - 1) * sizeof(pthread_t));
    int started = 0;
    if (th) {
        for (int i = 0; i < threads - 1; i++) {
            if (pthread_create(&th[i], NULL, worker_entry, &p) == 0) started++;
            else break;
        }
    }
    run_worker(&p);
    for (int i = 0; i < started; i++) pthread_join(th[i], NULL);
    free(th);
    pthread_mutex_destroy(&p.lock);
#endif
}
