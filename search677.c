// Erdős Problem #677 search.
// M(n,k) = lcm(n+1, ..., n+k). Look for M(n,k) == M(m,k) with m >= n+k.
// Checks every pair n < m with m+k <= X, for KMIN <= k <= KMAX (KMAX <= 100).
//
// Method:
//  1. For each y <= X, h(y) = sum over primes p > 100 (with multiplicity) of r(p),
//     where r is a fixed 64-bit pseudo-random function. B = prefix sums of h.
//     For p > 100 >= k, a window of k consecutive integers has at most one multiple
//     of p, so G(x,k) = B[x+k]-B[x] depends only on the "p > 100 part" of M(x,k).
//     A solution forces G(n,k) == G(m,k): no false negatives.
//  2. Sort all G(x,k), find equal values, list the pairs with m-n >= k.
//  3. Necessary size test (rigorous): (m+1)*C(m+k,k-1) divides M(m,k), so
//     prod(m+j)/(k-1)! <= M(m,k) = M(n,k) <= prod(n+i).
//  4. Exact check of the full lcm via prime factorisation (all primes).
//     Survivors are printed as SOLUTION and re-checked in Python with exact integers.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <omp.h>

typedef uint64_t u64;
typedef uint32_t u32;

#define BIGP 100

static u32 *spf;   // smallest prime factor (spf[p] = p for primes)
static u64 *B;     // prefix sums of h
static long X;

static inline u64 rnd(u64 p) {           // splitmix64
    u64 z = p * 0x9E3779B97F4A7C15ULL + 0x632BE59BD9B4E019ULL;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

// LSD radix sort, 4 passes of 16 bits
static void radix_sort(u64 *a, u64 *t, long n) {
    long (*cnt)[65536] = calloc(4, sizeof(*cnt));
    for (long i = 0; i < n; i++) {
        u64 v = a[i];
        cnt[0][v & 0xFFFF]++; cnt[1][(v >> 16) & 0xFFFF]++;
        cnt[2][(v >> 32) & 0xFFFF]++; cnt[3][v >> 48]++;
    }
    for (int p = 0; p < 4; p++) {
        long s = 0;
        for (int b = 0; b < 65536; b++) { long c = cnt[p][b]; cnt[p][b] = s; s += c; }
        int sh = 16 * p;
        for (long i = 0; i < n; i++) { u64 v = a[i]; t[cnt[p][(v >> sh) & 0xFFFF]++] = v; }
        u64 *tmp = a; a = t; t = tmp;
    }
    free(cnt);
    // 4 passes: result is back in the original array
}

typedef struct { u64 key; long x; } KX;
static int cmp_kx(const void *A, const void *Bv) {
    const KX *a = A, *b = Bv;
    if (a->key != b->key) return a->key < b->key ? -1 : 1;
    return (a->x > b->x) - (a->x < b->x);
}

typedef struct { u32 p; int e; } PE;
static int cmp_pe(const void *A, const void *Bv) {
    const PE *a = A, *b = Bv;
    if (a->p != b->p) return a->p < b->p ? -1 : 1;
    return b->e - a->e;  // larger exponent first
}

// Exact fingerprint of lcm(s+1..s+k) using ALL primes (max exponent per prime).
static u64 full_fp(long s, int k, PE *buf) {
    int c = 0;
    for (int i = 1; i <= k; i++) {
        u32 y = (u32)(s + i);
        while (y > 1) {
            u32 p = spf[y]; int e = 0;
            while (y % p == 0) { y /= p; e++; }
            buf[c].p = p; buf[c].e = e; c++;
        }
    }
    qsort(buf, c, sizeof(PE), cmp_pe);
    u64 f = 0;
    for (int i = 0; i < c; i++)
        if (i == 0 || buf[i].p != buf[i-1].p) f += (u64)buf[i].e * rnd(buf[i].p) * 0x2545F4914F6CDD1DULL;
    return f;
}

int main(int argc, char **argv) {
    X = argc > 1 ? atol(argv[1]) : 100000000L;
    int KMIN = argc > 2 ? atoi(argv[2]) : 3;
    int KMAX = argc > 3 ? atoi(argv[3]) : 100;
    if (KMAX > BIGP) { fprintf(stderr, "KMAX must be <= %d\n", BIGP); return 1; }
    double t0 = omp_get_wtime();

    spf = calloc(X + 1, sizeof(u32));
    for (long i = 2; i <= X; i++) {
        if (spf[i]) continue;
        spf[i] = (u32)i;
        if (i <= X / i) for (long j = i * i; j <= X; j += i) if (!spf[j]) spf[j] = (u32)i;
    }
    B = malloc((X + 1) * sizeof(u64));
    B[0] = 0;
    for (long y = 1; y <= X; y++) {
        u64 h = 0; u32 z = (u32)y;
        while (z > 1) { u32 p = spf[z]; z /= p; if (p > BIGP) h += rnd(p); }
        B[y] = B[y - 1] + h;
    }
    fprintf(stderr, "sieve+prefix done in %.1fs\n", omp_get_wtime() - t0);
    printf("X=%ld k=%d..%d\n", X, KMIN, KMAX); fflush(stdout);

    #pragma omp parallel for schedule(dynamic, 1)
    for (int k = KMIN; k <= KMAX; k++) {
        double tk = omp_get_wtime();
        long nw = X - k + 1;  // windows x = 0..X-k
        u64 *a = malloc(nw * sizeof(u64)), *t = malloc(nw * sizeof(u64));
        for (long x = 0; x < nw; x++) a[x] = B[x + k] - B[x];
        radix_sort(a, t, nw);
        // distinct duplicated keys -> hash set
        long nd = 0;
        for (long i = 1; i < nw; i++) if (a[i] == a[i-1] && (i == 1 || a[i-1] != a[i-2])) nd++;
        long hs = 1; while (hs < 4 * nd + 16) hs <<= 1;
        u64 *set = calloc(hs, sizeof(u64)); char *used = calloc(hs, 1);
        for (long i = 1; i < nw; i++) if (a[i] == a[i-1] && (i == 1 || a[i-1] != a[i-2])) {
            u64 v = a[i]; long h = (long)(rnd(v) & (hs - 1));
            while (used[h]) h = (h + 1) & (hs - 1);
            used[h] = 1; set[h] = v;
        }
        free(t); free(a);
        // collect positions with a duplicated key
        long cap = 1024, nl = 0; KX *L = malloc(cap * sizeof(KX));
        if (nd) for (long x = 0; x < nw; x++) {
            u64 v = B[x + k] - B[x]; long h = (long)(rnd(v) & (hs - 1));
            while (used[h]) {
                if (set[h] == v) {
                    if (nl == cap) { cap *= 2; L = realloc(L, cap * sizeof(KX)); }
                    L[nl].key = v; L[nl].x = x; nl++; break;
                }
                h = (h + 1) & (hs - 1);
            }
        }
        free(set); free(used);
        qsort(L, nl, sizeof(KX), cmp_kx);
        double lgk = lgamma((double)k);  // log((k-1)!)
        long gpairs = 0, sizeok = 0, sols = 0, nearshown = 0;
        PE *buf = malloc(sizeof(PE) * 64 * k);
        char line[256];
        for (long i = 0; i < nl; ) {
            long j = i; while (j < nl && L[j].key == L[i].key) j++;
            for (long u = i; u < j; u++) {
                // logs of window products
                long n = L[u].x;
                double ln = 0; for (int q = 1; q <= k; q++) ln += log((double)(n + q));
                for (long v = u + 1; v < j; v++) {
                    long m = L[v].x;
                    if (m - n < k) continue;
                    gpairs++;
                    double lm = 0; for (int q = 1; q <= k; q++) lm += log((double)(m + q));
                    if (lm - lgk > ln + 1e-6) continue;   // too big: impossible
                    sizeok++;
                    if (full_fp(n, k, buf) == full_fp(m, k, buf)) {
                        sols++;
                        snprintf(line, sizeof line, "SOLUTION k=%d n=%ld m=%ld\n", k, n, m);
                        #pragma omp critical
                        { fputs(line, stdout); fflush(stdout); }
                    } else if (nearshown < 5) {
                        nearshown++;
                        snprintf(line, sizeof line, "NEAR k=%d n=%ld m=%ld\n", k, n, m);
                        #pragma omp critical
                        { fputs(line, stdout); fflush(stdout); }
                    }
                }
            }
            i = j;
        }
        free(buf); free(L);
        snprintf(line, sizeof line, "k=%d windows=%ld dupkeys=%ld bigpart_pairs=%ld size_ok=%ld solutions=%ld (%.1fs)\n",
                 k, nw, nd, gpairs, sizeok, sols, omp_get_wtime() - tk);
        #pragma omp critical
        { fputs(line, stdout); fflush(stdout); }
    }
    fprintf(stderr, "total %.1fs\n", omp_get_wtime() - t0);
    return 0;
}
