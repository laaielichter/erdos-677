"""Independent slow check for small X: exact big-part pairs and exact lcm solutions.
Usage: python3 check_small.py X KMIN KMAX"""
import math, sys
from collections import defaultdict

X, KMIN, KMAX = (int(a) for a in sys.argv[1:4])

spf = list(range(X + 1))
for i in range(2, int(X ** 0.5) + 1):
    if spf[i] == i:
        for j in range(i * i, X + 1, i):
            if spf[j] == j:
                spf[j] = i

def big(y):  # part of y made of primes > 100
    r = 1
    while y > 1:
        p = spf[y]
        y //= p
        if p > 100:
            r *= p
    return r

bigs = [0] + [big(y) for y in range(1, X + 1)]
for k in range(KMIN, KMAX + 1):
    groups = defaultdict(list)
    for x in range(0, X - k + 1):
        groups[math.prod(bigs[x + 1:x + k + 1])].append(x)
    gp = sizeok = sols = 0
    lg = math.lgamma(k)
    for xs in groups.values():
        for a in range(len(xs)):
            n = xs[a]
            for b in range(a + 1, len(xs)):
                m = xs[b]
                if m - n < k:
                    continue
                gp += 1
                lm = sum(math.log(m + q) for q in range(1, k + 1))
                ln = sum(math.log(n + q) for q in range(1, k + 1))
                if lm - lg > ln + 1e-6:
                    continue
                sizeok += 1
                if math.lcm(*range(n + 1, n + k + 1)) == math.lcm(*range(m + 1, m + k + 1)):
                    sols += 1
                    print("SOLUTION", k, n, m)
    print(f"k={k} bigpart_pairs={gp} size_ok={sizeok} solutions={sols}")
