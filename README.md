# Erdős Problem #677: computer search and the case k = 3

[Erdős Problem #677](https://www.erdosproblems.com/677): let $M(n,k)=\mathrm{lcm}(n+1,\dots,n+k)$. Is it true that $M(n,k)\ne M(m,k)$ for all $m\ge n+k$?

## Results

1. **Search.** For every $3\le k\le 100$, there is no pair $0\le n<m$ with $m\ge n+k$, $m+k\le 10^8$ and $M(n,k)=M(m,k)$. (For $k\le 2$ the statement is trivial.)
2. **The case $k=3$ is settled.** $M(n,3)\ne M(m,3)$ for all $m\ge n+3$. The case reduces to $b(b+1)(b+2)=2a(a+1)(a+2)$. Its only positive solution is $(a,b)=(3,4)$ (T. H. Chan, [arXiv:2408.01306](https://arxiv.org/abs/2408.01306), Period. Math. Hungar. 2025). That solution gives $M(2,3)=M(3,3)=60$, which has $m-n=1$. A short independent proof, reducing to the Thue equations $p^3-2q^3\in\{-1,-2,-3,-6\}$, is in [`forum_post.md`](forum_post.md) and is checked by [`k3_proof_check.py`](k3_proof_check.py).

## Files

| File | What it does |
| --- | --- |
| `search677.c` | Fast search (C + OpenMP). The method is explained at the top of the file. |
| `check_small.py` | Slow, independent Python check for small ranges, used to cross-check the C program. |
| `k3_proof_check.py` | Solves the Thue equations with PARI (unconditional) and finishes the $k=3$ proof. |
| `run_1e8.txt` | Output of the full run ($m+k\le 10^8$, $3\le k\le 100$). |
| `forum_post.md` | Write-up posted on the erdosproblems.com forum. |

## Reproduce

```sh
gcc -O3 -march=native -fopenmp -o search677 search677.c -lm
./search677 100000000 3 100 > run_1e8.txt     # ~9 min on 2 cores, ~4.5 GB RAM
python3 check_small.py 30000 3 12             # compare with: ./search677 30000 3 12
pip install cypari2 && python3 k3_proof_check.py
```

## Method in short

For a prime $p>100\ge k$, a window of $k$ consecutive integers contains at most one multiple of $p$. So the $p>100$ part of $M(n,k)$ is the product of the $p>100$ parts of the window's elements. A random completely additive hash over primes $p>100$, stored as prefix sums, gives a necessary condition with no false negatives. The pairs that survive are filtered with the rigorous bound $\prod_j (m+j)/(k-1)! \le M(m,k)$, which holds because $(m+1)\binom{m+k}{k-1}$ divides $\mathrm{lcm}(m+1,\dots,m+k)$. Then they are compared exactly using full factorisations. In total there were 84,082 exact comparisons and none were equal.

Made with help from Claude (Anthropic).

## License

MIT
