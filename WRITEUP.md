**Computer search, and a complete answer for k = 3**

Write $M(n,k)=\mathrm{lcm}(n+1,\dots,n+k)$.

**1. Search.** For every $3\le k\le 100$ I checked all pairs $0\le n<m$ with $m\ge n+k$ and $m+k\le 10^8$. There is no pair with $M(n,k)=M(m,k)$.

Method: for a prime $p>100\ge k$, a window of $k$ consecutive integers contains at most one multiple of $p$. So the "$p>100$ part" of $M(n,k)$ is the product of the $p>100$ parts of $n+1,\dots,n+k$. I hashed this with prefix sums of a random completely additive function supported on primes $p>100$. Equal $M$ forces equal hashes, so this filter has no false negatives. Every pair with equal hash and $m\ge n+k$ was then tested against the size bound
$$\frac{(m+1)\cdots(m+k)}{(k-1)!}\le M(m,k)=M(n,k)\le (n+1)\cdots(n+k),$$
which holds because $(m+1)\binom{m+k}{k-1}$ divides $\mathrm{lcm}(m+1,\dots,m+k)$ (partial fractions). The survivors (84,082 pairs over all $k$) were compared exactly via full factorisations, and none were equal. The code agrees with an independent slow Python implementation for $m\le 3\cdot 10^4$, $k\le 12$. Runtime was about 9 minutes on 2 cores.

**2. The case $k=3$ is settled: $M(n,3)\ne M(m,3)$ for all $m\ge n+3$.**

Put $a=n+1$. Then $\mathrm{lcm}(a,a+1,a+2)$ equals $a(a+1)(a+2)$ if $a$ is odd and $a(a+1)(a+2)/2$ if $a$ is even. Since $m>n$, an equality $M(n,3)=M(m,3)$ forces $b(b+1)(b+2)=2a(a+1)(a+2)$ with $b=m+1$. By Theorem 1 of T. H. Chan (arXiv:2408.01306, Period. Math. Hungar. 2025), the only positive solution is $(a,b)=(3,4)$. That is $M(2,3)=M(3,3)=60$, which has $m-n=1<3$.

A short independent argument: let $X=a+1$, $Y=b+1$, so $Y^3-Y=2(X^3-X)$ and $X<Y<2X$. Write $g=\gcd(X,Y)$, $X=gq$, $Y=gp$. Then $g^2(p^3-2q^3)=p-2q$. With $N=p^3-2q^3\ne 0$ we get $p\equiv 2q \pmod N$, so $6q^3\equiv 0\pmod N$. Since $\gcd(q,N)=1$, this gives $N\mid 6$, and $N<0$ because $p<2q$. So $p^3-2q^3\in\{-1,-2,-3,-6\}$. PARI's `thue` (unconditional mode) gives all solutions: $(-1,0),(1,1)$; $(0,1)$; $(-1,1),(5,4)$; $(-2,-1)$. Only $(p,q)=(5,4)$ has $q<p<2q$. It gives $g=1$ and $(X,Y)=(4,5)$, the overlapping example above.

The search and this write-up were done with help from Claude (Anthropic). Code and output: https://github.com/laaielichter/erdos-677
