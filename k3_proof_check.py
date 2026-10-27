"""Rigorous computer part of the k=3 proof.
Step A: solve the Thue equations p^3 - 2 q^3 = N for N in {-1,-2,-3,-6} with PARI (unconditional mode).
Step B: from each (p,q) with q < p < 2q, gcd=1, build g^2 = (2q-p)/(2q^3-p^3) and the solutions (X,Y)=(gq,gp).
Step C: brute-force sanity check of Y^3-Y = 2(X^3-X) for X up to 10^6."""
import cypari2, math
pari = cypari2.Pari()
tnf = pari('thueinit(x^3 - 2, 1)')   # flag 1 = unconditional (no GRH)
found = []
for N in (-1, -2, -3, -6):
    sols = pari.thue(tnf, N)
    sols = [(int(s[0]), int(s[1])) for s in sols]
    print(f"p^3 - 2q^3 = {N}: all integer solutions {sols}")
    for p, q in sols:
        if q >= 1 and q < p < 2 * q and math.gcd(p, q) == 1:
            num, den = 2 * q - p, 2 * q**3 - p**3
            if num % den == 0 and math.isqrt(num // den) ** 2 == num // den:
                g = math.isqrt(num // den)
                X, Y = g * q, g * p
                assert Y**3 - Y == 2 * (X**3 - X)
                found.append((X, Y))
print("solutions (X,Y) of Y^3-Y = 2(X^3-X) with X>=2:", found)
for X, Y in found:
    a, b = X - 1, Y - 1
    print(f"  -> a={a}, b={b}: {b}*{b+1}*{b+2} = 2*{a}*{a+1}*{a+2}; windows n={a-1}, m={b-1}, m-n={b-a}")
# Step C
brute = [(X, Y) for X in range(2, 10**6) for Y in [round(2 ** (1/3) * X) + d for d in (-1, 0, 1)]
         if Y**3 - Y == 2 * (X**3 - X)]
print("brute force X < 10^6:", brute)
