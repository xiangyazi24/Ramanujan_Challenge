#!/usr/bin/env python3
"""Q6400 verification: Lucas collapse/bound, central binomial calibration,
and the Lemma-A (finite-offset apparition surrogate) falsification test on Apery hits."""
import math, time
from sympy import primerange, isprime

# ---------- V1: Lucas exact collapse (8.3) + Fibonacci bound (8.4) ----------
def fib_pair(m, p):
    """(F_m, F_{m+1}) mod p by fast doubling."""
    if m == 0: return (0, 1)
    a, b = fib_pair(m >> 1, p)
    c = a * ((2*b - a) % p) % p
    d = (a*a + b*b) % p
    return (d, (c + d) % p) if (m & 1) else (c, d)

def rho_F(p):
    """rank of apparition of Fibonacci mod p (p != 5): least d>0, p | F_d; rho | p - (5|p)."""
    chi = pow(5, (p-1)//2, p)
    chi = 1 if chi == 1 else -1
    N = p - chi
    # rho divides N: factor N, strip primes
    d = N
    n_ = N; f = {}
    x = n_
    for q in range(2, int(n_**0.5)+1):
        while x % q == 0: f[q] = f.get(q,0)+1; x //= q
    if x > 1: f[x] = f.get(x,0)+1
    for q in f:
        while d % q == 0 and fib_pair(d//q, p)[0] == 0:
            d //= q
    return d, chi

print("== V1: Lucas collapse + Fibonacci ==")
bad = 0
ps = [p for p in primerange(7, 2000) if p != 5]
for p in ps:
    rho, chi = rho_F(p)
    assert (p - chi) % rho == 0, (p, rho, chi)
n = 3001  # test the collapse: p | F_{n-p} <=> rho | n-chi
viol = 0
hits = 0
for p in primerange(n//2+1, n+1):
    if p == 5: continue
    rho, chi = rho_F(p)
    lhs = fib_pair(n-p, p)[0] == 0 if n-p > 0 else True
    rhs = ((n - chi) % rho == 0)
    if lhs != rhs: viol += 1
    if lhs: hits += 1
print(f"  rho | p-chi: {len(ps)}/{len(ps)} OK; collapse at n={n}: violations={viol}, H_F(n)={hits}")
# empirical H_F vs bound over a range of n
worst = (0, 0, 0)
for n in range(1000, 4001, 111):
    h = 0
    for p in primerange(n//2+1, n+1):
        if p == 5: continue
        if n - p > 0 and fib_pair(n-p, p)[0] == 0: h += 1
    # bound (5.13)
    def tau(m):
        t = 1; x = m
        for q in range(2, int(m**0.5)+1):
            e = 0
            while x % q == 0: e += 1; x //= q
            t *= e+1
        if x > 1: t *= 2
        return t
    L = math.log(n/2); cU = math.log((1+5**0.5)/2)
    bound = (tau(n-1)+tau(n+1))*(math.sqrt(2*n*cU/L)+1)
    if h > worst[1]: worst = (n, h, bound)
    assert h <= bound, (n, h, bound)
print(f"  H_F(n) <= bound(5.13) for n=1000..4000 step 111; worst case n={worst[0]}: H={worst[1]} vs bound={worst[2]:.0f}")

# ---------- V2: central binomial claim (8.2) ----------
print("== V2: central binomial ==")
for n in (500, 1234, 3000):
    lhs = set()
    for p in primerange(n//2+1, n+1):
        r = n - p
        # p | C(2r, r) iff carry adding r+r base p iff (2r >= p) for r < p
        c = math.comb(2*r, r) if r <= 600 else None
        kummer = (2*r >= p)
        if c is not None:
            assert (c % p == 0) == kummer, (n, p)
        if kummer: lhs.add(p)
    lo, hi = n/2, 2*n/3
    pred = sum(1 for p in primerange(int(lo)+1, n+1) if p <= hi)
    print(f"  n={n}: hits={len(lhs)} vs primes in (n/2,2n/3]={pred} (direct comb check where r<=600: OK)")

# ---------- V3: Lemma A falsification on Apery hits ----------
print("== V3: Lemma A (finite-offset apparition surrogate) falsification ==")
def apery_row(p):
    Z = []
    b = [0]*(p)
    b[0] = 1 % p
    if p > 1: b[1] = 5 % p
    inv = [0]*(p)
    inv[1] = 1
    for i in range(2, p): inv[i] = (p - (p//i)*inv[p % i]) % p
    for n_ in range(1, p-1):
        num = ((34*n_**3 + 51*n_**2 + 27*n_ + 5)*b[n_] - n_**3*b[n_-1]) % p
        b[n_+1] = num * pow(inv[n_+1], 3, p) % p
    for r in range(p):
        if b[r] == 0: Z.append(r)
    return Z, b
E = [-1, 0, 1]
uncaught = []
total_hits = 0
t0 = time.time()
for p in primerange(7, 4000):
    Z, b = apery_row(p)
    for r in Z:
        total_hits += 1
        # certificate: exists d>1, d|r, d|p-eps (eps in E), p | b_d  (carrier C_d = b_d)
        found = False
        if r == 0: continue  # r=0 -> b_0=1, never a hit anyway
        for eps in E:
            m = p - eps
            for d in range(2, r+1):
                if r % d == 0 and m % d == 0 and d < p and b[d] == 0:
                    found = True; break
            if found: break
        if not found:
            uncaught.append((p, r))
print(f"  Apery hits p<4000: total={total_hits}, uncaught by (d|r, d|p-eps, p|b_d, eps in {{0,±1}}): {len(uncaught)}")
print(f"  first uncaught examples: {uncaught[:8]}")
print(f"  => Lemma A with carrier C_d=b_d and E={{0,±1}}: {'FALSIFIED' if uncaught else 'survives (!!)'}  ({time.time()-t0:.0f}s)")
