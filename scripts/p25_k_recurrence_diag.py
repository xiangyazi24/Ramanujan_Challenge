#!/usr/bin/env python3
"""
P2.5 diagnostic companion to p25_k_recurrence.py.

1. Factor the denominators of f(k): holonomic rational sequences have
   denominators built from finitely many hypergeometric-type factors
   (values of the recurrence's leading coefficient). Irregular prime
   content = non-holonomic fingerprint.
2. Search recurrences for gauge rescalings:
     u(k) = 2^k C(2k,k) f(k)      (Legendre-gauge coefficient)
     and variants with (2k+1), 4^k, C(2k,k)^2 factors.
"""
from fractions import Fraction as F
from math import comb, gcd, lcm, isqrt
from functools import reduce
import sys, time

sys.setrecursionlimit(10000)

KMAX = 140

# --- copy of CMF and decomposition (same as p25_k_recurrence.py) ---
def M_entries(n):
    m11 = (-2*n-5)*(n+3)**2*(136*n**4+1424*n**3+5548*n**2+9551*n+6141)
    m12 = 384*n**6+6384*n**5+44168*n**4+162698*n**3+336377*n**2+369933*n+169011
    m13 = -(480*n**4+4980*n**3+19210*n**2+32690*n+20730)
    m21 = (n+2)**2*(n+3)**2*(4*n+10)*(48*n**3+386*n**2+1017*n+879)
    m22 = (n+2)**2*(-272*n**5-3848*n**4-21732*n**3-61184*n**2-85761*n-47808)
    m23 = (n+2)**2*(320*n**3+2540*n**2+6610*n+5640)
    m31 = (-4*n-10)*(n+2)**2*(n+3)**2*(32*n**4+302*n**3+1037*n**2+1530*n+813)
    m32 = (n+2)**2*(192*n**6+2984*n**5+19116*n**4+64452*n**3+120256*n**2+117279*n+46476)
    m33 = (n+2)**2*(-16*n**5-408*n**4-2912*n**3-8884*n**2-12254*n-6240)
    return [[m11,m12,m13],[m21,m22,m23],[m31,m32,m33]]

def delta_H(n):
    return -2*(n+2)**2*(n+3)**2*(2*n+5)*(2*n+7)**2

def delannoy_summand(N, k):
    if k < 0 or k > N:
        return F(0)
    return F(2**k * comb(2*k, k) * comb(N, k) * comb(N+k, k))

t0 = time.time()
rows = {'e1': [F(1), F(0), F(0)]}
history = {key: [list(v)] for key, v in rows.items()}
for N in range(KMAX):
    M = M_entries(N)
    d = F(delta_H(N))
    MH = [[F(M[i][j]) / d for j in range(3)] for i in range(3)]
    for key in rows:
        r = rows[key]
        rows[key] = [sum(r[i]*MH[i][k] for i in range(3)) for k in range(3)]
        history[key].append(list(rows[key]))
Qe1 = [history['e1'][N][0] for N in range(KMAX+1)]

def decompose(vals, NMAX):
    coeffs = []
    for K in range(NMAX + 1):
        rhs = vals[K]
        for k in range(K):
            rhs -= coeffs[k] * delannoy_summand(K, k)
        coeffs.append(rhs / delannoy_summand(K, K))
    return coeffs

f = decompose(Qe1, KMAX)
assert f[1] == F(5749, 3136)
print(f"data ready ({time.time()-t0:.1f}s)\n")

# ----------------------------------------------------------------------
# 1. Denominator structure
# ----------------------------------------------------------------------
def factorize(n, bound=200000):
    """Trial division by primes < bound; returns (factors, cofactor)."""
    fac = {}
    n = abs(n)
    d = 2
    while d < bound and d * d <= n:
        while n % d == 0:
            fac[d] = fac.get(d, 0) + 1
            n //= d
        d += 1 if d == 2 else 2
    return fac, n  # n = unfactored smooth-free cofactor (1 if fully factored)

print("=== Denominator factorization of f(k) (primes < 2e5) ===")
for k in range(0, 17):
    den = f[k].denominator
    fac, cof = factorize(den)
    s = " * ".join(f"{p}^{e}" if e > 1 else f"{p}" for p, e in sorted(fac.items()))
    tail = "" if cof == 1 else f" * [cofactor {cof}]"
    print(f"  k={k:2d}: den = {s if s else 1}{tail}")
print()

print("=== Largest small-prime factor of den(f(k)) vs k ===")
for k in range(2, 33, 2):
    den = f[k].denominator
    fac, cof = factorize(den)
    big = max(fac) if fac else 1
    print(f"  k={k:2d}: largest small prime {big}, cofactor{'=1' if cof==1 else ' big'}"
          f"  (2k+1={2*k+1}, 2k+3={2*k+3}, 2k+5={2*k+5})")
print()

# check: does den(f(k)) divide (2^k C(2k,k))^a * something small?
print("=== u(k) = 2^k C(2k,k) f(k): denominator ===")
for k in range(0, 17):
    u = f[k] * (2**k * comb(2*k, k))
    fac, cof = factorize(u.denominator)
    s = " * ".join(f"{p}^{e}" if e > 1 else f"{p}" for p, e in sorted(fac.items()))
    tail = "" if cof == 1 else f" * [cofactor {cof}]"
    print(f"  k={k:2d}: den(u) = {s if s else 1}{tail}")
print()
sys.stdout.flush()

# ----------------------------------------------------------------------
# 2. Recurrence search on rescalings (mod-p prescreen only: hit/no-hit)
# ----------------------------------------------------------------------
P0 = (1 << 61) - 1

def frac_modp(x, p):
    return (x.numerator % p) * pow(x.denominator % p, p - 2, p) % p

def nullspace_dim_modp(A, ncols, p):
    A = [row[:] for row in A]
    m = len(A)
    prow = 0
    for col in range(ncols):
        piv = next((r for r in range(prow, m) if A[r][col] % p != 0), None)
        if piv is None:
            continue
        A[prow], A[piv] = A[piv], A[prow]
        inv = pow(A[prow][col], p - 2, p)
        A[prow] = [(x * inv) % p for x in A[prow]]
        for r in range(m):
            if r != prow and A[r][col] % p != 0:
                fac = A[r][col]
                A[r] = [(A[r][j] - fac * A[prow][j]) % p for j in range(ncols)]
        prow += 1
        if prow == m:
            break
    return ncols - prow

def scan(seq, name, orders=range(2, 9), degrees=range(1, 25), kmin=5):
    K = len(seq) - 1
    sm = [frac_modp(x, P0) for x in seq]
    grid = sorted(((r, d) for r in orders for d in degrees),
                  key=lambda t: ((t[0]+1)*(t[1]+1), t[0]))
    for r, d in grid:
        nunk = (r + 1) * (d + 1)
        kmax_row = K - r
        if kmax_row - kmin + 1 < nunk + 8:
            continue
        A = []
        for k in range(kmin, kmax_row + 1):
            row = []
            for j in range(r + 1):
                fk = sm[k + j]
                kp = 1
                for _ in range(d + 1):
                    row.append(fk * kp % P0)
                    kp *= k
            A.append(row)
        dim = nullspace_dim_modp(A, nunk, P0)
        if dim > 0:
            print(f"  [{name}] HIT: r={r}, d={d}, nullity={dim}")
            return (r, d)
    print(f"  [{name}] no recurrence (r<=8, d<=24)")
    return None

print("=== Recurrence scan of rescaled sequences (mod p prescreen) ===")
variants = {
    "u = 2^k C(2k,k) f":        [f[k] * (2**k * comb(2*k, k)) for k in range(KMAX+1)],
    "C(2k,k) f":                [f[k] * comb(2*k, k) for k in range(KMAX+1)],
    "4^k f":                    [f[k] * 4**k for k in range(KMAX+1)],
    "(2k+1) 2^k C(2k,k) f":     [f[k] * ((2*k+1) * 2**k * comb(2*k, k)) for k in range(KMAX+1)],
    "(2^k C(2k,k))^2 f":        [f[k] * (2**k * comb(2*k, k))**2 for k in range(KMAX+1)],
    "f / (2^k C(2k,k))":        [f[k] / (2**k * comb(2*k, k)) for k in range(KMAX+1)],
}
for nm, seq in variants.items():
    scan(seq, nm)
    sys.stdout.flush()

print(f"\nTotal {time.time()-t0:.1f}s")
