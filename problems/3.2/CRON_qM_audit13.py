#!/usr/bin/env python3
"""CRON_qM_audit13.py — Q6509 AUDIT-1 与 AUDIT-3 机器审计（零失败=定理级入账）.

AUDIT-1 depth3_lucas_companion: p∈{13,17,29}, 0≤a≤3(a≥1), 0≤v<p:
  (i)  p³·c_{ap+v} ≡ c_a·b_v (mod p)   [c_n 为伴随解, 有理数; p³c 应 p-整]
  (ii) W_{ap+v}/p^{3(a−1)} ≡ (−1)^a·W_a·(v!)³·b_v (mod p),  W_n=(n!)³c_n ∈ Z
AUDIT-3 bernoulli_scalar: κ_p := ((p³c_p − 1)/p³) mod p ≡ −(2/3)·B_{p−3} (mod p);
  且 (W_p − ((p−1)!)³)/p³ ≡ −κ_p (mod p).
"""
from fractions import Fraction
from math import factorial
from sympy import bernoulli, Rational

def solve_pair(N):
    b = [Fraction(1), Fraction(5)]
    c = [Fraction(0), Fraction(1)]
    for n in range(2, N + 1):
        A = 34*n**3 - 51*n**2 + 27*n - 5
        B = (n - 1)**3
        b.append((A*b[n-1] - B*b[n-2]) / n**3)
        c.append((A*c[n-1] - B*c[n-2]) / n**3)
    return b, c

def vp(fr, p):
    if fr == 0: return 10**9
    n, d = fr.numerator, fr.denominator
    v = 0
    while n % p == 0: n //= p; v += 1
    while d % p == 0: d //= p; v -= 1
    return v

def res_mod(fr, p):
    n, d = fr.numerator, fr.denominator
    assert d % p != 0, "non-p-integral"
    return n * pow(d, -1, p) % p

fails = 0
def chk(cond, msg):
    global fails
    if not cond:
        fails += 1
        print(f"  [FAIL] {msg}", flush=True)

for p in (13, 17, 29):
    N = 3*p + p
    b, c = solve_pair(N)
    W = [c[n] * Fraction(factorial(n))**3 for n in range(N+1)]
    assert all(w.denominator == 1 for w in W), "W not integral"
    n_checked = 0
    for a in (1, 2, 3):
        ca_res = res_mod(c[a], p) if a >= 1 else None
        for v in range(p):
            n = a*p + v
            lhs = c[n] * p**3
            bv = res_mod(b[v], p)
            # (i) p³c_{ap+v} ≡ c_a·b_v mod p (两边可同为 0 当 p|b_v)
            if vp(lhs, p) < 0:
                chk(False, f"p={p} a={a} v={v}: p³c not p-integral (vp={vp(lhs,p)})")
                continue
            lhs_res = res_mod(lhs, p) if vp(lhs, p) == 0 else 0
            rhs_res = ca_res * bv % p
            chk(lhs_res == rhs_res, f"p={p} a={a} v={v}: (i) {lhs_res} != {rhs_res}")
            # (ii) W 版
            Wn = W[n].numerator
            q, r3 = divmod(Wn, p**(3*(a-1)))
            chk(r3 == 0, f"p={p} a={a} v={v}: (ii) valuation short")
            rhs2 = ((-1)**a * W[a].numerator * factorial(v)**3 * b[v].numerator
                    * pow(b[v].denominator, -1, p)) % p
            chk(q % p == rhs2 % p, f"p={p} a={a} v={v}: (ii) {q%p} != {rhs2%p}")
            n_checked += 1
    # AUDIT-3
    kap_fr = (c[p]*p**3 - 1) / p**3
    chk(vp(kap_fr, p) >= 0, f"p={p}: kappa not integral")
    kappa = res_mod(kap_fr, p)
    B = Rational(bernoulli(p - 3))
    Bres = int(B.p) * pow(int(B.q), -1, p) % p
    target = (-Fraction(2, 3))
    tres = (p - 2) * pow(3, -1, p) % p * Bres % p
    chk(kappa == tres, f"p={p}: AUDIT-3 kappa={kappa} != -(2/3)B_(p-3)={tres}")
    w_side = (W[p].numerator - factorial(p-1)**3)
    q3, r3 = divmod(w_side, p**3)
    chk(r3 == 0 and q3 % p == (p - kappa) % p, f"p={p}: W-side kappa mismatch")
    print(f"p={p}: AUDIT-1 {n_checked} 组 (a,v) 全过 + AUDIT-3 κ_p={kappa} = −(2/3)B_{{p−3}} mod p ✓"
          if fails == 0 else f"p={p}: FAILURES so far {fails}", flush=True)

print(f"\nRESULT: {'ALL PASS' if fails == 0 else f'{fails} FAILURES'}", flush=True)
