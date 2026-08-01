#!/usr/bin/env python3
"""Spot-check cron's s_p = truncation-of-fixed-series claim at p=13 (chi=-1) and p=29 (chi=+1).
F(x) = sum b_k x^k; sqrt(F) and sqrt(F/q), q = 1-34x+x^2, as power series mod p.
A_p(x) = sum_{k<p} b_k x^k mod p; compare its square-root factor with the truncation."""
import sympy as sp
x = sp.symbols('x')

def apery(N):
    b = [1, 5]
    for n in range(1, N):
        num = (2*n+1)*(17*n*n+17*n+5)*b[n] - n**3*b[n-1]
        qq, r = divmod(num, (n+1)**3); assert r == 0
        b.append(qq)
    return b

def series_sqrt(coeffs, N, p):
    # sqrt of power series with c0=1 mod p
    s = [1] + [0]*(N-1)
    inv2 = pow(2, p-2, p)
    for n in range(1, N):
        acc = coeffs[n] if n < len(coeffs) else 0
        t = sum(s[i]*s[n-i] for i in range(1, n)) % p
        s[n] = (acc - t) % p * inv2 % p
    return s

def series_div(a, bden, N, p):
    # a/bden as series mod p
    binv0 = pow(bden[0], p-2, p)
    c = [0]*N
    for n in range(N):
        t = a[n] if n < len(a) else 0
        t = (t - sum(bden[i]*c[n-i] for i in range(1, min(n, len(bden)-1)+1))) % p
        c[n] = t*binv0 % p
    return c

for p, chi in ((13, -1), (29, +1)):
    b = apery(p+2)
    bm = [v % p for v in b[:p]]
    Ap = sp.Poly(sum(bm[k]*x**k for k in range(p)), x, modulus=p)
    fac = sp.factor_list(Ap.as_expr(), modulus=p)
    # extract square part: multiply factors with even exponent halves
    sq = sp.Poly(1, x, modulus=p); rest = sp.Poly(1, x, modulus=p)
    for base, e in fac[1]:
        pb = sp.Poly(base, x, modulus=p)
        if e >= 2:
            sq = sq * pb**(e//2)
        if e % 2:
            rest = rest * pb
    print(f"p={p} (chi={chi}): deg A_p={Ap.degree()}, sqrt-part degree={sq.degree()}, odd-part={sp.factor(rest.as_expr(), modulus=p)}")
    # compare sq (up to scalar) with truncation of sqrt(F) or sqrt(F/q)
    if chi == +1:
        target = series_sqrt(bm, (p-1)//2 + 1, p)
    else:
        q = [1, (-34) % p, 1]
        Fq = series_div(bm, q, (p-1)//2 + 1, p)
        target = series_sqrt(Fq, (p-3)//2 + 1, p)
    tr = sp.Poly(sum(target[k]*x**k for k in range(len(target))), x, modulus=p)
    # compare projectively: sq ~ c*tr?
    if sq.degree() == tr.degree():
        c = sq.all_coeffs()[0] * pow(tr.all_coeffs()[0], p-2, p) % p
        match = all((a - c*bb) % p == 0 for a, bb in zip(sq.all_coeffs(), tr.all_coeffs()))
    else:
        match = False
    print(f"   truncated fixed-series matches square-root factor (up to scalar): {match}")
