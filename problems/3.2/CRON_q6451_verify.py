#!/usr/bin/env python3
"""CRON_q6451_verify.py — machine audit of ANSWER Q6451 (GARQI-1' gap-polynomial route audit).

Independently reconstructs the gap continuant and verifies the answer's spine:
  Derivation (independent): clearing denominators of the Apery recurrence with
  w_d = b_{r+d} * ((r+1)...(r+d))^3 gives w_{d+1} = A_d w_d - B_d w_{d-1},
  A_d = P(r+d) = 34(r+d)^3+51(r+d)^2+27(r+d)+5,  B_d = (r+d)^6.
  With b_r = 0: w_h = N_h(r) * w_1, N_1 = 1, N_2 = A_1, N_{d+1} = A_d N_d - B_d N_{d-1}.
  Hence (b_r ≡ 0 and b_{r+h} ≡ 0 mod p, all indices < p-1) ==> p | N_h(r).

Checks:
 1. deg N_h = 3(h-1) for h <= 10 (claimed exact degree).
 2. Reflection: N_h(-h-1-X) = ±N_h(X) (even for odd h, odd for even h after centering);
    for even h the forced linear factor (2X+h+1) divides N_h over Z.
 3. Real-data double zeros: for p in {17, 181, 379} compute Z_p, and for EVERY ordered
    pair r < r+h in Z_p verify p | N_h(r); for mirror pairs with even h verify
    2r+h+1 ≡ 0 mod p (forced factor = p exactly).
 4. Integer vanishing locus empty on a strip: N_h(r) != 0 for all integers
    r in [-h-5, 25], h <= 10; positivity for r >= 0.
 5. Quadratic degree sum: sum_{h<=H} deg N_h = 3*H*(H-1)/2 (the H^2 obstruction).
"""
import sys
from sympy import Poly, symbols, ZZ, factor, div

X = symbols('X')

def P_at(e):  # P(X+e) as sympy expr
    t = X + e
    return 34*t**3 + 51*t**2 + 27*t + 5

def build_N(hmax):
    N = {1: Poly(1, X, domain=ZZ)}
    if hmax >= 2:
        N[2] = Poly(P_at(1), X, domain=ZZ)
    for d in range(2, hmax):
        A = Poly(P_at(d), X, domain=ZZ)
        B = Poly((X + d)**6, X, domain=ZZ)
        N[d+1] = A * N[d] - B * N[d-1]
    return N

def apery_mod(p):
    b = [0]*p
    b[0] = 1 % p; b[1] = 5 % p
    inv = [0]*(p+1); inv[1] = 1
    for i in range(2, p+1):
        inv[i] = (p - (p//i)*inv[p % i]) % p if i < p else 0
    for n in range(1, p-1):
        Pn = (34*n**3 + 51*n**2 + 27*n + 5) % p
        iv = inv[(n+1) % p]
        b[n+1] = ((Pn*b[n] - pow(n, 3, p)*b[n-1]) * pow(iv, 3, p)) % p
    return [r for r in range(p-1) if b[r] == 0]

fails = 0
def chk(c, m):
    global fails
    print(("  [OK ] " if c else "  [FAIL] ") + m)
    if not c: fails += 1

HMAX = 12
N = build_N(HMAX)

print("== 1. exact degree 3(h-1) ==")
chk(all(N[h].degree() == 3*(h-1) for h in range(1, HMAX+1)),
    f"deg N_h = 3(h-1) for h=1..{HMAX}")

print("== 2. reflection + forced mirror factor ==")
ok_refl = True
for h in range(2, HMAX+1):
    refl = Poly(N[h].as_expr().subs(X, -h-1-X), X, domain=ZZ)
    sgn = 1 if (h % 2 == 1) else -1
    if refl != sgn * N[h]: ok_refl = False
chk(ok_refl, "N_h(-h-1-X) = (-1)^{h-1} N_h(X) for h=2..%d" % HMAX)
ok_fac = True
for h in range(2, HMAX+1, 2):
    q, r = div(N[h], Poly(2*X + h + 1, X, domain=ZZ), domain='QQ')
    if not r.is_zero: ok_fac = False
chk(ok_fac, "even h: (2X+h+1) | N_h over Q (forced mirror factor)")

print("== 3. real-data double zeros ==")
for p in (17, 181, 379):
    Z = apery_mod(p)
    pairs = [(r, s) for i, r in enumerate(Z) for s in Z[i+1:]]
    allok = True; mirror_ok = True; tested = 0
    for r, s in pairs:
        h = s - r
        if h > HMAX:
            # rebuild lazily for large gaps
            NN = build_N(h)
            val = NN[h].eval(r)
        else:
            val = N[h].eval(r)
        if val % p != 0: allok = False
        tested += 1
        if s == p - 1 - r:  # mirror pair
            if h % 2 != 0 or (2*r + h + 1) % p != 0: mirror_ok = False
    chk(allok, f"p={p}: p | N_h(r) for ALL {tested} zero pairs (Z_p={Z})")
    chk(mirror_ok, f"p={p}: every mirror pair has even h and forced factor 2r+h+1 ≡ 0 (mod p)")

print("== 4. integer vanishing locus empty on strip ==")
ok_nz = True; ok_pos = True
for h in range(1, HMAX+1):
    for r in range(-h-5, 26):
        v = N[h].eval(r)
        if v == 0: ok_nz = False
        if r >= 0 and v <= 0: ok_pos = False
chk(ok_nz, f"N_h(r) != 0 for h<=%d, r in [-h-5, 25]" % HMAX)
chk(ok_pos, "N_h(r) > 0 for r >= 0")

print("== 5. quadratic degree sum ==")
H = HMAX
chk(sum(N[h].degree() for h in range(1, H+1)) == 3*H*(H-1)//2,
    f"sum_h<=H deg N_h = 3H(H-1)/2 (quadratic obstruction), H={H}")

print()
print("ALL PASS" if fails == 0 else f"{fails} FAILURES")
sys.exit(0 if fails else 1)
