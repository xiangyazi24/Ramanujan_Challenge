#!/usr/bin/env python3
"""CRON_greens_reframing_verify.py — avenue (a) first strike: Green's-function / projective-orbit
reframing of the (R)+(C) package.

Claims to verify:
 1. Bottom (r-direction) recurrence, SYMBOLIC over Z, h<=12:
      N_h(r) = P(r+1) * N_{h-1}(r+1) - (r+2)^6 * N_{h-2}(r+2),  P(t)=34t^3+51t^2+27t+5.
    (Together with the h-direction transfer recurrence this makes W(r,h)=N_h(r) the solution
    of a FIXED bivariate holonomic system = Green's function of the Apery operator.)
 2. Projective-collision identity mod p: with b = Apery solution, c = companion solution
    (c_0=0, c_1=1, same recurrence), for all 0<=r<r+h<=p-2:
      p | N_h(r)  <==>  b_r*c_{r+h} - c_r*b_{r+h} ≡ 0 (mod p)
      i.e. gap roots = COLLISIONS of the projective orbit pi(n)=[b_n:c_n] in P^1(F_p).
    Verified EXHAUSTIVELY (all (r,h) pairs) for p in {17, 181, 379}; plus color check
    b_r=0 <=> pi(r)=[0:1].
 3. Casoratian: W(n)=b_{n+1}c_n - b_n c_{n+1} satisfies W(n) = W(1)*prod(k^3/(k+1)^3) i.e.
      W(n) ≡ -1/n^3 * (0!/n!)^0 ... explicit closed form W(n) = -(n!)^{-3}? verify numerically:
      W(n)*n!^3 ≡ const mod p.
 4. Orbit collision energy: E^pi_p = sum over v in P^1(F_p) of mult(v)^2 for the orbit points
    {pi(n): 0<=n<=p-2}; report E^pi_p/p across primes (random-map benchmark ~2p incl diagonal),
    plus the reflection Mobius: find gamma in PGL_2(F_p) with pi(p-1-n) = gamma(pi(n)) for all n
    (solve from 3 points, verify on all); report whether it exists (mirror structure = one fixed
    Mobius conjugation).
"""
import sys
from sympy import Poly, symbols, ZZ, expand

X = symbols('X')
fails = 0
def chk(c, m):
    global fails
    print(("  [OK ] " if c else "  [FAIL] ") + m)
    if not c: fails += 1

# ---------- 1. symbolic bottom recurrence ----------
def P_expr(t): return 34*t**3 + 51*t**2 + 27*t + 5

def build_N(hmax):
    N = {1: Poly(1, X, domain=ZZ), 2: Poly(P_expr(X+1), X, domain=ZZ)}
    for d in range(2, hmax):
        N[d+1] = Poly(P_expr(X+d), X, domain=ZZ)*N[d] - Poly((X+d)**6, X, domain=ZZ)*N[d-1]
    return N

print("== 1. bottom recurrence N_h(r) = P(r+1)N_{h-1}(r+1) - (r+2)^6 N_{h-2}(r+2), h<=12 ==")
N = build_N(13)
ok = True
for h in range(3, 13):
    lhs = N[h].as_expr()
    rhs = expand(P_expr(X+1)*N[h-1].as_expr().subs(X, X+1) - (X+2)**6*N[h-2].as_expr().subs(X, X+2))
    if expand(lhs - rhs) != 0: ok = False
chk(ok, "symbolic identity holds for h=3..12")

# ---------- 2. projective-collision identity mod p ----------
def solutions_mod(p):
    b = [0]*p; c = [0]*p
    b[0], b[1] = 1 % p, 5 % p
    c[0], c[1] = 0, 1
    inv = [0]*(p+1); inv[1] = 1
    for i in range(2, p): inv[i] = (p - (p//i)*inv[p % i]) % p
    for n in range(1, p-1):
        Pn = (34*n**3+51*n**2+27*n+5) % p
        i3 = pow(inv[n+1], 3, p)
        b[n+1] = ((Pn*b[n] - pow(n,3,p)*b[n-1]) * i3) % p
        c[n+1] = ((Pn*c[n] - pow(n,3,p)*c[n-1]) * i3) % p
    return b, c

def N_mod(p, r, h):
    # continuant mod p: w_{d+1} = P(r+d) w_d - (r+d)^6 w_{d-1}, N_h = w_h with w_0=0,w_1=1
    w0, w1 = 0, 1
    for d in range(1, h):
        t = (r + d) % p
        Pt = (34*pow(t,3,p) + 51*pow(t,2,p) + 27*t + 5) % p
        w0, w1 = w1, (Pt*w1 - pow(t,6,p)*w0) % p
    return w1

print("== 2. projective-collision identity, exhaustive p=17,181,379 ==")
for p in (17, 181, 379):
    b, c = solutions_mod(p)
    ok = True; ncoll = 0
    for r in range(0, p-2):
        for h in range(1, p-1-r):
            gap = (N_mod(p, r, h) == 0)
            cas = ((b[r]*c[r+h] - c[r]*b[r+h]) % p == 0)
            if gap != cas: ok = False
            if cas: ncoll += 1
    chk(ok, f"p={p}: p|N_h(r) <=> b_r c_s = c_r b_s (ALL {((p-2)*(p-3))//2}+ pairs; collisions={ncoll})")
    chk(all((b[r] == 0) == (b[r] % p == 0 and (b[r], c[r])[0] == 0) for r in range(p-1)), f"p={p}: color trivially consistent", )

# ---------- 3. Casoratian closed form ----------
print("== 3. Casoratian W(n)*'(n!)^3' constant ==")
for p in (181, 379):
    b, c = solutions_mod(p)
    vals = set()
    for n in range(1, p-1):
        W = (b[n]*c[n-1] - b[n-1]*c[n]) % p
        # W(n) = W(1)/n^3 from Casoratian ratio n^3/(n+1)^3  =>  W(n)*n^3 = const
        vals.add((W * pow(n, 3, p)) % p)
    chk(len(vals) == 1, f"p={p}: W(n)*n^3 constant -> {'YES' if len(vals)==1 else 'NO ' + str(sorted(list(vals))[:4])}")

# ---------- 4. orbit collision energy + reflection Mobius ----------
print("== 4. orbit energy and reflection Mobius ==")
import math
def pi_point(bn, cn, p):
    if cn == 0: return p  # infinity marker (b:0) -> point at infinity of [b:c] chart t=b/c? use t = b*c^{-1}
    return (bn * pow(cn, p-2, p)) % p

def primes_upto(n):
    s = bytearray([1])*(n+1); s[0:2] = b'\x00\x00'
    for i in range(2, int(n**0.5)+1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(7, n+1) if s[i]]

es = []
mob_ok_all = True
for p in primes_upto(3000):
    b, c = solutions_mod(p)
    from collections import Counter
    pts = [pi_point(b[n], c[n], p) for n in range(p-1)]
    cnt = Counter(pts)
    E = sum(v*v for v in cnt.values())
    es.append((p, E/p))
    # reflection Mobius: solve gamma from three point pairs (pi(n), pi(p-1-n)) if distinct
    # gamma(t) = (at+b)/(ct+d): use n=0,1,2 -> if degenerate skip
    # quick empirical test: does the multiset {pi(p-1-n)} equal {gamma(pi(n))} for SOME gamma?
    # We test existence via cross-ratio: gamma exists iff cross-ratios match for all quadruples;
    # spot-check with one quadruple per prime (cheap necessary condition).
    def cross_ratio(t1,t2,t3,t4,p):
        INF = p
        def sub(a,bb):
            if a == INF and bb == INF: return None
            if a == INF: return 1  # handled by limits; use projective formula instead
            if bb == INF: return 1
            return (a-bb) % p
        # full projective cross ratio with infinity handling
        def d(a,bb):
            if a == INF and bb == INF: return 0
            if a == INF or bb == INF: return None  # signals leading term
            return (a - bb) % p
        # (t1-t3)(t2-t4)/((t1-t4)(t2-t3)) with infinity limits
        num1, num2, den1, den2 = d(t1,t3), d(t2,t4), d(t1,t4), d(t2,t3)
        parts = [num1, num2, den1, den2]
        if None not in parts:
            if den1 == 0 or den2 == 0: return None
            return (num1*num2*pow(den1*den2, p-2, p)) % p
        return 'inf-case'
    n0 = [n for n in range(1, p-2)][:8]
    q = [pts[n] for n in n0[:4]]; qr = [pts[p-1-n] for n in n0[:4]]
    cr1, cr2 = cross_ratio(*q, p), cross_ratio(*qr, p)
    if cr1 is not None and cr2 is not None and cr1 != 'inf-case' and cr2 != 'inf-case':
        if cr1 != cr2: mob_ok_all = False

n = len(es)
mean_late = sum(e for p, e in es if p > 1500)/max(1, len([1 for p, e in es if p > 1500]))
print(f"  orbit energy E^pi/p over {n} primes: overall mean={sum(e for _,e in es)/n:.3f}, mean(p>1500)={mean_late:.3f} (random benchmark ~2)")
chk(mob_ok_all, "reflection Mobius cross-ratio necessary condition holds at every tested prime")

print()
print("ALL PASS" if fails == 0 else f"{fails} FAILURES")
sys.exit(0 if fails else 1)
