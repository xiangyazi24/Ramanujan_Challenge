#!/usr/bin/env python3
"""Go/no-go: measure N_coinc(H) at H ~ p^{2/3} (and other scales).
K_emp = (N_coinc - #S^2/p)/#S  must stay bounded for [COINC] to be true.
Also: verify (D1) Delta identity on a subsample; report n(0)=sum C_h vs H+sqrt(pH)."""
import sys, time
from collections import Counter

def run(p, exps=(0.25, 0.33, 0.5, 0.66)):
    # b, c mod p on regular window n=0..p-1
    b = [1, 5 % p]; c = [0, 1]
    inv = [0]*p
    inv[1] = 1
    for i in range(2, p): inv[i] = (p - (p//i)*inv[p % i]) % p
    for n in range(1, p-1):
        i3 = pow(inv[(n+1) % p], 3, p)
        Pn = ((2*n+1)*(17*n*n+17*n+5)) % p
        b.append(i3*(Pn*b[n] - pow(n,3,p)*b[n-1]) % p)
        c.append(i3*(Pn*c[n] - pow(n,3,p)*c[n-1]) % p)
    M = p-2
    out = []
    for e in exps:
        H = max(2, round(p**e))
        if H > M-1: H = M-1
        t0 = time.time()
        vals = Counter()
        nS = 0
        for h in range(1, H+1):
            for r in range(1, M-h+1):
                d = (b[r]*c[r+h] - b[r+h]*c[r]) % p
                vals[d] += 1
                nS += 1
        Nc = sum(v*v for v in vals.values())
        n0 = vals.get(0, 0)
        Kemp = (Nc - nS*nS/p)/nS
        # zero-zero share of the excess
        zz_excess = (n0*n0 - (nS/p)**2*0)  # raw n0^2 (its own centered excess ~ n0^2 since nS/p tiny per-cell? no: cell budget nS/p)
        n0_centered_sq = (n0 - nS/p)**2
        out.append((H, e, nS, Nc, n0, Kemp, n0_centered_sq/nS, time.time()-t0))
    # D1 spot-check
    bad = 0; tot = 0
    N_prev, N_cur = None, None
    import random
    random.seed(1)
    for _ in range(300):
        h = random.randint(1, 40); r = random.randint(1, M-h)
        if any((r+j) % p == 0 for j in range(1, h+1)): continue
        # N_h(r) via recurrence
        Nm1, N0 = 0, 1  # N_0 := 0? convention: N_1=1, N_2=P(r+1)
        N1 = 1; N2 = ((2*(r+1)+1)*(17*(r+1)**2+17*(r+1)+5)) % p
        if h == 1: Nh = N1
        else:
            a_, b_ = N1, N2
            for m in range(2, h):
                Pm = ((2*(r+m)+1)*(17*(r+m)**2+17*(r+m)+5)) % p
                a_, b_ = b_, (Pm*b_ - pow(r+m, 6, p)*a_) % p
            Nh = b_
        pr = 1
        for j in range(1, h+1): pr = pr*pow(r+j, 3, p) % p
        lhs = (b[r]*c[r+h] - b[r+h]*c[r]) % p
        rhs = Nh * pow(pr, p-2, p) % p
        tot += 1
        if lhs != rhs: bad += 1
    print(f"p={p}: D1 check {tot-bad}/{tot} pass")
    for H, e, nS, Nc, n0, Kemp, n0c, dt in out:
        import math
        print(f"  H=p^{e:.2f}={H:5d} #S={nS:9d} N_coinc={Nc:12d} K_emp={Kemp:8.3f} "
              f"n0={n0:7d} (H+sqrt(pH)={H+math.isqrt(p*H):6d}) n0^2/#S={n0c:8.3f} [{dt:.1f}s]")

for p in [499, 997, 1999, 4001, 7919]:
    run(p)
