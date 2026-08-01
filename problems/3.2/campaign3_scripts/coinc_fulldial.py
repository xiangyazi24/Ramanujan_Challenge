#!/usr/bin/env python3
"""K_emp across the full beta dial, up to H=M-1 (E1 scale)."""
import math
from collections import Counter
def run(p, exps):
    b=[1,5%p]; c=[0,1]
    inv=[0]*p; inv[1]=1
    for i in range(2,p): inv[i]=(p-(p//i)*inv[p%i])%p
    for n in range(1,p-1):
        i3=pow(inv[(n+1)%p],3,p); Pn=((2*n+1)*(17*n*n+17*n+5))%p
        b.append(i3*(Pn*b[n]-pow(n,3,p)*b[n-1])%p)
        c.append(i3*(Pn*c[n]-pow(n,3,p)*c[n-1])%p)
    M=p-2
    for e in exps:
        H=M-1 if e>=1 else min(M-1,round(p**e))
        vals=Counter(); nS=0
        for h in range(1,H+1):
            bh=b[h:]; ch=c[h:]
            for r in range(1,M-h+1):
                vals[(b[r]*c[r+h]-b[r+h]*c[r])%p]+=1
            nS+=M-h
        Nc=sum(v*v for v in vals.values()); n0=vals.get(0,0)
        K=(Nc-nS*nS/p)/nS
        zz=n0*(n0-1)/nS
        print(f"p={p} H={'M-1' if e>=1 else f'p^{e}'}={H:5d} #S={nS:9d} K_emp={K:7.3f} n0={n0:6d} zz/#S={zz:6.3f} n0 vs H+sqrt(pH): {n0}/{H+math.isqrt(p*H)}")
for p in [1999, 4001]:
    run(p,[0.75,0.9,1.0])
