#!/usr/bin/env python3
"""Verify [THM-BASE-RETURN-2/3]: d_D(r) <= min_Y ((3/2)Y(Y-1) + D//(Y+1)) for all r."""
import math
def run(p, D):
    b=[1,5%p]; c=[0,1]
    inv=[0]*p; inv[1]=1
    for i in range(2,p): inv[i]=(p-(p//i)*inv[p%i])%p
    for n in range(1,p-1):
        i3=pow(inv[(n+1)%p],3,p); Pn=((2*n+1)*(17*n*n+17*n+5))%p
        b.append(i3*(Pn*b[n]-pow(n,3,p)*b[n-1])%p)
        c.append(i3*(Pn*c[n]-pow(n,3,p)*c[n-1])%p)
    M=p-2
    bound=min(((3*Y*(Y-1))//2 + D//(Y+1)) for Y in range(2,D+1))
    worst=0; viol=0
    for r in range(1,M-1):
        Dr=min(D, M-r)
        cnt=0
        for d in range(1,Dr+1):
            if (b[r]*c[r+d]-b[r+d]*c[r])%p==0: cnt+=1
        if cnt>worst: worst=cnt
        bnd_r=min(((3*Y*(Y-1))//2 + Dr//(Y+1)) for Y in range(2,max(3,Dr+1)))
        if cnt>bnd_r: viol+=1
    print(f"p={p} D={D}: max_r d_D(r)={worst}, theorem bound={bound}, violations={viol}")
for p,D in [(499,60),(997,95),(1999,150)]:
    run(p,D)
