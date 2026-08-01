#!/usr/bin/env python3
"""max_t |S_h(t)|/sqrt(p) vs candidate constants 4h+1 / 4h-1 / trivial 6h-3."""
import cmath, math
def run(p, hs):
    b=[1,5%p]; c=[0,1]
    inv=[0]*p; inv[1]=1
    for i in range(2,p): inv[i]=(p-(p//i)*inv[p%i])%p
    for n in range(1,p-1):
        i3=pow(inv[(n+1)%p],3,p); Pn=((2*n+1)*(17*n*n+17*n+5))%p
        b.append(i3*(Pn*b[n]-pow(n,3,p)*b[n-1])%p)
        c.append(i3*(Pn*c[n]-pow(n,3,p)*c[n-1])%p)
    M=p-2
    for h in hs:
        # complete sum over r in F_p minus poles: delta = N_h(r)/prod (r+j)^3; use full F_p domain
        # compute N_h(r) via recurrence for each r (h small: O(ph))
        vals=[]
        for r in range(p):
            ok=True; pr=1
            for j in range(1,h+1):
                t=(r+j)%p
                if t==0: ok=False; break
                pr=pr*pow(t,3,p)%p
            if not ok: continue
            if h==1: Nh=1
            else:
                a_,b_=1,((2*(r+1)+1)*(17*(r+1)**2+17*(r+1)+5))%p
                for m in range(2,h):
                    Pm=((2*(r+m)+1)*(17*(r+m)**2+17*(r+m)+5))%p
                    a_,b_=b_,(Pm*b_-pow(r+m,6,p)*a_)%p
                Nh=b_
            vals.append(Nh*pow(pr,p-2,p)%p)
        hist=[0]*p
        for v in vals: hist[v]+=1
        best=0
        w=cmath.exp(2j*cmath.pi/p)
        # DFT via direct eval per t (O(p^2) small p ok)
        roots=[cmath.exp(2j*cmath.pi*k/p) for k in range(p)]
        for t in range(1,p):
            s=sum(hist[a]*roots[(t*a)%p] for a in range(p) if hist[a])
            m=abs(s)
            if m>best: best=m
        print(f"p={p} h={h}: max|S|/sqrt(p) = {best/math.sqrt(p):.2f}  vs 4h+1={4*h+1} 4h-1={4*h-1} deg-triv={6*h-3}")
for p in [199, 499]:
    run(p,[2,3,5,8])
