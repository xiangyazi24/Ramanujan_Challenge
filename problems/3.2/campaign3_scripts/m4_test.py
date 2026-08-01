#!/usr/bin/env python3
"""Empirical 4th moment M4 = (1/p^3) sum_{t!=0} |S_h(t)|^4 (normalized so Sp-Haar -> 3).
S_h over complete regular domain."""
import cmath, math
def run(p, hs):
    for h in hs:
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
        roots=[cmath.exp(2j*cmath.pi*k/p) for k in range(p)]
        s4=0.0
        for t in range(1,p):
            s=sum(hist[a]*roots[(t*a)%p] for a in range(p) if hist[a])
            s4+=abs(s)**4
        # normalization: for pure weight-1 traces, E|S|^2 ~ p * rank-variance; Haar-Sp M4: E|Tr|^4 = 3
        # normalize |S|^2/p as |Tr theta|^2; M4_emp = mean over t of |S/sqrt p|^4
        m4 = s4/((p-1)*p*p)
        print(f"p={p} h={h}: M4_emp = {m4:.3f}  (Sp-Haar target 3, O-Haar 3 too at 2nd... GL 2)")
for p in [499, 997]:
    run(p,[2,3,4,6,8])
