#!/usr/bin/env python3
"""End-to-end check: every fiber element with two forward lags roots Psi_{d1,d2} (or B_{d1}/B_{d2} for type-II legs), p=101,199."""
def apery_mod(p):
    b=[1,5]
    for n in range(1,p):
        b.append(((2*n+1)*(17*n*n+17*n+5)%p*b[n]-pow(n,3,p)*b[n-1])*pow(pow(n+1,3,p),-1,p)%p if n+1<p else 0)
    # last step n+1=p singular: recompute exact for safety
    if len(b)<p: b.append(0)
    return b
def apery_exact_mod(p):
    b=[1,5]
    for n in range(1,p+1):
        num=(2*n+1)*(17*n*n+17*n+5)*b[n]-n**3*b[n-1]
        q,rm=divmod(num,(n+1)**3); assert rm==0
        b.append(q)
    return [x%p for x in b[:p]]
def AB(d, rr, p):
    def co(x):
        d3=pow((x+1)%p,3,p)
        return (2*x+1)*(17*x*x+17*x+5)%p*pow(d3,-1,p)%p, (-(x**3))%p*pow(d3,-1,p)%p
    Ap,Bp=1,0
    a1,b1=co(rr%p); Ac,Bc=a1,b1
    for k in range(1,d):
        ak,bk=co((rr+k)%p)
        Ac,Ap=(ak*Ac+bk*Ap)%p,Ac
        Bc,Bp=(ak*Bc+bk*Bp)%p,Bc
    return Ac,Bc
for p in (101,199):
    bm=apery_exact_mod(p)
    from collections import defaultdict
    fib=defaultdict(list)
    for r in range(p): fib[bm[r]].append(r)
    total=checked=okc=0; fails=[]
    for c,rs in fib.items():
        if c==0 or len(rs)<3: continue
        for i in range(len(rs)-2):
            r,r1,r2=rs[i],rs[i+1],rs[i+2]
            d1,d2=r1-r,r2-r
            if r2>p-2 or r<1: continue
            total+=1
            A1,B1=AB(d1,r,p); A2,B2=AB(d2,r,p)
            if B1%p==0:
                ok=(A1%p==1)   # type II at d1
            elif B2%p==0:
                ok=(A2%p==1)
            else:
                psi=((1-A1)*B2-(1-A2)*B1)%p
                ok=(psi==0)    # r roots Psi
            checked+=1; okc+=ok
            if not ok: fails.append((p,c,r,d1,d2))
    print(f"p={p}: elements-with-two-lags checked {checked}, mechanism verified {okc}, FAILS {len(fails)} {fails[:5]}")
