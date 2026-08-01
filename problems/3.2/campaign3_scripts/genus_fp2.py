#!/usr/bin/env python3
"""Count X_{2,3} points over F_{p^2}: dev2=(cnt-p^2)/p bounded by 2g. g=19 -> typical ~6, g<=2 -> <=4."""
import sympy
from collections import Counter
def fp2_test(p):
    d=2
    while pow(d,(p-1)//2,p)!=p-1: d+=1
    # elements (a,b) = a+b*sqrt(d); arithmetic
    def mul(x,y): return ((x[0]*y[0]+d*x[1]*y[1])%p,(x[0]*y[1]+x[1]*y[0])%p)
    def add(x,y): return ((x[0]+y[0])%p,(x[1]+y[1])%p)
    def scal(k,x): return ((k*x[0])%p,(k*x[1])%p)
    def inv(x):
        n=(x[0]*x[0]-d*x[1]*x[1])%p
        ni=pow(n,p-2,p)
        return ((x[0]*ni)%p,((-x[1])*ni)%p)
    def P(u): return add(add(scal(34,mul(mul(u,u),u)),scal(51,mul(u,u))),add(scal(27,u),(5,0)))
    ONE=(1,0)
    c2=Counter(); c3=Counter()
    for a in range(p):
        for bb in range(p):
            r=(a,bb)
            r1=add(r,ONE); r2=add(r,(2,0)); r3=add(r,(3,0))
            # N2 = P(r+1); D2 = ((r+1)(r+2))^3
            d2=mul(r1,r2)
            if d2!=(0,0):
                v=mul(P(r1),inv(mul(mul(d2,d2),d2)))
                c2[v]+=1
            # N3 = P(r+2)P(r+1) - (r+2)^6 ; D3=((r1)(r2)(r3))^3
            d3=mul(mul(r1,r2),r3)
            if d3!=(0,0):
                r2sq=mul(r2,r2); r2c=mul(r2sq,r2)
                n3=add(mul(P(add(r,(2,0))),P(r1)),scal(p-1,mul(r2c,r2c)))
                v=mul(n3,inv(mul(mul(d3,d3),d3)))
                c3[v]+=1
    cnt=sum(c2[a]*c3[a] for a in c2 if a in c3)
    return (cnt-p*p)/p
devs=[]
for p in [q for q in range(53,140) if sympy.isprime(q)]:
    dv=fp2_test(p); devs.append(dv)
    print(f"p={p}: dev2={dv:+.2f}")
import statistics
print(f"mean={statistics.mean(devs):.2f} std={statistics.stdev(devs):.2f} max|.|={max(abs(x) for x in devs):.2f}  (2g bound: g=19->38, g=2->4)")
