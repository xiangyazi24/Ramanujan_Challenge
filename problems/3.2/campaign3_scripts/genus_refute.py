#!/usr/bin/env python3
"""Two checks:
(1) Verify N_h(-j) = (-1)^(j-1) * b_{j-1} * b_{h-j} * ((j-1)!)^3 * ((h-j)!)^3 over Q (exact ints).
(2) Point-count fluctuations of X_{2,3}: N_2(r)*D_3(r') = N_3(r')*D_2(r) over many primes.
    max |#X - p| / sqrt(p) discriminates genus ~19 (my RH computation) vs low genus ~5."""
from math import isqrt, factorial

# exact N_h over Z via recurrence in r (symbolic via polynomial coeff lists)
def polmul(a,b):
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): c[i+j]+=x*y
    return c
def poladd(a,b):
    n=max(len(a),len(b)); return [(a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0) for i in range(n)]
def polscale(a,s): return [x*s for x in a]
def poleval(a,x): 
    v=0
    for c in reversed(a): v=v*x+c
    return v
def P_shift(m):  # P(r+m) as poly in r
    # P(u)=34u^3+51u^2+27u+5, u=r+m
    u=[m,1]
    u2=polmul(u,u); u3=polmul(u2,u)
    return poladd(poladd(polscale(u3,34),polscale(u2,51)),poladd(polscale(u,27),[5]))
def pow6_shift(m):
    u=[m,1]; u2=polmul(u,u); u3=polmul(u2,u); return polmul(u3,u3)

def N_polys(hmax):
    N={1:[1],2:P_shift(1)}
    for h in range(2,hmax):
        N[h+1]=poladd(polmul(P_shift(h),N[h]),polscale(polmul(pow6_shift(h),N[h-1]),-1))
    return N

# Apery numbers exact
def apery(n):
    b=[1,5]
    for k in range(1,n):
        b.append(( (34*k**3+51*k**2+27*k+5)*b[k] - k**3*b[k-1] )//(k+1)**3)
    return b

hmax=10
N=N_polys(hmax)
b=apery(hmax+2)
ok=fail=0
for h in range(2,hmax+1):
    for j in range(1,h+1):
        lhs=poleval(N[h],-j)
        rhs=(-1)**(j-1)*b[j-1]*b[h-j]*(factorial(j-1)**3)*(factorial(h-j)**3)
        if lhs==rhs: ok+=1
        else:
            fail+=1
            if fail<5: print(f"FAIL h={h} j={j}: {lhs} vs {rhs}")
print(f"N_h(-j) Apery-product formula: {ok} pass, {fail} fail (h<=10)")

# (2) point counts of X_{2,3}
import sympy
primes=[q for q in range(300,4000) if sympy.isprime(q)][:250]
N2=N[2]; N3=N[3]
devs=[]
for p in primes:
    n2=[c%p for c in N2]; n3=[c%p for c in N3]
    # D_2(r)=((r+1)(r+2))^3, D_3(r')=((r'+1)(r'+2)(r'+3))^3
    cnt=0
    # tabulate f2(r)=N2(r)/D2(r) for all r with D2!=0; f3 likewise; count pairs equal value via value buckets
    from collections import Counter
    c2=Counter(); c3=Counter()
    for r in range(p):
        d=( (r+1)*(r+2) )%p
        if d==0: continue
        c2[ poleval(n2,r)*pow(pow(d,3,p),p-2,p)%p ]+=1
    for r in range(p):
        d=( (r+1)*(r+2)*(r+3) )%p
        if d==0: continue
        c3[ poleval(n3,r)*pow(pow(d,3,p),p-2,p)%p ]+=1
    cnt=sum(c2[a]*c3[a] for a in c2 if a in c3)
    devs.append((cnt-p)/isqrt(p))
import statistics
print(f"X_2,3 affine point-count deviations over {len(primes)} primes:")
print(f"  mean={statistics.mean(devs):.2f} std={statistics.stdev(devs):.2f} max|dev|={max(abs(d) for d in devs):.2f}")
print("  low-genus(<=5) predicts max|dev|<=2*5=10; my RH genus ~19 predicts larger tail")
