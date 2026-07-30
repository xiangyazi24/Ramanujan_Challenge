#!/usr/bin/env python3
"""Q5716 finite-field audit of the exact Legendre Green/CD reduction."""
from math import comb
MODS=(1000000007,1000000009); MS=(199,271,299,320,754)
def trim(a):
    while len(a)>1 and a[-1]==0:a.pop()
    return a
def add(a,b,p,sgn=1):
    n=max(len(a),len(b)); c=[0]*n
    for i in range(n): c[i]=((a[i] if i<len(a) else 0)+sgn*(b[i] if i<len(b) else 0))%p
    return trim(c)
def scale(a,c,p): return trim([(c*x)%p for x in a])
def mul(a,b,p):
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        if x:
            for j,y in enumerate(b): c[i+j]=(c[i+j]+x*y)%p
    return trim(c)
def mul_1pz(a,c,p):
    o=[0]*(len(a)+1)
    for i,x in enumerate(a): o[i]=(o[i]+c*x)%p; o[i+1]=(o[i+1]+c*x)%p
    return trim(o)
def mul_1mz2(a,c,p):
    o=[0]*(len(a)+2)
    for i,x in enumerate(a):
        o[i]=(o[i]+c*x)%p; o[i+1]=(o[i+1]-2*c*x)%p; o[i+2]=(o[i+2]+c*x)%p
    return trim([x%p for x in o])
def divmod_poly(a,b,p):
    a=a[:]; trim(a); trim(b); inv=pow(b[-1],p-2,p)
    if len(a)<len(b): return [0],a
    q=[0]*(len(a)-len(b)+1)
    while len(a)>=len(b) and not(len(a)==1 and a[0]==0):
        d=len(a)-len(b); c=a[-1]*inv%p; q[d]=c
        for j,y in enumerate(b): a[d+j]=(a[d+j]-c*y)%p
        trim(a)
    return trim(q),trim(a)
def gcd_poly(a,b,p):
    a=trim(a[:]); b=trim(b[:])
    while not(len(b)==1 and b[0]==0): _,r=divmod_poly(a,b,p); a,b=b,r
    return scale(a,pow(a[-1],p-2,p),p)
def J_list(N,p):
    J=[[1],[1,1]]
    for n in range(1,N):
        rhs=add(mul_1pz(J[n],2*n+1,p),mul_1mz2(J[n-1],n,p),p,-1)
        J.append(scale(rhs,pow(n+1,p-2,p),p))
    return J[:N+1]
def audit(M,p):
    tm2=[1]; tm1=scale([1,1],-(2*M+1),p); em1=[1]
    for k in range(1,M+1):
        th=add(mul_1pz(tm1,-(2*M+2*k+1),p),mul_1mz2(tm2,(M+k)**2,p),p,-1)
        w=comb(M,k)%p; w=w*w%p
        eta=add(scale(tm1,w,p),scale(em1,M+k,p),p,-1)
        tm2,tm1,em1=tm1,th,eta
    E=tm1; NM=em1
    pp2=[1]; pp1=scale([1,1],-(4*M+1),p); zp1=[1]
    for k in range(M-1,-1,-1):
        ph=add(mul_1pz(pp1,-(2*M+2*k+1),p),mul_1mz2(pp2,(M+k+1)**2,p),p,-1)
        w=comb(M,k)%p; w=w*w%p
        ze=add(scale(pp1,w,p),mul_1mz2(zp1,M+k+1,p),p,-1)
        pp2,pp1,zp1=pp1,ph,ze
    assert pp1==E; N0=zp1
    J=J_list(2*M+1,p); S=[0]
    for k in range(M+1):
        w=comb(M,k)%p; w=w*w%p; S=add(S,scale(J[M+k],w,p),p)
    lhs=mul(E,S,p)
    lhs=add(lhs,scale(mul(NM,J[2*M+1],p),2*M+1,p),p)
    lhs=add(lhs,mul_1mz2(mul(N0,J[M-1],p),M,p),p)
    assert all(x%p==0 for x in lhs)
    g0=gcd_poly(E,N0,p); gm=gcd_poly(E,NM,p); ga=gcd_poly(g0,NM,p)
    e1=sum(E)%p; pred=1
    for k in range(M+1): pred=pred*(2*(2*M+2*k+1))%p
    assert e1==((-1)**(M+1)*pred)%p
    print('M',M,'p',p,'degE',len(E)-1,'g0',len(g0)-1,'gm',len(gm)-1,'gall',len(ga)-1)
for p in MODS:
    for M in MS: audit(M,p)
print('PASS')
