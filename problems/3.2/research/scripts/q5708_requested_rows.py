#!/usr/bin/env python3
"""Focused exact Q5708 audit for n=200,272,300,321,755.

This is deliberately small: it finds large-prime targets, takes adjacent
pairs in each fixed quotient cell, and computes the exact Newton residuals.
"""
from functools import lru_cache
from math import comb, gcd, isqrt
from sympy import factorint, isprime


def C(n,k): return comb(n,k) if 0 <= k <= n else 0


def primes_upto(n):
    s=bytearray(b'\x01')*(n+1); s[:2]=b'\x00\x00'
    for p in range(2,isqrt(n)+1):
        if s[p]: s[p*p:n+1:p]=b'\x00'*(((n-p*p)//p)+1)
    return [p for p in range(2,n+1) if s[p]]


def apery_mod(r,p):
    return sum((C(r,k)*C(r+k,k))**2 for k in range(r+1))%p


@lru_cache(None)
def shell(M,d):
    a=M//d
    out=0
    for t in range(M+1):
        X=sum(C(M,M-t+d*u) for u in range(-a,a+1))
        Z=sum(C(2*M-t,M-t+d*v) for v in range(-a,a+1))
        out += C(M,t)*X*Z*Z
    return out


def delta(M,d,k):
    return sum((-1)**(k-i)*C(k,i)*shell(M,d+i) for i in range(k+1))


def weight(d,L,i):
    return (-1)**i*C(d+i,i)*C(d+L+1,L-i)


def G(M,d,L):
    return sum(weight(d,L,i)*shell(M,d+i) for i in range(L+1))


def val(x,p):
    x=abs(x); e=0
    while x and x%p==0: x//=p; e+=1
    return e


def fac(n):
    n=abs(n)
    if n in (0,1): return str(n)
    ff=factorint(n, limit=10**6, use_ecm=True)
    parts=[]
    for p,e in sorted((int(p),int(e)) for p,e in ff.items()):
        tag='' if isprime(p) else '[C]'
        parts.append(f'{p}{tag}' if e==1 else f'{p}{tag}^{e}')
    return ' * '.join(parts)


def targets(n):
    out=[]
    for q in primes_upto(n):
        if q<=isqrt(n): continue
        a,r=divmod(n,q)
        if 1<=r<=q-2 and apery_mod(r,q)==0: out.append((q,a,r))
    return out


def pairs(n):
    ts=targets(n); out=[]
    for a in sorted(set(a for _,a,_ in ts)):
        qs=sorted(q for q,aa,_ in ts if aa==a)
        out += [(q,e,a) for q,e in zip(qs,qs[1:])]
    return ts,out


def analyze(n,q,e,a):
    M=n-a; d=q-1; L=e-q
    gd=G(M,d,L); g1=G(M,d+1,L); H=delta(M,d,L+1)
    B=C(d+L+1,L)
    assert gd-g1 == (-1)**(L+1)*B*H
    R=gcd(gd,H); GG=gcd(gd,g1)
    x=gd//R
    assert GG == R*gcd(x,B)
    gm=G(M,d,L-1); g1m=G(M,d+1,L-1)
    assert L*gd+q*g1m == e*gm
    G4=gcd(gcd(abs(gd),abs(g1)),gcd(abs(gm),abs(g1m)))
    U=delta(M,d,L); V=delta(M,d+1,L)
    rho=gcd(abs(gm),gcd(abs(U),abs(V)))
    A=C(d+L,L); Cc=C(d+L,L-1); gg=gcd(A,Cc)
    smith=(1,gg,B)
    print(f'PAIR n={n} q={q} ell={e} a={a} M={M} d={d} L={L}')
    print(' digits G,H,R,GG,G4,rho=',tuple(len(str(abs(z))) for z in (gd,H,R,GG,G4,rho)))
    print(' R=',fac(R))
    print(' GG=',fac(GG))
    print(' G4=',fac(G4))
    print(' rho=',fac(rho))
    print(' v_R(q,ell)=',(val(R,q),val(R,e)),'v_GG=',(val(GG,q),val(GG,e)),'v_G4=',(val(G4,q),val(G4,e)),'v_rho=',(val(rho,q),val(rho,e)))
    print(' B=',fac(B),'v_B=',(val(B,q),val(B,e)))
    print(' A,C,g=',fac(A),fac(Cc),fac(gg),'smith=',smith)
    print(' primitive Pascal coefficients=',(L,q,e))


def main():
    for n in (200,272,300,321,755):
        ts,ps=pairs(n)
        print('ROW',n,'targets',ts,'pairs',ps)
        for q,e,a in ps: analyze(n,q,e,a)

    # Wider cheap scan: target-pair geometry only, exact through n=2000.
    hist={}; rows=[]
    for n in range(20,2001):
        _,ps=pairs(n)
        for q,e,a in ps:
            L=e-q
            hist[L]=hist.get(L,0)+1
            if L<=40: rows.append((n,q,e,a,L))
    print('WIDE n<=2000 same-cell pair count=',sum(hist.values()))
    print('WIDE gap histogram=',sorted(hist.items()))
    print('WIDE pairs L<=40=',rows)

if __name__=='__main__': main()
