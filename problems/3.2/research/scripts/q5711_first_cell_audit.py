#!/usr/bin/env python3
from functools import lru_cache
from math import comb, gcd, isqrt
from random import Random

def C(n,k): return comb(n,k) if 0 <= k <= n else 0

@lru_cache(None)
def apery(n):
    if n == 0: return 1
    a,b=1,5
    if n == 1: return b
    for m in range(1,n):
        num=(34*m**3+51*m**2+27*m+5)*b-m**3*a; den=(m+1)**3
        assert num%den==0; a,b=b,num//den
    return b

def apery_mod(r,p):
    if r == 0: return 1
    a,b=1,5%p
    if r == 1: return b
    for m in range(1,r):
        num=((34*m**3+51*m**2+27*m+5)*b-m**3*a)%p
        den=pow(m+1,3,p); a,b=b,num*pow(den,-1,p)%p
    return b

def primes_upto(n):
    s=bytearray(b'\x01')*(n+1)
    if n>=0:s[0]=0
    if n>=1:s[1]=0
    for p in range(2,isqrt(n)+1):
        if s[p]: s[p*p:n+1:p]=b'\x00'*(((n-p*p)//p)+1)
    return [p for p in range(2,n+1) if s[p]]

@lru_cache(None)
def shell(M,d):
    if d>M:return apery(M)
    a=M//d; out=0
    for t in range(M+1):
        X=sum(C(M,t+d*u) for u in range(-a,a+1))
        N=2*M-t; Z=sum(C(N,M+d*v) for v in range(-a,a+1))
        out += C(M,t)*X*Z*Z
    return out

def parts(M,d):
    assert M//2<d<=M
    r=M-d; core=low=0
    for t in range(M+1):
        A=C(M,t); N=2*M-t; B=C(N,M); P=C(N,r)
        core += A*A*(2*B*P+P*P)
        if t<=r:
            Q=C(N,r-t); U=C(M,r-t)
            low += A*(A*(2*(B+P)*Q+Q*Q)+U*(B+P+Q)**2)
    high=0
    for k in range(r+1):
        A=C(M,k); U=C(M,r-k); B=C(M+k,k); P=C(M+k,r)
        high += A*U*(B+P)**2
    b=apery(M); return b,core,low,high,b+core+low+high

def zeta2(M): return sum(C(M,k)**2*C(M+k,k) for k in range(M+1))
def endpoint(M):
    B=C(2*M,M); return apery(M)+2*zeta2(M)+B*B+7*B+11

def dvals(v):
    v=v[:]
    while len(v)>1:v=[v[i+1]-v[i] for i in range(len(v)-1)]
    return v[0]
def delta(M,d,k): return dvals([shell(M,d+i) for i in range(k+1)])
def weight(d,L,i): return (-1)**i*C(d+i,i)*C(d+L+1,L-i)
def G(M,d,L): return sum(weight(d,L,i)*shell(M,d+i) for i in range(L+1))

def targets(n):
    return [(q,*divmod(n,q)) for q in primes_upto(n) if q>isqrt(n) and (lambda ar:1<=ar[1]<=q-2 and apery_mod(ar[1],q)==0)(divmod(n,q))]

def audit_formula(M,d): assert shell(M,d)==parts(M,d)[-1]
def audit_split(M,d,k):
    rows=[parts(M,d+i) for i in range(k+1)]
    assert delta(M,d,k)==sum(dvals([r[j] for r in rows]) for j in (1,2,3))

def analyze(n):
    M=n-1; ts=targets(n); qs=sorted(q for q,a,r in ts if a==1)
    print('ROW',n,'targets',ts)
    if len(qs)<2:return
    d=M//2+1; L=qs[-1]-1-d
    gd=G(M,d,L); g1=G(M,d+1,L); H=delta(M,d,L+1); B=C(d+L+1,L)
    assert gd-g1==(-1)**(L+1)*B*H
    R=gcd(abs(gd),abs(H)); GG=gcd(abs(gd),abs(g1))
    assert GG==R*gcd(abs(gd//R),B)
    w=min(5,M-d-L); vals=[G(M,d+s,L) for s in range(w+1)]
    cg=0
    for x in vals:cg=gcd(cg,abs(x))
    bg=cg
    for s in range(w+1):bg=gcd(bg,abs(delta(M,d+s,L+1)))
    tp=1
    for q in qs:
        if d<=q-1<=d+L: assert gd%q==0; tp*=q
    print(' COVER',M,d,L,'digits',tuple(len(str(abs(x))) for x in (gd,H,R,GG,cg,bg)))
    print(' R',R,'GG',GG,'carrier_gcd',cg,'boundary_gcd',bg,'target_product',tp,'carrier_q',cg//tp if cg%tp==0 else -1)

for M in range(1,35):
    assert shell(M,M)==endpoint(M)
    for d in range(M//2+1,M+1):
        audit_formula(M,d)
        for k in range(1,min(5,M-d)+1):audit_split(M,d,k)
for M in (199,271,299,320,754):
    ds=sorted({M//2+1,M//2+2,3*M//4,M-7,M-1,M})
    for d in ds:
        if M//2<d<=M:audit_formula(M,d)
    assert shell(M,M)==endpoint(M)
    for d in ds:
        if M//2<d and d+4<=M:audit_split(M,d,4)
    print('FORMULA_PASS',M,ds)
rng=Random(5711)
for _ in range(40):
    M=rng.randrange(20,180); d=rng.randrange(M//2+1,M+1); audit_formula(M,d)
    if d<M:audit_split(M,d,rng.randrange(1,min(7,M-d+1)))
print('RANDOM_PASS')
for n in (200,272,300,321,755):analyze(n)
print('ALL_Q5711_CHECKS_PASS')
