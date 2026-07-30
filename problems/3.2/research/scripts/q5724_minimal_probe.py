#!/usr/bin/env python3
"""Platform-independent modular probe for Q5724.

Uses the exact first-cell shell formula, vectorized finite-field RREF,
monotone binary search in Ore order/degree, held-out verification, SymPy
endpoint factoring, and exact augmented-state gcd scans. No floating-point
arithmetic enters the recurrence calculation.
"""
from __future__ import annotations
from math import comb, gcd
import numpy as np

P = 1000003
MAX_ORDER = 38
MAX_DEGREE = 10
MS = (1400, 1600, 1800)


def binom_table_mod(n: int, p: int):
    fact=[1]*(n+1)
    for i in range(1,n+1): fact[i]=fact[i-1]*i%p
    invfact=[1]*(n+1)
    invfact[n]=pow(fact[n],p-2,p)
    for i in range(n,0,-1): invfact[i-1]=invfact[i]*i%p
    def C(a,b):
        if b<0 or b>a:return 0
        return fact[a]*invfact[b]%p*invfact[a-b]%p
    return C


def F_sequence_mod(M: int, p: int) -> list[int]:
    C=binom_table_mod(2*M,p)
    b=0
    for k in range(M+1):
        z=C(M,k)*C(M+k,k)%p
        b=(b+z*z)%p
    out=[]
    for r in range((M-1)//2+1):
        d=M-r; total=0
        for t in range(M+1):
            A=C(M,t)
            X=(C(M,t-d)+A+C(M,t+d))%p
            N=2*M-t
            Z=(C(N,M-d)+C(N,M)+C(N,M+d))%p
            total=(total+A*X%p*Z%p*Z)%p
        out.append((total-b)%p)
    return out


def make_matrix(seq: list[int], order: int, degree: int, row_ids) -> np.ndarray:
    rows=[]
    for r in row_ids:
        pw=[1]
        for _ in range(degree):pw.append(pw[-1]*r%P)
        row=[]
        for i in range(order+1):
            y=seq[r+i]
            row.extend(y*z%P for z in pw)
        rows.append(row)
    return np.asarray(rows,dtype=np.int64)


def rref_nullity(A: np.ndarray, p: int, need_vector=False):
    if A.size==0:return 0,None
    a=A.copy()%p
    m,n=a.shape;piv=[];rr=0
    for c in range(n):
        nz=np.flatnonzero(a[rr:,c])
        if nz.size==0:continue
        pivot=rr+int(nz[0])
        if pivot!=rr:a[[rr,pivot]]=a[[pivot,rr]]
        inv=pow(int(a[rr,c]),p-2,p)
        ar=(a[rr,c:]*inv)%p
        a[rr,c:]=ar
        factors=a[:,c].copy();factors[rr]=0
        if np.any(factors):
            a[:,c:] = (a[:,c:] - factors[:,None]*ar[None,:]) % p
        piv.append(c);rr+=1
        if rr==m:break
    nullity=n-rr
    if not need_vector or nullity==0:return nullity,None
    pset=set(piv);free=next(c for c in range(n) if c not in pset)
    v=np.zeros(n,dtype=np.int64);v[free]=1
    # Matrix is in reduced form, so one free variable suffices.
    for i,c in enumerate(piv):v[c]=(-a[i,free])%p
    return nullity,[int(x) for x in v]


def has_recurrence(seq,order,degree):
    cols=(order+1)*(degree+1); neq=len(seq)-order
    if neq<cols+16:return False
    train=min(neq-16,cols+20)
    return rref_nullity(make_matrix(seq,order,degree,range(train)),P)[0]>0


def first_true(lo,hi,pred):
    if not pred(hi):return None
    while lo<hi:
        mid=(lo+hi)//2
        if pred(mid):hi=mid
        else:lo=mid+1
    return lo


def verify(seq,order,degree,v,start=0):
    for r in range(start,len(seq)-order):
        pw=[1]
        for _ in range(degree):pw.append(pw[-1]*r%P)
        z=0;c=0
        for i in range(order+1):
            for k in range(degree+1):
                z=(z+v[c]*pw[k]%P*seq[r+i])%P;c+=1
        if z:return False,r,z
    return True,None,None


def solve_minimal(seq):
    order=first_true(1,MAX_ORDER,lambda o:has_recurrence(seq,o,MAX_DEGREE))
    if order is None:return None
    degree=first_true(0,MAX_DEGREE,lambda d:has_recurrence(seq,order,d))
    cols=(order+1)*(degree+1); neq=len(seq)-order
    train=min(neq-24,cols+30)
    nul,v=rref_nullity(make_matrix(seq,order,degree,range(train)),P,True)
    if nul!=1:return ('NONCYCLIC',order,degree,nul)
    assert verify(seq,order,degree,v,train)[0]
    nul2,v2=rref_nullity(make_matrix(seq,order,degree,range(neq)),P,True)
    assert nul2==1 and verify(seq,order,degree,v2)[0]
    q=next(x for x in reversed(v2) if x);iq=pow(q,P-2,P)
    return order,degree,[x*iq%P for x in v2]


def poly_coeff(v,i,d):return [v[i*(d+1)+k] for k in range(d+1)]


def factor_endpoint(coeff):
    import sympy as sp
    x=sp.symbols('r')
    f=sum(int(c)*x**i for i,c in enumerate(coeff))
    return sp.factor_list(sp.Poly(f,x,modulus=P),modulus=P)


def probe(M):
    seq=F_sequence_mod(M,P);ans=solve_minimal(seq)
    print('PROBE',M,P,'terms',len(seq),'ans',None if ans is None else ans[:2],flush=True)
    if ans is None or ans[0]=='NONCYCLIC':return ans
    o,d,v=ans;c0=poly_coeff(v,0,d);co=poly_coeff(v,o,d)
    print('TRAIL',factor_endpoint(c0),flush=True)
    print('LEAD',factor_endpoint(co),flush=True)
    print('HELD',verify(seq,o,d,v,max(0,len(seq)-o-60)),flush=True)
    print('VECTOR',v,flush=True)
    return ans


def apery_exact(M):return sum((comb(M,k)*comb(M+k,k))**2 for k in range(M+1))


def F_exact(M):
    b=apery_exact(M);out=[]
    def C(n,k):return comb(n,k) if 0<=k<=n else 0
    for r in range((M-1)//2+1):
        d=M-r;total=0
        for t in range(M+1):
            A=C(M,t);X=C(M,t-d)+A+C(M,t+d);N=2*M-t
            Z=C(N,M-d)+C(N,M)+C(N,M+d);total+=A*X*Z*Z
        out.append(total-b)
    return b,out


def scan_aug(width):
    import sympy as sp
    for M in (126,146,147,148,149,150):
        b,s=F_exact(M);hits=[]
        for r in range(len(s)-width+1):
            g=abs(b)
            for z in s[r:r+width]:g=gcd(g,abs(z))
            if g>1:hits.append((r,g,sp.factorint(g)))
        print('AUG',M,'width',width,'hits',hits,flush=True)


def main():
    answers=[probe(M) for M in MS]
    od={(a[0],a[1]) for a in answers if a and a[0]!='NONCYCLIC'}
    print('SUMMARY',od,flush=True)
    if len(od)==1:
        o,d=next(iter(od));scan_aug(o);scan_aug(o+1)
    print('PASS',flush=True)

if __name__=='__main__':main()
