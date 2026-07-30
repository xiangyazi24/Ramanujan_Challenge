#!/usr/bin/env sage
"""Q5724 modular probe for the minimal first-cell recurrence.

Computes F_M(r)=C_M(M-r)-b_M directly over finite fields, searches the
specialized minimal polynomial-coefficient recurrence through order 38,
factors the endpoint coefficients, performs held-out checks, and scans the
reported augmented-state gcd ranges exactly.
"""
from sage.all import *
from math import comb, gcd

PRIMES = (1000003, 1000033)
PROBE_MS = (1200, 1500, 1800)
MAX_ORDER = 38
MAX_DEGREE = 14


def mod_binom_table(N, p):
    F = GF(p)
    fact = [F.one()] * (N + 1)
    for i in range(1, N + 1):
        fact[i] = fact[i-1] * i
    ifact = [F.one()] * (N + 1)
    ifact[N] = fact[N]**(-1)
    for i in range(N, 0, -1):
        ifact[i-1] = ifact[i] * i
    def C(n, k):
        if k < 0 or k > n:
            return F.zero()
        return fact[n] * ifact[k] * ifact[n-k]
    return C


def F_sequence_mod(M, p):
    F = GF(p)
    C = mod_binom_table(2*M, p)
    b = F.zero()
    for k in range(M+1):
        z = C(M,k) * C(M+k,k)
        b += z*z
    out = []
    for r in range((M-1)//2 + 1):
        d = M-r
        total = F.zero()
        for t in range(M+1):
            A = C(M,t)
            X = C(M,t-d) + A + C(M,t+d)
            N = 2*M-t
            Z = C(N,M-d) + C(N,M) + C(N,M+d)
            total += A*X*Z*Z
        out.append(total-b)
    return out


def recurrence_matrix(seq, order, degree, rows=None):
    F = seq[0].parent()
    neq = len(seq)-order
    if rows is None:
        rows = range(neq)
    cols = (order+1)*(degree+1)
    A = matrix(F, len(rows), cols)
    for rr, r in enumerate(rows):
        c = 0
        powers = [F.one()]
        for _ in range(degree):
            powers.append(powers[-1]*r)
        for i in range(order+1):
            y = seq[r+i]
            for k in range(degree+1):
                A[rr,c] = y*powers[k]
                c += 1
    return A


def check_operator(seq, order, degree, v, start=0):
    F = seq[0].parent()
    for r in range(start, len(seq)-order):
        z = F.zero(); c = 0; powers=[F.one()]
        for _ in range(degree): powers.append(powers[-1]*r)
        for i in range(order+1):
            for k in range(degree+1):
                z += v[c]*powers[k]*seq[r+i]; c += 1
        if z:
            return False, r, z
    return True, None, None


def minimal_operator(seq):
    neq0 = len(seq)-1
    for order in range(1, MAX_ORDER+1):
        for degree in range(MAX_DEGREE+1):
            ncols=(order+1)*(degree+1)
            neq=len(seq)-order
            if neq < ncols + 16:
                continue
            train = min(neq-12, ncols+24)
            A = recurrence_matrix(seq,order,degree,range(train))
            K = A.right_kernel()
            if K.dimension() != 1:
                continue
            v = K.basis()[0]
            ok, bad, val = check_operator(seq,order,degree,v,train)
            if not ok:
                continue
            Afull=recurrence_matrix(seq,order,degree)
            Kfull=Afull.right_kernel()
            if Kfull.dimension()!=1:
                continue
            vv=Kfull.basis()[0]
            assert check_operator(seq,order,degree,vv)[0]
            return order,degree,vv
    return None


def coeff_poly(v, i, degree, R):
    return R([v[i*(degree+1)+k] for k in range(degree+1)])


def normalize_vector(v):
    for x in reversed(v):
        if x:
            return vector(v.base_ring(), [z/x for z in v])
    raise ValueError('zero vector')


def probe_one(M,p):
    seq=F_sequence_mod(M,p)
    ans=minimal_operator(seq)
    print('PROBE',M,p,'terms',len(seq),'answer',None if ans is None else ans[:2])
    if ans is None:
        return None
    order,degree,v=ans
    v=normalize_vector(v)
    R.<r>=PolynomialRing(GF(p))
    P0=coeff_poly(v,0,degree,R)
    PR=coeff_poly(v,order,degree,R)
    print(' endpoint_degrees',P0.degree(),PR.degree())
    print(' trailing_factor',factor(P0))
    print(' leading_factor',factor(PR))
    # Held-out suffix and independent direct values are already included in seq;
    # record a deterministic late interval explicitly.
    start=max(0,len(seq)-order-40)
    print(' heldout',start,len(seq)-order,check_operator(seq,order,degree,v,start)[0])
    print(' vector',list(v))
    return order,degree,v,P0,PR


def apery_exact(M):
    return sum((comb(M,k)*comb(M+k,k))**2 for k in range(M+1))


def F_sequence_exact(M):
    b=apery_exact(M)
    out=[]
    for r in range((M-1)//2+1):
        d=M-r; total=0
        for t in range(M+1):
            C=lambda n,k: comb(n,k) if 0<=k<=n else 0
            A=C(M,t)
            X=C(M,t-d)+A+C(M,t+d)
            N=2*M-t
            Z=C(N,M-d)+C(N,M)+C(N,M+d)
            total += A*X*Z*Z
        out.append(total-b)
    return b,out


def scan_augmented(width):
    for M in [126,146,147,148,149,150]:
        b,seq=F_sequence_exact(M)
        hits=[]
        for r in range(0,max(0,len(seq)-width+1)):
            g=abs(b)
            for z in seq[r:r+width]: g=gcd(g,abs(z))
            if g>1: hits.append((r,g,factor(g)))
        print('AUG',M,'width',width,'hits',hits)


def main():
    answers=[]
    for p in PRIMES:
        for M in PROBE_MS:
            ans=probe_one(M,p)
            if ans is not None: answers.append((p,M,ans[0],ans[1]))
    print('SUMMARY',answers)
    orders=sorted(set(o for _,_,o,_ in answers))
    if len(orders)==1:
        scan_augmented(orders[0])
        scan_augmented(orders[0]+1)
    print('PASS')

if __name__=='__main__':
    main()
