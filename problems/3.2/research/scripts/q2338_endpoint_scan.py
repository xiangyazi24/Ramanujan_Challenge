#!/usr/bin/env python3
"""Exact audits for Q2338 using only Python's standard library.

The script checks three distinct issues.

1. Exact characteristic-zero gcd data for
       gcd(b_r,b_s,K(r,s)).
2. Every strict lower-half Apéry double-zero pair for primes p <= PMAX,
   including whether gcd(P_r,C_{r,s}) has positive degree.
3. The second-kind/Dirichlet endpoint formulas
       K = -(r^3 b_{r-1}/6) N_h^(r)(lambda_r) / (r+1...s)^3
         =  (s^3 b_{s-1}/6) N_h^(r)(lambda_s) / (r+1...s)^3.

Here C_{r,s}=S_{s-r-1}^{(r)} is the degree-(s-r-1) transfer block.
"""
from __future__ import annotations

import json
import math
import random
from itertools import combinations

PMAX = 5000
EXACT_MAX = 55


def primes_upto(n: int) -> list[int]:
    sieve = bytearray(b"\x01") * (n + 1)
    sieve[:2] = b"\x00\x00"
    for q in range(2, math.isqrt(n) + 1):
        if sieve[q]:
            sieve[q*q:n+1:q] = b"\x00" * (((n-q*q)//q)+1)
    return [q for q in range(2, n + 1) if sieve[q]]


# ---------- exact integers ----------

def D_exact(n: int) -> list[int]:
    row = [1]
    cur = 1
    for k in range(n):
        cur = cur * (n-k) * (n+k+1) // ((k+1)*(k+1))
        row.append(cur)
    return row


def exact_rows(nmax: int) -> list[list[int]]:
    return [D_exact(n) for n in range(nmax + 1)]


def K_exact(rows: list[list[int]], r: int, s: int) -> int:
    return sum(x*y for x, y in zip(rows[r], rows[s]))


# Pollard-rho is used only to report prime support of usually small gcds.
def is_probable_prime(n: int) -> bool:
    if n < 2:
        return False
    small = (2,3,5,7,11,13,17,19,23,29,31,37)
    for q in small:
        if n % q == 0:
            return n == q
    d, t = n-1, 0
    while d % 2 == 0:
        d //= 2
        t += 1
    bases = small if n.bit_length() > 64 else (2,325,9375,28178,450775,9780504,1795265022)
    for a in bases:
        a %= n
        if a in (0,1):
            continue
        x = pow(a, d, n)
        if x in (1, n-1):
            continue
        for _ in range(t-1):
            x = x*x % n
            if x == n-1:
                break
        else:
            return False
    return True


def pollard(n: int) -> int:
    if n % 2 == 0:
        return 2
    if n % 3 == 0:
        return 3
    while True:
        c = random.randrange(1, n-1)
        x = random.randrange(2, n-1)
        y = x
        d = 1
        while d == 1:
            x = (x*x+c) % n
            y = (y*y+c) % n
            y = (y*y+c) % n
            d = math.gcd(abs(x-y), n)
        if d != n:
            return d


def factor(n: int, out: list[int]) -> None:
    if n == 1:
        return
    if is_probable_prime(n):
        out.append(n)
        return
    d = pollard(n)
    factor(d, out)
    factor(n//d, out)


def exact_gcd_audit() -> dict[str, object]:
    rows = exact_rows(EXACT_MAX)
    b = [sum(x*x for x in row) for row in rows]
    anomalies = []
    support_records = []
    maximum_ratio = (0.0, None)
    for r in range(1, EXACT_MAX):
        for s in range(r+1, EXACT_MAX+1):
            kval = K_exact(rows, r, s)
            g = math.gcd(math.gcd(b[r], b[s]), kval)
            if g == 1:
                continue
            fs: list[int] = []
            factor(g, fs)
            primes = sorted(set(fs))
            largest = max(primes)
            ratio = largest/s
            if ratio > maximum_ratio[0]:
                maximum_ratio = (ratio, (r,s,g,primes))
            large = [q for q in primes if q > 2*s+1]
            if large:
                anomalies.append({"r":r,"s":s,"gcd":g,"prime_factors":primes,"large":large})
            if len(support_records) < 40 or large:
                support_records.append({"r":r,"s":s,"gcd":g,"prime_factors":primes})
    return {
        "exact_max": EXACT_MAX,
        "nontrivial_samples": support_records,
        "prime_factor_gt_2s_plus_1": anomalies,
        "max_largest_prime_over_s": maximum_ratio,
    }


# ---------- finite-field polynomials ----------

def trim(a: list[int], p: int) -> list[int]:
    a = [x % p for x in a] or [0]
    while len(a)>1 and a[-1]==0:
        a.pop()
    return a


def scale(a: list[int], c: int, p: int) -> list[int]:
    return trim([(c*x)%p for x in a], p)


def sub(a: list[int], b: list[int], p: int) -> list[int]:
    n=max(len(a),len(b)); z=[0]*n
    for i in range(n):
        z[i]=((a[i] if i<len(a) else 0)-(b[i] if i<len(b) else 0))%p
    return trim(z,p)


def mul_linear(a: list[int], c: int, p: int) -> list[int]:
    # (T+c)*a
    z=[0]*(len(a)+1)
    for i,x in enumerate(a):
        z[i]=(z[i]+c*x)%p
        z[i+1]=(z[i+1]+x)%p
    return trim(z,p)


def divmod_poly(a: list[int], b: list[int], p: int) -> tuple[list[int],list[int]]:
    a=trim(a,p); b=trim(b,p)
    q=[0]*max(1,len(a)-len(b)+1)
    ib=pow(b[-1],-1,p)
    while a != [0] and len(a)>=len(b):
        d=len(a)-len(b); c=a[-1]*ib%p; q[d]=c
        for j,y in enumerate(b):
            a[d+j]=(a[d+j]-c*y)%p
        a=trim(a,p)
    return trim(q,p),a


def gcd_poly(a: list[int], b: list[int], p: int) -> list[int]:
    while trim(b,p) != [0]:
        a,b=b,divmod_poly(a,b,p)[1]
    a=trim(a,p)
    return scale(a,pow(a[-1],-1,p),p)


def eval_poly(a: list[int], x: int, p: int) -> int:
    z=0
    for c in reversed(a):
        z=(z*x+c)%p
    return z


def aa(n: int,p: int)->int:
    return (n*n+n+1)*pow(2,-1,p)%p


def beta(n: int,p: int)->int:
    return pow(n,6,p)*pow(4*(4*n*n-1)%p,-1,p)%p


def build_P(nmax:int,p:int)->list[list[int]]:
    ps=[[1]]
    if nmax==0:return ps
    ps.append([aa(0,p),1])
    for n in range(1,nmax):
        ps.append(sub(mul_linear(ps[n],aa(n,p),p),scale(ps[n-1],beta(n,p),p),p))
    return ps


def transfer_C(r:int,s:int,p:int)->list[int]:
    j=s-r-1
    if j==0:return [1]
    prev=[1]; cur=[aa(r+1,p),1]
    for k in range(1,j):
        n=r+k+1
        prev,cur=cur,sub(mul_linear(cur,aa(n,p),p),scale(prev,beta(n,p),p),p)
    return cur


def apery_half(p:int)->list[int]:
    m=(p-1)//2
    vals=[1]
    if m==0:return vals
    vals.append(5%p)
    for n in range(1,m):
        middle=(34*n**3+51*n*n+27*n+5)%p
        vals.append((middle*vals[n]-pow(n,3,p)*vals[n-1])%p*pow(pow(n+1,3,p),-1,p)%p)
    return vals[:m+1]


def row_D_mod(n:int,p:int)->list[int]:
    row=[1]; cur=1
    for k in range(n):
        cur=cur*(n-k)*(n+k+1)%p*pow((k+1)*(k+1)%p,-1,p)%p
        row.append(cur)
    return row


def K_mod(r:int,s:int,p:int)->int:
    a=row_D_mod(r,p); b=row_D_mod(s,p)
    return sum(x*y for x,y in zip(a,b))%p


def lam(n:int,p:int)->int:
    return n*(n+1)%p


def gamma(n:int,T:int,p:int)->int:
    return (pow(n+1,3,p)+pow(n,3,p)+2*(2*n+1)*T)%p


def endpoint_N(r:int,s:int,T:int,p:int)->int:
    # N_0=0,N_1=1; return N_{s-r}^{(r)}(T)
    h=s-r
    prev,cur=0,1
    if h==1:return cur
    for j in range(1,h):
        n=r+j
        prev,cur=cur,(gamma(n,T,p)*cur-pow(n,6,p)*prev)%p
    return cur


def double_zero_audit()->dict[str,object]:
    marked=[]; positive_gcd=[]; degree_gt_one=[]; bad_K=[]
    primes=primes_upto(PMAX)
    for p in primes:
        if p<7:continue
        vals=apery_half(p); m=(p-1)//2
        zeros=[n for n in range(1,m) if vals[n]==0]
        if len(zeros)<2:continue
        ps=None
        for r,s in combinations(zeros,2):
            kval=K_mod(r,s,p)
            if ps is None:
                ps=build_P(max(zeros),p)
            c=transfer_C(r,s,p)
            g=gcd_poly(ps[r],c,p); gd=len(g)-1
            if kval==0:bad_K.append((p,r,s,gd,g))
            if gd>0:positive_gcd.append((p,r,s,kval,gd,g))
            if gd>1:degree_gt_one.append((p,r,s,kval,gd,g))

            nr=endpoint_N(r,s,lam(r,p),p)
            ns=endpoint_N(r,s,lam(s,p),p)
            den=1
            for j in range(r+1,s+1):den=den*pow(j,3,p)%p
            left=(-pow(r,3,p)*vals[r-1]*pow(6,-1,p)*nr*pow(den,-1,p))%p
            right=(pow(s,3,p)*vals[s-1]*pow(6,-1,p)*ns*pow(den,-1,p))%p
            assert left==right==kval
            product_identity=(
                kval*kval +
                pow(r,3,p)*pow(s,3,p)*vals[r-1]*vals[s-1]*pow(36,-1,p)
                *nr*ns*pow(den,-2,p)
            )%p
            assert product_identity==0
            marked.append({
                "p":p,"r":r,"s":s,"h":s-r,"K":kval,
                "gcd_degree":gd,"endpoint_N_lambda_r":nr,
                "endpoint_N_lambda_s":ns,
            })
    return {
        "pmax":PMAX,
        "marked_pair_count":len(marked),
        "marked_pairs":marked,
        "positive_associated_gcd":positive_gcd,
        "associated_gcd_degree_gt_one":degree_gt_one,
        "K_zero_pairs":bad_K,
    }


def main()->None:
    random.seed(2338)
    out={
        "exact_gcd_audit":exact_gcd_audit(),
        "double_zero_audit":double_zero_audit(),
    }
    print("ANSWER Q2338 c289a614")
    print(json.dumps(out,sort_keys=True,separators=(",",":")))


if __name__=="__main__":
    main()
