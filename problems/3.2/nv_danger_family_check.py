#!/usr/bin/env python3
"""Independent audit of thm:nv-range's danger case (written by the auditing
session, NOT by the theorem's author): for primes p <= 3000, construct all
rank-of-apparition danger candidates h = a*pi, d = b*pi, k = (a+b)*pi < p
(k <= 60), build the ACTUAL Delta_{h,k} mod p, center it, and check the
theorem's claims: L = C1 = 0 automatic; actual C2 = 5/256(hB+dD+kBD);
when C2 = 0: B = D = -1, actual C3 = 0, actual C4 = -75hdk/2048 != 0.
Result on 2026-07-19: 309 cases, 164 killed by C2, 145 full C4 checks,
0 problems."""
from math import comb

def polymul(a,b,p):
    r=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        if x:
            for j,y in enumerate(b): r[i+j]=(r[i+j]+x*y)%p
    return r
def polyadd(a,b,p):
    n=max(len(a),len(b)); r=[0]*n
    for i,x in enumerate(a): r[i]=x%p
    for i,y in enumerate(b): r[i]=(r[i]+y)%p
    return r
def polyscale(a,c,p): return [(x*c)%p for x in a]
def shift_poly(a,s,p):
    n=len(a); r=[0]*n
    for i,c in enumerate(a):
        if c:
            for j in range(i+1): r[j]=(r[j]+c*comb(i,j)*pow(s,i-j,p))%p
    return r
def build_N(K,p):
    N=[[0],[1]]
    for m in range(1,K):
        lin=[m%p,1]; l2=polymul(lin,lin,p); l3=polymul(l2,lin,p); l6=polymul(l3,l3,p)
        P=polyadd(polyadd(polyscale(l3,34,p),polyscale(l2,51,p),p),polyadd(polyscale(lin,27,p),[5],p),p)
        N.append(polyadd(polymul(P,N[m],p),polyscale(polymul(l6,N[m-1],p),p-1,p),p))
    return N
def Pi(m,p):
    r=[1]
    for j in range(1,m+1):
        lj=[j%p,1]; r=polymul(r,polymul(polymul(lj,lj,p),lj,p),p)
    return r
def primes_upto(n):
    s=[True]*(n+1); s[0]=s[1]=False
    for i in range(2,int(n**.5)+1):
        if s[i]:
            for j in range(i*i,n+1,i): s[j]=False
    return [i for i in range(2,n+1) if s[i]]

def main():
    cases=0; c2killed=0; c4checked=0; problems=[]
    for p in primes_upto(3000):
        if p<7: continue
        ell=[0,1]; pi=None
        for m in range(1,2*p+2):
            ell.append((34*ell[-1]-ell[-2])%p)
            if ell[m]==0 and m>0: pi=m; break
        if pi is None or pi<2: continue
        for a in range(1,5):
            for b in range(1,5):
                h,d=a*pi,b*pi; k=h+d
                if k>=p or k>60: continue
                while len(ell)<=k: ell.append((34*ell[-1]-ell[-2])%p)
                if ell[h]!=0 or ell[d]!=0: continue
                B=ell[h-1]; D=ell[d-1]
                cases+=1
                N=build_N(k+1,p)
                Delta=polyadd(polyadd(polymul(N[h],shift_poly(Pi(d,p),h,p),p),polyscale(N[k],p-1,p),p),polymul(Pi(h,p),shift_poly(N[d],h,p),p),p)
                n_deg=3*(k-1)
                c=shift_poly(Delta,(-(k+1)*pow(2,-1,p))%p,p)
                while len(c)<n_deg+1: c.append(0)
                L,C1,C2,C3,C4=c[n_deg],c[n_deg-1],c[n_deg-2],c[n_deg-3],c[n_deg-4]
                C2f=5*pow(256,-1,p)*((h*B+d*D+k*B*D)%p)%p
                if L!=0 or C1!=0 or C2!=C2f:
                    problems.append((p,h,k,'C2')); continue
                if C2!=0: c2killed+=1; continue
                if not (B==p-1 and D==p-1):
                    problems.append((p,h,k,'sign')); continue
                C4f=(-75*h*d*k)%p*pow(2048,-1,p)%p
                c4checked+=1
                if C3!=0 or C4!=C4f or C4==0:
                    problems.append((p,h,k,'C4'))
    status="PASS" if not problems else "FAIL"
    print(f"{status}: {cases} danger-family cases; {c2killed} killed by C2; {c4checked} full C4 checks; {len(problems)} problems")
    for pr in problems[:10]: print(pr)
    return 0 if not problems else 1

if __name__=="__main__":
    raise SystemExit(main())
