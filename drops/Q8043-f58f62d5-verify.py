#!/usr/bin/env python3
"""Q8043 exact finite audit; standard library only."""
from fractions import Fraction
from itertools import combinations
from collections import Counter, defaultdict
from pathlib import Path
from math import comb, gcd
import argparse, json, struct


def A(n): return 34*n**3+51*n**2+27*n+5

def apery(N):
    b=[1] if N==0 else [1,5]
    for n in range(1,N):
        u=A(n)*b[n]-n**3*b[n-1]; d=(n+1)**3
        assert u%d==0; b.append(u//d)
    return b

def companion(N):
    a=[Fraction(0)] if N==0 else [Fraction(0),Fraction(6)]
    for n in range(1,N): a.append((A(n)*a[n]-n**3*a[n-1])/(n+1)**3)
    return a

def franel(N): return [sum(comb(n,k)**3 for k in range(n+1)) for n in range(N+1)]
def primes(N):
    z=bytearray(b'\1')*(N+1); z[:2]=b'\0\0'
    for p in range(2,int(N**.5)+1):
        if z[p]: z[p*p:N+1:p]=b'\0'*(((N-p*p)//p)+1)
    return [i for i in range(2,N+1) if z[i]]
def row(n): return [comb(n,k)*comb(n+k,k) for k in range(n+1)]

def coeffs(n,J,L,F):
    C=[sum(L[i]*F[i] for i in range(J+1))]+[0]*n
    for d in range(1,n+1):
        lo=max(0,J+1-d); hi=min(J,n-d); s=0
        for i in range(lo,hi+1):
            t=L[i+d]*comb(i+d,i)*comb(d-1,J-i)*F[i]
            s += -t if (J-i)&1 else t
        C[d]=s
    return C

def content(C):
    g=0
    for c in C: g=gcd(g,c)
    return abs(g)
def prim(C):
    g=content(C); assert g and all(c%g==0 for c in C)
    return g,[c//g for c in C]
def val(C,x):
    s=0
    for c in reversed(C): s=s*x+c
    return s
def fall(n,r):
    s=1
    for j in range(r): s*=n-j
    return s
def jet(C,x,r): return sum(C[d]*fall(d,r)*x**(d-r) for d in range(r,len(C)))
def vp(n,p):
    n=abs(n); e=0
    while n%p==0: n//=p; e+=1
    return e

def prefix_obstruction():
    b=apery(6); assert b==[1,5,73,1445,33001,819005,21460825]
    P1=val(b,1); P2=val(b,2)
    D11=jet(b,1,1); D12=jet(b,1,2)
    D21=jet(b,2,1); D22=jet(b,2,2)
    dd=P2-P1; curv=P2-2*P1+b[0]
    W1=P1*D21-P2*D11; W2=P1*D22-P2*D12
    assert (P1,P2,curv)==(22312715,1441667719,1397042290)
    assert P1%5==0 and D11%5==0 and D12%5!=0
    assert dd%5!=0 and curv%5==0 and W1%5==0 and W2%5!=0
    return {'P6(1)':P1,"P6'(1)":D11,"P6''(1)":D12,'P6(2)':P2,
            'dd(1,2)':dd,'curvature':curv,'jet_det_order1':W1,
            'jet_det_order2':W2,
            'mod5':[P1%5,D11%5,D12%5,dd%5,curv%5,W1%5,W2%5]}

def gapN(i,h):
    if h==1:return 1
    x,y=1,A(i+1)
    if h==2:return y
    for k in range(2,h): x,y=y,A(i+k)*y-(i+k)**6*x
    return y
def gapD(i,j):
    d=1
    for k in range(i+1,j+1): d*=k**3
    return d

def quotient_audit():
    m=12678; hits=[(379,171,33),(443,274,28),(499,203,25)]
    b=apery(274); a=companion(33); qs=sorted(q for _,_,q in hits)
    for p,r,q in hits:
        assert m==q*p+r and b[r]%p==0 and b[q]%p and gcd(a[q].denominator,p)==1
    pairs=[]
    for i,j in combinations(qs,2):
        N=gapN(i,j-i); det=a[i]*b[j]-a[j]*b[i]
        assert det==Fraction(-6*N,gapD(i,j)); pairs.append((i,j,N,det))
    out=[]
    for p,_,_ in hits:
        R=[N%p for _,_,N,_ in pairs]; assert all(R)
        slopes=[a[q]/b[q] for q in qs]
        S=[s.numerator%p*pow(s.denominator%p,-1,p)%p for s in slopes]
        assert len(set(S))==3
        out.append({'p':p,'gapN_mod_p':R,'slopes_mod_p':S})
    return {'m':m,'hits':hits,'q_pairs':[(i,j,len(str(abs(N)))) for i,j,N,_ in pairs],
            'prime_units':out}

def shell_jets(Xmax):
    N=Xmax*Xmax-1; F=franel(N); b=apery(N); ps=primes(2*Xmax)
    st=Counter(); first=[]; fail_support=fail_jet=fail_W=None
    targets_by_X=defaultdict(dict)
    for X in range(3,Xmax+1):
        shell=[p for p in ps if X<p<=2*X]
        for n in range(X,X*X):
            T=[]
            for p in shell:
                r=n%p; j=min(r,p-1-r)
                if b[j]%p==0:T.append((p,r,j))
            if not T:continue
            targets_by_X[X][n]={p for p,_,_ in T}; st['rows']+=1; st['hits']+=len(T)
            L=row(n); P={}
            for J in (X-1,X):
                G,Q=prim(coeffs(n,J,L,F)); P[J]=(G,Q)
                for p,r,j in T:
                    bad=[d for d,c in enumerate(Q) if d%p and c%p]
                    jets=[jet(Q,1,k)%p for k in (1,2,3)]
                    rec={'X':X,'n':n,'J':J,'p':p,'r':r,'j':j,
                         'vGamma':vp(G,p),'vbj':vp(b[j],p),'P1':val(Q,1)%p,
                         'jets':jets,'bad_degree':bad[0] if bad else None}
                    st['tests']+=1; st['support_pass']+=not bad; st['jet_pass']+=not any(jets)
                    if len(first)<10:first.append(rec)
                    if bad and fail_support is None:fail_support=rec
                    if any(jets) and fail_jet is None:fail_jet=rec
            G0,Q0=P[X-1]; G1,Q1=P[X]
            v0,v1=val(Q0,1),val(Q1,1); d0,d1=jet(Q0,1,1),jet(Q1,1,1)
            W=v0*d1-v1*d0
            for p,r,j in T:
                st['Wtests']+=1; st['Wpass']+=W%p==0
                if W%p and fail_W is None:
                    fail_W={'X':X,'n':n,'p':p,'r':r,'j':j,'Wmodp':W%p,
                            'Pmods':[v0%p,v1%p],'Dmods':[d0%p,d1%p],
                            'Wbits':abs(W).bit_length()}
    common={'count':0}; union={'count':0}; c7=u7=0
    for X,R in targets_by_X.items():
        for (m0,s0),(m1,s1) in combinations(sorted(R.items()),2):
            C=sorted(s0&s1); U=sorted(s0|s1)
            if len(C)>common['count']:common={'count':len(C),'X':X,'m':[m0,m1],'p':C,'rows':[sorted(s0),sorted(s1)]}
            if len(U)>union['count']:union={'count':len(U),'X':X,'m':[m0,m1],'p':U,'rows':[sorted(s0),sorted(s1)]}
            c7+=len(C)>=7; u7+=len(U)>=7
    return {'range':[3,Xmax],'stats':dict(st),'first':first,
            'first_support_failure':fail_support,'first_jet_failure':fail_jet,
            'first_W_failure':fail_W,'max_common':common,'max_union':union,
            'pairs_common_ge7':c7,'pairs_union_ge7':u7}

def census(root):
    data=(root/'problems/3.2/data_zp_pairs.bin').read_bytes(); assert len(data)%8==0
    Z=defaultdict(set); H=defaultdict(list)
    for o in range(0,len(data),8):
        p,r=struct.unpack_from('<II',data,o); Z[p].add(r); H[p+r].append(p)
    mx=max(map(len,H.values())); M=sorted(n for n,v in H.items() if len(v)==mx); assert mx==3
    m=12678; hits=sorted((p,m%p,m//p) for p,z in Z.items() if 256<p<=512 and m%p in z)
    assert hits==[(379,171,33),(443,274,28),(499,203,25)]
    first=[{'n':n,'p':sorted(H[n])} for n in M[:2]]
    U=sorted(set(first[0]['p'])|set(first[1]['p']))
    return {'records':len(data)//8,'max_top_half':mx,'maximizers':len(M),
            'first_two':first,'union':U,'union_count':len(U),
            'two_rows_can_have_7':2*mx>=7,'m12678_hits':hits}
def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--shell-x-max',type=int,default=24); z=ap.parse_args()
    root=Path(__file__).resolve().parents[1]; assert (root/'problems/3.2').is_dir()
    R={'prefix':prefix_obstruction(),'quotient':quotient_audit(),
       'shell':shell_jets(z.shell_x_max),'census':census(root)}
    print('Q8043_EXACT_AUDIT'); print(json.dumps(R,sort_keys=True,indent=2)); print('Q8043_VERIFIER_PASS')
if __name__=='__main__':main()
