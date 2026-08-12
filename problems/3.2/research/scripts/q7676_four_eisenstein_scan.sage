from sage.all import *
from math import gcd as pygcd
from collections import defaultdict

PMAX = 20000

def P(n):
    n=ZZ(n); return 34*n^3+51*n^2+27*n+5

def canon(v):
    v=tuple(int(x) for x in v)
    if not any(v): return None
    g=0
    for x in v: g=pygcd(g,abs(x))
    v=tuple(x//g for x in v)
    for x in v:
        if x:
            if x<0: v=tuple(-y for y in v)
            break
    return v

def dotmod(c,v,p): return sum(int(c[i])*int(v[i]) for i in range(4))%int(p)

def companion_vectors_for_prime(p,b,targets,avec):
    p=ZZ(p); m=max(targets); Fp=GF(p)
    qcoef=[Fp(1)]
    if m>=1: qcoef.append(Fp(17))
    for n in range(1,m):
        qcoef.append((Fp(17*(2*n+1))*qcoef[n]-Fp(n)*qcoef[n-1])/Fp(n+1))
    R=PowerSeriesRing(Fp,'t',default_prec=m+1)
    Fser=R([Fp(x) for x in b[:m+1]]).add_bigoh(m+1)
    Qser=R(qcoef).add_bigoh(m+1)
    Cser=(Qser/(Fser*Fser)).add_bigoh(m+1)
    g=[Fp(Cser[n]) for n in range(m+1)]
    srcO1=[Fp(0)]*(m+1); srcO2=[Fp(0)]*(m+1); srcP1=[Fp(0)]*(m+1); srcP2=[Fp(0)]*(m+1)
    for n in range(1,m+1):
        qprev=qcoef[n-1]
        srcO1[n]=Fp(35)*g[n]-(Fp(5) if n==1 else Fp(0))
        srcO2[n]=Fp(5)*g[n]-(Fp(35) if n==1 else Fp(0))
        srcP1[n]=Fp(37)*qcoef[n]-Fp(130)*qprev-Fp(37)*g[n]
        srcP2[n]=Fp(13)*qcoef[n]-Fp(130)*qprev-Fp(13)*g[n]
    def solve(src):
        z=[Fp(0)]*(m+1)
        for n in range(1,m+1):
            rhs=Fp(P(n-1))*z[n-1]+src[n]
            if n>=2: rhs-=Fp((n-1)^3)*z[n-2]
            z[n]=rhs/Fp(n^3)
        return z
    o1,o2,p1,p2=map(solve,[srcO1,srcO2,srcP1,srcP2])
    out={}
    for r in targets:
        v=tuple(int(x) for x in ((o1[r]+p1[r])/Fp(2),(o2[r]+p2[r])/Fp(8),(p2[r]-o2[r])/Fp(18),(p1[r]-o1[r])/Fp(72)))
        sep=dotmod((1,-28,63,-36),v,p)
        assert sep==(40*int(avec[r]))%int(p)
        if p>=7: assert sep!=0
        out[r]=v
    return out

records=[]; prime_target_counts=[]
for p in prime_range(5,PMAX+1):
    Fp=GF(p); b=[Fp(1),Fp(5)]; a=[Fp(0),Fp(6)]
    for n in range(2,int(p)):
        b.append((Fp(P(n-1))*b[n-1]-Fp((n-1)^3)*b[n-2])/Fp(n^3))
        a.append((Fp(P(n-1))*a[n-1]-Fp((n-1)^3)*a[n-2])/Fp(n^3))
    targets=[r for r in range(1,int(p)) if b[r]==0]
    if not targets: continue
    vecs=companion_vectors_for_prime(p,b,targets,a)
    prime_target_counts.append((int(p),len(targets)))
    for r in targets: records.append((int(p),int(r),vecs[r]))
print('Q7676_SCAN_PMAX',PMAX,'TARGET_PRIMES',len(prime_target_counts),'TARGET_COUNT',len(records),'MAX_Z',max(z for p,z in prime_target_counts))
r1713=[z for z in records if z[0]==17 and z[1]==13][0]
print('EICHLER_FRAME_VECTOR_17_13',r1713[2],'CANONICAL_EXPECTED',(9,8,9,7))
TR=(-3,4,-9,108); SEP=(1,-28,63,-36)
forms={'transverse':TR,'separator':SEP,'e++':(1,4,9,36),'e--':(1,-4,-9,36),'e+-':(1,4,-9,-36),'e-+':(1,-4,9,-36),'O1':(1,0,0,-36),'O2':(0,4,-9,0),'P1':(1,0,0,36),'P2':(0,4,9,0)}
for name,c in forms.items():
    hits=[(p,r) for p,r,v in records if dotmod(c,v,p)==0]
    print('FORM_HITS',name,len(hits),hits[:20])
eigs=[forms[x] for x in ['e++','e--','e+-','e-+']]
uncovered=[(p,r,v) for p,r,v in records if all(dotmod(c,v,p) for c in eigs)]
print('EIGEN_UNION_COVERED',len(records)-len(uncovered),'OF',len(records),'UNCOVERED',len(uncovered),'FIRST',uncovered[:10])
chars=defaultdict(int)
for p,r,v in records: chars[int(kronecker(-6,p))]+=1
print('CHAR_COUNTS',dict(chars))

def two_plane_cover(odd=True,Hpair=30):
    masks={}; union=0
    pairs=[]
    for aa in range(-Hpair,Hpair+1):
        for bb in range(-Hpair,Hpair+1):
            if aa==bb==0: continue
            c2=canon((aa,bb))
            if c2 not in pairs: pairs.append(c2)
    pairs=list(set(pairs))
    for tid,(p,r,v) in enumerate(records):
        A=(v[0]-36*v[3])%p if odd else (v[0]+36*v[3])%p
        B=(4*v[1]-9*v[2])%p if odd else (4*v[1]+9*v[2])%p
        for c2 in pairs:
            if (c2[0]*A+c2[1]*B)%p==0:
                masks[c2]=masks.get(c2,0)|(1<<tid); union|=1<<tid
    top=sorted(((m.bit_count(),c) for c,m in masks.items()),reverse=True)[:10]
    print('PLANE','odd' if odd else 'even','H',Hpair,'CANDS',len(masks),'COVERED',union.bit_count(),'MAX',top[:10])

two_plane_cover(True,30); two_plane_cover(False,30)

def full_cover(H):
    masks={}; union=0; rng=range(-H,H+1)
    for tid,(p,r,v) in enumerate(records):
        left=defaultdict(list)
        for a0 in rng:
            for a1 in rng: left[(a0*v[0]+a1*v[1])%p].append((a0,a1))
        local=set()
        for a2 in rng:
            for a3 in rng:
                need=(-(a2*v[2]+a3*v[3]))%p
                for a0,a1 in left.get(need,()):
                    if a0==a1==a2==a3==0: continue
                    c=canon((a0,a1,a2,a3))
                    if c and max(map(abs,c))<=H: local.add(c)
        if local: union|=1<<tid
        for c in local: masks[c]=masks.get(c,0)|(1<<tid)
    top=sorted(((m.bit_count(),c,m) for c,m in masks.items()),reverse=True)[:15]
    print('FULL_H',H,'CANDS',len(masks),'COVERED',union.bit_count(),'OF',len(records),'MULTI',sum(m.bit_count()>=2 for m in masks.values()),'TOP',[(n,c) for n,c,m in top])
    uncovered=(1<<len(records))-1; chosen=[]; vals=list(masks.items())
    while uncovered:
        c,m=max(vals,key=lambda cm:(cm[1]&uncovered).bit_count())
        n=(m&uncovered).bit_count()
        if not n: break
        chosen.append((c,n)); uncovered&=~m
    print('GREEDY',H,'COUNT',len(chosen),'UNCOVERED',uncovered.bit_count(),'HEAD',chosen[:15])
    if top: print('COVER_LOWER_BY_MAX',H,(len(records)+top[0][0]-1)//top[0][0])
for H in [4,8,12]: full_cover(H)
print('DONE')
