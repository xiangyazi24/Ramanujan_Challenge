#!/usr/bin/env python3
"""Independent exact audit of the Q3573 Gamma_47 carrier.

Standalone: no project computation imports and no Fable script reuse.
Reconstructs Apéry numbers, sampled forward differences, Newton selectors,
primitive three-row C* carriers, legal q=6 windows, and Gamma_47 from the
published formulas.  All arithmetic is integer-exact; stdlib only.
"""

from collections import Counter, defaultdict
from math import comb, gcd, isqrt, prod

H = 0
DELTAS = (0, 1, 2)
FIXTURES = ((19,8),(97,25),(139,61),(181,19),(293,47))
COLLISION_M = 2932
COLLISION_PAIRS = ((439,298),(443,274))
CENSUS_MINP = 43
CENSUS_LIMIT = 1000
# Exact H=0 far-edge values, derived from the slope-2 formula.
CAMPAIGN_F2 = (0,7,2,3,4,5)  # indexed by m mod 6


def apery_numbers(limit):
    if limit == 0: return [1]
    b=[1,5]
    for n in range(1,limit):
        num=(34*n**3+51*n**2+27*n+5)*b[n]-n**3*b[n-1]
        den=(n+1)**3
        assert num%den==0
        b.append(num//den)
    return b[:limit+1]


def primes_upto(limit):
    a=bytearray(b"\x01")*(limit+1)
    a[0:2]=b"\x00\x00"
    for p in range(2,isqrt(limit)+1):
        if a[p]: a[p*p:limit+1:p]=b"\x00"*(((limit-p*p)//p)+1)
    return [p for p in range(2,limit+1) if a[p]]


def sigma(m,h,s):
    out=[x for x in range(h+1,h+s+1) if (m+1+x)%s==0]
    assert len(out)==1
    return out[0]


def windows_from_f2(m,h,F2):
    s4=sigma(m,h,4); s7=sigma(m,h,7)
    X4=(m+1+s4)//4; X7=(m+1+s7)//7
    Phi=(m-F2)//6
    A4=max(1,2*Phi-X4+1); B4=min(X7-1,1+(m-s4)//4)
    A7=max(1,Phi-X7+1); B7=min(X7-1,1+(m-s7)//7)
    return F2,s4,X4,s7,X7,Phi,A4,B4,A7,B7


def legal_windows(m,h,delta2):
    s2=sigma(m,h,2)
    F2=3*s2+2*h-3+2*delta2
    return windows_from_f2(m,h,F2)


def campaign_windows(m):
    assert H==0
    F2=CAMPAIGN_F2[m%6]
    w=windows_from_f2(m,0,F2)
    assert any(legal_windows(m,0,d)[0]==F2 for d in DELTAS)
    return w


def neg_binom(X,k):
    return (-1 if k&1 else 1)*comb(X+k-1,k)


def forward_differences(values):
    row=list(values); out=[]
    while row:
        out.append(row[0])
        row=[row[i+1]-row[i] for i in range(len(row)-1)]
    return out


def slope_carriers(m,s,sg,X,A,B,b):
    assert A<=B-2,(m,s,A,B)
    sampled=[b[sg+s*t] for t in range(B)]
    d=forward_differences(sampled)
    U=[neg_binom(X,k) for k in range(B+1)]
    G=0
    for k in range(A,B+1): G=gcd(G,abs(U[k]))
    e=[u//G for u in U]
    eg=0
    for k in range(A,B+1): eg=gcd(eg,abs(e[k]))
    assert eg==1
    T=[0]*(B+1)
    for k in range(B): T[k+1]=T[k]+U[k]*d[k]
    out=[]
    for k in range(A,B-1):
        Ak=X+2*k+1
        if Ak&1:
            aa=(Ak+1)//2; bb=-(Ak-1)//2; gs=1
        else:
            aa=-1; bb=1; gs=2
        assert aa*Ak+bb*(Ak+2)==gs
        Vk=(X+k)*T[k]+(k+1)*T[k+1]
        V1=(X+k+1)*T[k+1]+(k+2)*T[k+2]
        C=aa*Vk+bb*V1
        Z=((aa*(k+1)+bb*(Ak+2))*e[k]*d[k]
           +bb*(k+2)*e[k+1]*d[k+1])
        assert C==gs*T[k]+G*Z
        out.append((k,C))
    return G,out


def gamma_record(m,win,b):
    F2,s4,X4,s7,X7,Phi,A4,B4,A7,B7=win
    G4,C4=slope_carriers(m,4,s4,X4,A4,B4,b)
    G7,C7=slope_carriers(m,7,s7,X7,A7,B7,b)
    gamma=0
    for _,z in C4+C7: gamma=gcd(gamma,abs(z))
    return dict(m=m,F2=F2,Phi=Phi,s4=s4,X4=X4,A4=A4,B4=B4,G4=G4,
                s7=s7,X7=X7,A7=A7,B7=B7,G7=G7,
                count4=len(C4),count7=len(C7),gamma=gamma)


def gamma_campaign(m,b): return gamma_record(m,campaign_windows(m),b)
def gamma_delta(m,d,b):
    r=gamma_record(m,legal_windows(m,H,d),b); r['delta2']=d; return r


def target_pairs_at_level(m,b,primes):
    return [(p,m-6*p) for p in primes
            if 0<=m-6*p<p and m//p==6 and b[m-6*p]%p==0]


def vp(n,p):
    e=0
    while n and n%p==0: n//=p; e+=1
    return e


def fixture_audit():
    max_m=max([6*p+j for p,j in FIXTURES]+[COLLISION_M])
    b=apery_numbers(max_m)
    ps=primes_upto(max_m)
    print('Q3574 INDEPENDENT GAMMA47 AUDIT')
    print('implementation=standalone exact recurrence/Newton/Bezout; no project imports')
    for p,j in FIXTURES:
        m=6*p+j; targets=target_pairs_at_level(m,b,ps); R=prod(q for q,_ in targets)
        r=gamma_campaign(m,b); eps=r['gamma']//R
        assert r['gamma']%R==0
        print(f'FIXTURE p={p} j={j} m={m} b_j_mod_p={b[j]%p} vp_bj={vp(b[j],p)}')
        print(f'  targets={targets} target_product={R}')
        print('  CAMPAIGN F2={F2} I4=[{A4},{B4}] I7=[{A7},{B7}] count=({count4},{count7}) gamma={gamma} epsilon={eps}'.format(eps=eps,**r))
        for d in DELTAS:
            rr=gamma_delta(m,d,b); ee=rr['gamma']//R
            assert rr['gamma']%R==0
            print('  delta2={delta2} F2={F2} gamma={gamma} epsilon={eps}'.format(eps=ee,**rr))
    print(f'COLLISION m={COLLISION_M}')
    for p,j in COLLISION_PAIRS:
        print(f'  p={p} j={j} b_j_mod_p={b[j]%p} vp_bj={vp(b[j],p)}')
    targets=target_pairs_at_level(COLLISION_M,b,ps); R=prod(p for p,_ in targets)
    r=gamma_campaign(COLLISION_M,b); eps=r['gamma']//R
    print(f'  targets={targets} target_product={R}')
    print('  CAMPAIGN F2={F2} I4=[{A4},{B4}] I7=[{A7},{B7}] count=({count4},{count7}) gamma={gamma} epsilon={eps}'.format(eps=eps,**r))
    for d in DELTAS:
        rr=gamma_delta(COLLISION_M,d,b); ee=rr['gamma']//R
        print('  delta2={delta2} F2={F2} gamma={gamma} epsilon={eps}'.format(eps=ee,**rr))
    assert r['gamma']==5*439*443
    print(f'CLAIMED_COLLISION_VALUE 5*439*443={5*439*443}')


def census_audit():
    # Discover target pairs with base prime p<1000, p>=43, using exact core filters.
    base_ps=[p for p in primes_upto(CENSUS_LIMIT-1) if p>=CENSUS_MINP]
    smallb=apery_numbers(CENSUS_LIMIT-1)
    pairs=[]
    for p in base_ps:
        for j in range(p):
            if smallb[j]%p: continue
            m=6*p+j
            F2,s4,X4,s7,X7,Phi,A4,B4,A7,B7=campaign_windows(m)
            if X7<=p<=Phi and F2<=j<=p-1-s7:
                pairs.append((p,j,m))
    by_m=defaultdict(list)
    for p,j,m in pairs: by_m[m].append((p,j))
    levels=sorted(by_m)
    maxm=max(levels)
    b=apery_numbers(maxm)
    # Target scan must cover the entire q=6 interval, not only p<1000.
    allps=primes_upto(maxm//6+2)
    dist=Counter(); rows=[]
    for idx,m in enumerate(levels,1):
        r=gamma_campaign(m,b); gamma=r['gamma']
        targets=target_pairs_at_level(m,b,allps)
        R=prod(p for p,_ in targets)
        assert gamma%R==0,(m,gamma,R,targets)
        eps=gamma//R
        dist[eps]+=1
        rows.append((eps,m,tuple(targets),gamma,R,tuple(by_m[m])))
        if idx%25==0: print(f'  LEVEL_PROGRESS {idx}/{len(levels)}')
    over=[x for x in rows if x[0]>25]
    print(f'CENSUS_PAIR_EVENTS={len(pairs)}')
    print(f'CENSUS_DISTINCT_LEVELS={len(levels)}')
    print(f'CENSUS_COLLISION_LEVELS={[(m,tuple(v)) for m,v in sorted(by_m.items()) if len(v)>1]}')
    print(f'CENSUS_EPS_SET={sorted(dist)}')
    print(f'CENSUS_EPS_DIST={dict(sorted(dist.items()))}')
    print(f'CENSUS_MAX_EPS={max(dist)}')
    print(f'CENSUS_OVER_25={len(over)}')
    print(f'CENSUS_OVER_25_ROWS={over}')
    assert len(levels)==143


def main():
    fixture_audit()
    census_audit()

if __name__=='__main__': main()
