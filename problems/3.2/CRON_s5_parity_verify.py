import math, sys
sys.setrecursionlimit(10000)
from CRON_b1_crosscorr import orbit

def apery_b(p):
    # b_j mod p, j=0..~D
    bs=[1%p,5%p]
    for n in range(2,400):
        A=(34*n**3-51*n**2+27*n-5)%p; B=((n-1)**3)%p
        inv=pow(n**3%p,-1,p)
        bs.append((A*bs[-1]-B*bs[-2])*inv%p)
    return bs

for p in [1009,3001,10007]:
    N=p-1; D=int(math.isqrt(N)*math.log(N))
    pts=orbit(p); bs=apery_b(p)
    R={}; bases={}
    for h in range(1,D+1):
        bl=[r for r in range(0,p-1-h) if pts[r]==pts[r+h]]
        R[h]=len(bl); bases[h]=bl
    S=sum(R.values()); E2=sum(v*v for v in R.values())
    mx=max(R.values()); argmx=[h for h in R if R[h]==mx]
    # parity law check
    viol_even=[h for h in R if h%2==0 and R[h]%2==0]
    odd_odd=[h for h in R if h%2==1 and R[h]%2==1]
    odd_odd_expl=[h for h in odd_odd if bs[(h-1)//2]==0]  # p | b_{(h-1)/2}
    # center collision check for even h
    center_missing=[h for h in R if h%2==0 and ((p-1-h)//2 not in bases[h])]
    # mirror closure of base sets
    mirror_viol=0; mirror_viol_detail=[]
    for h in R:
        s=set(bases[h])
        for r in s:
            m=p-1-h-r
            if m!=r and m not in s and 0<=m<=p-2-h:
                mirror_viol+=1; mirror_viol_detail.append((h,r))
    # boundary cases: r=0 root or r=p-1-h root
    bnd=[(h,r) for h in R for r in bases[h] if r==0 or r==p-1-h]
    # C(t) support (same-lag pairs by base offset)
    from collections import Counter
    Ct=Counter()
    for h in R:
        bl=bases[h]
        for i in range(len(bl)):
            for j in range(len(bl)):
                if i!=j: Ct[(bl[j]-bl[i])%p]+=1
    # histogram by parity
    he=Counter(R[h] for h in R if h%2==0); ho=Counter(R[h] for h in R if h%2==1)
    print(f"p={p} D={D} S_D={S} E2={E2} E2/D={E2/D:.2f} max={mx} argmax={argmx[:6]} argmax_parity={[h%2 for h in argmx[:6]]}")
    print(f"  parity violations: even-h with R even: {viol_even[:10]}{'...' if len(viol_even)>10 else ''} (n={len(viol_even)})")
    print(f"  odd-h with R odd: n={len(odd_odd)}, of which explained by p|b_(h-1)/2: {len(odd_odd_expl)}; unexplained: {[h for h in odd_odd if h not in odd_odd_expl][:10]}")
    print(f"  center-collision missing (even h): {center_missing[:10]} (n={len(center_missing)})")
    print(f"  mirror-closure violations: {mirror_viol} detail {mirror_viol_detail[:6]}")
    print(f"  boundary roots (r=0 or r=p-1-h): {bnd[:6]} (n={len(bnd)})")
    print(f"  even-h hist {dict(sorted(he.items()))}")
    print(f"  odd-h hist  {dict(sorted(ho.items()))}")
    nz=[t for t,c in Ct.items() if c>0]
    print(f"  C(t): support size={len(nz)} of {p-1}, total offdiag={sum(Ct.values())}, top5={Ct.most_common(5)}")
