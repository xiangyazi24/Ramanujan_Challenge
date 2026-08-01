import math, sys
sys.path.insert(0,'.')
from CRON_b1_crosscorr import orbit
import sympy as sp
out=[]
ps=[461,997,1009,2003,3001,5003,7013,10007,15013,20011,30011,40009]
for p in ps:
    N=p-1; D=int(math.isqrt(N)*math.log(N))
    pts=orbit(p)
    R=[0]*(D+1); heavy=[]
    for h in range(1,D+1):
        c=0
        for r in range(0,p-1-h):
            if pts[r]==pts[r+h]: c+=1
        R[h]=c
    mx=max(R[1:]); arg=[h for h in range(1,D+1) if R[h]==mx]
    S=sum(R[1:]); E2=sum(v*v for v in R[1:])
    viol=[h for h in range(1,D+1) if (h%2==0)!=(R[h]%2==1)]
    T0=math.ceil(math.sqrt(N/D))
    mu=[(R[h]-(1 if h%2==0 else 0))//2 for h in range(1,D+1)]
    Wmu=max([t*sum(1 for m in mu if m>=t) for t in range(T0+1, max(mu)+2)] or [0])
    print(f"p={p:6d} D={D:5d} T0={T0:3d} maxR={mx:3d}@h={arg[:3]} maxMu={max(mu)} S_D={S:6d} S/D={S/D:.3f} E2/D={E2/D:.3f} parityViol={len(viol)} W^mu={Wmu}", flush=True)
