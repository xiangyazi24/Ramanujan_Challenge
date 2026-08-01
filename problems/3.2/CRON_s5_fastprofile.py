import math, sys
from collections import defaultdict
def orbit_vals(p):
    """projective orbit as ints: value = c/b mod p, or p for b=0 (infinity)."""
    vals=[]
    b0,b1,c0,c1 = 1%p,5%p,0,1
    def key(x,y): return (y*pow(x,p-2,p))%p if x else p
    vals.append(key(b0,c0)); vals.append(key(b1,c1))
    bm2,bm1,cm2,cm1=b0,b1,c0,c1
    for n in range(2,p-1):
        A=(34*n*n*n-51*n*n+27*n-5)%p
        B=((n-1)**3)%p
        inv=pow((n*n*n)%p,p-2,p)
        bn=(A*bm1-B*bm2)*inv%p
        cn=(A*cm1-B*cm2)*inv%p
        vals.append(key(bn,cn))
        bm2,bm1,cm2,cm1=bm1,bn,cm1,cn
    return vals
def profile(p):
    N=p-1; D=int(math.isqrt(N)*math.log(N))
    vals=orbit_vals(p)
    fib=defaultdict(list)
    for r,v in enumerate(vals): fib[v].append(r)
    R=[0]*(D+2)
    for v,pos in fib.items():
        m=len(pos)
        if m<2: continue
        for i in range(m):
            pi_=pos[i]
            for j in range(i+1,m):
                d=pos[j]-pi_
                if d>D: break
                if pos[j] <= N-1:   # nonwrapping: r+h <= p-2
                    R[d]+=1
    return N,D,R
print("p N D T0 maxR argmax maxMu S_D S/D E2/D parityViol Wmu", flush=True)
import sympy
ps=[10007,30011,100003,300007,1000003]
for p in ps:
    N,D,R=profile(p)
    T0=math.ceil(math.sqrt(N/D))
    mx=max(R[1:D+1]); arg=[h for h in range(1,D+1) if R[h]==mx][:3]
    S=sum(R[1:D+1]); E2=sum(x*x for x in R[1:D+1])
    viol=[h for h in range(1,D+1) if (h%2==0)!=(R[h]%2==1)]
    mu=[(R[h]-(1 if h%2==0 else 0))//2 for h in range(1,D+1)]
    mm=max(mu)
    Wmu=max([t*sum(1 for x in mu if x>=t) for t in range(T0+1,mm+2)] or [0])
    print(f"{p} {N} {D} {T0} maxR={mx}@{arg} maxMu={mm} S_D={S} S/D={S/D:.3f} E2/D={E2/D:.3f} viol={len(viol)} Wmu={Wmu}", flush=True)
