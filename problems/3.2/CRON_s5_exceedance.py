import math, sys
from collections import defaultdict, Counter
exec(open('CRON_s5_fastprofile.py').read().split('print("p N D T0')[0])
import sympy
ps=[p for p in [2003,5003,10007,20011,50021,100003,200003,300007,500009,700001,1000003]]
print("p D T0 maxMu  hist(mu)  E[maxMu|Pois(.48)]  expected#{mu>=t} for t=3,4,5,6", flush=True)
for p in ps:
    N,D,R=profile(p)
    T0=math.ceil(math.sqrt(N/D))
    mu=[(R[h]-(1 if h%2==0 else 0))//2 for h in range(1,D+1)]
    c=Counter(mu); mm=max(mu)
    mubar=sum(mu)/len(mu)
    # Poisson(mubar) predictions
    import math as m
    def ppois_ge(t,lam):
        s=sum(lam**k/m.factorial(k) for k in range(t)); return 1-m.exp(-lam)*s
    exp_counts={t: round(D*ppois_ge(t,mubar),2) for t in [3,4,5,6,7]}
    obs={t: sum(1 for x in mu if x>=t) for t in [3,4,5,6,7]}
    # E[max] approx: smallest t with D*P(>=t)<1
    tmax=min([t for t in range(1,40) if D*ppois_ge(t,mubar)<1])
    print(f"p={p:8d} D={D:6d} T0={T0:2d} maxMu={mm} mubar={mubar:.3f} predMax={tmax} obs{obs} pred{exp_counts}", flush=True)
