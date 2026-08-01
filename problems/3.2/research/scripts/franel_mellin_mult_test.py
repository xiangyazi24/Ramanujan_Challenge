# Verify: (1) fiber discriminant of phi(x)=t is q(t)=t^2-34t+1 (up to square factor)
# (2) sum_x H(x)^2 phi(x)^{-r} = sum_t (1+chi(q(t))) A_p(t) t^{-r}  over valid ranges
# (3) M(r) = sum_{t in F_p*} A_p(t) t^{-r} = -(b_r + [r==0-ish corrections])
from math import comb
for p in [13, 29, 37]:
    b=[1,5]
    for n in range(1,p):
        b.append(((34*n**3+51*n**2+27*n+5)*b[n]-n**3*b[n-1])*pow((n+1)**3,p-2,p)%p)
    A=b[:p]
    f=[sum(comb(n,k)**3 for k in range(n+1))%p for n in range(p)]
    chi=lambda a: pow(a%p,(p-1)//2,p) if a%p else 0
    chi_s=lambda a: {0:0,1:1,p-1:-1}[chi(a)]
    Aval=lambda t: sum(A[n]*pow(t,n,p) for n in range(p))%p
    Hval=lambda x: sum(f[n]*pow(x,n,p) for n in range(p))%p
    # (1) fiber count vs 1+chi(q(t))
    from collections import Counter
    fib=Counter()
    for x in range(p):
        if (1+x)%p==0: continue
        fib[x*(1-8*x)%p*pow(1+x,p-2,p)%p]+=1
    ok1=all(fib.get(t,0)==1+chi_s(t*t-34*t+1) for t in range(p) if (t*t-34*t+1)%p)
    # (2)+(3) for a few r
    res=[]
    for r in [1,2,3,(p-1)//4,(p-1)//2]:
        LHS=sum(Hval(x)*Hval(x)*pow(x*(1-8*x)%p*pow(1+x,p-2,p)%p,p-1-r,p) for x in range(p) if (1+x)%p and x*(1-8*x)%p)%p
        RHS=sum((1+chi_s(t*t-34*t+1))*Aval(t)*pow(t,p-1-r,p) for t in range(1,p))%p
        M=sum(Aval(t)*pow(t,p-1-r,p) for t in range(1,p))%p
        res.append((r, LHS==RHS, (M + b[r] + b[p-1+r-(p-1)] if False else (M + b[r]))%p ))
    # M(r) should be -(b_r + b_{r+p-1 term}) : n==r and n==r+(p-1)<p i.e. r=0 only; also n=r-(p-1)... 
    print(f"p={p}: fiberdisc==q: {ok1}; [(r, LHS==RHS, M+b_r mod p)]:", res)
