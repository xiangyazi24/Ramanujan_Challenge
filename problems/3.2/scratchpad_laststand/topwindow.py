from sympy import primerange
import math
def b_mod(k,p):
    # b_k mod p via recurrence (k < p so denominators invertible)
    if k==0: return 1%p
    bm1,bc=1%p,5%p
    for n in range(1,k):
        A=(34*n**3+51*n**2+27*n+5)%p; B=(n**3)%p; inv=pow((n+1)**3%p,p-2,p)
        bm1,bc=bc,((A*bc-B*bm1)*inv)%p
    return bc
for n in [500,1000,2000,4000,8000]:
    hits=[]; np_=0
    for p in primerange(n//2+1,n+1):
        np_+=1
        r=n-p
        if r>=0 and b_mod(r,p)%p==0: hits.append((p,r))
    print(f"n={n}: primes in (n/2,n]={np_}, hits p|b_(n-p): {len(hits)} {hits[:4]}  (need ~n/log n={n/math.log(n):.0f} to break o(n))")
