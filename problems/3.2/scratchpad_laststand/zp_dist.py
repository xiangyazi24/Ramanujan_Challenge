from sympy import primerange
from collections import Counter
import math
def Zp_count(p):
    b0,b1=1%p,5%p; cnt=(1 if b0==0 else 0)+(1 if b1==0 else 0)
    bm1,bc=b0,b1
    for n in range(1,p-1):
        A=(34*n**3+51*n**2+27*n+5)%p; B=(n**3)%p; inv=pow((n+1)**3%p,p-2,p)
        nb=((A*bc-B*bm1)*inv)%p
        if nb==0: cnt+=1
        bm1,bc=bc,nb
    return cnt
C=Counter(); tot=0; N=0
for p in primerange(5,8000):
    z=Zp_count(p); C[z]+=1; tot+=z; N+=1
print(f"primes={N} mean={tot/N:.4f}  var={sum((k-tot/N)**2*v for k,v in C.items())/N:.4f}")
print("k : observed  | Poisson(mean) predicted | S_inf fixed-pt (Poisson(1))")
lam=tot/N
for k in range(0,6):
    pois=N*math.exp(-lam)*lam**k/math.factorial(k)
    pois1=N*math.exp(-1)/math.factorial(k)
    print(f"{k} : {C.get(k,0):5d}    | {pois:8.1f}   | {pois1:8.1f}")
