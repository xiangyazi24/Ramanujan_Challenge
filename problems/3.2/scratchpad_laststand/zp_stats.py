from sympy import primerange
import math
def Zp_count(p):
    # count r in [0,p-1] with b_r == 0 mod p
    b0,b1=1%p,5%p; cnt=(1 if b0==0 else 0)+(1 if b1==0 else 0)
    bm1,bc=b0,b1
    for n in range(1,p-1):
        A=(34*n**3+51*n**2+27*n+5)%p; B=(n**3)%p; inv=pow((n+1)**3%p,p-2,p)
        nb=((A*bc-B*bm1)*inv)%p
        if nb==0: cnt+=1
        bm1,bc=bc,nb
    return cnt
tot=0; s_logp=0; np_=0
for p in primerange(5,4000):
    z=Zp_count(p); tot+=z; s_logp+=z*math.log(p); np_+=1
    X=4000
print(f"primes counted={np_}, mean |Z_p|={tot/np_:.3f}")
print(f"sum_p log p * |Z_p| = {s_logp:.0f}  vs X=4000  ratio={s_logp/4000:.3f}")
