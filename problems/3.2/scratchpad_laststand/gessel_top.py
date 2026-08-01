from sympy import primerange
from math import comb
def apery(n): return sum(comb(n,k)**2*comb(n+k,k)**2 for k in range(n+1))
def b_mod(k,p):
    if k==0: return 1%p
    bm1,bc=1%p,5%p
    for n in range(1,k):
        A=(34*n**3+51*n**2+27*n+5)%p; B=(n**3)%p; inv=pow((n+1)**3%p,p-2,p)
        bm1,bc=bc,((A*bc-B*bm1)*inv)%p
    return bc
ok=True
for n in range(40,140):
    bn=apery(n)
    for p in primerange(n//2+1,n+1):
        lhs = (b_mod(n-p,p)%p==0)          # p | b_{n-p}
        rhs = (bn % p == 0)                # p | b_n
        pred= (bn % p) == (5*b_mod(n-p,p))%p   # Gessel: b_n = b_1*b_{n-p}
        if lhs!=rhs or not pred:
            ok=False; print(f"MISMATCH n={n} p={p}: lhs={lhs} rhs={rhs} gessel_ok={pred}"); break
    if not ok: break
print("Gessel top-window identity b_n == 5*b_{n-p} mod p AND (p|b_{n-p} <=> p|b_n):", "HOLDS for all n in [40,140)" if ok else "FAILED")
