from sympy import primerange
import math
def Zp(p):
    b0,b1=1%p,5%p; Z=[]
    if b0==0: Z.append(0)
    if b1==0: Z.append(1)
    bm1,bc=b0,b1
    for n in range(1,p-1):
        A=(34*n**3+51*n**2+27*n+5)%p; B=(n**3)%p; inv=pow((n+1)**3%p,p-2,p)
        nb=((A*bc-B*bm1)*inv)%p
        if nb==0: Z.append(n+1)
        bm1,bc=bc,nb
    return Z
def b_mod(k,p):
    if k==0: return 1%p
    bm1,bc=1%p,5%p
    for n in range(1,k):
        A=(34*n**3+51*n**2+27*n+5)%p; B=(n**3)%p; inv=pow((n+1)**3%p,p-2,p)
        bm1,bc=bc,((A*bc-B*bm1)*inv)%p
    return bc
X=600
# LHS: sum_{n<=X} TOP(n) where TOP(n)=#{p in (n/2,n]: p|b_n}
lhs=0
for n in range(2,X+1):
    for p in primerange(n//2+1,n+1):
        r=n-p
        if b_mod(r,p)%p==0 and not (r==0):  # r=0 -> b_0=1 never 0
            lhs+=1
# RHS: sum_{p<=X} |Z_p ∩ [1,p-1]| restricted so that n=p+r <= X
rhs=0; tot=0
for p in primerange(2,X+1):
    Z=Zp(p); tot+=len(Z)
    rhs+=sum(1 for r in Z if 1<=r and p+r<=X)
print(f"X={X}: sum_n TOP(n) = {lhs}   sum_p #{{r in Z_p : p+r<=X}} = {rhs}   -> identity {'HOLDS' if lhs==rhs else 'FAILS'}")
print(f"sum_{{p<=X}} |Z_p| = {tot}  vs pi(X)={len(list(primerange(2,X+1)))}  vs X={X}  vs record bound sum 3p^(2/3) = {sum(3*p**(2/3) for p in primerange(2,X+1)):.0f}")
