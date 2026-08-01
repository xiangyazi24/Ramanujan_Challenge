#!/usr/bin/env python3
"""Independent check: b_p - 5 ≡ -7 p^2 H2_{p-1} (mod p^5), H2_n = sum 1/j^2."""
from fractions import Fraction as F
def apery(N):
    b=[1,5]
    for n in range(1,N):
        num=(2*n+1)*(17*n*n+17*n+5)*b[n]-n**3*b[n-1]
        q,r=divmod(num,(n+1)**3); assert r==0
        b.append(q)
    return b
def sieve(n):
    s=bytearray([1])*(n+1); s[0]=s[1]=0
    for i in range(2,int(n**.5)+1):
        if s[i]: s[i*i::i]=bytearray(len(s[i*i::i]))
    return [i for i in range(n+1) if s[i]]
primes=[p for p in sieve(90) if p>=7]
b=apery(max(primes)+2)
ok=bad=0
for p in primes:
    p5=p**5
    H2=F(0)
    for j in range(1,p):
        H2+=F(1,j*j)
    # H2 is p-integral (j<p); reduce mod p^3 enough (multiplied by p^2)
    num,den=H2.numerator,H2.denominator
    H2m = num*pow(den,-1,p5)%p5
    lhs=(b[p]-5)%p5
    rhs=(-7)*p*p%p5*H2m%p5
    if lhs==rhs: ok+=1
    else:
        bad+=1
        if bad<=3: print("FAIL p=",p,lhs,rhs)
print(f"b_p - 5 = -7 p^2 H2_(p-1) mod p^5: {ok} verified, {bad} failed (7<=p<=89)")
