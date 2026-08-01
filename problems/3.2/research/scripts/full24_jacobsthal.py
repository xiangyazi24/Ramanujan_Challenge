#!/usr/bin/env python3
"""Complete the mod-24 table for quarter-point behavior of tau=sqrt(F), sigma=sqrt(F/q),
and hunt Jacobsthal value laws: value at/near quarter point vs (x,y) with p = x^2+6y^2
(classes 1,7 mod 24) or p = 2x^2+3y^2 (classes 5,11)."""
from fractions import Fraction as F

def apery(N):
    b=[1,5]
    for n in range(1,N):
        num=(2*n+1)*(17*n*n+17*n+5)*b[n]-n**3*b[n-1]
        q,r=divmod(num,(n+1)**3); assert r==0
        b.append(q)
    return b

def sqrt_mod(coeffs,N,p):
    s=[1]+[0]*(N-1); i2=pow(2,p-2,p)
    for n in range(1,N):
        acc=coeffs[n] if n<len(coeffs) else 0
        t=sum(s[i]*s[n-i] for i in range(1,n))%p
        s[n]=(acc-t)%p*i2%p
    return s

def div_mod(a,den,N,p):
    c=[0]*N; d0=pow(den[0],p-2,p)
    for n in range(N):
        t=a[n] if n<len(a) else 0
        t=(t-sum(den[i]*c[n-i] for i in range(1,min(n,len(den)-1)+1)))%p
        c[n]=t*d0%p
    return c

def rep(p, A, B):
    # p = A x^2 + B y^2
    x=1
    while A*x*x<=p:
        rem=p-A*x*x
        if rem%B==0:
            y2=rem//B
            y=int(y2**.5)
            for yy in (y-1,y,y+1):
                if yy>0 and yy*yy==y2: return x,yy
        x+=1
    return None

def sieve(n):
    s=bytearray([1])*(n+1); s[0]=s[1]=0
    for i in range(2,int(n**.5)+1):
        if s[i]: s[i*i::i]=bytearray(len(s[i*i::i]))
    return [i for i in range(n+1) if s[i]]

b=apery(1100)
from collections import defaultdict
table=defaultdict(lambda: [0,0])
vals=defaultdict(list)
for p in sieve(2000):
    if p<7: continue
    m=p%24
    bm=[v%p for v in b[:p]]
    J=(p-1)//4 if p%4==1 else (p-3)//4
    N=J+3
    tau=sqrt_mod(bm,N,p)
    q=[1,(-34)%p,1]
    sig=sqrt_mod(div_mod(bm,q,N,p),N,p)
    zt=tau[J]==0; zs=sig[J]==0
    table[m][0]+= zt; table[m][1]+= zs
    # value hunt: principal classes 1,7: p=x^2+6y^2; classes 5,11: p=2x^2+3y^2
    if m in (1,7):
        r=rep(p,1,6)
        if r: vals[m].append((p,tau[J],sig[J],r))
    elif m in (5,11):
        r=rep(p,2,3)
        if r: vals[m].append((p,tau[J],sig[J],r))
cnt=defaultdict(int)
for p in sieve(2000):
    if p>=7: cnt[p%24]+=1
print("class: #primes  tau-quarter-zeros  sigma-quarter-zeros")
for m in (1,5,7,11,13,17,19,23):
    print(f"  {m:>2}: {cnt[m]:>3}   {table[m][0]:>3}   {table[m][1]:>3}")
# Jacobsthal value test for class 1 (tau nonzero): tau_J vs 2x mod p? try ratios
print("\nvalue-law hunt (class, p, tau_J, sigma_J, (x,y)) + candidate ratios tau_J/x, tau_J/y mod p:")
for m in (1,5):
    for (p,tv,sv,(x,y)) in vals[m][:6]:
        rx = tv*pow(x,p-2,p)%p if x%p else '-'
        ry = tv*pow(y,p-2,p)%p if y%p else '-'
        print(f"  m={m} p={p}: tau_J={tv} sigma_J={sv} (x,y)=({x},{y})  tau/x={rx} tau/y={ry}")
