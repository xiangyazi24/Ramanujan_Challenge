#!/usr/bin/env python3
"""Verify codex-high's branch reflection law: for the truncated square-root polynomial
(deg D), coefficients satisfy a_{D-j} = (-2|p) a_j mod p. Test tau (chi=+1 classes) and
sigma across classes, p < 500."""
def legendre(a, p): return pow(a % p, (p-1)//2, p)
def apery(N):
    b=[1,5]
    for n in range(1,N):
        num=(2*n+1)*(17*n*n+17*n+5)*b[n]-n**3*b[n-1]
        q,r=divmod(num,(n+1)**3); assert r==0
        b.append(q)
    return b
def sqrt_mod(c,N,p):
    s=[1]+[0]*(N-1); i2=pow(2,p-2,p)
    for n in range(1,N):
        acc=c[n] if n<len(c) else 0
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
def sieve(n):
    s=bytearray([1])*(n+1); s[0]=s[1]=0
    for i in range(2,int(n**.5)+1):
        if s[i]: s[i*i::i]=bytearray(len(s[i*i::i]))
    return [i for i in range(n+1) if s[i]]
b=apery(520)
ok=bad=0
for p in sieve(500):
    if p<7: continue
    bm=[v%p for v in b[:p]]
    chi=legendre(-6,p)
    eps=legendre(-2,p)
    if chi==1:
        D=(p-1)//2
        t=sqrt_mod(bm,D+1,p)
    else:
        D=(p-3)//2
        q=[1,(-34)%p,1]
        t=sqrt_mod(div_mod(bm,q,D+1,p),D+1,p)
    good=all((t[D-j]-eps*t[j])%p==0 for j in range(D+1))
    ok+=good; bad+=not good
    if not good and bad<=3:
        j0=next(j for j in range(D+1) if (t[D-j]-eps*t[j])%p)
        print("FAIL p=",p,"first bad j=",j0, t[D-j0]%p, eps*t[j0]%p)
print(f"branch reflection a_(D-j) = (-2|p) a_j: {ok} primes pass, {bad} fail")
