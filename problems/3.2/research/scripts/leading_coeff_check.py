#!/usr/bin/env python3
"""Independent check: leading coefficient of the truncated sqrt branch = (-2|p) mod p.
tau_{(p-1)/2} for chi(-6|p)=+1... per cron: tau_{(p-1)/2}: 63/63 (all p?); sigma_{(p-3)/2} all four classes."""
def legendre(a,p): 
    v = pow(a % p, (p-1)//2, p)
    return v if v <= 1 else v - p
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
b=apery(1050)
okt=badt=oks=bads=0
for p in sieve(1000):
    if p<7: continue
    bm=[v%p for v in b[:p]]
    eps=legendre(-2,p) % p
    tau=sqrt_mod(bm,(p-1)//2+1,p)
    if tau[(p-1)//2]==eps: okt+=1
    else: badt+=1; print("tau FAIL",p,tau[(p-1)//2],eps) if badt<4 else None
    q=[1,(-34)%p,1]
    sig=sqrt_mod(div_mod(bm,q,(p-3)//2+1,p),(p-3)//2+1,p)
    if sig[(p-3)//2]==eps: oks+=1
    else: bads+=1; print("sig FAIL",p,sig[(p-3)//2],eps) if bads<4 else None
print(f"tau_((p-1)/2) = (-2|p): {okt} ok, {badt} fail;  sigma_((p-3)/2) = (-2|p): {oks} ok, {bads} fail (p<1000)")
