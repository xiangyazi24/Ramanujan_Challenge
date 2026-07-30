#!/usr/bin/env python3
"""Q5716: finite-field audit of the exact Legendre Green/CD reduction.

Standard library only.  For the long squared-binomial core
  S_M(z)=sum_{k=0}^M C(M,k)^2 J_{M+k}(z),
J_n(z)=sum_r C(n,r)^2 z^r,
construct the tridiagonal adjoint continuant E and the two Cramer endpoint
numerators N0,NM.  Verify
  E*S + (2M+1)*NM*J_{2M+1} + M*(1-z)^2*N0*J_{M-1}=0
coefficientwise modulo two large primes, and report gcd degrees.
"""
from math import comb

MODS=(1000000007,1000000009)
MS=(199,271,299,320,754)

def trim(a):
    while len(a)>1 and a[-1]==0:a.pop()
    return a

def add(a,b,p,sgn=1):
    n=max(len(a),len(b)); c=[0]*n
    for i in range(n): c[i]=((a[i] if i<len(a) else 0)+sgn*(b[i] if i<len(b) else 0))%p
    return trim(c)

def scale(a,c,p): return trim([(c*x)%p for x in a])

def mul(a,b,p):
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        if x:
            for j,y in enumerate(b):
                c[i+j]=(c[i+j]+x*y)%p
    return trim(c)

def mul_1pz(a,c,p):
    # c*(1+z)*a
    out=[0]*(len(a)+1)
    for i,x in enumerate(a):
        out[i]=(out[i]+c*x)%p; out[i+1]=(out[i+1]+c*x)%p
    return trim(out)

def mul_1mz2(a,c,p):
    # c*(1-z)^2*a
    out=[0]*(len(a)+2)
    for i,x in enumerate(a):
        out[i]=(out[i]+c*x)%p
        out[i+1]=(out[i+1]-2*c*x)%p
        out[i+2]=(out[i+2]+c*x)%p
    return trim([x%p for x in out])

def divmod_poly(a,b,p):
    a=a[:]; trim(a); trim(b); inv=pow(b[-1],p-2,p)
    if len(a)<len(b): return [0],a
    q=[0]*(len(a)-len(b)+1)
    while len(a)>=len(b) and not(len(a)==1 and a[0]==0):
        d=len(a)-len(b); c=a[-1]*inv%p; q[d]=c
        for j,y in enumerate(b): a[d+j]=(a[d+j]-c*y)%p
        trim(a)
    return trim(q),trim(a)

def gcd_poly(a,b,p):
    a=trim(a[:]); b=trim(b[:])
    while not(len(b)==1 and b[0]==0):
        _,r=divmod_poly(a,b,p); a,b=b,r
    inv=pow(a[-1],p-2,p); return scale(a,inv,p)

def J_list(N,p):
    J=[[1]]
    if N==0:return J
    J.append([1,1])
    one_mz2=[1,p-2,1]
    for n in range(1,N):
        # (n+1)J_{n+1}=(2n+1)(1+z)J_n-n(1-z)^2J_{n-1}
        rhs=add(mul_1pz(J[n],2*n+1,p),mul_1mz2(J[n-1],n,p),p,-1)
        J.append(scale(rhs,pow(n+1,p-2,p),p))
    return J

def audit(M,p):
    # left continuants theta and forward transformed RHS eta
    theta_m2=[1]
    theta_m1=scale([1,1],-(2*M+1),p)
    eta_m1=[1]  # w_0=1
    thetas=[theta_m1]
    etas=[eta_m1]
    for k in range(1,M+1):
        bpart=mul_1pz(theta_m1,-(2*M+2*k+1),p)
        cross=mul_1mz2(theta_m2,(M+k)*(M+k),p)
        theta=add(bpart,cross,p,-1)
        wk=comb(M,k)%p; wk=wk*wk%p
        eta=add(scale(theta_m1,wk,p),scale(eta_m1,M+k,p),p,-1)
        theta_m2,theta_m1=theta_m1,theta
        eta_m1=eta
        thetas.append(theta); etas.append(eta)
    E=theta_m1; NM=eta_m1

    # right continuants phi and reverse transformed RHS zeta
    phi_p2=[1]
    phi_p1=scale([1,1],-(4*M+1),p)  # b_M
    zeta_p1=[1]                     # w_M=1
    for k in range(M-1,-1,-1):
        bpart=mul_1pz(phi_p1,-(2*M+2*k+1),p)
        cross=mul_1mz2(phi_p2,(M+k+1)*(M+k+1),p)
        phi=add(bpart,cross,p,-1)
        wk=comb(M,k)%p; wk=wk*wk%p
        zeta=add(scale(phi_p1,wk,p),mul_1mz2(zeta_p1,M+k+1,p),p,-1)
        phi_p2,phi_p1=phi_p1,phi
        zeta_p1=zeta
    assert phi_p1==E
    N0=zeta_p1

    J=J_list(2*M+1,p)
    S=[0]
    for k in range(M+1):
        w=comb(M,k)%p; w=w*w%p
        S=add(S,scale(J[M+k],w,p),p)
    lhs=mul(E,S,p)
    lhs=add(lhs,scale(mul(NM,J[2*M+1],p),2*M+1,p),p)
    lhs=add(lhs,mul_1mz2(mul(N0,J[M-1],p),M,p),p)
    assert all(x%p==0 for x in lhs), (M,p,len(lhs),next((i,x) for i,x in enumerate(lhs) if x%p))

    g0=gcd_poly(E,N0,p); gm=gcd_poly(E,NM,p); gall=gcd_poly(g0,NM,p)
    # held-out coefficients: identity was checked fully; print five spread positions
    held=sorted(set([0,1,M//7,M//3,M//2,2*M-1,2*M]))
    print('M',M,'p',p,'degE',len(E)-1,'gcd(E,N0)',len(g0)-1,'gcd(E,NM)',len(gm)-1,'gcd_all',len(gall)-1,'held',held)
    # exact value at z=1 predicted by continuant
    e1=sum(E)%p
    pred=1
    for k in range(M+1): pred=pred*(2*(2*M+2*k+1))%p
    assert e1==((-1)**(M+1)*pred)%p

for p in MODS:
    for M in MS:audit(M,p)
print('PASS')
