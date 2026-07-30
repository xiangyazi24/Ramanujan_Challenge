#!/usr/bin/env python3
"""Exact/modular audit for Q5723.

- computes the actual gcd content of the 14 y/z ray classes (21 rays);
- computes content across short local state windows;
- guesses the specialized minimal P-recurrence of the full first-cell
  correction F_M(r) over two large finite fields, with held-out checks.

Standard library only.
"""
from math import comb, gcd

MODS=(1000000007,1000000009)
CLASSES=[]
for x, vals in [(-1,(-1,0,1)),(0,(-1,0,1)),(1,(0,1))]:
    for y in vals:
        for z in vals:
            if y>z: continue
            if (x,y,z)==(0,0,0): continue
            CLASSES.append((x,y,z,1 if y==z else 2))
assert len(CLASSES)==14 and sum(m for *_,m in CLASSES)==21

def C(n,k): return comb(n,k) if 0<=k<=n else 0

def ray_exact(M,r,kappa):
    d=M-r; x,y,z=kappa
    s=0
    for t in range(M+1):
        s += C(M,t)*C(M,t-d*x)*C(2*M-t,M-d*y)*C(2*M-t,M-d*z)
    return s

def state_content(M,r,width=0):
    g=0
    for u in range(width+1):
        rr=r+u
        if not (0<=rr and 2*rr<M): break
        for x,y,z,m in CLASSES:
            g=gcd(g,ray_exact(M,rr,(x,y,z)))
    return abs(g)

def cyclic_factor(M,r): return M//gcd(M,r) if r else 1

class ModBinom:
    def __init__(self,N,p):
        self.p=p; self.fact=[1]*(N+1); self.ifact=[1]*(N+1)
        for i in range(1,N+1): self.fact[i]=self.fact[i-1]*i%p
        self.ifact[N]=pow(self.fact[N],p-2,p)
        for i in range(N,0,-1): self.ifact[i-1]=self.ifact[i]*i%p
    def C(self,n,k):
        if k<0 or k>n:return 0
        return self.fact[n]*self.ifact[k]%self.p*self.ifact[n-k]%self.p

def F_sequence(M,p):
    B=ModBinom(2*M,p); Cp=B.C
    b=0
    for k in range(M+1):
        q=Cp(M,k)*Cp(M+k,k)%p
        b=(b+q*q)%p
    out=[]
    for r in range((M-1)//2+1):
        d=M-r; total=0
        for t in range(M+1):
            A=Cp(M,t)
            X=(Cp(M,t-d)+A+Cp(M,t+d))%p
            N=2*M-t
            Z=(Cp(N,M-d)+Cp(N,M)+Cp(N,M+d))%p
            total=(total+A*X%p*Z%p*Z)%p
        out.append((total-b)%p)
    return out

def null_vector(mat,p):
    if not mat:return None
    A=[row[:] for row in mat]; m=len(A); n=len(A[0]); piv=[]; rr=0
    for c in range(n):
        q=next((i for i in range(rr,m) if A[i][c]%p),None)
        if q is None: continue
        A[rr],A[q]=A[q],A[rr]
        inv=pow(A[rr][c]%p,p-2,p)
        A[rr]=[(x*inv)%p for x in A[rr]]
        for i in range(m):
            if i!=rr and A[i][c]%p:
                z=A[i][c]%p
                A[i]=[(x-z*y)%p for x,y in zip(A[i],A[rr])]
        piv.append(c); rr+=1
        if rr==m: break
    free=[c for c in range(n) if c not in set(piv)]
    if not free:return None
    f=free[-1]; v=[0]*n; v[f]=1
    for i in range(len(piv)-1,-1,-1):
        c=piv[i]; v[c]=(-sum(A[i][j]*v[j] for j in free))%p
    return v

def check_relation(seq,o,d,v,p,start=0):
    for r in range(start,len(seq)-o):
        s=0; q=0
        for i in range(o+1):
            rp=1
            for k in range(d+1):
                s=(s+v[q]*rp%p*seq[r+i])%p; q+=1; rp=rp*r%p
        if s%p:return False,r,s%p
    return True,None,None

def guess(seq,p,max_order=38,max_degree=10):
    E=len(seq)-1
    for o in range(1,max_order+1):
        for d in range(max_degree+1):
            n=(o+1)*(d+1); neq=len(seq)-o
            if neq<n+12: continue
            train=max(n+6,int(neq*.72)); train=min(train,neq-6)
            mat=[]
            for r in range(train):
                row=[]
                for i in range(o+1):
                    rp=1
                    for k in range(d+1): row.append(seq[r+i]*rp%p); rp=rp*r%p
                mat.append(row)
            v=null_vector(mat,p)
            if v is None:continue
            ok,_,_=check_relation(seq,o,d,v,p,train)
            if ok:
                # insist one-dimensional nullspace heuristically by adding all equations
                full=[]
                for r in range(neq):
                    row=[]
                    for i in range(o+1):
                        rp=1
                        for k in range(d+1):row.append(seq[r+i]*rp%p);rp=rp*r%p
                    full.append(row)
                vv=null_vector(full,p)
                if vv is not None and check_relation(seq,o,d,vv,p)[0]:return o,d,vv
    return None

def summarize_vector(o,d,v,p):
    # leading/trailing coefficient roots at small affine candidates via evaluation
    def poly(i,r):
        s=0
        for k in range(d,-1,-1):s=(s*r+v[i*(d+1)+k])%p
        return s
    roots0=[a for a in range(-20,21) if poly(0,a)%p==0]
    rootsR=[a for a in range(-20,21) if poly(o,a)%p==0]
    return roots0,rootsR

def main():
    print('CLASSES',CLASSES)
    cases=[(199,90),(271,120),(299,140),(320,159),(320,142),(320,128),(754,350)]
    for M,r in cases:
        if 2*r>=M:continue
        vals=[ray_exact(M,r,(x,y,z)) for x,y,z,m in CLASSES]
        g=0
        for a in vals:g=gcd(g,a)
        print('CONTENT',M,r,'digits',len(str(g)),'g',g,'cyclic',cyclic_factor(M,r),'quotient',g//cyclic_factor(M,r))
        for w in (1,2,4,8):
            gg=state_content(M,r,w)
            print(' WINDOW',w,'digits',len(str(gg)),'g',gg)
    for p in MODS:
        for M in (754,1200,1500):
            seq=F_sequence(M,p)
            ans=guess(seq,p)
            if ans is None:
                print('REC',p,M,'NONE')
            else:
                o,d,v=ans
                roots=summarize_vector(o,d,v,p)
                print('REC',p,M,'order',o,'degree',d,'roots',roots,'vector',v)
    print('PASS')
if __name__=='__main__':main()
