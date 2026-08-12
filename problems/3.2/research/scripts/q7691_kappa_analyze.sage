from sage.all import *
import glob,time

DEG=8
files=sorted(glob.glob('artifacts/**/*.tsv',recursive=True))
print('KAPPA_FILES',files)
samples=[]; seen=set()
for fn in files:
    with open(fn) as fh:
        assert fh.readline().strip().split('\t')==['p','r','bprev','W1','W2','W3','W6']
        for line in fh:
            if not line.strip(): continue
            s=tuple(map(ZZ,line.strip().split('\t')))
            assert (s[0],s[1]) not in seen
            seen.add((s[0],s[1])); samples.append(s)
samples.sort()

# This is the literal four-kappa affine class:
# L=P0(r)+sum_d Pd(r) kappa^(d), deg each <=8.
# At a target W_d=u*kappa_d with u=r^3*bprev a p-unit, so evaluate kappa=W/u.
def feature_kappa(s):
    p,r,bp,w1,w2,w3,w6=map(ZZ,s); Fp=GF(p)
    u=Fp(r)^3*Fp(bp); assert u != 0
    base=[Fp(1),Fp(w1)/u,Fp(w2)/u,Fp(w3)/u,Fp(w6)/u]
    out=[]
    for x in base:
        rr=Fp(1)
        for k in range(DEG+1): out.append(x*rr); rr*=Fp(r)
    return vector(Fp,out)

# Also the equivalent denominator-cleared class: u P0(r)+sum Pd(r) Wd.
def constraint_lattice(train):
    n=5*(DEG+1); B=identity_matrix(ZZ,n); used=0
    for s in train:
        p=ZZ(s[0]); Fp=GF(p); f=feature_kappa(s)
        h=[Fp(sum(B[i,j]*ZZ(f[j]) for j in range(n))) for i in range(n)]
        piv=next((i for i,x in enumerate(h) if x),None)
        if piv is None: continue
        inv=1/h[piv]; K=zero_matrix(ZZ,n,n); row=0
        for i in range(n):
            if i==piv: continue
            K[row,i]=1; z=ZZ((-h[i]*inv).lift())
            if z>p//2: z-=p
            K[row,piv]=z; row+=1
        K[row,piv]=p; B=K*B; used+=1
        if used%25==0: B=B.LLL(delta=0.75)
    return B.LLL(delta=0.99),used

def primitive(c):
    cc=[ZZ(x) for x in c]; g=gcd([abs(x) for x in cc if x] or [ZZ(1)]); cc=[x//g for x in cc]
    for x in cc:
        if x:
            if x<0: cc=[-y for y in cc]
            break
    return vector(ZZ,cc)

def ev(c,s):
    Fp=GF(ZZ(s[0])); f=feature_kappa(s)
    return sum(Fp(ZZ(c[j]))*f[j] for j in range(len(c)))

def rec(c,tr,ho):
    c=primitive(c)
    return (sum(ev(c,s)==0 for s in tr),sum(ev(c,s)==0 for s in ho),
            max([abs(x).nbits() for x in c if x] or [0]),sum(x*x for x in c),list(c))

tr=[s for s in samples if s[0]<=2500]; ho=[s for s in samples if s[0]>2500]
print('KAPPA_TRAIN',len(tr),'HOLD',len(ho),'NV',45)
t0=time.time(); B,used=constraint_lattice(tr)
rs=[rec(B.row(i),tr,ho) for i in range(min(12,B.nrows()))]; rs.sort(key=lambda z:(-z[1],z[3]))
print('KAPPA_SECONDS',time.time()-t0,'USED',used)
for i,z in enumerate(rs[:5]): print('KAPPA_CAND',i,'TRAIN',z[0],'HOLD',z[1],'MAXBITS',z[2],'COEFFS',z[4])

forms=[]
for cls in [1,5,7,11,13,17,19,23]:
    tt=[s for s in tr if s[0]%24==cls]; hh=[s for s in ho if s[0]%24==cls]
    if not tt: continue
    C,u=constraint_lattice(tt); c=primitive(C.row(0)); z=rec(c,tt,hh)
    print('KAPPA_STRATUM',cls,'TRAIN_N',len(tt),'HOLD_N',len(hh),'TRAIN_HITS',z[0],'HOLD_HITS',z[1],'MAXBITS',z[2])
    forms.append(c)
miss=[(s[0],s[1]) for s in ho if not any(ev(c,s)==0 for c in forms)]
hit=[(s[0],s[1]) for s in ho if any(ev(c,s)==0 for c in forms)]
print('KAPPA_STRATA_UNION_HITS',len(hit),'MISSES',len(miss),'FIRST20_MISS',miss[:20])
print('DONE_KAPPA_ANALYSIS')
