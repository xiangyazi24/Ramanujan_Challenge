from sage.all import *
import glob, os, collections, time

DEG=8
files=sorted(glob.glob(os.environ.get('SAMPLE_GLOB','artifacts/**/*.tsv'),recursive=True))
if not files:
    files=sorted(glob.glob('*.tsv'))
print('ANALYZE_FILES',files)
rows=[]
seen=set()
for fn in files:
    with open(fn) as fh:
        header=fh.readline().strip().split('\t')
        assert header==['p','r','bprev','W1','W2','W3','W6'],(fn,header)
        for line in fh:
            if not line.strip(): continue
            row=tuple(map(ZZ,line.strip().split('\t')))
            key=(row[0],row[1])
            assert key not in seen,key
            seen.add(key); rows.append(row)
rows.sort()
samples=rows
print('MERGED_TARGET_PAIRS',len(samples),'TARGET_PRIMES',len(set(s[0] for s in samples)))

cnt=collections.Counter(s[0] for s in samples)
hist=collections.Counter(cnt.values())
print('MERGED_MAX_TARGETS',max(cnt.values() or [0]),'HIST',sorted(hist.items()))


def vals(s):
    p,r,bp,w1,w2,w3,w6=s; Fp=GF(p)
    w1,w2,w3,w6=map(Fp,[w1,w2,w3,w6])
    Aodd=w1-36*w6
    Bodd=4*w2-9*w3
    assert Aodd-7*Bodd==Fp(240)
    Wtr=(-3*Aodd+Bodd)/Fp(20)
    Ceven=w1+36*w6
    Deven=4*w2+9*w3
    return {'W1':w1,'W2':w2,'W3':w3,'W6':w6,
            'Aodd':Aodd,'Bodd':Bodd,'Wtr':Wtr,
            'Ceven':Ceven,'Deven':Deven}

names=['W1','W2','W3','W6','Aodd','Bodd','Wtr','Ceven','Deven']
zero_counts={name:sum(1 for s in samples if vals(s)[name]==0) for name in names}
uncovered=[(s[0],s[1]) for s in samples if not any(vals(s)[name]==0 for name in names)]
print('MERGED_ZERO_COUNTS',zero_counts)
print('MERGED_NATURAL_UNCOVERED',len(uncovered),'FIRST20',uncovered[:20])

for lock in [(11,5),(17,13),(19,8),(2237,492)]:
    ss=[s for s in samples if (s[0],s[1])==lock]
    print('MERGED_LOCK',lock,'FOUND',len(ss))
    for s in ss:
        print('MERGED_LOCK_STATE',lock,'BPREV',s[2],'W',s[3:],
              'DERIVED',{k:ZZ(v) for k,v in vals(s).items()})

assert len([s for s in samples if (s[0],s[1])==(17,13)])==1
assert vals([s for s in samples if (s[0],s[1])==(17,13)][0])['Wtr']==0
assert len([s for s in samples if (s[0],s[1])==(2237,492)])==1
assert vals([s for s in samples if (s[0],s[1])==(2237,492)][0])['Wtr']==0

# Exact affine degree<=8 feature space after quotienting the unit identity.
# A form is sum_{X in {1,bprev,W2,W3,W6}} P_X(r) X, deg P_X<=8.
def feature_vector(s):
    p,r,bp,w1,w2,w3,w6=map(ZZ,s); Fp=GF(p)
    base=[Fp(1),Fp(bp),Fp(w2),Fp(w3),Fp(w6)]
    out=[]
    for x in base:
        rr=Fp(1)
        for k in range(DEG+1):
            out.append(x*rr); rr*=Fp(r)
    return vector(Fp,out)


def constraint_lattice(train):
    n=5*(DEG+1)
    B=identity_matrix(ZZ,n); used=0
    for s in train:
        p=ZZ(s[0]); Fp=GF(p); f=feature_vector(s)
        h=[Fp(sum(B[i,j]*ZZ(f[j]) for j in range(n))) for i in range(n)]
        pivot=next((i for i,x in enumerate(h) if x),None)
        if pivot is None: continue
        inv=1/h[pivot]
        K=zero_matrix(ZZ,n,n); row=0
        for i in range(n):
            if i==pivot: continue
            K[row,i]=1
            z=ZZ((-h[i]*inv).lift())
            if z>p//2: z-=p
            K[row,pivot]=z; row+=1
        K[row,pivot]=p
        B=K*B; used+=1
        if used%25==0: B=B.LLL(delta=0.75)
    return B.LLL(delta=0.99),used


def primitive(c):
    cc=[ZZ(x) for x in c]
    g=gcd([abs(x) for x in cc if x] or [ZZ(1)])
    cc=[x//g for x in cc]
    for x in cc:
        if x:
            if x<0: cc=[-y for y in cc]
            break
    return vector(ZZ,cc)


def eval_vec(c,s):
    Fp=GF(ZZ(s[0])); f=feature_vector(s)
    return sum(Fp(ZZ(c[j]))*f[j] for j in range(len(c)))


def rec(c,tr,ho):
    c=primitive(c)
    return (sum(1 for s in tr if eval_vec(c,s)==0),
            sum(1 for s in ho if eval_vec(c,s)==0),
            max([abs(x).nbits() for x in c if x] or [0]),
            sum(x*x for x in c),list(c))

train=[s for s in samples if s[0]<=2500]
hold=[s for s in samples if s[0]>2500]
print('MERGED_HYP_TRAIN',len(train),'HOLD',len(hold),'NV',45)
t0=time.time(); B,used=constraint_lattice(train)
records=[rec(B.row(i),train,hold) for i in range(min(12,B.nrows()))]
records.sort(key=lambda z:(-z[1],z[3]))
print('MERGED_HYP_SECONDS',time.time()-t0,'USED',used)
for j,z in enumerate(records[:5]):
    print('MERGED_HYP_CAND',j,'TRAIN',z[0],'HOLD',z[1],'MAXBITS',z[2],'COEFFS',z[4])

# Natural eight-class p mod 24 family, one shortest form learned per class on train.
forms=[]
for cls in [1,5,7,11,13,17,19,23]:
    tr=[s for s in train if s[0]%24==cls]
    ho=[s for s in hold if s[0]%24==cls]
    if not tr:
        print('MERGED_STRATUM',cls,'NO_TRAIN','HOLD',len(ho)); continue
    Bs,us=constraint_lattice(tr)
    c=primitive(Bs.row(0)); z=rec(c,tr,ho)
    print('MERGED_STRATUM',cls,'TRAIN_N',len(tr),'HOLD_N',len(ho),'USED',us,
          'TRAIN_HITS',z[0],'HOLD_HITS',z[1],'MAXBITS',z[2],'COEFFS',z[4])
    forms.append((cls,c))

miss=[]; hits=[]
for s in hold:
    if any(eval_vec(c,s)==0 for _,c in forms): hits.append((s[0],s[1]))
    else: miss.append((s[0],s[1]))
print('MERGED_STRATA_FIXED_UNION_HITS',len(hits),'MISSES',len(miss),'FIRST20_MISS',miss[:20])
print('DONE_MERGED_ANALYSIS')
