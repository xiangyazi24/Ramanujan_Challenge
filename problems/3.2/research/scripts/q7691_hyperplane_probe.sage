from sage.all import *

# Load the P_MAX=500 prototype state scan from this branch.
load('problems/3.2/research/scripts/q7691_state_scan.sage')

# Recover b_{r-1} mod p for the enlarged affine feature set.
bprev_cache = {}
for p in sorted(set(s[0] for s in samples)):
    bb, tt = apery_mod_targets(ZZ(p))
    for r in tt:
        bprev_cache[(ZZ(p),ZZ(r))] = ZZ(bb[r-1])


def base_features(sample, mode):
    p,r,w1,w2,w3,w6 = map(ZZ,sample)
    Fp = GF(p)
    bp = Fp(bprev_cache[(p,r)])
    if mode == 'W':
        # eliminate W1 by the exact unit identity; constant coordinate represents 240 too.
        return [Fp(1), bp, Fp(w2), Fp(w3), Fp(w6)]
    if mode == 'K':
        u = Fp(r)^3 * bp
        assert u != 0
        return [Fp(1), bp, Fp(w1)/u, Fp(w2)/u, Fp(w3)/u, Fp(w6)/u]
    raise ValueError(mode)


def feature_vector(sample, mode, deg):
    p,r = ZZ(sample[0]),ZZ(sample[1])
    Fp=GF(p)
    base=base_features(sample,mode)
    out=[]
    for x in base:
        rr=Fp(1)
        for k in range(deg+1):
            out.append(x*rr)
            rr *= Fp(r)
    return vector(Fp,out)


def constraint_lattice(train, mode, deg):
    n = len(base_features(train[0],mode))*(deg+1)
    B = identity_matrix(ZZ,n)
    used=0
    for idx,s in enumerate(train):
        p=ZZ(s[0]); Fp=GF(p)
        f=feature_vector(s,mode,deg)
        # residues of the current row basis against this modular constraint
        h=[Fp(sum(B[i,j]*ZZ(f[j]) for j in range(n))) for i in range(n)]
        pivot=next((i for i,x in enumerate(h) if x != 0),None)
        if pivot is None:
            continue
        inv=1/h[pivot]
        K=zero_matrix(ZZ,n,n)
        row=0
        for i in range(n):
            if i==pivot: continue
            K[row,i]=1
            K[row,pivot]=ZZ((-h[i]*inv).lift())
            row += 1
        K[row,pivot]=p
        B=K*B
        used += 1
        if used % 5 == 0:
            B=B.LLL(delta=0.75)
    B=B.LLL(delta=0.99)
    return B,used


def eval_vec(c,s,mode,deg):
    p=ZZ(s[0]); Fp=GF(p)
    f=feature_vector(s,mode,deg)
    return sum(Fp(ZZ(c[j]))*f[j] for j in range(len(c)))


def primitive(c):
    cc=[ZZ(x) for x in c]
    g=gcd([abs(x) for x in cc if x] or [ZZ(1)])
    cc=[x//g for x in cc]
    for x in cc:
        if x:
            if x<0: cc=[-y for y in cc]
            break
    return vector(ZZ,cc)

train=[s for s in samples if s[0] <= 251]
hold=[s for s in samples if s[0] > 251]
print('HYP_TRAIN',len(train),'HOLD',len(hold))

for mode in ['W','K']:
    for deg in [0,1,2,4,8]:
        if not train: continue
        B,used=constraint_lattice(train,mode,deg)
        cand=[]
        for i in range(min(12,B.nrows())):
            c=primitive(B.row(i))
            train_hits=sum(1 for s in train if eval_vec(c,s,mode,deg)==0)
            hold_hits=sum(1 for s in hold if eval_vec(c,s,mode,deg)==0)
            maxbits=max([abs(x).nbits() for x in c if x] or [0])
            cand.append((train_hits,hold_hits,maxbits,list(c)))
        cand.sort(key=lambda z:(-z[1],z[2]))
        print('HYP',mode,'DEG',deg,'NV',B.ncols(),'USED',used,'BEST',cand[:3])

# Atkin-Lehner residue strata: one form per p mod 24, tested only when enough data.
for mode in ['W','K']:
    deg=2
    print('STRATA_MODE',mode,'DEG',deg)
    for cls in [1,5,7,11,13,17,19,23]:
        tr=[s for s in train if s[0] % 24 == cls]
        ho=[s for s in hold if s[0] % 24 == cls]
        if len(tr)<3 or not ho:
            print('STRATUM',cls,'SKIP',len(tr),len(ho)); continue
        B,used=constraint_lattice(tr,mode,deg)
        best=None
        for i in range(min(12,B.nrows())):
            c=primitive(B.row(i))
            hh=sum(1 for s in ho if eval_vec(c,s,mode,deg)==0)
            mb=max([abs(x).nbits() for x in c if x] or [0])
            rec=(hh,mb,list(c))
            if best is None or (-hh,mb)<(-best[0],best[1]): best=rec
        print('STRATUM',cls,'TRAIN',len(tr),'HOLD',len(ho),'BEST',best)

print('DONE_HYP_PROBE')
