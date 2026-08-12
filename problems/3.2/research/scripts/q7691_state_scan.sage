from sage.all import *
import time

# Q7691: exact finite-field target-state scan and bounded-degree hyperplane audit.
P_MAX = 5000
DS = [1,2,3,6]


def P(n):
    n = ZZ(n)
    return 34*n^3 + 51*n^2 + 27*n + 5

# sigma_3 sieve, exact over ZZ
sig3 = [ZZ(0)]*(P_MAX+1)
for d in range(1, P_MAX+1):
    dd = ZZ(d)^3
    for m in range(d, P_MAX+1, d):
        sig3[m] += dd


def apery_mod_targets(p):
    Fp = GF(p)
    b = [Fp(1), Fp(5)]
    for n in range(1, p-1):
        nn = Fp(n)
        nxt = (Fp(P(n))*b[n] - nn^3*b[n-1]) / Fp(n+1)^3
        b.append(nxt)
    assert len(b) == p
    # Endpoint audit.  The theorem-level reason is reflection
    # b_{p-1-r}=b_r mod p, hence b_{p-1}=b_0=1.
    assert b[0] == 1 and b[p-1] == 1
    return b, [r for r in range(1,p-1) if b[r] == 0]


def q_objects_mod(p, prec):
    Fp = GF(p)
    PS = PowerSeriesRing(Fp, 'q', default_prec=prec)
    q = PS.gen()
    A = PS(1)
    for a in range(1, prec):
        if a % 3 != 0:
            A *= (1 + q^a)^12
            A = A.add_bigoh(prec)
    E = PS(1)
    for m in range(1, prec):
        e = -5
        if m % 2 == 0: e += 7
        if m % 3 == 0: e += 7
        if m % 6 == 0: e -= 5
        if e > 0:
            E *= (1-q^m)^e
        elif e < 0:
            E /= (1-q^m)^(-e)
        E = E.add_bigoh(prec)
    H = (1 - q*A.derivative()/A).add_bigoh(prec)
    EH = (E*H).add_bigoh(prec)
    return PS, q, A, E, H, EH


def state_for_target(p, r, bprev, objs):
    # At a target b_r=0, W_r^(d)=r^3*b_{r-1}*kappa_r^(d).
    # K^(d)(t(q))=E(q)*(240/d^3) S(d tau), and coefficient extraction gives
    # kappa_r^(d)=[q^r] U^(d)*E*H*A^r.
    PS,q,A,E,H,EH = objs
    Fp = GF(p)
    base = (EH * (A^r)).add_bigoh(r+1)
    W = {}
    for d in DS:
        kap = Fp(0)
        scale = Fp(240) / Fp(d)^3
        for m in range(1, r//d + 1):
            coeff = scale * Fp(sig3[m]) / Fp(m)^3
            kap += coeff * base[r-d*m]
        W[d] = Fp(r)^3 * bprev * kap
    # Exact unit relation; p>=7 makes 240 a unit.
    assert W[1] - 28*W[2] + 63*W[3] - 36*W[6] == Fp(240)
    return tuple(ZZ(W[d]) for d in DS)


t0 = time.time()
samples = []   # (p,r,b_{r-1},W1,W2,W3,W6), all representatives in [0,p)
prime_target_counts = []
for p in prime_range(7, P_MAX+1):
    bmod, targets = apery_mod_targets(ZZ(p))
    if not targets:
        continue
    maxr = max(targets)
    objs = q_objects_mod(ZZ(p), maxr+1)
    prime_target_counts.append((ZZ(p),len(targets)))
    for r in targets:
        W = state_for_target(ZZ(p), ZZ(r), bmod[r-1], objs)
        samples.append((ZZ(p),ZZ(r),ZZ(bmod[r-1])) + W)

scan_seconds = time.time()-t0
print('P_MAX',P_MAX)
print('TARGET_PAIRS',len(samples),'TARGET_PRIMES',len(prime_target_counts))
print('MAX_TARGETS_ONE_PRIME',max([c for _,c in prime_target_counts] or [0]))
print('TARGET_COUNT_HIST', sorted({c:sum(1 for _,cc in prime_target_counts if cc==c) for c in set(cc for _,cc in prime_target_counts)}.items()))
print('SCAN_SECONDS', scan_seconds)


def vals(sample):
    p,r,bp,w1,w2,w3,w6 = sample
    Fp=GF(p)
    w1,w2,w3,w6 = map(Fp,[w1,w2,w3,w6])
    Aodd = w1 - 36*w6
    Bodd = 4*w2 - 9*w3
    assert Aodd - 7*Bodd == Fp(240)
    Wtr = (-3*Aodd + Bodd)/Fp(20)
    Ceven = w1 + 36*w6
    Deven = 4*w2 + 9*w3
    return {
        'W1':w1,'W2':w2,'W3':w3,'W6':w6,
        'Aodd':Aodd,'Bodd':Bodd,'Wtr':Wtr,
        'Ceven':Ceven,'Deven':Deven,
    }

names=['W1','W2','W3','W6','Aodd','Bodd','Wtr','Ceven','Deven']
cover={name:[] for name in names}
for s in samples:
    vv=vals(s)
    for name in names:
        if vv[name]==0:
            cover[name].append((s[0],s[1]))
uncovered=[(s[0],s[1]) for s in samples if not any(vals(s)[name]==0 for name in names)]
print('ZERO_COUNTS', {name:len(cover[name]) for name in names})
print('NATURAL_UNION_UNCOVERED_COUNT',len(uncovered),'FIRST20',uncovered[:20])

for lock in [(11,5),(17,13),(19,8),(2237,492)]:
    ss=[s for s in samples if (s[0],s[1])==lock]
    print('LOCK',lock,'FOUND',len(ss))
    for s in ss:
        vv=vals(s)
        print('LOCK_STATE',lock,'BPREV',s[2],'W',s[3:],'DERIVED',{k:ZZ(vv[k]) for k in vv})

ss17=[s for s in samples if (s[0],s[1])==(17,13)]
assert len(ss17)==1
v17=vals(ss17[0])
assert all(v17[name] != 0 for name in ['W1','W2','W3','W6'])
assert v17['Wtr'] == 0

ss2237=[s for s in samples if (s[0],s[1])==(2237,492)]
assert len(ss2237)==1
assert vals(ss2237[0])['Wtr'] == 0

# ------------------------------------------------------------------
# Degree <= 8 congruence-lattice search.
# We quotient the universal identity by eliminating W1.  The enlarged affine
# feature set is (1,b_{r-1},W2,W3,W6), with a polynomial coefficient of degree
# <= 8 on each coordinate.  Thus a candidate has 5*9=45 integer coefficients.
# Any fixed W-linear form with polynomial coefficients and optional b_{r-1}/unit
# coordinate occurs in this class after substituting
# W1=240+28 W2-63 W3+36 W6.
# ------------------------------------------------------------------
DEG = 8


def feature_vector(sample):
    p,r,bp,w1,w2,w3,w6 = map(ZZ,sample)
    Fp=GF(p)
    base=[Fp(1),Fp(bp),Fp(w2),Fp(w3),Fp(w6)]
    out=[]
    for x in base:
        rr=Fp(1)
        for k in range(DEG+1):
            out.append(x*rr)
            rr *= Fp(r)
    return vector(Fp,out)


def constraint_lattice(train):
    n=5*(DEG+1)
    B=identity_matrix(ZZ,n)
    used=0
    for s in train:
        p=ZZ(s[0]); Fp=GF(p); f=feature_vector(s)
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
            # centered lift keeps intermediate bases a little smaller
            z=ZZ((-h[i]*inv).lift())
            if z > p//2: z -= p
            K[row,pivot]=z
            row += 1
        K[row,pivot]=p
        B=K*B
        used += 1
        if used % 25 == 0:
            B=B.LLL(delta=0.75)
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
    p=ZZ(s[0]); Fp=GF(p); f=feature_vector(s)
    return sum(Fp(ZZ(c[j]))*f[j] for j in range(len(c)))


def candidate_record(c,train,hold):
    c=primitive(c)
    return {
        'train_hits':sum(1 for s in train if eval_vec(c,s)==0),
        'hold_hits':sum(1 for s in hold if eval_vec(c,s)==0),
        'maxbits':max([abs(x).nbits() for x in c if x] or [0]),
        'norm2':sum(x*x for x in c),
        'coeffs':list(c),
    }

train=[s for s in samples if s[0] <= 2500]
hold=[s for s in samples if s[0] > 2500]
print('HYP_FULL_TRAIN',len(train),'HOLD',len(hold),'NV',5*(DEG+1))
B,used=constraint_lattice(train)
records=[candidate_record(B.row(i),train,hold) for i in range(min(12,B.nrows()))]
records.sort(key=lambda z:(-z['hold_hits'],z['norm2']))
print('HYP_FULL_USED',used)
for j,rec in enumerate(records[:5]):
    print('HYP_FULL_CAND',j,'TRAIN_HITS',rec['train_hits'],'HOLD_HITS',rec['hold_hits'],
          'MAXBITS',rec['maxbits'],'COEFFS',rec['coeffs'])

# Natural Atkin-Lehner prime strata.  Fit one shortest degree-8 form in each
# p mod 24 class on the training half, then test both its own class and the
# union of the eight fixed forms on the holdout half.
strata_forms=[]
for cls in [1,5,7,11,13,17,19,23]:
    tr=[s for s in train if s[0] % 24 == cls]
    ho=[s for s in hold if s[0] % 24 == cls]
    if not tr:
        print('HYP_STRATUM',cls,'SKIP_NO_TRAIN','HOLD',len(ho))
        continue
    Bs,us=constraint_lattice(tr)
    c=primitive(Bs.row(0))
    rec=candidate_record(c,tr,ho)
    print('HYP_STRATUM',cls,'TRAIN',len(tr),'HOLD',len(ho),'USED',us,
          'TRAIN_HITS',rec['train_hits'],'HOLD_HITS',rec['hold_hits'],
          'MAXBITS',rec['maxbits'],'COEFFS',rec['coeffs'])
    strata_forms.append((cls,c))

if strata_forms:
    union_hits=[]
    union_miss=[]
    for s in hold:
        if any(eval_vec(c,s)==0 for _,c in strata_forms):
            union_hits.append((s[0],s[1]))
        else:
            union_miss.append((s[0],s[1]))
    print('HYP_STRATA_FIXED_UNION_HITS',len(union_hits),'MISSES',len(union_miss),
          'FIRST_MISS',union_miss[:10])

print('DONE_STATE_SCAN_AND_HYPERPLANE_AUDIT')
