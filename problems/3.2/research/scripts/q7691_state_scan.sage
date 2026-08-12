from sage.all import *
import time

# Exact finite-field scan for Q7691.
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
    # Endpoint audit.  This is the finite scan check; the theorem-level reason
    # is Apéry reflection b_{p-1-r}=b_r mod p.
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
    assert W[1] - 28*W[2] + 63*W[3] - 36*W[6] == Fp(240)
    return tuple(ZZ(W[d]) for d in DS)


t0 = time.time()
samples = []
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
        samples.append((ZZ(p),ZZ(r)) + W)

scan_seconds = time.time()-t0
print('P_MAX',P_MAX)
print('TARGET_PAIRS',len(samples),'TARGET_PRIMES',len(prime_target_counts))
print('MAX_TARGETS_ONE_PRIME',max([c for _,c in prime_target_counts] or [0]))
print('TARGET_COUNT_HIST', sorted({c:sum(1 for _,cc in prime_target_counts if cc==c) for c in set(cc for _,cc in prime_target_counts)}.items()))
print('SCAN_SECONDS', scan_seconds)


def vals(sample):
    p,r,w1,w2,w3,w6 = sample
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

for lock in [(11,5),(17,13),(2237,492)]:
    ss=[s for s in samples if (s[0],s[1])==lock]
    print('LOCK',lock,'FOUND',len(ss))
    for s in ss:
        vv=vals(s)
        print('LOCK_STATE',lock,'W',s[2:],'DERIVED',{k:ZZ(vv[k]) for k in vv})

ss17=[s for s in samples if (s[0],s[1])==(17,13)]
assert len(ss17)==1
v17=vals(ss17[0])
assert all(v17[name] != 0 for name in ['W1','W2','W3','W6'])
assert v17['Wtr'] == 0

ss2237=[s for s in samples if (s[0],s[1])==(2237,492)]
assert len(ss2237)==1
assert vals(ss2237[0])['Wtr'] == 0

print('DONE_STATE_SCAN')
