from sage.all import *
from math import gcd as pygcd
from collections import defaultdict

PMAX = 20000


def P(n):
    n = ZZ(n)
    return 34*n^3 + 51*n^2 + 27*n + 5


def canon(v):
    v = tuple(int(x) for x in v)
    if not any(v):
        return None
    g = 0
    for x in v:
        g = pygcd(g, abs(x))
    v = tuple(x//g for x in v)
    for x in v:
        if x:
            if x < 0:
                v = tuple(-y for y in v)
            break
    return v


def dotmod(c, v, p):
    return sum(int(c[i])*int(v[i]) for i in range(4)) % int(p)


def companion_vectors_for_prime(p, b, targets, avec):
    p = ZZ(p)
    m = max(targets)
    Fp = GF(p)

    # Q(t)=1/sqrt(1-34t+t^2)=sum P_n(17)t^n.
    qcoef = [Fp(1)]
    if m >= 1:
        qcoef.append(Fp(17))
    for n in range(1, m):
        qn1 = (Fp(17*(2*n+1))*qcoef[n] - Fp(n)*qcoef[n-1]) / Fp(n+1)
        qcoef.append(qn1)

    R = PowerSeriesRing(Fp, 't', default_prec=m+1)
    Fser = R([Fp(x) for x in b[:m+1]]).add_bigoh(m+1)
    Qser = R(qcoef).add_bigoh(m+1)
    Cser = (Qser / (Fser*Fser)).add_bigoh(m+1)
    g = [Fp(Cser[n]) for n in range(m+1)]

    # Fricke basis of Eisenstein sources:
    # O1=E4-36E4(6), O2=4E4(2)-9E4(3),
    # P1=E4+36E4(6), P2=4E4(2)+9E4(3).
    # Exact source identities after multiplying by C=1/(EH):
    # L KO1=35(C-1)-5t, L KO2=5(C-1)-35t,
    # L KP1=(37-130t)/sqrt(Delta)-37C,
    # L KP2=(13-130t)/sqrt(Delta)-13C.
    srcO1 = [Fp(0)]*(m+1)
    srcO2 = [Fp(0)]*(m+1)
    srcP1 = [Fp(0)]*(m+1)
    srcP2 = [Fp(0)]*(m+1)
    for n in range(1, m+1):
        qprev = qcoef[n-1]
        srcO1[n] = Fp(35)*g[n] - (Fp(5) if n == 1 else Fp(0))
        srcO2[n] = Fp(5)*g[n] - (Fp(35) if n == 1 else Fp(0))
        srcP1[n] = Fp(37)*qcoef[n] - Fp(130)*qprev - Fp(37)*g[n]
        srcP2[n] = Fp(13)*qcoef[n] - Fp(130)*qprev - Fp(13)*g[n]

    def solve(src):
        z = [Fp(0)]*(m+1)
        for n in range(1, m+1):
            rhs = Fp(P(n-1))*z[n-1] + src[n]
            if n >= 2:
                rhs -= Fp((n-1)^3)*z[n-2]
            z[n] = rhs / Fp(n^3)
        return z

    o1 = solve(srcO1)
    o2 = solve(srcO2)
    p1 = solve(srcP1)
    p2 = solve(srcP2)

    out = {}
    for r in targets:
        k1 = (o1[r] + p1[r]) / Fp(2)
        k6 = (p1[r] - o1[r]) / Fp(72)
        k2 = (o2[r] + p2[r]) / Fp(8)
        k3 = (p2[r] - o2[r]) / Fp(18)
        v = tuple(int(x) for x in (k1,k2,k3,k6))
        # Exact unit-separator identity K1-28K2+63K3-36K6=40 a.
        sep = dotmod((1,-28,63,-36), v, p)
        assert sep == (40*int(avec[r])) % int(p), (p,r,v,sep,avec[r])
        if p >= 7:
            assert sep != 0, ("separator vanished",p,r,v)
        out[r] = v
    return out


records = []
prime_target_counts = []
for p in prime_range(5, PMAX+1):
    Fp = GF(p)
    b = [Fp(1), Fp(5)]
    a = [Fp(0), Fp(6)]
    for n in range(2, int(p)):
        b.append((Fp(P(n-1))*b[n-1] - Fp((n-1)^3)*b[n-2]) / Fp(n^3))
        a.append((Fp(P(n-1))*a[n-1] - Fp((n-1)^3)*a[n-2]) / Fp(n^3))
    targets = [r for r in range(1,int(p)) if b[r] == 0]
    if not targets:
        continue
    vecs = companion_vectors_for_prime(p,b,targets,a)
    prime_target_counts.append((int(p),len(targets)))
    for r in targets:
        records.append((int(p),int(r),vecs[r]))

print('Q7676_SCAN_PMAX', PMAX)
print('TARGET_PRIMES', len(prime_target_counts))
print('TARGET_COUNT', len(records))
print('MAX_Z', max(z for p,z in prime_target_counts))
print('FIRST_RECORDS', records[:12])

# Required normalization regression.
r1713 = [z for z in records if z[0] == 17 and z[1] == 13]
assert len(r1713) == 1
print('VECTOR_17_13', r1713[0][2])
assert r1713[0][2] == (9,8,9,7), r1713[0]

TR = (-3,4,-9,108)
SEP = (1,-28,63,-36)
O1 = (1,0,0,-36)
O2 = (0,4,-9,0)
P1 = (1,0,0,36)
P2 = (0,4,9,0)

for name,c in [('transverse',TR),('separator',SEP),('O1',O1),('O2',O2),('P1',P1),('P2',P2)]:
    hits = [(p,r) for p,r,v in records if dotmod(c,v,p) == 0]
    print('FORM_HITS', name, len(hits), hits[:30])

# Character split diagnostics; chi=(-6/p), p>=5.
chars = defaultdict(list)
for i,(p,r,v) in enumerate(records):
    chars[int(kronecker(-6,p))].append(i)
print('CHAR_COUNTS', {k:len(v) for k,v in chars.items()})

# Natural Fricke odd/even two-parameter families. Enumerate primitive (a,b)
# with |a|,|b| <= Hpair and accumulate exact target masks.
def two_plane_cover(odd=True,Hpair=300):
    masks = {}
    union = 0
    for tid,(p,r,v) in enumerate(records):
        if odd:
            A = (v[0] - 36*v[3]) % p
            B = (4*v[1] - 9*v[2]) % p
        else:
            A = (v[0] + 36*v[3]) % p
            B = (4*v[1] + 9*v[2]) % p
        local = set()
        for aa in range(-Hpair,Hpair+1):
            if aa == 0 and B == 0:
                # all b work; handled by loop below
                pass
            for bb in range(-Hpair,Hpair+1):
                if aa == 0 and bb == 0:
                    continue
                if (aa*A + bb*B) % p:
                    continue
                c2 = canon((aa,bb))
                if c2 is None or max(abs(c2[0]),abs(c2[1])) > Hpair:
                    continue
                local.add(c2)
        if local:
            union |= (1<<tid)
        for c2 in local:
            masks[c2] = masks.get(c2,0) | (1<<tid)
    top = sorted(((mask.bit_count(),c,mask) for c,mask in masks.items()), reverse=True)[:12]
    print('PLANE', 'odd' if odd else 'even', 'H',Hpair,'CANDS',len(masks),'COVERED',union.bit_count(),'MAX',top[0][0] if top else 0)
    print('PLANE_TOP', 'odd' if odd else 'even', [(n,c) for n,c,m in top])
    return masks,union

odd_masks, odd_union = two_plane_cover(True,300)
even_masks, even_union = two_plane_cover(False,300)

# Full 4D bounded-height family by meet-in-the-middle.
def full_cover(H):
    masks = {}
    union = 0
    rng = range(-H,H+1)
    for tid,(p,r,v) in enumerate(records):
        left = defaultdict(list)
        for a0 in rng:
            for a1 in rng:
                left[(a0*v[0]+a1*v[1]) % p].append((a0,a1))
        local = set()
        for a2 in rng:
            for a3 in rng:
                need = (-(a2*v[2]+a3*v[3])) % p
                for a0,a1 in left.get(need,()):
                    if a0==a1==a2==a3==0:
                        continue
                    c = canon((a0,a1,a2,a3))
                    if c is None or max(map(abs,c)) > H:
                        continue
                    if dotmod(c,v,p) == 0:
                        local.add(c)
        if local:
            union |= (1<<tid)
        for c in local:
            masks[c] = masks.get(c,0) | (1<<tid)
    top = sorted(((mask.bit_count(),c,mask) for c,mask in masks.items()), reverse=True)[:15]
    print('FULL_H',H,'CANDS',len(masks),'COVERED',union.bit_count(),'OF',len(records),'MULTI',sum(1 for m in masks.values() if m.bit_count()>=2),'MAX',top[0][0] if top else 0)
    print('FULL_TOP',H,[(n,c) for n,c,m in top])

    # Greedy set cover on the bounded candidate pool.
    uncovered = (1<<len(records))-1
    chosen=[]
    vals=[(c,m) for c,m in masks.items()]
    while uncovered:
        best_c=None; best_m=0; best_n=0
        for c,m in vals:
            n=(m & uncovered).bit_count()
            if n>best_n:
                best_n=n; best_c=c; best_m=m
        if best_n==0:
            break
        chosen.append((best_c,best_n))
        uncovered &= ~best_m
    print('GREEDY',H,'COUNT',len(chosen),'UNCOVERED',uncovered.bit_count(),'HEAD',chosen[:20])
    if top:
        lb = (len(records)+top[0][0]-1)//top[0][0]
        print('COVER_LOWER_BY_MAX',H,lb)
    return masks,union,chosen

for H in [4,8,12]:
    full_cover(H)

print('DONE')
