from sage.all import *

PMAX = 3000


def P(n):
    n = ZZ(n)
    return 34*n^3 + 51*n^2 + 27*n + 5


def apery_mod(p):
    Fp = GF(p)
    b = [Fp(1), Fp(5)]
    for n in range(2, p):
        b.append((Fp(P(n-1))*b[n-1] - Fp((n-1)^3)*b[n-2]) / Fp(n^3))
    return b


def g_mod(p, b):
    Fp = GF(p)
    R = PowerSeriesRing(Fp, 't', default_prec=p)
    q = [Fp(1)]
    if p > 1:
        q.append(Fp(17))
    for n in range(1, p-1):
        q.append((Fp(17*(2*n+1))*q[n] - Fp(n)*q[n-1]) / Fp(n+1))
    F = R(b).add_bigoh(p)
    Q = R(q).add_bigoh(p)
    G = (Q/(F*F)).add_bigoh(p)
    return [Fp(G[n]) for n in range(p)]


def kappa_mod(p, g):
    Fp = GF(p)
    k = [Fp(0), Fp(-36)]
    for n in range(2, p):
        k.append((Fp(P(n-1))*k[n-1] - Fp((n-1)^3)*k[n-2] - Fp(5)*g[n]) / Fp(n^3))
    return k


def xi_mod(p, b, g):
    Fp = GF(p)
    xi = [Fp(-1)]
    for n in range(1, p):
        xi.append(xi[-1] - Fp(5)*g[n]*b[n-1])
    return xi


def force_common_sources(p, b, g, rows):
    Fp = GF(p)
    out = list(g)
    deltas = []
    for r in rows:
        k = kappa_mod(p, out)
        assert b[r] == 0
        delta = Fp(r^3)*k[r]/Fp(5)
        out[r] += delta
        deltas.append(int(delta))
        assert kappa_mod(p, out)[r] == 0
    k = kappa_mod(p, out)
    xi = xi_mod(p, b, out)
    assert all(k[r] == 0 and xi[r] == 0 for r in rows)
    return out, deltas


def N_values(p, r, h):
    Fp = GF(p)
    if h == 1:
        return [Fp(1)]
    N = [None, Fp(1), Fp(P(r+1))]
    for j in range(2, h):
        N.append(Fp(P(r+j))*N[j] - Fp((r+j)^6)*N[j-1])
    return N


def marked_M(p, r, h, g):
    Fp = GF(p)
    N = N_values(p, r, h)
    total = Fp(0)
    for j in range(1, h):
        tail = Fp(1)
        for q in range(j+1, h):
            tail *= Fp((r+q)^3)
        total += g[r+j+1] * N[j] * tail
    return total, N

# A. Actual-source scan.  Finite evidence only; also search for collisions of
# Xi values on the Hasse-zero set, a stronger property than zero-fiber uniqueness.
collision_records = []
zero_common_records = []
repeated_common_primes = []
max_z = 0
for pz in prime_range(7, PMAX+1):
    p = int(pz)
    b = apery_mod(p)
    targets = [r for r in range(1,p) if b[r] == 0]
    max_z = max(max_z, len(targets))
    if not targets:
        continue
    g = g_mod(p,b)
    xi = xi_mod(p,b,g)
    seen = {}
    collisions = []
    commons = []
    for r in targets:
        key = int(xi[r])
        if key in seen:
            collisions.append((seen[key],r,key))
        else:
            seen[key] = r
        if xi[r] == 0:
            commons.append(r)
            zero_common_records.append((p,r))
    if collisions:
        collision_records.append((p,collisions[:5]))
    if len(commons) >= 2:
        repeated_common_primes.append((p,commons))

print('Q7684_SCAN_PMAX', PMAX)
print('MAX_HASSE_ZERO_COUNT', max_z)
print('XI_COLLISION_PRIMES', len(collision_records))
print('FIRST_XI_COLLISIONS', collision_records[:20])
print('COMMON_RECORDS', zero_common_records[:30])
print('REPEATED_COMMON_PRIMES', repeated_common_primes)

# B. Direct dyadic recurrence-compatible countermodel: p=19, R=5,
# Hasse zeros 8,10 lie in (R,2R], and p>2R.  Keep the canonical initial
# kappa line but change only the inhomogeneous source at the target rows.
p = 19
Fp = GF(p)
b19 = apery_mod(p)
g19 = g_mod(p,b19)
assert b19[8] == 0 and b19[10] == 0
forced19, deltas19 = force_common_sources(p,b19,g19,[8,10])
print('SOURCE_FLEX_DYADIC_19', 'R',5,'rows',(8,10),'deltas',tuple(deltas19))

# C. p=17 gives a clean marked-return example at rows 3 and 13.
p = 17
Fp = GF(p)
b = apery_mod(p)
g = g_mod(p,b)
k = kappa_mod(p,g)
xi = xi_mod(p,b,g)
assert b[3] == 0 and b[13] == 0
assert k[13] == 0 and k[3] != 0
assert xi[13] == 0 and xi[3] != 0

g2, deltas17 = force_common_sources(p,b,g,[3,13])
print('SOURCE_FLEX_17', 'deltas',tuple(deltas17))

r = 3
h = 10
M0,N0 = marked_M(p,r,h,g)
M2,N2 = marked_M(p,r,h,g2)
assert N0[h] == 0
assert N0[h-1] != 0
assert M0 != 0
assert M2 == 0
print('RETURN_17_3_13', 'Nh', int(N0[h]), 'Nhminus1', int(N0[h-1]),
      'Mcanonical', int(M0), 'Mperturbed', int(M2))

print('Q7684_AUDIT=PASS')
