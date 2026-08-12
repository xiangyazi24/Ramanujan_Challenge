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
    # Q=1/sqrt(1-34t+t^2)=sum P_n(17)t^n.
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

# A. Actual-source scan: stronger Green-injectivity on the Hasse-zero set.
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

# B. Exact source-flexibility countermodel at p=17, target rows 3 and 13.
p = 17
Fp = GF(p)
b = apery_mod(p)
g = g_mod(p,b)
k = kappa_mod(p,g)
xi = xi_mod(p,b,g)
assert b[3] == 0 and b[13] == 0
assert k[13] == 0 and k[3] != 0
assert xi[13] == 0 and xi[3] != 0

# Perturb g_3.  Since kappa_3 depends on g_3 with coefficient -5/3^3,
# this unique delta makes row 3 common while preserving the initial line.
g1 = list(g)
d3 = Fp(3^3) * k[3] / Fp(5)
g1[3] += d3
k1 = kappa_mod(p,g1)
assert k1[3] == 0
# Later source changes do not affect row 3.  Now use g_13 to force row 13.
g2 = list(g1)
d13 = Fp(13^3) * k1[13] / Fp(5)
g2[13] += d13
k2 = kappa_mod(p,g2)
xi2 = xi_mod(p,b,g2)
assert k2[3] == 0 and k2[13] == 0
assert xi2[3] == 0 and xi2[13] == 0
print('SOURCE_FLEX_17', 'delta3', int(d3), 'delta13', int(d13),
      'k3', int(k2[3]), 'k13', int(k2[13]))

# C. Verify the known continuant/marked-Green reduction on the same pair.
r = 3
h = 10
M0,N0 = marked_M(p,r,h,g)
M2,N2 = marked_M(p,r,h,g2)
assert N0[h] == 0
assert N0[h-1] != 0
assert M0 != 0
assert M2 == 0
# Terminal source coefficient is exactly N_{h-1}, hence a p-unit.
print('RETURN_17_3_13', 'Nh', int(N0[h]), 'Nhminus1', int(N0[h-1]),
      'Mcanonical', int(M0), 'Mperturbed', int(M2))

print('Q7684_AUDIT=PASS')
