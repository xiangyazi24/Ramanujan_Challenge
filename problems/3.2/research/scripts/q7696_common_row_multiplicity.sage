from sage.all import *
import os
from collections import Counter

PMAX = int(os.environ.get('Q7696_PMAX', '20000'))


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
    """Canonical level-six Green/Eichler source g=1/(F^2*sqrt(1-34t+t^2))."""
    Fp = GF(p)
    R = PowerSeriesRing(Fp, 't', default_prec=p)
    # q_n=[t^n](1-34t+t^2)^(-1/2), i.e. Legendre P_n(17).
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
    """Canonical inhomogeneous companion: n^3 k_n-P(n-1)k_{n-1}+(n-1)^3k_{n-2}=-5g_n."""
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


def check_green_casoratian(p, b, g, k, xi):
    Fp = GF(p)
    # Xi_n = n^3(b_{n-1} k_n-b_n k_{n-1}); this is the exact
    # Duhamel/Green bridge, not a second independent endpoint equation.
    for n in range(1, p):
        rhs = Fp(n^3) * (b[n-1]*k[n] - b[n]*k[n-1])
        assert xi[n] == rhs, (p, n, xi[n], rhs)


def check_hasse_reflection(p, b):
    # Apéry reciprocity A_{p-1-r}=A_r (mod p).
    for r in range(p):
        assert b[r] == b[p-1-r], (p, r, b[r], b[p-1-r])


def check_two_endpoint_shooting(p, b, g, k, xi, r, s):
    Fp = GF(p)
    assert 0 < r < s < p
    assert b[r] == 0 and b[s] == 0
    # Exact telescoped Green law between the two Hasse-zero endpoints.
    period = sum((g[m]*b[m-1] for m in range(r+1, s+1)), Fp(0))
    lhs = Fp(s^3)*b[s-1]*k[s] - Fp(r^3)*b[r-1]*k[r]
    assert lhs == -Fp(5)*period
    assert xi[s] - xi[r] == -Fp(5)*period
    # If r is common, the single shooting condition for s is period=0.
    if xi[r] == 0:
        assert k[r] == 0
        assert (xi[s] == 0) == (period == 0)
        assert (k[s] == 0) == (period == 0)


common_records = []
common_by_prime = {}
xi_collision_records = []
hasse_zero_count_max = 0
hzero_count_argmax = []
mult_hist = Counter()
primes_with_hasse_zeros = 0

for pz in prime_range(7, PMAX+1):
    p = int(pz)
    b = apery_mod(p)
    check_hasse_reflection(p, b)
    hzeros = [r for r in range(1, p) if b[r] == 0]
    hz = len(hzeros)
    if hz > hasse_zero_count_max:
        hasse_zero_count_max = hz
        hzero_count_argmax = [p]
    elif hz == hasse_zero_count_max and hz > 0:
        hzero_count_argmax.append(p)
    if not hzeros:
        mult_hist[0] += 1
        continue

    primes_with_hasse_zeros += 1
    g = g_mod(p, b)
    xi = xi_mod(p, b, g)
    k = kappa_mod(p, g)
    check_green_casoratian(p, b, g, k, xi)

    # Hasse-zero commonness is exactly kappa=0 because r^3*b_{r-1} is a unit.
    commons = []
    seen = {}
    collisions = []
    for r in hzeros:
        assert b[r-1] != 0, (p, r, 'consecutive Hasse zeros')
        assert (xi[r] == 0) == (k[r] == 0), (p, r)
        value = int(xi[r])
        if value in seen:
            collisions.append((seen[value], r, value))
        else:
            seen[value] = r
        if xi[r] == 0:
            commons.append(r)
            common_records.append((p, r))

    # Check the exact two-endpoint Green/Duhamel law for every Hasse-zero pair.
    for i in range(len(hzeros)):
        for j in range(i+1, len(hzeros)):
            check_two_endpoint_shooting(p, b, g, k, xi, hzeros[i], hzeros[j])

    if collisions:
        xi_collision_records.append((p, collisions))
    if commons:
        common_by_prime[p] = tuple(commons)
    mult_hist[len(commons)] += 1

# Mandatory guard: Xi is NOT injective on the Hasse-zero set.
p = 41
b41 = apery_mod(p)
g41 = g_mod(p, b41)
xi41 = xi_mod(p, b41, g41)
assert b41[10] == 0 and b41[30] == 0
assert xi41[10] == 7 and xi41[30] == 7
assert xi41[10] != 0

repeated = [(p, rows) for p, rows in sorted(common_by_prime.items()) if len(rows) >= 2]
max_common = max([len(rows) for rows in common_by_prime.values()] + [0])

print('Q7696_SCAN_PMAX', PMAX)
print('PRIMES_WITH_HASSE_ZEROS', primes_with_hasse_zeros)
print('MAX_HASSE_ZERO_COUNT', hasse_zero_count_max)
print('MAX_HASSE_ZERO_COUNT_PRIMES_FIRST20', hzero_count_argmax[:20])
print('COMMON_RECORDS', common_records)
print('COMMON_BY_PRIME', sorted(common_by_prime.items()))
print('COMMON_MULTIPLICITY_HISTOGRAM', sorted(mult_hist.items()))
print('MAX_COMMON_MULTIPLICITY', max_common)
print('REPEATED_COMMON_PRIMES', repeated)
print('XI_COLLISION_PRIME_COUNT', len(xi_collision_records))
print('FIRST_XI_COLLISIONS', xi_collision_records[:20])
print('GUARD_P41', (10, int(xi41[10])), (30, int(xi41[30])))
print('Q7696_EXACT_SCAN=PASS')
