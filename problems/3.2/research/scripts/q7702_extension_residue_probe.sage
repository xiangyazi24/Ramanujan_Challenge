from sage.all import *
import os

PMAX = int(os.environ.get('Q7702_RESIDUE_PMAX', '5000'))
ABOVE_COUNT = int(os.environ.get('Q7702_ABOVE_COUNT', '12'))


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


def kappa_mod(p, source):
    Fp = GF(p)
    k = [Fp(0), Fp(-36)]
    for n in range(2, p):
        k.append((Fp(P(n-1))*k[n-1] - Fp((n-1)^3)*k[n-2] - Fp(5)*source[n]) / Fp(n^3))
    return k


def xi_mod(p, b, source):
    Fp = GF(p)
    xi = [Fp(-1)]
    for n in range(1,p):
        xi.append(xi[-1] - Fp(5)*source[n]*b[n-1])
    return xi


def u_mod(p):
    Fp = GF(p)
    u = [Fp(0), Fp(1)]
    for n in range(2,p):
        u.append((Fp(P(n-1))*u[n-1] - Fp((n-1)^3)*u[n-2]) / Fp(n^3))
    return u


def phi_mod(p, source, u):
    Fp = GF(p)
    phi = [Fp(0)]
    for n in range(1,p):
        phi.append(phi[-1] + Fp(5)*source[n]*u[n-1])
    return phi


def audit_prime(p, full_defect=True):
    Fp = GF(p)
    b = apery_mod(p)
    g = g_mod(p,b)
    k = kappa_mod(p,g)
    xi = xi_mod(p,b,g)
    u = u_mod(p)
    phi = phi_mod(p,g,u)

    assert all(b[p-1-r] == b[r] for r in range(p))
    assert all(u[p-1-r] == u[r] for r in range(p))

    for n in range(1,p):
        assert Fp(n^3)*(b[n-1]*u[n]-b[n]*u[n-1]) == 1
        assert k[n] == xi[n]*u[n] + phi[n]*b[n]

    # Exact truncated extension-class residue:
    # L K_<p = -5(G_<p-1) - t - Xi_{p-1} t^p.
    if full_defect:
        for n in range(0,p+2):
            lhs = Fp(0)
            if n < p:
                lhs += Fp(n^3)*k[n]
            if 1 <= n and n-1 < p:
                lhs -= Fp(P(n-1))*k[n-1]
            if 2 <= n and n-2 < p:
                lhs += Fp((n-1)^3)*k[n-2]
            if n == 0:
                rhs = Fp(0)
            elif n == 1:
                rhs = -Fp(5)*g[1] - Fp(1)
            elif 2 <= n < p:
                rhs = -Fp(5)*g[n]
            elif n == p:
                rhs = -xi[p-1]
            else:
                rhs = Fp(0)
            assert lhs == rhs, (p,n,lhs,rhs)

    hzeros = [r for r in range(1,p) if b[r] == 0]
    commons = [r for r in hzeros if xi[r] == 0]
    return b,g,k,xi,u,phi,hzeros,commons


def force_rows_and_restore_terminal_periods(p, rows):
    """Formal-source no-go: force selected Hasse rows common while restoring
    both terminal variation periods (Xi_{p-1}, Phi_{p-1}) exactly.
    The canonical initial line kappa_0=0,kappa_1=-36 is never changed.
    """
    Fp = GF(p)
    b = apery_mod(p)
    u = u_mod(p)
    source = g_mod(p,b)
    base_xi = xi_mod(p,b,source)
    base_phi = phi_mod(p,source,u)
    assert all(b[r] == 0 for r in rows)
    assert max(rows) <= p-3  # b_{p-2}=5, b_{p-1}=1 for p != 5.

    forced = []
    for r in rows:
        k = kappa_mod(p,source)
        delta = Fp(r^3)*k[r]/Fp(5)
        source[r] += delta
        forced.append((r,int(delta)))
        assert kappa_mod(p,source)[r] == 0

    xi = xi_mod(p,b,source)
    phi = phi_mod(p,source,u)
    a, c = p-2, p-1
    response = matrix(Fp, [
        [-Fp(5)*b[a-1], -Fp(5)*b[c-1]],
        [ Fp(5)*u[a-1],  Fp(5)*u[c-1]],
    ])
    # det = 25/8 by the unit Casoratian at n=p-2, hence a p-unit for p>=7.
    assert response.det() == Fp(25)/Fp(8)
    rhs = vector(Fp, [base_xi[p-1]-xi[p-1], base_phi[p-1]-phi[p-1]])
    tail = response.solve_right(rhs)
    source[a] += tail[0]
    source[c] += tail[1]

    k2 = kappa_mod(p,source)
    xi2 = xi_mod(p,b,source)
    phi2 = phi_mod(p,source,u)
    assert all(k2[r] == 0 for r in rows)
    assert xi2[p-1] == base_xi[p-1]
    assert phi2[p-1] == base_phi[p-1]
    return forced, (int(tail[0]),int(tail[1])), int(base_xi[p-1]), int(base_phi[p-1])


residue_zero = []
residue_one = []
residue_minus_one = []
residue_with_common = []
for pz in prime_range(7, PMAX+1):
    p = int(pz)
    b,g,k,xi,u,phi,hz,commons = audit_prime(p, full_defect=True)
    e = int(xi[p-1])
    if e == 0 and len(residue_zero) < 20:
        residue_zero.append((p, tuple(commons), int(phi[p-1])))
    if e == 1 and len(residue_one) < 20:
        residue_one.append((p, tuple(commons), int(phi[p-1])))
    if e == p-1 and len(residue_minus_one) < 20:
        residue_minus_one.append((p, tuple(commons), int(phi[p-1])))
    if commons and len(residue_with_common) < 20:
        residue_with_common.append((p,tuple(commons),e,int(phi[p-1])))

# p=41: actual-source b-weighted period is zero, independent u-period is not.
p=41
b,g,k,xi,u,phi,hz,commons = audit_prime(p)
assert b[10] == b[30] == 0 and xi[10] == xi[30] == 7
b_period = sum((g[m]*b[m-1] for m in range(11,31)), GF(p)(0))
u_period = sum((g[m]*u[m-1] for m in range(11,31)), GF(p)(0))
phi_diff41 = phi[30]-phi[10]
assert b_period == 0
assert phi_diff41 == GF(p)(5)*u_period
assert u_period != 0

# Exact terminal-compensation example: p=19 has Hasse zeros 8,10.  Force both
# common, then restore BOTH global period coordinates with rows 17,18.
comp19 = force_rows_and_restore_terminal_periods(19, (8,10))

# Targeted primes strictly above the old p<=20000 uniqueness scan.  These test
# the mechanism (frame + extension defect), not a statistical uniqueness claim.
above = []
pz = next_prime(20000)
for _ in range(ABOVE_COUNT):
    p = int(pz)
    b,g,k,xi,u,phi,hz,commons = audit_prime(p, full_defect=True)
    above.append((p, len(hz), tuple(commons), int(xi[p-1]), int(phi[p-1])))
    pz = next_prime(p)

print('Q7702_RESIDUE_PMAX', PMAX)
print('RESIDUE_ZERO_FIRST20', residue_zero)
print('RESIDUE_ONE_FIRST20', residue_one)
print('RESIDUE_MINUS_ONE_FIRST20', residue_minus_one)
print('RESIDUE_AT_COMMON_PRIMES', residue_with_common)
print('P41_B_PERIOD', int(b_period), 'P41_U_PERIOD', int(u_period), 'P41_PHI_DIFF', int(phi_diff41))
print('P19_FORCE_AND_RESTORE', comp19)
print('ABOVE_20000_MECHANISM_PRIMES', above)
print('Q7702_EXTENSION_RESIDUE=PASS')
