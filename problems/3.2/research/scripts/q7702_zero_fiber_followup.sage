from sage.all import *
import os

PMAX = int(os.environ.get('Q7702_PMAX', '5000'))
ABOVE_COUNT = int(os.environ.get('Q7702_ABOVE_COUNT', '12'))


def P(n):
    n = ZZ(n)
    return 34*n^3 + 51*n^2 + 27*n + 5


def apery_mod(p):
    Fp = GF(p)
    b = [Fp(1), Fp(5)]
    for n in range(2,p):
        b.append((Fp(P(n-1))*b[n-1]-Fp((n-1)^3)*b[n-2])/Fp(n^3))
    return b


def source_mod(p,b):
    Fp = GF(p)
    R = PowerSeriesRing(Fp,'t',default_prec=p)
    q = [Fp(1)]
    if p > 1:
        q.append(Fp(17))
    for n in range(1,p-1):
        q.append((Fp(17*(2*n+1))*q[n]-Fp(n)*q[n-1])/Fp(n+1))
    F = R(b).add_bigoh(p)
    Q = R(q).add_bigoh(p)
    G = (Q/(F*F)).add_bigoh(p)
    return [Fp(G[n]) for n in range(p)]


def kappa_mod(p,source):
    Fp = GF(p)
    k = [Fp(0),Fp(-36)]
    for n in range(2,p):
        k.append((Fp(P(n-1))*k[n-1]-Fp((n-1)^3)*k[n-2]-Fp(5)*source[n])/Fp(n^3))
    return k


def xi_mod(p,b,source):
    Fp = GF(p)
    xi = [Fp(-1)]
    for n in range(1,p):
        xi.append(xi[-1]-Fp(5)*source[n]*b[n-1])
    return xi


def u_mod(p):
    Fp = GF(p)
    u = [Fp(0),Fp(1)]
    for n in range(2,p):
        u.append((Fp(P(n-1))*u[n-1]-Fp((n-1)^3)*u[n-2])/Fp(n^3))
    return u


def phi_mod(p,source,u):
    Fp = GF(p)
    phi = [Fp(0)]
    for n in range(1,p):
        phi.append(phi[-1]+Fp(5)*source[n]*u[n-1])
    return phi


def audit_prime(p, check_defect=True):
    Fp = GF(p)
    b = apery_mod(p)
    g = source_mod(p,b)
    k = kappa_mod(p,g)
    xi = xi_mod(p,b,g)
    u = u_mod(p)
    phi = phi_mod(p,g,u)

    # Full homogeneous Hasse reciprocity in the normalized basis.
    assert all(b[p-1-r] == b[r] for r in range(p))
    assert all(u[p-1-r] == u[r] for r in range(p))

    for n in range(1,p):
        # Unit Casoratian and variation-of-parameters frame.
        assert Fp(n^3)*(b[n-1]*u[n]-b[n]*u[n-1]) == 1
        assert k[n] == xi[n]*u[n]+phi[n]*b[n]

    if check_defect:
        # L K_<p = -5(G_<p-1)-t-Xi_{p-1}t^p.
        for n in range(p+2):
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
                rhs = -Fp(5)*g[1]-Fp(1)
            elif n < p:
                rhs = -Fp(5)*g[n]
            elif n == p:
                rhs = -xi[p-1]
            else:
                rhs = Fp(0)
            assert lhs == rhs, (p,n,lhs,rhs)

    hz = [r for r in range(1,p) if b[r] == 0]
    commons = [r for r in hz if xi[r] == 0]
    return b,g,k,xi,u,phi,hz,commons


def force_and_restore(p, rows):
    """Formal-source obstruction preserving canonical initial line and both
    global variation periods while forcing the selected Hasse rows common."""
    Fp = GF(p)
    b = apery_mod(p)
    u = u_mod(p)
    source = source_mod(p,b)
    base_xi = xi_mod(p,b,source)
    base_phi = phi_mod(p,source,u)
    assert all(b[r] == 0 for r in rows)
    assert max(rows) <= p-3

    forced = []
    for r in rows:
        k = kappa_mod(p,source)
        delta = Fp(r^3)*k[r]/Fp(5)
        source[r] += delta
        forced.append((r,int(delta)))
        assert kappa_mod(p,source)[r] == 0

    xi = xi_mod(p,b,source)
    phi = phi_mod(p,source,u)
    a,c = p-2,p-1
    M = matrix(Fp,[
        [-Fp(5)*b[a-1], -Fp(5)*b[c-1]],
        [ Fp(5)*u[a-1],  Fp(5)*u[c-1]],
    ])
    assert M.det() == Fp(25)/Fp(8)
    rhs = vector(Fp,[base_xi[p-1]-xi[p-1],base_phi[p-1]-phi[p-1]])
    tail = M.solve_right(rhs)
    source[a] += tail[0]
    source[c] += tail[1]

    k2 = kappa_mod(p,source)
    xi2 = xi_mod(p,b,source)
    phi2 = phi_mod(p,source,u)
    assert all(k2[r] == 0 for r in rows)
    assert xi2[p-1] == base_xi[p-1]
    assert phi2[p-1] == base_phi[p-1]
    return forced,(int(tail[0]),int(tail[1])),int(base_xi[p-1]),int(base_phi[p-1])


# p=5 endpoint case: Xi_n=-1 identically mod 5, hence no common row.
assert GF(5)(-1) != 0

q_partner_fail = []
xi_partner_fail = []
reflection_pair_examples = []
residue_zero = []
residue_one = []
residue_common = []
for pz in prime_range(7,PMAX+1):
    p = int(pz)
    b,g,k,xi,u,phi,hz,commons = audit_prime(p,True)

    # Candidate reciprocal coefficient law is false for the actual source.
    bad_q = [m for m in range(1,p) if g[p-m]*b[m] != g[m]*b[m-1]]
    if bad_q:
        m = bad_q[0]
        q_partner_fail.append((p,m,int(g[p-m]*b[m]),int(g[m]*b[m-1])))

    C = xi[p-1]-1
    bad_x = [r for r in range(p) if xi[r]+xi[p-1-r] != C]
    if bad_x:
        r = bad_x[0]
        xi_partner_fail.append((p,r,int(xi[r]),int(xi[p-1-r]),int(C)))

    if len(hz) >= 4 and len(reflection_pair_examples) < 12:
        pairs=[]
        for r in hz:
            s=p-1-r
            if r<s:
                pairs.append((r,s,int(xi[r]),int(xi[s]),int(xi[r]+xi[s])))
        reflection_pair_examples.append((p,pairs))

    e=int(xi[p-1])
    if e==0 and len(residue_zero)<20:
        residue_zero.append((p,tuple(commons),int(phi[p-1])))
    if e==1 and len(residue_one)<20:
        residue_one.append((p,tuple(commons),int(phi[p-1])))
    if commons and len(residue_common)<20:
        residue_common.append((p,tuple(commons),e,int(phi[p-1])))

# Exact actual-source obstruction to "second period also vanishes": p=41.
p=41
b,g,k,xi,u,phi,hz,commons=audit_prime(p,True)
assert b[10]==b[30]==0 and xi[10]==xi[30]==7
B41=sum((g[m]*b[m-1] for m in range(11,31)),GF(p)(0))
U41=sum((g[m]*u[m-1] for m in range(11,31)),GF(p)(0))
assert B41==0 and U41!=0
assert phi[30]-phi[10]==GF(p)(5)*U41
phi41=phi[30]-phi[10]

# Same-prime exact obstruction to a constant reflection-pair Xi sum.
p=181
b,g,k,xi,u,phi,hz,commons=audit_prime(p,True)
assert b[19]==b[161]==b[47]==b[133]==0
pair_sums_181=(int(xi[19]+xi[161]),int(xi[47]+xi[133]))
assert pair_sums_181[0] != pair_sums_181[1]

# Formal-source obstruction preserving both terminal periods.
comp19=force_and_restore(19,(8,10))

# Targeted mechanism test strictly beyond old p<=20000 uniqueness scan.
above=[]
pz=next_prime(20000)
for _ in range(ABOVE_COUNT):
    p=int(pz)
    b,g,k,xi,u,phi,hz,commons=audit_prime(p,True)
    above.append((p,len(hz),tuple(commons),int(xi[p-1]),int(phi[p-1])))
    pz=next_prime(p)

print('Q7702_PMAX',PMAX)
print('RECIPROCAL_Q_FAILURE_COUNT',len(q_partner_fail))
print('FIRST_RECIPROCAL_Q_FAILURES',q_partner_fail[:10])
print('AFFINE_XI_PARTNER_FAILURE_COUNT',len(xi_partner_fail))
print('FIRST_AFFINE_XI_PARTNER_FAILURES',xi_partner_fail[:10])
print('REFLECTION_PAIR_EXAMPLES',reflection_pair_examples[:3])
print('P181_PAIR_SUMS',pair_sums_181)
print('P41_B_PERIOD',int(B41),'P41_U_PERIOD',int(U41),'P41_PHI_DIFF',int(phi41))
print('RESIDUE_ZERO_FIRST20',residue_zero)
print('RESIDUE_ONE_FIRST20',residue_one)
print('RESIDUE_AT_COMMON_PRIMES',residue_common)
print('P19_FORCE_AND_RESTORE',comp19)
print('ABOVE_20000_MECHANISM_PRIMES',above)
print('Q7702_ZERO_FIBER_FOLLOWUP=PASS')
