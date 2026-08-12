from sage.all import *
import os

PMAX = int(os.environ.get('Q7702_PMAX', '5000'))


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


def xi_mod(p, b, g):
    Fp = GF(p)
    xi = [Fp(-1)]
    for n in range(1, p):
        xi.append(xi[-1] - Fp(5)*g[n]*b[n-1])
    return xi


def u_mod(p):
    Fp = GF(p)
    u = [Fp(0), Fp(1)]
    for n in range(2, p):
        u.append((Fp(P(n-1))*u[n-1] - Fp((n-1)^3)*u[n-2]) / Fp(n^3))
    return u


q_sym_fail = []
xi_affine_fail = []
u_reflect_fail = []
examples = []
for pz in prime_range(7, PMAX+1):
    p = int(pz)
    Fp = GF(p)
    b = apery_mod(p)
    g = g_mod(p, b)
    xi = xi_mod(p, b, g)
    u = u_mod(p)
    # Homogeneous reflection is exact: J u = A_p b - u.
    A = u[p-1]
    ok_u = all(u[p-1-r] == A*b[r] - u[r] for r in range(p))
    if not ok_u:
        u_reflect_fail.append(p)
        break

    # Candidate actual-source reciprocal increment law.
    ok_q = all(g[p-m]*b[m] == g[m]*b[m-1] for m in range(1,p))
    if not ok_q:
        first = next(m for m in range(1,p) if g[p-m]*b[m] != g[m]*b[m-1])
        q_sym_fail.append((p, first, int(g[p-first]*b[first]), int(g[first]*b[first-1])))
    # Consequent affine partner law, tested independently.
    C = xi[p-1] - 1
    ok_xi = all(xi[r] + xi[p-1-r] == C for r in range(p))
    if not ok_xi:
        first = next(r for r in range(p) if xi[r] + xi[p-1-r] != C)
        xi_affine_fail.append((p, first, int(xi[r]), int(xi[p-1-r]), int(C)))

    hz = [r for r in range(1,p) if b[r] == 0]
    if len(hz) >= 4 and len(examples) < 12:
        pairs = []
        for r in hz:
            s = p-1-r
            if r < s:
                pairs.append((r,s,int(xi[r]),int(xi[s]),int(xi[r]+xi[s])))
        examples.append((p, int(xi[p-1]), int(C), pairs))

print('Q7702_PMAX', PMAX)
print('U_REFLECTION_FAILURES', u_reflect_fail[:5])
print('Q_INCREMENT_SYMMETRY_FAILURE_COUNT', len(q_sym_fail))
print('FIRST_Q_INCREMENT_SYMMETRY_FAILURES', q_sym_fail[:10])
print('XI_AFFINE_PARTNER_FAILURE_COUNT', len(xi_affine_fail))
print('FIRST_XI_AFFINE_PARTNER_FAILURES', xi_affine_fail[:10])
print('MULTI_HASSE_PARTNER_EXAMPLES', examples)

# Locked guards.
p=41
b=apery_mod(p); g=g_mod(p,b); xi=xi_mod(p,b,g)
print('P41', int(xi[10]), int(xi[30]), int(xi[40]), int(xi[10]+xi[30]))
p=17
b=apery_mod(p); g=g_mod(p,b); xi=xi_mod(p,b,g)
print('P17_HZ', [(r,int(xi[r])) for r in range(1,p) if b[r]==0], 'XI_LAST', int(xi[p-1]))
print('Q7702_PARTNER_PROBE=PASS')
