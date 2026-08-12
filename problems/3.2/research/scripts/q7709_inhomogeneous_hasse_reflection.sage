#!/usr/bin/env sage -python
"""Q7709: exact inhomogeneous Hasse-reflection audit for the actual Apéry source.

This is ordinary Python syntax, run with Sage's Python:
  sage -python problems/3.2/research/scripts/q7709_inhomogeneous_hasse_reflection.sage

It verifies, over F_p for selected primes, all of the following directly from
repo definitions:
  * homogeneous anti-commutation L J = - J L and b_{p-1-r}=b_r;
  * g(t)=F(t)^(-2)(1-34t+t^2)^(-1/2) through degree p-1;
  * the canonical kappa boundary (0,-36) and row forcing -5 g_{r+1};
  * Xi = W(b,kappa), hence common-row equivalence;
  * C=J kappa-kappa has L C = 5(g_{p-r}+g_{r+1}) in the interior;
  * exact Hasse-rational reciprocity with its Taylor tail T_p;
  * an anti-invariant homogeneous variation coordinate Phi=W(a,kappa).

No reflected-depth equivalence is used anywhere.
"""

from sage.all import GF, PolynomialRing, PowerSeriesRing, ZZ

PRIMES = (17, 41, 181, 2237)


def A(r):
    return 34*r**3 + 51*r**2 + 27*r + 5


def make_b(p):
    K = GF(p)
    N = p - 1
    b = [K(0)] * (N + 1)
    b[0], b[1] = K(1), K(5)
    for r in range(1, N):
        b[r+1] = (K(A(r))*b[r] - K(r**3)*b[r-1]) / K((r+1)**3)
    return b


def make_homogeneous_v(p):
    """Independent interior homogeneous solution, v_0=0, v_1=1."""
    K = GF(p)
    N = p - 1
    v = [K(0)] * (N + 1)
    v[0], v[1] = K(0), K(1)
    for r in range(1, N):
        v[r+1] = (K(A(r))*v[r] - K(r**3)*v[r-1]) / K((r+1)**3)
    return v


def make_q(p):
    """Q=(1-34t+t^2)^(-1/2), via the repo recurrence."""
    K = GF(p)
    N = p - 1
    q = [K(0)] * (N + 1)
    q[0], q[1] = K(1), K(17)
    for n in range(1, N):
        q[n+1] = (K(17*(2*n+1))*q[n] - K(n)*q[n-1]) / K(n+1)
    return q


def make_g(p, b, q):
    """g=Q/F^2 in F_p[[t]]/(t^p)."""
    K = GF(p)
    N = p - 1
    PS = PowerSeriesRing(K, 'z', default_prec=p)
    z = PS.gen()
    F = PS(sum(b[j] * z**j for j in range(N+1)))
    Q = PS(sum(q[j] * z**j for j in range(N+1)))
    G = Q / (F*F)
    return [K(G[j]) for j in range(N+1)]


def make_kappa(p, g):
    """Canonical kappa: kappa_0=0, kappa_1=-36, L kappa=-5 g_{r+1}."""
    K = GF(p)
    N = p - 1
    k = [K(0)] * (N + 1)
    k[0], k[1] = K(0), K(-36)
    for r in range(1, N):
        k[r+1] = (
            K(A(r))*k[r] - K(r**3)*k[r-1] - K(5)*g[r+1]
        ) / K((r+1)**3)
    return k


def Lrow(p, y, r):
    K = GF(p)
    N = p - 1
    if r == 0:
        return y[1] - K(5)*y[0]
    if r == N:
        return K(5)*y[N] - y[N-1]
    return K((r+1)**3)*y[r+1] - K(A(r))*y[r] + K(r**3)*y[r-1]


def J(y):
    return list(reversed(y))


def wronskian(p, h, k, r):
    K = GF(p)
    if r == 0:
        raise ValueError('W_r is defined for r>=1')
    return K(r**3) * (h[r-1]*k[r] - h[r]*k[r-1])


def reverse_degree(poly, d, t):
    return sum(poly[i] * t**(d-i) for i in range(d+1))


def poly_support(poly):
    return [i for i in range(poly.degree()+1) if poly[i] != 0] if poly else []


def audit_prime(p):
    K = GF(p)
    N = p - 1
    b = make_b(p)
    v = make_homogeneous_v(p)
    q = make_q(p)
    g = make_g(p, b, q)
    k = make_kappa(p, g)

    # Q7702 homogeneous reflection, including the exact operator sign.
    assert A(-1) == -5
    for r in range(1, N):
        assert K(A(-r-1)) == -K(A(r))
        # test L J = - J L on an arbitrary independent homogeneous solution
        assert Lrow(p, J(v), r) == -Lrow(p, v, N-r)
        assert Lrow(p, b, r) == 0
        assert Lrow(p, v, r) == 0
    assert all(b[N-r] == b[r] for r in range(N+1))
    assert b[N] == 1

    # Actual source normalization: g0=1, g1=7 and left endpoint defect -1.
    assert g[0] == 1 and g[1] == 7
    assert k[0] == 0 and k[1] == K(-36)
    assert Lrow(p, k, 0) == K(-36)
    assert Lrow(p, k, 0) - K(-5)*g[1] == K(-1)
    for r in range(1, N):
        assert Lrow(p, k, r) == K(-5)*g[r+1]

    # Xi prefix and Wronskian identity.
    Xi = [K(0)] * (N + 1)
    Xi[0] = K(-1)
    for r in range(1, N+1):
        Xi[r] = Xi[r-1] - K(5)*g[r]*b[r-1]
    assert Xi[1] == K(-36)
    for r in range(1, N+1):
        assert Xi[r] == wronskian(p, b, k, r)

    common = []
    for r in range(1, N):
        if b[r] == 0:
            # b_{r-1} cannot also vanish for a nonzero second-order solution.
            assert b[r-1] != 0
            assert (k[r] == 0) == (Xi[r] == 0)
            if k[r] == 0:
                common.append(r)

    # Inhomogeneous cocycle C=Jk-k.  It is anti-invariant, its forcing is J-even.
    C = [J(k)[r] - k[r] for r in range(N+1)]
    assert all(C[N-r] == -C[r] for r in range(N+1))
    for r in range(1, N):
        rhs = K(5) * (g[p-r] + g[r+1])
        assert Lrow(p, C, r) == rhs
        assert rhs == K(5) * (g[p-(N-r)] + g[(N-r)+1])
    boundary_period = C[1] - K(5)*C[0]
    assert Lrow(p, C, 0) == boundary_period
    assert Lrow(p, C, N) == boundary_period

    # Actual rational reciprocity and its exact Taylor-tail correction.
    R = PolynomialRing(K, 't')
    t = R.gen()
    B = sum(b[j]*t**j for j in range(N+1))
    D = 1 - K(34)*t + t**2
    H = D**(N//2)
    G = sum(g[j]*t**j for j in range(N+1))
    assert reverse_degree(B, N, t) == B
    assert reverse_degree(H, N, t) == H

    numerator = H - B**2 * G
    # Taylor agreement at zero through t^N: the defect starts at t^(N+1)=t^p.
    assert all(numerator[i] == 0 for i in range(N+1))
    T, rem = numerator.quo_rem(t**(N+1))
    assert rem == 0
    Tdag = reverse_degree(T, 2*N-1, t)
    JG = reverse_degree(G, N, t)
    assert B**2 * JG == t**(2*N)*H - Tdag
    Sfull = (G - g[0]) // t + t*JG
    assert t*B**2*Sfull == (
        (1+t**(2*N+2))*H - g[0]*B**2 - t**(N+1)*T - t**2*Tdag
    )
    for r in range(1, N):
        assert Sfull[r] == g[r+1] + g[p-r]

    source = [g[r+1] + g[p-r] for r in range(1, N)]
    source_support = [r for r in range(1, N) if source[r-1] != 0]
    tail_support = poly_support(T)

    # Anti-invariant homogeneous coordinate a=v-Jv.  If this choice vanishes,
    # the + combination is independent of b and one can re-choose v; for these
    # primes the anti-invariant choice is nonzero and is checked explicitly.
    a = [v[r] - v[N-r] for r in range(N+1)]
    assert any(x != 0 for x in a)
    assert all(a[N-r] == -a[r] for r in range(N+1))
    assert all(Lrow(p, a, r) == 0 for r in range(1, N))
    Phi = [None] + [wronskian(p, a, k, r) for r in range(1, N+1)]
    for r in range(1, N):
        assert Phi[r+1] - Phi[r] == a[r] * K(-5) * g[r+1]

    common_rows = []
    for r in common:
        # At a common row Xi=0 but Phi records the free transverse slope.
        assert Xi[r] == 0 and b[r] == 0 and k[r] == 0
        assert a[r] != 0
        assert Phi[r] == -K(r**3)*a[r]*k[r-1]
        rr = N-r
        common_rows.append({
            'r': r,
            'reflected_r': rr,
            'kappa_reflected': int(k[rr]),
            'C_at_r': int(C[r]),
            'Phi_at_r': int(Phi[r]),
            'kappa_prev': int(k[r-1]),
        })

    # If two common rows existed, only the b-Wronskian has zero endpoints.
    # Check that the anti-invariant identity retains Phi endpoint data.
    pair_checks = []
    for i in range(len(common)):
        for j in range(i+1, len(common)):
            r, s = common[i], common[j]
            green_segment = sum(b[u] * K(-5) * g[u+1] for u in range(r, s))
            phi_segment = sum(a[u] * K(-5) * g[u+1] for u in range(r, s))
            assert green_segment == 0
            assert phi_segment == Phi[s] - Phi[r]
            pair_checks.append((r, s, int(Phi[r]), int(Phi[s])))

    print('P', p)
    print('  common', common_rows)
    print('  common_count', len(common))
    print('  boundary_period', int(boundary_period))
    print('  source_support', len(source_support), 'of', N-1)
    print('  tail_degree', T.degree(), 'tail_support', len(tail_support))
    print('  first_source_nonzero', source_support[:8])
    print('  pair_checks', pair_checks)
    print('  b_zero_count', sum(1 for r in range(1, N) if b[r] == 0))

    return {
        'p': p,
        'common': common,
        'boundary_period': int(boundary_period),
        'source_support': len(source_support),
        'tail_degree': int(T.degree()),
        'tail_support': len(tail_support),
    }


def main():
    rows = [audit_prime(p) for p in PRIMES]
    # The known actual common pairs must appear; this assertion uses direct
    # recurrence computation only and no reflected-depth equivalence.
    assert 13 in rows[0]['common']
    assert 492 in rows[-1]['common']
    print('SUMMARY', rows)
    print('Q7709_SAGE_VERIFY PASS')


if __name__ == '__main__':
    main()
