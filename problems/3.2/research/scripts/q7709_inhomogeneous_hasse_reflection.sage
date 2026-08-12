#!/usr/bin/env sage -python
"""Q7709 exact inhomogeneous Hasse-reflection audit.

Run with:
  sage -python problems/3.2/research/scripts/q7709_inhomogeneous_hasse_reflection.sage

The source is the actual repo source g=F^{-2}(1-34t+t^2)^{-1/2}.
No reflected-depth equivalence is used.
"""

from sage.all import GF, PolynomialRing, PowerSeriesRing

PRIMES = (17, 41, 181, 2237)


def A(r):
    return 34*r**3 + 51*r**2 + 27*r + 5


def homogeneous(p, y0, y1):
    K, N = GF(p), p-1
    y = [K(0)]*(N+1)
    y[0], y[1] = K(y0), K(y1)
    for r in range(1, N):
        y[r+1] = (K(A(r))*y[r] - K(r**3)*y[r-1]) / K((r+1)**3)
    return y


def make_q(p):
    K, N = GF(p), p-1
    q = [K(0)]*(N+1)
    q[0], q[1] = K(1), K(17)
    for n in range(1, N):
        q[n+1] = (K(17*(2*n+1))*q[n] - K(n)*q[n-1]) / K(n+1)
    return q


def make_g(p, b, q):
    K, N = GF(p), p-1
    PS = PowerSeriesRing(K, 'z', default_prec=p)
    z = PS.gen()
    F = PS(sum(b[j]*z**j for j in range(N+1)))
    Q = PS(sum(q[j]*z**j for j in range(N+1)))
    G = Q/(F*F)
    return [K(G[j]) for j in range(N+1)]


def make_kappa(p, g):
    K, N = GF(p), p-1
    k = [K(0)]*(N+1)
    k[0], k[1] = K(0), K(-36)
    for r in range(1, N):
        k[r+1] = (
            K(A(r))*k[r] - K(r**3)*k[r-1] - K(5)*g[r+1]
        ) / K((r+1)**3)
    return k


def J(y):
    return list(reversed(y))


def Lrow(p, y, r):
    K, N = GF(p), p-1
    if r == 0:
        return y[1] - K(5)*y[0]
    if r == N:
        return K(5)*y[N] - y[N-1]
    return K((r+1)**3)*y[r+1] - K(A(r))*y[r] + K(r**3)*y[r-1]


def W(p, h, y, r):
    K = GF(p)
    return K(r**3)*(h[r-1]*y[r] - h[r]*y[r-1])


def rev(poly, d, t):
    return sum(poly[i]*t**(d-i) for i in range(d+1))


def support(poly):
    if not poly:
        return []
    return [i for i in range(poly.degree()+1) if poly[i] != 0]


def audit(p):
    K, N = GF(p), p-1
    b = homogeneous(p, 1, 5)
    v = homogeneous(p, 0, 1)
    q = make_q(p)
    g = make_g(p, b, q)
    k = make_kappa(p, g)

    # Exact operator reflection: A(-r-1)=-A(r), hence LJ=-JL on interior.
    probe = [K((j+1)*(j+3)+7) for j in range(N+1)]
    for r in range(1, N):
        assert K(A(-r-1)) == -K(A(r))
        assert Lrow(p, J(probe), r) == -Lrow(p, probe, N-r)
        assert Lrow(p, b, r) == 0 and Lrow(p, v, r) == 0

    # Q7702 is pointwise fixedness of the full 2-dimensional kernel.
    # b,v are independent at rows 0,1 and both are J-fixed.
    assert all(b[N-r] == b[r] for r in range(N+1))
    assert all(v[N-r] == v[r] for r in range(N+1))
    assert b[N] == 1
    assert b[0]*v[1] - b[1]*v[0] == 1
    # Reflection preserves the constant Wronskian, so det(J|ker L)=+1.
    for r in range(1, N):
        assert W(p, J(b), J(v), r) == W(p, b, v, N-r+1)
    # There is therefore no anti-invariant homogeneous/adjoint direction.

    # Actual source normalization and left endpoint defect.
    assert g[0] == 1 and g[1] == 7
    assert k[0] == 0 and k[1] == K(-36)
    assert Lrow(p, k, 0) == K(-36)
    assert Lrow(p, k, 0) - K(-5)*g[1] == K(-1)
    for r in range(1, N):
        assert Lrow(p, k, r) == K(-5)*g[r+1]

    # Xi is exactly the b-Wronskian.
    Xi = [K(0)]*(N+1)
    Xi[0] = K(-1)
    for r in range(1, N+1):
        Xi[r] = Xi[r-1] - K(5)*g[r]*b[r-1]
        assert Xi[r] == W(p, b, k, r)
    assert Xi[1] == K(-36)

    common = []
    for r in range(1, N):
        if b[r] == 0:
            assert b[r-1] != 0
            assert (k[r] == 0) == (Xi[r] == 0)
            if k[r] == 0:
                common.append(r)

    # C=Jk-k: anti-invariant C, J-even forcing.
    C = [k[N-r]-k[r] for r in range(N+1)]
    assert all(C[N-r] == -C[r] for r in range(N+1))
    for r in range(1, N):
        rhs = K(5)*(g[p-r] + g[r+1])
        assert Lrow(p, C, r) == rhs
        assert rhs == K(5)*(g[p-(N-r)] + g[N-r+1])
    boundary_period = C[1] - K(5)*C[0]
    assert Lrow(p, C, 0) == boundary_period
    assert Lrow(p, C, N) == boundary_period

    # Exact rational Hasse reciprocity and the Taylor tail.
    R = PolynomialRing(K, 't')
    t = R.gen()
    B = sum(b[j]*t**j for j in range(N+1))
    H = (1-K(34)*t+t**2)**(N//2)
    G = sum(g[j]*t**j for j in range(N+1))
    assert rev(B, N, t) == B and rev(H, N, t) == H
    num = H - B**2*G
    assert all(num[i] == 0 for i in range(N+1))
    T, rem = num.quo_rem(t**(N+1))
    assert rem == 0
    Tdag = rev(T, 2*N-1, t)
    JG = rev(G, N, t)
    assert B**2*JG == t**(2*N)*H - Tdag
    Sfull = (G-g[0])//t + t*JG
    assert t*B**2*Sfull == (
        (1+t**(2*N+2))*H - g[0]*B**2 - t**(N+1)*T - t**2*Tdag
    )
    for r in range(1, N):
        assert Sfull[r] == g[r+1] + g[p-r]
    src_support = [r for r in range(1, N) if Sfull[r] != 0]
    tail_support = support(T)

    # Rank-two variation coordinates (Xi,Phi), using independent v.
    # v is also J-fixed; it is not an anti-invariant escape hatch.
    Phi = [None] + [W(p, v, k, r) for r in range(1, N+1)]
    for r in range(1, N):
        assert Phi[r+1]-Phi[r] == v[r]*K(-5)*g[r+1]

    rows = []
    for r in common:
        assert Xi[r] == 0 and b[r] == 0 and k[r] == 0
        # Since b_r=0 and b,v are independent, v_r cannot vanish.
        assert v[r] != 0
        assert Phi[r] == -K(r**3)*v[r]*k[r-1]
        rr = N-r
        rows.append({
            'r': r,
            'reflected_r': rr,
            'kappa_reflected': int(k[rr]),
            'C_at_r': int(C[r]),
            'Phi_at_r': int(Phi[r]),
            'kappa_prev': int(k[r-1]),
        })

    # If two common rows existed, Xi endpoints vanish and give Green/Duhamel;
    # Phi endpoints remain free slope data rather than a second return equation.
    pairs = []
    for i in range(len(common)):
        for j in range(i+1, len(common)):
            r, s = common[i], common[j]
            green = sum(b[u]*K(-5)*g[u+1] for u in range(r, s))
            transverse = sum(v[u]*K(-5)*g[u+1] for u in range(r, s))
            assert green == Xi[s]-Xi[r] == 0
            assert transverse == Phi[s]-Phi[r]
            pairs.append((r, s, int(Phi[r]), int(Phi[s])))

    result = {
        'p': p,
        'common': common,
        'boundary_period': int(boundary_period),
        'source_support': len(src_support),
        'source_rows': N-1,
        'tail_degree': int(T.degree()),
        'tail_support': len(tail_support),
        'b_zero_count': sum(1 for r in range(1, N) if b[r] == 0),
    }
    print('P', p)
    print('  common_rows', rows)
    print('  boundary_period', result['boundary_period'])
    print('  source_support', result['source_support'], 'of', result['source_rows'])
    print('  tail_degree', result['tail_degree'], 'tail_support', result['tail_support'])
    print('  b_zero_count', result['b_zero_count'])
    print('  pair_checks', pairs)
    return result


def main():
    out = [audit(p) for p in PRIMES]
    # Directly recomputed known common rows; no reflected-depth implication.
    assert 13 in out[0]['common']
    assert 492 in out[-1]['common']
    print('SUMMARY', out)
    print('Q7709_SAGE_VERIFY PASS')


if __name__ == '__main__':
    main()
