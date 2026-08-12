#!/usr/bin/env sage -python
"""Q7736 exact residual-affine transversality verifier.

For p > r the Apéry recurrence has companion matrices

    M_j = [[P(j)/(j+1)^3, -j^3/(j+1)^3], [1,0]].

Let T_1=I and T_{r+1}=M_r T_r.  The primitive Eichler coordinate satisfies

    K_{r+1}=M_r K_r + (-5*g_{r+1}/(r+1)^3,0)^t,

with K_1=(-36,0)^t.  Writing the accumulated translation as v_r gives

    b_r     = e1^t T_r (5,1)^t,
    kappa_r = e1^t (T_r (-36,0)^t + v_r),
    det T_r = r^(-3).

Thus (v_r,T_r) lies in the affine determinant coset
F_p^2 semidirect {T in GL_2(F_p): det T=r^(-3)}.

The ambient common-zero equations are, for T=[[a,b],[c,d]], v=(x,y),

    5*a+b = 0,
    -36*a+x = 0.

For every p >= 7 and every nonzero determinant d0, their intersection in the
5-dimensional affine determinant coset is smooth codimension two and has
exactly p^2(p-1) F_p-points, a fraction 1/(p(p+1)) of the coset.

This script checks the exact recurrence/matrix-coefficient identities on actual
Apéry rows, the locked pairs (17,13), (2237,492), the nontransverse control
(11,5), and exhaustively checks the ambient count/Jacobian rank for small p.
It does NOT assert that the actual nonautonomous row samples equidistribute in
the ambient affine group.
"""

from sage.all import GF, Matrix, PowerSeriesRing, vector


LOCKED = ((11, 5, False), (17, 13, True), (2237, 492, True))
SMALL_AMBIENT = (7, 11, 13)


def P(n):
    return 34*n**3 + 51*n**2 + 27*n + 5


def row_data(p, rmax):
    assert p >= 7 and 1 <= rmax < p
    K = GF(p)

    b = [K(0)] * (rmax + 1)
    b[0] = K(1)
    if rmax >= 1:
        b[1] = K(5)
    for n in range(1, rmax):
        b[n+1] = (
            K(P(n))*b[n] - K(n**3)*b[n-1]
        ) / K((n+1)**3)

    # Q=(1-34t+t^2)^(-1/2), (n+1)q_{n+1}=(34n+17)q_n-nq_{n-1}.
    q = [K(0)] * (rmax + 1)
    q[0] = K(1)
    if rmax >= 1:
        q[1] = K(17)
    for n in range(1, rmax):
        q[n+1] = (
            K(34*n + 17)*q[n] - K(n)*q[n-1]
        ) / K(n+1)

    PS = PowerSeriesRing(K, 't', default_prec=rmax + 1)
    t = PS.gen()
    F = PS(sum(b[n]*t**n for n in range(rmax + 1)))
    Q = PS(sum(q[n]*t**n for n in range(rmax + 1)))
    G = Q/(F*F)
    g = [K(G[n]) for n in range(rmax + 1)]

    kappa = [K(0)] * (rmax + 1)
    kappa[0] = K(0)
    if rmax >= 1:
        kappa[1] = K(-36)
    for r in range(1, rmax):
        kappa[r+1] = (
            K(P(r))*kappa[r] - K(r**3)*kappa[r-1] - K(5)*g[r+1]
        ) / K((r+1)**3)

    return K, b, g, kappa


def affine_path(p, rmax):
    K, b, g, kappa = row_data(p, rmax)
    B = vector(K, [5, 1])
    K0 = vector(K, [-36, 0])
    e1 = vector(K, [1, 0])

    T = Matrix.identity(K, 2)
    v = vector(K, [0, 0])

    rows = {}
    for r in range(1, rmax + 1):
        if r > 1:
            j = r - 1
            M = Matrix(K, [
                [K(P(j))/K((j+1)**3), -K(j**3)/K((j+1)**3)],
                [K(1), K(0)],
            ])
            source = vector(K, [-K(5)*g[r]/K(r**3), K(0)])
            v = M*v + source
            T = M*T

        assert T.det() == K(r)**(-3)
        bcoef = (e1*T*B)
        kcoef = e1*(T*K0 + v)
        assert bcoef == b[r]
        assert kcoef == kappa[r]

        a, bb = T[0,0], T[0,1]
        x = v[0]
        assert bcoef == K(5)*a + bb
        assert kcoef == K(-36)*a + x

        if bcoef == 0:
            # The first row of an invertible T cannot vanish. Since bb=-5a,
            # this forces a != 0.  At a common zero the translation is pinned.
            assert a != 0
            if kcoef == 0:
                assert x == K(36)*a

        rows[r] = (T, v, bcoef, kcoef)

    return K, b, g, kappa, rows


def locked_checks():
    for p, r, should_common in LOCKED:
        K, b, g, kappa, rows = affine_path(p, r)
        hb = (b[r] == 0)
        hk = (kappa[r] == 0)
        if p == 11:
            assert hb and not hk
        else:
            assert hb and hk == should_common
        T, v, bcoef, kcoef = rows[r]
        print('LOCKED', p, r,
              'b0', int(bcoef == 0),
              'k0', int(kcoef == 0),
              'det', int(T.det()),
              'a', int(T[0,0]),
              'v0', int(v[0]))


def jacobian_rank(K, a, b, c, d, x, y):
    # Equations det(T)-d0, 5a+b, -36a+x.
    J = Matrix(K, [
        [d, -c, -b, a, 0, 0],
        [5, 1, 0, 0, 0, 0],
        [-36, 0, 0, 0, 1, 0],
    ])
    return J.rank()


def ambient_exhaustive(p):
    K = GF(p)
    d0 = K(3)  # any nonzero determinant coset has the same count
    total = 0
    first_zero = 0
    common = 0
    smooth = 0

    els = list(K)
    for a in els:
        for b in els:
            for c in els:
                for d in els:
                    if a*d - b*c != d0:
                        continue
                    for x in els:
                        for y in els:
                            total += 1
                            f = K(5)*a + b
                            h = K(-36)*a + x
                            if f == 0:
                                first_zero += 1
                            if f == 0 and h == 0:
                                common += 1
                                assert a != 0
                                assert jacobian_rank(K, a, b, c, d, x, y) == 3
                                smooth += 1

    expected_total = p**3 * (p**2 - 1)
    expected_first = p**3 * (p - 1)
    expected_common = p**2 * (p - 1)
    assert total == expected_total
    assert first_zero == expected_first
    assert common == expected_common
    assert smooth == common
    # Inside the determinant coset, the two conditions have exact density
    # 1/(p(p+1)).
    assert total == common * p * (p + 1)
    print('AMBIENT', p,
          'total', total,
          'first_zero', first_zero,
          'common', common,
          'ratio_den', p*(p+1),
          'smooth', smooth)


def order_of_teich_exponent(p, r):
    from sage.all import gcd
    return (p - 1) // gcd(p - 1, r)


def bounded_order_identity_check(R, D):
    # Exact combinatorial identity behind the bounded-order sector:
    # ord(omega_p^r)=d <= D implies p=d*g+1 for g | r.
    from sage.all import is_prime, divisors
    checked = 0
    for r in range(R + 1, 2*R + 1):
        for d in range(1, D + 1):
            for g in divisors(r):
                p = d*g + 1
                if p <= 2*R or not is_prime(p) or r >= p:
                    continue
                order = order_of_teich_exponent(p, r)
                if order <= D:
                    assert (p - 1) % g == 0
                    gg = (p - 1) // order
                    assert r % gg == 0
                    assert p == order*gg + 1
                    checked += 1
    print('BOUNDED_ORDER', 'R', R, 'D', D, 'checked', checked)


def main():
    locked_checks()
    for p in SMALL_AMBIENT:
        ambient_exhaustive(p)
    bounded_order_identity_check(80, 8)
    print('Q7736_AFFINE_TRANSVERSALITY PASS')


if __name__ == '__main__':
    main()
