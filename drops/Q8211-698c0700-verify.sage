#!/usr/bin/env sage
"""Exact QQ/ZZ/GF(p) verifier for Q8211; no floating arithmetic.

Run:
    sage drops/Q8211-698c0700-verify.sage
Finite loops are regression checks only. The accompanying report proves the
identities for all parameters.
"""
from argparse import ArgumentParser
from sage.all import (GF, QQ, ZZ, PolynomialRing, binomial, diagonal_matrix,
                      factorial, identity_matrix, matrix, prime_range, prod,
                      vector)


def la(j):
    return ZZ(j) * ZZ(j + 1)


def u(n, k):
    if k < 0 or k > n:
        return ZZ(0)
    return ZZ(binomial(n, k) * binomial(n + k, k))


def invu(n, k):
    if k < 0 or k > n:
        return QQ(0)
    return (QQ((-1) ** (n-k) * (2*k+1) * factorial(n) ** 2)
            / QQ(factorial(n-k) * factorial(n+k+1)))


def build(nmax, base=QQ):
    R = PolynomialRing(base, 'Y')
    Y = R.gen()
    phi = [R.one()]
    for k in range(1, nmax+1):
        phi.append(phi[-1] * (Y-base((k-1)*k)) / base(k*k))
    rac = [sum(base(u(n, k))*phi[k] for k in range(n+1))
           for n in range(nmax+1)]
    return R, Y, phi, rac


def lift(f, z, S):
    ans = S.zero()
    power = S.one()
    for c in f.list():
        ans += S(c)*power
        power *= z
    return ans


def modq(q, K, p):
    q = QQ(q)
    den = ZZ(q.denominator()) % p
    assert den != 0
    return K(ZZ(q.numerator()) % p) / K(den)


def polynomial_checks(nmax):
    R, Y, phi, rac = build(nmax+1)
    for j in range(nmax+1):
        for k in range(nmax+1):
            assert phi[k](la(j)) == (u(j, k) if k <= j else 0)
        for n in range(nmax+1):
            value = sum(u(n, k)*u(j, k) for k in range(min(n, j)+1))
            assert rac[n](la(j)) == value == rac[j](la(n))
    for j in range(nmax+1):
        prev = rac[j-1] if j else R.zero()
        assert ((j+1)**3*rac[j+1]
                == (j**3+(j+1)**3+2*(2*j+1)*Y)*rac[j]-j**3*prev)

    S = PolynomialRing(QQ, names=('x', 'y'))
    x, y = S.gens()
    for L in range(1, nmax+1):
        rx = [lift(rac[j], x, S) for j in range(L+1)]
        ry = [lift(rac[j], y, S) for j in range(L+1)]
        Kxy = sum((2*j+1)*rx[j]*ry[j] for j in range(L))
        assert (2*(x-y)*Kxy
                == L**3*(rx[L]*ry[L-1]-rx[L-1]*ry[L]))
        Kxx = sum(QQ(2*j+1)*rac[j]**2 for j in range(L))
        assert (Kxx == QQ(L**3)/QQ(2)
                * (rac[L-1]*rac[L].derivative()
                   - rac[L]*rac[L-1].derivative()))


def matrix_checks(nmax):
    _, _, phi, rac = build(nmax)
    for N in range(nmax+1):
        d = N+1
        U = matrix(QQ, d, d, lambda i, j: QQ(u(i, j)))
        B = matrix(QQ, d, d, lambda i, j: invu(i, j))
        I = identity_matrix(QQ, d)
        assert U*B == I == B*U
        W = diagonal_matrix(QQ, [2*j+1 for j in range(d)])
        H = matrix(QQ, d, d,
                   lambda i, j: QQ((-1)**(i+j))/QQ(i+j+1))
        A = U.transpose()*W*U
        assert B*W.inverse()*B.transpose() == H
        assert A*H == I
        Aclosed = matrix(
            QQ, d, d,
            lambda i, j: QQ((i+j+1)
                * binomial(N+i+1, N-j)
                * binomial(N+j+1, N-i)
                * binomial(i+j, i)**2))
        assert A == Aclosed
        M = U*U.transpose()
        G = M*W*M
        assert G.inverse() == B.transpose()*H*B
        detU = prod(binomial(2*j, j) for j in range(d))
        assert U.det() == detU
        assert G.det() == detU**4 * prod(2*j+1 for j in range(d))
        for i in range(d):
            for j in range(d):
                assert M[i, j] == rac[i](la(j))
        v = vector(QQ, [1]*N + [N+1])
        expected = vector(QQ, [QQ((2*N+1)*(k+1)*u(N, k))/QQ(2*k+1)
                               for k in range(d)])
        assert U.transpose()*v == expected
        assert (sum(v[n]*rac[n] for n in range(d))
                == sum(expected[k]*phi[k] for k in range(d)))


def midpoint_check(p):
    p = ZZ(p)
    N = (p-1)//2
    d = N+1
    K = GF(p)
    U = matrix(K, d, d, lambda i, j: K(u(i, j)))
    assert U.det() != 0
    M = U*U.transpose()
    W = diagonal_matrix(K, [K(2*j+1) for j in range(d)])
    G = M*W*M
    v = vector(K, [1]*N + [N+1])
    assert G*v == 0 and G.rank() == N
    _, _, phi, rac = build(N, K)
    assert (sum(v[n]*rac[n] for n in range(d))
            == K(N+1)*K(u(N, N))*phi[N])
    last = vector(K, [modq(invu(N, k), K, p) for k in range(d)])
    assert last == K((-1)**N)/K(N+1)*v
    assert len({K(la(j)) for j in range(d)}) == d
    detG = (prod(binomial(2*j, j) for j in range(d))**4
            * prod(2*j+1 for j in range(d)))
    assert ZZ(detG).valuation(p) == 1


def obstruction_191():
    _, Y, _, rac = build(3)
    P = 10*Y**3 + 55*Y**2 + 66*Y + 18
    assert 18*rac[3] == P
    assert P.discriminant() == 584460 == 191*3060
    K = GF(191)
    y = K(148)
    q = P.change_ring(K)
    assert q(y) == 0
    assert q.derivative()(y) == 0
    assert q.derivative().derivative()(y) == K(13) != 0
    assert K(la(81)) == y and 191 > 2*3+1
    vals = [rac[j].change_ring(K)(y) for j in range(3)]
    assert vals == [K(1), K(106), K(67)]
    assert sum(K(2*j+1)*vals[j]**2 for j in range(3)) == 0
    assert rac[3](la(3)) == 1445


def main():
    ap = ArgumentParser()
    ap.add_argument('--max-N', type=int, default=9)
    ap.add_argument('--prime-bound', type=int, default=80)
    a = ap.parse_args()
    assert a.max_N >= 3
    polynomial_checks(a.max_N)
    matrix_checks(a.max_N)
    for p in prime_range(3, a.prime_bound):
        midpoint_check(p)
    obstruction_191()
    print('Q8211_EXACT_VERIFY PASS')
    print('Finite loops are regressions only; all formulas are proved in the report.')


if __name__ == '__main__':
    main()
