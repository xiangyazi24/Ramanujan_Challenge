#!/usr/bin/env sage
"""Exact verifier for Q8211.

This file uses QQ/ZZ/GF(p) arithmetic only.  Its finite loops are regression
checks for the general identities proved in the accompanying report; they are
not used as finite extrapolation.

Run from the repository root:

    sage drops/Q8211-698c0700-verify.sage

Optional larger exact regressions:

    sage drops/Q8211-698c0700-verify.sage --max-N 12 --prime-bound 120
"""

from argparse import ArgumentParser
from sage.all import (
    GF,
    QQ,
    ZZ,
    PolynomialRing,
    binomial,
    diagonal_matrix,
    factorial,
    identity_matrix,
    matrix,
    prime_range,
    prod,
    vector,
)


def lam(j):
    return ZZ(j) * ZZ(j + 1)


def U_entry(n, k):
    if k < 0 or k > n:
        return ZZ(0)
    return ZZ(binomial(n, k) * binomial(n + k, k))


def B_entry(n, k):
    """The reduced rational entry of U^{-1}."""
    if k < 0 or k > n:
        return QQ(0)
    return (
        QQ((-1) ** (n - k) * (2 * k + 1) * factorial(n) ** 2)
        / QQ(factorial(n - k) * factorial(n + k + 1))
    )


def build_polynomials(maximum, base=QQ):
    ring = PolynomialRing(base, "Y")
    Y = ring.gen()
    phi = [ring.one()]
    for k in range(1, maximum + 1):
        phi.append(phi[-1] * (Y - base((k - 1) * k)) / base(k * k))
    racah = []
    for n in range(maximum + 1):
        racah.append(
            sum(base(U_entry(n, k)) * phi[k] for k in range(n + 1))
        )
    return ring, Y, phi, racah


def lift_polynomial(polynomial, value, target_ring):
    answer = target_ring.zero()
    power = target_ring.one()
    for coefficient in polynomial.list():
        answer += target_ring(coefficient) * power
        power *= value
    return answer


def qreduce(value, field, prime):
    value = QQ(value)
    numerator = ZZ(value.numerator()) % prime
    denominator = ZZ(value.denominator()) % prime
    assert denominator != 0, (value, prime)
    return field(numerator) / field(denominator)


def check_newton_and_recurrence(maximum):
    ring, Y, phi, racah = build_polynomials(maximum + 1)

    # Newton evaluation and self-duality.
    for j in range(maximum + 1):
        for k in range(maximum + 1):
            expected = U_entry(j, k) if k <= j else 0
            assert phi[k](lam(j)) == expected, (j, k)
        for n in range(maximum + 1):
            expected = sum(
                U_entry(n, k) * U_entry(j, k)
                for k in range(min(n, j) + 1)
            )
            assert racah[n](lam(j)) == expected
            assert racah[n](lam(j)) == racah[j](lam(n))

    # Denominator-free three-term recurrence in the polynomial variable.
    for j in range(maximum + 1):
        previous = racah[j - 1] if j > 0 else ring.zero()
        lhs = (j + 1) ** 3 * racah[j + 1]
        rhs = (
            (j ** 3 + (j + 1) ** 3 + 2 * (2 * j + 1) * Y) * racah[j]
            - j ** 3 * previous
        )
        assert lhs == rhs, j

    # The displayed finite Green equation at physical spectral nodes.
    for n in range(maximum + 1):
        for j in range(maximum):
            rjm = racah[n](lam(j - 1)) if j > 0 else QQ(0)
            rj = racah[n](lam(j))
            rjp = racah[n](lam(j + 1))
            lhs = (j + 1) ** 3 * (rjp - rj) + j ** 3 * (rjm - rj)
            rhs = 2 * lam(n) * (2 * j + 1) * rj
            assert lhs == rhs


def check_christoffel_darboux(maximum):
    _, _, _, racah = build_polynomials(maximum)
    bivariate = PolynomialRing(QQ, names=("x", "y"))
    x, y = bivariate.gens()

    for L in range(1, maximum + 1):
        Rx = [lift_polynomial(racah[j], x, bivariate) for j in range(L + 1)]
        Ry = [lift_polynomial(racah[j], y, bivariate) for j in range(L + 1)]
        kernel = sum((2 * j + 1) * Rx[j] * Ry[j] for j in range(L))
        lhs = 2 * (x - y) * kernel
        rhs = L ** 3 * (Rx[L] * Ry[L - 1] - Rx[L - 1] * Ry[L])
        assert lhs == rhs, L

        diagonal_kernel = sum(
            QQ(2 * j + 1) * racah[j] ** 2 for j in range(L)
        )
        diagonal_rhs = QQ(L ** 3) / 2 * (
            racah[L - 1] * racah[L].derivative()
            - racah[L] * racah[L - 1].derivative()
        )
        assert diagonal_kernel == diagonal_rhs, L


def check_matrix_formulas(maximum):
    _, _, phi, racah = build_polynomials(maximum)

    for N in range(maximum + 1):
        size = N + 1
        U = matrix(QQ, size, size, lambda i, j: QQ(U_entry(i, j)))
        B = matrix(QQ, size, size, lambda i, j: B_entry(i, j))
        identity = identity_matrix(QQ, size)
        assert U * B == identity
        assert B * U == identity

        W = diagonal_matrix(QQ, [2 * j + 1 for j in range(size)])
        signed_hilbert = matrix(
            QQ,
            size,
            size,
            lambda i, j: QQ((-1) ** (i + j), i + j + 1),
        )
        coefficient_gram = U.transpose() * W * U

        # Shifted-Legendre moment identity and explicit inverse Hilbert matrix.
        assert B * W.inverse() * B.transpose() == signed_hilbert
        assert coefficient_gram * signed_hilbert == identity
        closed_coefficient_gram = matrix(
            QQ,
            size,
            size,
            lambda i, j: QQ(
                (i + j + 1)
                * binomial(N + i + 1, N - j)
                * binomial(N + j + 1, N - i)
                * binomial(i + j, i) ** 2
            ),
        )
        assert coefficient_gram == closed_coefficient_gram

        M = U * U.transpose()
        physical_gram = M * W * M
        assert physical_gram.inverse() == B.transpose() * signed_hilbert * B

        det_u = prod(binomial(2 * j, j) for j in range(size))
        assert U.det() == det_u
        assert M.det() == det_u ** 2
        assert physical_gram.det() == (
            det_u ** 4 * prod(2 * j + 1 for j in range(size))
        )

        # M is exactly the physical evaluation matrix and G is the CD kernel.
        for n in range(size):
            for j in range(size):
                assert M[n, j] == racah[n](lam(j))
            for m in range(size):
                direct = sum(
                    (2 * j + 1)
                    * racah[n](lam(j))
                    * racah[m](lam(j))
                    for j in range(size)
                )
                assert physical_gram[n, m] == direct

        # Exact characteristic-zero precursor of the midpoint kernel vector.
        p_symbol = 2 * N + 1
        v = vector(QQ, [1] * N + [N + 1])
        coefficients = U.transpose() * v
        expected_coefficients = vector(
            QQ,
            [
                QQ(p_symbol * (k + 1) * U_entry(N, k), 2 * k + 1)
                for k in range(size)
            ],
        )
        assert coefficients == expected_coefficients
        precursor = sum(v[n] * racah[n] for n in range(size))
        precursor_newton = sum(
            expected_coefficients[k] * phi[k] for k in range(size)
        )
        assert precursor == precursor_newton


def check_midpoint_prime(prime):
    prime = ZZ(prime)
    assert prime.is_prime() and prime % 2 == 1
    N = (prime - 1) // 2
    size = N + 1
    field = GF(prime)

    Uq = matrix(field, size, size, lambda i, j: field(U_entry(i, j)))
    assert Uq.det() != 0
    Mq = Uq * Uq.transpose()
    Wq = diagonal_matrix(field, [field(2 * j + 1) for j in range(size)])
    Gq = Mq * Wq * Mq
    vq = vector(field, [1] * N + [N + 1])
    assert vq != 0
    assert Gq * vq == 0
    assert Gq.rank() == N

    # The radical polynomial is the top Newton factor.
    _, _, phi, racah = build_polynomials(N, field)
    radical_polynomial = sum(vq[n] * racah[n] for n in range(size))
    expected = field(N + 1) * field(U_entry(N, N)) * phi[N]
    assert radical_polynomial == expected

    # The last row of U^{-1} is the same kernel vector, with all QQ
    # fractions reduced only after exact cancellation.
    row = vector(
        field,
        [qreduce(B_entry(N, k), field, prime) for k in range(size)],
    )
    scalar = field((-1) ** N) / field(N + 1)
    assert row == scalar * vq

    # Distinct nodes and the unique p-factor in the characteristic-zero Gram.
    assert len({field(lam(j)) for j in range(size)}) == size
    Uz = matrix(ZZ, size, size, lambda i, j: U_entry(i, j))
    Mz = Uz * Uz.transpose()
    Wz = diagonal_matrix(ZZ, [2 * j + 1 for j in range(size)])
    Gz = Mz * Wz * Mz
    assert ZZ(Gz.det()).valuation(prime) == 1


def check_exact_191_obstruction():
    ring, Y, _, racah = build_polynomials(3)
    integral_r3 = 10 * Y ** 3 + 55 * Y ** 2 + 66 * Y + 18
    assert 18 * racah[3] == integral_r3
    assert integral_r3.discriminant() == 584460 == 191 * 3060

    prime = ZZ(191)
    field = GF(prime)
    y = field(148)
    polynomial = integral_r3.change_ring(field)
    assert polynomial(y) == 0
    assert polynomial.derivative()(y) == 0
    assert polynomial.derivative().derivative()(y) == field(13) != 0
    assert field(lam(81)) == y
    assert prime > 2 * 3 + 1

    values = [racah[j].change_ring(field)(y) for j in range(3)]
    assert values == [field(1), field(106), field(67)]
    assert sum(field(2 * j + 1) * values[j] ** 2 for j in range(3)) == 0
    assert racah[3](lam(3)) == 1445


def main():
    parser = ArgumentParser()
    parser.add_argument("--max-N", type=int, default=9)
    parser.add_argument("--prime-bound", type=int, default=80)
    args = parser.parse_args()
    assert args.max_N >= 3

    check_newton_and_recurrence(args.max_N)
    check_christoffel_darboux(args.max_N)
    check_matrix_formulas(args.max_N)
    for prime in prime_range(3, args.prime_bound):
        check_midpoint_prime(prime)
    check_exact_191_obstruction()

    print("Q8211_EXACT_VERIFY PASS")
    print("Checked: Newton self-duality, recurrence, finite Green/CD,")
    print("         shifted-Hilbert inverses, exact midpoint radical,")
    print("         and the p=191 double physical-grid root.")
    print("Finite loops are regressions only; the report proves the formulas generally.")


if __name__ == "__main__":
    main()
