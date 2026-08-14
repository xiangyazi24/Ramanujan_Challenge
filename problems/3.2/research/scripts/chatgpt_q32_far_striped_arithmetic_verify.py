#!/usr/bin/env python3
"""Exact standard-library verifier for Q8372.

This file audits the algebra used in

    problems/3.2/ORACLE_COMM/chatgpt_q32_far_striped_arithmetic.md

It deliberately does not claim a new Apéry asymptotic.  It verifies:

1. the actual mixed-characteristic far edge (m,q,ell)=(321,179,193);
2. shell-prime partner uniqueness in the Q8345 physical range;
3. the exact layerwise labelled-radical identity;
4. freedom of the two uncontrolled 8 x 8 cross-layer blocks;
5. CRT saturation for sixteen distinct characteristics;
6. invertibility of the Boolean/Walsh transform on the 4-cube;
7. failure of deleting the physical quotient after an affine Boolean identity;
8. the support factorization of a natural mixed-characteristic product norm.

Only Python's standard library is used.  No probabilistic or floating-point
step enters any assertion.
"""

from fractions import Fraction
from itertools import product
from math import gcd, isqrt, prod


def P(n: int) -> int:
    return 34 * n**3 + 51 * n**2 + 27 * n + 5


def fold_signed(s: int) -> int:
    return s if s >= 0 else -s - 1


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    for d in range(3, isqrt(n) + 1, 2):
        if n % d == 0:
            return False
    return True


def primes_between(lo: int, hi: int) -> list[int]:
    return [p for p in range(lo + 1, hi + 1) if is_prime(p)]


def first_primes_after(lo: int, count: int) -> list[int]:
    out: list[int] = []
    n = lo + 1
    while len(out) < count:
        if is_prime(n):
            out.append(n)
        n += 1
    return out


def prime_support(n: int) -> set[int]:
    n = abs(n)
    out: set[int] = set()
    d = 2
    while d * d <= n:
        if n % d == 0:
            out.add(d)
            while n % d == 0:
                n //= d
        d = 3 if d == 2 else d + 2
    if n > 1:
        out.add(n)
    return out


def largest_prime_factor(n: int) -> int:
    support = prime_support(n)
    assert support
    return max(support)


def apery_exact(nmax: int) -> list[int]:
    b = [0] * (nmax + 1)
    b[0] = 1
    if nmax >= 1:
        b[1] = 5
    for n in range(1, nmax):
        numerator = P(n) * b[n] - n**3 * b[n - 1]
        denominator = (n + 1) ** 3
        assert numerator % denominator == 0
        b[n + 1] = numerator // denominator
    return b


def apery_mod(p: int, nmax: int) -> list[int]:
    assert is_prime(p) and nmax < p
    b = [0] * (nmax + 1)
    b[0] = 1 % p
    if nmax >= 1:
        b[1] = 5 % p
    for n in range(1, nmax):
        denominator = (n + 1) ** 3 % p
        assert denominator != 0
        numerator = (P(n) * b[n] - n**3 * b[n - 1]) % p
        b[n + 1] = numerator * pow(denominator, -1, p) % p
    return b


def shell_radical_gcd(X: int, physical: int, apery_value: int) -> int:
    common = gcd(abs(physical), abs(apery_value))
    return prod(p for p in primes_between(X, 2 * X) if common % p == 0)


def lcm(a: int, b: int) -> int:
    return abs(a // gcd(a, b) * b)


def egcd(a: int, b: int) -> tuple[int, int, int]:
    if b == 0:
        return abs(a), 1 if a >= 0 else -1, 0
    g, x1, y1 = egcd(b, a % b)
    return g, y1, x1 - (a // b) * y1


def identity(n: int) -> list[list[int]]:
    return [[int(i == j) for j in range(n)] for i in range(n)]


def zero_matrix(rows: int, cols: int) -> list[list[int]]:
    return [[0 for _ in range(cols)] for _ in range(rows)]


def mat_sub(A: list[list[int]], B: list[list[int]], p: int) -> list[list[int]]:
    return [[(A[i][j] - B[i][j]) % p for j in range(len(A[0]))]
            for i in range(len(A))]


def matmul(A: list[list[int]], B: list[list[int]], p: int) -> list[list[int]]:
    assert len(A[0]) == len(B)
    return [[sum(A[i][k] * B[k][j] for k in range(len(B))) % p
             for j in range(len(B[0]))]
            for i in range(len(A))]


def block_matrix(A: list[list[int]], B: list[list[int]],
                 C: list[list[int]], D: list[list[int]]) -> list[list[int]]:
    assert len(A) == len(B) and len(C) == len(D)
    return [A[i] + B[i] for i in range(len(A))] + [
        C[i] + D[i] for i in range(len(C))
    ]


def rank_mod(A: list[list[int]], p: int) -> int:
    M = [[x % p for x in row] for row in A]
    rows = len(M)
    cols = len(M[0]) if rows else 0
    rank = 0
    for col in range(cols):
        pivot = next((r for r in range(rank, rows) if M[r][col] % p), None)
        if pivot is None:
            continue
        M[rank], M[pivot] = M[pivot], M[rank]
        inv = pow(M[rank][col], -1, p)
        M[rank] = [(x * inv) % p for x in M[rank]]
        for r in range(rows):
            if r == rank:
                continue
            factor = M[r][col] % p
            if factor:
                M[r] = [
                    (M[r][c] - factor * M[rank][c]) % p
                    for c in range(cols)
                ]
        rank += 1
        if rank == rows:
            break
    return rank


def det_mod(A: list[list[int]], p: int) -> int:
    n = len(A)
    assert all(len(row) == n for row in A)
    M = [[x % p for x in row] for row in A]
    determinant = 1
    for col in range(n):
        pivot = next((r for r in range(col, n) if M[r][col] % p), None)
        if pivot is None:
            return 0
        if pivot != col:
            M[col], M[pivot] = M[pivot], M[col]
            determinant = -determinant
        pivot_value = M[col][col] % p
        determinant = determinant * pivot_value % p
        inv = pow(pivot_value, -1, p)
        for r in range(col + 1, n):
            factor = M[r][col] * inv % p
            if factor:
                for c in range(col, n):
                    M[r][c] = (M[r][c] - factor * M[col][c]) % p
    return determinant % p


def det_bareiss(A: list[list[int]]) -> int:
    n = len(A)
    assert all(len(row) == n for row in A)
    if n == 0:
        return 1
    M = [row[:] for row in A]
    sign = 1
    previous = 1
    for k in range(n - 1):
        if M[k][k] == 0:
            pivot = next((r for r in range(k + 1, n) if M[r][k] != 0), None)
            if pivot is None:
                return 0
            M[k], M[pivot] = M[pivot], M[k]
            sign = -sign
        pivot_value = M[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = M[i][j] * pivot_value - M[i][k] * M[k][j]
                assert numerator % previous == 0
                M[i][j] = numerator // previous
        previous = pivot_value
        for i in range(k + 1, n):
            M[i][k] = 0
        for j in range(k + 1, n):
            M[k][j] = 0
    return sign * M[n - 1][n - 1]


def vertices(dimension: int) -> list[tuple[int, ...]]:
    return list(product((0, 1), repeat=dimension))


def walsh_matrix(dimension: int) -> list[list[int]]:
    V = vertices(dimension)
    return [[(-1) ** sum(a * v for a, v in zip(A, x)) for x in V]
            for A in V]


def alternating_sum(values: list[int], dimension: int) -> int:
    V = vertices(dimension)
    assert len(values) == len(V)
    return sum(((-1) ** sum(v)) * value for v, value in zip(V, values))


def actual_mixed_edge_audit() -> None:
    # Q8345's exact physical edge.  Both signed rows are reflected (RR).
    X = 128
    m = 321
    g = 28
    lower_prime = 179
    upper_prime = 193
    lower_signed = -37
    upper_signed = -65

    assert upper_signed == lower_signed - g
    assert fold_signed(lower_signed) == 36
    assert fold_signed(upper_signed) == 64
    assert m - lower_signed == 2 * lower_prime
    assert m - upper_signed == 2 * upper_prime

    b179 = apery_mod(lower_prime, 64)
    b193 = apery_mod(upper_prime, 64)
    assert b179[36] == 0
    assert b193[64] == 0
    # The far edge is mixed characteristic; it is not a same-179 return.
    assert b179[64] != 0

    b = apery_exact(64)
    assert shell_radical_gcd(X, m - lower_signed, b[36]) == lower_prime
    assert shell_radical_gcd(X, m - upper_signed, b[64]) == upper_prime
    assert shell_radical_gcd(X, m - lower_signed, b[64]) == 1

    lower_order = (lower_prime - 1) // gcd(lower_prime - 1, 36)
    upper_order = (upper_prime - 1) // gcd(upper_prime - 1, 64)
    assert (lower_order, upper_order) == (89, 3)

    print("ACTUAL_MIXED_EDGE", m, lower_prime, upper_prime,
          fold_signed(lower_signed), fold_signed(upper_signed),
          lower_order, upper_order)


def shell_partner_uniqueness_audit() -> None:
    for X in (20, 40, 80, 128):
        shell = primes_between(X, 2 * X)
        for i, p in enumerate(shell):
            for q in shell[i:]:
                # This is the Q8345 product argument for physical integers
                # strictly below X^2+X.
                assert p * q > X * X + X
    print("SHELL_PARTNER_UNIQUENESS PASS")


def layerwise_radical_audit() -> None:
    # An exact algebra model for the theorem: each physical integer has one
    # shell-prime candidate and the upper-layer directed cross values are units.
    X = 100
    q = primes_between(X, 2 * X)[:8]
    assert len(q) == 8
    physical = [(i + 2) * q[i] for i in range(8)]

    # Proxy upper Apéry values: q_i divides only its own value among this layer.
    # Arbitrary small and opposite-layer factors are harmless here.
    upper_values = [q[i] * (2 * i + 3) for i in range(8)]
    R = shell_radical_gcd(X, prod(physical), prod(upper_values))
    assert R == prod(q)

    # Delete one corresponding target.  The radical deletes exactly its label.
    missing = 3
    upper_missing = upper_values[:]
    upper_missing[missing] //= q[missing]
    R_missing = shell_radical_gcd(X, prod(physical), prod(upper_missing))
    assert R_missing == prod(q[i] for i in range(8) if i != missing)

    print("LAYERWISE_LABELLED_RADICAL PASS", R.bit_length())


def cross_block_freedom_audit() -> None:
    # The clean within-layer zero pattern gives the two identity blocks.  The
    # 128 directed cross-layer values are two arbitrary 8 x 8 blocks X,Y.
    p = 101
    n = 8
    I = identity(n)
    O = zero_matrix(n, n)

    complete_unit = block_matrix(I, O, O, I)
    assert rank_mod(complete_unit, p) == 16
    assert det_mod(complete_unit, p) == 1

    # Even with Boolean cross blocks, every corank from 0 through 8 occurs.
    # Take X=I and Y=I-A with A diagonal of rank r.  Row elimination gives
    # rank [[I,I],[Y,I]] = 8 + rank(A), and determinant = det(A).
    for r in range(9):
        A = [[int(i == j and i < r) for j in range(n)] for i in range(n)]
        Y = mat_sub(I, A, p)
        Z = block_matrix(I, I, Y, I)
        assert rank_mod(Z, p) == n + r
        expected_det = 1 if r == n else 0
        assert det_mod(Z, p) == expected_det
        I_minus_YX = mat_sub(I, matmul(Y, I, p), p)
        assert det_mod(Z, p) == det_mod(I_minus_YX, p)

    print("CROSS_BLOCK_FREEDOM coranks=0..8 PASS")


def crt_saturation_audit() -> None:
    labels = first_primes_after(200, 16)
    modulus = 1
    for p in labels:
        modulus = lcm(modulus, p)
    assert modulus == prod(labels)
    for i, p in enumerate(labels):
        for q in labels[i + 1:]:
            g, a, b = egcd(p, q)
            assert g == 1 and a * p + b * q == 1
    print("CRT_SATURATION", len(labels), modulus.bit_length(), "bits")


def walsh_audit() -> None:
    H3 = walsh_matrix(3)
    H4 = walsh_matrix(4)
    det3 = det_bareiss(H3)
    det4 = det_bareiss(H4)
    assert abs(det3) == 8**4
    assert abs(det4) == 16**8
    print("WALSH_INVERTIBLE", det3, det4)


def post_radical_affine_audit() -> None:
    # A physical 4-cube is affine, so its fourth Boolean difference is zero.
    # Selecting one prime factor at each vertex destroys that identity; the
    # varying physical quotients are exactly what restores it.
    V = vertices(4)
    gaps = (3, 7, 11, 19)
    witness = None
    for base in range(200, 5000):
        values = [base - sum(e * d for e, d in zip(v, gaps)) for v in V]
        if min(values) <= 1:
            continue
        selected = [largest_prime_factor(value) for value in values]
        if alternating_sum(selected, 4) != 0:
            witness = (base, values, selected)
            break
    assert witness is not None
    base, values, selected = witness
    quotients = [value // prime for value, prime in zip(values, selected)]
    assert alternating_sum(values, 4) == 0
    assert alternating_sum(
        [q * prime for q, prime in zip(quotients, selected)], 4
    ) == 0
    assert alternating_sum(selected, 4) != 0
    print("POST_RADICAL_BOOLEAN_FAILURE", base,
          alternating_sum(selected, 4))


def mixed_norm_support_audit() -> None:
    # For a product of local elements in a common compositum K, norm
    # transitivity gives prod Norm_i(alpha_i)^[K:K_i].  Prime support is only
    # the union of the local supports; no cross-characteristic cancellation
    # occurs in this natural norm.
    local_degrees = [2, 3, 5]
    total_degree = 30
    local_norms = [2 * 3, 5, 7 * 11]
    tensor_norm = prod(
        norm ** (total_degree // degree)
        for norm, degree in zip(local_norms, local_degrees)
    )
    expected_support = set().union(*(prime_support(n) for n in local_norms))
    assert prime_support(tensor_norm) == expected_support
    print("MIXED_NORM_SUPPORT", sorted(expected_support))


def main() -> None:
    actual_mixed_edge_audit()
    shell_partner_uniqueness_audit()
    layerwise_radical_audit()
    cross_block_freedom_audit()
    crt_saturation_audit()
    walsh_audit()
    post_radical_affine_audit()
    mixed_norm_support_audit()
    print("CHATGPT_Q32_FAR_STRIPED_ARITHMETIC PASS")


if __name__ == "__main__":
    main()
