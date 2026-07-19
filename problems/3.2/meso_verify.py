#!/usr/bin/env python3
"""Independent exact checks for ``meso_result.tex`` and its data report.

The checks are finite symbolic cross-checks; the proofs are in
``meso_result.tex``.  Launching with ordinary Python automatically switches
to ``sage -python``.
"""

from __future__ import annotations

from collections import Counter
from math import isqrt, log
import os
from pathlib import Path
import shutil
import sys


SAGE_CACHE = "/tmp/ramanujan-meso-sage"
os.environ.setdefault("DOT_SAGE", SAGE_CACHE)
Path(os.environ["DOT_SAGE"]).mkdir(parents=True, exist_ok=True)

try:
    from sage.all import (  # type: ignore[import-not-found]
        GF,
        PolynomialRing,
        QQ,
        ZZ,
        factorial,
    )
except ModuleNotFoundError:
    sage = shutil.which("sage")
    if sage is None:
        raise SystemExit("SageMath is required (the `sage` executable was not found)")
    environment = os.environ.copy()
    environment["DOT_SAGE"] = SAGE_CACHE
    os.execvpe(sage, [sage, "-python", __file__, *sys.argv[1:]], environment)

sys.path.insert(0, str(Path(__file__).resolve().parent))
import meso_explore as me  # noqa: E402


HEIGHT = 24
MAX_H = HEIGHT - 2
EXPECTED_DIGEST = "c50c84589bc8426e978b8299e82aa415f4fb047de283083fe2549810f52d672f"


def check_gap_polynomials(n_polynomials: list) -> None:
    ell = me.pell_values(MAX_H)
    x = me.X
    for h in range(1, MAX_H + 1):
        polynomial = n_polynomials[h]
        assert polynomial.degree() == 3 * (h - 1)
        assert polynomial.leading_coefficient() == ell[h]
        reflected = polynomial(-x - h - 1)
        assert reflected == (-1) ** (h - 1) * polynomial
    print("GAP_RECURRENCE_REFLECTION PASS")


def check_centered_normalisation(
    n_polynomials: list, resultants: dict[tuple[int, int], int]
) -> None:
    qy = PolynomialRing(QQ, "Y")
    y = qy.gen()
    centered_polynomials = [qy.zero() for _ in range(MAX_H + 1)]
    for h in range(1, MAX_H + 1):
        degree = 3 * (h - 1)
        centered = qy(
            2**degree * n_polynomials[h]((y - h - 1) / 2)
        )
        centered_polynomials[h] = centered
        assert all(coefficient.denominator() == 1 for coefficient in centered)
        assert centered.leading_coefficient() == me.pell_values(h)[h]
        assert centered(-y) == (-1) ** (h - 1) * centered(y)
    for d, r in me.canonical_pairs(HEIGHT):
        degree_d = 3 * (d - 1)
        degree_r = 3 * (r - 1)
        actual = centered_polynomials[d].resultant(
            centered_polynomials[r](y + d + r)
        )
        assert actual == 2 ** (degree_d * degree_r) * resultants[d, r]
    print("CENTERED_INTEGER_NORMALISATION PASS")


def check_resultants(
    n_polynomials: list, resultants: dict[tuple[int, int], int]
) -> None:
    assert len(resultants) == 121
    assert len(me.ordered_pairs(HEIGHT)) == 231
    assert me.exact_digest(resultants) == EXPECTED_DIGEST
    for d, r in me.canonical_pairs(HEIGHT):
        forward = n_polynomials[d].resultant(me.shifted(n_polynomials[r], d))
        reverse = n_polynomials[r].resultant(me.shifted(n_polynomials[d], r))
        assert forward > 0
        assert abs(forward) == abs(reverse) == resultants[d, r]
    assert len(str(resultants[12, 12])) == 2246
    print("ALL_121_RESULTANTS_AND_SYMMETRY PASS")


def check_complete_small_factorisations(
    resultants: dict[tuple[int, int], int]
) -> None:
    factors: dict[tuple[int, int], tuple[tuple[int, int], ...]] = {
        (2, 2): ((2, 5), (5, 4), (11, 2), (17, 3), (71, 1)),
        (2, 3): ((5, 3), (23, 1), (103, 1), (129763, 1), (685784151628061, 1)),
        (2, 4): (
            (2, 5), (3, 6), (5, 7), (11, 1), (13, 1), (17, 3),
            (41, 1), (653, 1), (52385933, 1), (6024613671641, 1),
        ),
        (3, 3): (
            (3, 17), (5, 5), (7, 8), (11, 3), (89153, 1),
            (110629, 2), (315735829, 2),
        ),
        (2, 5): (
            (3, 6), (5, 8), (43, 1), (113, 1), (397, 1),
            (11119, 1), (30757, 1), (34283, 1), (233591, 1),
            (21869834584034319215521, 1),
        ),
        (3, 4): (
            (2, 14), (5, 5), (13, 1), (73, 1), (131, 1),
            (57131, 1), (75577, 1),
            (15260305013831376290046094059172675684932975108889335505111957, 1),
        ),
        (2, 6): (
            (2, 8), (3, 6), (5, 7), (17, 3), (41, 1), (1091, 1),
            (3919541, 1), (13365198871, 1), (5998798267201, 1),
            (511088437530213575388589, 1),
        ),
        (3, 5): (
            (2, 18), (5, 7), (59, 1), (129763, 1), (5784767, 1),
            (98608507, 1), (521339821, 1), (170262135457, 1),
            (111600466182360778561621, 1),
            (22456093031506928723038485585657593428847, 1),
        ),
        (4, 4): (
            (2, 42), (5, 11), (13, 1), (17, 3), (79, 1), (101, 2),
            (577, 3), (4787, 2), (8431, 2), (56003, 1),
            (271030933, 2), (45032941417, 2), (16802972635975249, 2),
        ),
    }
    for pair, factorisation in factors.items():
        product = ZZ.one()
        for prime, exponent in factorisation:
            assert ZZ(prime).is_prime(proof=True)
            product *= ZZ(prime) ** exponent
        assert product == resultants[pair]
    print("COMPLETE_FACTORISATIONS_D_PLUS_R_LE_8 PASS")


def check_adjacent_and_bezout(n_polynomials: list) -> None:
    apery = me.apery_values(12)
    for h in range(2, 12):
        signed_actual = n_polynomials[h].resultant(n_polynomials[h + 1])
        signed_previous = n_polynomials[h - 1].resultant(n_polynomials[h])
        assert signed_actual == (
            (-1) ** (h - 1)
            * n_polynomials[h](-h) ** 6
            * signed_previous
        )
        actual = abs(signed_actual)
        expected = ZZ.prod(
            (factorial(j) ** 3 * apery[j]) ** 6 for j in range(1, h)
        )
        assert actual == expected

    x = me.X
    for d in range(2, 11):
        cut_edge = ZZ.one()
        for j in range(2, d + 1):
            cut_edge *= x + j
        for e in range(d + 1, 13):
            left = (
                n_polynomials[e - 1](x + 1) * n_polynomials[d]
                - n_polynomials[d - 1](x + 1) * n_polynomials[e]
            )
            right = cut_edge**6 * n_polynomials[e - d](x + d)
            assert left == right

    finite_ring = PolynomialRing(GF(73), "x")
    xp = finite_ring.gen()
    cut_edge_gcd = finite_ring(n_polynomials[3]).gcd(
        finite_ring(n_polynomials[4])
    ).monic()
    assert cut_edge_gcd == xp + 3
    assert me.apery_values(2)[2] == 73
    for h in range(3, MAX_H + 1):
        assert n_polynomials[h](-3) % 73 == 0
    print("ADJACENT_RESULTANT_AND_BEZOUT PASS")


def parity_parts(poly, variable, target_ring):
    even = target_ring.zero()
    odd = target_ring.zero()
    u = target_ring.gen()
    for exponent in range(poly.degree() + 1):
        coefficient = poly[exponent]
        if exponent % 2:
            odd += coefficient * u ** ((exponent - 1) // 2)
        else:
            even += coefficient * u ** (exponent // 2)
    return even, odd


def quadratic_invariant_parts(poly, shift: int, target_ring):
    """Write poly(X)=A(V)+X*B(V), where V=X*(X+shift)."""

    variable = target_ring.gen()
    power_even = target_ring.one()
    power_odd = target_ring.zero()
    invariant_even = target_ring.zero()
    invariant_odd = target_ring.zero()
    for exponent in range(poly.degree() + 1):
        invariant_even += poly[exponent] * power_even
        invariant_odd += poly[exponent] * power_odd
        power_even, power_odd = (
            variable * power_odd,
            power_even - shift * power_odd,
        )
    return invariant_even, invariant_odd


def check_centered_norm_and_diagonal(
    n_polynomials: list, resultants: dict[tuple[int, int], int]
) -> None:
    qt = PolynomialRing(QQ, "T")
    qu = PolynomialRing(QQ, "U")
    t = qt.gen()

    for d in range(2, 13):
        n = 3 * (d - 1)
        f = qt(n_polynomials[d](t - QQ(2 * d + 1) / 2))
        even, odd = parity_parts(f, t, qu)
        invariant_even, invariant_odd = quadratic_invariant_parts(
            n_polynomials[d], 2 * d + 1, qu
        )
        shift = QQ((2 * d + 1) ** 2) / 4
        shifted_invariant_even = invariant_even(qu.gen() - shift)
        shifted_invariant_odd = invariant_odd(qu.gen() - shift)
        assert odd == shifted_invariant_odd
        assert even == (
            shifted_invariant_even
            - QQ(2 * d + 1) / 2 * shifted_invariant_odd
        )
        degree_gap = even.degree() - invariant_even.degree()
        assert degree_gap == 0
        if n % 2:
            assert invariant_odd.leading_coefficient() == n_polynomials[d].leading_coefficient()
            assert invariant_even.leading_coefficient() == (
                QQ(2 * d + 1 - n * d)
                * n_polynomials[d].leading_coefficient()
                / 2
            )
        integral_square_root = (
            invariant_odd.leading_coefficient() ** degree_gap
            * invariant_odd.resultant(invariant_even)
        )
        assert integral_square_root in ZZ
        assert odd.resultant(even) == integral_square_root
        rhs = (
            (-1) ** (d - 1)
            * 2**n
            * n_polynomials[d].leading_coefficient()
            * f(0)
            * odd.resultant(even) ** 2
        )
        direct = n_polynomials[d].resultant(me.shifted(n_polynomials[d], d))
        assert rhs.denominator() == 1
        assert ZZ(rhs) == direct == resultants[d, d]

    for d in range(2, 13):
        n = 3 * (d - 1)
        central_integer = ZZ(2**n * n_polynomials[d](-QQ.one() / 2))
        denominator = (
            n_polynomials[d].leading_coefficient() * abs(central_integer)
        )
        assert resultants[d, d] % denominator == 0
        square_part = ZZ(resultants[d, d] // denominator)
        assert square_part.is_square()

    qy = PolynomialRing(QQ, "Y")
    y = qy.gen()
    for d, r in me.ordered_pairs(HEIGHT):
        centered = qy(n_polynomials[d](y - QQ(d + 1) / 2))
        even, odd = parity_parts(centered, y, qu)
        if d % 2:
            assert odd.is_zero()
            a_poly = even
            center_factor = QQ.one()
        else:
            assert even.is_zero()
            a_poly = odd
            center_factor = n_polynomials[r](QQ(d - 1) / 2)
        g = qy(n_polynomials[r](y + QQ(d - 1) / 2))
        paired = g(y) * g(-y)
        assert all(paired[exponent] == 0 for exponent in range(1, paired.degree() + 1, 2))
        paired_u = qu(
            sum(paired[2 * j] * qu.gen() ** j for j in range(paired.degree() // 2 + 1))
        )
        norm_value = center_factor * a_poly.resultant(paired_u)
        assert abs(norm_value) == resultants[min(d, r), max(d, r)]
    print("CENTERED_NORM_AND_DIAGONAL_SQUARE_LAW PASS")


def odd_part(value: int) -> int:
    value = abs(ZZ(value))
    while value and value % 2 == 0:
        value //= 2
    return value


def center_value(n_polynomials: list, a: int, b: int) -> int:
    value = 2 ** (3 * (b - 1)) * n_polynomials[b](QQ(a - 1) / 2)
    assert value.denominator() == 1
    return ZZ(value)


def check_center_recurrence_lattice(
    n_polynomials: list, resultants: dict[tuple[int, int], int]
) -> None:
    recurrence_cases = 0
    divisibility_cases = 0
    for a in range(2, 21, 2):
        values = [ZZ.zero(), ZZ.one()]
        for b in range(1, MAX_H):
            if b >= 2:
                assert values[b] == center_value(n_polynomials, a, b)
            z = a + 2 * b - 1
            coefficient = 34 * z**3 + 102 * z**2 + 108 * z + 40
            values.append(coefficient * values[b] - z**6 * values[b - 1])
            recurrence_cases += 1

        for b in range(2, MAX_H + 1):
            if a + 2 * b > HEIGHT:
                continue
            central = center_value(n_polynomials, a, b)
            first = me.lookup_resultant(resultants, a, b)
            second = me.lookup_resultant(resultants, b, a + b)
            assert first % abs(central) == 0
            assert second % odd_part(central) == 0
            divisibility_cases += 1
    assert recurrence_cases == 210
    assert divisibility_cases == 55

    for d, r in me.canonical_pairs(HEIGHT):
        if d % 2 and r % 2:
            structural_factor = ZZ.one()
        elif d % 2 == 0 and r % 2:
            structural_factor = abs(center_value(n_polynomials, d, r))
        elif d % 2 and r % 2 == 0:
            structural_factor = abs(center_value(n_polynomials, r, d))
        else:
            numerator = abs(
                center_value(n_polynomials, d, r)
                * center_value(n_polynomials, r, d)
            )
            assert numerator % (2 * (d + r)) == 0
            structural_factor = numerator // (2 * (d + r))
        assert resultants[d, r] % structural_factor == 0
    print("CENTER_RECURRENCE_AND_RESULTANT_LATTICE PASS")


def check_resultant_validity_counterexamples(
    n_polynomials: list, resultants: dict[tuple[int, int], int]
) -> None:
    ring17 = PolynomialRing(GF(17), "x")
    x17 = ring17.gen()
    first17 = ring17(n_polynomials[2])
    second17 = ring17(n_polynomials[2])(x17 + 2)
    assert resultants[2, 2] % 17 == 0
    assert first17 == 5 * (2 * x17 + 3)
    assert second17 == 5 * (2 * x17 + 7)
    assert first17.gcd(second17).degree() == 0

    prime = 110629
    assert ZZ(prime).is_prime(proof=True)
    finite_ring = PolynomialRing(GF(prime), "x")
    xp = finite_ring.gen()
    gcd_poly = finite_ring(n_polynomials[3]).gcd(
        finite_ring(n_polynomials[3])(xp + 3)
    ).monic()
    assert gcd_poly == xp**2 + 7 * xp - 15495
    assert not GF(prime)(62029).is_square()
    assert resultants[3, 3] % prime == 0

    prime = 4787
    assert ZZ(prime).is_prime(proof=True)
    finite_ring = PolynomialRing(GF(prime), "x")
    xp = finite_ring.gen()
    gcd_poly = finite_ring(n_polynomials[4]).gcd(
        finite_ring(n_polynomials[4])(xp + 4)
    ).monic()
    assert gcd_poly == xp**2 + 9 * xp + 2367
    assert not gcd_poly.discriminant().is_square()
    assert resultants[4, 4] % prime == 0

    for prime, rank in ((577, 4), (1153, 6)):
        finite_ring = PolynomialRing(GF(prime), "x")
        for h in range(rank, MAX_H + 1, rank):
            reduced = finite_ring(n_polynomials[h])
            assert reduced.degree() == 3 * (h - 1) - 2
        for d, r in me.canonical_pairs(HEIGHT):
            if d % rank == 0 and r % rank == 0:
                assert resultants[d, r] % (prime**2) == 0
    print("RESULTANT_VALIDITY_COUNTEREXAMPLES PASS")


def check_centered_pell_transfer_failure(n_polynomials: list) -> None:
    prime = 71
    finite_ring = PolynomialRing(GF(prime), "x")
    x = finite_ring.gen()
    n2 = finite_ring(n_polynomials[2])
    x0 = -GF(prime)(5) / 2
    assert n2(x0) == 0
    assert n2(x0 + 2) == 0
    assert x0 + GF(prime)(3) / 2 == -1
    assert x0 + 2 + GF(prime)(3) / 2 == 1
    assert me.pell_values(2)[2] % prime == 34
    print("POINTWISE_CENTERED_PELL_OBSTRUCTION PASS")


def check_prime_statistics(
    n_polynomials: list, resultants: dict[tuple[int, int], int]
) -> None:
    records = me.prime_statistics(HEIGHT, 5000, n_polynomials, resultants)
    mesoscopic = [record for record in records if record.prime >= HEIGHT**2]
    assert len(records) == 669
    assert len(mesoscopic) == 564
    expected = {
        "resultant_support": {0: 506, 1: 3, 2: 38, 3: 2, 4: 10, 6: 3, 15: 1, 16: 1},
        "algebraic_support": {0: 508, 1: 3, 2: 38, 3: 2, 4: 10, 6: 2, 16: 1},
        "fp_support": {0: 509, 1: 2, 2: 39, 3: 1, 4: 10, 6: 2, 16: 1},
        "energy": {0: 509, 1: 2, 2: 39, 4: 11, 6: 2, 16: 1},
    }
    for field, target in expected.items():
        assert dict(Counter(getattr(record, field) for record in mesoscopic)) == target
    by_prime = {record.prime: record for record in mesoscopic}
    assert by_prime[653].resultant_support == 16
    assert by_prime[653].fp_support == 16
    assert by_prime[653].energy == 16
    assert by_prime[797].algebraic_support == 3
    assert by_prime[797].fp_support == 2
    assert by_prime[3109].fp_support == 3
    assert by_prime[3109].energy == 4
    assert (
        by_prime[577].resultant_support,
        by_prime[577].algebraic_support,
        by_prime[577].fp_support,
        by_prime[577].energy,
    ) == (15, 0, 0, 0)
    assert (
        by_prime[1153].resultant_support,
        by_prime[1153].algebraic_support,
        by_prime[1153].fp_support,
        by_prime[1153].energy,
    ) == (6, 0, 0, 0)
    assert (
        by_prime[4787].resultant_support,
        by_prime[4787].algebraic_support,
        by_prime[4787].fp_support,
        by_prime[4787].energy,
    ) == (1, 1, 0, 0)
    assert [
        record.prime
        for record in mesoscopic
        if record.resultant_support and not record.fp_support
    ] == [577, 1153, 4787]

    _, roots = me.polynomial_data_mod_prime(n_polynomials, MAX_H, 653)
    support = []
    for d, r in me.ordered_pairs(HEIGHT):
        common = sorted(x for x in roots[d] if (x + d) % 653 in roots[r])
        support.extend((d, r, x) for x in common)
    assert support == [
        (2, 4, 313), (2, 18, 313), (2, 22, 313), (4, 2, 333),
        (4, 14, 315), (4, 18, 315), (4, 20, 315), (6, 14, 313),
        (6, 18, 313), (14, 4, 319), (14, 6, 319), (18, 2, 319),
        (18, 4, 315), (18, 6, 315), (20, 4, 313), (22, 2, 315),
    ]
    finite_ring = PolynomialRing(GF(653), "x")
    xp = finite_ring.gen()
    for d, r, root in support:
        common_gcd = finite_ring(n_polynomials[d]).gcd(
            finite_ring(n_polynomials[r])(xp + d)
        ).monic()
        assert common_gcd == xp - root

    finite_ring = PolynomialRing(GF(3109), "x")
    xp = finite_ring.gen()
    common_gcd = finite_ring(n_polynomials[8]).gcd(
        finite_ring(n_polynomials[8])(xp + 8)
    ).monic()
    assert common_gcd.degree() == 2
    assert len(common_gcd.roots(multiplicities=False)) == 2

    assert max(record.resultant_support for record in records) == 231
    assert max(record.energy for record in records) == 747
    assert [record.prime for record in records if record.energy == 747] == [5]
    print("PRIME_STATISTICS_TO_5000 PASS")


def valuation(value: int, prime: int) -> int:
    exponent = 0
    while value % prime == 0:
        value //= prime
        exponent += 1
    return exponent


def check_report_factor_data(
    n_polynomials: list, resultants: dict[tuple[int, int], int]
) -> None:
    primes = me.prime_sieve(5000)
    resultant_support = {
        prime
        for prime in primes
        if any(value % prime == 0 for value in resultants.values())
    }
    apery = me.apery_values(MAX_H)
    pell = me.pell_values(MAX_H)
    apery_support_small = {
        prime for prime in primes if any(value % prime == 0 for value in apery[1:])
    }
    pell_support_small = {
        prime for prime in primes if any(value % prime == 0 for value in pell[1:])
    }
    assert len(resultant_support) == 132
    assert len(resultant_support & apery_support_small) == 14
    assert len(resultant_support & pell_support_small) == 37
    assert len(resultant_support - apery_support_small - pell_support_small) == 90

    expected_frequencies = {
        5: (231, 146),
        2: (210, 618),
        3: (203, 260),
        11: (151, 38),
        7: (140, 30),
        17: (120, 8),
        19: (84, 8),
        13: (79, 9),
    }
    ordered_values = [
        me.lookup_resultant(resultants, d, r)
        for d, r in me.ordered_pairs(HEIGHT)
    ]
    for prime, expected in expected_frequencies.items():
        actual = (
            sum(value % prime == 0 for value in ordered_values),
            max(valuation(value, prime) for value in ordered_values),
        )
        assert actual == expected

    def full_prime_support(values: list[int]) -> set[int]:
        return {
            int(prime)
            for value in values
            for prime, _ in ZZ(abs(value)).factor(proof=True)
        }

    apery_support = full_prime_support(apery[1:])
    pell_support = full_prime_support(pell[1:])
    known_support = apery_support | pell_support
    assert (len(apery_support), len(pell_support), len(known_support)) == (49, 69, 109)

    total_log_height = 0.0
    remaining_log_height = 0.0
    remaining_cofactors = []
    for value in resultants.values():
        remaining = ZZ(value)
        for prime in known_support:
            while remaining % prime == 0:
                remaining //= prime
        remaining_cofactors.append(remaining)
        total_log_height += log(int(value))
        remaining_log_height += log(int(remaining))
    explained_fraction = 1.0 - remaining_log_height / total_log_height
    assert abs(explained_fraction - 0.174972326547) < 1e-12
    assert all(value > 1 for value in remaining_cofactors)
    assert sum(value.is_prime(proof=True) for value in remaining_cofactors) == 1
    assert sum(not value.is_prime(proof=True) for value in remaining_cofactors) == 120
    assert max(len(str(value)) for value in remaining_cofactors) == 1796

    factors, cofactor = me.trial_factor(resultants[3, 6], primes)
    assert factors == (
        (2, 18), (3, 23), (5, 10), (7, 3), (11, 3), (23, 2), (53, 1)
    )
    for prime in (21187, 3605629):
        assert ZZ(prime).is_prime(proof=True)
        assert cofactor % prime == 0
        cofactor //= prime
    hard_composite = ZZ(
        "82331328057159957094082641924652164496892111499574525894473111366919585965389051882753610906780816656199"
    )
    assert cofactor == hard_composite
    assert len(str(hard_composite)) == 104
    assert not hard_composite.is_prime(proof=True)

    center_prime = ZZ(30681778082168266499711406058345639)
    assert center_prime.is_prime(proof=True)
    central_value = center_value(n_polynomials, 4, 10)
    assert central_value % center_prime == 0
    assert resultants[4, 10] % center_prime == 0
    assert resultants[10, 14] % center_prime == 0
    print("REPORT_FACTOR_AND_REPETITION_DATA PASS")


def check_infinity_rank_example() -> None:
    prime = ZZ(665857)
    assert prime.is_prime(proof=True)
    previous, current = 0, 1
    first_zero = None
    for index in range(1, 9):
        if current % prime == 0:
            first_zero = index
            break
        previous, current = current, (34 * current - previous) % prime
    assert first_zero == 8
    height = isqrt(int(prime))
    assert height == 816
    multiples = height // first_zero
    assert multiples * (multiples - 1) // 2 == 5151
    print("PELL_INFINITY_RANK_EXAMPLE PASS")


def main() -> None:
    n_polynomials = me.build_gap_polynomials(MAX_H)
    resultants = me.exact_resultants(n_polynomials, HEIGHT)
    check_gap_polynomials(n_polynomials)
    check_centered_normalisation(n_polynomials, resultants)
    check_resultants(n_polynomials, resultants)
    check_complete_small_factorisations(resultants)
    check_adjacent_and_bezout(n_polynomials)
    check_centered_norm_and_diagonal(n_polynomials, resultants)
    check_center_recurrence_lattice(n_polynomials, resultants)
    check_resultant_validity_counterexamples(n_polynomials, resultants)
    check_centered_pell_transfer_failure(n_polynomials)
    check_prime_statistics(n_polynomials, resultants)
    check_report_factor_data(n_polynomials, resultants)
    check_infinity_rank_example()
    print("MESO_VERIFY PASS")


if __name__ == "__main__":
    main()
