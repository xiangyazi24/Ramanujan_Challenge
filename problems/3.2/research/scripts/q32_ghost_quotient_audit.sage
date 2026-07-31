#!/usr/bin/env sage
"""Exact Route-B audit with independent common-annihilator ghosts.

The doubled-period construction supplies one shell coordinate which is
projectively aliased to the original shell at every candidate prime.  If L is
the least common left multiple of the certified J-operator and the Apéry
operator, then

    L(J) = L(b) = 0.

Consequently each shifted shell

    X_k(d) = sum_j l_j(n+k)
              (C_{n+k+j}(d-1) + C_{n+k+j}(d+1))

vanishes at every safe candidate node d=p-1 modulo p.  Every Newton carrier
of X_k over the core is therefore divisible by the whole candidate
primorial.  After dividing that universal factor, its cross-minor with the
original shell is still divisible by each target prime.

This script tests whether several such quotient ghosts give genuinely
stronger evaluated exterior content than the single doubled-period lift.
All arithmetic is exact.
"""

from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
load(str(HERE / "q32_doubled_period_gauge_audit.sage"))


CASES = (
    (200, 128, 63),
    (272, 180, 63),
    (300, 180, 57),
    (321, 168, 53),
)
GHOST_SHIFTS = tuple(range(5))


def primitive_integral_coefficients(operator):
    coefficients = [K(operator[i]) for i in range(operator.order() + 1)]
    polynomial_denominator = lcm(value.denominator() for value in coefficients)
    polynomials = [R(polynomial_denominator * value) for value in coefficients]
    rational_denominator = lcm(
        coefficient.denominator()
        for polynomial in polynomials
        for coefficient in polynomial.list()
    )
    polynomials = [R(rational_denominator * value) for value in polynomials]
    content = gcd(
        [
            ZZ(coefficient)
            for polynomial in polynomials
            for coefficient in polynomial.list()
        ]
    )
    polynomials = [R(value / content) for value in polynomials]
    assert gcd(polynomials) == 1
    return tuple(polynomials)


common_operator = P.lclm(apery_operator)
assert common_operator.order() == 5
ghost_coefficients = primitive_integral_coefficients(common_operator)
assert all(coefficient.degree() == 35 for coefficient in ghost_coefficients)

# Exact distinguished-solution audit of the common annihilator.
for index in range(2, 20):
    for values in (j_values, b_values):
        assert sum(
            ghost_coefficients[shift](n=index) * values[index + shift]
            for shift in range(len(ghost_coefficients))
        ) == 0


def endpoint(values, d, length):
    return newton_carrier(values, d, length)


def margin_normalizer(D, N, margin):
    d = D - margin + 1
    length = N + margin - 2
    return gcd_many(
        binomial(d + shift + length + 1, length)
        for shift in range(margin - 1)
    )


def target_carrier(Y, ghost, D, N, margin, candidate_product):
    d = D - margin + 1
    length = N + margin - 2
    y_left = endpoint(Y, d, length)
    y_right = endpoint(Y, D, length)
    x_left = endpoint(ghost, d, length)
    x_right = endpoint(ghost, D, length)
    assert x_left % candidate_product == 0
    assert x_right % candidate_product == 0
    normalizer = margin_normalizer(D, N, margin)
    assert normalizer % candidate_product == 0
    reduced_normalizer = normalizer // candidate_product
    numerator = (
        y_left * (x_right // candidate_product)
        - y_right * (x_left // candidate_product)
    )
    assert numerator % reduced_normalizer == 0
    return (
        numerator // reduced_normalizer,
        (y_left, y_right),
        (x_left // candidate_product, x_right // candidate_product),
        reduced_normalizer,
    )


def factor_if_small(value):
    value = abs(ZZ(value))
    return factor(value) if value.nbits() <= 300 else None


def audit_case(index, D, N):
    moment = index - 1
    maximum_margin = min(
        D - moment // 2,
        moment - D - N + 2,
    )
    margins = (maximum_margin - 1, maximum_margin)
    assert margins[0] >= 2
    lower = D - maximum_margin + 1
    upper = D + N + maximum_margin - 2
    nodes = list(range(lower, upper + 1))
    neighbours = list(range(lower - 1, upper + 2))

    candidates = tuple(prime_range(D + 1, D + N + 1))
    candidate_product = prod(candidates)
    targets = tuple(
        prime
        for prime in candidates
        if ZZ(apery(index - prime)) % prime == 0
    )
    target_product = prod(targets)

    Y = shell_batch(moment, nodes)
    ghosts = {}
    for ghost_shift in GHOST_SHIFTS:
        shell_rows = {
            operator_shift: shell_batch(
                index + ghost_shift + operator_shift,
                neighbours,
            )
            for operator_shift in range(len(ghost_coefficients))
        }
        ghost = {
            node: sum(
                ZZ(
                    ghost_coefficients[operator_shift](
                        n=index + ghost_shift
                    )
                )
                * (
                    shell_rows[operator_shift][node - 1]
                    + shell_rows[operator_shift][node + 1]
                )
                for operator_shift in range(len(ghost_coefficients))
            )
            for node in nodes
        }

        # Check the nodewise vanishing before using CRT on Newton rows.
        assert all(ghost[prime - 1] % prime == 0 for prime in candidates)
        ghosts[ghost_shift] = ghost

    data = {
        (ghost_shift, margin): target_carrier(
            Y,
            ghosts[ghost_shift],
            D,
            N,
            margin,
            candidate_product,
        )
        for ghost_shift in GHOST_SHIFTS
        for margin in margins
    }
    carriers = {key: value[0] for key, value in data.items()}
    assert all(value % target_product == 0 for value in carriers.values())

    smith_records = []
    for margin in margins:
        A = data[0, margin][1]
        reduced_normalizer = data[0, margin][3]
        assert all(
            data[ghost_shift, margin][1] == A
            and data[ghost_shift, margin][3] == reduced_normalizer
            for ghost_shift in GHOST_SHIFTS
        )
        quotient_columns = [
            data[ghost_shift, margin][2]
            for ghost_shift in GHOST_SHIFTS
        ]
        lattice_divisor = gcd_many(
            left[0] * right[1] - left[1] * right[0]
            for left_index, left in enumerate(quotient_columns)
            for right in quotient_columns[left_index + 1 :]
        )
        raw_wedge_gcd = gcd_many(
            A[0] * column[1] - A[1] * column[0]
            for column in quotient_columns
        )
        carrier_gcd = gcd_many(
            carriers[ghost_shift, margin]
            for ghost_shift in GHOST_SHIFTS
        )
        assert raw_wedge_gcd == reduced_normalizer * carrier_gcd
        scalar_content = gcd_many(A)
        assert (scalar_content * lattice_divisor) % raw_wedge_gcd == 0
        smith_records.append(
            (
                margin,
                abs(scalar_content).nbits(),
                abs(lattice_divisor).nbits(),
                abs(carrier_gcd).nbits(),
                factor_if_small(lattice_divisor),
            )
        )

    running = ZZ(0)
    history = []
    for ghost_shift in GHOST_SHIFTS:
        for margin in margins:
            running = gcd(running, carriers[ghost_shift, margin])
        history.append(
            (
                ghost_shift,
                abs(running).nbits(),
                factor_if_small(running // target_product),
            )
        )

    print(
        "GHOST_QUOTIENT",
        {
            "n": index,
            "core": (D, N),
            "margins": margins,
            "targets": targets,
            "carrier_bits": tuple(
                abs(carriers[key]).nbits() for key in sorted(carriers)
            ),
            "gcd_history": tuple(history),
            "smith_reduction": tuple(smith_records),
            "final_factor": factor_if_small(running),
        },
    )


for case in CASES:
    audit_case(*case)

print("Q32_GHOST_QUOTIENT_AUDIT=PASS")
