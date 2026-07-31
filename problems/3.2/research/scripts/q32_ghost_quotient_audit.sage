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


HOSTILE_CASES = (
    (200, 128, 63),
    (272, 180, 63),
    (300, 180, 57),
    (321, 168, 53),
)
CENTRAL_CASES = (
    (80, 52, 16),
    (120, 78, 24),
    (160, 104, 32),
    (200, 130, 40),
)
if "--hybrid-dense" in sys.argv:
    CASES = tuple(
        (
            index,
            (13 * index) // 20,
            index // 5,
        )
        for index in range(60, 501, 20)
    )
elif "--hybrid-central" in sys.argv:
    CASES = CENTRAL_CASES + HOSTILE_CASES
else:
    CASES = HOSTILE_CASES
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
integral_common_operator = sum(
    ghost_coefficients[shift] * OAK(Sn)^shift
    for shift in range(len(ghost_coefficients))
)
twisted_apery_operator = sum(
    K(apery_operator[shift])
    / K(multiplier(n=n + shift))
    * OAK(Sn)^shift
    for shift in range(apery_operator.order() + 1)
)
origin_to_ghost_quotient, origin_to_ghost_remainder = (
    (twisted_apery_operator * integral_operator).quo_rem(
        integral_common_operator
    )
)
assert origin_to_ghost_remainder == 0
assert origin_to_ghost_quotient.order() == 1

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


def prime_to_set_part(value, primes):
    """Remove every prime factor in ``primes`` from ``value``."""

    answer = abs(ZZ(value))
    for prime in primes:
        while answer % prime == 0:
            answer //= prime
    return answer


def det2(left, right):
    return ZZ(left[0]) * ZZ(right[1]) - ZZ(left[1]) * ZZ(right[0])


def finite_difference(values, start, order):
    return ZZ(
        sum(
            (-1)^(order - shift)
            * binomial(order, shift)
            * values[start + shift]
            for shift in range(order + 1)
        )
    )


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
    # A target must lie in the endpoint-mask union for both retained
    # margins in order to survive their gcd.  Since these unions are
    # nested, their intersection is the union for the smaller margin.
    common_margin = min(margins)
    covered_lower = D - common_margin + 2
    covered_upper = D + N + common_margin - 1
    safe_covered_targets = tuple(
        prime
        for prime in prime_range(covered_lower, covered_upper + 1)
        if 2 * prime > index
        and index - prime + max(GHOST_SHIFTS) + common_operator.order()
        <= prime - 5
        and ZZ(apery(index - prime)) % prime == 0
    )
    safe_covered_target_product = prod(safe_covered_targets)

    Y = shell_batch(moment, nodes)
    maximum_time_shift = max(
        max(GHOST_SHIFTS) + common_operator.order(),
        5 + origin_operator.order(),
    )
    time_shell_rows = {
        time_index: shell_batch(time_index, neighbours)
        for time_index in range(index, index + maximum_time_shift + 1)
    }
    ghosts = {}
    for ghost_shift in GHOST_SHIFTS:
        shell_rows = {
            operator_shift: time_shell_rows[
                index + ghost_shift + operator_shift
            ]
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

    # Intersect the quotient-ghost family with the independently certified
    # maximal-boundary family from Section 68.11.  The k=0 common-operator
    # shell rows already contain every moment needed for the origin-cancelled
    # coordinate W, so no second shell computation is required.
    origin_lifts = {
        lift_shift: {
            node: sum(
                ZZ(
                    integer_coefficients[operator_shift](
                        n=index + lift_shift
                    )
                )
                * (
                    time_shell_rows[
                        index + lift_shift + operator_shift
                    ][node - 1]
                    + time_shell_rows[
                        index + lift_shift + operator_shift
                    ][node + 1]
                )
                for operator_shift in range(
                    origin_operator.order() + 1
                )
            )
            for node in nodes
        }
        for lift_shift in range(6)
    }
    origin_coordinate = {
        node: origin_lifts[0][node]
        - ZZ(multiplier(n=index)) * Y[node]
        for node in nodes
    }
    boundary_carriers = []
    boundary_by_margin = {}
    boundary_data = {}
    for margin in margins:
        d = D - margin + 1
        length = N + margin - 2
        y_left, y_right = data[0, margin][1]
        w_left = endpoint(origin_coordinate, d, length)
        w_right = endpoint(origin_coordinate, D, length)
        normalizer = margin_normalizer(D, N, margin)
        determinant = y_left * w_right - y_right * w_left
        assert determinant % normalizer == 0
        boundary = determinant // normalizer
        assert boundary % target_product == 0
        boundary_carriers.append(boundary)
        boundary_by_margin[margin] = boundary
        boundary_data[margin] = (
            d,
            length,
            y_left,
            y_right,
            normalizer,
        )

    # The exact Ore identity A_M U = Q L implies that, modulo a prime
    # annihilating the ghost rows, every lifted determinant below follows
    # the same order-two recurrence in ``lift_shift``.  Its two-margin
    # time Casoratians are therefore a new target-preserving family.
    lifted_boundaries = {}
    for margin in margins:
        d, length, y_left, y_right, normalizer = boundary_data[margin]
        for lift_shift in range(6):
            z_left = endpoint(origin_lifts[lift_shift], d, length)
            z_right = endpoint(origin_lifts[lift_shift], D, length)
            determinant = y_left * z_right - y_right * z_left
            assert determinant % normalizer == 0
            lifted_boundaries[lift_shift, margin] = (
                determinant // normalizer
            )
        assert lifted_boundaries[0, margin] == boundary_carriers[
            margins.index(margin)
        ]
        for lift_shift in range(4):
            assert sum(
                K(twisted_apery_operator[operator_shift])(
                    n=index + lift_shift
                )
                * lifted_boundaries[
                    lift_shift + operator_shift,
                    margin,
                ]
                for operator_shift in range(
                    twisted_apery_operator.order() + 1
                )
            ) == sum(
                K(origin_to_ghost_quotient[operator_shift])(
                    n=index + lift_shift
                )
                * carriers[lift_shift + operator_shift, margin]
                for operator_shift in range(
                    origin_to_ghost_quotient.order() + 1
                )
            )

    lift_crosses = {}
    for left_shift in range(1, 6):
        for right_shift in range(left_shift + 1, 6):
            cross = (
                lifted_boundaries[left_shift, margins[0]]
                * lifted_boundaries[right_shift, margins[1]]
                - lifted_boundaries[right_shift, margins[0]]
                * lifted_boundaries[left_shift, margins[1]]
            )
            lift_crosses[left_shift, right_shift] = cross
    lift_cross_gcd = gcd_many(lift_crosses.values())

    # After division by the candidate primorial, the origin-cancelled
    # boundary column is a sixth ghost quotient column.  The second
    # determinantal divisor of this 2-by-6 matrix gives the exact Smith
    # obstruction to eliminating a common carrier prime.
    six_column_records = []
    six_column_deltas = {}
    for margin in margins:
        d, length, y_left, y_right, normalizer = boundary_data[margin]
        reduced_normalizer = data[0, margin][3]
        assert normalizer == candidate_product * reduced_normalizer
        w_left = endpoint(origin_coordinate, d, length)
        w_right = endpoint(origin_coordinate, D, length)
        assert w_left % candidate_product == 0
        assert w_right % candidate_product == 0
        quotient_columns = [
            (
                w_left // candidate_product,
                w_right // candidate_product,
            )
        ] + [
            data[ghost_shift, margin][2]
            for ghost_shift in GHOST_SHIFTS
        ]
        delta = gcd_many(
            det2(quotient_columns[left], quotient_columns[right])
            for left in range(len(quotient_columns))
            for right in range(left + 1, len(quotient_columns))
        )
        gamma = gcd_many((y_left, y_right))
        kappa = gcd_many(
            [boundary_by_margin[margin]]
            + [
                carriers[ghost_shift, margin]
                for ghost_shift in GHOST_SHIFTS
            ]
        )
        assert (gamma * delta) % (reduced_normalizer * kappa) == 0
        six_column_deltas[margin] = delta
        six_column_records.append(
            (
                margin,
                abs(reduced_normalizer).nbits(),
                abs(gamma).nbits(),
                abs(delta).nbits(),
                abs(kappa).nbits(),
                factor_if_small(delta),
                factor_if_small(gcd(delta, candidate_product)),
            )
        )

    origin_content = gcd_many(
        [ZZ(multiplier(n=index))]
        + [
            ZZ(coefficient(n=index))
            for coefficient in integer_coefficients
        ]
    )
    safe_origin_content = prime_to_set_part(origin_content, candidates)
    assert all(value % safe_origin_content == 0 for value in boundary_carriers)
    boundary_gcd = gcd_many(boundary_carriers) // safe_origin_content
    assert boundary_gcd % target_product == 0
    hybrid_gcd = gcd(boundary_gcd, running)
    assert hybrid_gcd % target_product == 0
    assert hybrid_gcd % safe_covered_target_product == 0
    augmented_gcd = gcd(hybrid_gcd, lift_cross_gcd)

    terminal_margin = maximum_margin
    terminal_d = D - terminal_margin + 1
    terminal_length = N + terminal_margin - 2
    left_boundary_packet = finite_difference(
        Y,
        terminal_d,
        terminal_length,
    )
    right_boundary_packet = finite_difference(
        Y,
        D,
        terminal_length,
    )
    left_pascal_multiplier = binomial(D + N, terminal_length)
    right_pascal_multiplier = binomial(
        D + terminal_length,
        terminal_length,
    )
    previous_margin = terminal_margin - 1
    previous_y_left, previous_y_right = data[0, previous_margin][1]
    terminal_y_left, terminal_y_right = data[0, terminal_margin][1]
    assert (
        terminal_y_left - previous_y_left
        == (-1)^terminal_length
        * left_pascal_multiplier
        * left_boundary_packet
    )
    assert (
        terminal_y_right - previous_y_right
        == (-1)^terminal_length
        * right_pascal_multiplier
        * right_boundary_packet
    )
    boundary_packet_gcd = gcd(
        left_boundary_packet,
        right_boundary_packet,
    )
    margin_content_product = prod(
        data[0, margin][3] for margin in margins
    )
    six_column_product = prod(
        six_column_deltas[margin] for margin in margins
    )
    pascal_product = left_pascal_multiplier * right_pascal_multiplier
    support_certificate = (
        safe_origin_content
        * margin_content_product
        * six_column_product
        * pascal_product
        * boundary_packet_gcd
    )
    hybrid_prime_support = tuple(
        prime for prime, _ in factor(abs(hybrid_gcd))
    )
    assert all(
        support_certificate % prime == 0
        for prime in hybrid_prime_support
    )

    print(
        "GHOST_QUOTIENT",
        {
            "n": index,
            "core": (D, N),
            "margins": margins,
            "targets": targets,
            "covered_prime_range": (covered_lower, covered_upper),
            "safe_covered_targets": safe_covered_targets,
            "carrier_bits": tuple(
                abs(carriers[key]).nbits() for key in sorted(carriers)
            ),
            "gcd_history": tuple(history),
            "smith_reduction": tuple(smith_records),
            "final_factor": factor_if_small(running),
            "boundary_factor": factor_if_small(boundary_gcd),
            "hybrid_bits": abs(hybrid_gcd).nbits(),
            "hybrid_factor": factor_if_small(hybrid_gcd),
            "hybrid_residual_factor": factor_if_small(
                hybrid_gcd // target_product
            ),
            "hybrid_safe_residual_factor": factor_if_small(
                hybrid_gcd // safe_covered_target_product
            ),
            "lift_cross_records": tuple(
                (
                    shifts,
                    abs(value).nbits(),
                    value % safe_covered_target_product == 0,
                    factor_if_small(value),
                )
                for shifts, value in sorted(lift_crosses.items())
            ),
            "lift_cross_gcd_bits": abs(lift_cross_gcd).nbits(),
            "lift_cross_gcd_factor": factor_if_small(lift_cross_gcd),
            "augmented_gcd_bits": abs(augmented_gcd).nbits(),
            "augmented_gcd_factor": factor_if_small(augmented_gcd),
            "six_column_records": tuple(six_column_records),
            "boundary_packet_bits": (
                abs(left_boundary_packet).nbits(),
                abs(right_boundary_packet).nbits(),
                abs(boundary_packet_gcd).nbits(),
            ),
            "boundary_packet_gcd_factor": factor_if_small(
                boundary_packet_gcd
            ),
            "hybrid_support_certificate_bits": abs(
                support_certificate
            ).nbits(),
            "hybrid_support_ledger": tuple(
                (
                    prime,
                    tuple(
                        label
                        for label, value in (
                            ("origin", safe_origin_content),
                            ("margin-content", margin_content_product),
                            ("six-column", six_column_product),
                            ("pascal", pascal_product),
                            ("boundary", boundary_packet_gcd),
                        )
                        if value % prime == 0
                    ),
                )
                for prime in hybrid_prime_support
            ),
        },
    )


for case in CASES:
    audit_case(*case)

print("Q32_GHOST_QUOTIENT_AUDIT=PASS")
