#!/usr/bin/env python3
"""Exact finite-field checks for toric_fiber_k3.tex."""

from itertools import product
from math import comb


def primes_through(limit: int):
    return [
        n
        for n in range(2, limit + 1)
        if all(n % d for d in range(2, int(n**0.5) + 1))
    ]


def projective_line(prime: int):
    return [(1, value) for value in range(prime)] + [(0, 1)]


def projective_value(point, prime: int):
    if point[0] == 0:
        return None
    return point[1] * pow(point[0], prime - 2, prime) % prime


def homogeneous_value(u, v, w, parameter: int, prime: int) -> int:
    u0, u1 = u
    v0, v1 = v
    w0, w1 = w
    return (
        -parameter * u0**2 * v0**2 * w0 * w1
        + parameter * u0**2 * v0 * v1 * w0 * w1
        + parameter * u0 * u1 * v0**2 * w0 * w1
        - u0 * u1 * v0 * v1 * w0**2
        + (1 - parameter) * u0 * u1 * v0 * v1 * w0 * w1
        + u1**2 * v1**2 * w0 * w1
        - u1**2 * v1**2 * w1**2
    ) % prime


def homogeneous_terms(parameter: int):
    return [
        (-parameter, (2, 0, 2, 0, 1, 1)),
        (parameter, (2, 0, 1, 1, 1, 1)),
        (parameter, (1, 1, 2, 0, 1, 1)),
        (-1, (1, 1, 1, 1, 2, 0)),
        (1 - parameter, (1, 1, 1, 1, 1, 1)),
        (1, (0, 2, 0, 2, 1, 1)),
        (-1, (0, 2, 0, 2, 0, 2)),
    ]


def local_gradient(u, v, w, parameter: int, prime: int):
    coordinates = u + v + w
    local_indices = [1 if point[0] else 0 for point in (u, v, w)]
    result = []
    for coordinate_index in local_indices:
        # Shift the coordinate index into the appropriate P1 factor.
        axis = len(result)
        differentiated_index = 2 * axis + coordinate_index
        derivative = 0
        for coefficient, exponents in homogeneous_terms(parameter):
            exponent = exponents[differentiated_index]
            if exponent == 0:
                continue
            term = coefficient * exponent
            for index, (value, power) in enumerate(zip(coordinates, exponents)):
                if index == differentiated_index:
                    power -= 1
                term *= value**power
            derivative += term
        result.append(derivative % prime)
    return result


def is_open(u, v, w, prime: int) -> bool:
    values = [projective_value(point, prime) for point in (u, v, w)]
    if any(value is None or value in (0, 1) for value in values):
        return False
    return values[0] * values[1] * values[2] % prime != 1


def equals(point, value, prime: int) -> bool:
    return projective_value(point, prime) == value


def boundary_memberships(u, v, w, prime: int):
    u0, u1 = u
    v0, v1 = v
    w0, w1 = w
    infinity = None
    return [
        equals(u, 0, prime) and equals(v, infinity, prime),
        equals(u, 0, prime) and equals(w, infinity, prime),
        equals(u, 0, prime) and equals(w, 0, prime),
        equals(u, 0, prime) and equals(v, 1, prime),
        equals(v, 0, prime) and equals(u, infinity, prime),
        equals(v, 0, prime) and equals(w, infinity, prime),
        equals(v, 0, prime) and equals(w, 0, prime),
        equals(v, 0, prime) and equals(u, 1, prime),
        equals(u, infinity, prime) and equals(w, 0, prime),
        equals(u, infinity, prime) and equals(w, 1, prime),
        equals(v, infinity, prime) and equals(w, 0, prime),
        equals(v, infinity, prime) and equals(w, 1, prime),
        equals(u, 1, prime) and equals(w, 1, prime),
        equals(u, 1, prime) and (v1 * w1 - v0 * w0) % prime == 0,
        equals(v, 1, prime) and equals(w, 1, prime),
        equals(v, 1, prime) and (u1 * w1 - u0 * w0) % prime == 0,
    ]


def polynomial_multiply(left, right, prime: int):
    result = {}
    for left_degree, left_value in left.items():
        for right_degree, right_value in right.items():
            degree = tuple(a + b for a, b in zip(left_degree, right_degree))
            result[degree] = (
                result.get(degree, 0) + left_value * right_value
            ) % prime
    return result


def univariate_power(constant: int, slope: int, exponent: int, axis: int,
                     prime: int):
    result = {}
    for degree in range(exponent + 1):
        key = [0, 0, 0]
        key[axis] = degree
        result[tuple(key)] = (
            comb(exponent, degree)
            * pow(constant, exponent - degree, prime)
            * pow(slope, degree, prime)
        ) % prime
    return result


def local_coordinate_pair(point):
    if point[0] == 0:
        return ((0, 1), (1, 0))
    return ((1, 0), (point[1], 1))


def local_polynomial(u, v, w, parameter: int, prime: int):
    coordinate_pairs = (
        local_coordinate_pair(u)
        + local_coordinate_pair(v)
        + local_coordinate_pair(w)
    )
    result = {}
    for coefficient, exponents in homogeneous_terms(parameter):
        term = {(0, 0, 0): coefficient % prime}
        for coordinate, exponent in enumerate(exponents):
            axis = coordinate // 2
            constant, slope = coordinate_pairs[coordinate]
            term = polynomial_multiply(
                term,
                univariate_power(constant, slope, exponent, axis, prime),
                prime,
            )
        for degree, value in term.items():
            result[degree] = (result.get(degree, 0) + value) % prime
    return result


def matrix_rank(matrix, prime: int) -> int:
    work = [[entry % prime for entry in row] for row in matrix]
    rank = 0
    for column in range(len(work[0])):
        pivot = next(
            (row for row in range(rank, len(work)) if work[row][column]),
            None,
        )
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        inverse = pow(work[rank][column], prime - 2, prime)
        work[rank] = [entry * inverse % prime for entry in work[rank]]
        for row in range(len(work)):
            if row == rank:
                continue
            multiplier = work[row][column]
            work[row] = [
                (entry - multiplier * pivot_entry) % prime
                for entry, pivot_entry in zip(work[row], work[rank])
            ]
        rank += 1
    return rank


def hessian_rank(polynomial, prime: int) -> int:
    coefficient = lambda degree: polynomial.get(degree, 0)
    hessian = [
        [2 * coefficient((2, 0, 0)), coefficient((1, 1, 0)),
         coefficient((1, 0, 1))],
        [coefficient((1, 1, 0)), 2 * coefficient((0, 2, 0)),
         coefficient((0, 1, 1))],
        [coefficient((1, 0, 1)), coefficient((0, 1, 1)),
         2 * coefficient((0, 0, 2))],
    ]
    return matrix_rank(hessian, prime)


def cubic_on_direction(polynomial, direction, prime: int) -> int:
    return sum(
        value
        * pow(direction[0], degree[0], prime)
        * pow(direction[1], degree[1], prime)
        * pow(direction[2], degree[2], prime)
        for degree, value in polynomial.items()
        if sum(degree) == 3
    ) % prime


def singular(polynomial, prime: int) -> bool:
    return all(
        polynomial.get(degree, 0) % prime == 0
        for degree in ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    )


def lambda_value(x: int, y: int, z: int, prime: int) -> int:
    numerator = (
        (1 + x)
        * (1 + y)
        * (1 + z)
        * ((1 + y) * (1 + z) + x * y * z)
    ) % prime
    denominator = x * y * z % prime
    return numerator * pow(denominator, prime - 2, prime) % prime


def direct_fibers(prime: int):
    result = [0] * prime
    for x, y, z in product(range(1, prime), repeat=3):
        result[lambda_value(x, y, z, prime)] += 1
    return result


def point(value, prime: int):
    if value is None:
        return (0, 1)
    return (1, value % prime)


def main() -> None:
    boundary_checks = 0
    fiber_checks = 0
    singularity_checks = 0
    trace_checks = 0

    for prime in primes_through(23):
        if prime < 5:
            continue
        line = projective_line(prime)
        multiplicities = {}
        component_sizes = [0] * 16
        boundary = set()
        for u, v, w in product(line, repeat=3):
            memberships = boundary_memberships(u, v, w, prime)
            multiplicity = sum(memberships)
            if multiplicity:
                boundary.add((u, v, w))
                multiplicities[multiplicity] = (
                    multiplicities.get(multiplicity, 0) + 1
                )
                for index, belongs in enumerate(memberships):
                    component_sizes[index] += int(belongs)
        assert component_sizes == [prime + 1] * 16
        assert multiplicities == {
            1: 16 * prime - 32,
            2: 16,
            3: 4,
            4: 1,
        }
        assert len(boundary) == 16 * prime - 11
        boundary_checks += 1

        direct = direct_fibers(prime)
        for parameter in range(1, prime):
            projective_count = 0
            open_count = 0
            singular_points = []
            for u, v, w in product(line, repeat=3):
                if homogeneous_value(u, v, w, parameter, prime):
                    continue
                projective_count += 1
                if is_open(u, v, w, prime):
                    open_count += 1
                if local_gradient(u, v, w, parameter, prime) == [0, 0, 0]:
                    polynomial = local_polynomial(u, v, w, parameter, prime)
                    assert singular(polynomial, prime)
                    singular_points.append((u, v, w, polynomial))

            assert open_count == direct[parameter]
            assert projective_count - open_count == 16 * prime - 11
            assert all(
                not homogeneous_value(u, v, w, parameter, prime)
                for u, v, w in boundary
            )
            fiber_checks += 1

            fixed_a1 = {
                (point(0, prime), point(1, prime), point(None, prime)),
                (point(0, prime), point(None, prime), point(None, prime)),
                (point(1, prime), point(0, prime), point(None, prime)),
                (point(None, prime), point(0, prime), point(None, prime)),
            }
            fixed_a2 = {
                (point(0, prime), point(None, prime), point(0, prime)),
                (point(None, prime), point(0, prime), point(0, prime)),
            }
            central = (point(1, prime), point(1, prime), point(1, prime))
            expected = fixed_a1 | fixed_a2 | {central}
            delta = int((parameter * parameter - 34 * parameter + 1) % prime == 0)
            if delta:
                inverse_twelve = pow(12, prime - 2, prime)
                critical_u = (parameter - 5) * inverse_twelve % prime
                critical_w = (29 - parameter) * inverse_twelve % prime
                expected.add(
                    (point(critical_u, prime), point(critical_u, prime),
                     point(critical_w, prime))
                )
            observed = {(u, v, w) for u, v, w, _ in singular_points}
            assert observed == expected

            for u, v, w, polynomial in singular_points:
                location = (u, v, w)
                rank = hessian_rank(polynomial, prime)
                if location in fixed_a1:
                    assert rank == 3
                elif location in fixed_a2:
                    assert rank == 2
                    direction = (
                        (parameter, 0, 1)
                        if location == (point(0, prime), point(None, prime),
                                        point(0, prime))
                        else (0, parameter, 1)
                    )
                    assert cubic_on_direction(polynomial, direction, prime)
                elif location == central:
                    if parameter == 1:
                        assert rank == 2
                        assert cubic_on_direction(
                            polynomial, (-1, -1, 1), prime
                        )
                    else:
                        assert rank == 3
                else:
                    assert delta and rank == 3
                singularity_checks += 1

            epsilon = int(parameter == 1)
            resolution_count = (
                projective_count + (9 + epsilon + delta) * prime
            )
            trace = resolution_count - 1 - prime * prime
            character_sum = open_count - (prime - 2) ** 2
            assert character_sum == (
                trace - (21 + epsilon + delta) * prime + 8
            )
            assert abs(trace) <= 22 * prime
            deflated_trace = trace - (12 + epsilon + delta) * prime
            assert character_sum == deflated_trace - 9 * prime + 8
            assert abs(deflated_trace) <= (10 - epsilon - delta) * prime
            assert abs(character_sum) <= 19 * prime + 8
            assert not (epsilon and delta)
            trace_checks += 1

        for parameter in range(1, prime):
            assert direct[parameter] == direct[pow(parameter, prime - 2, prime)]

    print("toric_fiber_k3_verify: PASS")
    print(f"boundary_checks={boundary_checks}")
    print(f"fiber_checks={fiber_checks}")
    print(f"singularity_checks={singularity_checks}")
    print(f"trace_checks={trace_checks}")


if __name__ == "__main__":
    main()
