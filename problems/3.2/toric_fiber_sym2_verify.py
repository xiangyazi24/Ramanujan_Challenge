#!/usr/bin/env python3
"""Exact checks for toric_fiber_sym2.tex.

The symbolic checks certify the Beukers--Peters change of variables, the
quadratic deck law, the two symmetric j-invariant identities, and the CM
specializations.  The finite-field check covers every nonzero fiber for every
prime 5 <= p <= 101, using the exact quadratic character formula for mu_p(a).
"""

from __future__ import annotations

from hashlib import sha256

import sympy as sp


def primes_through(limit: int) -> list[int]:
    return [
        n
        for n in range(2, limit + 1)
        if all(n % d for d in range(2, int(n**0.5) + 1))
    ]


def legendre(value: int, prime: int) -> int:
    value %= prime
    if value == 0:
        return 0
    return 1 if pow(value, (prime - 1) // 2, prime) == 1 else -1


def symbolic_checks() -> None:
    x, y, z, a = sp.symbols("x y z a")
    cap_x = y / (1 + y)
    cap_y = z / (1 + z)
    cap_z = x / (1 + x)
    bp_ratio = (
        1 - (1 - cap_x * cap_y) * cap_z
    ) / (
        cap_x
        * cap_y
        * cap_z
        * (1 - cap_x)
        * (1 - cap_y)
        * (1 - cap_z)
    )
    laurent = (
        (1 + x)
        * (1 + y)
        * (1 + z)
        * ((1 + y) * (1 + z) + x * y * z)
        / (x * y * z)
    )
    assert sp.factor(bp_ratio - laurent) == 0

    sigma = (1 - 8 * x) / (8 * (1 + x))
    phi = x * (1 - 8 * x) / (1 + x)
    assert sp.factor(sigma.subs(x, sigma) - x) == 0
    assert sp.factor(phi.subs(x, sigma) - phi) == 0

    def elliptic_j(parameter: sp.Expr) -> sp.Expr:
        a1 = 1 - 2 * parameter
        a3 = parameter**2
        b2 = a1**2
        b4 = a1 * a3
        b6 = a3**2
        c4 = b2**2 - 24 * b4
        discriminant = -8 * b4**3 - 27 * b6**2 + 9 * b2 * b4 * b6
        return sp.factor(c4**3 / discriminant)

    jx = elliptic_j(x)
    js = elliptic_j(sigma)
    cover = 8 * x**2 + (a - 1) * x + a

    def reduce_on_cover(expression: sp.Expr) -> sp.Expr:
        numerator, denominator = sp.cancel(expression).as_numer_denom()
        numerator = sp.rem(sp.Poly(numerator, x), sp.Poly(cover, x)).as_expr()
        denominator = sp.rem(
            sp.Poly(denominator, x), sp.Poly(cover, x)
        ).as_expr()
        inverse = sp.invert(sp.Poly(denominator, x), sp.Poly(cover, x)).as_expr()
        return sp.factor(
            sp.rem(sp.Poly(numerator * inverse, x), sp.Poly(cover, x)).as_expr()
        )

    expected_product = (
        (a + 1) ** 3 * (a**3 - 21 * a**2 + 219 * a + 1) ** 3 / a**7
    )
    expected_shifted = (
        (a**2 + 14 * a + 1) ** 2
        * (a**4 - 44 * a**3 + 198 * a**2 - 548 * a + 1) ** 2
        / a**7
    )
    assert sp.factor(reduce_on_cover(jx * js) - expected_product) == 0
    assert sp.factor(
        reduce_on_cover((jx - 1728) * (js - 1728)) - expected_shifted
    ) == 0

    # Bertin--Lecacheux's parameter may be chosen with s^2=a^{-1}.
    bl_product = (
        (a**-1 + 1) ** 3
        * (a**-3 + 219 * a**-2 - 21 * a**-1 + 1) ** 3
        / a**-5
    )
    bl_shifted = (
        (a**-2 + 14 * a**-1 + 1) ** 2
        * (a**-4 - 548 * a**-3 + 198 * a**-2 - 44 * a**-1 + 1) ** 2
        / a**-5
    )
    assert sp.factor(bl_product - expected_product) == 0
    assert sp.factor(bl_shifted - expected_shifted) == 0

    branch = 8 * x**2 + 16 * x - 1
    branch_j_numerator, branch_j_denominator = sp.cancel(jx).as_numer_denom()
    hilbert_numerator = sp.expand(
        branch_j_numerator**2
        - 4834944 * branch_j_numerator * branch_j_denominator
        + 14670139392 * branch_j_denominator**2
    )
    assert sp.rem(sp.Poly(hilbert_numerator, x), sp.Poly(branch, x)) == 0
    assert sp.rem(
        sp.Poly(16 * (1 + x) ** 2 - 18, x), sp.Poly(branch, x)
    ) == 0

    center = 8 * x**2 + 1
    center_identity = sp.together(jx - 8000)
    center_numerator = center_identity.as_numer_denom()[0]
    assert sp.rem(sp.Poly(center_numerator, x), sp.Poly(center, x)) == 0


def elliptic_trace(parameter: int, prime: int) -> int:
    a1 = (1 - 2 * parameter) % prime
    a3 = parameter * parameter % prime
    character_sum = 0
    for x_value in range(prime):
        linear = (a1 * x_value + a3) % prime
        discriminant = (linear * linear + 4 * x_value**3) % prime
        character_sum += legendre(discriminant, prime)
    return -character_sum


def elliptic_trace_quadratic(parameter: int, prime: int) -> int:
    """Trace over F_{p^2} at a root of 8x^2+(a-1)x+a."""

    inverse_eight = pow(8, -1, prime)
    relation_constant = -parameter * inverse_eight % prime
    relation_linear = (1 - parameter) * inverse_eight % prime

    def add(left: tuple[int, int], right: tuple[int, int]) -> tuple[int, int]:
        return ((left[0] + right[0]) % prime, (left[1] + right[1]) % prime)

    def scale(scalar: int, value: tuple[int, int]) -> tuple[int, int]:
        return (scalar * value[0] % prime, scalar * value[1] % prime)

    def multiply(
        left: tuple[int, int], right: tuple[int, int]
    ) -> tuple[int, int]:
        constant = (
            left[0] * right[0]
            + left[1] * right[1] * relation_constant
        ) % prime
        linear = (
            left[0] * right[1]
            + left[1] * right[0]
            + left[1] * right[1] * relation_linear
        ) % prime
        return constant, linear

    def quadratic_character(value: tuple[int, int]) -> int:
        if value == (0, 0):
            return 0
        # Norm(c+d*r)=c^2+c*d*(r+r')+d^2*r*r'.
        norm = (
            value[0] * value[0]
            + value[0] * value[1] * relation_linear
            - value[1] * value[1] * relation_constant
        ) % prime
        return legendre(norm, prime)

    root = (0, 1)
    root_square = (relation_constant, relation_linear)
    a1 = (1, -2 % prime)
    a3 = root_square
    character_sum = 0
    for constant in range(prime):
        for linear in range(prime):
            x_value = (constant, linear)
            x_square = multiply(x_value, x_value)
            x_cube = multiply(x_square, x_value)
            linear_y = add(multiply(a1, x_value), a3)
            discriminant = add(
                multiply(linear_y, linear_y), scale(4, x_cube)
            )
            character_sum += quadratic_character(discriminant)
    return -character_sum


def toric_mu(parameter: int, prime: int) -> int:
    inverses = [0] + [pow(value, -1, prime) for value in range(1, prime)]
    character_sum = 0
    for product_uv in range(1, prime):
        for w_value in range(1, prime):
            if w_value == 1 or product_uv * w_value % prime == 1:
                continue
            correction = (
                product_uv
                * (w_value - 1)
                * (1 - product_uv * w_value)
                * inverses[parameter]
                * inverses[w_value]
            ) % prime
            sum_uv = (product_uv + 1 - correction) % prime
            character_sum += legendre(
                sum_uv * sum_uv - 4 * product_uv, prime
            )
    return (prime - 2) ** 2 + character_sum


def finite_field_checks(limit: int = 101) -> None:
    counts = {"split": 0, "inert": 0, "branch": 0, "a1": 0}
    digest = sha256()
    total = 0

    for prime in primes_through(limit):
        if prime < 5:
            continue
        twist = legendre(-3, prime)
        for parameter in range(1, prime):
            discriminant = (parameter * parameter - 34 * parameter + 1) % prime
            epsilon = int(parameter == 1)
            counts["a1"] += epsilon
            actual = toric_mu(parameter, prime) - (
                prime * prime - 6 * prime + 12
            )

            if discriminant == 0:
                kind = "branch"
                root = -(parameter - 1) * pow(16, -1, prime) % prime
                trace = elliptic_trace(root, prime)
                predicted = (
                    trace * trace
                    - prime
                    - legendre(-6, prime) * prime
                )
            elif legendre(discriminant, prime) == 1:
                kind = "split"
                roots = [
                    value
                    for value in range(prime)
                    if (8 * value * value + (parameter - 1) * value + parameter)
                    % prime
                    == 0
                ]
                assert len(roots) == 2
                traces = [elliptic_trace(root, prime) for root in roots]
                assert traces[1] == twist * traces[0]
                predicted = traces[0] ** 2 - prime - epsilon * prime
            else:
                kind = "inert"
                trace = elliptic_trace_quadratic(parameter, prime)
                predicted = twist * trace - prime - epsilon * prime

            assert actual == predicted, (
                prime,
                parameter,
                kind,
                actual,
                predicted,
            )
            counts[kind] += 1
            total += 1
            digest.update(
                f"{prime}:{parameter}:{kind}:{actual}:{predicted}\n".encode()
            )

    assert counts == {"split": 531, "inert": 579, "branch": 22, "a1": 24}
    assert total == 1132
    expected_digest = (
        "1ce582dc5591366b87ecbf348148b70e5288d646f112e5d29a3fc2ef0e24b4fa"
    )
    assert digest.hexdigest() == expected_digest
    print("toric_fiber_sym2_verify: PASS")
    print(f"counts={counts} total={total}")
    print(f"sha256={digest.hexdigest()}")


if __name__ == "__main__":
    symbolic_checks()
    finite_field_checks()
