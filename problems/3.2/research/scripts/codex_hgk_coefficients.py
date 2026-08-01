#!/usr/bin/env python3
"""Verify the rank-two finite-field pullback and Jacobi convolution.

The checks are deliberately independent at the endpoints:

* ordinary binomial coefficients are reconstructed from Jacobi sums and from
  Morita-Gamma reductions (where Gross--Koblitz has no carry);
* the Franel Hasse polynomial is reconstructed from the 2F1 Hasse polynomial;
* the algebraic pullback is evaluated on every F_p-point, using F_{p^2} when
  the quadratic inverse does not split over F_p;
* Lagrange inversion recovers every relevant tau/sigma coefficient; and
* the resulting quadratic Mellin convolution recovers every b_r, 0 <= r < p.

Only Python's standard library is used.  The requested primes 13 and 29 cover
the sigma and tau branches respectively.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import comb


PRIMES = (13, 29)


def legendre(value: int, prime: int) -> int:
    value %= prime
    if value == 0:
        return 0
    answer = pow(value, (prime - 1) // 2, prime)
    return -1 if answer == prime - 1 else answer


def least_nonsquare(prime: int) -> int:
    return next(value for value in range(2, prime) if legendre(value, prime) == -1)


def square_root_mod(value: int, prime: int) -> int:
    value %= prime
    return next(root for root in range(prime) if root * root % prime == value)


@dataclass(frozen=True)
class Fp2:
    """Element a+b*w of F_p[w]/(w^2-nu), with nu a fixed nonsquare."""

    a: int
    b: int
    p: int
    nu: int

    def __post_init__(self) -> None:
        object.__setattr__(self, "a", self.a % self.p)
        object.__setattr__(self, "b", self.b % self.p)

    def coerce(self, other: int | Fp2) -> Fp2:
        if isinstance(other, Fp2):
            assert (other.p, other.nu) == (self.p, self.nu)
            return other
        return Fp2(other, 0, self.p, self.nu)

    def __add__(self, other: int | Fp2) -> Fp2:
        other = self.coerce(other)
        return Fp2(self.a + other.a, self.b + other.b, self.p, self.nu)

    __radd__ = __add__

    def __neg__(self) -> Fp2:
        return Fp2(-self.a, -self.b, self.p, self.nu)

    def __sub__(self, other: int | Fp2) -> Fp2:
        return self + (-self.coerce(other))

    def __rsub__(self, other: int | Fp2) -> Fp2:
        return self.coerce(other) - self

    def __mul__(self, other: int | Fp2) -> Fp2:
        other = self.coerce(other)
        return Fp2(
            self.a * other.a + self.b * other.b * self.nu,
            self.a * other.b + self.b * other.a,
            self.p,
            self.nu,
        )

    __rmul__ = __mul__

    def inverse(self) -> Fp2:
        norm = (self.a * self.a - self.nu * self.b * self.b) % self.p
        assert norm
        inverse_norm = pow(norm, -1, self.p)
        return Fp2(self.a * inverse_norm, -self.b * inverse_norm, self.p, self.nu)

    def __truediv__(self, other: int | Fp2) -> Fp2:
        return self * self.coerce(other).inverse()

    def __rtruediv__(self, other: int | Fp2) -> Fp2:
        return self.coerce(other) / self

    def __pow__(self, exponent: int) -> Fp2:
        if exponent < 0:
            return self.inverse() ** (-exponent)
        result = Fp2(1, 0, self.p, self.nu)
        base = self
        while exponent:
            if exponent & 1:
                result *= base
            base *= base
            exponent //= 2
        return result

    def is_base_field(self) -> bool:
        return self.b == 0


def fp2_sqrt_of_base(value: int, prime: int, nonsquare: int) -> Fp2:
    value %= prime
    if value == 0:
        return Fp2(0, 0, prime, nonsquare)
    if legendre(value, prime) == 1:
        return Fp2(square_root_mod(value, prime), 0, prime, nonsquare)
    quotient = value * pow(nonsquare, -1, prime) % prime
    return Fp2(0, square_root_mod(quotient, prime), prime, nonsquare)


def character_value(prime: int, exponent: int, value: int) -> int:
    """Reduction of omega^exponent, with every character extended by 0 at 0."""

    value %= prime
    if value == 0:
        return 0
    return pow(value, exponent % (prime - 1), prime)


def jacobi_sum_mod(prime: int, first: int, second: int) -> int:
    """J(bar(omega)^first, bar(omega)^second), reduced modulo p."""

    return sum(
        character_value(prime, -first, value)
        * character_value(prime, -second, 1 - value)
        for value in range(prime)
    ) % prime


def jacobi_binomial(prime: int, upper: int, lower: int) -> int:
    """A Jacobi-sum lift whose reduction is binom(upper, lower).

    At upper=p-1 and an endpoint, J(epsilon,epsilon) has the wrong endpoint
    normalization, so the two tautological endpoint values are separated.
    """

    assert 0 <= lower <= upper <= prime - 1
    if upper == 0 or (upper == prime - 1 and lower in (0, upper)):
        return 1
    return -jacobi_sum_mod(prime, lower, upper - lower) % prime


def gamma_at_fraction_mod(prime: int, numerator: int) -> int:
    """Gamma_p(numerator/(p-1)) modulo p, 0 <= numerator <= p-2."""

    assert 0 <= numerator <= prime - 2
    if numerator == 0:
        return 1
    representative = prime - numerator
    factorial = 1
    for value in range(1, representative):
        factorial = factorial * value % prime
    return (-1 if representative & 1 else 1) * factorial % prime


def gamma_binomial(prime: int, upper: int, lower: int) -> int:
    """The carry-free Gross--Koblitz Gamma ratio for upper <= p-2."""

    assert 1 <= upper <= prime - 2 and 0 <= lower <= upper
    numerator = (
        gamma_at_fraction_mod(prime, lower)
        * gamma_at_fraction_mod(prime, upper - lower)
    ) % prime
    return numerator * pow(gamma_at_fraction_mod(prime, upper), -1, prime) % prime


def hypergeometric_parameters(prime: int) -> tuple[int, int, int]:
    inverse_three = pow(3, -1, prime)
    first = -inverse_three % prime
    second = -2 * inverse_three % prime
    assert first + second == prime - 1
    return first, second, min(first, second)


def franel_coefficients(prime: int) -> list[int]:
    return [sum(comb(n, k) ** 3 for k in range(n + 1)) % prime for n in range(prime)]


def hypergeometric_hasse_coefficients(prime: int) -> list[int]:
    """Coefficients of the homogenized pulled-back 2F1 Hasse polynomial."""

    first, second, cutoff = hypergeometric_parameters(prime)
    result = [0] * prime
    for k in range(cutoff + 1):
        scalar = (
            jacobi_binomial(prime, first, k)
            * jacobi_binomial(prime, second, k)
            * pow(27, k, prime)
        ) % prime
        exponent = prime - 1 - 3 * k
        for offset in range(exponent + 1):
            degree = 2 * k + offset
            result[degree] = (
                result[degree]
                + scalar
                * jacobi_binomial(prime, exponent, offset)
                * pow(-2, offset, prime)
            ) % prime
    return result


def polynomial_value(coefficients: list[int], value: int, prime: int) -> int:
    result = 0
    for coefficient in reversed(coefficients):
        result = (result * value + coefficient) % prime
    return result


def polynomial_value_fp2(coefficients: list[int], value: Fp2) -> Fp2:
    result = value.coerce(0)
    for coefficient in reversed(coefficients):
        result = result * value + coefficient
    return result


def euler_2f1_value(prime: int, argument: Fp2) -> Fp2:
    """The single Kummer/Euler sum for P_p(z)."""

    first, second, _ = hypergeometric_parameters(prime)
    sign = -1 if (second + 1) & 1 else 1
    result = argument.coerce(0)
    for y in range(prime):
        result += (
            (argument.coerce(y) ** first)
            * (argument.coerce(1 - y) ** second)
            * ((1 - argument * y) ** first)
        )
    return sign * result


def hypergeometric_2f1_value(prime: int, argument: Fp2) -> Fp2:
    first, second, cutoff = hypergeometric_parameters(prime)
    result = argument.coerce(0)
    for k in range(cutoff + 1):
        result += (
            jacobi_binomial(prime, first, k)
            * jacobi_binomial(prime, second, k)
            * (argument**k)
        )
    return result


def homogenized_hasse_value(prime: int, x: Fp2) -> Fp2:
    first, second, cutoff = hypergeometric_parameters(prime)
    result = x.coerce(0)
    for k in range(cutoff + 1):
        result += (
            jacobi_binomial(prime, first, k)
            * jacobi_binomial(prime, second, k)
            * pow(27, k, prime)
            * (x ** (2 * k))
            * ((1 - 2 * x) ** (prime - 1 - 3 * k))
        )
    return result


def branch_coefficients(prime: int, branch: str) -> list[int]:
    half = (prime - 1) // 2
    degree = half if branch == "tau" else half - 1
    if branch == "tau":
        values = [1, 5 * pow(2, -1, prime) % prime]
    else:
        values = [1, 39 * pow(2, -1, prime) % prime]
    for n in range(1, degree):
        if branch == "tau":
            numerator = (
                2 * (68 * n * n + 34 * n + 5) * values[n]
                - (2 * n - 1) ** 2 * values[n - 1]
            )
        else:
            numerator = (
                2 * (68 * n * n + 102 * n + 39) * values[n]
                - (2 * n + 1) ** 2 * values[n - 1]
            )
        values.append(numerator * pow(4 * (n + 1) ** 2, -1, prime) % prime)
    return values[: degree + 1]


def branch_for_prime(prime: int) -> str:
    return "tau" if legendre(-6, prime) == 1 else "sigma"


def inverse_pullback_roots(prime: int, t: int) -> tuple[Fp2, Fp2]:
    nonsquare = least_nonsquare(prime)
    discriminant = (t * t - 34 * t + 1) % prime
    root_discriminant = fp2_sqrt_of_base(discriminant, prime, nonsquare)
    inverse_sixteen = pow(16, -1, prime)
    first = (1 - t + root_discriminant) * inverse_sixteen
    second = (1 - t - root_discriminant) * inverse_sixteen
    assert 8 * first * first + (t - 1) * first + t == first.coerce(0)
    assert 8 * second * second + (t - 1) * second + t == second.coerce(0)
    return first, second


def branch_pullback_value(prime: int, branch: str, x: Fp2) -> Fp2:
    half = (prime - 1) // 2
    hasse = homogenized_hasse_value(prime, x)
    if branch == "tau":
        return hasse / ((1 + x) ** half)
    correction = 1 - 16 * x - 8 * x * x
    return hasse / (correction * ((1 + x) ** (half - 1)))


def branch_values_from_pullback(prime: int, branch: str) -> list[int]:
    coefficients = branch_coefficients(prime, branch)
    values = [polynomial_value(coefficients, t, prime) for t in range(prime)]
    for t in range(prime):
        first, second = inverse_pullback_roots(prime, t)
        pulled_first = branch_pullback_value(prime, branch, first)
        pulled_second = branch_pullback_value(prime, branch, second)
        assert pulled_first == pulled_second
        assert pulled_first.is_base_field()
        assert pulled_first.a == values[t]

        # The Euler character sum is independently compared with the 2F1
        # Jacobi polynomial whenever the affine z-coordinate is defined.
        denominator = 1 - 2 * first
        if denominator != first.coerce(0):
            z = 27 * first * first / (denominator**3)
            assert euler_2f1_value(prime, z) == hypergeometric_2f1_value(prime, z)
            assert (denominator ** (prime - 1)) * euler_2f1_value(
                prime, z
            ) == homogenized_hasse_value(prime, first)
    return values


def phi_coefficient(prime: int, n: int, degree: int) -> int:
    return sum(
        jacobi_binomial(prime, n, index)
        * jacobi_binomial(prime, n + degree - index - 1, degree - index)
        * pow(8, degree - index, prime)
        for index in range(min(n, degree) + 1)
    ) % prime


def reciprocal_quadratic_coefficients(prime: int, limit: int) -> list[int]:
    values = [1]
    if limit:
        values.append(16 % prime)
    for degree in range(2, limit + 1):
        values.append((16 * values[-1] + 8 * values[-2]) % prime)
    return values


def pullback_series_coefficient(prime: int, branch: str, degree: int) -> int:
    """[x^degree] of H_p/(1+x)^e or H_p/(R(1+x)^(e-1))."""

    half = (prime - 1) // 2
    first, second, cutoff = hypergeometric_parameters(prime)
    reciprocal = reciprocal_quadratic_coefficients(prime, degree)
    result = 0
    for k in range(min(cutoff, degree // 2) + 1):
        scalar = (
            jacobi_binomial(prime, first, k)
            * jacobi_binomial(prime, second, k)
            * pow(27, k, prime)
        ) % prime
        exponent = prime - 1 - 3 * k
        for offset in range(min(exponent, degree - 2 * k) + 1):
            scalar_offset = (
                scalar
                * jacobi_binomial(prime, exponent, offset)
                * pow(-2, offset, prime)
            ) % prime
            remainder = degree - 2 * k - offset
            if branch == "tau":
                result += scalar_offset * jacobi_binomial(
                    prime, half + 1, remainder
                )
            else:
                result += scalar_offset * sum(
                    jacobi_binomial(prime, half + 2, power)
                    * reciprocal[remainder - power]
                    for power in range(remainder + 1)
                )
    return result % prime


def lagrange_branch_coefficients(prime: int, branch: str) -> list[int]:
    degree = (prime - 1) // 2 - (branch == "sigma")
    pullback = [
        pullback_series_coefficient(prime, branch, index)
        for index in range(degree + 1)
    ]
    result = [1]
    for n in range(1, degree + 1):
        result.append(
            sum(
                m * pullback[m] * phi_coefficient(prime, n, n - m)
                for m in range(1, n + 1)
            )
            * pow(n, -1, prime)
            % prime
        )
    return result


def apery_coefficients(prime: int) -> list[int]:
    return [
        sum(comb(n, k) ** 2 * comb(n + k, k) ** 2 for k in range(n + 1))
        % prime
        for n in range(prime)
    ]


def mellin_coefficients(prime: int, values: list[int], degree: int) -> list[int]:
    return [
        -sum(values[t] * pow(t, -index, prime) for t in range(1, prime)) % prime
        for index in range(degree + 1)
    ]


def mellin_kernel(prime: int, degree: int, r: int, first: int, second: int) -> int:
    lower = max(0, r - degree)
    upper = min(degree, r)
    return sum(
        pow(first, -index, prime) * pow(second, -(r - index), prime)
        for index in range(lower, upper + 1)
    ) % prime


def verify_quadratic_convolution(
    prime: int, branch: str, branch_values: list[int], coefficients: list[int]
) -> None:
    degree = len(coefficients) - 1
    apery = apery_coefficients(prime)
    convolution = []
    double_mellin = []
    for r in range(prime):
        lower = max(0, r - degree)
        upper = min(degree, r)
        convolution.append(
            sum(coefficients[i] * coefficients[r - i] for i in range(lower, upper + 1))
            % prime
        )
        double_mellin.append(
            sum(
                branch_values[t]
                * branch_values[u]
                * mellin_kernel(prime, degree, r, t, u)
                for t in range(1, prime)
                for u in range(1, prime)
            )
            % prime
        )
    assert convolution == double_mellin
    if branch == "tau":
        assert convolution == apery
    else:
        filtered = []
        for r, value in enumerate(convolution):
            corrected = value
            if r:
                corrected -= 34 * convolution[r - 1]
            if r >= 2:
                corrected += convolution[r - 2]
            filtered.append(corrected % prime)
        assert filtered == apery


def verify_prime(prime: int) -> None:
    # Jacobi and Gross--Koblitz reductions are checked independently against
    # Python's integer binomial coefficients.
    for upper in range(prime):
        for lower in range(upper + 1):
            expected = comb(upper, lower) % prime
            assert jacobi_binomial(prime, upper, lower) == expected
            if 1 <= upper <= prime - 2:
                assert gamma_binomial(prime, upper, lower) == expected

    franel = franel_coefficients(prime)
    hypergeometric = hypergeometric_hasse_coefficients(prime)
    assert hypergeometric == franel

    branch = branch_for_prime(prime)
    direct = branch_coefficients(prime, branch)
    pullback_values = branch_values_from_pullback(prime, branch)
    lagrange = lagrange_branch_coefficients(prime, branch)
    assert lagrange == direct
    mellin = mellin_coefficients(prime, pullback_values, len(direct) - 1)
    assert mellin == direct
    verify_quadratic_convolution(prime, branch, pullback_values, direct)

    print(
        f"p={prime}: branch={branch}, all {len(direct)} coefficients; "
        "Jacobi/Gamma, 2F1 pullback, Mellin inversion, and b_r convolution VERIFIED"
    )


def main() -> None:
    for prime in PRIMES:
        verify_prime(prime)


if __name__ == "__main__":
    main()
