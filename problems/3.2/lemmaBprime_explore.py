#!/usr/bin/env python3
"""Exact exploration for Lemma B' (Hodge/Newton/eigenvalue gate).

This file deliberately keeps two decompositions separate.

* ``basechanged_a1_trace_padic`` evaluates the *full* q-1 term A1 Mellin
  skeleton over F_{p^r}, compressed by cyclic convolution.  The q-1 terms are
  not Frobenius eigenvalues.
* ``hodge_data`` and ``frobenius_slopes`` concern the rank-three, cancelled
  middle-extension Frobenius factor.  Its Hodge slopes are 0,1,3.

The portable path uses only the Python standard library and includes a minimal
implementation of F_{p^r}.  The optional production scan imports Sage for fast
cyclic polynomial products.  The long-range table is a theorem-driven
enumeration from the cubic factor, with p-adic base-field A1 checks at every
Hasse zero.
"""

from __future__ import annotations

import argparse
import itertools
import random
from dataclasses import dataclass
from fractions import Fraction
from math import comb


# ---------------------------------------------------------------------------
# Elementary data: primes and Apéry residues


def primes_up_to(limit: int) -> list[int]:
    if limit < 2:
        return []
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[0:2] = b"\x00\x00"
    for q in range(2, int(limit**0.5) + 1):
        if sieve[q]:
            sieve[q * q : limit + 1 : q] = b"\x00" * (
                (limit - q * q) // q + 1
            )
    return [q for q in range(2, limit + 1) if sieve[q]]


def apery_residues(p: int) -> list[int]:
    """b_j modulo p for 0 <= j <= p-2, from the Apéry recurrence."""
    values = [0] * (p - 1)
    values[0] = 1
    values[1] = 5 % p
    for n in range(1, p - 2):
        coefficient = (34 * n**3 + 51 * n**2 + 27 * n + 5) % p
        rhs = (coefficient * values[n] - n**3 * values[n - 1]) % p
        values[n + 1] = rhs * pow((n + 1) ** 3, -1, p) % p
    return values


def apery_binomial_mod(n: int, p: int) -> int:
    return sum(comb(n, k) ** 2 * comb(n + k, k) ** 2 for k in range(n + 1)) % p


# ---------------------------------------------------------------------------
# Part 1: corrected AS/GKZ fractional-part computation


@dataclass(frozen=True)
class HodgeData:
    a: Fraction
    residues: tuple[Fraction, ...]
    rho: tuple[int, ...]
    generic_slopes: tuple[int, ...]
    conifold_slopes: tuple[int, ...]
    slope_zero_multiplicity: int


def fractional_part(x: Fraction) -> Fraction:
    return x - (x.numerator // x.denominator)


def hodge_data(p: int, j: int) -> HodgeData:
    """AS/Fedorov data for Hyp(1^4; A,A,A^-1,A^-1) at t=1.

    The upper exponents are alpha=(0,0,0,0).  The lower fractional exponents
    are beta=(a,a,1-a,1-a).  Since every beta is positive, the sorted-height
    statistic is rho(k)=#{i: alpha_i < beta_k}-k=4-k.  The generic circuit has
    slopes 0,1,2,3.  At the balanced conifold t=1, the middle-extension stalk
    is ker(N) and the Gr^W_4 slope-2 slot is removed.
    """
    if not (p >= 5 and 1 <= j <= p - 2):
        raise ValueError("hodge_data requires p>=5 and 1<=j<=p-2")
    a = Fraction(j, p - 1)
    residues = tuple(sorted((a, a, fractional_part(-a), fractional_part(-a))))
    rho = tuple(4 - k for k in range(1, 5))
    generic = tuple(sorted(rho))
    conifold = (0, 1, 3)
    return HodgeData(a, residues, rho, generic, conifold, 1)


def direct_tlambda_polytope_volume() -> Fraction:
    """Normalized volume of conv(0, Newt(Lambda) x {1}).

    The x-cross-section of Newt(Lambda) has area 4 on [-1,0] and (2-x)^2 on
    [0,1].  Thus 3! vol(Newt(Lambda))=38, also the normalized 4-volume of the
    pyramid for t*Lambda.  This is the *wrong* AS rank because t*Lambda is
    Newton-degenerate.
    """
    volume = Fraction(4) + (
        Fraction(4) - Fraction(2) + Fraction(1, 3)
    )  # integral_0^1 (2-x)^2 dx = 7/3
    assert volume == Fraction(19, 3)
    return 6 * volume


def tlambda_degenerate_witness() -> tuple[int, int, int, int]:
    """A toric point where t*Lambda and every logarithmic derivative vanish.

    At x=y=-1 two distinct simple factors of Lambda vanish.  Every first
    derivative of their product therefore vanishes.  z=t=1 avoids the toric
    boundary.
    """
    return (-1, -1, 1, 1)


# ---------------------------------------------------------------------------
# Part 2a: the extension-trace reduction and the cubic Newton factor


def extension_exponent(p: int, j: int, r: int) -> int:
    """Exponent j(1+p+...+p^(r-1)) of Lambda after reducing Teichmüller."""
    return j * (p**r - 1) // (p - 1)


def extension_trace_mod_p(p: int, j: int, r: int) -> int:
    """S_j^(r) modulo p, without enumerating F_{p^r}.

    Torus orthogonality gives CT Lambda^n with n=j(1+...+p^(r-1)); its base-p
    digits are r copies of j.  The Apéry Lucas congruence therefore gives
    S_r == b_n == b_j^r (mod p).  This is an extension-field trace statement,
    not a count of unit Mellin summands.
    """
    if not (1 <= j <= p - 2 and r >= 1):
        raise ValueError("invalid (p,j,r)")
    b = apery_residues(p)[j]
    return pow(b, r, p)


def newton_coefficients_mod_p(p: int, j: int) -> tuple[int, int, int]:
    """e_1,e_2,e_3 mod p from S_1,S_2,S_3 and Newton identities."""
    s1, s2, s3 = (extension_trace_mod_p(p, j, r) for r in (1, 2, 3))
    e1 = s1
    e2 = (s1 * s1 - s2) * pow(2, -1, p) % p
    e3 = (s1**3 - 3 * s1 * s2 + 2 * s3) * pow(6, -1, p) % p
    return e1, e2, e3


def slope_zero_count(p: int, j: int) -> int:
    """Slope-zero multiplicity of the reduced rank-three factor."""
    return int(apery_residues(p)[j] != 0)


def frobenius_slopes(trace_parameter_valuation: int) -> tuple[Fraction, ...]:
    """Slopes of (X-p)(X^2-aX+p^3) from v_p(a)."""
    m = trace_parameter_valuation
    if m < 0:
        raise ValueError("valuation must be nonnegative")
    if m == 0:
        pair = (Fraction(0), Fraction(3))
    elif m == 1:
        pair = (Fraction(1), Fraction(2))
    else:
        pair = (Fraction(3, 2), Fraction(3, 2))
    return tuple(sorted((Fraction(1),) + pair))


def factor_coefficients(p: int, s1):
    """Coefficients e1,e2,e3 of the cubic factor from its exact trace."""
    return s1, p * s1 + p**3 - p**2, p**4


def predicted_power_traces(p: int, s1: int, rmax: int) -> list[int]:
    """Power traces of (X-p)(X^2-(s1-p)X+p^3)."""
    a = s1 - p
    pair = [2, a]
    for r in range(2, rmax + 1):
        pair.append(a * pair[-1] - p**3 * pair[-2])
    return [p**r + pair[r] for r in range(1, rmax + 1)]


# ---------------------------------------------------------------------------
# Part 2b: base-field A1 in Z/p^M and valuations at exceptional points


def teichmuller_lift(value: int, p: int, precision: int) -> int:
    """The Teichmüller lift modulo p^precision of value in F_p."""
    modulus = p**precision
    x = value % modulus
    for _ in range(precision + 1):
        x = pow(x, p, modulus)
    assert x % p == value % p and pow(x, p, modulus) == x
    return x


def base_a1_trace_padic(p: int, j: int, precision: int = 4) -> int:
    """The exact A1 skeleton over F_p, evaluated modulo p^precision."""
    if not 1 <= j <= p - 2:
        raise ValueError("base_a1_trace_padic requires nontrivial j")
    h = p - 1
    modulus = p**precision
    roots = [0] + [teichmuller_lift(t, p, precision) for t in range(1, p)]
    chars = [[0] * p for _ in range(h)]
    for t in range(1, p):
        value = 1
        for exponent in range(h):
            chars[exponent][t] = value
            value = value * roots[t] % modulus

    cache: dict[tuple[int, int], int] = {}

    def jacobi(a: int, b: int) -> int:
        key = (a % h, b % h)
        if key not in cache:
            aa, bb = key
            cache[key] = sum(
                chars[aa][t] * chars[bb][(1 - t) % p] for t in range(p)
            ) % modulus
        return cache[key]

    total = 0
    for k in range(h):
        first = jacobi(-k, j)
        second = jacobi(-k, j + k)
        total = (total + first * first * second * second) % modulus
    return -total * pow(h, -1, modulus) % modulus


def valuation_mod(value: int, p: int, precision: int) -> int:
    value %= p**precision
    if value == 0:
        return precision
    answer = 0
    while value % p == 0:
        answer += 1
        value //= p
    return answer


def trace_parameter_valuation(p: int, j: int, precision: int = 4) -> int:
    """v_p(a_{p,j}) for a=S_1-p, to the requested precision."""
    b = apery_residues(p)[j]
    if b:
        return 0
    s1 = base_a1_trace_padic(p, j, precision)
    return valuation_mod(s1 - p, p, precision)


# ---------------------------------------------------------------------------
# Small finite fields and the full extension-field engines


def _digits(value: int, p: int, degree: int) -> list[int]:
    answer = []
    for _ in range(degree):
        answer.append(value % p)
        value //= p
    return answer


def _encode(coefficients: list[int], p: int) -> int:
    answer = 0
    for coefficient in reversed(coefficients):
        answer = answer * p + coefficient % p
    return answer


class FiniteField:
    """A small F_{p^r}, sufficient for the exact probes with r<=3."""

    def __init__(self, p: int, degree: int):
        if degree not in (1, 2, 3):
            raise ValueError("the small exact engine supports degrees 1,2,3")
        self.p = p
        self.degree = degree
        self.q = p**degree
        self.modulus = self._find_irreducible()
        self.generator = self._find_generator()
        self.exp = [0] * (self.q - 1)
        self.log = [-1] * self.q
        value = 1
        for n in range(self.q - 1):
            self.exp[n] = value
            self.log[value] = n
            value = self._raw_mul(value, self.generator)
        assert value == 1 and all(x >= 0 for x in self.log[1:])

    def _eval_poly(self, coefficients: tuple[int, ...], x: int) -> int:
        value = 1
        answer = 0
        for coefficient in coefficients:
            answer = (answer + coefficient * value) % self.p
            value = value * x % self.p
        return answer

    def _find_irreducible(self) -> tuple[int, ...]:
        if self.degree == 1:
            return (0,)
        for coefficients in itertools.product(range(self.p), repeat=self.degree):
            if coefficients[0] == 0:
                continue
            if all(self._eval_poly(coefficients + (1,), x) for x in range(self.p)):
                # For degrees 2 and 3, no root is equivalent to irreducibility.
                return tuple(coefficients)
        raise AssertionError("no irreducible polynomial found")

    def add(self, left: int, right: int) -> int:
        a = _digits(left, self.p, self.degree)
        b = _digits(right, self.p, self.degree)
        return _encode([(x + y) % self.p for x, y in zip(a, b)], self.p)

    def neg(self, value: int) -> int:
        return _encode([(-x) % self.p for x in _digits(value, self.p, self.degree)], self.p)

    def sub(self, left: int, right: int) -> int:
        return self.add(left, self.neg(right))

    def _raw_mul(self, left: int, right: int) -> int:
        if self.degree == 1:
            return left * right % self.p
        a = _digits(left, self.p, self.degree)
        b = _digits(right, self.p, self.degree)
        product = [0] * (2 * self.degree - 1)
        for i, x in enumerate(a):
            for k, y in enumerate(b):
                product[i + k] = (product[i + k] + x * y) % self.p
        for power in range(2 * self.degree - 2, self.degree - 1, -1):
            coefficient = product[power]
            if coefficient:
                for k, relation in enumerate(self.modulus):
                    product[power - self.degree + k] -= coefficient * relation
                    product[power - self.degree + k] %= self.p
        return _encode(product[: self.degree], self.p)

    def _raw_pow(self, value: int, exponent: int) -> int:
        answer = 1
        while exponent:
            if exponent & 1:
                answer = self._raw_mul(answer, value)
            value = self._raw_mul(value, value)
            exponent >>= 1
        return answer

    @staticmethod
    def _prime_factors(value: int) -> list[int]:
        factors = []
        divisor = 2
        while divisor * divisor <= value:
            if value % divisor == 0:
                factors.append(divisor)
                while value % divisor == 0:
                    value //= divisor
            divisor += 1
        if value > 1:
            factors.append(value)
        return factors

    def _find_generator(self) -> int:
        order = self.q - 1
        factors = self._prime_factors(order)
        for candidate in range(2, self.q):
            if all(self._raw_pow(candidate, order // ell) != 1 for ell in factors):
                return candidate
        raise AssertionError("no primitive generator found")

    def mul(self, left: int, right: int) -> int:
        if left == 0 or right == 0:
            return 0
        return self.exp[(self.log[left] + self.log[right]) % (self.q - 1)]

    def inv(self, value: int) -> int:
        if value == 0:
            raise ZeroDivisionError
        return self.exp[(-self.log[value]) % (self.q - 1)]

    def div(self, left: int, right: int) -> int:
        return self.mul(left, self.inv(right))

    def pow(self, value: int, exponent: int) -> int:
        if value == 0:
            return int(exponent == 0)
        return self.exp[(self.log[value] * exponent) % (self.q - 1)]

    def prime_subfield_value(self, value: int) -> int:
        coefficients = _digits(value, self.p, self.degree)
        assert all(coefficient == 0 for coefficient in coefficients[1:])
        return coefficients[0]


_FIELD_CACHE: dict[tuple[int, int], FiniteField] = {}


def small_field(p: int, r: int) -> FiniteField:
    key = (p, r)
    if key not in _FIELD_CACHE:
        _FIELD_CACHE[key] = FiniteField(p, r)
    return _FIELD_CACHE[key]


def exact_extension_trace_vector(p: int, j: int, r: int) -> list[int]:
    """Exact S_r as coefficients of 1,zeta_h,...,zeta_h^(h-1).

    This is an independent O(q^2) count in the arrangement coordinates.  It
    avoids floating point and is used to certify the small extension probes.
    """
    field = small_field(p, r)
    q = field.q
    h = p - 1
    nonzero = field.exp
    one = 1

    def character_exponent(value: int) -> int | None:
        if value == 0:
            return None
        return j * field.log[value] % h

    kernel: list[list[int]] = [[0] * h for _ in range(q)]
    for c in nonzero:
        counts = kernel[c]
        for w in nonzero:
            value = field.div(
                field.mul(field.sub(one, field.mul(c, w)), field.sub(w, one)),
                w,
            )
            exponent = character_exponent(value)
            if exponent is not None:
                counts[exponent] += 1

    total = [0] * h
    for u in nonzero:
        if u == one:
            continue
        for v in nonzero:
            if v == one:
                continue
            uv = field.mul(u, v)
            outside = field.div(uv, field.mul(field.sub(u, one), field.sub(v, one)))
            exponent = character_exponent(outside)
            assert exponent is not None
            for e, count in enumerate(kernel[uv]):
                total[(exponent + e) % h] -= count  # S_r has the outer minus sign
    return total


def evaluate_trace_vector_padic(
    coefficients: list[int], p: int, root: int, precision: int
) -> int:
    modulus = p**precision
    return sum(
        coefficient * pow(root, exponent, modulus)
        for exponent, coefficient in enumerate(coefficients)
    ) % modulus


def cyclic_square(values: list[int], modulus: int) -> list[int]:
    length = len(values)
    answer = [0] * length
    for i, left in enumerate(values):
        if left == 0:
            continue
        for k, right in enumerate(values):
            if right:
                index = (i + k) % length
                answer[index] = (answer[index] + left * right) % modulus
    return answer


def basechanged_a1_trace_padic(
    p: int, j: int, r: int, precision: int = 5
) -> int:
    """Full q-1 A1 skeleton over F_{p^r}, via exact cyclic convolution.

    Only p-1 of the q-1 Mellin characters are Hasse--Davenport lifts from
    F_p; this function keeps every extension character.  Its convolution is
    algebraically identical to the full Jacobi sum.
    """
    field = small_field(p, r)
    q = field.q
    length = q - 1
    h = p - 1
    modulus = p**precision
    base_generator = field.pow(field.generator, length // h)
    residue = field.prime_subfield_value(base_generator)
    zeta = teichmuller_lift(residue, p, precision)
    zeta_powers = [pow(zeta, exponent, modulus) for exponent in range(h)]

    first: list[int] = []
    second: list[int] = []
    for n, x in enumerate(field.exp):
        one_minus = field.sub(1, x)
        if one_minus == 0:
            first.append(0)
        else:
            first.append(zeta_powers[j * (field.log[one_minus] % h) % h])

        one_plus = field.add(1, x)
        if one_plus == 0:
            second.append(0)
        else:
            exponent = (n - field.log[one_plus]) % h
            second.append(zeta_powers[j * exponent % h])

    first_square = cyclic_square(first, modulus)
    second_square = cyclic_square(second, modulus)
    return -sum(
        left * right for left, right in zip(first_square, second_square)
    ) % modulus


def _newton_slopes_from_valuations(
    coefficient_valuations: tuple[int, ...],
) -> tuple[Fraction, ...]:
    """Lower-hull slopes for 1-e1*T+e2*T^2-... from v_p(e_i)."""
    points = [(0, 0)] + [
        (index, value) for index, value in enumerate(coefficient_valuations, 1)
    ]
    hull: list[tuple[int, int]] = []
    for point in points:
        while len(hull) >= 2:
            x0, y0 = hull[-2]
            x1, y1 = hull[-1]
            x2, y2 = point
            if (y1 - y0) * (x2 - x1) >= (y2 - y1) * (x1 - x0):
                hull.pop()
            else:
                break
        hull.append(point)
    slopes: list[Fraction] = []
    for (x0, y0), (x1, y1) in zip(hull, hull[1:]):
        slopes.extend([Fraction(y1 - y0, x1 - x0)] * (x1 - x0))
    return tuple(slopes)


def sage_full_extension_scan(bound: int = 59, precision: int = 6) -> dict[str, object]:
    """Independently scan S_1,S_2,S_3 using every F_{p^r} A1 character.

    Run this optional production path under ``sage -python``.  For each r it
    constructs F_{p^r}, retains all q-1 discrete-log positions, and squares the
    two cyclic-convolution polynomials over Z/p^precision.  Neither the Apéry
    recurrence nor the cubic factor is used to construct the traces.  The
    factor identities are checked only after Newton's identities have recovered
    the three coefficients.
    """
    try:
        import sage.all as sage  # type: ignore[import-not-found]
    except ImportError as error:
        raise RuntimeError(
            "the full range scan requires: sage -python lemmaBprime_explore.py "
            "--full-small-scan"
        ) from error

    def extension_exponents(p: int, r: int):
        q = p**r
        length = q - 1
        h = p - 1
        field = sage.GF(q, "a")
        generator = field.multiplicative_generator()
        logs = {}
        value = field.one()
        for n in range(length):
            logs[value] = n
            value *= generator
        first = []
        second = []
        value = field.one()
        for n in range(length):
            difference = field.one() - value
            first.append(None if difference == 0 else logs[difference] % h)
            total = field.one() + value
            second.append(None if total == 0 else (n - logs[total]) % h)
            value *= generator
        base_generator = int(generator ** (length // h))
        return first, second, base_generator

    def cyclic_square_sage(values, polynomial_ring, ring, length: int):
        square = polynomial_ring(values) ** 2
        answer = [ring(0)] * length
        for index, coefficient in enumerate(square):
            answer[index % length] += coefficient
        return answer

    def traces_for_degree(p: int, r: int) -> dict[int, int]:
        h = p - 1
        q = p**r
        length = q - 1
        modulus = p**precision
        ring = sage.Integers(modulus)
        polynomial_ring = sage.PolynomialRing(ring, "X")
        first_exponents, second_exponents, residue = extension_exponents(p, r)
        zeta = ring(teichmuller_lift(residue, p, precision))
        zeta_powers = [zeta**exponent for exponent in range(h)]
        answer = {}
        for j in range(1, h):
            first = [
                ring(0) if exponent is None else zeta_powers[j * exponent % h]
                for exponent in first_exponents
            ]
            second = [
                ring(0) if exponent is None else zeta_powers[j * exponent % h]
                for exponent in second_exponents
            ]
            left = cyclic_square_sage(first, polynomial_ring, ring, length)
            right = cyclic_square_sage(second, polynomial_ring, ring, length)
            answer[j] = int(-sum(
                (left[index] * right[index] for index in range(length)), ring(0)
            ))
        return answer

    records = []
    per_prime = {}
    for p in (prime for prime in primes_up_to(bound) if prime >= 5):
        modulus = p**precision
        traces = [traces_for_degree(p, r) for r in (1, 2, 3)]
        inverse_two = pow(2, -1, modulus)
        inverse_six = pow(6, -1, modulus)
        prime_counts: dict[tuple[Fraction, ...], int] = {}
        for j in range(1, p - 1):
            s1, s2, s3 = (traces[index][j] for index in range(3))
            e1 = s1 % modulus
            e2 = (s1 * s1 - s2) * inverse_two % modulus
            e3 = (s1**3 - 3 * s1 * s2 + 2 * s3) * inverse_six % modulus
            # These are post-computation checks, not inputs to the traces.
            if e2 != (p * e1 + p**3 - p**2) % modulus or e3 != p**4 % modulus:
                raise AssertionError((p, j, "cubic factor", e1, e2, e3))
            valuations = tuple(
                valuation_mod(coefficient, p, precision)
                for coefficient in (e1, e2, e3)
            )
            slopes = _newton_slopes_from_valuations(valuations)
            prime_counts[slopes] = prime_counts.get(slopes, 0) + 1
            records.append((p, j, s1 % p, valuations, slopes))
        per_prime[p] = prime_counts
    return {
        "bound": bound,
        "precision": precision,
        "records": records,
        "pairs": len(records),
        "ordinary": sum(slopes == (Fraction(0), Fraction(1), Fraction(3))
                        for _, _, _, _, slopes in records),
        "nonordinary": sum(slopes != (Fraction(0), Fraction(1), Fraction(3))
                           for _, _, _, _, slopes in records),
        "per_prime": per_prime,
    }


def reduce_zeta4(coefficients: list[int]) -> tuple[int, int]:
    if len(coefficients) != 4:
        raise ValueError
    return coefficients[0] - coefficients[2], coefficients[1] - coefficients[3]


def reduce_zeta6(coefficients: list[int]) -> tuple[int, int]:
    if len(coefficients) != 6:
        raise ValueError
    c = coefficients
    return c[0] - c[2] - c[3] + c[5], c[1] + c[2] - c[4] - c[5]


def _polynomial_div_exact(numerator: list[int], denominator: list[int]) -> list[int]:
    """Exact division in Z[X], coefficient lists in ascending order."""
    work = numerator[:]
    while len(work) > 1 and work[-1] == 0:
        work.pop()
    quotient = [0] * max(1, len(work) - len(denominator) + 1)
    while len(work) >= len(denominator):
        if work[-1] % denominator[-1]:
            raise AssertionError("nonexact polynomial division")
        coefficient = work[-1] // denominator[-1]
        shift = len(work) - len(denominator)
        quotient[shift] = coefficient
        for i, value in enumerate(denominator):
            work[shift + i] -= coefficient * value
        while len(work) > 1 and work[-1] == 0:
            work.pop()
    if any(work):
        raise AssertionError((numerator, denominator, work))
    return quotient


_CYCLOTOMIC_CACHE: dict[int, list[int]] = {1: [-1, 1]}


def cyclotomic_polynomial(order: int) -> list[int]:
    """Phi_order in Z[X], coefficients in ascending order."""
    if order not in _CYCLOTOMIC_CACHE:
        polynomial = [-1] + [0] * (order - 1) + [1]
        for divisor in range(1, order):
            if order % divisor == 0:
                polynomial = _polynomial_div_exact(
                    polynomial, cyclotomic_polynomial(divisor)
                )
        _CYCLOTOMIC_CACHE[order] = polynomial
    return _CYCLOTOMIC_CACHE[order][:]


def reduce_cyclotomic(coefficients: list[int]) -> list[int]:
    """Canonical coefficient vector in Z[X]/Phi_h for h=len(coefficients)."""
    phi = cyclotomic_polynomial(len(coefficients))
    degree = len(phi) - 1
    work = coefficients[:] + [0]
    for power in range(len(work) - 1, degree - 1, -1):
        coefficient = work[power]
        if coefficient:
            shift = power - degree
            for i, value in enumerate(phi):
                work[shift + i] -= coefficient * value
    if any(work[degree:]):
        raise AssertionError((coefficients, phi, work))
    return work[:degree]


EXACT_TRACE_TABLE: dict[tuple[int, int], tuple[int, ...]] = {
    (5, 1): (-5, -125, 2875, -8125),
    (5, 2): (3, -221, 867, 29891),
    (7, 1): (-2, -556, 8875, 133128),
    (7, 2): (-18, -12, 10443, -229176),
    (7, 3): (31, -61, -10529, -220797),
    (11, 5): (-33, -605, 91839),
    (13, 1): (5, -4161, 54413),
}


def exact_probe(p: int, j: int, rmax: int = 3, precision: int = 5) -> list[int]:
    """Run both exact-count and full-A1 engines and return rational traces.

    The built-in probe pairs have rational traces.  Each coefficient vector is
    reduced exactly in Q(zeta_{p-1}), then independently compared with the
    p-adic full extension A1 skeleton.
    """
    traces = []
    for r in range(1, rmax + 1):
        vector = exact_extension_trace_vector(p, j, r)
        reduced = reduce_cyclotomic(vector)
        if any(reduced[1:]):
            raise AssertionError((p, j, r, "non-rational trace", reduced))
        trace = reduced[0]
        a1 = basechanged_a1_trace_padic(p, j, r, precision)
        field = small_field(p, r)
        base_generator = field.pow(field.generator, (field.q - 1) // (p - 1))
        root = teichmuller_lift(field.prime_subfield_value(base_generator), p, precision)
        direct_padic = evaluate_trace_vector_padic(vector, p, root, precision)
        modulus = p**precision
        if a1 != direct_padic or a1 != trace % modulus:
            raise AssertionError((p, j, r, a1, direct_padic, trace % modulus))
        traces.append(trace)
    return traces


# ---------------------------------------------------------------------------
# Part 3: Delta critical values and exact branch tests


def quadratic_pair_mul(left: tuple[int, int], right: tuple[int, int]) -> tuple[int, int]:
    """Multiply a+b*sqrt(2) pairs."""
    a, b = left
    c, d = right
    return a * c + 2 * b * d, a * d + b * c


def quadratic_pair_pow(value: tuple[int, int], exponent: int) -> tuple[int, int]:
    answer = (1, 0)
    while exponent:
        if exponent & 1:
            answer = quadratic_pair_mul(answer, value)
        value = quadratic_pair_mul(value, value)
        exponent >>= 1
    return answer


def delta_critical_values() -> tuple[tuple[int, int], tuple[int, int]]:
    values = (quadratic_pair_pow((1, 1), 4), quadratic_pair_pow((1, -1), 4))
    assert values == ((17, 12), (17, -12))
    for a, b in values:
        square = quadratic_pair_mul((a, b), (a, b))
        # (a+b*s)^2 - 34(a+b*s) + 1 = 0, s^2=2.
        assert (square[0] - 34 * a + 1, square[1] - 34 * b) == (0, 0)
    return values


def split_prime_candidate_test_p7() -> dict[str, object]:
    """Exact p=7 rejection of Delta, Jacobi-field, and square candidates."""
    # Q_j=X^2-a_j X+7^3 for j=1,2,3; Delta roots modulo 7 are 2 and 4.
    parameters = {1: -9, 2: -25, 3: 24}
    delta_residues = {1: {2, 4}, 2: {2, 4}, 3: {1}}
    unit_residues = {j: a % 7 for j, a in parameters.items()}
    assert all(unit_residues[j] not in delta_residues[j] for j in parameters)
    # For j=1,2 the candidate minimal polynomial is Phi_3=X^2+X+1.
    remainders = {
        j: (-1 - parameters[j], 7**3 - 1) for j in (1, 2)
    }  # coefficient of X, constant
    assert remainders == {1: (8, 342), 2: (24, 342)}
    assert 1 - parameters[3] + 7**3 == 320
    assert all(pow(residue, 3, 7) == -1 % 7 for residue in unit_residues.values())

    def squarefree_part(value: int) -> int:
        sign = -1 if value < 0 else 1
        work = abs(value)
        answer = 1
        divisor = 2
        while divisor * divisor <= work:
            parity = 0
            while work % divisor == 0:
                work //= divisor
                parity ^= 1
            if parity:
                answer *= divisor
            divisor += 1
        if work > 1:
            answer *= work
        return sign * answer

    raw_discriminants = {
        j: parameter * parameter - 4 * 7**3
        for j, parameter in parameters.items()
    }
    discriminant_fields = tuple(
        squarefree_part(raw_discriminants[j]) for j in sorted(parameters)
    )
    assert discriminant_fields == (-1291, -83, -199)
    # Q(sqrt(d))=Q(sqrt(-3)) exactly when the squarefree parts agree.
    assert all(value != -3 for value in discriminant_fields)

    # In Q_7 a square has even valuation, and a unit square has square residue.
    # The three unit roots have nonsquare residues; the Tate and complementary
    # roots have odd valuations 1 and 3.
    eigenvalue_square_obstructions = {
        "unit_residue_nonsquare": all(
            pow(residue, 3, 7) == -1 % 7 for residue in unit_residues.values()
        ),
        "other_valuations_odd": (1 % 2 == 1 and 3 % 2 == 1),
    }
    all_local_eigenvalues_nonsquare = all(eigenvalue_square_obstructions.values())
    return {
        "delta_roots_mod_7": (2, 4),
        "unit_residues": unit_residues,
        "phi3_remainders": remainders,
        "Q3_at_1": 320,
        "unit_root_fields": discriminant_fields,
        "square_obstructions": eigenvalue_square_obstructions,
        "all_local_eigenvalues_nonsquare": all_local_eigenvalues_nonsquare,
    }


# ---------------------------------------------------------------------------
# Required ranges and command-line report


def small_range_report(bound: int = 59) -> dict[str, object]:
    records = []
    ordinary = 0
    nonordinary = 0
    exceptional: list[tuple[int, int]] = []
    for p in primes_up_to(bound):
        if p < 5:
            continue
        values = apery_residues(p)
        for j in range(1, p - 1):
            m = trace_parameter_valuation(p, j) if values[j] == 0 else 0
            slopes = frobenius_slopes(m)
            record = (p, j, values[j], slope_zero_count(p, j), slopes)
            records.append(record)
            if slopes == (Fraction(0), Fraction(1), Fraction(3)):
                ordinary += 1
            else:
                nonordinary += 1
                exceptional.append((p, j))
    return {
        "bound": bound,
        "records": records,
        "pairs": len(records),
        "ordinary": ordinary,
        "nonordinary": nonordinary,
        "exceptional": exceptional,
    }


def stretch_report(bound: int = 300, seed: int = 3_202_026) -> dict[str, object]:
    rng = random.Random(seed)
    zero_records = []
    random_records = []
    for p in primes_up_to(bound):
        if p < 5:
            continue
        values = apery_residues(p)
        zeros = [j for j in range(1, p - 1) if values[j] == 0]
        all_twists = list(range(1, p - 1))
        for j in zeros:
            m = trace_parameter_valuation(p, j)
            zero_records.append((p, j, m, frobenius_slopes(m)))
        for j in sorted(rng.sample(all_twists, min(20, len(all_twists)))):
            m = trace_parameter_valuation(p, j) if values[j] == 0 else 0
            random_records.append((p, j, m, frobenius_slopes(m)))
    return {
        "bound": bound,
        "seed": seed,
        "zero_records": zero_records,
        "random_records": random_records,
        "zeros": len(zero_records),
        "random": len(random_records),
    }


def format_slopes(slopes: tuple[Fraction, ...]) -> str:
    return "(" + ",".join(str(slope) for slope in slopes) + ")"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--small-bound", type=int, default=59)
    parser.add_argument("--stretch-bound", type=int, default=300)
    parser.add_argument("--seed", type=int, default=3_202_026)
    parser.add_argument(
        "--exact-probes", action="store_true", help="run portable F_{p^r} probes"
    )
    parser.add_argument(
        "--full-small-scan",
        action="store_true",
        help="under sage -python, independently scan all r=1,2,3 traces",
    )
    args = parser.parse_args()

    sample = hodge_data(11, 5)
    print(
        "HODGE",
        f"a={sample.a}",
        f"beta={sample.residues}",
        f"rho={sample.rho}",
        f"generic={sample.generic_slopes}",
        f"t=1={sample.conifold_slopes}",
        f"h0={sample.slope_zero_multiplicity}",
    )
    print(
        "DIRECT_AS_OBSTRUCTION",
        f"normalized_volume={direct_tlambda_polytope_volume()}",
        f"degenerate_witness={tlambda_degenerate_witness()}",
    )

    small = small_range_report(args.small_bound)
    print(
        "SMALL_NEWTON",
        f"bound={small['bound']}",
        f"pairs={small['pairs']}",
        f"HP=NP={small['ordinary']}",
        f"nonordinary={small['nonordinary']}",
        "slope0_ge2=0",
    )
    print("SMALL_EXCEPTIONAL", small["exceptional"])

    stretch = stretch_report(args.stretch_bound, args.seed)
    print(
        "STRETCH_NEWTON",
        f"bound={stretch['bound']}",
        f"zeros={stretch['zeros']}",
        f"random={stretch['random']}",
        f"random_zero_hits={sum(m > 0 for _, _, m, _ in stretch['random_records'])}",
        "zero_slopes=(1,1,2)",
        "slope0_ge2=0",
    )

    print("DELTA", delta_critical_values(), "critical_values_not_eigenbranches")
    print("P7_BRANCH_TEST", split_prime_candidate_test_p7())
    print("EXACT_FACTOR_TABLE")
    for (p, j), traces in EXACT_TRACE_TABLE.items():
        predicted = predicted_power_traces(p, traces[0], len(traces))
        if list(traces) != predicted:
            raise AssertionError((p, j, traces, predicted))
        e1, e2, e3 = factor_coefficients(p, traces[0])
        m = 0 if traces[0] % p else valuation_mod(traces[0] - p, p, 8)
        print(
            f"  p={p} j={j} S={traces} e=({e1},{e2},{e3}) "
            f"slopes={format_slopes(frobenius_slopes(m))}"
        )

    if args.exact_probes:
        for p, j, rmax in (
            (5, 1, 3),
            (5, 2, 3),
            (7, 1, 3),
            (7, 2, 3),
            (7, 3, 3),
            (11, 5, 3),
            (13, 1, 1),
        ):
            traces = exact_probe(p, j, rmax)
            expected = list(EXACT_TRACE_TABLE[p, j][:rmax])
            if traces != expected:
                raise AssertionError((p, j, traces, expected))
            print(f"EXACT_EXTENSION_A1 p={p} j={j} traces={traces} PASS")

    if args.full_small_scan:
        full = sage_full_extension_scan(args.small_bound, 6)
        print(
            "FULL_EXTENSION_SCAN",
            f"bound={full['bound']}",
            f"precision={full['precision']}",
            f"pairs={full['pairs']}",
            f"ordinary={full['ordinary']}",
            f"nonordinary={full['nonordinary']}",
        )

    print("VERDICT G1")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
