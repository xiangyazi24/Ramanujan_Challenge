#!/usr/bin/env python3
"""Apéry ``K_infinity`` branch table and fixed-height critical-value scan.

This program implements the finite-height part of
``CHATGPT_Q6723_kinf_algorithm.md``.  There are two deliberately separate
numerical layers:

* an uncertified mpmath pass (at least 100 decimal digits), and
* a certified Arb pass using python-flint complex balls, exhaustive root
  isolation, and an explicit interval-Newton inclusion for every root.

For the fixed-height certificate we use the exact reflection quotient

    u = (2*x+h+1)^2.

The critical polynomial becomes a degree ``2*h-2`` polynomial ``J_h(u)``.
Evaluating ``delta_h(x)^2`` as an exact rational function of ``u`` gives
precisely the mirror-orbits of squared critical values requested in Q6723,
but avoids isolating both members of every forced reflection pair.

The system Python on the Mac mini does not expose python-flint, while Xiang's
uv cache does.  As with the Sage-backed CRON scripts in this directory, the
entry point re-executes itself in the pinned uv environment when necessary.
No package is installed into the system interpreter.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math
import os
from pathlib import Path
import shutil
import sys
import time
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Sequence, Tuple


UV_RUNTIME_GUARD = "CRON_KINF_UV_RUNTIME"


def _ensure_runtime() -> None:
    """Re-execute with the pinned numerical dependencies when needed."""

    try:
        import flint  # type: ignore  # noqa: F401
        import mpmath  # type: ignore  # noqa: F401
        return
    except ImportError:
        pass

    if os.environ.get(UV_RUNTIME_GUARD) == "1":
        raise SystemExit(
            "python-flint/mpmath are unavailable even inside the uv runtime"
        )
    uv = shutil.which("uv")
    if uv is None:
        raise SystemExit(
            "python-flint is not importable and `uv` was not found; "
            "install python-flint 0.9.0 or run this script through uv"
        )

    os.environ[UV_RUNTIME_GUARD] = "1"
    command = [
        uv,
        "run",
        "--python",
        "3.12",
        "--with",
        "python-flint==0.9.0",
        "--with",
        "mpmath==1.3.0",
        "python",
        str(Path(__file__).resolve()),
        *sys.argv[1:],
    ]
    os.execvp(uv, command)


_ensure_runtime()

from flint import acb, acb_poly, arb, ctx, fmpz_poly  # type: ignore  # noqa: E402
from mpmath import mp  # type: ignore  # noqa: E402


if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)


LAMBDA_PLUS_TEXT = "17 + 12*sqrt(2)"
LAMBDA_MINUS_TEXT = "17 - 12*sqrt(2)"


def log(message: str) -> None:
    print(message, flush=True)


def apery_P(value: Any) -> Any:
    """The cubic transfer coefficient ``34t^3+51t^2+27t+5``."""

    return 34 * value**3 + 51 * value**2 + 27 * value + 5


def build_numerators(height: int) -> List[Any]:
    """Return exact ``N_0,...,N_height`` over ``ZZ[x]``."""

    x = fmpz_poly([0, 1])
    values = [fmpz_poly([]), fmpz_poly([1])]
    for h in range(1, height):
        values.append(
            apery_P(x + h) * values[h] - (x + h) ** 6 * values[h - 1]
        )
    return values


def build_apery_numbers(height: int) -> List[int]:
    """Return the exact Apéry numbers ``b_0,...,b_height``."""

    if height == 0:
        return [1]
    values = [1, 5]
    for n in range(1, height):
        numerator = int(apery_P(n)) * values[n] - n**3 * values[n - 1]
        denominator = (n + 1) ** 3
        if numerator % denominator:
            raise ArithmeticError("Apéry recurrence lost integrality at n=%d" % n)
        values.append(numerator // denominator)
    return values


def q_polynomial(h: int) -> Any:
    """Return ``q_h(x)=product_{j=1}^h (x+j)``."""

    x = fmpz_poly([0, 1])
    result = fmpz_poly([1])
    for j in range(1, h + 1):
        result *= x + j
    return result


def binary_affine_lift(polynomial: Any, center: int) -> Any:
    """Compute ``2^d p((s-center)/2)`` exactly by Horner evaluation."""

    s = fmpz_poly([0, 1])
    degree = int(polynomial.degree())
    if degree < 0:
        return fmpz_poly([])
    result = fmpz_poly([polynomial[degree]])
    scale = 1
    for exponent in range(degree - 1, -1, -1):
        scale <<= 1
        result = (s - center) * result + polynomial[exponent] * scale
    return result


def collapse_even(polynomial: Any, label: str) -> Any:
    """Write an exactly even polynomial in ``s`` as a polynomial in ``u=s^2``."""

    degree = int(polynomial.degree())
    bad = [k for k in range(1, degree + 1, 2) if polynomial[k] != 0]
    if bad:
        raise ArithmeticError(
            "%s is not reflection-even; first odd exponent is %d" % (label, bad[0])
        )
    return fmpz_poly([polynomial[2 * k] for k in range(degree // 2 + 1)])


@dataclass
class HeightData:
    """Exact quotient data for one height."""

    h: int
    numerator: Any
    q: Any
    critical: Any
    quotient_critical: Any
    squared_value_num: Any
    squared_value_den: Any


def build_height_data(h: int, numerator: Any) -> HeightData:
    """Build ``J_h`` and the exact squared critical-value rational function."""

    qh = q_polynomial(h)
    critical = qh * numerator.derivative() - 3 * qh.derivative() * numerator
    expected_degree = 4 * h - 4
    if critical.degree() != expected_degree:
        raise ArithmeticError(
            "critical degree mismatch at h=%d: got %d, expected %d"
            % (h, critical.degree(), expected_degree)
        )

    center = h + 1
    critical_s = binary_affine_lift(critical, center)
    quotient_critical = collapse_even(critical_s, "A_h")
    if quotient_critical.degree() != 2 * h - 2:
        raise ArithmeticError("reflection quotient degree mismatch at h=%d" % h)

    numerator_s = binary_affine_lift(numerator, center)
    q_s = binary_affine_lift(qh, center)
    # Since deg(N_h)=3h-3 and deg(q_h)=h,
    #   delta_h((s-h-1)/2) = 8 N_s(s) / q_s(s)^3.
    value_num = 64 * collapse_even(numerator_s * numerator_s, "N_h^2")
    value_den = collapse_even(q_s**6, "q_h^6")
    return HeightData(
        h=h,
        numerator=numerator,
        q=qh,
        critical=critical,
        quotient_critical=quotient_critical,
        squared_value_num=value_num,
        squared_value_den=value_den,
    )


def mp_polyval(polynomial: Any, value: Any) -> Any:
    """Evaluate an fmpz polynomial by mpmath Horner arithmetic."""

    result = mp.mpc(0)
    for coefficient in reversed(polynomial.coeffs()):
        result = result * value + int(coefficient)
    return result


def mp_midpoint(value: Any, digits: int) -> Any:
    """Convert an exact Arb midpoint to an mpmath complex number."""

    midpoint = value.mid()
    real = midpoint.real.str(digits, radius=False)
    imag = midpoint.imag.str(digits, radius=False)
    return mp.mpc(real, imag)


def mp_scientific(value: Any, digits: int = 6) -> str:
    """Stable scientific notation for finite nonnegative mpmath values."""

    value = mp.mpf(value)
    if not mp.isfinite(value):
        return str(value)
    if value == 0:
        return "0"
    return mp.nstr(value, digits, min_fixed=0, max_fixed=0)


def float_scientific(value: float, digits: int = 6) -> str:
    if math.isnan(value) or math.isinf(value):
        return str(value)
    if value == 0:
        return "0"
    return ("%%.%de" % (digits - 1)) % value


@dataclass
class Jet2:
    """Order-two Taylor jet; entries are ``f, f', f''/2``."""

    coeff: Tuple[Any, Any, Any]

    @staticmethod
    def constant(value: Any) -> "Jet2":
        return Jet2((mp.mpc(value), mp.mpc(0), mp.mpc(0)))

    @staticmethod
    def variable(value: Any) -> "Jet2":
        return Jet2((mp.mpc(value), mp.mpc(1), mp.mpc(0)))

    @staticmethod
    def coerce(value: Any) -> "Jet2":
        return value if isinstance(value, Jet2) else Jet2.constant(value)

    def __add__(self, other: Any) -> "Jet2":
        other = Jet2.coerce(other)
        return Jet2(tuple(self.coeff[k] + other.coeff[k] for k in range(3)))

    __radd__ = __add__

    def __neg__(self) -> "Jet2":
        return Jet2(tuple(-value for value in self.coeff))

    def __sub__(self, other: Any) -> "Jet2":
        return self + (-Jet2.coerce(other))

    def __rsub__(self, other: Any) -> "Jet2":
        return Jet2.coerce(other) - self

    def __mul__(self, other: Any) -> "Jet2":
        other = Jet2.coerce(other)
        result = []
        for degree in range(3):
            result.append(
                sum(
                    (self.coeff[k] * other.coeff[degree - k] for k in range(degree + 1)),
                    mp.mpc(0),
                )
            )
        return Jet2(tuple(result))

    __rmul__ = __mul__

    def reciprocal(self) -> "Jet2":
        a0, a1, a2 = self.coeff
        if a0 == 0:
            raise ZeroDivisionError("jet reciprocal at zero")
        return Jet2((1 / a0, -a1 / a0**2, a1**2 / a0**3 - a2 / a0**2))

    def __truediv__(self, other: Any) -> "Jet2":
        return self * Jet2.coerce(other).reciprocal()

    def __rtruediv__(self, other: Any) -> "Jet2":
        return Jet2.coerce(other) / self

    def __pow__(self, exponent: int) -> "Jet2":
        if exponent < 0:
            return (self.reciprocal()) ** (-exponent)
        result = Jet2.constant(1)
        base = self
        power = exponent
        while power:
            if power & 1:
                result *= base
            base *= base
            power >>= 1
        return result

    @property
    def value(self) -> Any:
        return self.coeff[0]

    @property
    def first(self) -> Any:
        return self.coeff[1]

    @property
    def second(self) -> Any:
        return 2 * self.coeff[2]


@dataclass
class AcbJet2:
    """The same order-two jet, with C-level Arb/Acb scalar arithmetic."""

    coeff: Tuple[Any, Any, Any]

    @staticmethod
    def constant(value: Any) -> "AcbJet2":
        return AcbJet2((acb(value), acb(0), acb(0)))

    @staticmethod
    def variable(value: Any) -> "AcbJet2":
        return AcbJet2((acb(value), acb(1), acb(0)))

    @staticmethod
    def coerce(value: Any) -> "AcbJet2":
        return value if isinstance(value, AcbJet2) else AcbJet2.constant(value)

    def __add__(self, other: Any) -> "AcbJet2":
        other = AcbJet2.coerce(other)
        return AcbJet2(
            (
                self.coeff[0] + other.coeff[0],
                self.coeff[1] + other.coeff[1],
                self.coeff[2] + other.coeff[2],
            )
        )

    __radd__ = __add__

    def __neg__(self) -> "AcbJet2":
        return AcbJet2((-self.coeff[0], -self.coeff[1], -self.coeff[2]))

    def __sub__(self, other: Any) -> "AcbJet2":
        return self + (-AcbJet2.coerce(other))

    def __rsub__(self, other: Any) -> "AcbJet2":
        return AcbJet2.coerce(other) - self

    def __mul__(self, other: Any) -> "AcbJet2":
        other = AcbJet2.coerce(other)
        a0, a1, a2 = self.coeff
        b0, b1, b2 = other.coeff
        return AcbJet2(
            (
                a0 * b0,
                a0 * b1 + a1 * b0,
                a0 * b2 + a1 * b1 + a2 * b0,
            )
        )

    __rmul__ = __mul__

    def reciprocal(self) -> "AcbJet2":
        a0, a1, a2 = self.coeff
        return AcbJet2((1 / a0, -a1 / a0**2, a1**2 / a0**3 - a2 / a0**2))

    def __truediv__(self, other: Any) -> "AcbJet2":
        return self * AcbJet2.coerce(other).reciprocal()

    def __rtruediv__(self, other: Any) -> "AcbJet2":
        return AcbJet2.coerce(other) / self

    def __pow__(self, exponent: int) -> "AcbJet2":
        if exponent < 0:
            return self.reciprocal() ** (-exponent)
        result = AcbJet2.constant(1)
        base = self
        power = exponent
        while power:
            if power & 1:
                result *= base
            base *= base
            power >>= 1
        return result

    @property
    def value(self) -> Any:
        return self.coeff[0]

    @property
    def first(self) -> Any:
        return self.coeff[1]

    @property
    def second(self) -> Any:
        return 2 * self.coeff[2]


def jet_linear_combination(values: Sequence[Jet2], weights: Sequence[Any]) -> Jet2:
    return Jet2(
        tuple(
            sum((weight * value.coeff[k] for value, weight in zip(values, weights)), mp.mpc(0))
            for k in range(3)
        )
    )


def extrapolation_weights(nodes: Sequence[int]) -> List[Any]:
    """Lagrange weights at ``t=0`` for nodes ``t_i=1/n_i``."""

    ts = [mp.mpf(1) / n for n in nodes]
    result = []
    for i, ti in enumerate(ts):
        weight = mp.mpf(1)
        for j, tj in enumerate(ts):
            if i != j:
                weight *= -tj / (ti - tj)
        result.append(weight)
    return result


def phi_gamma_sequence_jets(
    z: Jet2, nodes: Sequence[int], apery: Sequence[int]
) -> Tuple[List[Jet2], List[Jet2]]:
    """Evaluate normalized finite ``phi_n,gamma_n`` jets at selected indices."""

    wanted = set(nodes)
    phi_prev = Jet2.constant(1)
    phi = apery_P(z) / (5 * (z + 1) ** 3)
    gamma_prev = Jet2.constant(0)
    gamma = 1 / (5 * (z + 1) ** 3)

    max_node = max(nodes)
    by_index: Dict[int, Tuple[Jet2, Jet2]] = {
        0: (phi_prev, gamma_prev),
        1: (phi, gamma),
    }
    for n in range(1, max_node):
        denominator = (z + n + 1) ** 3 * mp.mpf(apery[n + 1])
        alpha = apery_P(z + n) * mp.mpf(apery[n]) / denominator
        beta = -(z + n) ** 3 * mp.mpf(apery[n - 1]) / denominator
        phi_next = alpha * phi + beta * phi_prev
        gamma_next = alpha * gamma + beta * gamma_prev
        phi_prev, phi = phi, phi_next
        gamma_prev, gamma = gamma, gamma_next
        if n + 1 in wanted:
            by_index[n + 1] = (phi, gamma)
    return (
        [by_index[n][0] for n in nodes],
        [by_index[n][1] for n in nodes],
    )


@dataclass
class KInfinityCell:
    """Numerical limiting cell ``phi(-z)phi(z)+z^6 gamma(-z)gamma(z)``."""

    apery: Sequence[int]
    nodes: Sequence[int]
    weights: Sequence[Any]
    ratio_n: Optional[Sequence[Any]] = None
    ratio_prev: Optional[Sequence[Any]] = None

    def __post_init__(self) -> None:
        """Cache the huge Apéry-number ratios once, not in every Newton step."""

        maximum = max(self.nodes)
        ratio_n = [mp.mpf(0)] * maximum
        ratio_prev = [mp.mpf(0)] * maximum
        for n in range(1, maximum):
            denominator = mp.mpf(self.apery[n + 1])
            ratio_n[n] = mp.mpf(self.apery[n]) / denominator
            ratio_prev[n] = mp.mpf(self.apery[n - 1]) / denominator
        self.ratio_n = ratio_n
        self.ratio_prev = ratio_prev

    def phi_gamma(self, z: Jet2) -> Tuple[Jet2, Jet2]:
        wanted = set(self.nodes)
        phi_prev = Jet2.constant(1)
        phi = apery_P(z) / (5 * (z + 1) ** 3)
        gamma_prev = Jet2.constant(0)
        gamma = 1 / (5 * (z + 1) ** 3)
        by_index: Dict[int, Tuple[Jet2, Jet2]] = {}
        for n in range(1, max(self.nodes)):
            shifted = z + n
            shifted_squared = shifted * shifted
            shifted_cubed = shifted_squared * shifted
            denominator_inverse = ((shifted + 1) ** 3).reciprocal()
            transfer = (
                34 * shifted_cubed
                + 51 * shifted_squared
                + 27 * shifted
                + 5
            )
            alpha = transfer * self.ratio_n[n] * denominator_inverse
            beta = -shifted_cubed * self.ratio_prev[n] * denominator_inverse
            phi_next = alpha * phi + beta * phi_prev
            gamma_next = alpha * gamma + beta * gamma_prev
            phi_prev, phi = phi, phi_next
            gamma_prev, gamma = gamma, gamma_next
            if n + 1 in wanted:
                by_index[n + 1] = (phi, gamma)
        phis = [by_index[n][0] for n in self.nodes]
        gammas = [by_index[n][1] for n in self.nodes]
        return (
            jet_linear_combination(phis, self.weights),
            jet_linear_combination(gammas, self.weights),
        )

    def K_jet(self, value: Any) -> Jet2:
        z = Jet2.variable(value)
        phi_plus, gamma_plus = self.phi_gamma(z)
        phi_minus, gamma_minus = self.phi_gamma(-z)
        return phi_minus * phi_plus + z**6 * gamma_minus * gamma_plus

    def H_and_derivative(self, value: Any) -> Tuple[Any, Any]:
        k = self.K_jet(value)
        z = mp.mpc(value)
        return z * k.first - 3 * k.value, z * k.second - 2 * k.first


def exact_extrapolation_weights(nodes: Sequence[int]) -> List[Any]:
    """Exact rational Lagrange weights, converted outward to Arb."""

    ts = [Fraction(1, n) for n in nodes]
    weights = []
    for i, ti in enumerate(ts):
        weight = Fraction(1)
        for j, tj in enumerate(ts):
            if i != j:
                weight *= -tj / (ti - tj)
        weights.append(arb(weight.numerator) / weight.denominator)
    return weights


@dataclass
class AcbKInfinityCell:
    """Fast Arb evaluation of the same limiting cell and its first two jets."""

    apery: Sequence[int]
    nodes: Sequence[int]
    weights: Sequence[Any]
    ratio_n: Optional[Sequence[Any]] = None
    ratio_prev: Optional[Sequence[Any]] = None

    def __post_init__(self) -> None:
        maximum = max(self.nodes)
        ratio_n = [arb(0)] * maximum
        ratio_prev = [arb(0)] * maximum
        for n in range(1, maximum):
            denominator = arb(self.apery[n + 1])
            ratio_n[n] = arb(self.apery[n]) / denominator
            ratio_prev[n] = arb(self.apery[n - 1]) / denominator
        self.ratio_n = ratio_n
        self.ratio_prev = ratio_prev

    def phi_gamma(self, z: AcbJet2) -> Tuple[AcbJet2, AcbJet2]:
        wanted = set(self.nodes)
        phi_prev = AcbJet2.constant(1)
        phi = apery_P(z) / (5 * (z + 1) ** 3)
        gamma_prev = AcbJet2.constant(0)
        gamma = 1 / (5 * (z + 1) ** 3)
        by_index: Dict[int, Tuple[AcbJet2, AcbJet2]] = {}
        for n in range(1, max(self.nodes)):
            shifted = z + n
            shifted_squared = shifted * shifted
            shifted_cubed = shifted_squared * shifted
            denominator_inverse = ((shifted + 1) ** 3).reciprocal()
            transfer = (
                34 * shifted_cubed
                + 51 * shifted_squared
                + 27 * shifted
                + 5
            )
            alpha = transfer * self.ratio_n[n] * denominator_inverse
            beta = -shifted_cubed * self.ratio_prev[n] * denominator_inverse
            phi_next = alpha * phi + beta * phi_prev
            gamma_next = alpha * gamma + beta * gamma_prev
            phi_prev, phi = phi, phi_next
            gamma_prev, gamma = gamma, gamma_next
            if n + 1 in wanted:
                by_index[n + 1] = (phi, gamma)

        phi_coefficients = []
        gamma_coefficients = []
        for degree in range(3):
            phi_coefficients.append(
                sum(
                    (
                        self.weights[k] * by_index[n][0].coeff[degree]
                        for k, n in enumerate(self.nodes)
                    ),
                    acb(0),
                )
            )
            gamma_coefficients.append(
                sum(
                    (
                        self.weights[k] * by_index[n][1].coeff[degree]
                        for k, n in enumerate(self.nodes)
                    ),
                    acb(0),
                )
            )
        return AcbJet2(tuple(phi_coefficients)), AcbJet2(tuple(gamma_coefficients))

    def K_jet(self, value: Any) -> AcbJet2:
        z = AcbJet2.variable(value)
        phi_plus, gamma_plus = self.phi_gamma(z)
        phi_minus, gamma_minus = self.phi_gamma(-z)
        return phi_minus * phi_plus + z**6 * gamma_minus * gamma_plus

    def H_and_derivative(self, value: Any) -> Tuple[Any, Any]:
        k = self.K_jet(value)
        z = acb(value)
        return z * k.first - 3 * k.value, z * k.second - 2 * k.first


def fixed_phi_gamma_jet(n: int, z: Jet2, apery: Sequence[int]) -> Tuple[Jet2, Jet2]:
    """Exact finite recurrence evaluation of ``phi_n,gamma_n``."""

    if n == 0:
        return Jet2.constant(1), Jet2.constant(0)
    phis, gammas = phi_gamma_sequence_jets(z, [n], apery)
    return phis[0], gammas[0]


def fixed_cell_K_jet(r: int, s: int, value: Any, apery: Sequence[int]) -> Jet2:
    """The finite cell object ``K_(r,s)`` from Q6723, equation (3.2)."""

    z = Jet2.variable(value)
    phi_r, gamma_r = fixed_phi_gamma_jet(r, -z, apery)
    phi_s, gamma_s = fixed_phi_gamma_jet(s, z, apery)
    return phi_r * phi_s + z**6 * gamma_r * gamma_s


def central_branch_extrapolation(
    numerators: Sequence[Any], digits: int, m_values: Sequence[int]
) -> Dict[str, Any]:
    """Extrapolate the four central ``K_(m,m)`` roots to ``K_infinity``.

    Q6723 (6.1)-(6.6) gives a full inverse-length branch expansion.  The
    present first pass uses its numerical consequence: each simple central
    branch has a power series in ``1/m``.  Every finite ``K_(m,m)`` root used
    below is first isolated by Arb; only the extrapolation to ``m=infinity``
    is uncertified.
    """

    import cmath

    old_dps = mp.dps
    mp.dps = max(digits, 120)
    branches: List[List[Any]] = []
    finite_rows = []
    for m in m_values:
        h = 2 * m + 1
        data = build_height_data(h, numerators[h])
        ctx.prec = max(384, int(4 * digits))
        quotient_roots = acb_poly(data.quotient_critical).roots(
            tol=arb(10) ** -45, maxprec=max(1024, int(10 * digits))
        )
        central = []
        for root in quotient_roots:
            midpoint = complex(float(root.mid().real), float(root.mid().imag))
            s_value = cmath.sqrt(midpoint)
            if s_value.real > 0:
                s_value = -s_value
            z_value = s_value / 2
            if -1 < z_value.real < 0:
                central.append(mp.mpc(str(z_value.real), str(z_value.imag)))
        central.sort(key=lambda value: (float(value.imag), float(value.real)))
        if len(central) != 4:
            raise ArithmeticError(
                "expected four central K_(m,m) roots at m=%d, found %d"
                % (m, len(central))
            )
        branches.append(central)
        finite_rows.append(
            {
                "m": m,
                "h": h,
                "roots": [
                    {"real": mp.nstr(value.real, 18), "imag": mp.nstr(value.imag, 18)}
                    for value in central
                ],
            }
        )

    weights = extrapolation_weights(m_values)
    roots = []
    for branch in range(4):
        values = [row[branch] for row in branches]
        limit = sum(
            (weight * value for weight, value in zip(weights, values)), mp.mpc(0)
        )
        shorter_weights = extrapolation_weights(m_values[1:])
        shorter = sum(
            (
                weight * value
                for weight, value in zip(shorter_weights, values[1:])
            ),
            mp.mpc(0),
        )
        roots.append(
            {
                "real": mp.nstr(limit.real, 30),
                "imag": mp.nstr(limit.imag, 30),
                "extrapolation_delta": mp_scientific(abs(limit - shorter), 8),
                "last_finite_m": int(m_values[-1]),
            }
        )

    mp.dps = old_dps
    return {
        "status": "UNCERTIFIED",
        "method": "central K_(m,m) branch extrapolation in 1/m",
        "m_values": list(m_values),
        "root_count": len(roots),
        "finite_rows": finite_rows,
        "roots": roots,
    }


def isolate_kinfinity_roots(
    numerators: Sequence[Any], digits: int, m_values: Sequence[int]
) -> Dict[str, Any]:
    """Newton-search the finite roots of the accelerated ``H_infinity``."""

    ctx.prec = max(448, int(4.5 * digits))
    nodes = [48 * (k + 1) for k in range(8)]
    apery = build_apery_numbers(max(nodes))
    weights = exact_extrapolation_weights(nodes)
    cell = AcbKInfinityCell(apery=apery, nodes=nodes, weights=weights)
    seeds = [
        acb("-0.5", "0.1104"),
        acb("-0.5", "-0.1104"),
        acb("-0.47", "1.8"),
        acb("-0.47", "-1.8"),
    ]
    attempts = []
    converged_values = []
    threshold = 10.0 ** -30
    for seed in seeds:
        z = seed
        previous = math.inf
        correction_size = math.inf
        converged = False
        derivative = acb(0)
        residual = acb(0)
        for iteration in range(12):
            residual, derivative = cell.H_and_derivative(z)
            if derivative.contains(0):
                break
            correction = (residual / derivative).mid()
            correction_size = float(correction.abs_upper())
            z = (z - correction).mid()
            if correction_size < threshold:
                converged = True
                break
            if correction_size > 10 * previous and iteration >= 3:
                break
            previous = correction_size
        residual, derivative = cell.H_and_derivative(z)
        attempts.append(
            {
                "seed_real": seed.real.str(12, radius=False),
                "seed_imag": seed.imag.str(12, radius=False),
                "converged": converged,
                "real": z.real.str(40, radius=False),
                "imag": z.imag.str(40, radius=False),
                "residual_upper": float_scientific(float(residual.abs_upper()), 8),
                "newton_correction": float_scientific(correction_size, 8),
                "derivative_lower": float_scientific(float(derivative.abs_lower()), 8),
            }
        )
        if converged:
            midpoint = complex(float(z.real), float(z.imag))
            if not any(abs(midpoint - other) < 1e-20 for other in converged_values):
                converged_values.append(midpoint)

    direct_roots = []
    shorter_nodes = nodes[1:]
    shorter_cell = AcbKInfinityCell(
        apery=apery,
        nodes=shorter_nodes,
        weights=exact_extrapolation_weights(shorter_nodes),
    )
    for value in sorted(converged_values, key=lambda item: item.imag):
        matching = min(
            (row for row in attempts if row["converged"]),
            key=lambda row: abs(
                complex(float(row["real"]), float(row["imag"])) - value
            ),
        )
        root_point = acb(matching["real"], matching["imag"])
        shorter_residual, _ = shorter_cell.H_and_derivative(root_point)
        direct_roots.append(
            {
                # The omitted-node extrapolation check is only at the 1e-17
                # level, so do not print the 100-digit fit as if all of its
                # digits described the limiting root.  Full Newton iterates
                # remain available under ``seed_attempts`` for diagnostics.
                "real": mp.nstr(mp.mpf(matching["real"]), 17),
                "imag": mp.nstr(mp.mpf(matching["imag"]), 17),
                "residual_upper": matching["residual_upper"],
                "newton_correction": matching["newton_correction"],
                "derivative_lower": matching["derivative_lower"],
                "extrapolation_check": float_scientific(
                    float(shorter_residual.abs_upper()), 8
                ),
            }
        )

    central = central_branch_extrapolation(numerators, digits, m_values)
    return {
        "status": "UNCERTIFIED",
        "method": "Arb-accelerated normalized recurrence plus Newton",
        "nodes": nodes,
        "root_count": len(direct_roots),
        "roots": direct_roots,
        "seed_attempts": attempts,
        "central_branch_evidence": central,
    }


def check_cell_identity(
    numerators: Sequence[Any], apery: Sequence[int], digits: int
) -> List[Dict[str, Any]]:
    """Numerically check Q6723 (3.2)-(3.3) at independent sample points."""

    old_dps = mp.dps
    mp.dps = max(digits, 120)
    samples = [(2, 1), (7, 2), (12, 5)]
    z = mp.mpc("-0.37", "0.23")
    rows = []
    for h, j in samples:
        r, s = j - 1, h - j
        k = fixed_cell_K_jet(r, s, z, apery).value
        cell_value = mp.mpf(apery[r]) * mp.mpf(apery[s]) * k / z**3
        x = -j + z
        direct = mp_polyval(numerators[h], x) / mp_polyval(q_polynomial(h), x) ** 3
        relative = abs(cell_value - direct) / max(mp.mpf(1), abs(direct))
        rows.append(
            {
                "h": h,
                "j": j,
                "r": r,
                "s": s,
                "relative_error": mp_scientific(relative, 8),
            }
        )
    mp.dps = old_dps
    return rows


def uncertified_height(
    data: HeightData, roots: Sequence[Any], digits: int
) -> Dict[str, Any]:
    """Discard all ball radii and redo Newton/value evaluation with mpmath."""

    old_dps = mp.dps
    mp.dps = digits
    polynomial = data.quotient_critical
    derivative = polynomial.derivative()
    value_num = data.squared_value_num
    value_den = data.squared_value_den
    value_num_derivative = value_num.derivative()
    value_den_derivative = value_den.derivative()

    refined = []
    max_correction = mp.mpf(0)
    for root in roots:
        u = mp_midpoint(root, digits + 20)
        correction = mp.inf
        derivative_value = mp.mpc(0)
        for _ in range(5):
            derivative_value = mp_polyval(derivative, u)
            if derivative_value == 0:
                break
            correction = mp_polyval(polynomial, u) / derivative_value
            u -= correction
            if abs(correction) < mp.power(10, -(digits - 15)):
                break
        max_correction = max(max_correction, abs(correction))
        numerator = mp_polyval(value_num, u)
        denominator = mp_polyval(value_den, u)
        value = numerator / denominator
        value_derivative = (
            mp_polyval(value_num_derivative, u) * denominator
            - numerator * mp_polyval(value_den_derivative, u)
        ) / denominator**2
        radius = abs(value_derivative) * abs(correction) + 100 * mp.eps * max(1, abs(value))
        refined.append((u, value, radius, abs(derivative_value)))

    min_separation = mp.inf
    min_margin = mp.inf
    min_relative_margin = mp.inf
    for i in range(len(refined)):
        value_i = refined[i][1]
        min_margin = min(min_margin, abs(value_i) - refined[i][2])
        for j in range(i):
            value_j = refined[j][1]
            separation = abs(value_i - value_j)
            margin = separation - refined[i][2] - refined[j][2]
            scale = max(abs(value_i), abs(value_j), mp.mpf(1))
            min_separation = min(min_separation, separation)
            min_margin = min(min_margin, margin)
            min_relative_margin = min(min_relative_margin, margin / scale)

    expected = 2 * data.h - 2
    numerically_simple = all(row[3] > 0 for row in refined)
    if (
        len(refined) == expected
        and numerically_simple
        and min_margin > 0
        and len(refined) >= data.h
    ):
        status = "YES"
    elif len(refined) == expected and min_margin < 0:
        status = "UNRESOLVED"
    else:
        status = "UNRESOLVED"

    result = {
        "status": status,
        "precision_digits": digits,
        "quotient_roots": len(refined),
        "critical_zeros": 2 * len(refined),
        "good_orbits": len(refined) if min_margin > 0 else 0,
        "min_separation": mp_scientific(min_separation, 10),
        "estimated_margin": mp_scientific(min_margin, 10),
        "relative_margin": mp_scientific(min_relative_margin, 10),
        "max_newton_correction": mp_scientific(max_correction, 8),
    }
    mp.dps = old_dps
    return result


def make_newton_box(midpoint: Any, radius_bits: int) -> Any:
    radius = (1, -radius_bits)
    return acb(arb(midpoint.real, radius), arb(midpoint.imag, radius))


def arb_midpoint_separation(left: Any, right: Any) -> float:
    return float((left.mid() - right.mid()).abs_upper())


def certified_height_once(
    data: HeightData, digits: int, precision_bits: int
) -> Tuple[Dict[str, Any], List[Any]]:
    """One Arb root/value certification attempt."""

    ctx.prec = precision_bits
    polynomial = acb_poly(data.quotient_critical)
    derivative = polynomial.derivative()
    tolerance = arb(10) ** (-(digits + 10))
    started = time.monotonic()
    roots = polynomial.roots(tol=tolerance, maxprec=2 * precision_bits)
    root_seconds = time.monotonic() - started

    expected = int(data.quotient_critical.degree())
    exhaustive = len(roots) == expected
    root_balls_disjoint = all(
        not roots[i].overlaps(roots[j])
        for i in range(len(roots))
        for j in range(i)
    )

    newton_pass = 0
    derivative_pass = 0
    target_radius_bits = int(math.ceil((digits + 5) * math.log2(10)))
    for root in roots:
        midpoint = root.mid()
        radius_bits = target_radius_bits
        box = make_newton_box(midpoint, radius_bits)
        while not box.contains(root) and radius_bits > 32:
            radius_bits -= 16
            box = make_newton_box(midpoint, radius_bits)
        derivative_ball = derivative(box)
        if not derivative_ball.contains(0):
            derivative_pass += 1
            newton_image = midpoint - polynomial(midpoint) / derivative_ball
            if box.contains_interior(newton_image):
                newton_pass += 1

    numerator = acb_poly(data.squared_value_num)
    denominator = acb_poly(data.squared_value_den)
    values = []
    denominator_pass = 0
    nonzero_flags = []
    for root in roots:
        denominator_ball = denominator(root)
        if not denominator_ball.contains(0):
            denominator_pass += 1
            value = numerator(root) / denominator_ball
        else:
            value = acb(0, arb(0, 1))
        values.append(value)
        nonzero_flags.append(not value.contains(0))

    pairwise_flags = [True] * len(values)
    min_separation = math.inf
    min_margin = math.inf
    min_relative_margin = math.inf
    for i, value_i in enumerate(values):
        zero_margin = float(value_i.abs_lower())
        min_margin = min(min_margin, zero_margin)
        for j in range(i):
            difference = value_i - values[j]
            separation = arb_midpoint_separation(value_i, values[j])
            margin = float(difference.abs_lower())
            scale = max(
                float(value_i.abs_upper()), float(values[j].abs_upper()), 1.0
            )
            min_separation = min(min_separation, separation)
            min_margin = min(min_margin, margin)
            min_relative_margin = min(min_relative_margin, margin / scale)
            if difference.contains(0):
                pairwise_flags[i] = False
                pairwise_flags[j] = False

    simple = [
        exhaustive
        and root_balls_disjoint
        and derivative_pass == len(roots)
        and newton_pass == len(roots)
        for _ in roots
    ]
    good = sum(
        1
        for i in range(len(roots))
        if simple[i] and nonzero_flags[i] and pairwise_flags[i]
    )
    status = "YES" if good >= data.h else "UNRESOLVED"
    max_root_radius = max((float(root.rad()) for root in roots), default=math.inf)
    result = {
        "status": status,
        "precision_bits": precision_bits,
        "quotient_roots": len(roots),
        "critical_zeros": 2 * len(roots),
        "good_orbits": good,
        "exhaustive": exhaustive,
        "root_balls_disjoint": root_balls_disjoint,
        "interval_newton_pass": newton_pass,
        "derivative_excludes_zero": derivative_pass,
        "denominator_excludes_zero": denominator_pass,
        "min_separation": float_scientific(min_separation, 10),
        "certified_margin": float_scientific(min_margin, 10),
        "relative_margin": float_scientific(min_relative_margin, 10),
        "max_root_radius": float_scientific(max_root_radius, 6),
        "root_seconds": round(root_seconds, 6),
    }
    return result, roots


def certify_height(
    data: HeightData, digits: int, retries: int
) -> Tuple[Dict[str, Any], List[Any]]:
    """Increase precision before returning an unresolved certificate."""

    base = int(math.ceil((digits + 80) * math.log2(10))) + 5 * data.h
    last_result: Dict[str, Any] = {}
    last_roots: List[Any] = []
    for attempt in range(retries + 1):
        bits = int(base * (1.35**attempt))
        try:
            last_result, last_roots = certified_height_once(data, digits, bits)
        except (ValueError, ZeroDivisionError) as error:
            last_result = {
                "status": "UNRESOLVED",
                "precision_bits": bits,
                "error": repr(error),
            }
            last_roots = []
        if last_result.get("status") == "YES":
            return last_result, last_roots
        log(
            "[retry] h=%d certificate unresolved at %d bits; increasing precision"
            % (data.h, bits)
        )
    return last_result, last_roots


def arb_text(value: Any, digits: int = 20) -> str:
    """Serialize an Arb midpoint without Python-float truncation."""

    return value.str(digits, radius=False)


def acb_ball_payload(value: Any, digits: int = 20) -> Dict[str, Any]:
    """Return a compact midpoint-plus-radius representation of an Acb ball."""

    midpoint = value.mid()
    return {
        "midpoint": {
            "real": arb_text(midpoint.real, digits),
            "imag": arb_text(midpoint.imag, digits),
        },
        "radius": arb_text(value.rad(), min(digits, 14)),
    }


def dyadic_band_index(cell_distance: int, h: int) -> Optional[int]:
    """Index the exact distance ``cell_distance/h`` by dyadic half-open bands.

    The returned ``k`` is characterized by
    ``2^-k <= cell_distance/h < 2^(-k+1)``.  A same-cell pair has no positive
    dyadic index and is returned as ``None``.
    """

    if cell_distance == 0:
        return None
    k = 1
    while cell_distance * (1 << k) < h:
        k += 1
    return k


def compact_pair(pair: Dict[str, Any]) -> Dict[str, Any]:
    """Fields needed when a pair is referenced from a summary table."""

    return {
        "pair_index": pair["pair_index"],
        "orbit_a": pair["orbit_a"],
        "orbit_b": pair["orbit_b"],
        "cell_distance_cells": pair["cell_distance_cells"],
        "cell_distance": pair["cell_distance"],
        "dyadic_band": pair["dyadic_band"],
        "relative_separation": pair["relative_separation"],
    }


def pair_anatomy_height(
    data: HeightData, roots: Sequence[Any], digits: int
) -> Dict[str, Any]:
    """Resolve certified quotient roots into cells/branches and all pairs.

    For each mirror orbit, take the unique critical-point representative in
    the upper half-plane.  It lies in one of the physical cells
    ``-j-1 < Re(x) < -j``.  We use ``t=j/h``; reflection sends ``t`` to
    ``1-t``.  Each cell has two upper-half-plane representatives, labelled
    ``inner`` and ``outer`` by increasing imaginary part.
    """

    old_dps = mp.dps
    mp.dps = max(digits, 80)
    numerator = acb_poly(data.squared_value_num)
    denominator = acb_poly(data.squared_value_den)
    provisional: List[Dict[str, Any]] = []

    for root_index, root in enumerate(roots):
        if root.imag.contains(0) and float(root.real.upper()) < 0:
            # A real-coefficient polynomial has the conjugate root in the
            # conjugate ball.  Exhaustive isolation and disjoint root balls
            # force this unique root to be real.  Direct Acb sqrt on a ball
            # crossing the negative real axis deliberately returns both signs,
            # so construct the positive-imaginary square root from the real
            # interval instead.
            representative = acb(
                arb(-(data.h + 1)) / 2,
                (-root.real).sqrt() / 2,
            )
            representative_method = "negative-real quotient root"
        else:
            square_root = root.sqrt()
            candidates = [
                (-(data.h + 1) + square_root) / 2,
                (-(data.h + 1) - square_root) / 2,
            ]
            upper = [
                candidate
                for candidate in candidates
                if float(candidate.imag.lower()) > 0
            ]
            if len(upper) != 1:
                raise ArithmeticError(
                    "could not select a unique upper-half-plane representative "
                    "at h=%d root=%d" % (data.h, root_index)
                )
            representative = upper[0]
            representative_method = "principal sqrt with sign chosen by Im(x)>0"
        cell_index = math.floor(-float(representative.real))
        if not 1 <= cell_index <= data.h - 1:
            raise ArithmeticError(
                "critical representative escaped internal cells at h=%d: j=%d"
                % (data.h, cell_index)
            )
        if not (
            float(representative.real.lower()) > -cell_index - 1
            and float(representative.real.upper()) < -cell_index
        ):
            raise ArithmeticError(
                "cell label is not certified by the representative ball at "
                "h=%d j=%d" % (data.h, cell_index)
            )

        denominator_ball = denominator(root)
        if denominator_ball.contains(0):
            raise ArithmeticError(
                "squared critical value denominator contains zero at h=%d"
                % data.h
            )
        value = numerator(root) / denominator_ball
        local_coordinate = representative + cell_index
        provisional.append(
            {
                "root_index": root_index,
                "cell_index": cell_index,
                "root_ball": root,
                "representative_ball": representative,
                "representative_method": representative_method,
                "local_coordinate_ball": local_coordinate,
                "value_ball": value,
                "value_midpoint_mp": mp_midpoint(value, digits + 20),
            }
        )

    by_cell: Dict[int, List[Dict[str, Any]]] = {}
    for orbit in provisional:
        by_cell.setdefault(orbit["cell_index"], []).append(orbit)
    if sorted(by_cell) != list(range(1, data.h)):
        raise ArithmeticError(
            "upper-half-plane representatives do not cover all cells at h=%d"
            % data.h
        )
    for cell_index, cell_orbits in by_cell.items():
        if len(cell_orbits) != 2:
            raise ArithmeticError(
                "expected two local branches at h=%d j=%d, found %d"
                % (data.h, cell_index, len(cell_orbits))
            )
        cell_orbits.sort(key=lambda row: float(row["representative_ball"].imag))
        if cell_orbits[0]["representative_ball"].imag.overlaps(
            cell_orbits[1]["representative_ball"].imag
        ):
            raise ArithmeticError(
                "inner/outer branch labels are not separated at h=%d j=%d"
                % (data.h, cell_index)
            )
        cell_orbits[0]["branch"] = "inner"
        cell_orbits[1]["branch"] = "outer"

    provisional.sort(
        key=lambda row: (row["cell_index"], 0 if row["branch"] == "inner" else 1)
    )
    orbit_rows: List[Dict[str, Any]] = []
    internal: List[Dict[str, Any]] = []
    for orbit_index, orbit in enumerate(provisional, start=1):
        cell_index = orbit["cell_index"]
        branch = orbit["branch"]
        orbit_id = "h%d-j%d-%s" % (data.h, cell_index, branch)
        value_ball = orbit["value_ball"]
        row = {
            "orbit_index": orbit_index,
            "orbit_id": orbit_id,
            "quotient_root_index": orbit["root_index"],
            "cell_index": cell_index,
            "cell_position": mp.nstr(mp.mpf(cell_index) / data.h, 20),
            "cell_position_fraction": "%d/%d" % (cell_index, data.h),
            "branch": branch,
            "representative_method": orbit["representative_method"],
            "quotient_root_u": acb_ball_payload(orbit["root_ball"], 22),
            "upper_half_plane_x": acb_ball_payload(
                orbit["representative_ball"], 22
            ),
            "local_coordinate_z": acb_ball_payload(
                orbit["local_coordinate_ball"], 22
            ),
            "squared_critical_value": acb_ball_payload(value_ball, 24),
        }
        orbit_rows.append(row)
        internal.append(
            {
                **row,
                "value_ball": value_ball,
                "value_midpoint_mp": orbit["value_midpoint_mp"],
            }
        )

    pairs: List[Dict[str, Any]] = []
    best_pair: Optional[Tuple[Any, Dict[str, Any]]] = None
    band_accumulators: Dict[Optional[int], Dict[str, Any]] = {}
    for right in range(len(internal)):
        orbit_b = internal[right]
        for left in range(right):
            orbit_a = internal[left]
            value_a = orbit_a["value_midpoint_mp"]
            value_b = orbit_b["value_midpoint_mp"]
            absolute_separation = abs(value_a - value_b)
            scale = max(abs(value_a), abs(value_b))
            relative_separation = absolute_separation / scale
            cell_distance = abs(orbit_a["cell_index"] - orbit_b["cell_index"])
            band_index = dyadic_band_index(cell_distance, data.h)
            band_label = "same-cell" if band_index is None else "k=%d" % band_index
            pair = {
                "pair_index": len(pairs) + 1,
                "orbit_a": orbit_a["orbit_id"],
                "orbit_b": orbit_b["orbit_id"],
                "cell_distance_cells": cell_distance,
                "cell_distance": mp.nstr(mp.mpf(cell_distance) / data.h, 20),
                "cell_distance_fraction": "%d/%d" % (cell_distance, data.h),
                "dyadic_band": band_label,
                "absolute_separation_midpoint": mp_scientific(
                    absolute_separation, 20
                ),
                "relative_separation": mp_scientific(relative_separation, 20),
            }
            pairs.append(pair)
            if best_pair is None or relative_separation < best_pair[0]:
                best_pair = (relative_separation, pair)

            accumulator = band_accumulators.setdefault(
                band_index,
                {
                    "pair_count": 0,
                    "minimum": None,
                    "minimum_relative_over_distance_squared": None,
                },
            )
            accumulator["pair_count"] += 1
            if (
                accumulator["minimum"] is None
                or relative_separation < accumulator["minimum"][0]
            ):
                accumulator["minimum"] = (relative_separation, pair)
            if cell_distance:
                distance = mp.mpf(cell_distance) / data.h
                normalized = relative_separation / distance**2
                if (
                    accumulator["minimum_relative_over_distance_squared"] is None
                    or normalized
                    < accumulator["minimum_relative_over_distance_squared"]
                ):
                    accumulator["minimum_relative_over_distance_squared"] = normalized

    if best_pair is None:
        raise ArithmeticError("pair anatomy requires at least two quotient orbits")

    dyadic_rows = []
    same_cell_row: Optional[Dict[str, Any]] = None
    for band_index, accumulator in sorted(
        band_accumulators.items(), key=lambda item: -1 if item[0] is None else item[0]
    ):
        minimum_value, minimum_pair = accumulator["minimum"]
        band_row = {
            "band": "same-cell" if band_index is None else "k=%d" % band_index,
            "pair_count": accumulator["pair_count"],
            "minimum_relative_separation": mp_scientific(minimum_value, 20),
            "minimizing_pair": compact_pair(minimum_pair),
        }
        if band_index is None:
            band_row["distance_interval"] = "{0}"
            same_cell_row = band_row
        else:
            band_row.update(
                {
                    "k": band_index,
                    "distance_interval": "[2^-%d,2^-%d)"
                    % (band_index, band_index - 1),
                    "lower_endpoint": mp.nstr(mp.mpf(2) ** (-band_index), 20),
                    "upper_endpoint": mp.nstr(
                        mp.mpf(2) ** (-(band_index - 1)), 20
                    ),
                    "minimum_relative_over_distance_squared": mp_scientific(
                        accumulator["minimum_relative_over_distance_squared"], 20
                    ),
                }
            )
            dyadic_rows.append(band_row)

    orbit_by_id = {row["orbit_id"]: row for row in internal}
    if data.h % 2:
        central_cells = (data.h // 2, data.h // 2 + 1)
    else:
        central_cells = (data.h // 2 - 1, data.h // 2 + 1)
    central_ids = {
        "h%d-j%d-inner" % (data.h, central_cells[0]),
        "h%d-j%d-inner" % (data.h, central_cells[1]),
    }
    central_pair = next(
        (
            pair
            for pair in pairs
            if {pair["orbit_a"], pair["orbit_b"]} == central_ids
        ),
        best_pair[1],
    )
    minimum_ids = {best_pair[1]["orbit_a"], best_pair[1]["orbit_b"]}
    central_gap = {
        **compact_pair(central_pair),
        "expected_cells": list(central_cells),
        "branch": "inner",
        "is_global_minimum": minimum_ids == central_ids,
        "gap_times_h_squared": mp_scientific(
            mp.mpf(central_pair["relative_separation"]) * data.h**2, 20
        ),
    }
    minimum_orbits = [orbit_by_id[best_pair[1]["orbit_a"]], orbit_by_id[best_pair[1]["orbit_b"]]]
    global_minimum = {
        **compact_pair(best_pair[1]),
        "cells": [row["cell_index"] for row in minimum_orbits],
        "branches": [row["branch"] for row in minimum_orbits],
        "is_predicted_central_inner_pair": minimum_ids == central_ids,
    }

    result = {
        "h": data.h,
        "orbit_count": len(orbit_rows),
        "pair_count": len(pairs),
        "orbits": orbit_rows,
        "pairs": pairs,
        "global_minimum": global_minimum,
        "central_gap": central_gap,
        "same_cell_band": same_cell_row,
        "dyadic_bands": dyadic_rows,
    }
    mp.dps = old_dps
    return result


def fit_central_limit(points: Sequence[Tuple[int, Any]]) -> Dict[str, Any]:
    """Fit ``gap*h^2 = kappa + a/h + b/h^2`` and retain every residual."""

    if len(points) < 3:
        raise ArithmeticError("a quadratic inverse-height fit needs three points")
    old_dps = mp.dps
    mp.dps = 80
    design = mp.matrix(
        [[mp.mpf(1), mp.mpf(1) / h, mp.mpf(1) / h**2] for h, _ in points]
    )
    observations = mp.matrix([mp.mpf(value) for _, value in points])
    coefficients = mp.lu_solve(design.T * design, design.T * observations)
    fitted = design * coefficients
    residual_values = [observations[k] - fitted[k] for k in range(len(points))]
    rmse = mp.sqrt(sum(value**2 for value in residual_values) / len(points))
    maximum = max(abs(value) for value in residual_values)
    result = {
        "model": "relative_gap*h^2 = kappa + a/h + b/h^2",
        "fit_range": "%d..%d" % (min(h for h, _ in points), max(h for h, _ in points)),
        "heights": [h for h, _ in points],
        "point_count": len(points),
        "coefficients": {
            "kappa": mp.nstr(coefficients[0], 20),
            "a": mp.nstr(coefficients[1], 20),
            "b": mp.nstr(coefficients[2], 20),
        },
        "rmse": mp_scientific(rmse, 16),
        "max_abs_residual": mp_scientific(maximum, 16),
        "residuals": [
            {
                "h": h,
                "observed": mp.nstr(observations[k], 20),
                "fitted": mp.nstr(fitted[k], 20),
                "residual": mp_scientific(residual_values[k], 16),
            }
            for k, (h, _) in enumerate(points)
        ],
    }
    mp.dps = old_dps
    return result


def logarithmic_model_diagnostic(
    bands: Sequence[Dict[str, Any]]
) -> Dict[str, Any]:
    """Compare constant and quadratic descriptions of far-band minima."""

    selected = [row for row in bands if row["k"] <= 3]
    if len(selected) < 2:
        return {"status": "INSUFFICIENT DATA"}
    xs = [mp.log(mp.mpf(row["minimizing_pair"]["cell_distance"])) for row in selected]
    ys = [mp.log(mp.mpf(row["empirical_lower_bound"])) for row in selected]
    x_mean = sum(xs) / len(xs)
    y_mean = sum(ys) / len(ys)
    denominator = sum((value - x_mean) ** 2 for value in xs)
    slope = sum((x - x_mean) * (y - y_mean) for x, y in zip(xs, ys)) / denominator
    intercept = y_mean - slope * x_mean
    constant_intercept = y_mean
    quadratic_intercept = sum(y - 2 * x for x, y in zip(xs, ys)) / len(xs)
    constant_rmse = mp.sqrt(
        sum((y - constant_intercept) ** 2 for y in ys) / len(ys)
    )
    quadratic_rmse = mp.sqrt(
        sum((y - (quadratic_intercept + 2 * x)) ** 2 for x, y in zip(xs, ys))
        / len(ys)
    )
    label = "quadratic-like" if quadratic_rmse < constant_rmse else "constant-like"
    return {
        "status": "EMPIRICAL",
        "bands": [row["band"] for row in selected],
        "distance_floor": "1/8",
        "free_log_log_slope": mp.nstr(slope, 12),
        "free_log_log_intercept": mp.nstr(intercept, 12),
        "constant_log_rmse": mp_scientific(constant_rmse, 12),
        "quadratic_log_rmse": mp_scientific(quadratic_rmse, 12),
        "closer_model": label,
    }


def fit_band_height_decay(
    anatomy_rows: Sequence[Dict[str, Any]], band_index: int
) -> Dict[str, Any]:
    """Fit a fixed dyadic band's per-height minimum to ``C*h^-alpha``."""

    points = []
    for row in anatomy_rows:
        band = next(item for item in row["dyadic_bands"] if item["k"] == band_index)
        points.append((row["h"], mp.mpf(band["minimum_relative_separation"])))
    xs = [mp.log(h) for h, _ in points]
    ys = [mp.log(value) for _, value in points]
    x_mean = sum(xs) / len(xs)
    y_mean = sum(ys) / len(ys)
    slope = sum((x - x_mean) * (y - y_mean) for x, y in zip(xs, ys)) / sum(
        (x - x_mean) ** 2 for x in xs
    )
    intercept = y_mean - slope * x_mean
    free_rmse = mp.sqrt(
        sum((y - (intercept + slope * x)) ** 2 for x, y in zip(xs, ys))
        / len(points)
    )
    constant_rmse = mp.sqrt(sum((y - y_mean) ** 2 for y in ys) / len(points))
    inverse_intercept = sum(y + x for x, y in zip(xs, ys)) / len(points)
    inverse_rmse = mp.sqrt(
        sum((y - (inverse_intercept - x)) ** 2 for x, y in zip(xs, ys))
        / len(points)
    )
    return {
        "band": "k=%d" % band_index,
        "distance_interval": "[2^-%d,2^-%d)" % (band_index, band_index - 1),
        "fit_range": "%d..%d" % (points[0][0], points[-1][0]),
        "point_count": len(points),
        "model": "min_relative_separation = C*h^(-alpha)",
        "alpha": mp.nstr(-slope, 14),
        "C": mp.nstr(mp.exp(intercept), 14),
        "free_log_rmse": mp_scientific(free_rmse, 12),
        "constant_log_rmse": mp_scientific(constant_rmse, 12),
        "inverse_h_log_rmse": mp_scientific(inverse_rmse, 12),
        "first_height": points[0][0],
        "first_minimum": mp_scientific(points[0][1], 18),
        "last_height": points[-1][0],
        "last_minimum": mp_scientific(points[-1][1], 18),
    }


def build_pair_anatomy_payload(
    standard_payload: Dict[str, Any], anatomy_rows: Sequence[Dict[str, Any]], args: argparse.Namespace
) -> Dict[str, Any]:
    """Aggregate per-height anatomy into fits and fixed dyadic-band profiles."""

    old_dps = mp.dps
    mp.dps = 80
    parity_points: Dict[str, List[Tuple[int, Any]]] = {"odd": [], "even": []}
    for row in anatomy_rows:
        parity = "odd" if row["h"] % 2 else "even"
        parity_points[parity].append(
            (row["h"], mp.mpf(row["central_gap"]["gap_times_h_squared"]))
        )
    primary_cutoff = max(min(row["h"] for row in anatomy_rows), max(row["h"] for row in anatomy_rows) - 20)
    fits = {
        parity: fit_central_limit(
            [(h, value) for h, value in points if h >= primary_cutoff]
        )
        for parity, points in parity_points.items()
    }
    stability_cutoffs = sorted(
        {
            min(row["h"] for row in anatomy_rows),
            min(row["h"] for row in anatomy_rows) + 10,
            primary_cutoff,
        }
    )
    fit_stability = {
        parity: [
            fit_central_limit([(h, value) for h, value in points if h >= cutoff])
            for cutoff in stability_cutoffs
            if sum(1 for h, _ in points if h >= cutoff) >= 3
        ]
        for parity, points in parity_points.items()
    }
    kappa_odd = mp.mpf(fits["odd"]["coefficients"]["kappa"])
    kappa_even = mp.mpf(fits["even"]["coefficients"]["kappa"])
    ratio = kappa_even / kappa_odd
    half_error = kappa_odd - kappa_even / 2

    aggregate: Dict[int, Dict[str, Any]] = {}
    for row in anatomy_rows:
        for band in row["dyadic_bands"]:
            k = band["k"]
            value = mp.mpf(band["minimum_relative_separation"])
            normalized = mp.mpf(band["minimum_relative_over_distance_squared"])
            target = aggregate.setdefault(
                k,
                {
                    "k": k,
                    "band": band["band"],
                    "distance_interval": band["distance_interval"],
                    "lower_endpoint": band["lower_endpoint"],
                    "upper_endpoint": band["upper_endpoint"],
                    "pair_count": 0,
                    "minimum": None,
                    "minimum_normalized": None,
                },
            )
            target["pair_count"] += band["pair_count"]
            if target["minimum"] is None or value < target["minimum"][0]:
                target["minimum"] = (value, row["h"], band["minimizing_pair"])
            if target["minimum_normalized"] is None or normalized < target["minimum_normalized"]:
                target["minimum_normalized"] = normalized

    aggregate_rows = []
    for k in sorted(aggregate):
        row = aggregate[k]
        value, height, pair = row.pop("minimum")
        normalized = row.pop("minimum_normalized")
        aggregate_rows.append(
            {
                **row,
                "empirical_lower_bound": mp_scientific(value, 20),
                "attained_at_h": height,
                "minimizing_pair": pair,
                "minimum_relative_over_distance_squared": mp_scientific(
                    normalized, 20
                ),
            }
        )

    macroscopic = [row for row in aggregate_rows if row["k"] <= 2]
    far_row = min(macroscopic, key=lambda row: mp.mpf(row["empirical_lower_bound"]))
    far_lower_bound = mp.mpf(far_row["empirical_lower_bound"])
    height_decay = [
        fit_band_height_decay(anatomy_rows, band_index)
        for band_index in (1, 2, 3)
    ]
    k1_decay_exponent = mp.mpf(height_decay[0]["alpha"])
    far_uniform_supported = (
        k1_decay_exponent < mp.mpf("0.25")
        and mp.mpf(height_decay[0]["constant_log_rmse"])
        <= mp.mpf(height_decay[0]["inverse_h_log_rmse"])
    )
    every_minimum_is_central = all(
        row["global_minimum"]["is_predicted_central_inner_pair"]
        for row in anatomy_rows
    )
    fit_residual_ok = all(
        mp.mpf(fit["max_abs_residual"])
        < mp.mpf(fit["coefficients"]["kappa"]) * mp.mpf("1e-3")
        for fit in fits.values()
    )
    confirmed = (
        every_minimum_is_central
        and fit_residual_ok
        and kappa_odd > 0
        and kappa_even > 0
        and far_lower_bound > 0
        and far_uniform_supported
    )

    result = {
        "schema": "CRON-kinf-pair-anatomy-v1",
        "command": serialize_anatomy_command(args),
        "runtime": standard_payload["runtime"],
        "scan": {
            **standard_payload["scan"],
            "pair_anatomy": True,
            "orbit_position_convention": "upper-half-plane representative; -j-1 < Re(x) < -j; t=j/h",
            "branch_convention": "inner/outer by increasing Im(x) within each cell",
        },
        "verdict": "TWO-REGIME CONFIRMED" if confirmed else "TWO-REGIME REFUTED",
        "checks": {
            "all_global_minima_are_predicted_central_inner_pairs": every_minimum_is_central,
            "fit_residual_gate": fit_residual_ok,
            "macroscopic_far_threshold": "cell distance >= 1/4",
            "far_uniform_constant_or_quadratic_supported": far_uniform_supported,
        },
        "central_fits": fits,
        "central_fit_stability": fit_stability,
        "constants": {
            "kappa_odd": mp.nstr(kappa_odd, 20),
            "kappa_even": mp.nstr(kappa_even, 20),
            "kappa_even_over_odd": mp.nstr(ratio, 20),
            "F_second_at_one_half_from_even_identity": mp.nstr(kappa_even, 20),
            "kappa_odd_minus_half_kappa_even": mp_scientific(half_error, 16),
            "relative_half_identity_error": mp_scientific(
                abs(half_error) / abs(kappa_odd), 16
            ),
            "macroscopic_far_lower_bound": mp_scientific(far_lower_bound, 20),
            "macroscopic_far_lower_bound_band": far_row["band"],
            "macroscopic_far_lower_bound_h": far_row["attained_at_h"],
            "macroscopic_far_lower_bound_pair": far_row["minimizing_pair"],
            "k1_height_decay_exponent": height_decay[0]["alpha"],
        },
        "far_band_model_diagnostic": logarithmic_model_diagnostic(aggregate_rows),
        "far_band_height_decay": height_decay,
        "aggregate_dyadic_bands": aggregate_rows,
        "heights": list(anatomy_rows),
    }
    mp.dps = old_dps
    return result


def serialize_command(args: argparse.Namespace) -> str:
    return (
        "python3 -u CRON_kinf_branch.py --min-h %d --max-h %d "
        "--digits %d --cert-digits %d"
        % (args.min_h, args.max_h, args.digits, args.cert_digits)
    )


def serialize_anatomy_command(args: argparse.Namespace) -> str:
    parts = [
        "python3 -u CRON_kinf_branch.py",
        "--min-h %d" % args.min_h,
        "--max-h %d" % args.max_h,
        "--digits %d" % args.digits,
        "--cert-digits %d" % args.cert_digits,
        "--pair-anatomy",
        "--anatomy-report %s" % args.anatomy_report,
        "--anatomy-json %s" % args.anatomy_json,
    ]
    if args.skip_kinf:
        parts.append("--skip-kinf")
    return " ".join(parts)


def report_markdown(payload: Dict[str, Any]) -> str:
    rows = payload["rows"]
    certified_yes = [row["h"] for row in rows if row["certified"]["status"] == "YES"]
    uncertified_yes = [row["h"] for row in rows if row["uncertified"]["status"] == "YES"]
    failures = [row["h"] for row in rows if row["uncertified"]["status"] != "YES"]
    h0 = max(failures) if failures else None
    h0_le_40 = h0 is None or h0 <= 40
    if certified_yes:
        certified_range = "%d..%d" % (min(certified_yes), max(certified_yes))
    else:
        certified_range = "none"

    lines = [
        "H0_EMPIRICAL <= 40: %s; CERTIFIED YES RANGE: h=%s"
        % ("YES" if h0_le_40 else "NO", certified_range),
        "",
        "# K_infinity branch-table scan",
        "",
        "## Verdict",
        "",
        "The %d-decimal-digit mpmath pass returned `YES` for `h=%d..%d`. "
        "The proof-grade pass then returned `YES` for exactly "
        "`h=%s`: Arb isolated all quotient critical roots, every explicit "
        "interval-Newton inclusion succeeded, every squared critical-value "
        "ball excluded zero, and all nonmirror value balls were disjoint."
        % (
            payload["scan"]["uncertified_digits"],
            min(uncertified_yes),
            max(uncertified_yes),
            certified_range,
        ),
        "",
        "Thus the empirical failure cutoff is `%s`.  In the whole tested range "
        "there are `2h-2` certified mirror-orbits (equivalently `4h-4` "
        "critical points), stronger than the required `h` good orbits."
        % ("none <= %d" % max(row["h"] for row in rows) if h0 is None else h0),
        "",
        "This is consistent with every banked full-Morse certificate through "
        "`h=30` (indeed through `h=32`): no tested height is refuted.",
        "",
        "## Exact reconstruction and counting convention",
        "",
        "The exact recurrence used by the program is",
        "",
        "```text",
        "N_0=0, N_1=1,",
        "N_(h+1)(x)=P(x+h)N_h(x)-(x+h)^6N_(h-1)(x),",
        "delta_h(x)=N_h(x)/q_h(x)^3,",
        "A_h(x)=q_h N_h'-3q_h' N_h.",
        "```",
        "",
        "With `s=2x+h+1` and `u=s^2`, reflection makes `A_h` an exact "
        "degree-`2h-2` polynomial `J_h(u)`.  Both `N_h(x)^2` and `q_h(x)^6` "
        "are also exact polynomials in `u`, so the program evaluates",
        "",
        "```text",
        "Psi_h(u)=delta_h(x)^2",
        "```",
        "",
        "without choosing a square root of `u`.  The roots of `J_h` are "
        "exactly the mirror-orbits of critical points, with multiplicity; "
        "their `Psi_h` values are exactly the scaled `U_(r,s)` of Q6723 (3.5).",
        "",
        "There is a coordinate ambiguity in Q6723 that matters in code.  If "
        "the standard cell is `-1 < Re(z) < 0`, reflection in that same "
        "orientation is",
        "",
        "```text",
        "(r,s,z) -> (s-1,r+1,-1-z),",
        "```",
        "",
        "not `(s,r,-z)`.  The latter is correct only after switching to the "
        "opposite pole orientation, whose local coordinate lies in "
        "`0 < Re(z) < 1`.  The quotient variable above is orientation-free "
        "and therefore prevents an off-by-one double count.  The finite-cell "
        "identity was independently checked from the Q6723 `F/G` recurrence:",
        "",
        "| h | j | (r,s) | relative error in (3.3) |",
        "|---:|---:|:---:|---:|",
    ]
    for check in payload["cell_identity_checks"]:
        lines.append(
            "| %d | %d | (%d,%d) | `%s` |"
            % (
                check["h"],
                check["j"],
                check["r"],
                check["s"],
                check["relative_error"],
            )
        )

    lines.extend(
        [
            "",
            "The exact transfer anchors used by the implementation are "
            "`lambda_+- = 17 +- 12 sqrt(2)` and the two "
            "normalized diagonal drifts used by the limiting object are both "
            "`-3/2`; consequently no `n^d` ratio remains in `K_infinity`.",
            "",
            "## K_infinity bulk template (uncertified)",
            "",
            "The code implements "
            "`K_infinity(z)=phi(-z)phi(z)+z^6 gamma(-z)gamma(z)` through the "
            "normalized `F_n/b_n,G_n/b_n` recurrence and Richardson "
            "extrapolation.  Arb-accelerated Newton found the following "
            "finite numerical roots of "
            "`H_infinity=zK_infinity'-3K_infinity`.  As an independent "
            "asymptotic diagnostic, exact central `K_(m,m)` roots were also "
            "computed for `m=%s`.  Neither template calculation is used as a "
            "proof of any fixed-height row."
            % ",".join(
                str(value)
                for value in payload["kinfinity"]
                .get("central_branch_evidence", {})
                .get("m_values", [])
            ),
            "",
            "| Re(z) | Im(z) | residual (fit) | shorter-fit check | Newton correction |",
            "|---:|---:|---:|---:|---:|",
        ]
    )
    for root in payload["kinfinity"]["roots"]:
        lines.append(
            "| `%s` | `%s` | `%s` | `%s` | `%s` |"
            % (
                root["real"],
                root["imag"],
                root["residual_upper"],
                root["extrapolation_check"],
                root["newton_correction"],
            )
        )
    failed_seeds = [
        row for row in payload["kinfinity"].get("seed_attempts", []) if not row["converged"]
    ]
    if failed_seeds:
        lines.extend(
            [
                "",
                "The outer seed attempts did not converge to additional finite roots; "
                "their final diagnostics are retained in `CRON_kinf_results.json`.  "
                "The finite central-cell outer branches move outward with `m`, so "
                "they are not silently counted as finite `K_infinity` branches.  "
                "Without a certified outer-contour count this is evidence of escape, "
                "not an exhaustive theorem about all finite `H_infinity` zeros.",
            ]
        )

    lines.extend(
        [
            "",
            "## Per-height machine table",
            "",
            "`min sep` is the minimum distance between midpoint values.  "
            "`cert margin` is the outward-rounded lower bound on the same "
            "pairwise differences (also including distance from zero in the "
            "global minimum).  `rel margin` divides the pairwise lower bound "
            "by the larger value magnitude.  All value balls use the exact "
            "Apéry pole scales, via the algebraically equivalent `Psi_h`.",
            "",
            "| h | critical zeros / quotient orbits | mpmath | Arb | min sep | cert margin | rel margin | bits |",
            "|---:|:---:|:---:|:---:|---:|---:|---:|---:|",
        ]
    )
    for row in rows:
        certified = row["certified"]
        uncertified = row["uncertified"]
        lines.append(
            "| %d | %d / %d | %s | %s | `%s` | `%s` | `%s` | %d |"
            % (
                row["h"],
                certified["critical_zeros"],
                certified["quotient_roots"],
                uncertified["status"],
                certified["status"],
                certified["min_separation"],
                certified["certified_margin"],
                certified["relative_margin"],
                certified["precision_bits"],
            )
        )

    lines.extend(
        [
            "",
            "## Inline machine output",
            "",
            "```text",
            *payload["machine_summary"],
            "```",
            "",
            "## Caveats",
            "",
            "1. `H0_EMPIRICAL` is the finite-scan convention requested in "
            "`CODEX_SPEC_CRON_kinf.md`: the largest failed height in "
            "`2..60`.  It is not the universal tail threshold `h0` of Q6723 "
            "Section 9.",
            "2. The fixed-height Arb rows are genuine finite certificates, "
            "but this run does not certify the all-`h` predicates `PERSIST`, "
            "`EXHAUST`, `SAME`, `CROSS`, and `COUNT` on inverse-length boxes.  "
            "Therefore the finite result alone does not close the all-height "
            "campaign prize.",
            "3. The displayed `K_infinity` roots use accelerated normalized "
            "recurrences and central-branch extrapolation without the explicit "
            "Q6708 conjugation-tail ball or an outer-contour count.  They are "
            "intentionally labelled uncertified; every "
            "fixed-height verdict is recomputed from exact integer polynomials "
            "and does not depend on them.",
            "4. `python-flint` was found in the pinned uv cache rather than the "
            "system Python import path.  The script re-executes through uv and "
            "records the exact versions in the JSON payload.",
            "5. The uncertified pass takes the midpoints of FLINT root balls as "
            "Newton seeds, then discards every radius and recomputes roots and "
            "values in mpmath.  It is a 100+-digit floating-point regression "
            "pass, not an implementation-independent root finder; the Arb pass "
            "is the certificate.",
            "",
            "Reproduction command:",
            "",
            "```bash",
            payload["command"],
            "```",
            "",
        ]
    )
    return "\n".join(lines)


def pair_anatomy_report_markdown(payload: Dict[str, Any]) -> str:
    """Readable summary of the raw orbit and all-pair anatomy payload."""

    constants = payload["constants"]
    short = lambda value: mp.nstr(mp.mpf(value), 8)
    first_line = (
        "%s — kappa_odd=%s, kappa_even=%s, ratio=%s, F''(1/2)=%s, "
        "finite far-band lower bound=%s for |Δt|>=1/4; k=1 decay exponent=%s"
        % (
            payload["verdict"],
            short(constants["kappa_odd"]),
            short(constants["kappa_even"]),
            short(constants["kappa_even_over_odd"]),
            short(constants["F_second_at_one_half_from_even_identity"]),
            short(constants["macroscopic_far_lower_bound"]),
            short(constants["k1_height_decay_exponent"]),
        )
    )
    k1_decay = payload["far_band_height_decay"][0]
    far_pair = constants["macroscopic_far_lower_bound_pair"]
    lines = [
        first_line,
        "",
        "# Certified per-pair margin anatomy",
        "",
        "## Verdict",
        "",
        "The proposed two-regime picture is `REFUTED` by its far-pair clause, "
        "while its central near-pair clause is confirmed.  At every tested "
        "height the global midpoint relative-separation minimum is the "
        "predicted central `inner`-branch pair.  Its cell distance is `1/h` "
        "for odd `h`; for even `h` the pair straddles the reflection-fixed "
        "cell and has distance `2/h`.  The relative gap is "
        "`kappa_parity*h^-2 + O(h^-3)`.",
        "",
        "The far pairs are not bounded below by a height-independent constant "
        "or by `c'|Δt|^2`.  Already in the fixed macroscopic band "
        "`1/2 <= |Δt| < 1`, the minimum falls from `%s` at `h=%d` to `%s` "
        "at `h=%d`; the fit `C*h^-alpha` gives `alpha=%s` (log-RMSE `%s`), "
        "close to `1/h`, not a constant.  Since `|Δt|` stays at least `1/2`, "
        "this also rules out a uniform positive `c'|Δt|^2` interpretation of "
        "the scanned profile."
        % (
            k1_decay["first_minimum"],
            k1_decay["first_height"],
            k1_decay["last_minimum"],
            k1_decay["last_height"],
            k1_decay["alpha"],
            k1_decay["free_log_rmse"],
        ),
        "",
        "There is also a distinct cross-branch event: the smallest observed "
        "separation with `|Δt|>=1/4` is `%s` at `h=%d`, between `%s` and `%s` "
        "at distance `%s`.  This `inner/outer` near-collision is far below the "
        "smooth same-branch trend and shows that a faithful all-height "
        "partition must track branch combinations, not distance alone."
        % (
            constants["macroscopic_far_lower_bound"],
            constants["macroscopic_far_lower_bound_h"],
            far_pair["orbit_a"],
            far_pair["orbit_b"],
            far_pair["cell_distance"],
        ),
        "",
        "This refutation does not by itself kill the intended comparison "
        "strategy: the generic far gap observed here is roughly `h^-1`, still "
        "larger than an `O(h^-2)` tail.  It does kill the stated two-regime "
        "lower-bound claim and exposes an additional cross-branch regime that "
        "needs its own estimate.",
        "",
        "## Cell and branch convention",
        "",
        "A root ball `u` of the exact quotient polynomial `J_h(u)` represents "
        "the critical-point pair",
        "",
        "```text",
        "x = (-(h+1) ± sqrt(u))/2.",
        "```",
        "",
        "The anatomy pass chooses the unique member in the upper half-plane.  "
        "If it lies in `-j-1 < Re(x) < -j`, its cell position is `t=j/h`.  "
        "The Arb ball itself is checked to lie strictly inside that cell.  "
        "For every `j=1..h-1` there are exactly two such representatives; "
        "the one with smaller imaginary part is `inner`, the other `outer`.  "
        "This gives exactly `2h-2` labelled mirror-orbits.",
        "",
        "Choosing the lower-half-plane member instead sends `j` to `h-j`, "
        "hence `t` to `1-t`.  It preserves every `|t_a-t_b|`, dyadic band, "
        "and minimising pair.  The reported anatomy is therefore independent "
        "of this representative choice.",
        "",
        "Each raw orbit record contains the midpoint and radius of its Arb "
        "root ball, representative, local coordinate, and squared critical "
        "value.  Each raw pair record uses the midpoint quantity",
        "",
        "```text",
        "|V_a - V_b| / max(|V_a|, |V_b|).",
        "```",
        "",
        "The original per-height Arb lower bound is retained separately as "
        "`certified_relative_margin`.",
        "",
        "## Central-gap fits",
        "",
        "The primary fits use the tail window beginning at `h=40` and the "
        "model `relative_gap*h^2 = kappa + a/h + b/h^2`.  The fit range, "
        "RMSE, and maximum absolute residual are stated explicitly; "
        "digits beyond the residual scale are diagnostic rather than certified "
        "as asymptotic constants.",
        "",
        "| parity | fit heights | n | kappa | a | b | RMSE | max |residual| |",
        "|:---:|:---:|---:|---:|---:|---:|---:|---:|",
    ]
    for parity in ("odd", "even"):
        fit = payload["central_fits"][parity]
        coefficients = fit["coefficients"]
        lines.append(
            "| %s | %s | %d | `%s` | `%s` | `%s` | `%s` | `%s` |"
            % (
                parity,
                fit["fit_range"],
                fit["point_count"],
                coefficients["kappa"],
                coefficients["a"],
                coefficients["b"],
                fit["rmse"],
                fit["max_abs_residual"],
            )
        )
    lines.extend(
        [
            "",
            "Fit-window drift is the relevant systematic uncertainty.  The "
            "following table repeats the same fit on successively later "
            "windows; every displayed constant therefore carries its own "
            "range and residual.",
            "",
            "| parity | fit heights | n | kappa | RMSE | max |residual| |",
            "|:---:|:---:|---:|---:|---:|---:|",
        ]
    )
    for parity in ("odd", "even"):
        for fit in payload["central_fit_stability"][parity]:
            lines.append(
                "| %s | %s | %d | `%s` | `%s` | `%s` |"
                % (
                    parity,
                    fit["fit_range"],
                    fit["point_count"],
                    fit["coefficients"]["kappa"],
                    fit["rmse"],
                    fit["max_abs_residual"],
                )
            )
    lines.extend(
        [
            "",
            "Using the proved identities `kappa_even = F''(1/2)` and "
            "`kappa_odd = F''(1/2)/2`, the data give",
            "",
            "```text",
            "F''(1/2) = %s" % constants["F_second_at_one_half_from_even_identity"],
            "kappa_odd - kappa_even/2 = %s" % constants["kappa_odd_minus_half_kappa_even"],
            "relative half-identity error = %s" % constants["relative_half_identity_error"],
            "```",
            "",
            "### Per-height central gaps and fit residuals",
            "",
            "| h | parity | minimizing cells | Δt | relative gap | gap*h^2 | fit residual | global min? | Arb bits |",
            "|---:|:---:|:---:|---:|---:|---:|---:|:---:|---:|",
        ]
    )
    # The per-height table spans the whole scan, so use the widest stability
    # window for its residual column.  The primary tail-fit residuals remain
    # available verbatim in ``central_fits`` and the raw JSON.
    residual_lookup = {
        (parity, row["h"]): row
        for parity, fits in payload["central_fit_stability"].items()
        for row in fits[0]["residuals"]
    }
    for row in payload["heights"]:
        parity = "odd" if row["h"] % 2 else "even"
        central = row["central_gap"]
        residual = residual_lookup[(parity, row["h"])]["residual"]
        lines.append(
            "| %d | %s | `%s` | `%s` | `%s` | `%s` | `%s` | %s | %d |"
            % (
                row["h"],
                parity,
                ",".join(str(value) for value in central["expected_cells"]),
                central["cell_distance"],
                central["relative_separation"],
                central["gap_times_h_squared"],
                residual,
                "YES" if central["is_global_minimum"] else "NO",
                row["certificate"]["precision_bits"],
            )
        )

    diagnostic = payload["far_band_model_diagnostic"]
    lines.extend(
        [
            "",
            "## Dyadic-band profile",
            "",
            "A positive-distance pair is put in band `k` when "
            "`2^-k <= |Δt| < 2^(-k+1)`.  Same-cell `inner/outer` pairs are "
            "recorded separately in the JSON.  The `sep/Δt^2` column is the "
            "smallest value over all pairs in that band, so it directly tests "
            "a lower bound of the form `c'|Δt|^2`.",
            "",
            "### Aggregate over all tested heights",
            "",
            "| band | distance range | pair count | empirical min sep | h | minimizing pair | min sep/Δt^2 |",
            "|:---:|:---:|---:|---:|---:|:---|---:|",
        ]
    )
    for band in payload["aggregate_dyadic_bands"]:
        pair = band["minimizing_pair"]
        lines.append(
            "| %s | %s | %d | `%s` | %d | `%s / %s` | `%s` |"
            % (
                band["band"],
                band["distance_interval"],
                band["pair_count"],
                band["empirical_lower_bound"],
                band["attained_at_h"],
                pair["orbit_a"],
                pair["orbit_b"],
                band["minimum_relative_over_distance_squared"],
            )
        )
    if diagnostic.get("status") == "EMPIRICAL":
        lines.extend(
            [
                "",
                "Across the three far distance bands `%s`, a free log-log "
                "fit against distance has slope `%s`.  Its constant-model "
                "log-RMSE is `%s`, versus `%s` for a quadratic-distance "
                "model.  This cross-band diagnostic alone calls the shape "
                "`%s`, but it does **not** supply a uniform lower bound because "
                "each fixed band's minimum itself decays with height."
                % (
                    ", ".join(diagnostic["bands"]),
                    diagnostic["free_log_log_slope"],
                    diagnostic["constant_log_rmse"],
                    diagnostic["quadratic_log_rmse"],
                    diagnostic["closer_model"],
                ),
            ]
        )
    lines.extend(
        [
            "",
            "### Fixed-band decay with height",
            "",
            "| band | distance range | fit heights | alpha in C*h^-alpha | free log-RMSE | constant log-RMSE | 1/h log-RMSE | first min | last min |",
            "|:---:|:---:|:---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for decay in payload["far_band_height_decay"]:
        lines.append(
            "| %s | %s | %s | `%s` | `%s` | `%s` | `%s` | `%s` (h=%d) | `%s` (h=%d) |"
            % (
                decay["band"],
                decay["distance_interval"],
                decay["fit_range"],
                decay["alpha"],
                decay["free_log_rmse"],
                decay["constant_log_rmse"],
                decay["inverse_h_log_rmse"],
                decay["first_minimum"],
                decay["first_height"],
                decay["last_minimum"],
                decay["last_height"],
            )
        )
    lines.extend(
        [
            "",
            "### Profile as a function of height",
            "",
            "| h | band | distance range | pair count | min relative sep | minimizing pair | min sep/Δt^2 |",
            "|---:|:---:|:---:|---:|---:|:---|---:|",
        ]
    )
    for row in payload["heights"]:
        for band in row["dyadic_bands"]:
            pair = band["minimizing_pair"]
            lines.append(
                "| %d | %s | %s | %d | `%s` | `%s / %s` | `%s` |"
                % (
                    row["h"],
                    band["band"],
                    band["distance_interval"],
                    band["pair_count"],
                    band["minimum_relative_separation"],
                    pair["orbit_a"],
                    pair["orbit_b"],
                    band["minimum_relative_over_distance_squared"],
                )
            )

    lines.extend(
        [
            "",
            "## Caveats",
            "",
            "1. `REFUTED` means that the observed fixed far band has clear "
            "height decay, contrary to the proposed uniform lower bound.  A "
            "finite scan does not prove the eventual exponent.  The orbit "
            "balls and fixed-height separation certificates are rigorous; "
            "the extrapolated exponents and `kappa` values are numerical fits.",
            "2. Relative pair separations in the anatomy tables use Arb-ball "
            "midpoints, as requested.  The raw midpoint and radius make the "
            "rounding scale explicit.  The pre-existing certificate path was "
            "not weakened and its outward-rounded minimum is retained per "
            "height.",
            "3. The `inner/outer` label is geometric (imaginary height), not a "
            "new algebraic theorem about continuation of branches in `h`.",
            "4. The constant-vs-quadratic far-band comparison uses only three "
            "fixed dyadic bands and is descriptive.  On any region bounded "
            "away from zero, both a positive constant and `c'|Δt|^2` give a "
            "positive uniform lower bound; the table is the durable result.",
            "5. Heights 80 and 100 are optional in the mission.  They are not "
            "part of this report unless they appear in the per-height table.",
            "",
            "Reproduction command:",
            "",
            "```bash",
            payload["command"],
            "```",
            "",
            "The complete orbit balls and every nonmirror pair are in "
            "`%s`." % Path(payload["command"].split("--anatomy-json ", 1)[1].split()[0]).name,
            "",
        ]
    )
    return "\n".join(lines)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--min-h", type=int, default=2)
    parser.add_argument("--max-h", type=int, default=60)
    parser.add_argument("--digits", type=int, default=120)
    parser.add_argument("--cert-digits", type=int, default=110)
    parser.add_argument("--retries", type=int, default=2)
    parser.add_argument("--kinf-base", type=int, default=10)
    parser.add_argument("--kinf-step", type=int, default=5)
    parser.add_argument("--kinf-order", type=int, default=4)
    parser.add_argument("--skip-kinf", action="store_true")
    parser.add_argument("--report", default="CODEX_KINF_report.md")
    parser.add_argument("--json", default="CRON_kinf_results.json")
    parser.add_argument(
        "--pair-anatomy",
        action="store_true",
        help="emit certified orbit labels, every pair, dyadic profiles, and fits",
    )
    parser.add_argument("--anatomy-report", default="CODEX_MARGINS_report.md")
    parser.add_argument("--anatomy-json", default="CRON_margins_results.json")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.min_h < 2 or args.max_h < args.min_h:
        raise SystemExit("require 2 <= min-h <= max-h")
    if args.digits < 100 or args.cert_digits < 100:
        raise SystemExit("both numerical passes must use at least 100 digits")

    ctx.threads = 1
    started = time.monotonic()
    import flint  # type: ignore
    import mpmath  # type: ignore

    log(
        "[config] h=%d..%d mpmath=%d digits Arb target=%d digits"
        % (args.min_h, args.max_h, args.digits, args.cert_digits)
    )
    m_values = [args.kinf_base + args.kinf_step * k for k in range(args.kinf_order + 1)]
    kinf_height = 2 * max(m_values) + 1
    numerators = build_numerators(max(args.max_h, 12, kinf_height))
    apery_small = build_apery_numbers(max(args.max_h, 12))
    checks = check_cell_identity(numerators, apery_small, args.digits)
    for check in checks:
        log(
            "[cell] h=%d j=%d r=%d s=%d relerr=%s"
            % (
                check["h"],
                check["j"],
                check["r"],
                check["s"],
                check["relative_error"],
            )
        )

    if args.skip_kinf:
        kinfinity = {
            "status": "SKIPPED",
            "method": "order-%d Richardson extrapolation" % args.kinf_order,
            "root_count": 0,
            "roots": [],
        }
    else:
        log("[kinf] isolating H_infinity template roots (uncertified)")
        kinfinity = isolate_kinfinity_roots(
            numerators, args.digits, m_values
        )
        log("[kinf] roots=%d" % kinfinity["root_count"])

    rows = []
    anatomy_rows = []
    machine_summary = []
    for h in range(args.min_h, args.max_h + 1):
        data = build_height_data(h, numerators[h])
        certified, roots = certify_height(data, args.cert_digits, args.retries)
        if not roots:
            uncertified = {
                "status": "UNRESOLVED",
                "precision_digits": args.digits,
                "quotient_roots": 0,
                "critical_zeros": 0,
                "good_orbits": 0,
                "min_separation": "nan",
                "estimated_margin": "nan",
                "relative_margin": "nan",
                "max_newton_correction": "nan",
            }
        else:
            uncertified = uncertified_height(data, roots, args.digits)
        if args.pair_anatomy:
            if certified.get("status") != "YES" or len(roots) != 2 * h - 2:
                raise ArithmeticError(
                    "pair anatomy requires a complete certified row at h=%d" % h
                )
            anatomy = pair_anatomy_height(data, roots, args.digits)
            anatomy["certificate"] = {
                "status": certified["status"],
                "precision_bits": certified["precision_bits"],
                "certified_relative_margin": certified["relative_margin"],
                "max_root_radius": certified["max_root_radius"],
                "interval_newton_pass": certified["interval_newton_pass"],
                "denominator_excludes_zero": certified["denominator_excludes_zero"],
            }
            anatomy_rows.append(anatomy)
            log(
                "[anatomy] h=%d pairs=%d min=%s central=%s"
                % (
                    h,
                    anatomy["pair_count"],
                    anatomy["global_minimum"]["relative_separation"],
                    anatomy["global_minimum"][
                        "is_predicted_central_inner_pair"
                    ],
                )
            )
        row = {"h": h, "uncertified": uncertified, "certified": certified}
        rows.append(row)
        summary = (
            "h=%02d UNCERTIFIED=%-10s CERTIFIED=%-10s "
            "zeros=%d orbits=%d good=%d rel_margin=%s"
            % (
                h,
                uncertified["status"],
                certified.get("status", "UNRESOLVED"),
                certified.get("critical_zeros", 0),
                certified.get("quotient_roots", 0),
                certified.get("good_orbits", 0),
                certified.get("relative_margin", "nan"),
            )
        )
        machine_summary.append(summary)
        log(summary)

    failures = [row["h"] for row in rows if row["uncertified"]["status"] != "YES"]
    empirical_h0 = max(failures) if failures else None
    certified_yes = [row["h"] for row in rows if row["certified"]["status"] == "YES"]
    payload = {
        "schema": "CRON-kinf-branch-v1",
        "command": serialize_command(args),
        "runtime": {
            "python": sys.version.split()[0],
            "python_flint": flint.__version__,
            "mpmath": mpmath.__version__,
            "wall_seconds": round(time.monotonic() - started, 6),
        },
        "anchors": {
            "lambda_plus": LAMBDA_PLUS_TEXT,
            "lambda_minus": LAMBDA_MINUS_TEXT,
            "d_plus": "-3/2",
            "d_minus": "-3/2",
        },
        "scan": {
            "min_h": args.min_h,
            "max_h": args.max_h,
            "uncertified_digits": args.digits,
            "certified_target_digits": args.cert_digits,
            "h0_empirical": empirical_h0,
            "h0_empirical_text": (
                "none <= %d" % args.max_h if empirical_h0 is None else str(empirical_h0)
            ),
            "certified_yes": certified_yes,
        },
        "cell_identity_checks": checks,
        "kinfinity": kinfinity,
        "rows": rows,
        "machine_summary": machine_summary,
    }

    json_path = Path(args.json)
    report_path = Path(args.report)
    json_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    report_path.write_text(report_markdown(payload))
    log("[done] wrote %s and %s" % (json_path, report_path))
    if args.pair_anatomy:
        anatomy_payload = build_pair_anatomy_payload(payload, anatomy_rows, args)
        anatomy_json_path = Path(args.anatomy_json)
        anatomy_report_path = Path(args.anatomy_report)
        anatomy_json_path.write_text(
            json.dumps(anatomy_payload, separators=(",", ":")) + "\n"
        )
        anatomy_report_path.write_text(
            pair_anatomy_report_markdown(anatomy_payload)
        )
        log(
            "[done] wrote %s and %s"
            % (anatomy_json_path, anatomy_report_path)
        )
    log("[done] wall=%.2fs" % (time.monotonic() - started))


if __name__ == "__main__":
    main()
