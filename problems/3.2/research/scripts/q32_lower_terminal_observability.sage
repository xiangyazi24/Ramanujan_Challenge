#!/usr/bin/env sage
"""Exact lower-packet observability after the top-half Freshman reduction.

For ``n = p + r`` and the terminal nodes

    d_t = n - 1 - t = p + (r - 1 - t),

the unique Freshman decomposition (away from the bounded borrow strip)
reduces the actual ``Y/W`` raw vector modulo ``p`` to the vector
constructed below.  This script computes the gcd of its 2-by-2 minors.

The calculation is diagnostic until the unique-decomposition lemma and
the resulting formulas are written out in the accompanying notes.
"""

from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

load(str(HERE / "q32_doubled_period_gauge_audit.sage"))

from q32_cartier_packet_audit import LAMBDA, coefficient


def weighted_correlation(moment, exponent, modulus=None):
    """Return sum_k [k]Lambda * [exponent*k]Lambda^moment."""

    total = 0
    for point, weight in LAMBDA.items():
        total += weight * coefficient(
            moment,
            exponent * point[0],
            exponent * point[1],
            exponent * point[2],
            modulus,
        )
        if modulus:
            total %= modulus
    return ZZ(total if modulus is None else total % modulus)


def lower_raw_vector(residue, terminal_offset, modulus=None):
    """The reduced terminal vector at d=p+r-1-terminal_offset."""

    exponent = residue - 1 - terminal_offset
    y_value = weighted_correlation(
        residue - 1, exponent, modulus
    )
    z_value = 0
    for shift in range(origin_operator.order() + 1):
        scalar = ZZ(
            integer_coefficients[shift](n=residue)
        )
        z_value += scalar * (
            weighted_correlation(
                residue + shift, exponent - 1, modulus
            )
            + weighted_correlation(
                residue + shift, exponent + 1, modulus
            )
        )
        if modulus:
            z_value %= modulus
    w_value = z_value - ZZ(multiplier(n=residue)) * y_value
    if modulus:
        w_value %= modulus
    return vector(
        ZZ,
        (
            y_value if modulus is None else y_value % modulus,
            w_value if modulus is None else w_value % modulus,
        ),
    )


def observability_content(residue, width=7):
    """Return the raw vectors, all minors, and their exact gcd."""

    vectors = [
        lower_raw_vector(residue, offset)
        for offset in range(width + 1)
    ]
    minors = {
        (left, right):
        vectors[left][0] * vectors[right][1]
        - vectors[left][1] * vectors[right][0]
        for left in range(width + 1)
        for right in range(left + 1, width + 1)
    }
    return vectors, minors, gcd(minors.values())


def origin_content(index):
    return gcd(
        [
            ZZ(multiplier(n=index)),
            *[
                ZZ(polynomial(n=index))
                for polynomial in integer_coefficients
            ],
        ]
    )


def origin_bezout_constant():
    """Return a fixed integer containing every pointwise origin content."""

    polynomials = [*integer_coefficients, multiplier]
    bezout = [R.one()] + [R.zero()] * (len(polynomials) - 1)
    running_gcd = polynomials[0]
    for index in range(1, len(polynomials)):
        new_gcd, left, right = running_gcd.xgcd(
            polynomials[index]
        )
        bezout = [left * entry for entry in bezout]
        bezout[index] += right
        running_gcd = new_gcd
    assert running_gcd == 1
    assert sum(
        bezout[index] * polynomials[index]
        for index in range(len(polynomials))
    ) == 1
    denominator = lcm(
        coefficient.denominator()
        for polynomial in bezout
        for coefficient in polynomial.list()
    )
    return ZZ(denominator), bezout


if "--origin-bezout" in sys.argv:
    constant, _ = origin_bezout_constant()
    print("ORIGIN_BEZOUT_CONSTANT_BITS", constant.nbits())
    if "--factor-origin-bezout" in sys.argv:
        print("ORIGIN_BEZOUT_CONSTANT_FACTOR", factor(constant))


residues = tuple(range(12, 61))
if "--short" in sys.argv:
    residues = (12, 13, 17, 19, 23, 29, 31, 37)

for residue in residues:
    _, minors, content = observability_content(residue)
    local_origin_content = origin_content(residue)
    quotient = content // gcd(content, local_origin_content)
    predicted_support = ZZ(2 * 3 * 5 * 7) * prod(
        residue + shift for shift in range(2, 8)
    )
    unsupported_factors = tuple(
        (ZZ(prime), exponent)
        for prime, exponent in factor(quotient)
        if predicted_support % prime
    )
    large_factors = tuple(
        (ZZ(prime), exponent)
        for prime, exponent in factor(quotient)
        if prime > residue + 11
    )
    smallest = sorted(
        minors,
        key=lambda pair: abs(minors[pair]).nbits(),
    )[:3]
    print(
        "LOWER_OBSERVABILITY",
        "r", residue,
        "bits", abs(content).nbits(),
        "factor", factor(content),
        "origin_content", factor(local_origin_content),
        "quotient_factor", factor(quotient),
        "unsupported_factors", unsupported_factors,
        "prime_above_unique_range", large_factors,
        "smallest_minors", tuple(
            (pair, abs(minors[pair]).nbits())
            for pair in smallest
        ),
    )
    assert not unsupported_factors

print("Q32_LOWER_TERMINAL_OBSERVABILITY=PASS")
