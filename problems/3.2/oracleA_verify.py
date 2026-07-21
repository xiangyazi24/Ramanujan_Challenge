#!/usr/bin/env python3
"""Exact finite checks for the analytic-oracle identities.

The checks here separate four facts which are easy to conflate.

* With ``e_p(t) = exp(2*pi*i*t/p)`` and the plus-sign transform
  ``F_p(a) = sum_{r in Z_p} e_p(a*r)``, reflection under
  ``r -> p-1-r`` gives

      conjugate(F_p(a)) = e_p(a) F_p(a),
      e_p(a/2) F_p(a) is real,

  where ``1/2`` in ``e_p(a/2)`` is the inverse of 2 modulo p.  In
  ordinary real arguments, the real rotation is therefore
  ``e(-a*(p-1)/(2*p))``, not the same expression with a plus sign.

* If ``Z_p = {r, p-1-r}``, with ``r < (p-1)/2`` and
  ``h = p-1-2*r``, then

      F_p(a) = 2 e_p(a*(p-1)/2) cos(pi*a*h/p)
             = 2 e_p(a*(3*p-1)/2) cos(pi*a*h/p).

* For the top prime block ``N < p <= 2*N``, put ``s=p-1-r`` and
  ``L=2*N-p``.  A residue u has its unique hit in ``I_N=(N,2*N]`` at
  ``u+p`` if ``u <= L``, no hit if ``L < u <= N``, and at ``u`` if
  ``N < u < p``.  Consequently the two positions with centre
  ``(3*p-1)/2`` are the actual doublet hits exactly in the low--low
  case ``s <= L``.  Other doublets need boundary terms.

* A circular short interval of length ``delta <= 1/2`` has the exact
  pair expansion

      integral |sum_{x<alpha_j<=x+delta} c_j|^2 dx
        = sum_{j,k} c_j conjugate(c_k)
          max(delta-||alpha_j-alpha_k||, 0).

  A finite spectral sum, together with an explicit positive tail bound,
  independently checks the full Fejer expansion, including the integer
  mode ``k=0`` and the removable definition of ``gamma_M(0)``.

* For distinct primes, the script exhausts the least-residue bijection
  ``(a,b) <-> k=<a*q-b*p>_(pq)`` and compares the two Fourier products as
  exact multisets of ``pq``-th roots of unity on actual nonempty zero sets.

* Small reflected anchored families are checked directly for reflection,
  distinctness, nonadjacency, their common hit at ``2*N``, both ``D(2*N)``
  and ``D(0)``, finite Fourier Parseval for ``Q``, and
  ``gamma_(4*N)(2*N)``.

The recurrence checks use both the divided Apéry recurrence and the
division-free recurrence for ``(n!)^3 b_n`` and compare them with the raw
``data_zp_pairs.bin`` bank.  Numerical evaluations of roots of unity are
secondary guards; the phase, lift, Fourier-inversion, and kernel checks are
performed with integer counters or ``Fraction`` arithmetic.

This script verifies identities and transcription only.  It does not treat a
finite computation as a proof of the missing cross-prime large-sieve bound.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction
from math import cos, fsum, gcd, isqrt, log, pi, sin
from pathlib import Path
import cmath
import struct
import sys


HERE = Path(__file__).resolve().parent
PAIR_BANK = HERE / "data_zp_pairs.bin"
TOP_N = 1024


class CheckFailure(RuntimeError):
    """A check failed even when Python is run with optimization enabled."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise CheckFailure(message)


def primes_up_to(limit: int) -> list[int]:
    sieve = bytearray(b"\x01") * (limit + 1)
    if limit >= 0:
        sieve[0] = 0
    if limit >= 1:
        sieve[1] = 0
    for prime in range(2, isqrt(limit) + 1):
        if sieve[prime]:
            sieve[prime * prime : limit + 1 : prime] = b"\x00" * (
                (limit - prime * prime) // prime + 1
            )
    return [number for number in range(2, limit + 1) if sieve[number]]


def apery_zeros_divided(prime: int) -> tuple[int, ...]:
    """Compute b_0,...,b_{p-1} modulo p from the original recurrence."""

    values = [1 % prime, 5 % prime]
    for index in range(1, prime - 1):
        coefficient = (
            34 * index**3 + 51 * index**2 + 27 * index + 5
        ) % prime
        numerator = (
            coefficient * values[index] - index**3 * values[index - 1]
        ) % prime
        denominator = pow(index + 1, 3, prime)
        values.append(numerator * pow(denominator, -1, prime) % prime)
    return tuple(index for index, value in enumerate(values) if value == 0)


def apery_zeros_cleared(prime: int) -> tuple[int, ...]:
    """Compute the same zero set using Y_n=(n!)^3 b_n."""

    previous, current = 1 % prime, 5 % prime
    zeros: list[int] = []
    if previous == 0:
        zeros.append(0)
    if current == 0:
        zeros.append(1)
    for index in range(1, prime - 1):
        coefficient = (
            34 * index**3 + 51 * index**2 + 27 * index + 5
        ) % prime
        following = (
            coefficient * current - pow(index, 6, prime) * previous
        ) % prime
        previous, current = current, following
        if current == 0:
            zeros.append(index + 1)
    return tuple(zeros)


def load_pair_bank() -> dict[int, tuple[int, ...]]:
    raw = PAIR_BANK.read_bytes()
    require(len(raw) % 8 == 0, "pair bank length is not divisible by 8")
    pairs = list(struct.iter_unpack("<II", raw))
    require(pairs == sorted(pairs), "pair bank is not sorted by (p,r)")
    require(len(pairs) == len(set(pairs)), "pair bank contains duplicates")
    grouped: dict[int, list[int]] = defaultdict(list)
    for prime, residue in pairs:
        require(7 <= prime <= 2_000_000, f"bank prime out of range: {prime}")
        require(residue < prime, f"bank residue out of range: ({prime},{residue})")
        grouped[prime].append(residue)
    for prime, residues in grouped.items():
        residue_set = set(residues)
        require(
            all(prime - 1 - residue in residue_set for residue in residues),
            f"bank reflection failure at p={prime}",
        )
    return {prime: tuple(residues) for prime, residues in grouped.items()}


def check_recurrence_and_bank() -> str:
    bank = load_pair_bank()
    primes = [prime for prime in primes_up_to(2 * TOP_N) if prime >= 7]
    active = 0
    zero_count = 0
    zero_sets: dict[int, tuple[int, ...]] = {}
    for prime in primes:
        divided = apery_zeros_divided(prime)
        cleared = apery_zeros_cleared(prime)
        require(divided == cleared, f"recurrence mismatch at p={prime}")
        require(cleared == bank.get(prime, ()), f"bank mismatch at p={prime}")
        require(
            all(prime - 1 - residue in cleared for residue in cleared),
            f"recurrence reflection failure at p={prime}",
        )
        active += bool(cleared)
        zero_count += len(cleared)
        zero_sets[prime] = cleared

    expected = {
        7: (),
        11: (5,),
        17: (3, 13),
        181: (19, 47, 133, 161),
        379: (85, 152, 171, 207, 226, 293),
        1087: (77, 1009),
        1097: (200, 896),
        1213: (58, 1154),
        1559: (755, 803),
        1657: (586, 1070),
    }
    for prime, wanted in expected.items():
        require(zero_sets[prime] == wanted, f"fixed zero-set mismatch at p={prime}")

    # Save data for the later checks without making them trust the binary bank.
    check_recurrence_and_bank.zero_sets = zero_sets  # type: ignore[attr-defined]
    return (
        f"two recurrences and pair bank agree for {len(primes)} primes "
        f"through {2 * TOP_N} ({active} active, {zero_count} zeros)"
    )


def ep(prime: int, exponent: int) -> complex:
    return cmath.exp(2j * pi * (exponent % prime) / prime)


def fourier_value(prime: int, zeros: tuple[int, ...], frequency: int) -> complex:
    return sum(ep(prime, frequency * residue) for residue in zeros)


def actual_zero_sets() -> dict[int, tuple[int, ...]]:
    zero_sets = getattr(check_recurrence_and_bank, "zero_sets", None)
    require(zero_sets is not None, "recurrence check must run first")
    return zero_sets


def check_palindromic_phase() -> str:
    zero_sets = actual_zero_sets()
    frequency_checks = 0
    for prime, zeros in zero_sets.items():
        inverse_two = (prime + 1) // 2
        centre = (prime - 1) // 2
        require((2 * inverse_two) % prime == 1, f"bad inverse of 2 at p={prime}")
        require(
            inverse_two % prime == (-centre) % prime,
            f"modular/ordinary phase mismatch at p={prime}",
        )
        for frequency in range(prime):
            # Exact termwise form of conjugate(F)=e_p(a)F.
            conjugate_exponents = Counter(
                (-frequency * residue) % prime for residue in zeros
            )
            shifted_exponents = Counter(
                (frequency + frequency * residue) % prime for residue in zeros
            )
            require(
                conjugate_exponents == shifted_exponents,
                f"conjugation phase failure at p={prime}, a={frequency}",
            )

            # After multiplying by e_p(a/2), the exponent multiset is closed
            # under negation, so the corresponding root-of-unity sum is real.
            rotated = Counter(
                (frequency * (residue + inverse_two)) % prime
                for residue in zeros
            )
            require(
                rotated == Counter({(-key) % prime: value for key, value in rotated.items()}),
                f"real-phase exponent failure at p={prime}, a={frequency}",
            )

            value = ep(prime, frequency * inverse_two) * fourier_value(
                prime, zeros, frequency
            )
            require(
                abs(value.imag) <= 5e-11 * max(1, len(zeros)),
                f"numerical phase guard failed at p={prime}, a={frequency}",
            )
            frequency_checks += 1

    # The plus ordinary phase printed earlier in proof.tex is not the real
    # rotation for the plus-sign Fourier convention.  This actual doublet is
    # an explicit guard against silently changing the sign.
    prime, frequency = 17, 1
    zeros = zero_sets[prime]
    centre = (prime - 1) // 2
    wrong_rotation = ep(prime, frequency * centre) * fourier_value(
        prime, zeros, frequency
    )
    require(
        abs(wrong_rotation.imag) > 0.1,
        "the wrong-sign phase counterexample at (p,a)=(17,1) disappeared",
    )
    return (
        f"exact phase checked at {frequency_checks} (p,a) pairs; "
        "p=17 rejects the wrong ordinary plus sign"
    )


def check_doublet_fourier_formula() -> str:
    zero_sets = actual_zero_sets()
    doublets = 0
    frequency_checks = 0
    for prime, zeros in zero_sets.items():
        if len(zeros) != 2:
            continue
        lower, upper = zeros
        centre = (prime - 1) // 2
        inverse_two = (prime + 1) // 2
        gap = prime - 1 - 2 * lower
        phase = (3 * prime - 1) // 2
        require(lower < centre < upper, f"unordered doublet at p={prime}")
        require(upper == prime - 1 - lower, f"bad reflected doublet at p={prime}")
        require(gap == upper - lower > 0 and gap % 2 == 0, f"bad gap at p={prime}")
        require(phase % prime == centre, f"3p phase mismatch at p={prime}")

        for frequency in range(prime):
            lhs = Counter(
                ((frequency * lower) % prime, (frequency * upper) % prime)
            )
            rhs = Counter(
                (
                    (frequency * (phase - gap // 2)) % prime,
                    (frequency * (phase + gap // 2)) % prime,
                )
            )
            require(lhs == rhs, f"doublet exponent identity failed at p={prime}, a={frequency}")

            rotated = Counter(
                (frequency * (residue + inverse_two)) % prime
                for residue in zeros
            )
            wanted_rotated = Counter(
                (
                    (-frequency * (gap // 2)) % prime,
                    (frequency * (gap // 2)) % prime,
                )
            )
            require(
                rotated == wanted_rotated,
                f"rotated cosine exponents failed at p={prime}, a={frequency}",
            )

            direct = fourier_value(prime, zeros, frequency)
            factored = (
                2
                * ep(prime, frequency * phase)
                * cos(pi * frequency * gap / prime)
            )
            require(
                abs(direct - factored) <= 1e-10,
                f"numerical doublet formula failed at p={prime}, a={frequency}",
            )
            frequency_checks += 1
        doublets += 1
    require(doublets > 0, "no actual doublets were checked")
    return (
        f"doublet cosine/phase identity checked for {doublets} primes "
        f"and {frequency_checks} frequencies"
    )


def residue_hits_top_block(scale: int, prime: int, residue: int) -> tuple[int, ...]:
    """All m in (scale,2*scale] congruent to residue modulo prime."""

    first_multiplier = max(0, (scale - residue) // prime + 1)
    last_multiplier = (2 * scale - residue) // prime
    return tuple(
        residue + multiplier * prime
        for multiplier in range(first_multiplier, last_multiplier + 1)
    )


def lift_region(scale: int, prime: int, residue: int) -> str:
    boundary = 2 * scale - prime
    if residue <= boundary:
        return "low"
    if residue <= scale:
        return "middle"
    return "high"


def check_top_block_hit_positions() -> str:
    zero_sets = actual_zero_sets()
    scale = TOP_N
    top_primes = [
        prime for prime in zero_sets if scale < prime <= 2 * scale
    ]

    # First check the complete, residue-by-residue lift partition.  This is
    # independent of which residues happen to be Apéry zeros.
    residue_checks = 0
    for prime in top_primes:
        boundary = 2 * scale - prime
        for residue in range(prime):
            region = lift_region(scale, prime, residue)
            if region == "low":
                expected = (residue + prime,)
                require(residue <= boundary, "low-region transcription error")
            elif region == "middle":
                expected = ()
                require(boundary < residue <= scale, "middle-region transcription error")
            else:
                expected = (residue,)
                require(scale < residue < prime, "high-region transcription error")
            require(
                residue_hits_top_block(scale, prime, residue) == expected,
                f"top-block lift failure at N={scale}, p={prime}, r={residue}",
            )
            residue_checks += 1

    cases: Counter[tuple[str, str]] = Counter()
    low_low_doublets = 0
    for prime in top_primes:
        zeros = zero_sets[prime]
        if len(zeros) != 2:
            continue
        lower, upper = zeros
        gap = prime - 1 - 2 * lower
        lower_region = lift_region(scale, prime, lower)
        upper_region = lift_region(scale, prime, upper)
        case = (lower_region, upper_region)
        cases[case] += 1
        actual_hits = tuple(
            sorted(
                residue_hits_top_block(scale, prime, lower)
                + residue_hits_top_block(scale, prime, upper)
            )
        )
        proposed = (prime + lower, prime + upper)

        if case == ("low", "low"):
            require(actual_hits == proposed, f"low-low hits failed at p={prime}")
            require(sum(actual_hits) == 3 * prime - 1, f"low-low centre failed at p={prime}")
            require(actual_hits[1] - actual_hits[0] == gap, f"low-low gap failed at p={prime}")
            low_low_doublets += 1
        elif case == ("low", "middle"):
            require(actual_hits == (prime + lower,), f"low-middle hits failed at p={prime}")
        elif case == ("low", "high"):
            require(actual_hits == (upper, prime + lower), f"low-high hits failed at p={prime}")
            require(sum(actual_hits) == 2 * prime - 1, f"low-high centre failed at p={prime}")
            require(
                actual_hits[1] - actual_hits[0] == 2 * lower + 1 == prime - gap,
                f"low-high gap failed at p={prime}",
            )
        elif case == ("middle", "middle"):
            require(actual_hits == (), f"middle-middle hits failed at p={prime}")
        elif case == ("middle", "high"):
            require(actual_hits == (upper,), f"middle-high hits failed at p={prime}")
        else:
            raise CheckFailure(f"unexpected doublet lift case {case} at p={prime}")

        require(
            (actual_hits == proposed) == (case == ("low", "low")),
            f"claimed iff-domain for (3p-1)/2 centre failed at p={prime}",
        )

    expected_cases = Counter(
        {
            ("low", "low"): 5,
            ("low", "middle"): 4,
            ("low", "high"): 15,
            ("middle", "middle"): 2,
            ("middle", "high"): 10,
        }
    )
    require(cases == expected_cases, f"actual N={scale} case census changed: {cases}")

    # Actual recurrence counterexamples to the unqualified hit-position claim.
    require(
        tuple(
            sorted(
                residue_hits_top_block(scale, 1087, 77)
                + residue_hits_top_block(scale, 1087, 1009)
            )
        )
        == (1164,),
        "p=1087 boundary counterexample changed",
    )
    require(
        tuple(
            sorted(
                residue_hits_top_block(scale, 1213, 58)
                + residue_hits_top_block(scale, 1213, 1154)
            )
        )
        == (1154, 1271),
        "p=1213 mixed-lift counterexample changed",
    )
    return (
        f"top-block lift partition checked at {residue_checks} residues; "
        f"only {low_low_doublets}/{sum(cases.values())} actual doublets have centre (3p-1)/2"
    )


GaussianRational = tuple[Fraction, Fraction]


def circular_distance(first: Fraction, second: Fraction) -> Fraction:
    difference = (first - second) % 1
    return min(difference, 1 - difference)


def short_arc_integral_direct(
    frequencies: list[Fraction],
    coefficients: list[GaussianRational],
    delta: Fraction,
) -> Fraction:
    """Integrate the step function directly between all arc endpoints."""

    endpoints = {Fraction(0), Fraction(1)}
    for frequency in frequencies:
        endpoints.add(frequency % 1)
        endpoints.add((frequency - delta) % 1)
    ordered = sorted(endpoints)
    integral = Fraction(0)
    for left, right in zip(ordered, ordered[1:]):
        midpoint = (left + right) / 2
        real = Fraction(0)
        imaginary = Fraction(0)
        for frequency, coefficient in zip(frequencies, coefficients):
            forward = (frequency - midpoint) % 1
            if 0 < forward <= delta:
                real += coefficient[0]
                imaginary += coefficient[1]
        integral += (right - left) * (real * real + imaginary * imaginary)
    return integral


def short_arc_integral_pairs(
    frequencies: list[Fraction],
    coefficients: list[GaussianRational],
    delta: Fraction,
) -> Fraction:
    answer = Fraction(0)
    for first_frequency, first_coefficient in zip(frequencies, coefficients):
        for second_frequency, second_coefficient in zip(frequencies, coefficients):
            kernel = max(
                Fraction(0),
                delta - circular_distance(first_frequency, second_frequency),
            )
            real_inner_product = (
                first_coefficient[0] * second_coefficient[0]
                + first_coefficient[1] * second_coefficient[1]
            )
            answer += kernel * real_inner_product
    return answer


def check_short_arc_pair_kernel() -> str:
    frequencies = [
        Fraction(1, 17),
        Fraction(3, 19),
        Fraction(7, 23),
        Fraction(13, 29),
        Fraction(18, 19),
        Fraction(22, 23),
        Fraction(28, 29),
    ]
    coefficients: list[GaussianRational] = [
        (Fraction(2), Fraction(1)),
        (Fraction(-3, 2), Fraction(4, 3)),
        (Fraction(5, 7), Fraction(-2)),
        (Fraction(0), Fraction(3, 5)),
        (Fraction(-4), Fraction(-1, 2)),
        (Fraction(9, 4), Fraction(7, 6)),
        (Fraction(-1, 3), Fraction(5, 2)),
    ]
    checks = 0
    for denominator in (2, 3, 5, 8, 13):
        delta = Fraction(1, denominator)
        direct = short_arc_integral_direct(frequencies, coefficients, delta)
        paired = short_arc_integral_pairs(frequencies, coefficients, delta)
        require(direct == paired, f"pair-kernel mismatch at delta={delta}")

        diagonal = delta * sum(
            real * real + imaginary * imaginary
            for real, imaginary in coefficients
        )
        off_diagonal = paired - diagonal
        independently_grouped = Fraction(0)
        for first in range(len(frequencies)):
            for second in range(first + 1, len(frequencies)):
                kernel = max(
                    Fraction(0),
                    delta
                    - circular_distance(frequencies[first], frequencies[second]),
                )
                dot = (
                    coefficients[first][0] * coefficients[second][0]
                    + coefficients[first][1] * coefficients[second][1]
                )
                independently_grouped += 2 * kernel * dot
        require(
            off_diagonal == independently_grouped,
            f"diagonal/off-diagonal pair grouping failed at delta={delta}",
        )
        checks += 1
    return f"circular triangular pair kernel checked exactly for {checks} arc lengths"


def gamma_fejer(scale: int, mode: int) -> float:
    """The full Fejer weight, with its removable value at mode zero."""

    if mode == 0:
        return 1.0 / (scale * scale)
    if mode % scale == 0:
        return 0.0
    # Reduce the sine argument before floating evaluation.  This also makes
    # the exact zeroes at nonzero multiples of scale evaluate as zero.
    sine_residue = mode % (2 * scale)
    numerator = sin(pi * sine_residue / scale)
    return (numerator / (pi * mode)) ** 2


def spectral_polynomial_value(
    frequencies: list[Fraction],
    coefficients: list[GaussianRational],
    mode: int,
) -> complex:
    answer = 0j
    for frequency, (real, imaginary) in zip(frequencies, coefficients):
        phase_numerator = (-mode * frequency.numerator) % frequency.denominator
        phase = phase_numerator / frequency.denominator
        answer += complex(float(real), float(imaginary)) * cmath.exp(2j * pi * phase)
    return answer


def check_full_fejer_spectral_identity() -> str:
    """Independently enclose the infinite spectral sum by a finite one."""

    frequencies = [
        Fraction(1, 5),
        Fraction(2, 7),
        Fraction(4, 11),
        Fraction(6, 13),
        Fraction(12, 13),
    ]
    coefficients: list[GaussianRational] = [
        (Fraction(1), Fraction(0)),
        (Fraction(-1, 2), Fraction(1, 3)),
        (Fraction(2, 5), Fraction(-3, 7)),
        (Fraction(1, 4), Fraction(1, 2)),
        (Fraction(-2, 3), Fraction(-1, 5)),
    ]
    cutoff = 20_000
    coefficient_bound = sum(
        float(abs(real) + abs(imaginary)) for real, imaginary in coefficients
    )
    checks = 0
    for scale in (3, 5, 8):
        delta = Fraction(1, scale)
        exact = float(short_arc_integral_pairs(frequencies, coefficients, delta))

        summed_real = sum(real for real, _ in coefficients)
        summed_imaginary = sum(imaginary for _, imaginary in coefficients)
        zero_mode_exact = Fraction(
            summed_real * summed_real + summed_imaginary * summed_imaginary,
            scale * scale,
        )
        zero_mode_spectral = gamma_fejer(scale, 0) * abs(
            spectral_polynomial_value(frequencies, coefficients, 0)
        ) ** 2
        require(
            abs(zero_mode_spectral - float(zero_mode_exact)) <= 2e-15,
            f"Fejer k=0 term failed at M={scale}",
        )
        require(gamma_fejer(scale, scale) == 0.0, f"gamma_M(M) failed at M={scale}")
        require(
            gamma_fejer(scale, 2 * scale) == 0.0,
            f"gamma_M(2M) failed at M={scale}",
        )

        terms = []
        for mode in range(-cutoff, cutoff + 1):
            value = spectral_polynomial_value(frequencies, coefficients, mode)
            terms.append(gamma_fejer(scale, mode) * abs(value) ** 2)
        partial = fsum(terms)

        # gamma_M(k) <= 1/(pi*k)^2 and |D(k)| <= sum_j |c_j|.
        # The latter is bounded above by sum_j (|Re c_j|+|Im c_j|), so
        # the omitted two-sided tail is at most the following quantity.
        tail_bound = 2 * coefficient_bound**2 / (pi**2 * cutoff)
        difference = exact - partial
        require(
            difference >= -2e-12,
            f"finite Fejer sum exceeds the exact integral at M={scale}",
        )
        require(
            difference <= tail_bound + 2e-12,
            f"Fejer spectral tail enclosure failed at M={scale}",
        )
        checks += 1
    return (
        f"full Fejer identity (including k=0 and gamma_M) enclosed at "
        f"{checks} scales with |k|<={cutoff}"
    )


def least_absolute_residue(value: int, modulus: int) -> int:
    """The unique least-absolute residue for the odd moduli used here."""

    require(modulus % 2 == 1, "least-residue check expects an odd modulus")
    residue = value % modulus
    if residue > modulus // 2:
        residue -= modulus
    return residue


def check_cross_prime_reciprocal_bijection() -> str:
    zero_sets = actual_zero_sets()
    cases = [
        (17, 19, 32),
        (19, 17, 32),
        (37, 41, 64),
        (41, 37, 64),
        (181, 191, 256),
        (191, 181, 256),
    ]
    near_pair_checks = 0
    for prime, other_prime, scale in cases:
        modulus = prime * other_prime
        inverse_other = pow(other_prime, -1, prime)
        inverse_prime = pow(prime, -1, other_prime)
        near_pairs: dict[int, tuple[int, int]] = {}

        for first_frequency in range(1, prime):
            for second_frequency in range(1, other_prime):
                determinant = (
                    first_frequency * other_prime
                    - second_frequency * prime
                )
                least = least_absolute_residue(determinant, modulus)
                require(least != 0, f"zero determinant for p={prime}, q={other_prime}")
                require(
                    gcd(least, modulus) == 1,
                    f"non-coprime least residue for p={prime}, q={other_prime}",
                )
                if scale * abs(least) >= modulus:
                    continue
                require(
                    least not in near_pairs,
                    f"reciprocal parametrization is not injective at k={least}",
                )
                near_pairs[least] = (first_frequency, second_frequency)

                recovered_first = (least * inverse_other) % prime
                recovered_second = (-least * inverse_prime) % other_prime
                require(
                    (recovered_first, recovered_second)
                    == (first_frequency, second_frequency),
                    f"reciprocal recovery failed at k={least}",
                )

                distance = circular_distance(
                    Fraction(first_frequency, prime),
                    Fraction(second_frequency, other_prime),
                )
                require(
                    distance == Fraction(abs(least), modulus),
                    f"least-residue distance failed at k={least}",
                )

        largest = (modulus - 1) // scale
        expected = {
            least
            for least in range(-largest, largest + 1)
            if least != 0
            and scale * abs(least) < modulus
            and gcd(least, modulus) == 1
        }
        require(
            set(near_pairs) == expected,
            f"reciprocal parametrization is not surjective for p={prime}, q={other_prime}",
        )

        first_zeros = zero_sets[prime]
        second_zeros = zero_sets[other_prime]
        require(
            first_zeros and second_zeros,
            f"Fourier-product guard would be vacuous for p={prime}, q={other_prime}",
        )
        for least, (first_frequency, second_frequency) in near_pairs.items():
            reciprocal_first = (least * inverse_other) % prime
            reciprocal_second = (least * inverse_prime) % other_prime
            require(reciprocal_first == first_frequency, "first reciprocal frequency changed")
            require(
                reciprocal_second == (-second_frequency) % other_prime,
                "second reciprocal frequency has the wrong sign",
            )

            # Compare the products as exact multisets of pq-th roots of unity.
            original_exponents = Counter(
                (
                    first_frequency * first_residue * other_prime
                    - second_frequency * second_residue * prime
                )
                % modulus
                for first_residue in first_zeros
                for second_residue in second_zeros
            )
            reciprocal_exponents = Counter(
                (
                    reciprocal_first * first_residue * other_prime
                    + reciprocal_second * second_residue * prime
                )
                % modulus
                for first_residue in first_zeros
                for second_residue in second_zeros
            )
            require(
                original_exponents == reciprocal_exponents,
                f"exact Fourier-product identity failed at k={least}",
            )

            original_value = fourier_value(
                prime, first_zeros, first_frequency
            ) * fourier_value(
                other_prime, second_zeros, second_frequency
            ).conjugate()
            reciprocal_value = fourier_value(
                prime, first_zeros, reciprocal_first
            ) * fourier_value(
                other_prime, second_zeros, reciprocal_second
            )
            require(
                abs(original_value - reciprocal_value) <= 2e-10,
                f"numerical Fourier-product guard failed at k={least}",
            )
            near_pair_checks += 1
    return (
        f"least-residue (a,b)<->k bijection and Fourier product checked at "
        f"{near_pair_checks} ordered cross-prime pairs"
    )


def nonzero_fourier_inversion_exact(
    prime: int, zeros: tuple[int, ...], integer: int
) -> Fraction:
    """Collapse the a-sum using exact root-of-unity orthogonality."""

    numerator = sum(
        prime - 1 if (residue - integer) % prime == 0 else -1
        for residue in zeros
    )
    return Fraction(numerator, prime)


def check_fourier_inversion_and_interval_twist() -> str:
    zero_sets = actual_zero_sets()
    inversion_checks = 0
    for prime in (11, 17, 181, 379, 1087, 1213):
        zeros = zero_sets[prime]
        for integer in range(-2, 2 * prime + 3, max(1, prime // 17)):
            expected = Fraction(int(integer % prime in zeros), 1) - Fraction(
                len(zeros), prime
            )
            exact = nonzero_fourier_inversion_exact(prime, zeros, integer)
            require(exact == expected, f"Fourier inversion failed at p={prime}, m={integer}")

            numerical = sum(
                fourier_value(prime, zeros, frequency)
                * ep(prime, -frequency * integer)
                for frequency in range(1, prime)
            ) / prime
            require(
                abs(numerical - float(expected)) <= 2e-10,
                f"numerical Fourier sign guard failed at p={prime}, m={integer}",
            )
            inversion_checks += 1

    # A plus sign in the inversion exponent recovers the wrong reflected
    # residue.  The actual p=17 zero at m=3 makes this visible.
    prime, integer = 17, 3
    zeros = zero_sets[prime]
    wrong_sign = sum(
        fourier_value(prime, zeros, frequency)
        * ep(prime, frequency * integer)
        for frequency in range(1, prime)
    ) / prime
    expected = 1 - len(zeros) / prime
    require(abs(wrong_sign - expected) > 0.5, "Fourier inversion sign guard disappeared")

    # Translating m=N+j to j in (0,N] twists every coefficient by
    # e_p(-a*N).  This is the exact input required by a shifted-interval
    # Gallagher application; it is not the untwisted coefficient family.
    shift_checks = 0
    for prime in (17, 181, 1087, 1213):
        for frequency in (1, 2, (prime - 1) // 2, prime - 1):
            for offset in (1, 2, TOP_N // 2, TOP_N):
                direct_exponent = (-frequency * (TOP_N + offset)) % prime
                translated_exponent = (
                    -frequency * TOP_N - frequency * offset
                ) % prime
                require(
                    direct_exponent == translated_exponent,
                    f"interval-twist exponent failed at p={prime}, a={frequency}, j={offset}",
                )
                shift_checks += 1
    return (
        f"nonzero Fourier inversion checked {inversion_checks} times with exponent -am; "
        f"the required I_N translation twist was checked {shift_checks} times"
    )


def reflected_anchored_primes(scale: int) -> list[int]:
    excluded_product = (4 * scale + 1) * (2 * scale + 1)
    return [
        prime
        for prime in primes_up_to(2 * scale)
        if scale < prime <= 2 * scale and excluded_product % prime != 0
    ]


def check_reflected_anchored_family() -> str:
    scales = (8, 16, 32, 64)
    column_checks = 0
    for scale in scales:
        primes = reflected_anchored_primes(scale)
        require(primes, f"anchored family is empty at N={scale}")
        d_at_anchor_fourier = 0j
        d_at_anchor_incidence = 0.0
        d_at_zero_fourier = 0j
        d_at_zero_incidence = 0.0
        q_fourier = 0.0
        q_parseval = 0.0

        for prime in primes:
            first = (2 * scale) % prime
            second = (-2 * scale - 1) % prime
            zeros = (first, second)
            require(first != second, f"anchored residues coincide at N={scale}, p={prime}")
            require(0 not in zeros, f"anchored family contains zero at N={scale}, p={prime}")
            require(
                (first - second) % prime not in (1, prime - 1),
                f"anchored residues are adjacent at N={scale}, p={prime}",
            )
            require(
                {prime - 1 - residue for residue in zeros} == set(zeros),
                f"anchored reflection failed at N={scale}, p={prime}",
            )
            require(
                (4 * scale + 1) % prime != 0,
                f"distinctness exclusion failed at N={scale}, p={prime}",
            )
            require(
                (2 * scale + 1) % prime != 0,
                f"zero exclusion failed at N={scale}, p={prime}",
            )

            hits = set(residue_hits_top_block(scale, prime, first)) | set(
                residue_hits_top_block(scale, prime, second)
            )
            require(2 * scale in hits, f"common anchor is missing at N={scale}, p={prime}")
            require(len(hits) <= 2, f"anchored column has too many hits at N={scale}, p={prime}")

            anchor_exact = nonzero_fourier_inversion_exact(
                prime, zeros, 2 * scale
            )
            zero_exact = nonzero_fourier_inversion_exact(prime, zeros, 0)
            require(
                anchor_exact == 1 - Fraction(2, prime),
                f"exact anchored D(2N) summand failed at N={scale}, p={prime}",
            )
            require(
                zero_exact == -Fraction(2, prime),
                f"exact anchored D(0) summand failed at N={scale}, p={prime}",
            )

            weight = log(prime)
            anchor_fourier_summand = sum(
                fourier_value(prime, zeros, frequency)
                * ep(prime, -frequency * 2 * scale)
                for frequency in range(1, prime)
            ) / prime
            zero_fourier_summand = sum(
                fourier_value(prime, zeros, frequency)
                for frequency in range(1, prime)
            ) / prime
            d_at_anchor_fourier += weight * anchor_fourier_summand
            d_at_anchor_incidence += weight * (1 - 2 / prime)
            d_at_zero_fourier += weight * zero_fourier_summand
            d_at_zero_incidence -= weight * 2 / prime

            # Exact finite Fourier Parseval for this two-point set:
            # sum_{a=1}^{p-1}|F_p(a)|^2 = 2p-4.
            exact_energy = sum(
                prime - 1 if (left - right) % prime == 0 else -1
                for left in zeros
                for right in zeros
            )
            require(
                exact_energy == 2 * prime - 4,
                f"anchored exact Q Parseval failed at N={scale}, p={prime}",
            )
            numerical_energy = fsum(
                abs(fourier_value(prime, zeros, frequency)) ** 2
                for frequency in range(1, prime)
            )
            require(
                abs(numerical_energy - exact_energy) <= 2e-9,
                f"anchored numerical Q guard failed at N={scale}, p={prime}",
            )
            q_fourier += weight**2 * numerical_energy / prime**2
            q_parseval += weight**2 * (2 / prime - 4 / prime**2)
            column_checks += 1

        require(
            abs(d_at_anchor_fourier - d_at_anchor_incidence) <= 2e-10,
            f"anchored aggregate D(2N) failed at N={scale}",
        )
        require(
            abs(d_at_zero_fourier - d_at_zero_incidence) <= 2e-10,
            f"anchored aggregate D(0) failed at N={scale}",
        )
        require(
            abs(q_fourier - q_parseval) <= 2e-10,
            f"anchored aggregate Q Parseval failed at N={scale}",
        )

        gamma = gamma_fejer(4 * scale, 2 * scale)
        expected_gamma = 1 / (4 * pi**2 * scale**2)
        require(
            abs(gamma - expected_gamma) <= 2e-18,
            f"gamma_(4N)(2N) identity failed at N={scale}",
        )
        require(
            abs(gamma * (4 * pi**2 * scale**2) - 1) <= 3e-15,
            f"normalized gamma_(4N)(2N) identity failed at N={scale}",
        )
    return (
        f"reflected anchored family checked at {column_checks} columns over "
        f"N={','.join(map(str, scales))}, including nonadjacency, D(2N), "
        f"D(0), Q, and gamma_(4N)(2N)"
    )


def main() -> int:
    checks = [
        ("recurrence/bank", check_recurrence_and_bank),
        ("palindromic phase", check_palindromic_phase),
        ("doublet Fourier", check_doublet_fourier_formula),
        ("top-block hits", check_top_block_hit_positions),
        ("short-arc kernel", check_short_arc_pair_kernel),
        ("full Fejer spectrum", check_full_fejer_spectral_identity),
        ("cross-prime reciprocal", check_cross_prime_reciprocal_bijection),
        ("Fourier/Gallagher transcription", check_fourier_inversion_and_interval_twist),
        ("reflected anchored family", check_reflected_anchored_family),
    ]
    failures = 0
    for name, check in checks:
        try:
            result = check()
        except Exception as error:  # Deliberately turn every failure into exit 1.
            failures += 1
            print(f"FAIL: {name}: {error}", file=sys.stderr)
        else:
            print(f"PASS: {result}")
    if failures:
        print(f"FAILED {failures}/{len(checks)} CHECKS", file=sys.stderr)
        return 1
    print(f"ALL {len(checks)} CHECKS PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
