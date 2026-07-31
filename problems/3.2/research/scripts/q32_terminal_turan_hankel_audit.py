#!/usr/bin/env python3
"""Exact audit of terminal Turan and Hankel target carriers.

For n=M+1, let

    f_L = G_{M-L,L}(C_M),  F_j = f_{L0-j},
    L0 = M - (floor(M/2)+1).

The script checks the normalized adjacent Turan minors

    E_j = (F_{j-1} F_{j+1} - F_j^2)
          / gcd(binomial(n,L0-j), binomial(n,L0-j+1))

and the first two Hankel determinants.  These use the scalar shell
Y_d=C_M(d), not the two-coordinate Y/W shell.
"""

from math import comb, gcd
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from q32_cartier_packet_audit import shell_batch  # noqa: E402
from q32_terminal_family_audit import newton, primes_up_to  # noqa: E402


def valuation(value, prime):
    exponent = 0
    while value and value % prime == 0:
        value //= prime
        exponent += 1
    return exponent


def factor_over_primes(value, bound):
    value = abs(value)
    answer = []
    for prime in primes_up_to(bound):
        exponent = valuation(value, prime)
        if exponent:
            answer.append((prime, exponent))
            value //= prime**exponent
    if value > 1:
        answer.append((value, 1))
    return tuple(answer)


def terminal_data(n):
    moment = n - 1
    start = moment // 2 + 1
    maximum = moment - start
    values = shell_batch(moment, range(start, moment + 1))
    f_values = [
        newton(values, moment - order, order)
        for order in range(maximum + 1)
    ]
    terminal = tuple(reversed(f_values))
    pascal = tuple(comb(n, maximum - shift) for shift in range(maximum + 1))
    return start, maximum, terminal, pascal


def turan_carrier(terminal, pascal, shift):
    """Return E_shift for 1 <= shift < L0."""

    numerator = (
        terminal[shift - 1] * terminal[shift + 1]
        - terminal[shift] ** 2
    )
    denominator = gcd(pascal[shift], pascal[shift - 1])
    assert numerator % denominator == 0
    return numerator // denominator


def hankel_two(terminal):
    return terminal[0] * terminal[2] - terminal[1] ** 2


def hankel_three(terminal):
    f0, f1, f2, f3, f4 = terminal[:5]
    return (
        f0 * f2 * f4
        + 2 * f1 * f2 * f3
        - f0 * f3**2
        - f4 * f1**2
        - f2**3
    )


def pascal_gcd(pascal, width):
    answer = 0
    for value in pascal[:width]:
        answer = gcd(answer, value)
    return answer


def target_primes(n, start, maximum, terminal, width):
    """Targets in the interval where F_0,...,F_width all select the node."""

    return tuple(
        prime
        for prime in primes_up_to(n)
        if start + width < prime <= n and terminal[0] % prime == 0
    )


def candidate_primes(n, start, width):
    """Primes for which F_0,...,F_width share the selected node p-1."""

    return tuple(
        prime
        for prime in primes_up_to(n)
        if start + width < prime <= n
    )


def assert_rank_one_alias(terminal, pascal, start, n, max_shift=7):
    """Check the once-divided rank-one formula on every safe candidate.

    This deliberately includes non-target candidates.  Targetness is the
    extra condition ``terminal[shift] == 0 (mod p)``; it is not used in
    the derivation of the alias formula.
    """

    checks = 0
    for shift in range(1, min(max_shift + 1, len(terminal) - 1)):
        for prime in candidate_primes(n, start, shift + 1):
            left_increment = terminal[shift - 1] - terminal[shift]
            right_increment = terminal[shift] - terminal[shift + 1]
            denominator = gcd(pascal[shift - 1], pascal[shift])
            assert left_increment % prime == 0
            assert right_increment % prime == 0
            assert valuation(denominator, prime) == 1
            unit = (denominator // prime) % prime
            assert unit
            alpha_difference = (
                left_increment // prime - right_increment // prime
            ) % prime
            carrier = turan_carrier(terminal, pascal, shift)
            assert (
                unit * carrier - terminal[shift] * alpha_difference
            ) % prime == 0
            checks += 1
    return checks


def audit_record(n):
    start, maximum, terminal, pascal = terminal_data(n)
    assert maximum >= 4

    e1 = turan_carrier(terminal, pascal, 1)
    e2 = turan_carrier(terminal, pascal, 2)
    turan_gcd = gcd(e1, e2)

    a = terminal[0] - terminal[1]
    b = terminal[1] - terminal[2]
    c = terminal[2] - terminal[3]
    q0 = gcd(pascal[0], pascal[1])
    q1 = gcd(pascal[1], pascal[2])
    assert (b - c) % q1 == 0
    assert (a - b) % q0 == 0
    assert b * (a * c - b * b) % (q0 * q1) == 0
    boundary_residual = b * (a * c - b * b) // (q0 * q1)
    assert (
        (b - c) // q1 * e1
        - (a - b) // q0 * e2
        == boundary_residual
    )

    d2 = hankel_two(terminal)
    q2 = pascal_gcd(pascal, 2)
    assert d2 % q2 == 0
    t2 = d2 // q2
    assert t2 == e1

    d3 = hankel_three(terminal)
    q3 = pascal_gcd(pascal, 3)
    q4 = pascal_gcd(pascal, 4)
    assert d3 % (q3 * q4) == 0
    t3 = d3 // (q3 * q4)

    # Desnanot--Jacobi in the indexing used by the report.
    a0 = terminal[0] * terminal[2] - terminal[1] ** 2
    a1 = terminal[1] * terminal[3] - terminal[2] ** 2
    a2 = terminal[2] * terminal[4] - terminal[3] ** 2
    assert terminal[2] * d3 == a0 * a2 - a1**2

    targets = target_primes(n, start, maximum, terminal, 3)
    for prime in targets:
        assert turan_gcd % prime == 0
        assert t2 % prime == 0
        assert t3 % prime == 0
        # The hostile examples have exactly the expected one surviving digit.
        assert valuation(t2, prime) >= 1
        assert valuation(t3, prime) >= 1

    # General adjacent elimination identity.
    for shift in range(1, min(8, maximum - 2)):
        a0 = terminal[shift - 1] - terminal[shift]
        a1 = terminal[shift] - terminal[shift + 1]
        a2 = terminal[shift + 1] - terminal[shift + 2]
        raw0 = (
            terminal[shift - 1] * terminal[shift + 1]
            - terminal[shift] ** 2
        )
        raw1 = (
            terminal[shift] * terminal[shift + 2]
            - terminal[shift + 1] ** 2
        )
        assert (
            (a1 - a2) * raw0 - (a0 - a1) * raw1
            == a1 * (a0 * a2 - a1 * a1)
        )

    rank_one_checks = assert_rank_one_alias(
        terminal, pascal, start, n
    )

    return {
        "n": n,
        "start": start,
        "maximum": maximum,
        "targets": targets,
        "e1_bits": abs(e1).bit_length(),
        "e2_bits": abs(e2).bit_length(),
        "turan_gcd": turan_gcd,
        "turan_gcd_factors": factor_over_primes(turan_gcd, n),
        "boundary_residual_bits": abs(boundary_residual).bit_length(),
        "t3_bits": abs(t3).bit_length(),
        "hankel_gcd": gcd(t2, t3),
        "hankel_gcd_factors": factor_over_primes(gcd(t2, t3), n),
        "rank_one_checks": rank_one_checks,
    }


def audit_rank_one_dense():
    """Exercise both target and non-target candidates on consecutive n."""

    checks = 0
    targets = 0
    nontargets = 0
    for n in range(12, 61):
        start, maximum, terminal, pascal = terminal_data(n)
        if maximum < 2:
            continue
        checks += assert_rank_one_alias(
            terminal, pascal, start, n, max_shift=5
        )
        for shift in range(1, min(6, len(terminal) - 1)):
            for prime in candidate_primes(n, start, shift + 1):
                if terminal[shift] % prime == 0:
                    targets += 1
                else:
                    nontargets += 1
    assert checks == targets + nontargets
    assert targets
    assert nontargets
    return checks, targets, nontargets


if __name__ == "__main__":
    for test_n in (80, 120, 160, 200, 272, 300, 321, 340, 380, 500):
        record = audit_record(test_n)
        print(
            "N",
            record["n"],
            "TARGET_INTERVAL",
            f"({record['start'] + 3},{record['n']}]",
            "TARGETS",
            record["targets"],
            "E_BITS",
            (record["e1_bits"], record["e2_bits"]),
            "TURAN_GCD",
            record["turan_gcd_factors"],
            "BOUNDARY_RESIDUAL_BITS",
            record["boundary_residual_bits"],
            "T3_BITS",
            record["t3_bits"],
            "HANKEL_GCD",
            record["hankel_gcd_factors"],
            "RANK_ONE_CHECKS",
            record["rank_one_checks"],
        )
    dense_checks, dense_targets, dense_nontargets = audit_rank_one_dense()
    print(
        "DENSE_RANK_ONE",
        dense_checks,
        "TARGET_CASES",
        dense_targets,
        "NONTARGET_CASES",
        dense_nontargets,
    )
    print("Q32_TERMINAL_TURAN_HANKEL_AUDIT=PASS")
