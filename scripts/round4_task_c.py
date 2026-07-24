#!/usr/bin/env python3
"""Round 4, Task C: independently verify the record Z(p) at p=159977.

For the zeta(3) Apery numbers

    b_j = sum_{k=0}^j binom(j,k)^2 binom(j+k,k)^2,

the script computes every residue b_0,...,b_{p-1} with the standard
three-term recurrence, lists all zero positions, and verifies their complete
reflection pairing.  It also performs two non-tautological reproducibility
checks:

* the entire residue table is checked for b_j = b_{p-1-j};
* every detected zero and deterministic control indices are recomputed from
  the binomial-sum definition using factorial tables, independently of the
  recurrence.

The default report path is /tmp/round4_task_c.txt.  Only the Python standard
library is required.
"""

from __future__ import annotations

import argparse
import hashlib
import math
import random
from pathlib import Path
from typing import Sequence


DEFAULT_PRIME = 159_977
DEFAULT_OUTPUT = Path("/tmp/round4_task_c.txt")
EXPECTED_ZERO_COUNT = 12
SAMPLE_SEED = 20260715


def is_prime(n: int) -> bool:
    """Deterministic trial-division primality test (ample for the target p)."""
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    for divisor in range(3, math.isqrt(n) + 1, 2):
        if n % divisor == 0:
            return False
    return True


def inverse_table(p: int) -> list[int]:
    """Return inverses of 1,...,p-1 in F_p in O(p) time."""
    inverses = [0] * p
    inverses[1] = 1
    for value in range(2, p):
        inverses[value] = (p - (p // value) * inverses[p % value] % p) % p
    return inverses


def apery_table_from_recurrence(p: int) -> list[int]:
    """Compute b_j mod p for every 0 <= j < p.

    Initial values and recurrence:

        b_0 = 1,  b_1 = 5,
        (n+1)^3 b_{n+1}
          = (34n^3+51n^2+27n+5)b_n - n^3 b_{n-1}.

    Only n=1,...,p-2 is used, so division by (n+1)^3 is valid in F_p.
    """
    if p < 5 or not is_prime(p):
        raise ValueError(f"p must be a prime >= 5, got {p}")

    inverses = inverse_table(p)
    residues = [0] * p
    residues[0] = 1
    residues[1] = 5 % p

    for n in range(1, p - 1):
        n2 = n * n % p
        n3 = n2 * n % p
        coefficient = (((34 * n + 51) * n + 27) * n + 5) % p
        inv = inverses[n + 1]
        inverse_cube = inv * inv % p * inv % p
        residues[n + 1] = (
            (coefficient * residues[n] - n3 * residues[n - 1])
            * inverse_cube
            % p
        )

    return residues


def recurrence_failures(residues: Sequence[int], p: int) -> list[int]:
    """Return indices where the recurrence identity fails after multiplication."""
    failures: list[int] = []
    for n in range(1, p - 1):
        coefficient = (34 * n**3 + 51 * n**2 + 27 * n + 5) % p
        lhs = (n + 1) ** 3 % p * residues[n + 1] % p
        rhs = (coefficient * residues[n] - n**3 * residues[n - 1]) % p
        if lhs != rhs:
            failures.append(n)
    return failures


def factorial_tables(p: int) -> tuple[list[int], list[int]]:
    """Build factorial and inverse-factorial tables in F_p."""
    factorial = [1] * p
    for n in range(1, p):
        factorial[n] = factorial[n - 1] * n % p

    inverse_factorial = [1] * p
    inverse_factorial[p - 1] = pow(factorial[p - 1], -1, p)
    for n in range(p - 1, 0, -1):
        inverse_factorial[n - 1] = inverse_factorial[n] * n % p
    return factorial, inverse_factorial


def apery_binomial_mod(
    j: int,
    p: int,
    factorial: Sequence[int],
    inverse_factorial: Sequence[int],
) -> int:
    """Compute b_j mod p directly from its binomial sum.

    When j+k >= p, Lucas' theorem makes binom(j+k,k) zero because 0 < j < p
    (the endpoint cases agree as well).  Thus it is sufficient to sum through
    min(j,p-j-1), keeping every factorial argument below p.
    """
    if not 0 <= j < p:
        raise ValueError(f"j must lie in [0,p), got j={j}, p={p}")

    last_k = min(j, p - j - 1)
    total = 0
    inv_factorial_j = inverse_factorial[j]
    factorial_j = factorial[j]
    for k in range(last_k + 1):
        choose_j_k = (
            factorial_j
            * inverse_factorial[k]
            % p
            * inverse_factorial[j - k]
            % p
        )
        choose_j_plus_k_k = (
            factorial[j + k]
            * inverse_factorial[k]
            % p
            * inv_factorial_j
            % p
        )
        product = choose_j_k * choose_j_plus_k_k % p
        total = (total + product * product) % p
    return total


def direct_check_indices(p: int, zeros: Sequence[int]) -> tuple[int, ...]:
    """Choose deterministic definition-level checks, including every zero."""
    indices = set(zeros)
    indices.update(
        {
            0,
            1,
            2,
            3,
            5,
            8,
            11,
            31,
            73,
            97,
            (p - 1) // 4,
            (p - 1) // 3,
            (p - 1) // 2,
            2 * (p - 1) // 3,
            3 * (p - 1) // 4,
            p - 4,
            p - 3,
            p - 2,
            p - 1,
        }
    )
    rng = random.Random(SAMPLE_SEED)
    indices.update(rng.randrange(p) for _ in range(12))
    return tuple(sorted(j for j in indices if 0 <= j < p))


def residue_digest(residues: Sequence[int]) -> str:
    """SHA-256 of the residues encoded as unsigned 4-byte big-endian words."""
    digest = hashlib.sha256()
    for residue in residues:
        digest.update(residue.to_bytes(4, byteorder="big", signed=False))
    return digest.hexdigest()


def build_report(p: int, residues: Sequence[int]) -> str:
    zeros = tuple(j for j, value in enumerate(residues) if value == 0)
    zero_set = set(zeros)
    center = (p - 1) // 2
    pairs = tuple((j, p - 1 - j) for j in zeros if j < center)

    recurrence_bad = recurrence_failures(residues, p)
    palindrome_bad = [
        j for j in range(p) if residues[j] != residues[p - 1 - j]
    ]
    unpaired = [j for j in zeros if p - 1 - j not in zero_set]

    factorial, inverse_factorial = factorial_tables(p)
    check_indices = direct_check_indices(p, zeros)
    direct_values = {
        j: apery_binomial_mod(j, p, factorial, inverse_factorial)
        for j in check_indices
    }
    direct_bad = [j for j in check_indices if direct_values[j] != residues[j]]

    # These assertions make the script a regression test, not just a report.
    assert len(residues) == p
    assert residues[0] == 1 and residues[1] == 5 % p
    assert not recurrence_bad, f"recurrence failures: {recurrence_bad[:10]}"
    assert not palindrome_bad, f"palindrome failures: {palindrome_bad[:10]}"
    assert not unpaired, f"unpaired zeros: {unpaired}"
    assert len(zeros) == EXPECTED_ZERO_COUNT, (
        f"expected Z({p})={EXPECTED_ZERO_COUNT}, got {len(zeros)}"
    )
    assert len(pairs) * 2 == len(zeros)
    assert all(left + right == p - 1 for left, right in pairs)
    assert not direct_bad, f"binomial-sum mismatches: {direct_bad}"
    assert all(direct_values[j] == 0 for j in zeros)

    direct_rows = ", ".join(
        f"{j}:{direct_values[j]}" for j in check_indices
    )
    lines = [
        "Round 4 - Task C: record Z(p) verification",
        "=" * 52,
        f"p = {p} (primality check: PASS)",
        f"computed range = 0 <= j <= {p - 1} ({len(residues)} residues)",
        "initial values = b_0=1, b_1=5",
        "recurrence = (n+1)^3 b_(n+1) - "
        "(34n^3+51n^2+27n+5)b_n + n^3 b_(n-1) = 0 (mod p)",
        f"full recurrence recheck = PASS ({p - 2} identities)",
        "",
        f"Z({p}) = {len(zeros)}",
        f"zero positions = {list(zeros)}",
        f"center j={(p - 1) // 2} is zero = {center in zero_set}",
        f"palindromic pairs = {list(pairs)}",
        f"pair sums = {[left + right for left, right in pairs]}",
        f"zero pairing j <-> p-1-j = PASS ({len(pairs)} pairs)",
        f"full-table palindrome b_j=b_(p-1-j) = PASS ({p} comparisons)",
        "",
        "Independent binomial-sum cross-check",
        "--------------------------------------",
        "definition = b_j=sum_k binom(j,k)^2 binom(j+k,k)^2 (mod p)",
        f"checked indices ({len(check_indices)}) = {list(check_indices)}",
        f"direct values j:b_j = {direct_rows}",
        "all recurrence/direct comparisons = PASS",
        "all 12 recurrence-detected zeros independently evaluate to 0 = PASS",
        "",
        "Reproducibility fingerprint",
        "---------------------------",
        "encoding = concatenated 4-byte unsigned big-endian residues b_0,...,b_(p-1)",
        f"SHA-256 = {residue_digest(residues)}",
        "",
        f"FINAL: verified exactly {len(zeros)} zeros at p={p}, in {len(pairs)} "
        "palindromic pairs.",
    ]
    return "\n".join(lines) + "\n"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prime", type=int, default=DEFAULT_PRIME)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    residues = apery_table_from_recurrence(args.prime)
    report = build_report(args.prime, residues)
    args.output.write_text(report, encoding="utf-8")
    print(report, end="")
    print(f"Report written to {args.output}")


if __name__ == "__main__":
    main()
