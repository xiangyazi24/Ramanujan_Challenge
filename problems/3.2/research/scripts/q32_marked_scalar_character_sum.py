#!/usr/bin/env python3
"""The marked scalar is an exponential sum over (F_p^*)^3.

Claim (exact, checked below):  for every prime p >= 5 and every moment M,

    C_M(p-1)  ==  - sum_{x,y,z in F_p^*} Lambda(x,y,z)^M   (mod p),

where Lambda(x,y,z) = (1+x)(1+y)(1+z)((1+y)(1+z)+xyz)/(xyz) is the Apery Laurent
polynomial and C_M(d) is its d-section shell.  Equivalently, with

    N_p(t) = #{ (x,y,z) in (F_p^*)^3 : Lambda(x,y,z) = t },

    C_M(p-1) == - sum_{t in F_p^*} t^{M mod (p-1)} N_p(t)   (mod p),

so the marked scalar is the r-th MOMENT of the point-count function of the Apery
family, r = M mod (p-1).  Since n = p + r gives M = n-1 == r (mod p-1), the target
condition p | F_0 is exactly the vanishing of that moment, in agreement with the
Apery-Lucas fact p | b_n <=> p | b_r.

Why it is worth recording: it rewrites the one marked mod-p scalar of the whole
programme in the standard language of point counts of the (modular, weight-4 level-8)
Apery family, where Hasse-Witt / unit-root technology is the natural tool -- and it
makes the rank-one phenomenon expected, since modulo p only the one-dimensional
unit-root datum of the family survives.
"""

from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from q32_cartier_packet_audit import apery, shell_fast  # noqa: E402
from q32_terminal_family_audit import primes_up_to  # noqa: E402


def lambda_value(x, y, z, prime):
    numerator = (1 + x) * (1 + y) % prime * (1 + z) % prime
    bracket = ((1 + y) * (1 + z) + x * y * z) % prime
    inverse = pow(x * y * z % prime, prime - 2, prime)
    return numerator * bracket % prime * inverse % prime


def exponential_sum(prime, moment):
    total = 0
    for x in range(1, prime):
        for y in range(1, prime):
            for z in range(1, prime):
                total += pow(lambda_value(x, y, z, prime), moment, prime)
        total %= prime
    return total % prime


def point_counts(prime):
    counts = [0] * prime
    for x in range(1, prime):
        for y in range(1, prime):
            for z in range(1, prime):
                counts[lambda_value(x, y, z, prime)] += 1
    return counts


def moment(prime, counts, residue):
    return sum(
        pow(t, residue, prime) * counts[t] for t in range(1, prime)
    ) % prime


def audit(prime_bound=23):
    shell_checks = 0
    moment_checks = 0
    apery_checks = 0
    for prime in primes_up_to(prime_bound):
        if prime < 5:
            continue
        counts = point_counts(prime)
        for moment_index in range(prime, 3 * prime, 2):
            left = shell_fast(moment_index, prime - 1, modulus=prime) % prime
            right = (-exponential_sum(prime, moment_index)) % prime
            assert left == right
            shell_checks += 1

            residue = moment_index % (prime - 1)
            assert left == (-moment(prime, counts, residue)) % prime
            moment_checks += 1

            if 0 <= moment_index + 1 - prime < prime - 1:
                assert left == apery(moment_index + 1 - prime) % prime
                apery_checks += 1
    return shell_checks, moment_checks, apery_checks


if __name__ == "__main__":
    shell_checks, moment_checks, apery_checks = audit()
    print("SHELL_VS_EXPONENTIAL_SUM", shell_checks)
    print("SHELL_VS_POINT_COUNT_MOMENT", moment_checks)
    print("SHELL_VS_APERY", apery_checks)
    print("Q32_MARKED_SCALAR_CHARACTER_SUM=PASS")
