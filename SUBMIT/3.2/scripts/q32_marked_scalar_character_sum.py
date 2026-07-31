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


def apery_mod_sequence(prime):
    """b_0..b_{prime-1} mod prime by the Apery recurrence."""
    inverses = [0, 1] + [0] * (prime - 2)
    for k in range(2, prime):
        inverses[k] = (prime - (prime // k) * inverses[prime % k] % prime) % prime
    values = [1 % prime, 5 % prime]
    for index in range(1, prime - 1):
        inverse = inverses[index + 1]
        lead = inverse * inverse % prime * inverse % prime
        values.append(
            ((34 * index ** 3 + 51 * index ** 2 + 27 * index + 5) % prime
             * values[index]
             - index ** 3 % prime * values[index - 1]) % prime * lead % prime
        )
    return values


def audit_reflection(prime_bound=60):
    """The reflection law, derived rather than observed.

    N_p(t) = N_p(1/t) for every t, hence sum_t t^r N_p(t) = sum_t t^{-r} N_p(t),
    hence b_{p-1-r} = b_r (mod p), hence Z_p is stable under r -> p-1-r and |Z_p|
    is even unless the fixed point (p-1)/2 lies in Z_p.
    """
    inversion_checks = 0
    palindrome_checks = 0
    involution_checks = 0
    odd_primes = []
    for prime in primes_up_to(prime_bound):
        if prime < 5:
            continue
        counts = point_counts(prime)
        for t in range(1, prime):
            assert counts[t] == counts[pow(t, prime - 2, prime)]
            inversion_checks += 1
        values = apery_mod_sequence(prime)
        for residue in range(prime - 1):
            assert values[prime - 1 - residue] == values[residue]
            palindrome_checks += 1
        zeros = {r for r in range(prime - 1) if values[r] == 0}
        for r in zeros:
            assert prime - 1 - r in zeros
            involution_checks += 1
        if len(zeros) % 2:
            assert (prime - 1) % 2 == 0 and (prime - 1) // 2 in zeros
            odd_primes.append(prime)
    return inversion_checks, palindrome_checks, involution_checks, odd_primes


if __name__ == "__main__":
    shell_checks, moment_checks, apery_checks = audit()
    print("SHELL_VS_EXPONENTIAL_SUM", shell_checks)
    print("SHELL_VS_POINT_COUNT_MOMENT", moment_checks)
    print("SHELL_VS_APERY", apery_checks)
    inv_c, pal_c, invol_c, odd_p = audit_reflection()
    print("POINT_COUNT_INVERSION_SYMMETRY", inv_c)
    print("APERY_PALINDROME_MOD_P", pal_c)
    print("ZERO_SET_INVOLUTION", invol_c)
    print("PRIMES_WITH_ODD_ZERO_COUNT", odd_p)
    print("Q32_MARKED_SCALAR_CHARACTER_SUM=PASS")
