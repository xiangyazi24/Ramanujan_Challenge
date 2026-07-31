#!/usr/bin/env python3
"""Two exact facts about the bad-prime criterion, in the shift variable.

(1) THE 27-TERM CRITERION.  Since p = 1 (mod p-1), the exponent is p-independent:
    for p in (n/2, n] and m = n-1,

        p | b_n  <=>  sum_{x,y,z in F_p^*} Lambda(x,y,z)^m = 0 (mod p)
                 <=>  p | V(n, p-1),
        V(n,s) := sum_{eps in {-1,0,1}^3} c_m(eps_1 s, eps_2 s, eps_3 s),

    the 27 terms being the only lattice points of mP in the sublattice sZ^3 once s > m/2.

(2) THE SHIFT OPERATOR.  C(s) = c_m(s,0,0) satisfies, exactly,

        q2(s) C(s+2) + q1(s) C(s+1) + q0(s) C(s) = 0,
        q0(s) = (s-m)^3,
        q1(s) = 2s^3 - 4(m-1)s^2 + (m^2-7m+3)s + (2m^3+4m^2-2m+1),
        q2(s) = s^3 - (m-4)s^2 - (m+5)(m-1)s + (m-1)^2(m+2).

    Order 2, degree 3 in s, with coefficients polynomial in m as well.  Note q0 vanishes to
    order 3 at s = m, the edge of the support (C(s) = 0 for s > m).

Why this matters: the criterion needs C at s = p-1, i.e. at s = -1 mod p.  A certificate
C_n divisible by every bad prime with log|C_n| = o(n) would close the conjecture, and the
transfer-matrix product of this operator over a full period mod p is exactly the p-curvature
object of a (globally nilpotent) Picard-Fuchs-type operator.

Both facts are checked below.
"""

import itertools
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from q32_cartier_packet_audit import apery, coefficient  # noqa: E402


def shell_sum(moment, dilation, modulus=None):
    total = 0
    for signs in itertools.product((-1, 0, 1), repeat=3):
        total += coefficient(
            moment, *(s * dilation for s in signs), modulus=modulus
        )
        if modulus:
            total %= modulus
    return total if modulus is None else total % modulus


def audit_criterion(cases):
    checks = 0
    for n, prime in cases:
        assert prime * 2 > n and prime <= n
        divides = apery(n) % prime == 0
        value = shell_sum(n - 1, prime - 1, modulus=prime)
        assert divides == (value == 0), (n, prime, divides, value)
        checks += 1
    return checks


def operator(moment, s):
    return (
        (s - moment) ** 3,
        2 * s ** 3 - 4 * (moment - 1) * s ** 2
        + (moment ** 2 - 7 * moment + 3) * s
        + (2 * moment ** 3 + 4 * moment ** 2 - 2 * moment + 1),
        s ** 3 - (moment - 4) * s ** 2
        - (moment + 5) * (moment - 1) * s
        + (moment - 1) ** 2 * (moment + 2),
    )


def audit_operator(moments=(20, 30, 41, 55), width=25):
    checks = 0
    for moment in moments:
        for s in range(1, width):
            q0, q1, q2 = operator(moment, s)
            lhs = (
                q2 * coefficient(moment, s + 2, 0, 0)
                + q1 * coefficient(moment, s + 1, 0, 0)
                + q0 * coefficient(moment, s, 0, 0)
            )
            assert lhs == 0, (moment, s, lhs)
            checks += 1
    return checks


if __name__ == "__main__":
    cases = (
        (200, 139), (200, 181), (200, 151), (200, 163),
        (321, 179), (321, 193), (321, 211), (321, 197),
        (60, 53), (60, 47), (100, 89), (100, 83),
    )
    print("CRITERION_CHECKS", audit_criterion(cases))
    print("SHIFT_OPERATOR_CHECKS", audit_operator())
    print("Q32_SHIFT_OPERATOR_AUDIT=PASS")
