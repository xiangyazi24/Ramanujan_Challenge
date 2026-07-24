#!/usr/bin/env python3
"""Test the first endpoint-asymptotic cancellation for fixed odd q.

With K=J+1, S=sum_{k<K} L(n,k)F_k, and T=L(n,K)F_K, the endpoint term ratio
has limit rho_q=32*q*(q+1).  Hence S/T tends to 1/(rho_q-1).  The integer
certificate

    E=(rho_q-1)S-T

cancels the leading constant but, empirically and asymptotically, only gains
a polynomial factor.  Every fixed-q bad prime still divides E because it
divides both S and the universal carrier T.
"""

from __future__ import annotations

from math import comb, log

from q32_strehl_gcd import franel_numbers


INDICES = (120, 240, 480, 960, 1200)


def main() -> None:
    franel = franel_numbers(max(INDICES))
    for quotient in (1, 3):
        limiting_ratio = 32 * quotient * (quotient + 1)
        print(f"q={quotient} rho={limiting_ratio}")
        for n in INDICES:
            cutoff = (n - quotient) // (2 * quotient + 1)
            boundary = cutoff + 1
            strehl = sum(
                comb(n, k) * comb(n + k, k) * franel[k]
                for k in range(boundary)
            )
            endpoint = (
                comb(n, boundary)
                * comb(n + boundary, boundary)
                * franel[boundary]
            )
            residual = (limiting_ratio - 1) * strehl - endpoint
            print(
                f"n={n} J={cutoff} "
                f"log_abs_E/n={log(abs(residual))/n:.9f} "
                f"log_T/n={log(endpoint)/n:.9f}"
            )


if __name__ == "__main__":
    main()
