#!/usr/bin/env sage -python
"""Exact verifier for the Q7703 fixed-exponent horizontal reduction.

For fixed n and a quotient cell q, primes in the cell satisfy
    n = q*p + r,  1 <= r <= p-2.
Put s=min(r,p-1-r) and
    A_n(k) = (binomial(n,k)*binomial(n+k,k))^2,
    S_n(s) = sum_{0<=k<=s} A_n(k).
The exact reductions verified here are
    b_r == S_n(s) == b_s  (mod p),
and
    p=(n-s)/q                 if 2r <= p-1,
    p=(n+s+1)/(q+1)           if 2r >  p-1.
Thus the defining-prime Mellin-zero condition is exactly one of two affine
prime-divisor rays for the fixed characteristic-zero Apéry sequence b_s.

Run, e.g.
    sage -python problems/3.2/research/scripts/q7703_fixed_exponent_prefix.sage \
        --n 321 --qmax 8
"""

from sage.all import ZZ, binomial, prime_range
import argparse


def apery(r):
    r = ZZ(r)
    return ZZ(sum((binomial(r, k) * binomial(r + k, k)) ** 2
                  for k in range(int(r) + 1)))


def fixed_term(n, k):
    n, k = ZZ(n), ZZ(k)
    return ZZ((binomial(n, k) * binomial(n + k, k)) ** 2)


def fixed_prefixes(n, smax):
    """Return S_n(0),...,S_n(smax) as exact integers."""
    n = ZZ(n)
    out = []
    total = ZZ(0)
    for k in range(smax + 1):
        total += fixed_term(n, k)
        out.append(total)
    return out


def check_prefix_recurrence(n, prefixes):
    """Verify the fixed degree-4 recurrence for S_n(s).

    If T_s=S_s-S_{s-1}, then
       s^4 T_s = ((n-s+1)(n+s))^2 T_{s-1}.
    """
    n = ZZ(n)
    if len(prefixes) < 3:
        return
    for s in range(2, len(prefixes)):
        S0 = ZZ(prefixes[s - 2])
        S1 = ZZ(prefixes[s - 1])
        S2 = ZZ(prefixes[s])
        R = ((n - s + 1) * (n + s)) ** 2
        assert ZZ(s) ** 4 * (S2 - S1) == R * (S1 - S0)


def cell_primes(n, q):
    """All primes p with floor(n/p)=q, endpoints included for filtering."""
    n, q = ZZ(n), ZZ(q)
    lo = n // (q + 1) + 1
    hi = n // q
    if hi < lo:
        return []
    return [ZZ(p) for p in prime_range(lo, hi + 1)]


def scan(n, qmax):
    n = ZZ(n)
    qmax = int(qmax)

    # s <= (p-1)/2 <= n/2 in every interior cell, so this one fixed prefix
    # table is enough for all requested q.
    smax = int(n // 2)
    prefixes = fixed_prefixes(n, smax)
    check_prefix_recurrence(n, prefixes)

    total_rows = 0
    total_bad = 0
    for q in range(1, qmax + 1):
        rows = []
        for p in cell_primes(n, q):
            r = n - q * p
            if not (1 <= r <= p - 2):
                continue
            s = min(r, p - 1 - r)
            S = ZZ(prefixes[int(s)])
            br = apery(r)
            bs = apery(s)

            # Fixed-exponent prefix theorem and folded Apéry reflection.
            assert (S - br) % p == 0
            assert (S - bs) % p == 0
            assert (br - bs) % p == 0

            if 2 * r <= p - 1:
                branch = "left"
                assert s == r
                assert (n - s) % q == 0
                assert (n - s) // q == p
                # qp | q*b_s iff p | b_s.
                assert ((q * bs) % (n - s) == 0) == (bs % p == 0)
                assert (2 * q + 1) * s <= n - q
            else:
                branch = "right"
                assert s == p - 1 - r
                assert (n + s + 1) % (q + 1) == 0
                assert (n + s + 1) // (q + 1) == p
                # (q+1)p | (q+1)b_s iff p | b_s.
                assert (((q + 1) * bs) % (n + s + 1) == 0) == (bs % p == 0)
                assert (2 * q + 1) * s <= n - 2 * q - 1

            # The exponent in the Mellin formulation is fixed inside the q-cell.
            m = n - q
            assert m % (p - 1) == r

            bad = int(bs % p == 0)
            rows.append((int(p), int(r), int(s), branch, bad))

        bad_rows = [row for row in rows if row[-1]]
        print(
            "CELL",
            q,
            "prime_rows",
            len(rows),
            "bad",
            len(bad_rows),
            "bad_data",
            bad_rows,
        )
        total_rows += len(rows)
        total_bad += len(bad_rows)

    print("TOTAL_ROWS", total_rows)
    print("TOTAL_BAD", total_bad)
    print("Q7703_FIXED_PREFIX_VERIFY PASS")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=321)
    ap.add_argument("--qmax", type=int, default=8)
    args = ap.parse_args()
    assert args.n >= 5
    assert args.qmax >= 1
    scan(args.n, args.qmax)


if __name__ == "__main__":
    main()
