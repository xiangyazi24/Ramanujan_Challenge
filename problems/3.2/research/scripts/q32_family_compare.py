#!/usr/bin/env python3
"""Compare R(n) = log rad_{p<=n}(u_n) across Apery-like sequences.

If the Apery zeta(3) numbers are special, this shows in which company they sit.
Control case: C(2n,n), which is built entirely from primes <= 2n, so its R(n) is ~ 2n --
the conjecture-analogue is FALSE there, and provably so.

Recurrences are handled division-free: for L(n) u_{n+1} = A(n) u_n - B(n) u_{n-1},
set U_n = u_n * prod_{j<n} L(j); then U_{n+1} = A(n) U_n - B(n) L(n-1) U_{n-1},
and for p > n, p | u_n iff p | U_n.
"""
import sys
from math import comb, log

SEQS = {
    # name: (L, A, B, u0, u1)   with u_{n+1} = (A(n) u_n - B(n) u_{n-1}) / L(n)
    "apery_zeta3": (lambda n: (n + 1) ** 3,
                    lambda n: 34 * n ** 3 + 51 * n ** 2 + 27 * n + 5,
                    lambda n: n ** 3, 1, 5),
    "apery_zeta2": (lambda n: (n + 1) ** 2,
                    lambda n: 11 * n ** 2 + 11 * n + 3,
                    lambda n: -n ** 2, 1, 3),
    "franel": (lambda n: (n + 1) ** 2,
               lambda n: 7 * n ** 2 + 7 * n + 2,
               lambda n: -8 * n ** 2, 1, 2),
    "domb": (lambda n: (n + 1) ** 3,
             lambda n: 2 * (2 * n + 1) * (5 * n ** 2 + 5 * n + 2),
             lambda n: -64 * n ** 3, 1, 4),
    "almkvist_zudilin": (lambda n: (n + 1) ** 3,
                         lambda n: (2 * n + 1) * (3 * n ** 2 + 3 * n + 1),
                         lambda n: -27 * n ** 3, 1, -3),
    "cooper_s7": (lambda n: (n + 1) ** 3,
                  lambda n: 26 * n ** 3 + 39 * n ** 2 + 21 * n + 4,
                  lambda n: -27 * n ** 3 + 3 * n, 1, 4),
    "cooper_s10": (lambda n: (n + 1) ** 3,
                   lambda n: 12 * n ** 3 + 18 * n ** 2 + 10 * n + 2,
                   lambda n: -64 * n ** 3 + 4 * n, 1, 2),
    "cooper_s18": (lambda n: (n + 1) ** 3,
                   lambda n: 28 * n ** 3 + 42 * n ** 2 + 26 * n + 6,
                   lambda n: 192 * n ** 3 - 12 * n, 1, 6),
}


def primes_upto(n):
    s = bytearray([1]) * (n + 1)
    s[0:2] = b'\x00\x00'
    for i in range(2, int(n ** .5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]


def zero_set(name, p):
    """Z_p = {r < p : p | u_r}, via the scaled recurrence (valid only for r < p)."""
    L, A, B, u0, u1 = SEQS[name]
    U0, U1 = u0 % p, (u1 * L(0)) % p
    Z = set()
    if U0 == 0:
        Z.add(0)
    if U1 == 0:
        Z.add(1)
    for n in range(1, p - 1):
        Un = (A(n) % p * U1 - B(n) % p * (L(n - 1) % p) % p * U0) % p
        if Un == 0:
            Z.add(n + 1)
        U0, U1 = U1, Un
    return Z


def divisibility_rows(name, N, ps):
    """p | u_n iff some base-p digit of n lies in Z_p (Lucas property, Malik-Straub)."""
    hits = {}
    for p in ps:
        Z = zero_set(name, p)
        s = set()
        if Z:
            for n in range(p, N + 1):
                m = n
                while m:
                    if m % p in Z:
                        s.add(n)
                        break
                    m //= p
        hits[p] = s
    return hits


def central_binomial_rows(N, ps):
    hits = {}
    for p in ps:
        s = set()
        for n in range(N + 1):
            # Kummer: p | C(2n,n) iff adding n+n in base p has a carry
            m, carry = n, False
            while m:
                if 2 * (m % p) >= p:
                    carry = True
                    break
                m //= p
            if carry:
                s.add(n)
        hits[p] = s
    return hits


def report(name, hits, N, ps):
    R = [0.0] * (N + 1)
    W = [0] * (N + 1)
    for p in ps:
        lp = log(p)
        for n in hits[p]:
            if n <= N and p <= n:
                R[n] += lp
                W[n] += 1
    out = []
    for lo, hi in ((100, 500), (500, 2000)):
        best = max(range(lo, min(hi, N) + 1), key=lambda n: R[n] / n)
        out.append(f"n in [{lo},{hi}]: max R/n = {R[best]/best:.4f} (n={best}, R={R[best]:.1f}), "
                   f"max #primes = {max(W[lo:min(hi,N)+1])}")
    print(f"{name:20s} " + " | ".join(out), flush=True)


if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
    ps = [p for p in primes_upto(N) if p >= 5]
    for name in SEQS:
        report(name, divisibility_rows(name, N, ps), N, ps)
    report("central_binomial", central_binomial_rows(N, ps), N, ps)
