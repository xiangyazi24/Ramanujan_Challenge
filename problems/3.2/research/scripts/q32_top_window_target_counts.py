#!/usr/bin/env python3
"""K(n) and T(n): the quantity that actually has to be bounded: for each n,

    T(n)  = sum over primes p in (n/2, n] with p | b_{n-p} of log p,
    K(n)  = the number of such primes ("targets at level n").

The Apery-Lucas step says p | b_n <=> p | b_{n-p} for p in that range, so T(n) is exactly
the top-half sum, and the open statement is T(n) = o(n).

Computes Z_p = {r < p : p | b_r} for every prime p <= N by iterating the Apery recurrence
mod p, then assembles T(n), K(n) for all n <= N.  Prints progress.
"""
import sys, time
from math import log


def primes_upto(n):
    s = bytearray([1]) * (n + 1)
    s[0:2] = b'\x00\x00'
    for i in range(2, int(n ** .5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]


def zeros_mod(p):
    """{ r < p : p | b_r } by iterating the Apery recurrence mod p."""
    out = []
    inv = [0, 1] + [0] * (p - 2)
    for k in range(2, p):
        inv[k] = (p - (p // k) * inv[p % k] % p) % p
    b0, b1 = 1 % p, 5 % p
    if b0 == 0:
        out.append(0)
    if b1 == 0:
        out.append(1)
    for j in range(1, p - 1):
        ik = inv[j + 1]
        lead = ik * ik % p * ik % p
        nxt = ((34 * j ** 3 + 51 * j ** 2 + 27 * j + 5) % p * b1
               - j ** 3 % p * b0) % p * lead % p
        if nxt == 0:
            out.append(j + 1)
        b0, b1 = b1, nxt
    return out


def main(N):
    t0 = time.time()
    ps = primes_upto(N)
    print(f"# {len(ps)} primes up to {N}", flush=True)
    T = [0.0] * (N + 1)
    K = [0] * (N + 1)
    done = 0
    for p in ps:
        if p < 5:
            continue
        Z = zeros_mod(p)
        lp = log(p)
        for r in Z:
            n = p + r
            if n <= N and p > n / 2:
                T[n] += lp
                K[n] += 1
        done += 1
        if done % 200 == 0:
            print(f"  {done}/{len(ps)} primes, {time.time()-t0:.0f}s", flush=True)
    print(f"# assembled [{time.time()-t0:.0f}s]", flush=True)

    best = max(range(2, N + 1), key=lambda n: K[n])
    bestT = max(range(2, N + 1), key=lambda n: T[n] / max(n, 1))
    print(f"max K(n) = {K[best]} at n = {best}", flush=True)
    print(f"max T(n)/n = {T[bestT]/bestT:.5f} at n = {bestT} (T = {T[bestT]:.2f})", flush=True)
    for cut in (2, 3, 4, 5, 6):
        c = sum(1 for n in range(2, N + 1) if K[n] >= cut)
        print(f"count n<={N} with K(n)>={cut}: {c}", flush=True)
    # growth of the running max
    run = 0
    marks = []
    for n in range(2, N + 1):
        if K[n] > run:
            run = K[n]
            marks.append((n, run))
    print("running max of K(n):", marks, flush=True)
    # mean
    print(f"mean K = {sum(K)/N:.4f}, mean T/n = {sum(T[n]/n for n in range(2,N+1))/N:.6f}", flush=True)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 20000)
