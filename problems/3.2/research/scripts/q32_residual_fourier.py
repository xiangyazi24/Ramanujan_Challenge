#!/usr/bin/env python3
"""Residual Fourier spectrum of Apery residues (Q6441 sec.5 experiment).

S_p(u) = sum_{r=0}^{p-2} e(u * (b_r mod p) / p),  u in F_p  — computed as the
FFT of the residue histogram. Exact identity: |Z_p| = (1/p) sum_u S_p(u).
Measured: max_{u!=0} |S_p(u)| / sqrt(p) (pointwise square-root cancellation
test for [GAP-RES-WEYL] / [RES-SQRT]), the identity check, and Poisson tally.
"""
import math

import numpy as np


def primes_in(lo, hi):
    sieve = bytearray([1]) * (hi + 1)
    sieve[0:2] = b"\x00\x00"
    for i in range(2, int(hi**0.5) + 1):
        if sieve[i]:
            sieve[i * i :: i] = bytearray(len(sieve[i * i :: i]))
    return [i for i in range(lo, hi + 1) if sieve[i]]


def apery_residues(p):
    """b_r mod p for r = 0..p-2 via the recurrence (exact integer arithmetic)."""
    res = [1, 5 % p]
    for n in range(1, p - 2):
        num = (((34 * n + 51) * n + 27) * n + 5) % p
        b1 = (num * res[n] - n**3 * res[n - 1]) * pow((n + 1) ** 3, p - 2, p) % p
        res.append(b1)
    return res[: p - 1]


def spectrum(p):
    res = apery_residues(p)
    hist = np.bincount(np.array(res, dtype=np.int64), minlength=p).astype(float)
    zp = int(hist[0])
    S = np.fft.fft(hist)  # S[u] = sum_c hist[c] e(-uc/p); |S| symmetric in sign convention
    ident = S.real.sum() / p
    mags = np.abs(S[1:])
    return zp, ident, mags.max() / math.sqrt(p), mags.mean() / math.sqrt(p)


def main():
    ps = primes_in(1000, 30000)
    print(f"testing {len(ps)} primes in [1000, 30000]")
    tally = {}
    bad = 0
    rows = []
    for p in ps:
        zp, ident, mxn, avgn = spectrum(p)
        if abs(ident - zp) > 1e-5:
            bad += 1
            print(f"IDENTITY FAIL p={p}: {ident} vs {zp}")
        tally[zp] = tally.get(zp, 0) + 1
        rows.append((p, zp, mxn, avgn))
    n = len(ps)
    mxs = [r[2] for r in rows]
    avgs = [r[3] for r in rows]
    print(f"identity failures: {bad}/{n}")
    print("Z_p tally:", dict(sorted(tally.items())))
    lam = sum(k * v for k, v in tally.items()) / n
    print(f"mean |Z_p| = {lam:.4f}")
    print(f"max|S|/sqrt(p): mean {sum(mxs)/n:.3f}, min {min(mxs):.3f}, max {max(mxs):.3f}")
    print(f"avg|S|/sqrt(p): mean {sum(avgs)/n:.4f}")
    # growth check: does max|S|/sqrt(p) grow with p? bucket by size
    buckets = [(1000, 4000), (4000, 10000), (10000, 20000), (20000, 30000)]
    for lo, hi in buckets:
        sel = [r[2] for r in rows if lo <= r[0] < hi]
        if sel:
            print(f"  p in [{lo},{hi}): mean max|S|/sqrt(p) = {sum(sel)/len(sel):.3f} (n={len(sel)})")
    # sqrt(log p) comparison (Gaussian max heuristic over p frequencies: ~sqrt(log p))
    for lo, hi in buckets:
        sel = [r[2] / math.sqrt(math.log(r[0])) for r in rows if lo <= r[0] < hi]
        if sel:
            print(f"  p in [{lo},{hi}): mean max|S|/sqrt(p log p) = {sum(sel)/len(sel):.3f}")
    print("VERIFIED" if bad == 0 else "FAILED")


if __name__ == "__main__":
    main()
