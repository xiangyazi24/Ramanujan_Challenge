#!/usr/bin/env python3
"""Fable 2026-07-31: decide the value-law of b_m mod p (m=0..p-1).

Competing predictions for D(p) = #distinct{b_m mod p}:
  (DS)    D(p) ~ c*sqrt(p)          ["p^{3/2} folding" story]
  (Fable) D(p) ~ (1-e^{-1/2})*p     [uniform Poisson conditioned on the
                                     reflection FE b_r = b_{p-1-r}, no other structure]
Also check, against the "random mod FE" model:
  E(p)   = sum_a N(a)^2            -> 3p + O(sqrt p), fluct stdev ~ 2*sqrt(2p)
  |C_p(1)|                          -> Rayleigh, mean = sqrt(pi/4 * 2p) = 1.2533*sqrt(p)
  max_h |C_p(h)|                    -> Gumbel scale sqrt(2p * ln((p-1)/2))
"""
import math, sys
import numpy as np

PMAX = 1000

def primes_upto(n):
    s = np.ones(n + 1, dtype=bool); s[:2] = False
    for i in range(2, int(n**0.5) + 1):
        if s[i]: s[i*i::i] = False
    return np.nonzero(s)[0].tolist()

# exact Apery b_n (zeta(3) numbers) via integer recurrence
def apery_table(N):
    b = [1, 5]
    for n in range(1, N):
        num = (2*n + 1)*(17*n*n + 17*n + 5)*b[n] - n**3 * b[n-1]
        q, r = divmod(num, (n + 1)**3)
        assert r == 0, n
        b.append(q)
    return b

def main():
    ps = [p for p in primes_upto(PMAX) if p >= 11]
    b = apery_table(PMAX)
    print(f"{'p':>5} {'FE':>3} {'D(p)':>5} {'0.393p':>7} {'D/p':>6} {'D/sqrt(p)':>9} "
          f"{'(E-3p)/sq':>9} {'|C1|/sq':>8} {'maxC/gum':>8} {'maxmult':>7}")
    Ds, dp_ratios, dsq_ratios, efluct, c1s, gums = [], [], [], [], [], []
    for p in ps:
        v = np.array([b[m] % p for m in range(p)], dtype=np.int64)
        fe_ok = bool(np.all(v == v[::-1]))          # b_r = b_{p-1-r}
        hist = np.bincount(v, minlength=p).astype(float)
        D = int(np.count_nonzero(hist))
        E = float(np.dot(hist, hist))
        C = np.fft.fft(hist)                        # C[h] = sum_r N(r) e(-hr/p)
        mag = np.abs(C); mag[0] = 0.0
        c1 = mag[1] / math.sqrt(p)
        gum = math.sqrt(2 * p * math.log((p - 1) / 2))
        maxc = float(mag.max()) / gum
        ef = (E - 3 * p) / math.sqrt(p)
        Ds.append(D); dp_ratios.append(D / p); dsq_ratios.append(D / math.sqrt(p))
        efluct.append(ef); c1s.append(c1); gums.append(maxc)
        if p in (101, 199, 307, 401, 503, 601, 701, 797, 887, 997) or not fe_ok:
            print(f"{p:>5} {'ok' if fe_ok else 'FAIL':>3} {D:>5} {0.393*p:>7.1f} "
                  f"{D/p:>6.3f} {D/math.sqrt(p):>9.2f} {ef:>9.2f} {c1:>8.2f} "
                  f"{maxc:>8.2f} {int(hist.max()):>7}")
    n = len(ps)
    print(f"\n{n} primes in [11,{PMAX}]")
    print(f"D/p          : mean {np.mean(dp_ratios):.4f}  std {np.std(dp_ratios):.4f}   "
          f"[Fable pred (1-e^-1/2) = {1-math.exp(-0.5):.4f}]")
    print(f"D/sqrt(p)    : first {dsq_ratios[0]:.2f} ... last {dsq_ratios[-1]:.2f}  "
          f"(constant iff DS sqrt-law; drifting iff linear law)")
    print(f"(E-3p)/sq(p) : mean {np.mean(efluct):.2f}  std {np.std(efluct):.2f}   "
          f"[random-mod-FE pred: mean ~0, std ~ 2*sqrt(2) = {2*math.sqrt(2):.2f}]")
    print(f"|C_p(1)|/sq  : mean {np.mean(c1s):.3f}   [Rayleigh pred sqrt(pi/2)/sq2*2 = "
          f"{math.sqrt(math.pi/4*2):.4f}]")
    print(f"max_h/gumbel : mean {np.mean(gums):.3f}  std {np.std(gums):.3f}   [pred ~1]")

if __name__ == "__main__":
    main()
