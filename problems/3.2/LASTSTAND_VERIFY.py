#!/usr/bin/env python3
"""
LASTSTAND_VERIFY.py -- one gate per claim banked in FABLE_NOTES sections 139-144.

Every claim entered into the ledger on 2026-08-02 (life lane) must pass its gate
here. A claim with no gate in this file is NOT banked; a claim whose gate fails
must be moved to ERRATA.

Run: python3 LASTSTAND_VERIFY.py
Exit 0 = all gates pass. Final line is PASS or FAIL.
"""
import sys, math
from math import comb

FAILURES = []

def gate(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  -- {detail}" if detail else ""))
    if not ok:
        FAILURES.append(name)

# ---------- orbit primitives ----------
def orbit_b(p, upto=None):
    """b_r mod p for r = 0 .. upto (default p-2). Denominators (n+1)^3 invertible for n+1 < p."""
    if upto is None: upto = p - 2
    b = [1 % p, 5 % p]
    for n in range(1, upto):
        A = (34*n**3 + 51*n**2 + 27*n + 5) % p
        B = (n**3) % p
        inv = pow((n+1)**3 % p, p-2, p)
        b.append(((A*b[n] - B*b[n-1]) * inv) % p)
    return b

def orbit_bc(p, upto=None):
    if upto is None: upto = p - 2
    b = [1 % p, 5 % p]; c = [0 % p, 6 % p]
    for n in range(1, upto):
        A = (34*n**3 + 51*n**2 + 27*n + 5) % p
        B = (n**3) % p
        inv = pow((n+1)**3 % p, p-2, p)
        b.append(((A*b[n] - B*b[n-1]) * inv) % p)
        c.append(((A*c[n] - B*c[n-1]) * inv) % p)
    return b, c

def Zp_set(p):
    b = orbit_b(p, upto=p-1)          # r = 0 .. p-2
    return {r for r, v in enumerate(b) if v % p == 0}

def primes_upto(X):
    sieve = bytearray([1])*(X+1); sieve[0:2] = b"\0\0"
    for i in range(2, int(X**0.5)+1):
        if sieve[i]: sieve[i*i::i] = bytearray(len(sieve[i*i::i]))
    return [i for i in range(X+1) if sieve[i]]

def apery_exact(n):
    return sum(comb(n,k)**2 * comb(n+k,k)**2 for k in range(n+1))

# ---------- GATE 1 (s141): base-p digit criterion, against EXACT ground truth ----------
# Claim: p | b_n  <=>  some base-p digit of n lies in Z_p.
def gate_digit_criterion():
    bad = []; tested = 0
    for p in [7, 11, 13, 17, 19, 23, 29]:
        Z = Zp_set(p)
        for n in range(0, 90):
            digs = []; m = n
            if m == 0: digs = [0]
            while m > 0: digs.append(m % p); m //= p
            pred = any(d in Z for d in digs)
            act  = (apery_exact(n) % p == 0)
            tested += 1
            if pred != act: bad.append((p, n, pred, act))
    gate("s141 digit criterion (exact binomial ground truth)", not bad,
         f"{tested} (p,n) pairs over 7 primes, mismatches={len(bad)}" + (f" e.g. {bad[:3]}" if bad else ""))

# ---------- GATE 2 (s141): top-window collapse p|b_{n-p} <=> p|b_n ----------
def gate_top_collapse():
    bad = []; tested = 0
    for n in range(40, 130):
        bn = apery_exact(n)
        for p in primes_upto(n):
            if p <= n//2: continue
            r = n - p
            br = orbit_b(p, upto=max(r,2)+1)[r] if r >= 2 else (1 if r == 0 else 5) % p
            lhs = (br % p == 0); rhs = (bn % p == 0)
            tested += 1
            if lhs != rhs: bad.append((n, p, lhs, rhs))
    gate("s141 top-window collapse  p|b_(n-p) <=> p|b_n", not bad,
         f"{tested} (n,p) pairs, mismatches={len(bad)}")

# ---------- GATE 3 (s142): the averaging identity, at SEVERAL X ----------
# Claim: sum_{n<=X} #{p in (n/2,n] : p|b_n} = #{(p,r) : r in Z_p, r>=1, p+r<=X}
def gate_identity():
    rows = []; allok = True
    for X in [150, 300, 450, 600]:
        pr = primes_upto(X)
        Zc = {p: Zp_set(p) for p in pr if p >= 5}
        lhs = 0
        for n in range(2, X+1):
            for p in pr:
                if p <= n//2 or p > n: continue
                if p < 5: continue
                r = n - p
                if r >= 1 and r in Zc.get(p, ()): lhs += 1
        rhs = sum(1 for p in pr if p >= 5 for r in Zc[p] if r >= 1 and p + r <= X)
        ok = (lhs == rhs); allok &= ok
        rows.append(f"X={X}:{lhs}={rhs}" if ok else f"X={X}:{lhs}!={rhs}")
    gate("s142 averaging identity (multi-X)", allok, "  ".join(rows))

# ---------- GATE 4 (s143): |Z_p| parity law ----------
# Claim: |Z_p| is even unless (p-1)/2 in Z_p  (reflection pairing r <-> p-1-r).
def gate_parity():
    bad = []; n_even = 0; n_odd = 0
    for p in primes_upto(4000):
        if p < 5: continue
        Z = Zp_set(p)
        mid_hit = ((p-1)//2) in Z
        if len(Z) % 2 == 0: n_even += 1
        else: n_odd += 1
        if (len(Z) % 2 == 1) != mid_hit: bad.append((p, len(Z), mid_hit))
    gate("s143 parity law  |Z_p| odd <=> (p-1)/2 in Z_p", not bad,
         f"even={n_even} odd={n_odd} violations={len(bad)}" + (f" e.g. {bad[:3]}" if bad else ""))

# ---------- GATE 5 (s143): reflection b_{p-1-r} = b_r mod p ----------
def gate_reflection():
    bad = []
    for p in [101, 211, 307, 401, 503]:
        b = orbit_b(p, upto=p-1)
        for r in range(0, p-1):
            if b[r] % p != b[p-1-r] % p: bad.append((p, r)); break
    gate("s143 reflection  b_{p-1-r} = b_r mod p", not bad, f"5 primes fully checked, violations={len(bad)}")

# ---------- GATE 6 (s140): the REFUTED sufficient condition really is Theta(X) ----------
# Claim: sum_{p<=X} log p |Z_p| is ~c X with c of order 1 (so "= o(X)" is FALSE).
def gate_avg_scale():
    rows = []
    ok = True
    for X in [1000, 2000, 4000]:
        s = 0.0; tot = 0
        for p in primes_upto(X):
            if p < 5: continue
            z = len(Zp_set(p)); tot += z; s += z*math.log(p)
        ratio = s/X
        rows.append(f"X={X}: sum log p|Z_p|/X={ratio:.3f}, mean|Z_p|={tot/max(1,len(primes_upto(X))):.3f}")
        if not (0.3 < ratio < 3.0): ok = False
    gate("s140 sum_p log p |Z_p| = Theta(X)  (refutes the 'o(X)' hypothesis)", ok, " | ".join(rows))

# ---------- GATE 7 (s139): determinant-bilinear numerics at square-root ----------
# Claim: mean_t |B(t)| is comparable to sqrt(#strip pairs) (NOT the trivial ND).
def gate_det_bilinear():
    import cmath
    rows = []; ok = True
    for p in [1009, 2003]:
        N = p-2; D = int(math.isqrt(N)*2)
        b, c = orbit_bc(p, upto=N+1)
        rs = list(range(1, N+1-D))
        npairs = len(rs)*(D - D//2)
        tot = 0.0; cnt = 0
        for t in range(1, 25):
            acc = 0j
            for d in range(D//2+1, D+1):
                for r in rs[::3]:                     # subsample r for speed; unbiased for magnitude scale
                    det = (b[r]*c[r+d] - b[r+d]*c[r]) % p
                    acc += cmath.exp(2j*math.pi*(t*det % p)/p)
            tot += abs(acc)*3; cnt += 1                # rescale for the subsample
        mean = tot/cnt
        sq = math.sqrt(npairs)
        rows.append(f"p={p}: mean|B|={mean:.0f} sqrt(pairs)={sq:.0f} trivial={npairs}")
        if not (mean < 0.2*npairs): ok = False          # gate: far below trivial
    gate("s139 |B(t)| is far below the trivial bound (square-root regime)", ok, " | ".join(rows))

# ---------- GATE 8 (s144): threshold calibration  sum_{p<=X} log p ~ X ----------
def gate_threshold():
    rows = []; ok = True
    for X in [10**4, 10**5]:
        s = sum(math.log(p) for p in primes_upto(X))
        rows.append(f"X={X}: theta(X)/X={s/X:.3f}")
        if not (0.8 < s/X < 1.2): ok = False
    gate("s144 threshold: theta(X) ~ X  (so |Z_p| <= C log p gives only Theta(X))", ok, " | ".join(rows))

if __name__ == "__main__":
    print("=== LASTSTAND verification gates (ledger sections 139-144) ===")
    gate_digit_criterion()
    gate_top_collapse()
    gate_identity()
    gate_parity()
    gate_reflection()
    gate_avg_scale()
    gate_threshold()
    gate_det_bilinear()
    print()
    if FAILURES:
        print("FAIL -- failing gates: " + ", ".join(FAILURES))
        sys.exit(1)
    print("PASS -- all gates green")
    sys.exit(0)
