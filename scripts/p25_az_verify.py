#!/usr/bin/env python3
"""
P2.5: Extract and verify the AZ certificate found by p25_az_correct.py.

AUDIT NOTE (2026-07-30)
-----------------------
This program verifies only the cleared differential identity

    Phi P_(n+1) - P_n M_H = D_u A + D_v B + D_t C.

It does *not* turn that identity into a Delannoy constant-term period.  The
covariant derivatives in p25_az_correct.py are integrated against the ordinary
form ``du dv dt``.  After the factor ``u v (1+t^2)`` in P cancels the carrier
denominator, small u- and v-contours extract the coefficient of
``u^-1 v^-1`` rather than the constant term.  At eps=0 the remaining t
integrand is constant, so its ordinary closed-contour integral is zero.  Thus
the previously proposed logarithmic-measure interpretation (``du/u``,
``dv/v``) changes the functional and is invalid.

An independent exact recurrence check also rules out the claimed carrier row:
the proposed r_0 = [1, 1, 1/4] gives
r_0 M_H(0) = [15/8, -11/8, 1/4], not
[81/4, 27/2, 9/4].  Do not cite this certificate as a connection coefficient
or a proof of the Catalan limit.

The search found P_n = (1/4) * uv * (1+t^2) * [(n+2)^2, 2(n+2), 1]
at n-degree 2, P support [0,2]^3, for ALL values of eps.

This script:
1. Re-runs the modular solve with eps=0 (concrete, simplest)
2. Extracts the full solution vector (P, A, B, C)
3. Reconstructs rational coefficients via CRT + rational reconstruction
4. Verifies the identity EXACTLY over Z at sufficiently many n values
5. Reports the full certificate in human-readable form
"""

import os
import sys
import time
import numpy as np
from fractions import Fraction
from math import gcd, isqrt

# Import the search module
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p25_az_correct import (
    solve_config, full_vector_for_P, build_terms, make_sup,
    Mmat, deltaH, PHIG, G, GVGT, GUGT, GUGV, GU, GV, GT,
    U1U2, V1V2, U2M2, V2M2, ONE_T2,
    lmul, ladd, lscale, lshift, tag,
    ratrec, crt,
    PRIME1, PRIME2,
)


def main():
    t0 = time.time()
    Pu, Pv, Pt = (0, 2), (0, 2), (0, 2)
    dP, dC = 2, 2
    eps_val = 0  # concrete eps=0

    print("=" * 78)
    print("P2.5 AZ certificate: extraction and exact verification")
    print(f"  dP={dP}, dC={dC}, eps={eps_val}")
    print(f"  P support: u in {Pu}, v in {Pv}, t in {Pt}")
    print("=" * 78)

    # Solve with two primes
    print("\n--- Solving mod p1 ---")
    res1 = solve_config(dP, dC, str(eps_val), PRIME1, seed=1,
                        Pu=Pu, Pv=Pv, Pt=Pt)
    print("\n--- Solving mod p2 ---")
    res2 = solve_config(dP, dC, str(eps_val), PRIME2, seed=1,
                        Pu=Pu, Pv=Pv, Pt=Pt)

    if res1['Pdim'] == 0:
        print("ERROR: no P found. Exiting.")
        return
    if res1['Pdim'] != res2['Pdim']:
        print("WARNING: P-dim mismatch between primes!")

    print(f"\nP-dim = {res1['Pdim']} (both primes agree: {res1['Pdim'] == res2['Pdim']})")

    # Extract first P basis vector, lift to full solution in both primes
    prow1 = res1['Pbasis'][0]
    prow2 = res2['Pbasis'][0]
    v1 = full_vector_for_P(res1, prow1)
    v2 = full_vector_for_P(res2, prow2)

    if v1 is None or v2 is None:
        print("ERROR: could not lift P to full certificate vector")
        return

    N = len(v1)
    p1, p2 = PRIME1, PRIME2
    pp = p1 * p2  # ~4.6 * 10^18
    bound = isqrt(pp // 2)

    print(f"\nRational reconstruction (CRT with p1*p2 = {pp:.2e}, bound = {bound:.2e})...")

    # CRT + rational reconstruction
    xr = []
    fail_count = 0
    for i in range(N):
        a = crt(int(v1[i]), int(v2[i]), p1, p2)
        r = ratrec(a, pp, bound)
        if r is None:
            fail_count += 1
            xr.append(None)
        else:
            xr.append(r)

    if fail_count > 0:
        print(f"  {fail_count} / {N} entries failed rational reconstruction!")
        # Try with third prime
        print("  Trying third prime for a wider CRT product...")
        PRIME3 = 2147483587
        res3 = solve_config(dP, dC, str(eps_val), PRIME3, seed=1,
                            Pu=Pu, Pv=Pv, Pt=Pt)
        prow3 = res3['Pbasis'][0]
        v3 = full_vector_for_P(res3, prow3)
        # Use all three primes
        ppp = p1 * p2 * PRIME3
        bound3 = isqrt(ppp // 2)
        xr = []
        for i in range(N):
            a12 = crt(int(v1[i]), int(v2[i]), p1, p2)
            a123 = crt(a12, int(v3[i]), pp, PRIME3)
            r = ratrec(a123, ppp, bound3)
            if r is None:
                fail_count += 1
                xr.append(None)
            else:
                xr.append(r)
        if any(x is None for x in xr):
            print(f"  STILL failed after 3 primes. Falling back to modular report.")
            report_modular(res1, prow1)
            return

    print(f"  All {N} entries reconstructed successfully!")

    # Clear denominators
    den = 1
    for r in xr:
        den = den * r.denominator // gcd(den, r.denominator)
    xint = [int(r * den) for r in xr]
    print(f"  Common denominator: {den}")
    maxabs = max(abs(x) for x in xint)
    print(f"  Max absolute coefficient: {maxabs}")

    # Exact verification
    unknowns = res1['unknowns']
    supP = res1['supP']
    supA, supB, supC = res1['supA'], res1['supB'], res1['supC']
    _, terms = build_terms(dP, dC, str(eps_val), eps_val, supP, supA, supB, supC)

    maxk = max(dP, dC)
    ndeg = 8 + maxk
    print(f"\nExact verification: checking identity at n = 1..{ndeg + 3}...")

    bad = 0
    for n0 in range(1, ndeg + 4):
        dH = deltaH(n0)
        Mn = Mmat(n0)
        sc = {}
        for k in range(maxk + 1):
            sc[('phi', k)] = dH * (n0 + 1) ** k
            sc[('c0', k)] = -dH * n0 ** k
            sc[('c1', k)] = -dH * n0 ** (k + 1)
            for i in range(3):
                for j in range(3):
                    sc[('m', k, i, j)] = -(n0 ** k) * Mn[i][j]
        acc = {}
        for col, tl in enumerate(terms):
            x = xint[col]
            if x == 0:
                continue
            for sid, d in tl:
                s = sc[sid] * x
                if s == 0:
                    continue
                for key, v in d.items():
                    acc[key] = acc.get(key, 0) + s * v
        residual = {k: v for k, v in acc.items() if v}
        if residual:
            bad += 1
            if bad <= 3:
                print(f"  n={n0}: NONZERO residual ({len(residual)} terms)")
                for k in sorted(residual.keys())[:5]:
                    print(f"    {k}: {residual[k]}")

    if bad == 0:
        print(f"  EXACT VERIFICATION PASSED at all {ndeg + 3} values of n!")
        print(f"  (Degree bound {ndeg} for n-identity => proven.)")
    else:
        print(f"  FAILED: {bad} nonzero residuals")
        return

    # Report the certificate
    print("\n" + "=" * 78)
    print("CERTIFICATE STRUCTURE")
    print("=" * 78)

    blocks = {}
    for idx, (blk, k, comp, mono) in enumerate(unknowns):
        if xint[idx] == 0:
            continue
        r = xr[idx]
        if blk not in blocks:
            blocks[blk] = {}
        if comp not in blocks[blk]:
            blocks[blk][comp] = []
        blocks[blk][comp].append((k, mono, r))

    for blk_name in ['P', 'A', 'B', 'C']:
        if blk_name not in blocks:
            print(f"\n{blk_name}: zero")
            continue
        print(f"\n{blk_name}:")
        for comp in range(3):
            if comp not in blocks[blk_name]:
                print(f"  [{comp}]: 0")
                continue
            terms_str = []
            for k, (a, b, c), r in sorted(blocks[blk_name][comp]):
                coeff = str(r) if r.denominator != 1 else str(r.numerator)
                mono_parts = []
                if k > 0:
                    mono_parts.append(f"n^{k}" if k > 1 else "n")
                for s, e in [("u", a), ("v", b), ("t", c)]:
                    if e != 0:
                        mono_parts.append(f"{s}^{e}" if abs(e) != 1 else
                                         (s if e == 1 else f"{s}^{e}"))
                mono_str = "*".join(mono_parts) if mono_parts else "1"
                terms_str.append(f"{coeff}*{mono_str}")
            print(f"  [{comp}]: {' + '.join(terms_str)}")

    # Summary
    print("\n" + "=" * 78)
    print("SUMMARY")
    print("=" * 78)
    print(f"  Phase: Phi(u,v) = phi(u)*phi(v), phi(u) = (1+u)(u+2)/u")
    print(f"  P_n = (1/4) * uv * (1+t^2) * [(n+2)^2, 2(n+2), 1]")
    print(f"  eps = {eps_val}")
    print(f"  n-degree of P: {dP}")
    print(f"  n-degree of certificates: {dC}")
    print(f"  Exact verification: PASSED")
    print(f"  Total time: {time.time()-t0:.1f}s")
    print("=" * 78)


def report_modular(res, prow):
    """Fallback: report the P structure modularly."""
    from p25_az_correct import describe_P
    print("\n--- Modular P basis (rational reconstruction partial) ---")
    for line in describe_P(res, prow):
        print(line)


if __name__ == '__main__':
    main()
