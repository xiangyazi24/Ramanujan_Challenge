#!/usr/bin/env python3
"""
Systematic Barnes intercept search for P2.7.

For each tuple (a1, a2, a3, beta) of half-integer parameters,
compute the "raw" residue sum S_n of
  R_n(t) = prod_j (t-n+a_j)_n / (t+beta)_{n+1}
and check if S_1/S_0 and S_2/S_0 match q̂_1/q̂_0 and q̂_2/q̂_0.

The normalization kappa(n) cancels in ratios.
"""
from fractions import Fraction as F
from itertools import combinations_with_replacement
import sys

def pochhammer(x, n):
    """(x)_n = x(x+1)...(x+n-1), exact rational."""
    result = F(1)
    for i in range(n):
        result *= (x + i)
    return result

def raw_residue_sum(n, a1, a2, a3, beta):
    """
    Sum of residues of prod_j (t-n+a_j)_n / (t+beta)_{n+1}
    at t = -beta, -beta-1, ..., -beta-n.

    Residue at t = -beta-k is:
      prod_j ((-beta-k-n+a_j)_n) * (-1)^k / (k! * (n-k)!)
    """
    total = F(0)
    for k in range(n + 1):
        # Numerator factor: product of three Pochhammer symbols
        num = F(1)
        for a in [a1, a2, a3]:
            # (-beta-k-n+a)_n = product_{m=0}^{n-1} (-beta-k-n+a+m)
            #                 = product_{m=0}^{n-1} (a - beta - k - n + m)
            poch_val = pochhammer(a - beta - k - n, n)
            num *= poch_val

        # Denominator: (-1)^k / (k! * (n-k)!)
        sign = (-1)**k
        factorial_k = F(1)
        for i in range(1, k + 1):
            factorial_k *= i
        factorial_nk = F(1)
        for i in range(1, n - k + 1):
            factorial_nk *= i

        term = num * sign / (factorial_k * factorial_nk)
        total += term

    return total

# Target ratios
q0 = F(-215040420000)
q1 = F(-167282265043404, 905)
q2 = F(-964185327658080, 6071)

qhat0 = q0  # 64^0 = 1
qhat1 = F(64) * q1
qhat2 = F(64)**2 * q2

target_r1 = qhat1 / qhat0  # q̂_1 / q̂_0
target_r2 = qhat2 / qhat0  # q̂_2 / q̂_0

print(f"Target ratios:")
print(f"  q̂_1/q̂_0 = {float(target_r1):.15e}")
print(f"  q̂_2/q̂_0 = {float(target_r2):.15e}")
print()

# Search over half-integer parameters
# a1, a2, a3 in {-2, -3/2, -1, -1/2, 0, 1/2, 1, 3/2, 2, 5/2, 3, 7/2, 4, 9/2, 5}
# beta in same range
# a1 <= a2 <= a3 (WLOG by symmetry)

half_ints = [F(i, 2) for i in range(-4, 11)]  # -2, -3/2, -1, ..., 5
# For beta, use a wider range
beta_range = [F(i, 2) for i in range(-4, 11)]

print(f"Search parameters: {len(half_ints)} values for a_j, {len(beta_range)} for beta")
print(f"Total ordered triples: {len(list(combinations_with_replacement(half_ints, 3)))}")

hits = []
count = 0
total = len(list(combinations_with_replacement(half_ints, 3))) * len(beta_range)

for a_triple in combinations_with_replacement(half_ints, 3):
    a1, a2, a3 = a_triple
    for beta in beta_range:
        count += 1

        # Skip if any Pochhammer would have zero factors that break things
        # (This happens when a_j - beta is a non-negative integer <= 2n)
        skip = False
        for a in [a1, a2, a3]:
            diff = a - beta
            if diff.denominator == 1:  # integer
                d = int(diff)
                if 0 <= d <= 10:  # could hit zeros
                    skip = True
                    break
        if skip:
            continue

        try:
            S0 = raw_residue_sum(0, a1, a2, a3, beta)
            if S0 == 0:
                continue
            S1 = raw_residue_sum(1, a1, a2, a3, beta)
            S2 = raw_residue_sum(2, a1, a2, a3, beta)

            r1 = S1 / S0
            r2 = S2 / S0

            # Check if ratios match
            if r1 == target_r1 and r2 == target_r2:
                print(f"\n*** EXACT MATCH! a=({a1},{a2},{a3}), beta={beta} ***")
                print(f"  S0={S0}, S1={S1}, S2={S2}")
                hits.append((a1, a2, a3, beta, 'exact'))
            elif abs(float(r1 - target_r1)) < 1e-10 and abs(float(r2 - target_r2)) < 1e-10:
                print(f"\n*** NEAR MATCH! a=({a1},{a2},{a3}), beta={beta} ***")
                print(f"  r1 err = {float(r1 - target_r1):.2e}")
                print(f"  r2 err = {float(r2 - target_r2):.2e}")
                hits.append((a1, a2, a3, beta, 'near'))
        except (ZeroDivisionError, ValueError):
            continue

    if count % 1000 == 0:
        print(f"  progress: {count}/{total} ({100*count/total:.1f}%)", file=sys.stderr)

print(f"\nSearch complete. {count} candidates tested. {len(hits)} hits.")

if not hits:
    print("\nNo matches found in the basic half-integer box.")
    print("Trying extended search with larger parameters...")

    # Extended: add more values
    ext_ints = [F(i, 2) for i in range(-8, 20)]
    ext_beta = [F(i, 2) for i in range(-8, 20)]

    count2 = 0
    for a_triple in combinations_with_replacement(ext_ints, 3):
        a1, a2, a3 = a_triple
        for beta in ext_beta:
            count2 += 1
            if count2 > 500000:  # cap the search
                break

            skip = False
            for a in [a1, a2, a3]:
                diff = a - beta
                if diff.denominator == 1:
                    d = int(diff)
                    if 0 <= d <= 10:
                        skip = True
                        break
            if skip:
                continue

            try:
                S0 = raw_residue_sum(0, a1, a2, a3, beta)
                if S0 == 0:
                    continue
                S1 = raw_residue_sum(1, a1, a2, a3, beta)
                r1 = S1 / S0
                if abs(float(r1 - target_r1)) > 0.1:
                    continue
                S2 = raw_residue_sum(2, a1, a2, a3, beta)
                r2 = S2 / S0
                if abs(float(r1 - target_r1)) < 1e-10 and abs(float(r2 - target_r2)) < 1e-10:
                    print(f"\n*** EXTENDED MATCH! a=({a1},{a2},{a3}), beta={beta} ***")
                    print(f"  S0={S0}, S1={S1}, S2={S2}")
                    hits.append((a1, a2, a3, beta, 'extended'))
            except:
                continue

        if count2 > 500000:
            break

    print(f"\nExtended search: {count2} candidates tested. {len(hits)} total hits.")

# Also try the Zudilin parameters as a sanity check
print("\n=== Sanity check: Zudilin parameters (0,0,0; 0) ===")
S0_z = raw_residue_sum(0, F(0), F(0), F(0), F(0))
S1_z = raw_residue_sum(1, F(0), F(0), F(0), F(0))
S2_z = raw_residue_sum(2, F(0), F(0), F(0), F(0))
print(f"  S0={S0_z}, S1={S1_z}, S2={S2_z}")
print(f"  b_0=1, b_1=7, b_2=163")
print(f"  Match: {S0_z == 1 and S1_z == 7 and S2_z == 163}")
