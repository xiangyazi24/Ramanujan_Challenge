#!/usr/bin/env python3
"""P2.7: Find the minimal recurrence of W_n.

W_n = (1/256^n) Σ_{j=0}^{2n} C(2n,j)(-2)^{2n-j} T_j
where T_k is Cooper's level-11 sequence.

We know W_n does NOT satisfy order 3. Find its actual minimal order.
"""
from fractions import Fraction
import math

def compute_T(N):
    T = [Fraction(1), Fraction(4), Fraction(28)]
    for k in range(2, N):
        num = (2*(2*k+1)*(5*k**2+5*k+2)*T[k]
               - 8*k*(7*k**2+1)*T[k-1]
               + 22*k*(2*k-1)*(k-1)*T[k-2])
        T.append(num / Fraction((k+1)**3))
    return T

print("Computing T_k and W_n...", flush=True)
NMAX = 120
T = compute_T(2 * NMAX + 10)

W = []
for n in range(NMAX + 1):
    m = 2 * n
    val = Fraction(0)
    for j in range(m + 1):
        val += math.comb(m, j) * Fraction(-2)**(m - j) * T[j]
    W.append(val / Fraction(256)**n)

print(f"  Computed {len(W)} W_n values", flush=True)
print(f"  W[0]={W[0]}, W[1]={W[1]}, W[2]={W[2]}", flush=True)

# Find minimal recurrence: Σ_{j=0}^{r} P_j(n) W_{n+j} = 0
# where P_j(n) is a polynomial of degree ≤ d
# For each (r, d), set up linear system and check

for r in range(1, 10):
    for d in range(30):
        n_unknowns = (r + 1) * (d + 1)
        n_eqs = n_unknowns + 10

        if n_eqs + r >= len(W):
            break

        # Build system
        A_rows = []
        for n in range(n_eqs):
            row = []
            for j in range(r + 1):
                for k in range(d + 1):
                    row.append(Fraction(n)**k * W[n + j])
            A_rows.append(row)

        # Gaussian elimination to find null space
        m = len(A_rows)
        aug = [list(row) for row in A_rows]
        pivot_cols = []
        row_idx = 0

        for col in range(n_unknowns):
            found = -1
            for rr in range(row_idx, m):
                if aug[rr][col] != 0:
                    found = rr
                    break
            if found == -1:
                continue
            aug[row_idx], aug[found] = aug[found], aug[row_idx]
            piv = aug[row_idx][col]
            for j2 in range(n_unknowns):
                aug[row_idx][j2] /= piv
            for rr in range(m):
                if rr == row_idx:
                    continue
                if aug[rr][col] == 0:
                    continue
                factor = aug[rr][col]
                for j2 in range(n_unknowns):
                    aug[rr][j2] -= factor * aug[row_idx][j2]
            pivot_cols.append(col)
            row_idx += 1

        rank = len(pivot_cols)
        nullity = n_unknowns - rank

        if nullity > 0 and d > 0:
            # Check if the nullity is genuine (not just from low degree)
            # Also check previous degree had 0 nullity
            prev_nullity = 0
            if d > 0:
                prev_unknowns = (r + 1) * d
                prev_rank = min(prev_unknowns, rank)  # approximate
                # actually need to recompute... skip for now

            # Find null vector
            free_cols = [c for c in range(n_unknowns) if c not in pivot_cols]
            x = [Fraction(0)] * n_unknowns
            x[free_cols[0]] = Fraction(1)
            for pi in range(rank - 1, -1, -1):
                pc = pivot_cols[pi]
                val = Fraction(0)
                for j2 in range(n_unknowns):
                    if j2 != pc:
                        val += aug[pi][j2] * x[j2]
                x[pc] = -val

            # Verify on extra values
            all_ok = True
            for n in range(n_eqs, min(len(W) - r, NMAX)):
                check = Fraction(0)
                for j in range(r + 1):
                    for k in range(d + 1):
                        check += x[j * (d + 1) + k] * Fraction(n)**k * W[n + j]
                if check != 0:
                    all_ok = False
                    break

            if all_ok:
                print(f"\n*** RECURRENCE FOUND: order r={r}, coeff degree d={d} ***", flush=True)
                print(f"  Unknowns={n_unknowns}, rank={rank}, nullity={nullity}")

                # Display polynomial coefficients
                for j in range(r + 1):
                    coeffs = [x[j * (d + 1) + k] for k in range(d + 1)]
                    nonzero = [(k, c) for k, c in enumerate(coeffs) if c != 0]
                    if nonzero:
                        terms = []
                        for k, c in nonzero:
                            if k == 0:
                                terms.append(f"{c}")
                            else:
                                terms.append(f"({c})n^{k}")
                        print(f"  P_{j}(n) = {' + '.join(terms)}")

                # Try to extract integer/rational polynomial
                # Clear denominators
                all_coeffs = []
                for j in range(r + 1):
                    for k in range(d + 1):
                        if x[j*(d+1)+k] != 0:
                            all_coeffs.append(x[j*(d+1)+k])

                if all_coeffs:
                    # Find LCM of denominators
                    from math import gcd
                    lcm = 1
                    for c in all_coeffs:
                        lcm = lcm * c.denominator // gcd(lcm, c.denominator)

                    print(f"\n  Cleared denominators (×{lcm}):")
                    for j in range(r + 1):
                        coeffs = [int(x[j*(d+1)+k] * lcm) for k in range(d + 1)]
                        nonzero = [(k, c) for k, c in enumerate(coeffs) if c != 0]
                        if nonzero:
                            terms = []
                            for k, c in nonzero:
                                if k == 0:
                                    terms.append(f"{c}")
                                elif k == 1:
                                    terms.append(f"{c}·n")
                                else:
                                    terms.append(f"{c}·n^{k}")
                            print(f"    P_{j}(n) = {' + '.join(terms)}")

                break  # found minimal for this order
        elif nullity == 0:
            if d == 0 or d % 10 == 0:
                print(f"  r={r}, d={d}: full rank ({rank}), no recurrence", flush=True)

    else:
        continue
    break

print("\nDone.")
