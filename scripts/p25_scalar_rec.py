#!/usr/bin/env python3
"""P2.5: Find the CMF scalar recurrence for column 0.

Uses modular arithmetic (mod large prime) for speed, then verifies
with exact Fraction arithmetic.
"""
from fractions import Fraction

P = 2**61 - 1  # Mersenne prime

def M_entries_mod(n, p):
    m11 = (-2*n-5)*(n+3)**2*(136*n**4+1424*n**3+5548*n**2+9551*n+6141) % p
    m12 = (384*n**6+6384*n**5+44168*n**4+162698*n**3+336377*n**2+369933*n+169011) % p
    m13 = (-(480*n**4+4980*n**3+19210*n**2+32690*n+20730)) % p
    m21 = ((n+2)**2*(n+3)**2*(4*n+10)*(48*n**3+386*n**2+1017*n+879)) % p
    m22 = ((n+2)**2*(-272*n**5-3848*n**4-21732*n**3-61184*n**2-85761*n-47808)) % p
    m23 = ((n+2)**2*(320*n**3+2540*n**2+6610*n+5640)) % p
    m31 = ((-4*n-10)*(n+2)**2*(n+3)**2*(32*n**4+302*n**3+1037*n**2+1530*n+813)) % p
    m32 = ((n+2)**2*(192*n**6+2984*n**5+19116*n**4+64452*n**3+120256*n**2+117279*n+46476)) % p
    m33 = ((n+2)**2*(-16*n**5-408*n**4-2912*n**3-8884*n**2-12254*n-6240)) % p
    return [[m11,m12,m13],[m21,m22,m23],[m31,m32,m33]]

def delta_H_mod(n, p):
    return (-2*(n+2)**2*(n+3)**2*(2*n+5)*(2*n+7)**2) % p

def modinv(a, p):
    return pow(a % p, p - 2, p)

def MH_mod(n, p):
    M = M_entries_mod(n, p)
    d = delta_H_mod(n, p)
    di = modinv(d, p)
    return [[(M[i][j] * di) % p for j in range(3)] for i in range(3)]

def mat_vec_mod(M, v, p):
    return [sum(M[i][j] * v[j] for j in range(3)) % p for i in range(3)]

def row_mat_mod(row, M, p):
    return [sum(row[i] * M[i][j] for i in range(3)) % p for j in range(3)]

# Compute CMF Q̂_n (column 0) mod p
print("Computing CMF column-0 values mod p...", flush=True)
NMAX = 200

A_init = [[30921 % P, (-32972) % P, 8240 % P],
          [33750 % P, (-36000) % P, 9000 % P]]

q_mod = []
q_row = [33750 % P, (-36000) % P, 9000 % P]
for N in range(NMAX):
    q_mod.append(q_row[0])
    MH = MH_mod(N, P)
    q_row = row_mat_mod(q_row, MH, P)

print(f"  {len(q_mod)} values computed", flush=True)

# Search for recurrence: Σ_{j=0}^3 ell_j(n) · q_{n+j} ≡ 0 (mod p)
# ell_j(n) = Σ_{k=0}^d c_{j,k} · n^k

print("\n=== Searching for order-3 recurrence ===", flush=True)

for d in range(25):
    n_unknowns = 4 * (d + 1)
    n_train = n_unknowns + 10
    n_holdout = 10

    if n_train + n_holdout + 3 >= len(q_mod):
        print(f"  d={d}: need {n_train+n_holdout+3} values, have {len(q_mod)}")
        break

    # Build system mod p
    rows = []
    for n_idx in range(n_train):
        row = []
        for j in range(4):
            for k in range(d + 1):
                row.append((pow(n_idx, k, P) * q_mod[n_idx + j]) % P)
        rows.append(row)

    # Gaussian elimination mod p
    aug = [list(row) for row in rows]
    m = n_train
    pivot_cols = []
    row_idx = 0

    for col in range(n_unknowns):
        found = -1
        for rr in range(row_idx, m):
            if aug[rr][col] % P != 0:
                found = rr
                break
        if found == -1:
            continue
        aug[row_idx], aug[found] = aug[found], aug[row_idx]
        piv_inv = modinv(aug[row_idx][col], P)
        for j2 in range(n_unknowns):
            aug[row_idx][j2] = (aug[row_idx][j2] * piv_inv) % P
        for rr in range(m):
            if rr == row_idx:
                continue
            if aug[rr][col] % P == 0:
                continue
            factor = aug[rr][col]
            for j2 in range(n_unknowns):
                aug[rr][j2] = (aug[rr][j2] - factor * aug[row_idx][j2]) % P
        pivot_cols.append(col)
        row_idx += 1

    rank = len(pivot_cols)
    nullity = n_unknowns - rank

    if nullity > 0:
        # Find null vector
        free_cols = [c for c in range(n_unknowns) if c not in pivot_cols]
        x = [0] * n_unknowns
        x[free_cols[0]] = 1
        for pi in range(rank - 1, -1, -1):
            pc = pivot_cols[pi]
            val = 0
            for j2 in range(n_unknowns):
                if j2 != pc:
                    val = (val + aug[pi][j2] * x[j2]) % P
            x[pc] = (-val) % P

        # Verify on holdout
        all_ok = True
        for n_idx in range(n_train, n_train + n_holdout):
            check = 0
            for j in range(4):
                for k in range(d + 1):
                    check = (check + x[j*(d+1)+k] * pow(n_idx, k, P) * q_mod[n_idx + j]) % P
            if check % P != 0:
                all_ok = False
                break

        if all_ok:
            print(f"\n  *** RECURRENCE FOUND: degree d={d} ***", flush=True)
            print(f"  Unknowns={n_unknowns}, rank={rank}, nullity={nullity}")

            # Now verify with exact arithmetic
            print("  Verifying with exact Fraction arithmetic...", flush=True)

            # Compute exact q values
            q_exact = []
            q_row_ex = [Fraction(33750), Fraction(-36000), Fraction(9000)]
            for N in range(min(NMAX, d * 3 + 30)):
                q_exact.append(q_row_ex[0])
                if N < NMAX - 1:
                    def M_entries_frac(n):
                        n = Fraction(n)
                        m11 = (-2*n-5)*(n+3)**2*(136*n**4+1424*n**3+5548*n**2+9551*n+6141)
                        m12 = 384*n**6+6384*n**5+44168*n**4+162698*n**3+336377*n**2+369933*n+169011
                        m13 = -(480*n**4+4980*n**3+19210*n**2+32690*n+20730)
                        m21 = (n+2)**2*(n+3)**2*(4*n+10)*(48*n**3+386*n**2+1017*n+879)
                        m22 = (n+2)**2*(-272*n**5-3848*n**4-21732*n**3-61184*n**2-85761*n-47808)
                        m23 = (n+2)**2*(320*n**3+2540*n**2+6610*n+5640)
                        m31 = (-4*n-10)*(n+2)**2*(n+3)**2*(32*n**4+302*n**3+1037*n**2+1530*n+813)
                        m32 = (n+2)**2*(192*n**6+2984*n**5+19116*n**4+64452*n**3+120256*n**2+117279*n+46476)
                        m33 = (n+2)**2*(-16*n**5-408*n**4-2912*n**3-8884*n**2-12254*n-6240)
                        return [[m11,m12,m13],[m21,m22,m23],[m31,m32,m33]]

                    def delta_H_frac(n):
                        n = Fraction(n)
                        return Fraction(-2)*(n+2)**2*(n+3)**2*(2*n+5)*(2*n+7)**2

                    MF = M_entries_frac(N)
                    df = delta_H_frac(N)
                    MH = [[MF[i][j] / df for j in range(3)] for i in range(3)]
                    q_row_ex = [sum(q_row_ex[i] * MH[i][j] for i in range(3)) for j in range(3)]

            # The mod-p solution x needs to be lifted to rational.
            # But we can just re-solve exactly with the known degree.
            # Actually, let's just verify a few: substitute the mod-p coefficients
            # into the exact recurrence and check.

            # Better: re-solve exactly with known d
            print(f"  Re-solving with exact arithmetic at degree {d}...", flush=True)
            n_exact = n_unknowns + 5
            A_rows = []
            for n_idx in range(n_exact):
                row = []
                for j in range(4):
                    for k in range(d + 1):
                        row.append(Fraction(n_idx)**k * q_exact[n_idx + j])
                A_rows.append(row)

            # Gaussian elimination
            aug_ex = [list(row) for row in A_rows]
            pivot_cols_ex = []
            row_idx = 0
            for col in range(n_unknowns):
                found = -1
                for rr in range(row_idx, n_exact):
                    if aug_ex[rr][col] != 0:
                        found = rr
                        break
                if found == -1:
                    continue
                aug_ex[row_idx], aug_ex[found] = aug_ex[found], aug_ex[row_idx]
                piv = aug_ex[row_idx][col]
                for j2 in range(n_unknowns):
                    aug_ex[row_idx][j2] /= piv
                for rr in range(n_exact):
                    if rr == row_idx:
                        continue
                    if aug_ex[rr][col] == 0:
                        continue
                    factor = aug_ex[rr][col]
                    for j2 in range(n_unknowns):
                        aug_ex[rr][j2] -= factor * aug_ex[row_idx][j2]
                pivot_cols_ex.append(col)
                row_idx += 1

            rank_ex = len(pivot_cols_ex)
            nullity_ex = n_unknowns - rank_ex
            print(f"  Exact: rank={rank_ex}, nullity={nullity_ex}")

            if nullity_ex > 0:
                free_cols_ex = [c for c in range(n_unknowns) if c not in pivot_cols_ex]
                x_ex = [Fraction(0)] * n_unknowns
                x_ex[free_cols_ex[0]] = Fraction(1)
                for pi in range(rank_ex - 1, -1, -1):
                    pc = pivot_cols_ex[pi]
                    val = Fraction(0)
                    for j2 in range(n_unknowns):
                        if j2 != pc:
                            val += aug_ex[pi][j2] * x_ex[j2]
                    x_ex[pc] = -val

                # Verify on extra values
                ok = True
                for n_idx in range(n_exact, len(q_exact) - 3):
                    check = Fraction(0)
                    for j in range(4):
                        for k in range(d + 1):
                            check += x_ex[j*(d+1)+k] * Fraction(n_idx)**k * q_exact[n_idx + j]
                    if check != 0:
                        ok = False
                        print(f"  FAIL at n={n_idx}")
                        break

                if ok:
                    print(f"  VERIFIED exactly!")

                    # Clear denominators and display
                    from math import gcd
                    all_c = [c for c in x_ex if c != 0]
                    lcm_d = 1
                    for c in all_c:
                        lcm_d = lcm_d * c.denominator // gcd(lcm_d, c.denominator)

                    print(f"\n  Recurrence coefficients (×{lcm_d}):")
                    for j in range(4):
                        coeffs = [int(x_ex[j*(d+1)+k] * lcm_d) for k in range(d+1)]
                        nonzero = [(k, c) for k, c in enumerate(coeffs) if c != 0]
                        if nonzero:
                            terms = []
                            for k, c in nonzero:
                                if k == 0:
                                    terms.append(str(c))
                                elif k == 1:
                                    terms.append(f"{c}*n")
                                else:
                                    terms.append(f"{c}*n^{k}")
                            print(f"    ell_{j}(n) = {' + '.join(terms)}")
                        else:
                            print(f"    ell_{j}(n) = 0")

                    # Store for later use
                    print(f"\n  Degree of ell_j:")
                    for j in range(4):
                        coeffs = [x_ex[j*(d+1)+k] for k in range(d+1)]
                        deg_j = max((k for k, c in enumerate(coeffs) if c != 0), default=-1)
                        print(f"    ell_{j}: degree {deg_j}")

            break

    if d % 5 == 0:
        print(f"  d={d}: rank={rank}, nullity={nullity}", flush=True)

print("\nDone.")
