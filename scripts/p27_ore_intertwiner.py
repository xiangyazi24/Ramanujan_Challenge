#!/usr/bin/env python3
"""P2.7: Search for the Ore intertwiner q_n = Σ P_j(n)·W_{n+j}.

We have exact values of q_n (from the P2.7 recurrence) and W_n
(from the level-11 binomial transform). Search for polynomial
coefficients P_0(n), ..., P_r(n) such that q_n = Σ_j P_j(n)·W_{n+j}.
"""
from fractions import Fraction
import math

# P2.7 recurrence coefficients (exact)
def A27(n):
    return 1024*(2*n+5)**4*(2*n+7)**3*(2*n+9)**3*(946*n**2+6407*n+10860)

def B27(n):
    P6 = 104060*n**6+1745370*n**5+12145238*n**4+44886481*n**3+92943995*n**2+102256019*n+46709052
    return 128*(2*n+7)**3*(2*n+9)**3*P6

def C27(n):
    P5 = 3784*n**5+57792*n**4+351019*n**3+1059230*n**2+1587211*n+944620
    return 16*(n+3)**4*(2*n+9)**3*P5

def D27(n):
    return (n+3)**4*(n+4)**6*(946*n**2+4515*n+5399)

# Compute q_n values
print("Computing q_n from P2.7 recurrence...", flush=True)
q = [Fraction(-215040420000),
     Fraction(-167282265043404, 905),
     Fraction(-964185327658080, 6071)]

NMAX = 80
for n in range(2, NMAX):
    q_next = (Fraction(B27(n)) * q[n] / Fraction(A27(n))
              - Fraction(C27(n-1)) * q[n-1] / Fraction(A27(n-1))
              + Fraction(D27(n-2)) * q[n-2] / Fraction(A27(n-2)))
    q.append(q_next)

print(f"  Computed {len(q)} q_n values", flush=True)

# Compute W_n = b_{2n}/256^n from Cooper's level-11
def compute_T(N):
    T = [Fraction(1), Fraction(4), Fraction(28)]
    for k in range(2, N):
        num = (2*(2*k+1)*(5*k**2+5*k+2)*T[k]
               - 8*k*(7*k**2+1)*T[k-1]
               + 22*k*(2*k-1)*(k-1)*T[k-2])
        T.append(num / Fraction((k+1)**3))
    return T

print("Computing W_n from level-11 binomial transform...", flush=True)
T = compute_T(200)

W = []
for n in range(NMAX + 10):
    m = 2 * n
    val = Fraction(0)
    for j in range(min(m + 1, len(T))):
        val += math.comb(m, j) * Fraction(-2)**(m - j) * T[j]
    W.append(val / Fraction(256)**n)

print(f"  Computed {len(W)} W_n values", flush=True)

# Verify a few
print(f"  q[0] = {float(q[0]):.6e}, W[0] = {float(W[0]):.6e}")
print(f"  q[1] = {float(q[1]):.6e}, W[1] = {float(W[1]):.6e}")

# Search for: q_n = Σ_{j=0}^{r} P_j(n) · W_{n+j}
# where P_j(n) = Σ_{k=0}^{d} c_{j,k} · n^k
# Unknowns: (r+1)*(d+1) coefficients c_{j,k}
# Equations: one per n value

print("\n=== Ore intertwiner search: q_n = Σ P_j(n)·W_{n+j} ===", flush=True)

for r in range(1, 7):  # order r = number of W terms minus 1
    for d in range(20):  # degree of P_j
        n_unknowns = (r + 1) * (d + 1)
        n_train = n_unknowns + 5
        n_holdout = 5

        if n_train + n_holdout + r >= len(q) or n_train + n_holdout + r >= len(W):
            break

        # Build system: for each n, q[n] = Σ_j Σ_k c_{j,k} · n^k · W[n+j]
        A_rows = []
        b_vec = []
        for n in range(n_train + n_holdout):
            row = []
            for j in range(r + 1):
                for k in range(d + 1):
                    row.append(Fraction(n)**k * W[n + j])
            A_rows.append(row)
            b_vec.append(q[n])

        # Solve train system
        m = n_train
        aug = [list(A_rows[i]) + [b_vec[i]] for i in range(m)]
        n_cols = n_unknowns

        pivot_cols = []
        row_idx = 0
        for col in range(n_cols):
            found = -1
            for rr in range(row_idx, m):
                if aug[rr][col] != 0:
                    found = rr
                    break
            if found == -1:
                continue
            aug[row_idx], aug[found] = aug[found], aug[row_idx]
            piv = aug[row_idx][col]
            for j2 in range(n_cols + 1):
                aug[row_idx][j2] /= piv
            for rr in range(m):
                if rr == row_idx:
                    continue
                if aug[rr][col] == 0:
                    continue
                factor = aug[rr][col]
                for j2 in range(n_cols + 1):
                    aug[rr][j2] -= factor * aug[row_idx][j2]
            pivot_cols.append(col)
            row_idx += 1

        rank = len(pivot_cols)
        if rank < n_unknowns:
            continue

        # Extract solution
        x = [Fraction(0)] * n_unknowns
        for pi, pc in enumerate(pivot_cols):
            x[pc] = aug[pi][n_cols]

        # Verify on holdout
        max_err = Fraction(0)
        all_exact = True
        for n in range(n_train, n_train + n_holdout):
            pred = Fraction(0)
            for j in range(r + 1):
                for k in range(d + 1):
                    pred += x[j * (d + 1) + k] * Fraction(n)**k * W[n + j]
            if pred != q[n]:
                all_exact = False
                if q[n] != 0:
                    err = abs(pred - q[n])
                    if err > max_err:
                        max_err = err
            # else: exact match

        if all_exact:
            print(f"\n  *** EXACT MATCH: order r={r}, degree d={d} ***", flush=True)
            print(f"  Unknowns: {n_unknowns}, train: {n_train}, holdout: {n_holdout}")

            # Display coefficients
            for j in range(r + 1):
                coeffs = [x[j * (d + 1) + k] for k in range(d + 1)]
                nonzero = [(k, c) for k, c in enumerate(coeffs) if c != 0]
                if nonzero:
                    terms = []
                    for k, c in nonzero:
                        if k == 0:
                            terms.append(f"{c}")
                        else:
                            terms.append(f"({c})·n^{k}")
                    print(f"  P_{j}(n) = {' + '.join(terms)}")
                else:
                    print(f"  P_{j}(n) = 0")

            # Extra verification at larger n
            extra_ok = True
            for n in range(n_train + n_holdout, min(len(q) - r, len(W) - r)):
                pred = Fraction(0)
                for j in range(r + 1):
                    for k in range(d + 1):
                        pred += x[j * (d + 1) + k] * Fraction(n)**k * W[n + j]
                if pred != q[n]:
                    extra_ok = False
                    print(f"  FAIL at n={n}")
                    break

            if extra_ok:
                print(f"  Verified for ALL n=0..{min(len(q)-r, len(W)-r)-1} ✓")
            else:
                print(f"  Some extra values fail")

            # Check if we can find the pattern
            break
        else:
            if d <= 3 or d % 5 == 0:
                print(f"  r={r}, d={d}: train rank={rank}, holdout NOT exact", flush=True)

    else:
        continue
    break

print("\nDone.")
