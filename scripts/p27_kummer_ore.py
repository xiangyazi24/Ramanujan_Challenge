#!/usr/bin/env python3
"""P2.7: Kummer-twisted Ore intertwiner search.

From Q4887: polynomial Ore is impossible (exponent gap 3/2).
The correct form is q_n = h_n · Σ u_j(n) W_{n+j}
where h_{n+1}/h_n = 1 + 3/(2n) + O(1/n²).

Try several natural twists h_n.
"""
from fractions import Fraction
import math

# P2.7 recurrence coefficients
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

print("Computing q_n...", flush=True)
q = [Fraction(-215040420000),
     Fraction(-167282265043404, 905),
     Fraction(-964185327658080, 6071)]

NMAX = 80
for n in range(2, NMAX):
    q_next = (Fraction(B27(n)) * q[n] / Fraction(A27(n))
              - Fraction(C27(n-1)) * q[n-1] / Fraction(A27(n-1))
              + Fraction(D27(n-2)) * q[n-2] / Fraction(A27(n-2)))
    q.append(q_next)
print(f"  {len(q)} values", flush=True)

# Cooper level-11
print("Computing W_n...", flush=True)
T = [Fraction(1), Fraction(4), Fraction(28)]
for k in range(2, 200):
    num = (2*(2*k+1)*(5*k**2+5*k+2)*T[k]
           - 8*k*(7*k**2+1)*T[k-1]
           + 22*k*(2*k-1)*(k-1)*T[k-2])
    T.append(num / Fraction((k+1)**3))

W = []
for n in range(NMAX + 10):
    m = 2 * n
    val = Fraction(0)
    for j in range(m + 1):
        val += math.comb(m, j) * Fraction(-2)**(m - j) * T[j]
    W.append(val / Fraction(256)**n)
print(f"  {len(W)} values", flush=True)

def rising_frac(a, k):
    result = Fraction(1)
    for i in range(k):
        result *= (a + i)
    return result

# Define twist candidates
twists = {}

# Twist 1: h_n = (5/2)_n / n!
h1 = [rising_frac(Fraction(5,2), n) / Fraction(math.factorial(n)) for n in range(NMAX+1)]
twists['(5/2)_n/n!'] = h1

# Twist 2: h_n = (3/2)_n / n!
h2 = [rising_frac(Fraction(3,2), n) / Fraction(math.factorial(n)) for n in range(NMAX+1)]
twists['(3/2)_n/n!'] = h2

# Twist 3: h_n = C(2n,n) / 4^n  (has ratio (2n+1)/(2n+2) = 1 + 1/(2n) - ... not 3/(2n))
# Actually this has the wrong leading term. Skip.

# Twist 4: h_n = (2n+3)!! / (2^n * n!)
# (2n+3)!! = 1·3·5·...·(2n+3)
h4 = [Fraction(1)]
for n in range(1, NMAX+1):
    h4.append(h4[-1] * Fraction(2*n+3, 2*n))
twists['(2n+3)!!/(2^n n!)'] = h4

# Twist 5: h_n = Γ(n+5/2) / (Γ(5/2)·Γ(n+1)) = (5/2)_n / n! (same as twist 1)
# Skip duplicate

# Twist 6: h_n = (3/2)_n (5/2)_n / (n!)^2
h6 = [rising_frac(Fraction(3,2), n) * rising_frac(Fraction(5,2), n)
      / Fraction(math.factorial(n))**2 for n in range(NMAX+1)]
twists['(3/2)(5/2)/n!^2'] = h6

# Twist 7: h_n = (7/2)_n / n!
h7 = [rising_frac(Fraction(7,2), n) / Fraction(math.factorial(n)) for n in range(NMAX+1)]
twists['(7/2)_n/n!'] = h7

# Twist 8: h_n = C(2n+2,n+1) / 4^(n+1) ≈ (n+1)^{-1/2} / √π
# ratio = (2n+4)(2n+3)/((n+2)·4) = (2n+3)/2 · (2n+4)/(4(n+2)) = (2n+3)(n+2)/(2(n+2)) = ...
# Not the right asymptotics. Skip.

# Twist 9: h_n = (2n+5)!! / (2^(n+2) * n!)
# Try this variant
h9 = [Fraction(1)]
for n in range(1, NMAX+1):
    h9.append(h9[-1] * Fraction(2*n+5, 2*n+2))
twists['(2n+5)!!_shift'] = h9

# Check ratios to verify asymptotics
print("\n=== Twist ratio h_{n+1}/h_n at n=100 ===", flush=True)
for name, h in twists.items():
    if len(h) > 101 and h[100] != 0:
        ratio = h[101] / h[100]
        approx = float(ratio)
        expected = 1 + 1.5/100  # 1 + 3/(2n)
        print(f"  {name}: ratio = {approx:.8f}, expected ~{expected:.8f}, diff = {approx-expected:.2e}")

def search_ore(name, h, max_r=6, max_d=20):
    """Search for q_n/h_n = Σ P_j(n) W_{n+j}"""
    print(f"\n=== Twist: {name} ===", flush=True)

    # Compute twisted sequence
    qh = []
    for n in range(min(len(q), len(h))):
        if h[n] == 0:
            return
        qh.append(q[n] / h[n])

    # Check ratio qh[n]/W[n] behavior
    for n in [0, 1, 5, 10, 20]:
        if n < len(qh) and n < len(W) and W[n] != 0:
            r = float(qh[n] / W[n])
            print(f"  qh[{n}]/W[{n}] = {r:.6e}", flush=True)

    for r in range(1, max_r + 1):
        for d in range(max_d + 1):
            n_unknowns = (r + 1) * (d + 1)
            n_train = n_unknowns + 5
            n_holdout = 5

            if n_train + n_holdout + r >= len(qh) or n_train + n_holdout + r >= len(W):
                break

            # Build system
            A_rows = []
            b_vec = []
            for n in range(n_train + n_holdout):
                row = []
                for j in range(r + 1):
                    for k in range(d + 1):
                        row.append(Fraction(n)**k * W[n + j])
                A_rows.append(row)
                b_vec.append(qh[n])

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

            x = [Fraction(0)] * n_unknowns
            for pi, pc in enumerate(pivot_cols):
                x[pc] = aug[pi][n_cols]

            # Verify holdout
            all_exact = True
            for n in range(n_train, n_train + n_holdout):
                pred = Fraction(0)
                for j in range(r + 1):
                    for k in range(d + 1):
                        pred += x[j * (d + 1) + k] * Fraction(n)**k * W[n + j]
                if pred != qh[n]:
                    all_exact = False
                    break

            if all_exact:
                print(f"\n  *** EXACT MATCH: r={r}, d={d} ***", flush=True)
                for j in range(r + 1):
                    coeffs = [x[j*(d+1)+k] for k in range(d+1)]
                    nonzero = [(k,c) for k,c in enumerate(coeffs) if c != 0]
                    if nonzero:
                        terms = [f"({c})·n^{k}" if k > 0 else f"{c}" for k, c in nonzero]
                        print(f"    u_{j}(n) = {' + '.join(terms)}")

                # Full verification
                ok = True
                for n in range(n_train + n_holdout, min(len(qh) - r, len(W) - r)):
                    pred = Fraction(0)
                    for j in range(r + 1):
                        for k in range(d + 1):
                            pred += x[j*(d+1)+k] * Fraction(n)**k * W[n + j]
                    if pred != qh[n]:
                        ok = False
                        print(f"    FAIL at n={n}")
                        break
                if ok:
                    print(f"    Verified ALL n=0..{min(len(qh)-r, len(W)-r)-1}")
                return True

            if d == 0 or d == 5 or d == 10 or d == 15:
                print(f"  r={r}, d={d}: no match", flush=True)

    print(f"  No match found", flush=True)
    return False

for name, h in twists.items():
    if search_ore(name, h):
        break

print("\nDone.")
