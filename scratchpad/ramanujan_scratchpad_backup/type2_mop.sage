"""
Test whether the type-II Hermite-Padé polynomials for measures
  mu_2 = -log(t) dt  (moments 1/(k+1)^2)
  mu_3 = (1/2) log^2(t) dt  (moments 1/(k+1)^3)
on (0,1) with step-line indices give the P2.7 recurrence.

The type-II polynomial Q_n(t) of degree 2n on the step line
(n,n) satisfies:
  int_0^1 t^j Q_n(t) (-log t) dt = 0  for j = 0,...,n-1
  int_0^1 t^j Q_n(t) (1/2 log^2 t) dt = 0  for j = 0,...,n-1

Equivalently, with moments:
  m2(k) = 1/(k+1)^2,  m3(k) = 1/(k+1)^3
The orthogonality conditions for monic Q_n(t) = t^(2n) + ... are:
  sum_{l=0}^{2n} q_l * m2(j+l) = 0  for j = 0,...,n-1
  sum_{l=0}^{2n} q_l * m3(j+l) = 0  for j = 0,...,n-1
"""
from sage.all import *

Q = QQ

def m2(k):
    return Q(1) / Q(k+1)^2

def m3(k):
    return Q(1) / Q(k+1)^3

def type2_mop(n):
    """Compute the monic type-II polynomial Q_n(t) of degree 2n."""
    if n == 0:
        return [Q(1)]  # Q_0 = 1

    deg = 2*n
    # Monic: coefficient of t^{2n} is 1
    # Unknowns: q_0, q_1, ..., q_{2n-1} (coefficients of t^0,...,t^{2n-1})
    # q_{2n} = 1

    num_unknowns = deg  # q_0, ..., q_{2n-1}
    num_equations = 2*n  # n from mu_2, n from mu_3

    # Build the linear system
    rows = []
    rhs = []

    # Orthogonality for mu_2
    for j in range(n):
        row = [m2(j+l) for l in range(deg)]
        rows.append(row)
        rhs.append(-m2(j + deg))  # from the monic term

    # Orthogonality for mu_3
    for j in range(n):
        row = [m3(j+l) for l in range(deg)]
        rows.append(row)
        rhs.append(-m3(j + deg))

    M = matrix(Q, rows)
    b = vector(Q, rhs)

    try:
        x = M.solve_right(b)
    except ValueError:
        print(f"  n={n}: system is singular!")
        return None

    coeffs = list(x) + [Q(1)]  # q_0, ..., q_{2n-1}, 1
    return coeffs

# Compute type-II MOPs for small n
print("Computing type-II MOPs for n=0,...,6")
Rt = PolynomialRing(Q, 't')
t_var = Rt.gen()

for n in range(7):
    coeffs = type2_mop(n)
    if coeffs is None:
        continue
    Q_poly = sum(c * t_var^i for i, c in enumerate(coeffs))
    val_at_1 = Q_poly(1)
    print(f"n={n}: Q_n(1) = {val_at_1}")
    if n <= 3:
        print(f"  Q_n(t) = {Q_poly}")

# P2.7 initial values for q_n
q_p27 = [
    Q(-215040420000),
    Q(-167282265043404) / Q(905),
    Q(-964185327658080) / Q(6071),
]

print(f"\nP2.7 initial values:")
print(f"  q_0 = {q_p27[0]}")
print(f"  q_1 = {q_p27[1]}")
print(f"  q_2 = {q_p27[2]}")

# Check if Q_n(1) matches P2.7 q_n (up to a common scaling)
print(f"\nChecking if Q_n(1) is proportional to P2.7 q_n:")
Q0_at_1 = type2_mop(0)
Q1_at_1 = sum(c for c in type2_mop(1))
Q2_at_1 = sum(c for c in type2_mop(2))

vals = [Q0_at_1, Q1_at_1, Q2_at_1]
if isinstance(vals[0], list):
    vals = [sum(c for c in v) for v in [type2_mop(i) for i in range(3)]]

print(f"  MOP Q_0(1) = {vals[0]}")
print(f"  MOP Q_1(1) = {vals[1]}")
print(f"  MOP Q_2(1) = {vals[2]}")

if vals[0] != 0 and q_p27[0] != 0:
    ratio_0 = q_p27[0] / vals[0]
    print(f"\n  ratio q_0/Q_0(1) = {ratio_0}")
    if vals[1] != 0:
        ratio_1 = q_p27[1] / vals[1]
        print(f"  ratio q_1/Q_1(1) = {ratio_1}")
        if ratio_0 == ratio_1:
            print("  *** RATIOS MATCH! Same proportionality constant. ***")
        else:
            print(f"  Ratios differ: {float(ratio_0):.6e} vs {float(ratio_1):.6e}")

# Also try the alternating step line: (1,0), (1,1), (2,1), (2,2), ...
print("\n\n=== Alternating step line ===")
def type2_mop_step(r, s):
    """Type-II MOP for multi-index (r, s), degree r+s."""
    deg = r + s
    if deg == 0:
        return [Q(1)]

    num_unknowns = deg
    rows = []
    rhs = []

    for j in range(r):
        row = [m2(j+l) for l in range(deg)]
        rows.append(row)
        rhs.append(-m2(j + deg))

    for j in range(s):
        row = [m3(j+l) for l in range(deg)]
        rows.append(row)
        rhs.append(-m3(j + deg))

    M = matrix(Q, rows)
    b = vector(Q, rhs)

    try:
        x = M.solve_right(b)
    except ValueError:
        print(f"  (r,s)=({r},{s}): system is singular!")
        return None

    return list(x) + [Q(1)]

# Step line: (0,0), (1,0), (1,1), (2,1), (2,2), (3,2), (3,3), ...
step_indices = [(0,0), (1,0), (1,1), (2,1), (2,2), (3,2), (3,3),
                (4,3), (4,4), (5,4), (5,5), (6,5), (6,6)]

print("Step line evaluations:")
step_vals = []
for r, s in step_indices:
    coeffs = type2_mop_step(r, s)
    if coeffs is None:
        step_vals.append(None)
        continue
    val = sum(coeffs)  # Q(1)
    step_vals.append(val)
    print(f"  ({r},{s}): Q(1) = {val}")

# Check 4-term recurrence on step line values
print("\n=== Checking for 4-term recurrence on diagonal Q_n(1) ===")
diag_vals = []
for n in range(7):
    coeffs = type2_mop(n)
    if coeffs is None:
        break
    diag_vals.append(sum(coeffs))

print("Diagonal Q_n(1):", diag_vals)

# Try: a(n)*v(n+1) + b(n)*v(n) + c(n)*v(n-1) + d(n)*v(n-2) = 0
# for n = 2, 3, 4 with polynomial a, b, c, d
if len(diag_vals) >= 6:
    # Just check the recurrence residual with P2.7 coefficients
    def AA(n):
        return Q(1024*(2*n+5)^4*(2*n+7)^3*(2*n+9)^3
                *(946*n^2 + 6407*n + 10860))
    def BB(n):
        return Q(128*(2*n+7)^3*(2*n+9)^3
                *(104060*n^6 + 1745370*n^5 + 12145238*n^4
                  + 44886481*n^3 + 92943995*n^2
                  + 102256019*n + 46709052))
    def CC(n):
        return Q(16*(n+3)^4*(2*n+9)^3
                *(3784*n^5 + 57792*n^4 + 351019*n^3
                  + 1059230*n^2 + 1587211*n + 944620))
    def DD(n):
        return Q((n+3)^4*(n+4)^6
                *(946*n^2 + 4515*n + 5399))

    print("\nChecking P2.7 recurrence on Q_n(1):")
    v = diag_vals
    for n in range(2, min(len(v)-1, 6)):
        # Monic: v[n+1] = B(n)/A(n)*v[n] - C(n-1)/A(n-1)*v[n-1] + D(n-2)/A(n-2)*v[n-2]
        predicted = BB(n)/AA(n)*v[n] - CC(n-1)/AA(n-1)*v[n-1] + DD(n-2)/AA(n-2)*v[n-2]
        residual = v[n+1] - predicted
        if residual == 0:
            print(f"  n={n}: P2.7 recurrence HOLDS")
        else:
            print(f"  n={n}: residual = {float(residual):.6e} (FAILS)")

print("\nDone.")
