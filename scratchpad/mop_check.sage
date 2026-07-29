#!/usr/bin/env sage
"""
Check whether the step-line type-II Multiple Orthogonal Polynomials (MOP)
for measures dμ₂ = -log(t) dt and dμ₃ = ½ log²(t) dt on [0,1]
produce the P2.7 sequences.

Moments:
  m₂(k) = ∫₀¹ t^k (-log t) dt = 1/(k+1)²
  m₃(k) = ∫₀¹ t^k (½ log²t) dt = 1/(k+1)³

Step-line multi-index: (n₁, n₂) alternating
  n=0: (0,0), n=1: (1,0), n=2: (1,1), n=3: (2,1), n=4: (2,2), ...

Type-II MOP Q_n(x) of degree n₁+n₂ satisfying:
  ∫₀¹ Q_n(t) t^j (-log t) dt = 0  for j = 0,...,n₁-1
  ∫₀¹ Q_n(t) t^j (½ log²t) dt = 0 for j = 0,...,n₂-1

Then:
  q_n = Q_n(1)  (the denominator)
  p_n = Σ c_{n,j} (H_j^{(2)} + H_j^{(3)})  (the numerator)
  e_n = p_n - (ζ(2)+ζ(3)) q_n = -∫∫ K(x,y) Q_n(xy) dx dy

If these satisfy P2.7, then c₀(e)=0 is proved!
"""
import mpmath
mpmath.mp.dps = 50

# Moments
def m2(k):
    return QQ(1) / QQ(k+1)**2

def m3(k):
    return QQ(1) / QQ(k+1)**3

# Step-line multi-index
def step_line_index(n):
    k = n // 2
    if n % 2 == 0:
        return (k, k)
    else:
        return (k+1, k)

# Compute type-II MOP Q_n as a polynomial
# Q_n(x) = Σ_{j=0}^{deg} c_j x^j, monic (c_{deg} = 1)
# Orthogonality: ∫ Q(t) t^l w_s(t) dt = 0 for appropriate l, s

def compute_mop(n):
    n1, n2 = step_line_index(n)
    deg = n1 + n2

    if deg == 0:
        return [QQ(1)]  # Q_0 = 1

    # Set up the linear system for c_0, ..., c_{deg-1} (c_{deg} = 1 monic)
    # Conditions from μ₂: ∫ Q(t) t^l (-log t) dt = 0 for l = 0,...,n1-1
    # This gives: Σ_j c_j m2(j+l) + m2(deg+l) = 0 for l = 0,...,n1-1

    # Conditions from μ₃: ∫ Q(t) t^l (½ log²t) dt = 0 for l = 0,...,n2-1
    # This gives: Σ_j c_j m3(j+l) + m3(deg+l) = 0 for l = 0,...,n2-1

    rows = []
    rhs = []

    # μ₂ conditions
    for l in range(n1):
        row = [m2(j+l) for j in range(deg)]
        rows.append(row)
        rhs.append(-m2(deg+l))

    # μ₃ conditions
    for l in range(n2):
        row = [m3(j+l) for j in range(deg)]
        rows.append(row)
        rhs.append(-m3(deg+l))

    M = matrix(QQ, rows)
    b = vector(QQ, rhs)

    try:
        c = M.solve_right(b)
        coeffs = list(c) + [QQ(1)]
        return coeffs
    except:
        print("  Singular system at n=%d" % n)
        return None

# Compute MOP for several n values
print("=== Step-line Type-II MOP ===")
N_max = 20

q_mop = []  # Q_n(1)
p_mop = []  # Σ c_j (H_j^{(2)} + H_j^{(3)})

def H(s, k):
    """Harmonic number H_k^{(s)} = Σ_{j=1}^k 1/j^s"""
    return sum(QQ(1)/QQ(j)**s for j in range(1, k+1))

for n in range(N_max):
    coeffs = compute_mop(n)
    if coeffs is None:
        q_mop.append(None)
        p_mop.append(None)
        continue

    deg = len(coeffs) - 1
    n1, n2 = step_line_index(n)

    # q = Q_n(1) = Σ c_j
    q = sum(coeffs)
    q_mop.append(q)

    # p = Σ c_j (H_j^{(2)} + H_j^{(3)})
    p = sum(coeffs[j] * (H(2,j) + H(3,j)) for j in range(len(coeffs)))
    p_mop.append(p)

    # Error e = p - L*q
    L = float(mpmath.zeta(2) + mpmath.zeta(3))
    e = float(p) - L * float(q)

    print("n=%2d (%d,%d): deg=%d, q=%.6g, p=%.6g, e=%.3e" %
          (n, n1, n2, deg, float(q), float(p), e))

# Check ratios q_{n+1}/q_n for Poincaré root
print("\n=== q ratio (Poincaré check) ===")
for n in range(1, min(15, len(q_mop))):
    if q_mop[n] and q_mop[n-1] and q_mop[n-1] != 0:
        print("n=%2d: q_{n}/q_{n-1} = %.6f" % (n, float(q_mop[n]/q_mop[n-1])))

# Check if the MOP sequences satisfy a 4-term recurrence
print("\n=== Guessing recurrence for q_mop ===")
from ore_algebra import OreAlgebra
Rn = PolynomialRing(QQ, 'n')
n_var = Rn.gen()
OS = OreAlgebra(Rn, 'Sn')
Sn = OS.gen()

try:
    from ore_algebra import guess
    q_list = [q_mop[i] for i in range(len(q_mop)) if q_mop[i] is not None]
    rec_q = guess(q_list, OS, order=3)
    print("Found order-3 recurrence for q_mop!")
    print("rec_q = %s" % rec_q)

    # Extract leading coefficient degree
    for j in range(4):
        print("  P_%d: degree %d" % (j, rec_q[j].degree()))
except Exception as e:
    print("Guessing failed: %s" % e)

# Also try p_mop
print("\n=== Guessing recurrence for p_mop ===")
try:
    p_list = [p_mop[i] for i in range(len(p_mop)) if p_mop[i] is not None]
    rec_p = guess(p_list, OS, order=3)
    print("Found order-3 recurrence for p_mop!")
    print("rec_p = %s" % rec_p)
except Exception as e:
    print("Guessing failed: %s" % e)

# Compare q_mop to P2.7 initial conditions
print("\n=== Compare to P2.7 ===")
# P2.7 q values
q_p27_0 = QQ(-215040420000)
q_p27_1 = QQ(-167282265043404) / QQ(905)
q_p27_2 = QQ(-964185327658080) / QQ(6071)

# Zudilin a values
a_vals = []
for i in range(20):
    s = QQ(0)
    for k in range(i+1):
        s += binomial(i,k)^2 * binomial(i+k,i) * binomial(i+2*k,i)
    a_vals.append(s)

print("Zudilin a_n:", a_vals[:5])
print("MOP q_n:", [q_mop[i] for i in range(5)])

# Check gauge: q_p27 / q_mop ?
if q_mop[0] and q_mop[0] != 0:
    g0 = q_p27_0 / q_mop[0]
    print("\nq_p27_0 / q_mop_0 = %s" % g0)
if q_mop[1] and q_mop[1] != 0:
    g1 = q_p27_1 / q_mop[1]
    print("q_p27_1 / q_mop_1 = %s" % g1)
if q_mop[2] and q_mop[2] != 0:
    g2 = q_p27_2 / q_mop[2]
    print("q_p27_2 / q_mop_2 = %s" % g2)

# Check gauge with a_n
g0_a = q_p27_0 / a_vals[0]
g1_a = q_p27_1 / a_vals[1]
g2_a = q_p27_2 / a_vals[2]
print("\nq_p27 / a_n:")
print("  n=0: %s" % g0_a)
print("  n=1: %s" % g1_a)
print("  n=2: %s" % g2_a)

# Factor q_p27_0
print("\n=== Factoring q_0 ===")
print("q_0 = %s" % q_p27_0)
print("factor: %s" % factor(ZZ(q_p27_0)))
print("q_1 num: %s" % factor(ZZ(-167282265043404)))
print("q_1 den: %s = %s" % (905, factor(905)))
print("q_2 num: %s" % factor(ZZ(-964185327658080)))
print("q_2 den: %s = %s" % (6071, factor(6071)))

# Compute h_n Pochhammer gauge
print("\n=== Pochhammer gauge h_n ===")
def pochhammer(a, n):
    return prod(a + j for j in range(n))

for i in range(5):
    h_i = QQ(2)**(-20*i) * pochhammer(QQ(3), i)**4 * pochhammer(QQ(4), i)**6 / \
          (pochhammer(QQ(5)/2, i)**4 * pochhammer(QQ(7)/2, i)**3 * pochhammer(QQ(9)/2, i)**3)
    print("h_%d = %s = %.6e" % (i, h_i, float(h_i)))
    if i < 3:
        print("  q_p27_%d / h_%d = %s" % (i, i, q_p27_0 / h_i if i == 0 else (q_p27_1 / h_i if i == 1 else q_p27_2 / h_i)))

# Also verify the Zudilin identity: a_n = (-1)^n Σ (-1)^k C(n,k) C(n+k,k)^3
print("\n=== Verify Zudilin identity ===")
for i in range(6):
    s1 = sum(binomial(i,k)^2 * binomial(i+k,i) * binomial(i+2*k,i) for k in range(i+1))
    s2 = (-1)^i * sum((-1)^k * binomial(i,k) * binomial(i+k,k)^3 for k in range(i+1))
    print("n=%d: Σ C² C C = %s, (-1)^n Σ (-1)^k C·C³ = %s, match: %s" % (i, s1, s2, s1 == s2))
