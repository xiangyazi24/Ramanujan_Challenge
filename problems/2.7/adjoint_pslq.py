#!/usr/bin/env python3
"""
Compute the adjoint minimal solution w^{(0)} for P2.7 to high precision,
then try PSLQ to identify w^{(0)}_0, w^{(0)}_1, w^{(0)}_2 as algebraic
combinations of known constants.

If w^{(0)} initial values have closed forms, then
  <w^{(0)}, X(e)> = 0
could be verified EXACTLY, proving c_0(e) = 0.
"""
from mpmath import mp, mpf, zeta, pi, euler, catalan, log, pslq, matrix, power
from fractions import Fraction as Q

mp.dps = 300

# === P2.7 recurrence coefficients ===
def A_c(n): return 1024*(2*n+5)**4*(2*n+7)**3*(2*n+9)**3*(946*n*n+6407*n+10860)
def B_c(n): return 128*(2*n+7)**3*(2*n+9)**3*(104060*n**6+1745370*n**5+12145238*n**4+44886481*n**3+92943995*n**2+102256019*n+46709052)
def C_c(n): return 16*(n+3)**4*(2*n+9)**3*(3784*n**5+57792*n**4+351019*n**3+1059230*n**2+1587211*n+944620)
def D_c(n): return (n+3)**4*(n+4)**6*(946*n*n+4515*n+5399)

# Monic operator: u_{n+3} + a2(n)u_{n+2} + a1(n)u_{n+1} + a0(n)u_n = 0
def a0(n): return mpf(-D_c(n)) / mpf(A_c(n))
def a1(n): return mpf(C_c(n+1)) / mpf(A_c(n+1))
def a2(n): return mpf(-B_c(n+2)) / mpf(A_c(n+2))

# === Compute adjoint minimal solution via Miller backward recurrence ===
# Adjoint: w_n = -a2(n+1)*w_{n+1} - a1(n+2)*w_{n+2} - a0(n+3)*w_{n+3}

M = 500  # backward start
w = [mpf(0)] * (M + 5)
w[M] = mpf(1)
w[M-1] = mpf(0)
w[M-2] = mpf(0)

for n in range(M-3, -4, -1):
    if n < -3:
        break
    w[n] = -a2(n+1)*w[n+1] - a1(n+2)*w[n+2] - a0(n+3)*w[n+3]

# Normalize: w[0] = 1
norm = w[0]
if norm != 0:
    for i in range(len(w)):
        w[i] /= norm

print("=== Adjoint minimal solution w^{(0)} (w[0]=1) ===")
for n in range(8):
    print(f"  w[{n}] = {mp.nstr(w[n], 50)}")

# Also do a second backward computation from M2=600 to check convergence
M2 = 600
w2 = [mpf(0)] * (M2 + 5)
w2[M2] = mpf(1)
w2[M2-1] = mpf(0)
w2[M2-2] = mpf(0)

for n in range(M2-3, -4, -1):
    if n < -3:
        break
    w2[n] = -a2(n+1)*w2[n+1] - a1(n+2)*w2[n+2] - a0(n+3)*w2[n+3]

norm2 = w2[0]
if norm2 != 0:
    for i in range(len(w2)):
        w2[i] /= norm2

print("\n=== Convergence check (M=500 vs M=600) ===")
for n in range(5):
    diff = abs(w[n] - w2[n])
    print(f"  w[{n}] diff = {mp.nstr(diff, 10)}")

# Use the M=600 result (more precise)
w = w2[:M+5] + [mpf(0)] * max(0, M+5 - len(w2))

# === Known constants for PSLQ ===
z2 = zeta(2)  # pi^2/6
z3 = zeta(3)
z23 = z2 + z3
pi2 = pi**2
pi4 = pi**4
ln2 = log(2)
G = catalan  # Catalan's constant
z4 = zeta(4)  # pi^4/90
z5 = zeta(5)

print("\n=== PSLQ identification of w[1], w[2] ===")
print(f"w[0] = 1 (normalized)")
print(f"w[1] = {mp.nstr(w[1], 80)}")
print(f"w[2] = {mp.nstr(w[2], 80)}")

# Try PSLQ: w[1] = c_1 + c_2*z2 + c_3*z3 + c_4*z2^2 + c_5*z2*z3 + c_6*z3^2
# with rational c_i (we look for integer relations)

# For w[1]:
print("\n--- PSLQ for w[1] against {1, ζ(2), ζ(3)} ---")
rel = pslq([w[1], mpf(1), z2, z3])
if rel:
    print(f"  Found: {rel[0]}*w[1] + {rel[1]} + {rel[2]}*ζ(2) + {rel[3]}*ζ(3) = 0")
    if rel[0] != 0:
        print(f"  => w[1] = {-Q(rel[1],rel[0])} + {-Q(rel[2],rel[0])}*ζ(2) + {-Q(rel[3],rel[0])}*ζ(3)")
else:
    print("  No relation found")

print("\n--- PSLQ for w[1] against {1, ζ(2), ζ(3), ζ(2)², ζ(2)ζ(3), ζ(3)²} ---")
rel = pslq([w[1], mpf(1), z2, z3, z2**2, z2*z3, z3**2])
if rel:
    print(f"  Found: {rel}")
else:
    print("  No relation found")

print("\n--- PSLQ for w[1] against {1, ζ(2), ζ(3), π², ln2} ---")
rel = pslq([w[1], mpf(1), z2, z3, pi2, ln2])
if rel:
    print(f"  Found: {rel}")
else:
    print("  No relation found")

# For w[2]:
print("\n--- PSLQ for w[2] against {1, ζ(2), ζ(3)} ---")
rel = pslq([w[2], mpf(1), z2, z3])
if rel:
    print(f"  Found: {rel[0]}*w[2] + {rel[1]} + {rel[2]}*ζ(2) + {rel[3]}*ζ(3) = 0")
    if rel[0] != 0:
        print(f"  => w[2] = {-Q(rel[1],rel[0])} + {-Q(rel[2],rel[0])}*ζ(2) + {-Q(rel[3],rel[0])}*ζ(3)")
else:
    print("  No relation found")

print("\n--- PSLQ for w[2] against {1, ζ(2), ζ(3), ζ(2)², ζ(2)ζ(3), ζ(3)²} ---")
rel = pslq([w[2], mpf(1), z2, z3, z2**2, z2*z3, z3**2])
if rel:
    print(f"  Found: {rel}")
else:
    print("  No relation found")

# For w[3]:
print("\n--- PSLQ for w[3] against {1, ζ(2), ζ(3)} ---")
rel = pslq([w[3], mpf(1), z2, z3])
if rel:
    print(f"  Found: {rel[0]}*w[3] + {rel[1]} + {rel[2]}*ζ(2) + {rel[3]}*ζ(3) = 0")
    if rel[0] != 0:
        print(f"  => w[3] = {-Q(rel[1],rel[0])} + {-Q(rel[2],rel[0])}*ζ(2) + {-Q(rel[3],rel[0])}*ζ(3)")
else:
    print("  No relation found")

# Now try a different normalization: instead of w[0]=1, try the UNMONIC version
# where we work with exact integers
print("\n\n=== Alternative: UNMONIC adjoint with integer/rational values ===")
# The unmonic recurrence is:
#   A(n) u_{n+1} - B(n) u_n + C(n-1) u_{n-1} - D(n-2) u_{n-2} = 0
#
# Adjoint of this: operate backward
# The formal adjoint of Σ c_k(n) u_{n+k} = 0 is Σ c_k(n-k) w_{n-k} = 0 (with sign conventions)
# For our recurrence L = A(n)S - B(n) + C(n-1)S^{-1} - D(n-2)S^{-2}
# Formally: L = A(n)E - B(n) + C(n-1)E^{-1} - D(n-2)E^{-2}
# Adjoint L* maps w_n to: A(n)w_{n+1} - B(n)w_n + C(n)w_{n-1} - D(n)w_{n-2}
# Wait, need to be careful.
#
# Actually for recurrence a(n)u_{n+1} + b(n)u_n + c(n)u_{n-1} + d(n)u_{n-2} = 0
# The adjoint recurrence is: a(n-1)w_{n-1} + b(n)w_n + c(n+1)w_{n+1} + d(n+2)w_{n+2} = 0
# (obtained by summing by parts / transposing the recurrence matrix)
#
# Our recurrence: A(n)u_{n+1} - B(n)u_n + C(n-1)u_{n-1} - D(n-2)u_{n-2} = 0
# So: a(n)=A(n), b(n)=-B(n), c(n)=C(n-1), d(n)=-D(n-2)
# Adjoint: A(n-1)w_{n-1} - B(n)w_n + C(n)w_{n+1} - D(n)w_{n+2} = 0
# Or: D(n)w_{n+2} = C(n)w_{n+1} - B(n)w_n + A(n-1)w_{n-1}
#
# For backward: w_{n-1} = [D(n)w_{n+2} - C(n)w_{n+1} + B(n)w_n] / A(n-1)

# Compute unmonic adjoint backward
M3 = 500
wu = [mpf(0)] * (M3 + 5)
wu[M3] = mpf(1)
wu[M3-1] = mpf(0)
wu[M3-2] = mpf(0)

for n in range(M3-1, 0, -1):
    # w_{n-1} = [D(n)*w_{n+2} - C(n)*w_{n+1} + B(n)*w_n] / A(n-1)
    wu[n-1] = (mpf(D_c(n))*wu[n+2] - mpf(C_c(n))*wu[n+1] + mpf(B_c(n))*wu[n]) / mpf(A_c(n-1))

print("Unmonic adjoint w (backward from M=500):")
for n in range(8):
    print(f"  wu[{n}] = {mp.nstr(wu[n], 50)}")

# Try PSLQ on unmonic adjoint
# First normalize so wu[0] = 1
norm_u = wu[0]
for i in range(len(wu)):
    wu[i] /= norm_u

print(f"\n  wu[0] = 1 (normalized)")
print(f"  wu[1] = {mp.nstr(wu[1], 80)}")
print(f"  wu[2] = {mp.nstr(wu[2], 80)}")

print("\n--- PSLQ for unmonic wu[1] against {1, ζ(2), ζ(3)} ---")
rel = pslq([wu[1], mpf(1), z2, z3])
if rel:
    print(f"  Found: {rel[0]}*wu[1] + {rel[1]} + {rel[2]}*ζ(2) + {rel[3]}*ζ(3) = 0")
else:
    print("  No relation found")

print("\n--- PSLQ for unmonic wu[2] against {1, ζ(2), ζ(3)} ---")
rel = pslq([wu[2], mpf(1), z2, z3])
if rel:
    print(f"  Found: {rel[0]}*wu[2] + {rel[1]} + {rel[2]}*ζ(2) + {rel[3]}*ζ(3) = 0")
else:
    print("  No relation found")

# Try with denominator clearing — maybe w values are small rationals times known constants
# Try PSLQ with larger basis
print("\n--- PSLQ for wu[1] against {1, ζ(2), ζ(3), ζ(4), ζ(5), ζ(2)ζ(3)} ---")
rel = pslq([wu[1], mpf(1), z2, z3, z4, z5, z2*z3])
if rel:
    print(f"  Found: {rel}")
else:
    print("  No relation found")

# Try recognizing w[1]/w[0] as a simple rational number
# (since the recurrence has rational coefficients, if the adjoint minimal
# solution is uniquely characterized, its ratios might be rational)
print("\n--- Checking if wu[1] is rational (PSLQ against {1}) ---")
# Check: is wu[1] close to a rational with small denominator?
for d in range(1, 200):
    frac = wu[1] * d
    rounded = round(float(frac))
    if abs(frac - rounded) < mpf(10)**(-250):
        print(f"  wu[1] ≈ {rounded}/{d}")
        break

# Check bilinear concomitant with unmonic version
print("\n\n=== Bilinear concomitant J(wu, e) ===")

# Forward solutions
N = 120
q = [mpf(0)] * (N + 10)
q[0] = mpf(-215040420000)
q[1] = mpf(Q(-167282265043404, 905))
q[2] = mpf(Q(-964185327658080, 6071))
for i in range(3, N + 5):
    n = i - 1
    q[i] = mpf(B_c(n))/mpf(A_c(n))*q[i-1] - mpf(C_c(n-1))/mpf(A_c(n-1))*q[i-2] + mpf(D_c(n-2))/mpf(A_c(n-2))*q[i-3]

p = [mpf(0)] * (N + 10)
p[0] = mpf(-612218384750)
p[1] = mpf(Q(-9525021973931919, 18100))
p[2] = mpf(Q(-29561828382772029, 65380))
for i in range(3, N + 5):
    n = i - 1
    p[i] = mpf(B_c(n))/mpf(A_c(n))*p[i-1] - mpf(C_c(n-1))/mpf(A_c(n-1))*p[i-2] + mpf(D_c(n-2))/mpf(A_c(n-2))*p[i-3]

e = [p[n] - z23*q[n] for n in range(N+5)]

# Bilinear concomitant for the monic operator
def J_monic(w, u, n):
    return (w[n-1]*u[n+2]
            + (w[n-2] + a2(n-1)*w[n-1])*u[n+1]
            + (w[n-3] + a2(n-2)*w[n-2] + a1(n-1)*w[n-1])*u[n])

# Use the monic adjoint w (not wu) for the monic bilinear form
# Re-compute monic adjoint from M=600
w_m = w2[:M+5] + [mpf(0)] * max(0, M+5 - len(w2))
norm_m = w_m[0]
for i in range(len(w_m)):
    w_m[i] /= norm_m

print("J(w, e) at various n:")
for n in [5, 10, 20, 30, 40, 50]:
    if n+2 < N+5 and n-3 >= 0 and n < len(w_m):
        val = J_monic(w_m, e, n)
        print(f"  n={n}: J = {mp.nstr(val, 15)}")

print("\nJ(w, q) at n=10:")
Jq = J_monic(w_m, q, 10)
print(f"  J(w,q) = {mp.nstr(Jq, 50)}")

print("\nJ(w, p) at n=10:")
Jp = J_monic(w_m, p, 10)
print(f"  J(w,p) = {mp.nstr(Jp, 50)}")

print(f"\nJ(w,p)/J(w,q) = {mp.nstr(Jp/Jq, 50)}")
print(f"ζ(2)+ζ(3)    = {mp.nstr(z23, 50)}")
print(f"Difference   = {mp.nstr(Jp/Jq - z23, 15)}")

# Try to PSLQ the concomitant values themselves
print("\n\n=== PSLQ on J(w, q) ===")
# J(w,q) should be some specific constant. Let's try to identify it.
# Try against common constants
print(f"J(w,q) = {mp.nstr(Jq, 80)}")

# It might be a product of Pochhammer-type rationals times pi powers etc.
# But more likely it's just a specific number determined by the normalization.
# Since we normalized w[0]=1, the value J(w,q) depends on the initial values.

# More interesting: the RATIO J(w,p)/J(w,q) = ζ(2)+ζ(3) is what matters.
# This is already verified to high precision. What we need is J(w,e) = 0.

# Let's compute |J(w,e)| / |J(w,q)| to see the relative magnitude
print(f"\n|J(w,e)| / |J(w,q)| = {mp.nstr(abs(J_monic(w_m, e, 10)) / abs(Jq), 15)}")
