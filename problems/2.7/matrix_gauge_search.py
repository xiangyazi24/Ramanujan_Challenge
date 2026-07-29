#!/usr/bin/env python3
"""
Search for matrix gauge G(n) ∈ GL_3(Q(n)) such that
  C_P^hat(n) · G(n) = G(n+1) · C_Z(n)

where:
- C_P^hat(n) is the companion matrix for q̂_n = 64^n q_n (P2.7 scaled)
- C_Z(n) is the companion matrix for b_n (Zudilin)

If such G exists, the P2.7 and Zudilin difference modules are isomorphic,
and the error bound transfers via the gauge.

Method: compute G(n) = C_P^hat(n)^{-1} · ... numerically at many n values,
then identify it as a matrix of rational functions by Padé approximation.
"""
from fractions import Fraction as F
import numpy as np

# === P2.7 coefficients ===
def A_c(n): return 1024*(2*n+5)**4*(2*n+7)**3*(2*n+9)**3*(946*n*n+6407*n+10860)
def B_c(n): return 128*(2*n+7)**3*(2*n+9)**3*(104060*n**6+1745370*n**5+12145238*n**4+44886481*n**3+92943995*n**2+102256019*n+46709052)
def C_c(n): return 16*(n+3)**4*(2*n+9)**3*(3784*n**5+57792*n**4+351019*n**3+1059230*n**2+1587211*n+944620)
def D_c(n): return (n+3)**4*(n+4)**6*(946*n*n+4515*n+5399)

# === Zudilin recurrence coefficients ===
# b_{n+1} = MZ(n)/QZ(n)/(2n+1)/(n+1)^3 * b_n
#          - n*NZ(n)/QZ(n)/(2n+1)/(n+1)^3 * b_{n-1}
#          + RZ(n+1)*n*(n-1)^3/[2*QZ(n)*(2n+1)*(n+1)^3] * b_{n-2}

# Zudilin recurrence in forward form b_{n+3} = α'_n b_{n+2} + β'_n b_{n+1} + γ'_n b_n
# From Zudilin eq 6.4 (shifted by n → n+2):
# 2(946(n+2)²-731(n+2)+153)(2(n+2)+1)((n+2)+1)³ b_{n+3}
# - 2·MZ(n+2) b_{n+2}
# + 2(n+2)·NZ(n+2) b_{n+1}
# - (946(n+2)²+1161(n+2)+368)(n+2)(n+1)³ b_n = 0

def QZ(n): return 946*n**2 - 731*n + 153
def MZ(n): return 104060*n**6+127710*n**5+12788*n**4-34525*n**3-8482*n**2+3298*n+1071
def NZ(n): return 3784*n**5-1032*n**4-1925*n**3+853*n**2+328*n-184
def RZ(n): return 946*n**2+1161*n+368

def zudilin_companion(n):
    """Companion matrix for Zudilin's recurrence for b_n (shifted to forward form).

    b_{n+3} = α' b_{n+2} + β' b_{n+1} + γ' b_n

    From: 2·QZ(n+2)·(2n+5)·(n+3)³ · b_{n+3}
        = 2·MZ(n+2) · b_{n+2}
        - 2·(n+2)·NZ(n+2) · b_{n+1}
        + RZ(n+2)·(n+2)·(n+1)³ · b_n
    """
    m = n + 2  # shift
    denom = F(2 * QZ(m) * (2*m+1) * (m+1)**3)
    alpha = F(2 * MZ(m), denom.denominator) / denom * denom.denominator if denom != 0 else F(0)

    # More carefully:
    d = 2 * QZ(m) * (2*m+1) * (m+1)**3
    a = F(2 * MZ(m), d)
    b = F(-2 * m * NZ(m), d)
    c = F(RZ(m) * m * (m-1)**3, d)

    return np.array([[float(a), float(b), float(c)],
                     [1.0, 0.0, 0.0],
                     [0.0, 1.0, 0.0]])

def p27_companion(n):
    """Companion matrix for P2.7 monic recurrence.

    u_{n+3} = (B(n+2)/A(n+2)) u_{n+2} - (C(n+1)/A(n+1)) u_{n+1} + (D(n)/A(n)) u_n
    """
    a = float(F(B_c(n+2), A_c(n+2)))
    b = float(F(-C_c(n+1), A_c(n+1)))
    c = float(F(D_c(n), A_c(n)))

    return np.array([[a, b, c],
                     [1.0, 0.0, 0.0],
                     [0.0, 1.0, 0.0]])

def p27_hat_companion(n):
    """Companion matrix for q̂_n = 64^n q_n.

    C_hat(n) = 64 · D · C_P(n) · D^{-1}
    where D = diag(64², 64, 1).
    """
    C = p27_companion(n)
    # Scale: multiply row 0 by 64, col 0 by 64², col 1 by 64, col 2 by 1
    # Then divide row 0 by 64², row 1 by 64, row 2 by 1

    # C_hat = 64 * D * C * D^{-1}
    # D = diag(4096, 64, 1), D^{-1} = diag(1/4096, 1/64, 1)
    D = np.diag([4096.0, 64.0, 1.0])
    Dinv = np.diag([1.0/4096, 1.0/64, 1.0])
    return 64.0 * D @ C @ Dinv

# Verify: compute q̂_n terms from both representations
print("=== Verify companion matrices ===")

# Zudilin b_n
b = [F(1), F(7), F(163)]
for n in range(2, 15):
    m = n
    d = 2 * QZ(m) * (2*m+1) * (m+1)**3
    bn3 = F(2*MZ(m), d) * b[n] + F(-2*m*NZ(m), d) * b[n-1] + F(RZ(m)*m*(m-1)**3, d) * b[n-2]
    b.append(bn3)

print(f"b[0:6] = {[float(x) for x in b[:6]]}")

# P2.7 q_n and q̂_n
q = [F(-215040420000), F(-167282265043404, 905), F(-964185327658080, 6071)]
for n in range(2, 15):
    qn1 = F(B_c(n), A_c(n)) * q[n] + F(-C_c(n-1), A_c(n-1)) * q[n-1] + F(D_c(n-2), A_c(n-2)) * q[n-2]
    q.append(qn1)

qhat = [F(64)**n * q[n] for n in range(15)]
print(f"q̂[0:4] = {[float(x) for x in qhat[:4]]}")

# Verify companion matrix multiplication
C_Z_0 = zudilin_companion(0)
state_Z = np.array([float(b[2]), float(b[1]), float(b[0])])
for n in range(3):
    C = zudilin_companion(n)
    new_state = C @ state_Z[:3] if n == 0 else C @ state_Z
    state_Z = new_state
    print(f"  Zudilin state after step {n}: b_{n+3} = {state_Z[0]:.6f}, actual = {float(b[n+3]):.6f}")

# === Compute G(n) at several values of n ===
# If G exists: G(n+1) = C_hat(n) · G(n) · C_Z(n)^{-1}
# Start with G(0) = identity (just to see what happens)

print("\n=== Computing G(n) chain ===")

# Actually, we need to find G(n) such that the SOLUTION VECTORS match.
# q̂_n = first component of G(n) · [b_{n+2}, b_{n+1}, b_n]^T
# So G(n) · state_Z(n) = state_P(n) (the scaled P2.7 state)

# Method: compute G(n) at each n from the state vectors.
# We need THREE solution pairs to determine the 3x3 matrix.
# Zudilin has three solutions: b_n, b̃_n, b̃̃_n
# P2.7 has three solutions: q_n, p_n, ???

# Zudilin companions
bt = [F(0), F(23,2), F(2145,8)]  # b̃
btt = [F(0), F(17,2), F(3135,16)]  # b̃̃

for n in range(2, 15):
    m = n
    d = 2 * QZ(m) * (2*m+1) * (m+1)**3
    bt.append(F(2*MZ(m), d) * bt[n] + F(-2*m*NZ(m), d) * bt[n-1] + F(RZ(m)*m*(m-1)**3, d) * bt[n-2])
    btt.append(F(2*MZ(m), d) * btt[n] + F(-2*m*NZ(m), d) * btt[n-1] + F(RZ(m)*m*(m-1)**3, d) * btt[n-2])

# P2.7 numerator p_n
p = [F(-612218384750), F(-9525021973931919, 18100), F(-29561828382772029, 65380)]
for n in range(2, 15):
    pn1 = F(B_c(n), A_c(n)) * p[n] + F(-C_c(n-1), A_c(n-1)) * p[n-1] + F(D_c(n-2), A_c(n-2)) * p[n-2]
    p.append(pn1)

phat = [F(64)**n * p[n] for n in range(15)]

# Third P2.7 solution: we can use e_n = p_n - L*q_n as a third solution
# But we'd need L = ζ(2)+ζ(3) which introduces irrationals.
# Alternative: use any linear combination, e.g. the trivial one [1, 0, 0, ...]
# extended to satisfy the recurrence.

# Actually, the third solution is determined by the initial conditions.
# Let's use: s_0 = 1, s_1 = 0, s_2 = 0 (gives a specific P2.7 solution)
s = [F(1), F(0), F(0)]
for n in range(2, 15):
    sn1 = F(B_c(n), A_c(n)) * s[n] + F(-C_c(n-1), A_c(n-1)) * s[n-1] + F(D_c(n-2), A_c(n-2)) * s[n-2]
    s.append(sn1)

shat = [F(64)**n * s[n] for n in range(15)]

# Now compute G(n) as the matrix mapping [b_{n+2}, b_{n+1}, b_n; bt_{n+2}, ...; btt_{n+2}, ...]
# to [qhat_{n+2}, qhat_{n+1}, qhat_n; phat_{n+2}, ...; shat_{n+2}, ...]

# At index n, the Zudilin state matrix (columns = solutions):
# Z(n) = [[b_{n+2}, bt_{n+2}, btt_{n+2}],
#          [b_{n+1}, bt_{n+1}, btt_{n+1}],
#          [b_n,     bt_n,     btt_n    ]]

# The P2.7 state matrix:
# P(n) = [[qhat_{n+2}, phat_{n+2}, shat_{n+2}],
#          [qhat_{n+1}, phat_{n+1}, shat_{n+1}],
#          [qhat_n,     phat_n,     shat_n    ]]

# Then G(n) = P(n) · Z(n)^{-1}

print("\n=== G(n) at n=0,1,...,8 ===")
for n in range(9):
    Z = np.array([[float(b[n+2]),  float(bt[n+2]),  float(btt[n+2])],
                  [float(b[n+1]),  float(bt[n+1]),  float(btt[n+1])],
                  [float(b[n]),    float(bt[n]),    float(btt[n])]])

    P = np.array([[float(qhat[n+2]), float(phat[n+2]), float(shat[n+2])],
                  [float(qhat[n+1]), float(phat[n+1]), float(shat[n+1])],
                  [float(qhat[n]),   float(phat[n]),   float(shat[n])]])

    try:
        G = P @ np.linalg.inv(Z)
        print(f"\nG({n}):")
        for i in range(3):
            print(f"  [{G[i,0]:20.6f}  {G[i,1]:20.6f}  {G[i,2]:20.6f}]")

        # Check: is G(n) = G(0) for all n? (Would mean constant gauge)
        if n == 0:
            G0 = G.copy()
        else:
            diff = np.max(np.abs(G - G0))
            rel = diff / (np.max(np.abs(G)) + 1e-30)
            print(f"  max |G({n})-G(0)| = {diff:.6e}, relative = {rel:.6e}")
    except np.linalg.LinAlgError:
        print(f"\nG({n}): Zudilin matrix singular")

# Check if G(n) entries are rational functions of n
print("\n\n=== Checking if G(n) entries are rational in n ===")
G_values = []
for n in range(12):
    Z = np.array([[float(b[n+2]),  float(bt[n+2]),  float(btt[n+2])],
                  [float(b[n+1]),  float(bt[n+1]),  float(btt[n+1])],
                  [float(b[n]),    float(bt[n]),    float(btt[n])]])
    P = np.array([[float(qhat[n+2]), float(phat[n+2]), float(shat[n+2])],
                  [float(qhat[n+1]), float(phat[n+1]), float(shat[n+1])],
                  [float(qhat[n]),   float(phat[n]),   float(shat[n])]])
    try:
        G = P @ np.linalg.inv(Z)
        G_values.append(G)
    except:
        G_values.append(None)

# For each entry (i,j), fit a rational function p(n)/q(n)
# with deg(p) <= d, deg(q) <= d
print("\nFitting rational functions to G(n) entries (deg=3):")
for i in range(3):
    for j in range(3):
        vals = [(n, G_values[n][i,j]) for n in range(len(G_values)) if G_values[n] is not None]
        ns = np.array([v[0] for v in vals])
        gs = np.array([v[1] for v in vals])

        # Try polynomial fit first
        for deg in range(1, 6):
            if len(vals) < deg + 2:
                break
            coeffs = np.polyfit(ns[:deg+2], gs[:deg+2], deg)
            pred = np.polyval(coeffs, ns)
            resid = np.max(np.abs(pred - gs))
            if resid < 1e-3:
                print(f"  G[{i},{j}]: polynomial deg {deg}, max resid = {resid:.2e}")
                print(f"    coeffs = {coeffs}")
                break
        else:
            # No good polynomial fit, try rational
            print(f"  G[{i},{j}]: no polynomial fit up to deg 5, values = {gs[:5]}")
