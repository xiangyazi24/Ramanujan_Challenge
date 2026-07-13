#!/usr/bin/env python3
"""Problem 2.5: Extract scalar recurrence from the 3x3 CMF matrix.
Compute det(M(n)) and trace patterns to identify the underlying structure."""
from mpmath import mp, mpf, matrix, catalan, det

mp.dps = 50

def M(n):
    """3x3 matrix M(n) from Problem 2.5."""
    m11 = (-2*n-5)*(n+3)**2*(136*n**4+1424*n**3+5548*n**2+9551*n+6141)
    m12 = 384*n**6+6384*n**5+44168*n**4+162698*n**3+336377*n**2+369933*n+169011
    m13 = -480*n**4-4980*n**3-19210*n**2-32690*n-20730
    m21 = (n+2)**2*(n+3)**2*(4*n+10)*(48*n**3+386*n**2+1017*n+879)
    m22 = (n+2)**2*(-272*n**5-3848*n**4-21732*n**3-61184*n**2-85761*n-47808)
    m23 = (n+2)**2*(320*n**3+2540*n**2+6610*n+5640)
    m31 = (-4*n-10)*(n+2)**2*(n+3)**2*(32*n**4+302*n**3+1037*n**2+1530*n+813)
    m32 = (n+2)**2*(192*n**6+2984*n**5+19116*n**4+64452*n**3+120256*n**2+117279*n+46476)
    m33 = (n+2)**2*(-16*n**5-408*n**4-2912*n**3-8884*n**2-12254*n-6240)
    return matrix([[mpf(m11),mpf(m12),mpf(m13)],
                   [mpf(m21),mpf(m22),mpf(m23)],
                   [mpf(m31),mpf(m32),mpf(m33)]])

# Compute det(M(n)) for several n
print("=== det(M(n)) for Problem 2.5 ===")
for n in range(10):
    d = det(M(n))
    print(f"  det(M({n})) = {mp.nstr(d, 15)}")

# Factor pattern: check if det(M(n)) = c * (n+2)^a * (n+3)^b * poly(n)
print("\n=== det(M(n)) / ((n+2)^6 * (n+3)^4) ===")
for n in range(10):
    d = det(M(n))
    normalized = d / ((n+2)**6 * (n+3)**4)
    print(f"  n={n}: {mp.nstr(normalized, 15)}")

# Compute P_{N,j}/Q_{N,j} for small N to extract the scalar recurrence
A = matrix([[mpf(30921), mpf(-32972), mpf(8240)],
            [mpf(33750), mpf(-36000), mpf(9000)]])

prod = matrix([[mpf(1),mpf(0),mpf(0)],
               [mpf(0),mpf(1),mpf(0)],
               [mpf(0),mpf(0),mpf(1)]])

P_vals = []  # P_{N,1} values
Q_vals = []  # Q_{N,1} values

for N in range(20):
    result = A * prod
    P_vals.append(result[0,0])
    Q_vals.append(result[1,0])
    prod = prod * M(N)

# Append N=20
result = A * prod
P_vals.append(result[0,0])
Q_vals.append(result[1,0])

print("\n=== Scalar recurrence test for Q_{N,1} ===")
print("Testing: a(N)*Q_N + b(N)*Q_{N-1} + c(N)*Q_{N-2} + d(N)*Q_{N-3} = 0")
# For an order-3 recurrence, we need at least 4 consecutive Q values
# to determine the coefficients (up to normalization)
for start in range(2, 15):
    # Solve for (a,b,c) in: a*Q[n+3] + b*Q[n+2] + c*Q[n+1] + Q[n] = 0
    # Three equations for three unknowns (a,b,c), normalizing d=1
    n = start
    M_sys = matrix([
        [Q_vals[n+3], Q_vals[n+2], Q_vals[n+1]],
        [Q_vals[n+4], Q_vals[n+3], Q_vals[n+2]],
        [Q_vals[n+5], Q_vals[n+4], Q_vals[n+3]],
    ])
    rhs_sys = matrix([[-Q_vals[n]], [-Q_vals[n+1]], [-Q_vals[n+2]]])
    try:
        sol = M_sys**(-1) * rhs_sys
        a, b, c = sol[0,0], sol[1,0], sol[2,0]
        # Verify with next value
        residual = a*Q_vals[n+6] + b*Q_vals[n+5] + c*Q_vals[n+4] + Q_vals[n+3]
        print(f"  n={n}: a={mp.nstr(a,8)}, b={mp.nstr(b,8)}, c={mp.nstr(c,8)}, residual={mp.nstr(residual,6)}")
    except:
        print(f"  n={n}: singular system")

# Also check ratios Q_{N+1}/Q_N for Poincaré root
print("\n=== Q_{N+1}/Q_N (Poincaré root) ===")
for N in range(5, 20):
    if Q_vals[N] != 0:
        ratio = Q_vals[N+1] / Q_vals[N]
        print(f"  N={N}: Q_{{N+1}}/Q_N = {mp.nstr(ratio/N**2, 10)}")

print(f"\nTarget G = {catalan}")
print(f"P_20/Q_20 = {P_vals[20]/Q_vals[20]}")
print(f"diff = {P_vals[20]/Q_vals[20] - catalan}")
