"""
Test whether ĥ_n satisfies the P2.7 recurrence.
If so, q_n and ĥ_n are solutions of the SAME operator, and the
intertwiner problem reduces to initial-value matching.
"""
from fractions import Fraction as Q
from math import comb, factorial

def pochhammer(a, n):
    result = Q(1)
    for i in range(n):
        result *= Q(a) + i
    return result

def A_c(n):
    n = Q(n)
    return (1024*(2*n+5)**4*(2*n+7)**3*(2*n+9)**3*(946*n*n+6407*n+10860))

def B_c(n):
    n = Q(n)
    return (128*(2*n+7)**3*(2*n+9)**3*(104060*n**6+1745370*n**5+12145238*n**4+44886481*n**3+92943995*n**2+102256019*n+46709052))

def C_c(n):
    n = Q(n)
    return (16*(n+3)**4*(2*n+9)**3*(3784*n**5+57792*n**4+351019*n**3+1059230*n**2+1587211*n+944620))

def D_c(n):
    n = Q(n)
    return ((n+3)**4*(n+4)**6*(946*n*n+4515*n+5399))

def compute_b(n):
    return sum(comb(n,k)**2 * comb(n+k,n) * comb(n+2*k,n) for k in range(n+1))

N = 40
b = [Q(compute_b(n)) for n in range(N)]

h_hat = []
for n in range(N):
    val = Q(1, 64**n) * pochhammer(Q(5,2), n) / Q(factorial(n)) * b[n]
    h_hat.append(val)

# Test: does ĥ_n satisfy u_{n+1} = (B(n)/A(n))u_n - (C(n-1)/A(n-1))u_{n-1} + (D(n-2)/A(n-2))u_{n-2}?
print("=== Testing if ĥ_n satisfies P2.7 recurrence ===")
for n in range(2, 35):
    pred = Q(B_c(n), A_c(n)) * h_hat[n] - Q(C_c(n-1), A_c(n-1)) * h_hat[n-1] + Q(D_c(n-2), A_c(n-2)) * h_hat[n-2]
    actual = h_hat[n+1]
    residual = pred - actual
    if residual == 0:
        print(f"  n={n}: EXACT ZERO residual")
    else:
        rel = float(abs(residual / actual)) if actual != 0 else float('inf')
        print(f"  n={n}: residual={float(residual):.6e}, relative={rel:.6e}")

# If not, compute the TWISTED Zudilin recurrence for ĥ_n
# ĥ_n = g(n) * b_n where g(n) = 64^{-n} * (5/2)_n / n!
# If b_n satisfies Zudilin recurrence with coefficients α_0(n)b_{n-2} + α_1(n)b_{n-1} + α_2(n)b_n + α_3(n)b_{n+1} = 0
# Then ĥ_n satisfies the gauged version with coefficients α_i(n) * g(n+i-2)/g(n) etc.

# Zudilin recurrence: from Q5175 eq (6.4)
# 2(946n²-731n+153)(2n+1)(n+1)³ b_{n+1}
# - 2(104060n⁶+127710n⁵+12788n⁴-34525n³-8482n²+3298n+1071) b_n
# + 2n(3784n⁵-1032n⁴-1925n³+853n²+328n-184) b_{n-1}
# - (946n²+1161n+368)n(n-1)³ b_{n-2} = 0

def Z3(n):
    """Coefficient of b_{n+1} in Zudilin recurrence"""
    n = Q(n)
    return 2*(946*n**2-731*n+153)*(2*n+1)*(n+1)**3

def Z2(n):
    """Coefficient of b_n"""
    n = Q(n)
    return -2*(104060*n**6+127710*n**5+12788*n**4-34525*n**3-8482*n**2+3298*n+1071)

def Z1(n):
    """Coefficient of b_{n-1}"""
    n = Q(n)
    return 2*n*(3784*n**5-1032*n**4-1925*n**3+853*n**2+328*n-184)

def Z0(n):
    """Coefficient of b_{n-2}"""
    n = Q(n)
    return -(946*n**2+1161*n+368)*n*(n-1)**3

# Verify Zudilin recurrence on b_n
print("\n=== Verifying Zudilin recurrence on b_n ===")
for n in range(2, 20):
    res = Z3(n)*b[n+1] + Z2(n)*b[n] + Z1(n)*b[n-1] + Z0(n)*b[n-2]
    print(f"  n={n}: residual = {res}")

# Now compute the recurrence for ĥ_n
# ĥ_n = g(n) * b_n, g(n) = 64^{-n} * (5/2)_n / n!
# g(n)/g(n-1) = 64^{-1} * (5/2+n-1)/n = (n+3/2) / (64n)
# g(n)/g(n-2) = g(n)/g(n-1) * g(n-1)/g(n-2) = (n+3/2)/(64n) * (n+1/2)/(64(n-1))
# = (n+3/2)(n+1/2) / (64² * n(n-1))
# g(n+1)/g(n-2) = g(n+1)/g(n) * g(n)/g(n-1) * g(n-1)/g(n-2)
# = (n+5/2)/(64(n+1)) * (n+3/2)/(64n) * (n+1/2)/(64(n-1))
# = (n+5/2)(n+3/2)(n+1/2) / (64³ * (n+1)*n*(n-1))

# The recurrence Z3(n)·g(n+1)·ĥ_{n+1}/g(n+1) + Z2(n)·g(n)·ĥ_n/g(n) + ... = 0
# becomes Z3(n)·(ĥ_{n+1}/g(n+1))·...
# Actually: b_k = ĥ_k/g(k), so substitute:
# Z3(n)·ĥ_{n+1}/g(n+1) + Z2(n)·ĥ_n/g(n) + Z1(n)·ĥ_{n-1}/g(n-1) + Z0(n)·ĥ_{n-2}/g(n-2) = 0
# Multiply through by g(n):
# Z3(n)·g(n)/g(n+1)·ĥ_{n+1} + Z2(n)·ĥ_n + Z1(n)·g(n)/g(n-1)·ĥ_{n-1} + Z0(n)·g(n)/g(n-2)·ĥ_{n-2} = 0
#
# g(n)/g(n+1) = 64(n+1)/(n+5/2) = 64(n+1)/(2n+5) * 2 = 128(n+1)/(2n+5)
# g(n)/g(n-1) = 64n/(n+3/2) = 64n/(2n+3) * 2 = 128n/(2n+3)
# g(n)/g(n-2) = 64²*n(n-1)/((n+3/2)(n+1/2)) = 4096*n(n-1)/((2n+3)(2n+1)) * 4 = 16384*n(n-1)/((2n+3)(2n+1))

print("\n=== Computing twisted recurrence coefficients ===")

def g_ratio_forward(n):
    """g(n)/g(n+1) = 64(n+1)/(n+3/2) where the +3/2 comes from (5/2)_n"""
    n = Q(n)
    return Q(64)*(n+1) / (n + Q(5,2))

def g_ratio_back1(n):
    """g(n)/g(n-1) = 64n/(n+3/2-1) = 64n/(n+1/2)... wait let me recompute"""
    # g(n) = 64^{-n} * (5/2)_n / n!
    # g(n-1) = 64^{-(n-1)} * (5/2)_{n-1} / (n-1)!
    # g(n)/g(n-1) = 64^{-1} * (5/2)_n/(5/2)_{n-1} * (n-1)!/n!
    #             = (1/64) * (5/2 + n - 1) * (1/n)
    #             = (n + 3/2) / (64n)
    n = Q(n)
    return (n + Q(3,2)) / (Q(64)*n)

def g_ratio_back2(n):
    """g(n)/g(n-2) = g(n)/g(n-1) * g(n-1)/g(n-2)"""
    n = Q(n)
    r1 = (n + Q(3,2)) / (Q(64)*n)
    r2 = (n - 1 + Q(3,2)) / (Q(64)*(n-1))
    return r1 * r2

# So twisted recurrence: ĥ_{n+1} coefficient, ĥ_n, ĥ_{n-1}, ĥ_{n-2}
# Z3(n) * g(n)/g(n+1) * ĥ_{n+1} + Z2(n) * ĥ_n + Z1(n) * g(n)/g(n-1)^{-1}...
# Wait, I need g(n)/g(n-1) vs g(n-1)/g(n). Let me redo.
#
# Original: Z3(n)*b_{n+1} + Z2(n)*b_n + Z1(n)*b_{n-1} + Z0(n)*b_{n-2} = 0
# b_k = ĥ_k / g(k)
# Z3(n)*ĥ_{n+1}/g(n+1) + Z2(n)*ĥ_n/g(n) + Z1(n)*ĥ_{n-1}/g(n-1) + Z0(n)*ĥ_{n-2}/g(n-2) = 0
# Multiply by g(n):
# Z3(n)*(g(n)/g(n+1))*ĥ_{n+1} + Z2(n)*ĥ_n + Z1(n)*(g(n)/g(n-1))*ĥ_{n-1} + Z0(n)*(g(n)/g(n-2))*ĥ_{n-2} = 0

# Verify this on ĥ_n
print("\n=== Verifying twisted Zudilin recurrence on ĥ_n ===")
for n in range(2, 20):
    gn_over_gnp1 = g_ratio_forward(n)  # = 64(n+1)/(n+5/2)
    # Wait, g(n)/g(n+1): g(n+1)/g(n) = (n+1+3/2)/(64*(n+1)) = (n+5/2)/(64(n+1))
    # So g(n)/g(n+1) = 64(n+1)/(n+5/2)

    # g(n)/g(n-1) = 1/[g(n-1)/g(n)] = 1/[(n+3/2)/(64n)]^{-1}...
    # g(n-1)/g(n) = 64n/(n+3/2)... no.
    # g(n)/g(n-1) = (1/64) * (n+3/2)/n  [from the definition]
    # Hmm wait: g(n) = 64^{-n} * (5/2)_n / n!
    # g(n)/g(n-1) = (64^{-n}/64^{-(n-1)}) * ((5/2)_n/(5/2)_{n-1}) * ((n-1)!/n!)
    #            = (1/64) * (5/2+n-1) * (1/n) = (n+3/2)/(64n)
    gn_over_gnm1 = (Q(n) + Q(3,2)) / (Q(64)*Q(n))
    gn_over_gnm2 = gn_over_gnm1 * (Q(n-1) + Q(3,2)) / (Q(64)*Q(n-1))
    gn_over_gnp1_v2 = Q(64)*(Q(n)+1) / (Q(n) + Q(5,2))

    res = (Z3(n) * gn_over_gnp1_v2 * h_hat[n+1]
           + Z2(n) * h_hat[n]
           + Z1(n) * gn_over_gnm1 * h_hat[n-1]
           + Z0(n) * gn_over_gnm2 * h_hat[n-2])
    print(f"  n={n}: residual = {res}")

# Now: compare the twisted Zudilin recurrence with the P2.7 recurrence
# Twisted Zudilin: coefficients for [ĥ_{n-2}, ĥ_{n-1}, ĥ_n, ĥ_{n+1}] = 0
# P2.7: ĥ_{n+1} = (B(n)/A(n))ĥ_n - (C(n-1)/A(n-1))ĥ_{n-1} + (D(n-2)/A(n-2))ĥ_{n-2}
# i.e., -D(n-2)/A(n-2)*ĥ_{n-2} + C(n-1)/A(n-1)*ĥ_{n-1} - B(n)/A(n)*ĥ_n + ĥ_{n+1} = 0

# Twisted Zudilin in monic form (divide by coeff of ĥ_{n+1}):
# ĥ_{n+1} = -Z2(n)/(Z3(n)*r_+1) * ĥ_n - Z1(n)*r_{-1}/(Z3(n)*r_+1) * ĥ_{n-1} - Z0(n)*r_{-2}/(Z3(n)*r_+1) * ĥ_{n-2}
# where r_+1 = g(n)/g(n+1), r_{-1} = g(n)/g(n-1), r_{-2} = g(n)/g(n-2)

print("\n=== Comparing P2.7 vs Twisted Zudilin monic coefficients ===")
for n in range(2, 15):
    gn_gnp1 = Q(64)*(Q(n)+1) / (Q(n) + Q(5,2))
    gn_gnm1 = (Q(n) + Q(3,2)) / (Q(64)*Q(n))
    gn_gnm2 = gn_gnm1 * (Q(n-1) + Q(3,2)) / (Q(64)*Q(n-1))

    # Twisted Zudilin monic coefficients
    tw_c2 = -Z2(n) / (Z3(n) * gn_gnp1)  # coeff of ĥ_n
    tw_c1 = -Z1(n) * gn_gnm1 / (Z3(n) * gn_gnp1)  # coeff of ĥ_{n-1}
    tw_c0 = -Z0(n) * gn_gnm2 / (Z3(n) * gn_gnp1)  # coeff of ĥ_{n-2}

    # P2.7 monic coefficients
    p27_c2 = Q(B_c(n), A_c(n))      # coeff of u_n (which is ĥ_n if ĥ satisfies P2.7)
    p27_c1 = -Q(C_c(n-1), A_c(n-1))  # coeff of u_{n-1}
    p27_c0 = Q(D_c(n-2), A_c(n-2))   # coeff of u_{n-2}

    print(f"n={n}:")
    print(f"  P2.7:    c0={float(p27_c0):.10e}, c1={float(p27_c1):.10e}, c2={float(p27_c2):.10e}")
    print(f"  Tw.Zud.: c0={float(tw_c0):.10e}, c1={float(tw_c1):.10e}, c2={float(tw_c2):.10e}")
    print(f"  Ratios:  c0={float(p27_c0/tw_c0):.10e}, c1={float(p27_c1/tw_c1):.10e}, c2={float(p27_c2/tw_c2):.10e}")
