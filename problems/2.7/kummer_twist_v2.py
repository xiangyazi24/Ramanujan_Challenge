"""
Deeper analysis of the gauge between P2.7 and twisted Zudilin.
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

q = [Q(-215040420000), Q(-167282265043404, 905), Q(-964185327658080, 6071)]
for i in range(3, N):
    n = i - 1
    new_q = Q(B_c(n), A_c(n)) * q[-1] - Q(C_c(n-1), A_c(n-1)) * q[-2] + Q(D_c(n-2), A_c(n-2)) * q[-3]
    q.append(new_q)

p = [Q(-612218384750), Q(-9525021973931919, 18100), Q(-29561828382772029, 65380)]
for i in range(3, N):
    n = i - 1
    new_p = Q(B_c(n), A_c(n)) * p[-1] - Q(C_c(n-1), A_c(n-1)) * p[-2] + Q(D_c(n-2), A_c(n-2)) * p[-3]
    p.append(new_p)

# Ratio q_n / ĥ_n
print("=== q_n / ĥ_n ===")
for n in range(20):
    if h_hat[n] != 0:
        r = q[n] / h_hat[n]
        print(f"  n={n}: {float(r):.6e}")

# Ratio p_n / ĥ_n
print("\n=== p_n / ĥ_n ===")
for n in range(15):
    if h_hat[n] != 0:
        r = p[n] / h_hat[n]
        print(f"  n={n}: {float(r):.6e}")

# Check: does q_n / ĥ_n grow as a polynomial in n?
print("\n=== (q_n / ĥ_n) / n^k for various k ===")
for k in [0, 1, 2, 3, 4]:
    print(f"  k={k}:")
    for n in range(1, 15):
        r = q[n] / h_hat[n] / Q(n)**k
        print(f"    n={n}: {float(r):.6e}")

# Look at successive ratios of (q_n/ĥ_n)
print("\n=== (q_{n+1}/ĥ_{n+1}) / (q_n/ĥ_n) ===")
for n in range(1, 20):
    r0 = q[n] / h_hat[n]
    r1 = q[n+1] / h_hat[n+1]
    print(f"  n={n}: {float(r1/r0):.10f}")

# Maybe the relationship involves C(2n,n): try q_n * C(2n,n)^{-1} vs ĥ_n
print("\n=== q_n / (C(2n,n) * ĥ_n) ===")
for n in range(20):
    if h_hat[n] != 0:
        cn = Q(comb(2*n, n))
        r = q[n] / (cn * h_hat[n])
        print(f"  n={n}: {float(r):.6e}")

# Relationship: A_n^{AESZ} = C(2n,n) * b_n. So C(2n,n)*ĥ_n = C(2n,n)*64^{-n}*(5/2)_n/n!*b_n
# = 64^{-n}*(5/2)_n/n! * A_n^{AESZ}
# Note: C(2n,n) * (5/2)_n / n! = C(2n,n) * Γ(5/2+n)/(Γ(5/2)*n!)
# By duplication: C(2n,n) = 4^n * (1/2)_n / n!
# So C(2n,n) * (5/2)_n / n! = 4^n * (1/2)_n * (5/2)_n / (n!)^2

# Actually, 64^{-n} * C(2n,n) = 64^{-n} * 4^n * (1/2)_n/n! = 4^{-n} * (1/2)_n / n!
# Hmm, let me just check what C(2n,n)*ĥ_n equals
print("\n=== C(2n,n) * ĥ_n = C(2n,n) * 64^{-n} * (5/2)_n/n! * b_n ===")
for n in range(10):
    val = Q(comb(2*n,n)) * h_hat[n]
    print(f"  n={n}: {float(val):.6e}")
    # Compare with AESZ: A_n = C(2n,n)*b_n, twisted: 64^{-n}*(5/2)_n/n!*A_n
    val2 = Q(comb(2*n,n)) * b[n] * Q(1,64**n) * pochhammer(Q(5,2),n) / Q(factorial(n))
    assert val == val2

# Try: q_n = g(n) * A_n^{AESZ,twisted} where A_n^{tw} = 64^{-n}*(5/2)_n/n!*C(2n,n)*b_n
# = C(2n,n) * ĥ_n
print("\n=== q_n / A_n^{AESZ,twisted} ===")
A_tw = [Q(comb(2*n,n)) * h_hat[n] for n in range(N)]
for n in range(20):
    if A_tw[n] != 0:
        r = q[n] / A_tw[n]
        print(f"  n={n}: {float(r):.6e}")

# Successive ratios
print("\n=== (q_{n+1}/A_{n+1}^{tw}) / (q_n/A_n^{tw}) ===")
for n in range(1, 20):
    r0 = q[n] / A_tw[n]
    r1 = q[n+1] / A_tw[n+1]
    print(f"  n={n}: {float(r1/r0):.10f}")
