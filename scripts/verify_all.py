#!/usr/bin/env python3
"""Numerical verification for all Ramanujan Challenge problems."""
from mpmath import mp, mpf, pi, euler, e, catalan, zeta, log, polylog, inf, binomial, fac
from fractions import Fraction

mp.dps = 50  # 50 decimal places

def verify_2_1():
    """Problem 2.1: Polynomial continued fraction for π."""
    print("=== Problem 2.1: PCF for π ===")
    target = 6 / (3 - pi)
    print(f"Target 6/(3-π) = {target}")

    def a(n):
        return -220*n**3 - 484*n**2 - 301*n - 42

    def b(n):
        return 4*n**2 * (2*n+1)**2 * (5*n-4) * (5*n+6)

    # Compute continued fraction from bottom up
    N = 200
    val = mpf(0)
    for n in range(N, 0, -1):
        val = b(n) / (a(n) + val)
    val = a(0) + val
    print(f"PCF (N={N})     = {val}")
    print(f"Difference      = {val - target}")
    print()

def verify_2_2():
    """Problem 2.2: Euler's γ as Apéry limit."""
    print("=== Problem 2.2: γ as Apéry limit ===")
    target = euler
    print(f"Target γ = {target}")

    # 4-term recurrence coefficients
    def c0(n): return -8*n**3 - 51*n**2 - 105*n - 68
    def c1(n): return 24*n**5 + 337*n**4 + 1833*n**3 + 4818*n**2 + 6092*n + 2928
    def c2(n): return -(n+2)*(n+3)*(24*n**5 + 273*n**4 + 1150*n**3 + 2154*n**2 + 1635*n + 268)
    def c3(n): return (n+1)*(n+2)**4*(n+3)*(8*n**3 + 75*n**2 + 231*n + 232)

    # p sequence: p_{-3}=0, p_{-2}=7, p_{-1}=179
    # q sequence: q_{-3}=1, q_{-2}=12, q_{-1}=306
    p = [mpf(0), mpf(7), mpf(179)]
    q = [mpf(1), mpf(12), mpf(306)]

    N = 100
    for n in range(0, N):
        # recurrence: c0(n)*u_n + c1(n)*u_{n-1} + c2(n)*u_{n-2} + c3(n)*u_{n-3} = 0
        # u_n = -(c1(n)*u_{n-1} + c2(n)*u_{n-2} + c3(n)*u_{n-3}) / c0(n)
        pn = -(c1(n)*p[-1] + c2(n)*p[-2] + c3(n)*p[-3]) / c0(n)
        qn = -(c1(n)*q[-1] + c2(n)*q[-2] + c3(n)*q[-3]) / c0(n)
        p.append(pn)
        q.append(qn)

    ratio = p[-1] / q[-1]
    print(f"p_n/q_n (N={N}) = {ratio}")
    print(f"Difference      = {ratio - target}")
    print()

def verify_2_3():
    """Problem 2.3: π + e as Apéry limit."""
    print("=== Problem 2.3: π + e as Apéry limit ===")
    target = pi + e
    print(f"Target π+e = {target}")

    # 5-term recurrence
    def c0(n): return -n**3 + 2*n**2 + 7*n + 3
    def c1(n): return (n+2)*(2*n**4 + n**3 - 26*n**2 - 48*n - 19)
    def c2(n): return (n+2)*(n**6 + 9*n**5 + 8*n**4 - 87*n**3 - 249*n**2 - 234*n - 68)
    def c3(n): return (n+1)**2*(n+2)*(2*n**5 + 3*n**4 - 13*n**3 - 21*n**2 + 4)
    def c4(n): return -n**3*(n+1)**2*(n+2)*(n**3 + n**2 - 8*n - 11)

    # p: p_{-3}=1, p_{-2}=1, p_{-1}=20, p_0=296
    # q: q_{-3}=1, q_{-2}=0, q_{-1}=4,  q_0=48
    p = [mpf(1), mpf(1), mpf(20), mpf(296)]
    q = [mpf(1), mpf(0), mpf(4), mpf(48)]

    N = 80
    for n in range(1, N):
        pn = -(c1(n)*p[-1] + c2(n)*p[-2] + c3(n)*p[-3] + c4(n)*p[-4]) / c0(n)
        qn = -(c1(n)*q[-1] + c2(n)*q[-2] + c3(n)*q[-3] + c4(n)*q[-4]) / c0(n)
        p.append(pn)
        q.append(qn)

    ratio = p[-1] / q[-1]
    print(f"p_n/q_n (N={N}) = {ratio}")
    print(f"Difference      = {ratio - target}")
    print()

def verify_2_4():
    """Problem 2.4: Harmonic/polylogarithm identity."""
    print("=== Problem 2.4: Harmonic + polylog + zeta ===")

    # RHS
    rhs = (20*polylog(4, mpf(1)/2) + mpf(5)/6 * log(2)**4 + 10*zeta(2)
           - mpf(65)/9 * zeta(2)**2 - log(2)**2 * (12 + 5*zeta(2))
           + mpf(1)/2 * zeta(3) + log(2) * (mpf(35)/2 * zeta(3) - 16))
    print(f"RHS = {rhs}")

    # LHS: double sum
    lhs = mpf(0)
    H = [mpf(0)]  # H_0 = 0
    for k in range(1, 201):
        H.append(H[-1] + mpf(1)/k)

    for m in range(0, 150):
        inner = mpf(0)
        for k in range(0, m+1):
            binom_mk = binomial(m, k)
            inner += binom_mk**2 * H[k]**2
        binom_2m_m = binomial(2*m, m)
        lhs += inner / ((m+1)**2 * binom_2m_m)

    print(f"LHS (M=150)     = {lhs}")
    print(f"Difference      = {lhs - rhs}")
    print()

def verify_2_6():
    """Problem 2.6: Series for ζ(2) + ζ(3)."""
    print("=== Problem 2.6: Series for ζ(2)+ζ(3) ===")
    target = zeta(2) + zeta(3)
    print(f"Target ζ(2)+ζ(3) = {target}")

    # u_1 = -93/4480, u_2 = -117/14000
    u = [mpf(0), mpf(-93)/4480, mpf(-117)/14000]

    N = 200
    for n in range(3, N):
        # 0 = -2(n+3)^3(2n+5)(3n+5) u_n + (n+2)^2(15n^3+85n^2+155n+93) u_{n-1}
        #     - (n+1)^3(n+2)(3n+8) u_{n-2}
        c_n = -2*(n+3)**3*(2*n+5)*(3*n+5)
        c_nm1 = (n+2)**2*(15*n**3 + 85*n**2 + 155*n + 93)
        c_nm2 = -(n+1)**3*(n+2)*(3*n+8)
        un = -(c_nm1*u[-1] + c_nm2*u[-2]) / c_n
        u.append(un)

    series_sum = mpf(2077)/720 + sum(u[1:])
    print(f"2077/720 + Σu_j = {series_sum}")
    print(f"Difference       = {series_sum - target}")
    print()

def verify_2_5():
    """Problem 2.5: Catalan's constant G via 3x3 matrix recurrence."""
    print("=== Problem 2.5: Catalan's constant G ===")
    target = catalan
    print(f"Target G = {target}")

    # Matrix entries m_{ij}(n) for the 3x3 matrix M(n)
    def M(n):
        m11 = (-2*n-5)*(n+3)**2*(136*n**4+1424*n**3+5548*n**2+9551*n+6141)
        m12 = 384*n**6+6384*n**5+44168*n**4+162698*n**3+336377*n**2+369933*n+169011
        m13 = -480*n**4-4980*n**3-19210*n**2-32690*n-20730
        m21 = (n+2)**2*(n+3)**2*(4*n+10)*(48*n**3+386*n**2+1017*n+879)
        m22 = (n+2)**2*(-272*n**5-3848*n**4-21732*n**3-61184*n**2-85761*n-47808)
        m23 = (n+2)**2*(320*n**3+2540*n**2+6610*n+5640)
        m31 = (-4*n-10)*(n+2)**2*(n+3)**2*(32*n**4+302*n**3+1037*n**2+1530*n+813)
        m32 = (n+2)**2*(192*n**6+2984*n**5+19116*n**4+64452*n**3+120256*n**2+117279*n+46476)
        m33 = (n+2)**2*(-16*n**5-408*n**4-2912*n**3-8884*n**2-12254*n-6240)
        return [[mpf(m11), mpf(m12), mpf(m13)],
                [mpf(m21), mpf(m22), mpf(m23)],
                [mpf(m31), mpf(m32), mpf(m33)]]

    def mat_mul(A, B):
        n = len(A)
        m = len(B[0])
        k = len(B)
        C = [[mpf(0)]*m for _ in range(n)]
        for i in range(n):
            for j in range(m):
                for l in range(k):
                    C[i][j] += A[i][l] * B[l][j]
        return C

    # A * M(0) * M(1) * ... * M(N-1) = [[P_{N,1}, P_{N,2}, P_{N,3}],
    #                                     [Q_{N,1}, Q_{N,2}, Q_{N,3}]]
    A = [[mpf(30921), mpf(-32972), mpf(8240)],
         [mpf(33750), mpf(-36000), mpf(9000)]]

    # Accumulate product M(0)*M(1)*...*M(N-1)
    prod = [[mpf(1),mpf(0),mpf(0)],[mpf(0),mpf(1),mpf(0)],[mpf(0),mpf(0),mpf(1)]]
    N = 60
    for n in range(N):
        prod = mat_mul(prod, M(n))

    result = mat_mul(A, prod)
    for j in range(3):
        ratio = result[0][j] / result[1][j]
        print(f"P_{{N,{j+1}}}/Q_{{N,{j+1}}} = {ratio}")
    print(f"Difference (j=1) = {result[0][0]/result[1][0] - target}")
    print()

if __name__ == "__main__":
    verify_2_1()
    verify_2_2()
    verify_2_3()
    verify_2_4()
    verify_2_5()
    verify_2_6()
    print("Done.")
