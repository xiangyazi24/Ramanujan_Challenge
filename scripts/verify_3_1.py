#!/usr/bin/env python3
"""High-precision verification of Ramanujan Challenge Problem 3.1.
Compute the integral of (log x dy/y - log y dx/x) over the A-polynomial
curve of knot 7_2 and compare with 4*pi^2/85.

Based on ChatGPT analysis (Q4654). Run on uisai2 for higher precision."""
from __future__ import annotations
import mpmath as mp

DPS = 200
mp.mp.dps = DPS

def _coefficients(M):
    M2 = M * M; M4 = M2 * M2; M6 = M4 * M2
    M8 = M4 * M4; M10 = M8 * M2; M12 = M10 * M2
    M14 = M12 * M2; M16 = M14 * M2; M18 = M16 * M2
    M20 = M18 * M2; M22 = M20 * M2
    c4 = M14 - M12 + 3*M4 + 4*M2 - 2
    c3 = -2*M18 + 5*M16 + M14 - 4*M12 + 6*M8 + 5*M6 + 2*M4 - 4*M2 + 1
    c2 = M22 - 4*M20 + 2*M18 + 5*M16 + 6*M14 - 4*M10 + M8 + 5*M6 - 2*M4
    c1 = -2*M22 + 4*M20 + 3*M18 - M10 + M8
    c0 = M22
    return c4, c3, c2, c1, c0

def A(M, L):
    c4, c3, c2, c1, c0 = _coefficients(M)
    return ((((L + c4)*L + c3)*L + c2)*L + c1)*L + c0

def A_L(M, L):
    c4, c3, c2, c1, _ = _coefficients(M)
    return (((5*L + 4*c4)*L + 3*c3)*L + 2*c2)*L + c1

def A_M(M, L):
    d4 = 14*M**13 - 12*M**11 + 12*M**3 + 8*M
    d3 = -36*M**17 + 80*M**15 + 14*M**13 - 48*M**11 + 48*M**7 + 30*M**5 + 8*M**3 - 8*M
    d2 = 22*M**21 - 80*M**19 + 36*M**17 + 80*M**15 + 84*M**13 - 40*M**9 + 8*M**7 + 30*M**5 - 8*M**3
    d1 = -44*M**21 + 80*M**19 + 54*M**17 - 10*M**9 + 8*M**7
    d0 = 22*M**21
    return (((d4*L + d3)*L + d2)*L + d1)*L + d0

def slope(x, y):
    return -A_M(x, y) / A_L(x, y)

def ode_rhs(x, state):
    y, integral = state
    yp = slope(x, y)
    omega = mp.log(x)*yp/y - mp.log(y)/x
    return [yp, omega]

def main():
    root_tol = mp.power(10, -(DPS - 20))
    alpha = mp.findroot(lambda x: A(x, mp.sqrt(x)), mp.mpf("0.349269"),
                        solver="newton", tol=root_tol)
    beta = mp.findroot(lambda x: A(x, x), mp.mpf("0.406813"),
                       solver="newton", tol=root_tol)

    y0 = mp.sqrt(alpha)
    print(f"precision: {DPS} digits")
    print(f"alpha = {mp.nstr(alpha, 40)}")
    print(f"beta  = {mp.nstr(beta, 40)}")
    print(f"A(alpha, sqrt(alpha)) = {mp.nstr(A(alpha, y0), 10)}")
    print(f"A(beta, beta)         = {mp.nstr(A(beta, beta), 10)}")

    ode_tol = mp.power(10, -(DPS - 25))
    solution = mp.odefun(ode_rhs, alpha, (y0, mp.mpf("0")),
                         tol=ode_tol, degree=50, method="taylor")
    y_beta, integral_val = solution(beta)

    target = 4*mp.pi**2 / 85
    print(f"\ny(beta) = {mp.nstr(y_beta, 30)}")
    print(f"y(beta) - beta = {mp.nstr(y_beta - beta, 15)}")
    print(f"\nIntegral      = {mp.nstr(integral_val, 40)}")
    print(f"4*pi^2/85     = {mp.nstr(target, 40)}")
    print(f"Difference    = {mp.nstr(integral_val - target, 15)}")

    if abs(integral_val - target) < mp.power(10, -30):
        print("\nVERIFIED: integral matches 4*pi^2/85!")
    else:
        print(f"\nWARNING: mismatch of {mp.nstr(abs(integral_val - target), 6)}")

if __name__ == "__main__":
    main()
