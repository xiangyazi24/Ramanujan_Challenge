#!/usr/bin/env python3
"""P2.5: Set up the double contour integral carrier for D_n^2 and verify.

Key representation:
  D_n = P_n(3) = [w^n] h(w)^n  where h(w) = (w+2)(w+4)/2

Proof: Via Schlafli's integral P_n(x) = (1/2pi i) oint (z^2-1)^n / [2^n (z-x)^{n+1}] dz
  With x=3, z=3+w: D_n = [w^n] ((8+6w+w^2)/2)^n = [w^n] ((w+2)(w+4)/2)^n

Therefore:
  D_n^2 = [u^n v^n] [h(u)h(v)]^n

Phase: phi(u,v) = h(u)h(v)/(uv) = [(u+2)(u+4)/(2u)] * [(v+2)(v+4)/(2v)]
Critical points: u = ±2sqrt(2), v = ±2sqrt(2)
Critical values: phi(2√2,2√2) = (3+2√2)^2 = 17+12√2  [dominant]
                 phi(-2√2,-2√2) = (3-2√2)^2 = 17-12√2  [recessive]
                 phi(±2√2,∓2√2) = 1  [intermediate]
"""
import mpmath as mp
mp.mp.dps = 100

# Verify the representation D_n = [w^n] h(w)^n
def h(w):
    return (w+2)*(w+4)/2

def D_from_coeff(n, prec=50):
    """Compute D_n via coefficient extraction from h(w)^n."""
    from sympy import symbols, Poly, Rational
    w = symbols('w')
    hw = (w+2)*(w+4)/Rational(1,2)
    p = Poly(hw**n, w)
    coeffs = p.all_coeffs()
    deg = p.degree()
    if n <= deg:
        return coeffs[deg - n]
    return 0

def D_legendre(n):
    """D_n = P_n(3) via Legendre recurrence."""
    if n == 0: return 1
    if n == 1: return 3
    a, b = 1, 3
    for k in range(1, n):
        a, b = b, ((2*k+1)*3*b - k*a) // (k+1)
    return b

print("=== Verifying D_n = [w^n] h(w)^n ===")
for n in range(8):
    d_leg = D_legendre(n)
    d_coeff = D_from_coeff(n)
    print(f"  D_{n} = {d_leg}, coeff = {d_coeff}, match = {d_leg == d_coeff}")

print("\n=== Critical values of the phase ===")
s2 = mp.sqrt(2)
u_plus = 2*s2
u_minus = -2*s2

def phi_u(u):
    return (u+2)*(u+4)/(2*u)

for label, u, v in [("(+,+)", u_plus, u_plus), ("(-,-)", u_minus, u_minus),
                     ("(+,-)", u_plus, u_minus), ("(-,+)", u_minus, u_plus)]:
    val = phi_u(u) * phi_u(v)
    print(f"  phi{label} = {mp.nstr(val, 30)}")

print(f"\n  17+12√2 = {mp.nstr(17+12*s2, 30)}")
print(f"  17-12√2 = {mp.nstr(17-12*s2, 30)}")
print(f"  Ratio dominant/intermed = {mp.nstr((17+12*s2)/1, 30)} (= 17+12√2)")

# CMF Poincare roots (with factor -16)
print("\n=== CMF Poincare roots (×(-16)) ===")
for label, val in [("c_+", -16*(17+12*s2)), ("c_0", -16*1), ("c_-", -16*(17-12*s2))]:
    print(f"  {label} = {mp.nstr(val, 20)}")

# Verify convergence rate
rho = 17 - 12*s2
print(f"\n  Convergence rate rho = c_0/c_+ = 1/(17+12√2) = {mp.nstr(rho, 20)}")
print(f"  Digits per step = -log10(rho) = {mp.nstr(-mp.log10(rho), 10)}")
