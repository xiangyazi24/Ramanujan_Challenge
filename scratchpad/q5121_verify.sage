#!/usr/bin/env sage
"""
Q5121 verification: Ore factorization L27 = L16 · g(n),
Mellin dual M† = M16 · E, and attempt to factor M16.
"""
from sage.all import *
from ore_algebra import OreAlgebra

Rn.<n> = PolynomialRing(QQ)
Kn = Rn.fraction_field()
OS.<S> = OreAlgebra(Kn, 'Sn')

Rt.<t> = PolynomialRing(QQ)
Kt = Rt.fraction_field()
OD.<Dt> = OreAlgebra(Kt, 'Dt')
theta = Kt(t)*Dt

Rx.<x> = PolynomialRing(QQ)

def eval_poly_at_op(f, T):
    out = OD.zero()
    for a in reversed(f.list()):
        out = out*T + Kt(a)
    return out

def primitive_coefficients(c):
    den = Rn.one()
    for a in c:
        den = lcm(den, Rn(a.denominator()))
    nums = [Rn(den*a) for a in c]
    gg = nums[0]
    for a in nums[1:]:
        gg = gcd(gg, a)
    nums = [a.quo_rem(gg)[0] for a in nums]
    scale = QQ(1)/QQ(nums[-1].leading_coefficient())
    return [Rn(scale*a) for a in nums]

def mellin_dual(coeffs, scale_t=QQ(1)):
    out = OD.zero()
    for j, pj in enumerate(coeffs):
        out += Kt((scale_t*t)**j)*eval_poly_at_op(Rn(pj), -theta-j-1)
    return out

def shifted_gcd(coeffs):
    fs = [Rx(Rn(coeffs[j])(x-j)) for j in range(len(coeffs))]
    gg = fs[0]
    for f in fs[1:]:
        gg = gcd(gg, f)
    return gg.monic()

# P2.7 recurrence coefficients
def A(z): return 1024*(2*z+5)**4*(2*z+7)**3*(2*z+9)**3*(946*z**2+6407*z+10860)
def B(z): return 128*(2*z+7)**3*(2*z+9)**3*(104060*z**6+1745370*z**5+12145238*z**4+44886481*z**3+92943995*z**2+102256019*z+46709052)
def C(z): return 16*(z+3)**4*(2*z+9)**3*(3784*z**5+57792*z**4+351019*z**3+1059230*z**2+1587211*z+944620)
def D(z): return (z+3)**4*(z+4)**6*(946*z**2+4515*z+5399)

crat = [Kn(-D(n)/A(n)), Kn(C(n+1)/A(n+1)), Kn(-B(n+2)/A(n+2)), Kn(1)]
p = primitive_coefficients(crat)
print(f"Primitive degrees: {[f.degree() for f in p]}")
assert [f.degree() for f in p] == [18,18,18,18]

L27 = sum((Kn(p[j])*S**j for j in range(4)), OS.zero())

# Q209 and scalar gauge g(n)
Q209 = 946*x**2 - 2623*x + 1830
assert Q209(61/22-x) == Q209, "Reflection symmetry failed"
print("Q209 reflection symmetry: OK")

g = Rn(n**2 + QQ(105)/22*n + QQ(5399)/946)
assert Rx(946*g(x)) == Q209(x+QQ(83)/22)
print(f"g(n) = {g}")

# Ore factorization: p_j(n) = r_j(n) * g(n+j)
r = []
for j in range(4):
    divisor = Rn(g(n=n+j))
    qj, rem = p[j].quo_rem(divisor)
    assert rem == 0, f"g(n+{j}) does not divide p_{j}!"
    assert qj.degree() == 16
    r.append(qj)
print("Ore factorization p_j = r_j * g(n+j): VERIFIED")

L16rec = sum((Kn(r[j])*S**j for j in range(4)), OS.zero())
assert L27 == L16rec*Kn(g), "L27 != L16 * g failed"
print("L27 = L16_rec * g(n): VERIFIED")

# No further Euler factor
assert shifted_gcd(r).degree() == 0
print("No further Euler factor in L16: VERIFIED")

# Mellin dual
print("\nBuilding Mellin adjoint M†...")
M = mellin_dual(p)
assert M.order() == 18
print(f"M† order: {M.order()}")

M16 = mellin_dual(r)
E = eval_poly_at_op(g, -theta-1)
E_expected = theta**2 - QQ(61)/22*theta + QQ(915)/473
assert E == E_expected
print(f"E = {E}")

assert M == M16*E, "M† != M16 * E failed"
print("M† = M16 * E: VERIFIED")
print(f"M16 order: {M16.order()}")

# Leading coefficient
chi = t**3 - QQ(55)/64*t**2 + QQ(1)/2048*t - QQ(1)/2**20
assert M16[16] == Kt(t**16*chi)
print(f"M16 leading coeff = t^16 * chi(t): VERIFIED")

# Indicial polynomials
Rrho.<rho> = PolynomialRing(QQ)
ind0_16 = Rrho(r[0](-rho-1))
indinf_16 = Rrho(r[3](rho-4))
print(f"\nM16 indicial at 0: {factor(ind0_16)}")
print(f"M16 indicial at inf: {factor(indinf_16)}")

# Print factored forms of r_j
print("\n--- Degree-16 recurrence coefficients ---")
for j in range(4):
    print(f"r_{j}(n) = {factor(r[j])}")

# Zudilin/AESZ comparison
print("\n--- AESZ inner Mellin dual comparison ---")
def QZ(z): return 946*z**2-731*z+153
def PZ(z): return 104060*z**6+127710*z**5+12788*z**4-34525*z**3-8482*z**2+3298*z+1071
def SZ(z): return 3784*z**5-1032*z**4-1925*z**3+853*z**2+328*z-184

zcoef = [-QZ(n+3)*(n+2)*(n+1)**3, 2*(n+2)*SZ(n+2), -2*PZ(n+2), 2*QZ(n+2)*(2*n+5)*(n+3)**3]
assert shifted_gcd(zcoef).degree() == 0
MZ = mellin_dual(zcoef)
assert MZ.order() == 6
print(f"AESZ inner Mellin: order {MZ.order()}")

# Check Q209 factor
Qop, remZ = M.quo_rem(E)
# Already verified M = M16*E above

_, remZ2 = MZ.quo_rem(E)
if remZ2 == 0:
    print("AESZ inner Mellin has Q209 factor: YES")
else:
    print("AESZ inner Mellin has Q209 factor: NO")

# Attempt M16 factorization
print("\n--- Attempting M16 factorization ---")
print("(This may take a long time for order 16...)")

# First try right_factors of small order
for ord in [1, 2, 3]:
    print(f"Checking order-{ord} right factors...", end=" ", flush=True)
    try:
        rfs = M16.right_factors(order=ord)
        if rfs:
            print(f"FOUND {len(rfs)}!")
            for f in rfs:
                print(f"  {f}")
        else:
            print("none")
    except Exception as e:
        print(f"error: {e}")

# Try the full factor() -- may be very slow
print("\nAttempting full factor(M16)...", flush=True)
import signal
def handler(signum, frame):
    raise TimeoutError("factorization timed out")

signal.signal(signal.SIGALRM, handler)
signal.alarm(300)  # 5 minute timeout
try:
    fac16 = M16.factor()
    print(f"factor(M16) = {fac16}")
except TimeoutError:
    print("Factorization timed out after 5 minutes")
except Exception as e:
    print(f"factor(M16) error: {e}")
finally:
    signal.alarm(0)

print("\nAll verifications completed.")
