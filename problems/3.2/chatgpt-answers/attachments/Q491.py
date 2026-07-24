from sage.all import *
from ore_algebra import OreAlgebra

R.<m> = PolynomialRing(QQ)
K = R.fraction_field()


def P(x):
    return 34*x^3 + 51*x^2 + 27*x + 5

# Even Apéry section X_m = A_{2m}
a0 = (2*m+2)^3 * (2*m+1)^3 * P(2*m-1)
a1 = (
    P(2*m+1) * (P(2*m)*P(2*m-1) - (2*m)^6)
    - (2*m+1)^6 * P(2*m-1)
)
a2 = P(2*m+1) * (2*m)^3 * (2*m-1)^3

# Odd Apéry section Y_m = A_{2m+1}
b0 = (2*m+3)^3 * (2*m+2)^3 * P(2*m)
b1 = (
    P(2*m+2) * (P(2*m+1)*P(2*m) - (2*m+1)^6)
    - (2*m+2)^6 * P(2*m)
)
b2 = P(2*m+2) * (2*m+1)^3 * (2*m)^3

# Hypergeometric ratios for B_{2m}, B_{2m+1}
h0 = 2*(2*m+1)/(m+1)
h0prev = 2*(2*m-1)/m
h1 = 2*(2*m+3)/(m+2)
h1prev = 2*(2*m+1)/(m+1)

rho0 = K(a0)*h0 - K(a1) + K(a2)/h0prev
rho1 = K(b0)*h1 - K(b1) + K(b2)/h1prev

R0 = R(2*(m+1)*(2*m-1)*rho0)
R1 = R(2*(m+2)*(2*m+1)*rho1)

assert R0.degree() == 11
assert R1.degree() == 11
assert R0.leading_coefficient() == -80059392
assert R1.leading_coefficient() == -80059392

# Generate exact Apéry numbers.
def apery(N):
    A = [ZZ(1), ZZ(5)]
    for n in range(1, N):
        num = P(ZZ(n))*A[n] - n^3*A[n-1]
        den = (n+1)^3
        assert num % den == 0
        A.append(num // den)
    return A

AA = apery(30)

for j in range(1, 12):
    assert a0(m=j)*AA[2*j+2] - a1(m=j)*AA[2*j] + a2(m=j)*AA[2*j-2] == 0
    assert b0(m=j)*AA[2*j+3] - b1(m=j)*AA[2*j+1] + b2(m=j)*AA[2*j-1] == 0

# Exact binomial recurrences.
for j in range(0, 12):
    E0 = binomial(2*j, j)
    E1 = binomial(2*j+2, j+1)
    assert (j+1)*E1 == 2*(2*j+1)*E0

    O0 = binomial(2*j+1, j)
    O1 = binomial(2*j+3, j+1)
    assert (j+2)*O1 == 2*(2*j+3)*O0

# Ore right gcds.
Ore.<S> = OreAlgebra(K)

# Shift m -> m+1 to put the parity recurrences in forward form.
def sh(f, k=1):
    return K(f(m=m+k))

LA0 = sh(a2) - sh(a1)*S + sh(a0)*S^2
LB0 = -2*(2*m+1) + (m+1)*S

LA1 = sh(b2) - sh(b1)*S + sh(b0)*S^2
LB1 = -2*(2*m+3) + (m+2)*S

assert LA0.gcrd(LB0).order() == 0
assert LA1.gcrd(LB1).order() == 0

print("even compatibility polynomial degree:", R0.degree())
print("odd compatibility polynomial degree:", R1.degree())
print("even/odd leading coefficient:", R0.leading_coefficient())
print("Ore right gcds are units")

for p in prime_range(7, 200):
    n = 2*p - 1
    assert binomial(n, (n-1)//2) % p != 0
    assert AA[n] % p == 5 if n < len(AA) else True
