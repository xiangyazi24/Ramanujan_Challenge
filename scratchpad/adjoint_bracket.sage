"""
Adjoint bracket computation for P2.7.
Verifies J(w⁰, p)/J(w⁰, q) = ζ(2)+ζ(3) using Miller backward iteration.
Based on Q5128's exact bracket formula.
"""
from sage.all import *

def A_Q(k):
    k = QQ(k)
    return QQ(1024)*(2*k+5)^4*(2*k+7)^3*(2*k+9)^3*(946*k^2+6407*k+10860)

def B_Q(k):
    k = QQ(k)
    return QQ(128)*(2*k+7)^3*(2*k+9)^3*(104060*k^6+1745370*k^5+12145238*k^4+44886481*k^3+92943995*k^2+102256019*k+46709052)

def C_Q(k):
    k = QQ(k)
    return QQ(16)*(k+3)^4*(2*k+9)^3*(3784*k^5+57792*k^4+351019*k^3+1059230*k^2+1587211*k+944620)

def D_Q(k):
    k = QQ(k)
    return (k+3)^4*(k+4)^6*(946*k^2+4515*k+5399)

def alpha_Q(k): return B_Q(k)/A_Q(k)
def beta_Q(k): return -C_Q(k-1)/A_Q(k-1)
def gamma_Q(k): return D_Q(k-2)/A_Q(k-2)

p_init = [QQ(-612218384750), QQ(-9525021973931919)/QQ(18100), QQ(-29561828382772029)/QQ(65380)]
q_init = [QQ(-215040420000), QQ(-167282265043404)/QQ(905), QQ(-964185327658080)/QQ(6071)]

# Exact bracket matrix
K0 = matrix(QQ, [
    [0, 0, 1],
    [gamma_Q(2), beta_Q(2), 0],
    [0, gamma_Q(3), 0],
])
print("K0 =")
print(K0)
print()
print("gamma_2 = D(0)/A(0) =", gamma_Q(2))
print("beta_2 = -C(1)/A(1) =", beta_Q(2))
print("gamma_3 = D(1)/A(1) =", gamma_Q(3))

# High precision
BITS = 2000
RF = RealField(BITS)

def alpha(k): return RF(alpha_Q(k))
def beta(k): return RF(beta_Q(k))
def gamma(k): return RF(gamma_Q(k))

def forward(init, N):
    u = [RF(0)]*(N+3)
    for j in range(3):
        u[j] = RF(init[j])
    for n in range(2, N+2):
        u[n+1] = alpha(n)*u[n] + beta(n)*u[n-1] + gamma(n)*u[n-2]
    return u

NTEST = 110
p = forward(p_init, NTEST)
q = forward(q_init, NTEST)

# Miller backward iteration for slow adjoint
# Adjoint: w_n = alpha(n+2)*w_{n+1} + beta(n+3)*w_{n+2} + gamma(n+4)*w_{n+3}
def miller(M, seed):
    w = [RF(0)]*(M+3)
    w[M] = RF(seed[0])
    w[M+1] = RF(seed[1])
    w[M+2] = RF(seed[2])
    for n in range(M-1, -1, -1):
        w[n] = alpha(n+2)*w[n+1] + beta(n+3)*w[n+2] + gamma(n+4)*w[n+3]
    return w

def J0(w, u0):
    return gamma(2)*w[1]*RF(u0[0]) + (beta(2)*w[1] + gamma(3)*w[2])*RF(u0[1]) + w[0]*RF(u0[2])

def Jm(w, u, m):
    return gamma(m+2)*w[m+1]*u[m] + (beta(m+2)*w[m+1] + gamma(m+3)*w[m+2])*u[m+1] + w[m]*u[m+2]

def normalized_miller(M, seed):
    w = miller(M, seed)
    scale = J0(w, q_init)
    assert scale != 0
    return [x/scale for x in w]

print("\n=== Miller backward iteration ===")
M = 220
wA = normalized_miller(M, (1, 0, 1))
wB = normalized_miller(M, (0, 1, 1))

seed_error = max(abs(wA[j]-wB[j]) for j in range(3))
print(f"Miller seed disagreement = {float(seed_error):.3e}")

w = wA

# Bracket ratio
Jq = J0(w, q_init)
Jp = J0(w, p_init)
ratio = Jp/Jq

Ltarget = RF.pi()^2/RF(6) + RF(zeta(3).n(prec=BITS))
error = abs(ratio - Ltarget)
matching_digits = -float(log(error, 10)) if error > 0 else float('inf')

print(f"\nJ(w,q) = 1 (by normalization)")
print(f"J(w,p) = {str(ratio)[:80]}...")
print(f"ζ(2)+ζ(3) = {str(Ltarget)[:80]}...")
print(f"absolute error = {float(error):.5e}")
print(f"matching decimal digits ≈ {matching_digits:.0f}")

# Constancy check
print("\n=== Bracket constancy ===")
for m in [0, 1, 2, 5, 10, 30, 60, 90]:
    if m == 0:
        jqm = J0(w, q_init)
        jpm = J0(w, p_init)
    else:
        jqm = Jm(w, q, m)
        jpm = Jm(w, p, m)
    err_q = float(abs(jqm - Jq))
    err_p = float(abs(jpm - Jp))
    print(f"  m={m:3d}: |Jq-Jq0| = {err_q:.3e}, |Jp-Jp0| = {err_p:.3e}")

# Generic directions give different ratios
print("\n=== Generic adjoint directions ===")
for W in [(1,0,0), (0,1,0), (0,0,1), (1,1,1)]:
    w0, w1, w2 = map(RF, W)
    num = w0*RF(p_init[2]) + w1*(gamma(2)*RF(p_init[0]) + beta(2)*RF(p_init[1])) + w2*gamma(3)*RF(p_init[1])
    den = w0*RF(q_init[2]) + w1*(gamma(2)*RF(q_init[0]) + beta(2)*RF(q_init[1])) + w2*gamma(3)*RF(q_init[1])
    if den != 0:
        r = num/den
        print(f"  W={W}: ratio = {float(r):.10f}")

# Characteristic scales
Rx = PolynomialRing(QQ, 'x')
x = Rx.gen()
P27 = 4*x^3 - 220*x^2 + 8*x - 1
mu0 = P27.roots(RF, multiplicities=False)[0]
print(f"\nμ₀ = {float(mu0):.15f}")
print(f"Poincaré root μ₀/64 = {float(mu0/64):.15f}")
print(f"Slow adjoint base 64/μ₀ = {float(RF(64)/mu0):.15f}")

print("\nDone.")
