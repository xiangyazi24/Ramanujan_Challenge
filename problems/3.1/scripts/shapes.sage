## Exact algebraic computation of the four tetrahedron shapes at both endpoints,
## from the (T,U,V,W) twist-knot deformation chart in the project write-up:
##
##   X = M^2,  u = (L + X^3)/(X(L+X)),  r = -(1+sqrt(1+4u^2))/(2u),  tau = 1-r^2
##   T = tau,  U = u,  V = u/X,  W = 1/(1-uX)
##
## alpha endpoint:  M = a^2, L = a,  a = sqrt(alpha) root of the degree-12 f_alpha
## beta  endpoint:  M = L = b,       b root of A(x,x) near 0.406813

R.<z> = PolynomialRing(QQ)

falpha = (z^12 - 3*z^11 + 4*z^10 - 5*z^9 + 6*z^8 - 7*z^7
          + 7*z^6 - 7*z^5 + 6*z^4 - 5*z^3 + 4*z^2 - 3*z + 1)

prec = 200
RRp = RealField(prec)

def shapes(Mval, Lval, label):
    X   = Mval^2
    u   = (Lval + X^3) / (X*(Lval + X))
    disc = 1 + 4*u^2
    sq  = sqrt(RRp(disc))
    r   = -(1 + sq) / (2*u)
    tau = 1 - r^2
    T = tau; U = u; V = u/X; W = 1/(1 - u*X)
    print("=== %s ===" % label)
    print("  X   =", RRp(X))
    print("  u   =", RRp(u))
    print("  r   =", RRp(r))
    print("  T   =", RRp(T))
    print("  U   =", RRp(U))
    print("  V   =", RRp(V))
    print("  W   =", RRp(W))
    return [RRp(T), RRp(U), RRp(V), RRp(W)]

# --- alpha ---
K.<aa> = NumberField(falpha, embedding=RRp("0.590989428670256"))
a = RRp(aa)
Sa = shapes(a^2, a, "alpha endpoint  (M=a^2, L=a),  a = %s" % RRp(a).n(30))
print("  expected  T=-0.15788  U=6.8158  V=55.872  W=5.9329")

# --- beta ---
R2.<x> = PolynomialRing(QQ)
Rb.<M,L> = PolynomialRing(QQ, 2)
A = (L^5
     + L^4*( M^14 - M^12 + 3*M^4 + 4*M^2 - 2 )
     + L^3*( -2*M^18 + 5*M^16 + M^14 - 4*M^12 + 6*M^8 + 5*M^6 + 2*M^4 - 4*M^2 + 1 )
     + L^2*( M^22 - 4*M^20 + 2*M^18 + 5*M^16 + 6*M^14 - 4*M^10 + M^8 + 5*M^6 - 2*M^4 )
     + L  *( -2*M^22 + 4*M^20 + 3*M^18 - M^10 + M^8 )
     + M^22)
Abeta = R2(A.subs({M: x, L: x}))
bfac = None
for f, e in Abeta.factor():
    rr = [t for t in f.roots(RR, multiplicities=False) if abs(t - 0.406813081336790) < 1e-9]
    if rr:
        bfac = f
print()
print("beta minimal polynomial degree:", bfac.degree(), " palindromic:",
      bfac.coefficients(sparse=False) == list(reversed(bfac.coefficients(sparse=False))))
print("  ", bfac)
Kb.<bb> = NumberField(bfac, embedding=RRp("0.406813081336790"))
b = RRp(bb)
print()
Sb = shapes(b, b, "beta endpoint  (M=L=b),  b = %s" % RRp(b).n(30))
print("  expected  T=-0.25829  U=4.3430  V=26.242  W=3.5555")

# --- Rogers dilogarithm check ---
print()
print("=== extended Rogers dilogarithm sum ===")
CCp = ComplexField(prec)
def Rhat(zv):
    zc = CCp(zv)
    return dilog(zc) + log(zc)*log(1-zc)/2 - CCp.pi()^2/6

Salpha = sum(Rhat(t) for t in Sa)
Sbeta  = sum(Rhat(t) for t in Sb)
D = (Sbeta - Salpha).real()
print("  Re[sum Rhat(beta) - sum Rhat(alpha)] =", RRp(D))
pi2 = RRp(CCp.pi()^2)
print("  divided by pi^2 =", RRp(D/pi2))
for num, den in [(-4,85), (4,85), (-16,85), (16,85), (74,15), (-74,15)]:
    print("     vs %4d/%-3d = %s" % (num, den, RRp(num)/den))
