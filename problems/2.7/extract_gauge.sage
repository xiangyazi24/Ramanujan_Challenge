# Extract the full rational gauge R(n) found by sage_gauge_search.sage
# Re-run the computation and print the full matrix

from sage.all import *

P.<n> = PolynomialRing(QQ)
K = P.fraction_field()

def sh(f, k=1):
    f = K(f)
    return K(f.numerator()(n+k)) / K(f.denominator()(n+k))

def shM(M, k=1):
    return matrix(K, M.nrows(), M.ncols(), [sh(x, k) for x in M.list()])

def ev(f, a):
    f = K(f)
    den = QQ(f.denominator()(a))
    if den == 0:
        raise ZeroDivisionError("pole at n=%s" % a)
    return QQ(f.numerator()(a)) / den

def evM(M, a, base=QQ):
    return matrix(base, M.nrows(), M.ncols(), [base(ev(x, a)) for x in M.list()])

def zero33(base=K):
    return zero_matrix(base, 3, 3)

# P2.7 companion
def A(x):
    return (1024*(2*x+5)^4*(2*x+7)^3*(2*x+9)^3*(946*x^2+6407*x+10860))

def B(x):
    return (128*(2*x+7)^3*(2*x+9)^3*(104060*x^6+1745370*x^5+12145238*x^4+44886481*x^3+92943995*x^2+102256019*x+46709052))

def C(x):
    return (16*(x+3)^4*(2*x+9)^3*(3784*x^5+57792*x^4+351019*x^3+1059230*x^2+1587211*x+944620))

def Dp(x):
    return (x+3)^4*(x+4)^6*(946*x^2+4515*x+5399)

alphaP = K(64*B(n+2)/A(n+2))
betaP  = K(-64^2*C(n+1)/A(n+1))
gammaP = K(64^3*Dp(n)/A(n))

CP = matrix(K, [[alphaP, betaP, gammaP], [1, 0, 0], [0, 1, 0]])

# Zudilin companion
def Qz(x): return 946*x^2 - 731*x + 153
def Pz(x): return (104060*x^6+127710*x^5+12788*x^4-34525*x^3-8482*x^2+3298*x+1071)
def Sz(x): return (3784*x^5-1032*x^4-1925*x^3+853*x^2+328*x-184)
def Rz(x): return 946*x^2+1161*x+368

m = n + 2
denZ = Qz(m)*(2*m+1)*(m+1)^3
alphaZ = K(Pz(m)/denZ)
betaZ  = K(-m*Sz(m)/denZ)
gammaZ = K(Rz(m)*m*(m-1)^3/(2*denZ))

CZ = matrix(K, [[alphaZ, betaZ, gammaZ], [1, 0, 0], [0, 1, 0]])

# Rank-one twist
def Q209(x): return 946*x^2 - 2623*x + 1830
r = K((n+4)^3 / ((n+QQ(5)/2)*(n+QQ(7)/2)*(n+QQ(9)/2)))
CZh = r * CZ

Delta = K((n+1)^3*(n+2)^4 / ((n+3)^3*Q209(n+QQ(83)/22)*Q209(n+3)))

def gauge_residual(R):
    return shM(R)*CZh - CP*R

# Initial states
q0 = QQ(-215040420000)
q1 = QQ(-167282265043404)/905
q2 = QQ(-964185327658080)/6071
p0 = QQ(-612218384750)
p1 = QQ(-9525021973931919)/18100
p2 = QQ(-29561828382772029)/65380

zb = vector(QQ, [163, 7, 1])
z2 = vector(QQ, [QQ(2145)/8, QQ(23)/2, 0])
z3 = vector(QQ, [QQ(3135)/16, QQ(17)/2, 0])
zm = z2 + z3

Z0 = matrix(QQ, 3, 3)
Z0.set_column(0, zb)
Z0.set_column(1, z2)
Z0.set_column(2, z3)

xq = vector(QQ, [64^2*q2, 64*q1, q0])
xp = vector(QQ, [64^2*p2, 64*p1, p0])

# === Run the search that found the gauge ===
def lcm_polys(polys):
    ans = P.one()
    for f in polys:
        ans = lcm(ans, P(f))
    return ans

def matrix_denominator(M):
    return lcm_polys([K(x).denominator() for x in M.list()])

Ddelta = P(Delta.denominator())
Dsing = lcm(matrix_denominator(CP), matrix_denominator(CZh))

candidates = [
    P.one(),
    Ddelta,
    lcm(Ddelta, prod(f for f, e in Dsing.factor())),
    lcm(Ddelta, Dsing),
]

# Use approach C (symbolic propagation) since it's faster
S.<s2, s1, s0> = PolynomialRing(QQ)

svec = vector(S, [s2, s1, s0])
P0 = matrix(S, 3, 3)
P0.set_column(0, vector(S, xq))
P0.set_column(1, svec)
P0.set_column(2, vector(S, xp) - svec)

R0_family = P0 * matrix(S, Z0.inverse())

def propagated_family(N):
    values = [R0_family]
    for k in range(N-1):
        CPk = matrix(S, evM(CP, k))
        CZk = matrix(S, evM(CZh, k))
        values.append(CPk*values[-1]*CZk.inverse())
    return values

def fit_propagated_family(Dcommon, degree, N=None):
    Dcommon = P(Dcommon)
    if N is None:
        N = degree + 8
    vals = propagated_family(N)

    n_unknown = 3 + 9*(degree+1)
    rows = []
    rhs = []

    def idx(i, j, ell):
        return 3 + (3*i+j)*(degree+1) + ell

    for k in range(N):
        if Dcommon(k) == 0:
            continue
        for i in range(3):
            for j in range(3):
                expr = S(QQ(Dcommon(k))*vals[k][i, j])
                row = [QQ(0)]*n_unknown
                row[0] = expr.monomial_coefficient(s2)
                row[1] = expr.monomial_coefficient(s1)
                row[2] = expr.monomial_coefficient(s0)
                for ell in range(degree+1):
                    row[idx(i, j, ell)] = -QQ(k)^ell
                rows.append(row)
                rhs.append(-expr.constant_coefficient())

    M = matrix(QQ, rows)
    b = vector(QQ, rhs)
    if M.rank() != M.augment(b).rank():
        return None

    sol = M.solve_right(b)
    s_found = vector(QQ, sol[:3])

    Nmat = zero_matrix(P, 3, 3)
    for i in range(3):
        for j in range(3):
            Nmat[i, j] = sum(sol[idx(i, j, ell)]*n^ell for ell in range(degree+1))

    Rfit = matrix(K, 3, 3, [K(x/Dcommon) for x in Nmat.list()])

    if any(x != 0 for x in gauge_residual(Rfit).list()):
        return None
    if evM(Rfit, 0)*zb != xq or evM(Rfit, 0)*zm != xp:
        return None

    return Rfit, s_found

# The winning candidate was lcm(Ddelta, Dsing), degree 25
Dcand = lcm(Ddelta, Dsing)
print("Denominator candidate:", Dcand.factor())
print("Degree:", Dcand.degree())

ans = fit_propagated_family(Dcand, 25)
if ans is not None:
    R, s_found = ans
    print("\n=== GAUGE MATRIX R(n) FOUND ===")
    print("\nThird solution split s =", s_found)
    print("\ny^(2) initial state:", evM(R, 0)*z2)
    print("y^(3) initial state:", evM(R, 0)*z3)

    # Verify
    print("\n=== VERIFICATION ===")
    res = gauge_residual(R)
    print("Gauge residual zero?", all(x == 0 for x in res.list()))
    print("R(0)*zb == xq?", evM(R, 0)*zb == xq)
    print("R(0)*zm == xp?", evM(R, 0)*zm == xp)
    print("det(R) != 0?", R.det() != 0)

    det_ratio = K(R.det()/Delta)
    print("det(R)/Delta =", det_ratio)
    print("Period-one?", sh(det_ratio) == det_ratio)

    # Print each entry of R(n)
    print("\n=== R(n) ENTRIES ===")
    for i in range(3):
        for j in range(3):
            entry = R[i,j]
            num = P(K(entry).numerator())
            den = P(K(entry).denominator())
            print(f"\nR[{i},{j}]:")
            print(f"  numerator degree: {num.degree()}")
            print(f"  denominator degree: {den.degree()}")
            print(f"  numerator: {num}")
            # Don't print full denominator (it's shared)

    print(f"\nCommon denominator degree: {Dcand.degree()}")
    print(f"Common denominator factored: {Dcand.factor()}")

    # Print R(0), R(1) to verify
    print("\n=== R(0) ===")
    print(evM(R, 0))
    print("\n=== R(1) ===")
    print(evM(R, 1))

    # Export the full gauge to a file
    with open("gauge_R.txt", "w") as f:
        f.write("# Rational gauge R(n) for P2.7 <-> h-twisted Zudilin\n")
        f.write(f"# Common denominator: {Dcand}\n")
        f.write(f"# Common denominator factored: {Dcand.factor()}\n")
        f.write(f"# det(R)/Delta = {det_ratio}\n\n")
        for i in range(3):
            for j in range(3):
                entry = R[i,j]
                num = P(K(entry).numerator())
                den = P(K(entry).denominator())
                f.write(f"R[{i},{j}] = ({num}) / ({den})\n\n")
        f.write(f"\ny^(2) initial: {evM(R, 0)*z2}\n")
        f.write(f"y^(3) initial: {evM(R, 0)*z3}\n")

    print("\nFull gauge written to gauge_R.txt")
else:
    print("NOT FOUND at this denominator/degree")
