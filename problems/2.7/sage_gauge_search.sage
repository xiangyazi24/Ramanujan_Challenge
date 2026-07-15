from sage.all import *

# ============================================================
# 0. Rational-function infrastructure
# ============================================================
P.<n> = PolynomialRing(QQ)
K = P.fraction_field()


def sh(f, k=1):
    """Apply n -> n+k to f in QQ(n)."""
    f = K(f)
    return K(f.numerator()(n+k)) / K(f.denominator()(n+k))


def shM(M, k=1):
    return matrix(K, M.nrows(), M.ncols(),
                  [sh(x, k) for x in M.list()])


def ev(f, a):
    f = K(f)
    den = QQ(f.denominator()(a))
    if den == 0:
        raise ZeroDivisionError("pole at n=%s" % a)
    return QQ(f.numerator()(a)) / den


def evM(M, a, base=QQ):
    return matrix(base, M.nrows(), M.ncols(),
                  [base(ev(x, a)) for x in M.list()])


def zero33(base=K):
    return zero_matrix(base, 3, 3)


# ============================================================
# 1. P2.7 companion matrix for qhat_n = 64^n q_n
# ============================================================
def A(x):
    return (1024*(2*x+5)^4*(2*x+7)^3*(2*x+9)^3
            *(946*x^2+6407*x+10860))


def B(x):
    return (128*(2*x+7)^3*(2*x+9)^3
            *(104060*x^6+1745370*x^5+12145238*x^4
              +44886481*x^3+92943995*x^2
              +102256019*x+46709052))


def C(x):
    return (16*(x+3)^4*(2*x+9)^3
            *(3784*x^5+57792*x^4+351019*x^3
              +1059230*x^2+1587211*x+944620))


def Dp(x):
    return (x+3)^4*(x+4)^6*(946*x^2+4515*x+5399)


alphaP = K(64*B(n+2)/A(n+2))
betaP  = K(-64^2*C(n+1)/A(n+1))
gammaP = K(64^3*Dp(n)/A(n))

CP = matrix(K, [
    [alphaP, betaP, gammaP],
    [1,      0,     0],
    [0,      1,     0],
])


# ============================================================
# 2. Zudilin companion matrix
# ============================================================
def Qz(x):
    return 946*x^2 - 731*x + 153


def Pz(x):
    return (104060*x^6+127710*x^5+12788*x^4-34525*x^3
            -8482*x^2+3298*x+1071)


def Sz(x):
    return (3784*x^5-1032*x^4-1925*x^3+853*x^2
            +328*x-184)


def Rz(x):
    return 946*x^2+1161*x+368


m = n + 2
denZ = Qz(m)*(2*m+1)*(m+1)^3
alphaZ = K(Pz(m)/denZ)
betaZ  = K(-m*Sz(m)/denZ)
gammaZ = K(Rz(m)*m*(m-1)^3/(2*denZ))

CZ = matrix(K, [
    [alphaZ, betaZ, gammaZ],
    [1,      0,     0],
    [0,      1,     0],
])


# ============================================================
# 3. Rank-one Kummer twist and determinant identity
# ============================================================
def Q209(x):
    return 946*x^2 - 2623*x + 1830


r = K((n+4)^3 / ((n+QQ(5)/2)*(n+QQ(7)/2)*(n+QQ(9)/2)))

# IMPORTANT: this is the rank-one system twist h_n I_3.
CZh = r * CZ

Delta = K((n+1)^3*(n+2)^4 /
          ((n+3)^3*Q209(n+QQ(83)/22)*Q209(n+3)))

assert CP.det()/CZh.det() == sh(Delta)/Delta


def gauge_residual(R):
    return shM(R)*CZh - CP*R


# Nine-dimensional first-order system for vec(R).
# For column-major vec, vec(A R B) = (B^T tensor A) vec(R).
A9 = CZh.inverse().transpose().tensor_product(CP)


# ============================================================
# 4. Marked initial states
# ============================================================
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
assert Z0.det() == QQ(825)/32

xq = vector(QQ, [64^2*q2, 64*q1, q0])
xp = vector(QQ, [64^2*p2, 64*p1, p0])

# ============================================================
# 5. Fixed-common-denominator rational gauge search
# ============================================================
def lcm_polys(polys):
    ans = P.one()
    for f in polys:
        ans = lcm(ans, P(f))
    return ans


def matrix_denominator(M):
    return lcm_polys([K(x).denominator() for x in M.list()])


def ansatz_basis(Dcommon, degree):
    Dcommon = P(Dcommon)
    basis = []
    for i in range(3):
        for j in range(3):
            for k in range(degree+1):
                E = zero33(K)
                E[i, j] = K(n^k/Dcommon)
                basis.append(E)
    return basis


def exact_gauge_basis(Dcommon, degree):
    """
    Return the full QQ-vector space of gauges R=N/Dcommon with
    deg N_ij <= degree satisfying sh(R)*CZh = CP*R.
    """
    basis = ansatz_basis(Dcommon, degree)
    residuals = [gauge_residual(E) for E in basis]

    # Bound the numerator degree of every possible residual entry.
    maxdeg = 0
    forbidden = P.one()
    for i in range(3):
        for j in range(3):
            dij = lcm_polys([R[i, j].denominator()
                              for R in residuals])
            forbidden = lcm(forbidden, dij)
            for R in residuals:
                f = K(R[i, j]*dij)
                assert f.denominator() == 1
                if f != 0:
                    maxdeg = max(maxdeg, P(f.numerator()).degree())

    # More than maxdeg regular points force each residual numerator to zero.
    points = []
    a = 0
    while len(points) < maxdeg + 1:
        if forbidden(a) != 0:
            points.append(a)
        a += 1

    rows = []
    for a in points:
        for i in range(3):
            for j in range(3):
                rows.append([ev(R[i, j], a) for R in residuals])

    M = matrix(QQ, rows)
    ker = M.right_kernel()

    gauges = []
    for v in ker.basis():
        G = zero33(K)
        for c, E in zip(v, basis):
            G += c*E
        assert all(x == 0 for x in gauge_residual(G).list())
        gauges.append(G)
    return gauges


def marked_affine_space(gauges):
    """
    Impose R(0) zb = xq and R(0) (z2+z3) = xp.
    Returns a particular marked gauge and homogeneous marked directions.
    """
    rdim = len(gauges)
    if rdim == 0:
        return None

    cols = []
    for G in gauges:
        G0 = evM(G, 0)
        cols.append(vector(QQ, list(G0*zb) + list(G0*zm)))
    M = matrix(QQ, 6, rdim)
    for j, col in enumerate(cols):
        M.set_column(j, col)

    rhs = vector(QQ, list(xq) + list(xp))
    if M.rank() != M.augment(rhs).rank():
        return None

    cpart = M.solve_right(rhs)
    cker = M.right_kernel().basis()

    Rpart = zero33(K)
    for c, G in zip(cpart, gauges):
        Rpart += c*G

    Rdirs = []
    for v in cker:
        H = zero33(K)
        for c, G in zip(v, gauges):
            H += c*G
        Rdirs.append(H)

    assert all(x == 0 for x in gauge_residual(Rpart).list())
    assert evM(Rpart, 0)*zb == xq
    assert evM(Rpart, 0)*zm == xp
    return Rpart, Rdirs


def verify_candidate(R):
    assert all(x == 0 for x in gauge_residual(R).list())
    assert evM(R, 0)*zb == xq
    assert evM(R, 0)*zm == xp
    assert R.det() != 0

    det_ratio = K(R.det()/Delta)
    assert sh(det_ratio) == det_ratio
    # A period-one element of QQ(n) is constant.
    print("det(R)/Delta =", det_ratio)

    y2_initial = evM(R, 0)*z2
    y3_initial = evM(R, 0)*z3
    assert y2_initial + y3_initial == xp
    print("X_0(y^(2)) =", y2_initial)
    print("X_0(y^(3)) =", y3_initial)


# First denominator candidates. Failure here is not a proof of nonexistence.
Ddelta = P(Delta.denominator())
Dsing = lcm(matrix_denominator(CP), matrix_denominator(CZh))

candidates = [
    P.one(),
    Ddelta,
    lcm(Ddelta, prod(f for f, e in Dsing.factor())),
    lcm(Ddelta, Dsing),
]

for Dcand in candidates:
    for degree in range(0, min(40, Dcand.degree()+13)):
        Hbasis = exact_gauge_basis(Dcand, degree)
        marked = marked_affine_space(Hbasis)
        if marked is None:
            continue
        Rpart, Rdirs = marked
        if Rpart.det() != 0:
            print("FOUND", Dcand.factor(), degree)
            verify_candidate(Rpart)
            raise SystemExit
        # If Rpart is singular but marked directions remain, test exact
        # small combinations or solve det(Rpart + sum t_i Rdirs[i])) != 0.

# ============================================================
# 6. Exact symbolic propagation of the unknown target split
# ============================================================
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
    """
    Find s2,s1,s0 and a matrix N(n) of degree <= degree such that
    propagated R(k;s) = N(k)/Dcommon for the sampled exact integers k.
    The returned candidate is then verified symbolically.
    """
    Dcommon = P(Dcommon)
    if N is None:
        N = degree + 8
    vals = propagated_family(N)

    # Unknown vector:
    # [s2,s1,s0, c_(0,0,0),...,c_(2,2,degree)]
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
    print("free dimension after fitting =", M.right_kernel().dimension())
    s_found = vector(QQ, sol[:3])

    Nmat = zero_matrix(P, 3, 3)
    for i in range(3):
        for j in range(3):
            Nmat[i, j] = sum(sol[idx(i, j, ell)]*n^ell
                              for ell in range(degree+1))

    Rfit = matrix(K, 3, 3,
                  [K(x/Dcommon) for x in Nmat.list()])

    # Definitive checks: finite interpolation alone is not accepted.
    if any(x != 0 for x in gauge_residual(Rfit).list()):
        return None
    if evM(Rfit, 0)*zb != xq or evM(Rfit, 0)*zm != xp:
        return None

    print("correct third state s =", s_found)
    verify_candidate(Rfit)
    return Rfit, s_found


for Dcand in candidates:
    for degree in range(0, min(40, Dcand.degree()+13)):
        ans = fit_propagated_family(Dcand, degree)
        if ans is not None:
            R, s_found = ans
            raise SystemExit
