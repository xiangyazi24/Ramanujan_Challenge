#!/usr/bin/env python3
from fractions import Fraction as F
from math import comb


# ================================================================
# Exact univariate polynomials over Q, coefficients in ascending order
# ================================================================

class Poly:
    def __init__(self, coeffs=(0,)):
        c = [F(x) for x in coeffs]
        while len(c) > 1 and c[-1] == 0:
            c.pop()
        self.c = tuple(c)

    @staticmethod
    def coerce(x):
        return x if isinstance(x, Poly) else Poly((x,))

    @property
    def degree(self):
        return len(self.c) - 1

    def __add__(self, other):
        other = Poly.coerce(other)
        m = max(len(self.c), len(other.c))
        return Poly([
            (self.c[i] if i < len(self.c) else 0)
            + (other.c[i] if i < len(other.c) else 0)
            for i in range(m)
        ])

    __radd__ = __add__

    def __neg__(self):
        return Poly([-x for x in self.c])

    def __sub__(self, other):
        return self + (-Poly.coerce(other))

    def __rsub__(self, other):
        return Poly.coerce(other) - self

    def __mul__(self, other):
        other = Poly.coerce(other)
        out = [F(0)] * (len(self.c) + len(other.c) - 1)
        for i, a in enumerate(self.c):
            for j, b in enumerate(other.c):
                out[i + j] += a * b
        return Poly(out)

    __rmul__ = __mul__

    def __pow__(self, e):
        if e < 0:
            raise ValueError("negative polynomial power")
        ans = Poly((1,))
        base = self
        while e:
            if e & 1:
                ans = ans * base
            base = base * base
            e >>= 1
        return ans

    def __call__(self, x):
        x = F(x)
        ans = F(0)
        for a in reversed(self.c):
            ans = ans * x + a
        return ans

    def shift(self, h=1):
        """Return p(n+h), exactly."""
        h = F(h)
        out = [F(0)] * len(self.c)
        for i, a in enumerate(self.c):
            for j in range(i + 1):
                out[j] += a * comb(i, j) * h ** (i - j)
        return Poly(out)

    def __eq__(self, other):
        return self.c == Poly.coerce(other).c

    def __repr__(self):
        return f"Poly({self.c})"


def poly(*ascending_coeffs):
    return Poly(ascending_coeffs)


def det3(A):
    return (
        A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1])
        - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0])
        + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0])
    )


N = Poly((0, 1))

# Printed CMF matrix, exact polynomial form.
m11 = (-2*N-5)*(N+3)**2*poly(6141, 9551, 5548, 1424, 136)
m12 = poly(169011, 369933, 336377, 162698, 44168, 6384, 384)
m13 = -poly(20730, 32690, 19210, 4980, 480)

m21 = (N+2)**2*(N+3)**2*(4*N+10)*poly(879, 1017, 386, 48)
m22 = (N+2)**2*(-poly(47808, 85761, 61184, 21732, 3848, 272))
m23 = (N+2)**2*poly(5640, 6610, 2540, 320)

m31 = (-4*N-10)*(N+2)**2*(N+3)**2*poly(813, 1530, 1037, 302, 32)
m32 = (N+2)**2*poly(46476, 117279, 120256, 64452, 19116, 2984, 192)
m33 = (N+2)**2*(-poly(6240, 12254, 8884, 2912, 408, 16))

Mpoly = [
    [m11, m12, m13],
    [m21, m22, m23],
    [m31, m32, m33],
]

detM = det3(Mpoly)
detM_expected = (
    -8*(N+1)*(N+2)**6*(N+3)**5
    *(2*N+3)**2*(2*N+5)**3*(2*N+7)**4
)
assert detM == detM_expected
print("det(M) factorization verified; degree =", detM.degree)

# det(C) = detC_num / detC_den.
detC_num = (2*N+5)*(N+1)**2
detC_den = (2*N+3)*(N+3)**2

# Proposed slope-only twist.
s0 = -16*(N+1)**7
ratio0_num = s0**3 * detC_num
ratio0_den = detM * detC_den

R0_num = 512*(N+1)**22
R0_den = (
    (N+2)**6*(N+3)**7*(2*N+3)**3
    *(2*N+5)**2*(2*N+7)**4
)
assert ratio0_num * R0_den == ratio0_den * R0_num
print("R0 verified exactly:")
print("  512*(n+1)^22 / ((n+2)^6*(n+3)^7*(2n+3)^3*(2n+5)^2*(2n+7)^4)")

# Shift-orbit sums for R0 after extracting powers of 2.
integer_orbit_sum = 22 - 6 - 7
half_integer_orbit_sum = -3 - 2 - 4
print("orbit totals:", integer_orbit_sum, half_integer_orbit_sum)
assert integer_orbit_sum == 9
assert half_integer_orbit_sum == -9
print("=> R0 is not f(n+1)/f(n) for any nonzero f in Q(n).")

# Corrected twist.
s1 = -2*(N+1)**4*(2*N+3)**3
ratio1_num = s1**3 * detC_num
ratio1_den = detM * detC_den

R1_num = (N+1)**13*(2*N+3)**6
R1_den = (
    (N+2)**6*(N+3)**7*(2*N+5)**2*(2*N+7)**4
)
assert ratio1_num * R1_den == ratio1_den * R1_num

# delta(n) = 1 / delta_den(n).
delta_den = (N+1)**13*(N+2)**7*(2*N+3)**6*(2*N+5)**4
assert R1_num * delta_den.shift(1) == R1_den * delta_den
print("corrected determinant ratio verified as delta(n+1)/delta(n)")


# ================================================================
# Exact matrix evaluation
# ================================================================

def M_at(n):
    return [[entry(n) for entry in row] for row in Mpoly]


def C_at(n):
    n = F(n)
    Q = 35*n*n + 140*n + 131
    a3 = (2*n+3)*(n+3)**2
    a2 = (2*n+5)*Q
    a1 = (2*n+3)*Q
    a0 = (2*n+5)*(n+1)**2
    return [
        [F(0), F(1), F(0)],
        [F(0), F(0), F(1)],
        [a0/a3, -a1/a3, a2/a3],
    ]


def matmul(A, B):
    return [[
        sum(A[i][k] * B[k][j] for k in range(3))
        for j in range(3)
    ] for i in range(3)]


def matsub(A, B):
    return [[A[i][j] - B[i][j] for j in range(3)] for i in range(3)]


def matscale(c, A):
    return [[c*A[i][j] for j in range(3)] for i in range(3)]


def is_zero_matrix(A):
    return all(x == 0 for row in A for x in row)


def det3_num(A):
    return (
        A[0][0]*(A[1][1]*A[2][2] - A[1][2]*A[2][1])
        - A[0][1]*(A[1][0]*A[2][2] - A[1][2]*A[2][0])
        + A[0][2]*(A[1][0]*A[2][1] - A[1][1]*A[2][0])
    )


# ================================================================
# Fraction-only Gaussian elimination and polynomial-map search
# ================================================================

def nullspace_basis(A):
    """Exact nullspace basis over Q, returned as Fraction vectors."""
    if not A:
        return []
    A = [[F(x) for x in row] for row in A]
    m, ncols = len(A), len(A[0])
    pivot_cols = []
    r = 0

    for c in range(ncols):
        pivot = next((i for i in range(r, m) if A[i][c] != 0), None)
        if pivot is None:
            continue
        A[r], A[pivot] = A[pivot], A[r]
        p = A[r][c]
        A[r] = [x / p for x in A[r]]

        for i in range(m):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [A[i][j] - f*A[r][j] for j in range(ncols)]

        pivot_cols.append(c)
        r += 1
        if r == m:
            break

    free_cols = [c for c in range(ncols) if c not in pivot_cols]
    basis = []
    for fc in free_cols:
        v = [F(0)] * ncols
        v[fc] = F(1)
        for row, pc in reversed(list(enumerate(pivot_cols))):
            v[pc] = -sum(A[row][j] * v[j] for j in free_cols)
        basis.append(v)
    return basis


def index(i, j, r, d):
    return ((3*i + j)*(d+1) + r)


def poly_matrix_at(v, d, n):
    n = F(n)
    P = [[F(0) for _ in range(3)] for _ in range(3)]
    for i in range(3):
        for j in range(3):
            P[i][j] = sum(v[index(i, j, r, d)] * n**r for r in range(d+1))
    return P


def build_system(d, s_fun, den_fun=lambda n: F(1), sample_count=None):
    """
    Search U(n)=P(n)/den(n), with every P_ij of degree <= d.
    After clearing scalar denominators:
        den(n) M(n) P(n+1)
        - s(n) den(n+1) P(n) C(n) = 0.
    """
    if sample_count is None:
        sample_count = 2*d + 30

    nvars = 9*(d+1)
    rows = []

    for nn in range(sample_count):
        M = M_at(nn)
        C = C_at(nn)
        s = F(s_fun(nn))
        q0 = F(den_fun(nn))
        q1 = F(den_fun(nn+1))

        for i in range(3):
            for j in range(3):
                row = [F(0)] * nvars

                # q(n) * M(n) * P(n+1)
                for k in range(3):
                    for r in range(d+1):
                        row[index(k, j, r, d)] += q0 * M[i][k] * F(nn+1)**r

                # -s(n) * q(n+1) * P(n) * C(n)
                for k in range(3):
                    for r in range(d+1):
                        row[index(i, k, r, d)] -= s * q1 * F(nn)**r * C[k][j]

                rows.append(row)

    return rows, sample_count


def U_at(v, d, n, den_fun=lambda n: F(1)):
    P = poly_matrix_at(v, d, n)
    q = F(den_fun(n))
    return [[P[i][j]/q for j in range(3)] for i in range(3)]


def residual(v, d, n, s_fun, den_fun=lambda n: F(1)):
    left = matmul(M_at(n), U_at(v, d, n+1, den_fun))
    right = matscale(F(s_fun(n)), matmul(U_at(v, d, n, den_fun), C_at(n)))
    return matsub(left, right)


def verify_candidate(v, d, s_fun, den_fun, start, count=25):
    return all(
        is_zero_matrix(residual(v, d, n, s_fun, den_fun))
        for n in range(start, start + count)
    )


def generically_invertible(v, d):
    # det(P(n)) has degree at most 3d. Testing 3d+1 points is decisive.
    return any(
        det3_num(poly_matrix_at(v, d, n)) != 0
        for n in range(3*d + 1)
    )


def vector_sum(a, b):
    return [x+y for x, y in zip(a, b)]


def search(max_degree, s_fun, den_fun=lambda n: F(1), label="search"):
    print("\n--", label, "--")
    for d in range(max_degree + 1):
        A, training_count = build_system(d, s_fun, den_fun)
        basis = nullspace_basis(A)
        print(f"degree {d}: nullity {len(basis)}")

        # Basis vectors and simple pairwise combinations.
        candidates = list(basis)
        for i in range(len(basis)):
            for j in range(i+1, len(basis)):
                candidates.append(vector_sum(basis[i], basis[j]))

        for v in candidates:
            if not verify_candidate(v, d, s_fun, den_fun, training_count):
                continue
            if generically_invertible(v, d):
                print("FOUND generically invertible candidate at degree", d)
                return d, v
            else:
                print("verified singular lower-rank candidate at degree", d)

    print("no invertible candidate found through degree", max_degree)
    return None


# Scalar multipliers and denominator ansatz.
def s_bad(n):
    return -16*(n+1)**7


def s_good(n):
    return -2*(n+1)**4*(2*n+3)**3


def q_den(n):
    return (n+1)**5*(n+2)**3*(2*n+3)**2*(2*n+5)**2


if __name__ == "__main__":
    print("\nThe determinant orbit test proves that the bad-twist search")
    print("cannot contain an invertible rational or polynomial solution.")

    # Optional: this can still detect singular lower-rank polynomial maps.
    SEARCH_SINGULAR_BAD = False
    if SEARCH_SINGULAR_BAD:
        search(8, s_bad, label="singular maps for slope-only twist")

    # Correct search: rational U=P/q with determinant-compatible twist.
    # Increase the bound if desired; the determinant test no longer rules it out.
    search(6, s_good, den_fun=q_den,
           label="corrected rational coboundary U=P/q")
