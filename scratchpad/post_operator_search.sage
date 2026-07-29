"""
Post-operator search based on Q5104 §5,8.
CORRECTED: uses proof.tex's monic convention (each coeff / own A).
"""
from sage.all import *

QQt = PolynomialRing(QQ, 't')
t = QQt.gen()

# P2.7 coefficient polynomials
def AA(n):
    n = QQ(n)
    return (1024*(2*n+5)**4*(2*n+7)**3*(2*n+9)**3
            *(946*n**2 + 6407*n + 10860))

def BB(n):
    n = QQ(n)
    return (128*(2*n+7)**3*(2*n+9)**3
            *(104060*n**6 + 1745370*n**5 + 12145238*n**4
              + 44886481*n**3 + 92943995*n**2
              + 102256019*n + 46709052))

def CC(n):
    n = QQ(n)
    return (16*(n+3)**4*(2*n+9)**3
            *(3784*n**5 + 57792*n**4 + 351019*n**3
              + 1059230*n**2 + 1587211*n + 944620))

def DD(n):
    n = QQ(n)
    return ((n+3)**4*(n+4)**6
            *(946*n**2 + 4515*n + 5399))

# CORRECT forward recurrence (proof.tex eq:rec)
# u_{n+1} = B(n)/A(n) u_n - C(n-1)/A(n-1) u_{n-1} + D(n-2)/A(n-2) u_{n-2}
def alpha(n): return BB(n) / AA(n)
def beta(n): return -CC(n-1) / AA(n-1)
def gamma(n): return DD(n-2) / AA(n-2)

def target_q(N):
    q = [
        QQ(-215040420000),
        QQ(-167282265043404) / QQ(905),
        QQ(-964185327658080) / QQ(6071),
    ]
    for n in range(2, N - 1):
        q.append(alpha(n)*q[n] + beta(n)*q[n-1] + gamma(n)*q[n-2])
    return q[:N]

# Harmonic sums
def H(k, r):
    return sum(QQ(1) / j**r for j in range(1, k + 1))

# AESZ #209 inner polynomial
def aesz_poly(n):
    return sum(
        binomial(n, k)**2
        * binomial(n + k, n)
        * binomial(n + 2*k, n)
        * t**k
        for k in range(n + 1)
    )

# Cauchy twist coefficients (tau=1)
def cauchy_twist_coeffs(N, tau=QQ(1)):
    c = [QQ(1)]
    for r in range(N - 1):
        c1 = c[r - 1] if r >= 1 else QQ(0)
        c2 = c[r - 2] if r >= 2 else QQ(0)
        num = ((440*r + 660*tau) * c[r]
               + (144*tau - 16*r + 16) * c1
               + (2*r - 4 - 39*tau) * c2)
        c.append(num / (8*(r + 1)))
    return c

def twisted_polys(N, tau=QQ(1)):
    A_list = [aesz_poly(n) for n in range(N)]
    c = cauchy_twist_coeffs(N, tau)
    return [
        sum(c[n-m] * A_list[m] for m in range(n + 1)) / QQ(64)**n
        for n in range(N)
    ]

# theta_t operator
def theta_t_op(P, j):
    if j == 0:
        return P
    return sum(P[k] * k**j * t**k for k in range(P.degree() + 1))

# Post-operator search
def search_post_operator(Bpolys, qvals,
                         dz=2, dt=2, zdeg=4,
                         nfit=20, nverify=35):
    triples = [(i, j, h)
               for i in range(dz + 1)
               for j in range(dt + 1)
               for h in range(zdeg + 1)]

    print(f"  Search space: {len(triples)} unknowns, {nfit} equations")

    rows = []
    rhs = []

    for n in range(nfit):
        row = []
        for i, j, h in triples:
            if n < h:
                row.append(QQ(0))
            else:
                row.append(QQ(n-h)**i * theta_t_op(Bpolys[n-h], j)(1))
        rows.append(row)
        rhs.append(qvals[n])

    M = matrix(QQ, rows)
    y = vector(QQ, rhs)

    print(f"  Matrix: {M.nrows()} x {M.ncols()}, rank {M.rank()}")

    try:
        sol = M.solve_right(y)
    except ValueError as e:
        print(f"  No solution: {e}")
        return None

    # Check if solution is unique or if there's a kernel
    ker = M.right_kernel()
    print(f"  Kernel dimension: {ker.dimension()}")

    def candidate(n):
        out = QQt.zero()
        for coeff, (i, j, h) in zip(sol, triples):
            if coeff != 0 and n >= h:
                out += coeff * QQ(n-h)**i * theta_t_op(Bpolys[n-h], j)
        return out

    # Verify unused q values
    success = True
    for n in range(nfit, min(nverify, len(qvals))):
        predicted = candidate(n)(1)
        actual = qvals[n]
        if predicted != actual:
            rel = float(abs(predicted - actual) / abs(actual)) if actual != 0 else float('inf')
            print(f"  FAIL at n={n}: rel error = {rel:.6e}")
            success = False
            break

    if success:
        print(f"  PASSED verification for n={nfit}..{min(nverify, len(qvals))-1}")

        # Check if polynomials satisfy the P2.7 recurrence coefficientwise
        print("  Checking recurrence on polynomials...")
        rec_ok = True
        for n in range(2, min(nverify - 1, len(Bpolys) - 1)):
            Qn1 = candidate(n+1)
            Qn = candidate(n)
            Qn_1 = candidate(n-1)
            Qn_2 = candidate(n-2)
            # Use CORRECT monic recurrence
            residual = Qn1 - alpha(n)*Qn - beta(n)*Qn_1 - gamma(n)*Qn_2
            if residual != 0:
                print(f"    Recurrence FAILS at n={n}")
                rec_ok = False
                break
        if rec_ok:
            print("    Recurrence HOLDS coefficientwise!")

    if not success:
        return None
    return triples, sol, candidate


if __name__ == '__main__':
    N = 45
    print("Computing twisted polynomials...")
    Bpolys = twisted_polys(N, tau=1)
    print("Computing target q values...")
    qvals = target_q(N)

    print(f'\nAESZ values:   {[aesz_poly(n)(1) for n in range(5)]}')
    print(f'Twisted at t=1: {[float(Bpolys[n](1)) for n in range(5)]}')
    print(f'Target q:       {[float(qvals[n]) for n in range(5)]}')
    print(f'q ratios:       {[float(qvals[n]/qvals[0]) for n in range(5)]}')

    # Search with increasing complexity
    for dz in range(3):
        for dt in range(3):
            for zdeg in range(6):
                nunk = (dz+1)*(dt+1)*(zdeg+1)
                if nunk < 3 or nunk > 60:
                    continue
                print(f"\n--- dz={dz}, dt={dt}, zdeg={zdeg} ({nunk} unknowns) ---")
                ans = search_post_operator(Bpolys, qvals,
                                           dz=dz, dt=dt, zdeg=zdeg,
                                           nfit=max(nunk+5, 20),
                                           nverify=min(N-2, 40))
                if ans is not None:
                    triples, sol, candidate = ans
                    print("\n*** POST-OPERATOR FOUND! ***")
                    for coeff, (i, j, h) in zip(sol, triples):
                        if coeff != 0:
                            print(f"  u[{i},{j},{h}] = {coeff}")
                    break
            else:
                continue
            break
        else:
            continue
        break
    else:
        print("\nNo post-operator found in the search space.")
