from sage.all import *


def apery_terms(N):
    """Return A_0,...,A_N exactly."""
    if N == 0:
        return [ZZ(1)]
    A = [ZZ(1), ZZ(5)]
    for n in range(1, N):
        P = 34*n^3 + 51*n^2 + 27*n + 5
        num = P*A[n] - n^3*A[n-1]
        den = (n+1)^3
        assert num % den == 0
        A.append(num // den)
    return A


def newton_coefficients(A):
    c = []
    for k in range(len(A)):
        c.append(sum((-1)^(k-j)*binomial(k,j)*A[j]
                     for j in range(k+1)))
    return c


def F_value(c, J, X):
    if J < 0:
        return ZZ(0)
    return sum(c[k]*binomial(X,k) for k in range(J+1))


N = 300
A = apery_terms(N)
c = newton_coefficients(A)
assert c[:8] == [1, 4, 64, 1240, 27640, 667744,
                 17013976, 450174736]

# Exact folded-support verification.
for n in range(7, N+1):
    J = (n-1)//3
    U = F_value(c, J, n)
    V = F_value(c, J, -n-1)
    for p in prime_range(n//2 + 1, n + 1):
        if p == 5:
            continue
        r = n-p
        j = min(r, p-1-r)
        assert 0 <= j <= J
        assert (A[n] - 5*A[j]) % p == 0
        if p == n-j:
            assert (U-A[j]) % p == 0
            assert (A[n]-5*U) % p == 0
            assert (A[n] % p == 0) == (U % p == 0)
        else:
            assert 2*p == n+1+j
            assert (V-A[j]) % p == 0
            assert (A[n]-5*V) % p == 0
            assert (A[n] % p == 0) == (V % p == 0)

print("folded congruences: PASS")

# Heights.
for n in [100, 300]:
    J = (n-1)//3
    U = abs(F_value(c, J, n))
    V = abs(F_value(c, J, -n-1))
    print(n, RR(log(U)/n), RR(log(V)/n))

lam = 17 + 12*sqrt(RR(2))
gam = lam - 1
alpha_U = log(27*gam/4)/3
alpha_V = log(256*gam/27)/3
print("limiting rates", alpha_U, alpha_V)


# Exact polynomial-recurrence guesser with hold-out verification.
def guess_recurrence(seq, max_order=16, max_degree=30):
    R = PolynomialRing(QQ, 'm')
    m = R.gen()
    M = len(seq)

    for order in range(1, max_order+1):
        for degree in range(max_degree+1):
            cols = (order+1)*(degree+1)
            train = min(M-order, 2*cols + 5)
            if train <= cols:
                continue

            rows = []
            for n in range(train):
                rows.append([
                    QQ(seq[n+i]) * QQ(n)^j
                    for i in range(order+1)
                    for j in range(degree+1)
                ])
            K = matrix(QQ, rows).right_kernel()
            if K.dimension() == 0:
                continue

            for vec in K.basis():
                den = lcm([x.denominator() for x in vec])
                iv = [ZZ(x*den) for x in vec]
                g = gcd(iv)
                iv = [x//g for x in iv]
                if next(x for x in reversed(iv) if x) < 0:
                    iv = [-x for x in iv]

                pol = []
                for i in range(order+1):
                    pol.append(sum(QQ(iv[i*(degree+1)+j])*m^j
                                   for j in range(degree+1)))

                ok = True
                for n in range(M-order):
                    if sum(pol[i](n)*seq[n+i]
                           for i in range(order+1)) != 0:
                        ok = False
                        break
                if ok:
                    return pol
    return None


# Build the three sections.
def U_section(s, M):
    out = []
    for m in range(M):
        n = 3*m+s
        J = m-1 if s == 0 else m
        out.append(F_value(c, J, n))
    return out


def V_section(s, M):
    out = []
    for m in range(M):
        n = 3*m+s
        J = m-1 if s == 0 else m
        out.append(F_value(c, J, -n-1))
    return out

# Generate more A,c terms before a serious recurrence search.
# Then, for example:
# polU = guess_recurrence(U_section(1, 120), 16, 30)
# polV = guess_recurrence(V_section(1, 120), 16, 30)
# print(polU)
# print(polV)
