from math import comb, gcd, log, prod
from sympy import primerange


def apery_table(N: int) -> list[int]:
    A = [0] * (N + 1)
    A[0] = 1
    if N >= 1:
        A[1] = 5
    for m in range(1, N):
        num = (34*m**3 + 51*m**2 + 27*m + 5) * A[m] \
              - m**3 * A[m-1]
        den = (m + 1)**3
        assert num % den == 0
        A[m+1] = num // den
    return A


def central(n: int) -> int:
    return comb(n, n // 2)


def top_bad_radical(n: int, An: int) -> int:
    return prod(p for p in primerange(n // 2 + 1, n + 1)
                if An % p == 0)


N = 3000
A = apery_table(N)

# Exact support comparison; n=9 is the finite midpoint exception.
for n in range(10, N + 1):
    g = gcd(A[n], central(n))
    lhs = top_bad_radical(n, A[n])
    rhs = prod(p for p in primerange(n // 2 + 1, n + 1)
               if g % p == 0)
    assert lhs == rhs

# Full-gcd dyadic maxima, as in the question.
for lo, hi in [(80,160), (160,320), (320,640),
               (640,1280), (1280,2560)]:
    val, arg = max((log(gcd(A[n], central(n))) / n, n)
                   for n in range(lo, min(hi, N + 1)))
    print((lo, hi), arg, val)

# The statistic actually needed for q=1.
for lo, hi in [(80,160), (320,640), (640,1280), (1280,2560)]:
    val, arg = max((log(top_bad_radical(n, A[n])) / n, n)
                   for n in range(lo, min(hi, N + 1)))
    print("top", (lo, hi), arg, val)

from math import comb


def D(n, k):
    if k < 0 or k > n:
        return 0
    return comb(n, k) * comb(n + k, k)


def franel(k):
    return sum(comb(k, j)**3 for j in range(k + 1))


def apery_direct(n):
    return sum(D(n, k)**2 for k in range(n + 1))


def apery_strehl(n):
    return sum(D(n, k) * franel(k) for k in range(n + 1))


for n in range(30):
    assert apery_direct(n) == apery_strehl(n)

for p in list(primerange(7, 100)):
    for r in range(p):
        n = p + r
        for s in range(r + 1):
            assert (D(n, s) - D(r, s)) % p == 0
            assert (D(n, p+s) - 2*D(r, s)) % p == 0
            assert (franel(p+s) - 2*franel(s)) % p == 0
        assert (apery_direct(n) - 5*apery_direct(r)) % p == 0
        assert (apery_strehl(n) - 5*apery_strehl(r)) % p == 0
