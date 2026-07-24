from math import comb, gcd
from functools import reduce


def gcd_list(values):
    return reduce(gcd, (abs(v) for v in values), 0)


def franel(k):
    return sum(comb(k, a) ** 3 for a in range(k + 1))


def L(n, k):
    if k < 0 or k > n:
        return 0
    return comb(n, k) * comb(n + k, k)


def coefficient_family(n):
    J = (n - 1) // 3
    m = J + 1
    F = [franel(i) for i in range(J + 1)]

    C = [0] * (n + 1)
    C[0] = sum(L(n, i) * F[i] for i in range(J + 1))

    for d in range(1, n + 1):
        lo = max(0, J - d + 1)
        hi = min(J, n - d)
        C[d] = sum(
            (-1) ** (J - i)
            * L(n, i + d)
            * comb(i + d, i)
            * comb(d - 1, J - i)
            * F[i]
            for i in range(lo, hi + 1)
        )

    a = [0] * (n + 1)
    a[0] = 1
    for d in range(m, n + 1):
        a[d] = (-1) ** J * L(n, d) * comb(d - 1, J)

    Bminus = comb(n, m)
    Bplus = comb(n + m, m)
    return J, m, F, C, a, Bminus, Bplus


def branch_gcds(n):
    J, m, F, C, a, Bminus, Bplus = coefficient_family(n)
    gamma_minus = gcd_list([Bminus] + C)
    gamma_plus = gcd_list([Bplus] + C)
    gamma_all = gcd_list(C)
    return gamma_minus, gamma_plus, gamma_all

assert M * V == block_before_row_reduction

expected = L(n, m) ** J
for d in range(1, J + 1):
    expected *= binomial(m, d)
assert abs(Mlow.det()) == expected

for d in range(1, n + 1):
    for i in range(1, J + 1):
        assert (M[d, i] - a[d] * M[0, i]) % Pi == 0
