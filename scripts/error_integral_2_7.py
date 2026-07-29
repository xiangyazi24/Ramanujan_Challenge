from __future__ import annotations

from fractions import Fraction as Q


def A(n: int) -> int:
    return (
        1024 * (2*n + 5)**4 * (2*n + 7)**3 * (2*n + 9)**3
        * (946*n*n + 6407*n + 10860)
    )


def B(n: int) -> int:
    return 128 * (2*n + 7)**3 * (2*n + 9)**3 * (
        104060*n**6 + 1745370*n**5 + 12145238*n**4
        + 44886481*n**3 + 92943995*n**2
        + 102256019*n + 46709052
    )


def C(n: int) -> int:
    return 16 * (n + 3)**4 * (2*n + 9)**3 * (
        3784*n**5 + 57792*n**4 + 351019*n**3
        + 1059230*n**2 + 1587211*n + 944620
    )


def D(n: int) -> int:
    return (n + 3)**4 * (n + 4)**6 * (
        946*n*n + 4515*n + 5399
    )


def next_term(n: int, u_nm2: Q, u_nm1: Q, u_n: Q) -> Q:
    """Verifier convention: compute u_{n+1} for n>=2."""
    return (
        Q(B(n), A(n)) * u_n
        - Q(C(n - 1), A(n - 1)) * u_nm1
        + Q(D(n - 2), A(n - 2)) * u_nm2
    )


p = [
    Q(-612218384750),
    Q(-9525021973931919, 18100),
    Q(-29561828382772029, 65380),
]
q = [
    Q(-215040420000),
    Q(-167282265043404, 905),
    Q(-964185327658080, 6071),
]

# A polynomial a+b*t is stored as the pair (a,b).
def kernel_poly(pn: Q, qn: Q) -> tuple[Q, Q]:
    return (pn / 2 - qn, -pn / 2)


R = [kernel_poly(p[j], q[j]) for j in range(3)]

assert R[0] == (Q(-91068772375), Q(306109192375))
assert R[1] == (
    Q(-2833731372195759, 36200),
    Q(9525021973931919, 36200),
)
assert R[2] == (
    Q(-114331877231773977, 1699880),
    Q(384303768976036377, 1699880),
)

for n in range(2, 20):
    p.append(next_term(n, p[-3], p[-2], p[-1]))
    q.append(next_term(n, q[-3], q[-2], q[-1]))

    a_next = next_term(n, R[-3][0], R[-2][0], R[-1][0])
    b_next = next_term(n, R[-3][1], R[-2][1], R[-1][1])
    R.append((a_next, b_next))

    # Exact identity R_n(t)=p_n(1-t)/2-q_n.
    assert R[-1] == kernel_poly(p[-1], q[-1])

print("Exact recurrence-defined kernel verified through n=20.")

import mpmath as mp

mp.mp.dps = 100


def I_of_s(n: int, s):
    pn = mp.mpf(p[n].numerator) / p[n].denominator
    qn = mp.mpf(q[n].numerator) / q[n].denominator
    return -qn * mp.zeta(2, s + 1) + pn / (2 * (s + 1)**2)


L = mp.zeta(2) + mp.zeta(3)
for n in range(3):
    mixed = I_of_s(n, 0) - mp.diff(lambda ss: I_of_s(n, ss), 0) / 2
    pn = mp.mpf(p[n].numerator) / p[n].denominator
    qn = mp.mpf(q[n].numerator) / q[n].denominator
    assert mp.almosteq(mixed, pn - L * qn)
