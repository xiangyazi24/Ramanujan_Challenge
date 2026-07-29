from fractions import Fraction
from math import comb


def franel(n: int) -> int:
    return sum(comb(n, r) ** 3 for r in range(n + 1))


def harmonic(n: int) -> Fraction:
    return sum((Fraction(1, r) for r in range(1, n + 1)), Fraction(0))


def d_coeff(J: int, m: int, F: list[int]) -> Fraction:
    """Coefficient [x^m] D_J."""
    if m <= J:
        return (
            (-1) ** m
            * comb(J, m)
            * F[m]
            * (harmonic(J) - harmonic(J - m))
        )
    return Fraction(
        (-1) ** (J + 1) * F[m],
        m * comb(m - 1, J),
    )


def M_coeff(J: int, n: int, F: list[int]) -> Fraction:
    """Coefficient [x^n] M_J from the universal convolution."""
    ans = Fraction(0)
    for b in range(min(J, n) + 1):
        p_b = (-1) ** b * comb(J, b) * F[b]
        ans += (n - 2 * b) * p_b * d_coeff(J, n - b, F)
    return ans


def tail_formula(J: int, n: int, F: list[int]) -> Fraction:
    assert n > 2 * J
    ans = Fraction(0)
    for b in range(J + 1):
        ans += Fraction(
            (-1) ** (J + 1 + b)
            * comb(J, b)
            * F[b]
            * F[n - b]
            * (n - 2 * b),
            (n - b) * comb(n - b - 1, J),
        )
    return ans


N = 120
F = [franel(n) for n in range(N + 1)]

for J in range(0, 20):
    for n in range(1, N + 1):
        c = M_coeff(J, n, F)
        assert c.denominator == 1
        if n > 2 * J:
            assert c == tail_formula(J, n, F)

print("all exact Wronskian checks passed")
