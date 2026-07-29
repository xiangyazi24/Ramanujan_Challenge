from fractions import Fraction
from mpmath import mp


def delannoy_basis(max_n: int):
    """Return D_n and E_n through index max_n+3 exactly."""
    D = [Fraction(1), Fraction(3)]
    E = [Fraction(0), Fraction(1)]

    # At step n, append the value with index n+1.
    for n in range(1, max_n + 3):
        D.append(
            (3 * (2 * n + 1) * D[n] - n * D[n - 1])
            / Fraction(n + 1)
        )
        E.append(
            (3 * (2 * n + 1) * E[n] - n * E[n - 1])
            / Fraction(n + 1)
        )
    return D, E


def sym2_residual(U, n: int):
    h = 35 * n * n + 140 * n + 131
    return (
        (2 * n + 3) * (n + 3) ** 2 * U[n + 3]
        - (2 * n + 5) * h * U[n + 2]
        + (2 * n + 3) * h * U[n + 1]
        - (2 * n + 5) * (n + 1) ** 2 * U[n]
    )


def exact_check():
    D, E = delannoy_basis(10)
    basis = {
        "D^2": [D[n] ** 2 for n in range(len(D))],
        "D*E": [D[n] * E[n] for n in range(len(D))],
        "E^2": [E[n] ** 2 for n in range(len(D))],
    }

    for name, U in basis.items():
        for n in range(11):
            r = sym2_residual(U, n)
            assert r == 0, (name, n, r)
        print(f"{name}: exact recurrence verified for n=0,...,10")

    print("D[0:7] =", D[:7])
    print("E[0:7] =", E[:7])


# High-precision integral checks.
mp.dps = 80
sqrt2 = mp.sqrt(2)
Lambda = 3 + 2 * sqrt2
Lambda_inv = 1 / Lambda


def P_integral(n: int):
    return mp.quad(
        lambda th: (3 + 2 * sqrt2 * mp.cos(th)) ** n,
        [0, mp.pi],
    ) / mp.pi


def Q_integral(n: int):
    # u=tanh(t/2) transforms the infinite interval to [0,1].
    return 2 * mp.quad(
        lambda u: (1 - u * u) ** n
        / (Lambda - Lambda_inv * u * u) ** (n + 1),
        [0, 1],
    )


def integral_family(n: int, s, t):
    P = P_integral(n)
    Q = Q_integral(n)
    return (P + s * Q) * (P + t * Q)


def numerical_integral_check():
    s = mp.mpf("0.3")
    t = mp.mpf("-0.7")
    U = [integral_family(n, s, t) for n in range(14)]

    max_residual = mp.mpf("0")
    for n in range(11):
        h = 35 * n * n + 140 * n + 131
        r = (
            (2 * n + 3) * (n + 3) ** 2 * U[n + 3]
            - (2 * n + 5) * h * U[n + 2]
            + (2 * n + 3) * h * U[n + 1]
            - (2 * n + 5) * (n + 1) ** 2 * U[n]
        )
        max_residual = max(max_residual, abs(r))

    print("max integral-family residual =", mp.nstr(max_residual, 12))
    print("Q_0(3) - log(2)/2 =",
          mp.nstr(Q_integral(0) - mp.log(2) / 2, 12))


# Verification of the mixed period by parameter differentiation.
def mixed_period_check():
    Phi = lambda s: mp.zeta(2, 1 + s)  # Hurwitz zeta
    value = Phi(0) - mp.diff(Phi, 0) / 2
    target = mp.zeta(2) + mp.zeta(3)
    print("mixed period =", mp.nstr(value, 60))
    print("difference    =", mp.nstr(value - target, 12))


if __name__ == "__main__":
    exact_check()
    numerical_integral_check()
    mixed_period_check()

D^2: exact recurrence verified for n=0,...,10
D*E: exact recurrence verified for n=0,...,10
E^2: exact recurrence verified for n=0,...,10
