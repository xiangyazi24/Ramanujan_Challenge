#!/usr/bin/env python3
from __future__ import annotations

from fractions import Fraction
from math import gcd, isqrt

Q_LIMIT = 200
SERIES_DEG = 4


def primes_upto(limit: int) -> list[int]:
    mark = bytearray(b"\x01") * (limit + 1)
    mark[:2] = b"\x00\x00"
    for p in range(2, isqrt(limit) + 1):
        if mark[p]:
            mark[p * p : limit + 1 : p] = b"\x00" * (
                (limit - p * p) // p + 1
            )
    return [p for p in range(2, limit + 1) if mark[p]]


def P(n: int) -> int:
    return 34 * n**3 + 51 * n**2 + 27 * n + 5


def exact_apery(limit: int) -> list[int]:
    if limit == 0:
        return [1]
    out = [1, 5]
    for n in range(1, limit):
        num = P(n) * out[n] - n**3 * out[n - 1]
        den = (n + 1) ** 3
        assert num % den == 0
        out.append(num // den)
    return out


def exact_companion(limit: int) -> list[Fraction]:
    if limit == 0:
        return [Fraction(0)]
    out = [Fraction(0), Fraction(6)]
    for n in range(1, limit):
        out.append((P(n) * out[n] - n**3 * out[n - 1]) / (n + 1) ** 3)
    return out


def add(a: list[Fraction], b: list[Fraction]) -> list[Fraction]:
    return [a[i] + b[i] for i in range(SERIES_DEG + 1)]


def scale(a: list[Fraction], c: Fraction) -> list[Fraction]:
    return [c * x for x in a]


def mul(a: list[Fraction], b: list[Fraction]) -> list[Fraction]:
    out = [Fraction(0) for _ in range(SERIES_DEG + 1)]
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            if i + j <= SERIES_DEG:
                out[i + j] += x * y
    return out


def div_series(num: list[Fraction], den: list[Fraction]) -> list[Fraction]:
    assert den[0] != 0
    out = [Fraction(0) for _ in range(SERIES_DEG + 1)]
    for k in range(SERIES_DEG + 1):
        rhs = num[k]
        for j in range(1, k + 1):
            rhs -= den[j] * out[k - j]
        out[k] = rhs / den[0]
    return out


def translated_poly(n: int) -> list[Fraction]:
    # P(n+x).
    return [
        Fraction(P(n)),
        Fraction(102 * n**2 + 102 * n + 27),
        Fraction(102 * n + 51),
        Fraction(34),
        Fraction(0),
    ]


def cube_poly(n: int) -> list[Fraction]:
    # (n+x)^3.
    return [
        Fraction(n**3),
        Fraction(3 * n**2),
        Fraction(3 * n),
        Fraction(1),
        Fraction(0),
    ]


def transfer_series(limit: int) -> tuple[list[list[Fraction]], list[list[Fraction]]]:
    zero = [Fraction(0) for _ in range(SERIES_DEG + 1)]
    one = zero.copy()
    one[0] = 1

    f0 = one
    f1 = div_series(translated_poly(0), cube_poly(1))
    g0 = zero
    gnum = zero.copy()
    gnum[3] = -1
    g1 = div_series(gnum, cube_poly(1))

    F = [f0, f1]
    G = [g0, g1]
    for r in range(1, limit):
        den = cube_poly(r + 1)
        nf = add(mul(translated_poly(r), F[r]), scale(mul(cube_poly(r), F[r - 1]), -1))
        ng = add(mul(translated_poly(r), G[r]), scale(mul(cube_poly(r), G[r - 1]), -1))
        F.append(div_series(nf, den))
        G.append(div_series(ng, den))
    return F, G


def mod_fraction(x: Fraction | int, modulus: int) -> int:
    x = Fraction(x)
    assert gcd(x.denominator, modulus) == 1, (x, modulus)
    return x.numerator % modulus * pow(x.denominator % modulus, -1, modulus) % modulus


def qadic_pair(q: int, limit: int, final_precision: int = 5) -> tuple[list[int], list[int]]:
    # Compute b_n and z_n=q^3 a_n modulo q^final_precision.
    # Every seam (n+1) divisible by q consumes exactly three q-adic digits.
    seams = limit // q
    precision = final_precision + 3 * seams
    modulus = q**precision
    q5 = q**final_precision

    bp, bc = 1 % modulus, 5 % modulus
    zp, zc = 0, (6 * q**3) % modulus
    bvals = [bp % q5, bc % q5]
    zvals = [zp % q5, zc % q5]

    for n in range(1, limit):
        d = n + 1
        nb_num = (P(n) * bc - n**3 * bp) % modulus
        nz_num = (P(n) * zc - n**3 * zp) % modulus
        if d % q:
            inv = pow(pow(d, 3, modulus), -1, modulus)
            nb = nb_num * inv % modulus
            nz = nz_num * inv % modulus
        else:
            t = d // q
            assert nb_num % q**3 == 0, ("b seam", q, d, precision)
            assert nz_num % q**3 == 0, ("z seam", q, d, precision)
            precision -= 3
            newmod = q**precision
            inv = pow(pow(t, 3, newmod), -1, newmod)
            nb = (nb_num // q**3) % newmod * inv % newmod
            nz = (nz_num // q**3) % newmod * inv % newmod
            modulus = newmod
            bp %= modulus
            bc %= modulus
            zp %= modulus
            zc %= modulus
        bp, bc = bc % modulus, nb
        zp, zc = zc % modulus, nz
        bvals.append(bc % q5)
        zvals.append(zc % q5)
    return bvals, zvals


def carry_digit(x: Fraction, digit0: int, q: int) -> int:
    return mod_fraction((x - digit0) / q, q)


def endpoint_digits_direct(
    q: int,
    a: int,
    bsmall: list[int],
    asmall: list[Fraction],
    bmod: list[int],
    zmod: list[int],
) -> tuple[int, int, int]:
    q5 = q**5
    aa = asmall[a]
    aa5 = mod_fraction(aa, q5)
    ba = bsmall[a]
    x = a * q

    eplus = (ba * zmod[x] - aa5 * bmod[x]) % q5
    assert eplus % q**2 == 0, ("D endpoint q2", q, a, eplus)
    p = eplus // q**2 % q
    m = eplus // q**3 % q

    q2 = q**2
    aa2 = mod_fraction(aa, q2)
    eminus = (ba * zmod[x - 1] - aa2 * bmod[x - 1]) % q2
    constant = mod_fraction(Fraction(-6, a**3), q2)
    diff = (eminus - constant) % q2
    assert diff % q == 0, ("D predecessor q", q, a, diff)
    lam = diff // q % q
    return p, m, lam


def endpoint_digits_reflected(
    q: int,
    a: int,
    bsmall: list[int],
    asmall: list[Fraction],
    bmod: list[int],
    zmod: list[int],
) -> tuple[int, int, int]:
    q5 = q**5
    aa = asmall[a]
    aa5 = mod_fraction(aa, q5)
    ba = bsmall[a]
    c = a + 1
    x = c * q - 1

    eplus = (ba * zmod[x] - aa5 * bmod[x]) % q5
    assert eplus % q**2 == 0, ("R endpoint q2", q, a, eplus)
    p = eplus // q**2 % q
    m = eplus // q**3 % q

    # This next-endpoint correction starts at q^2, so its q^1 digit is zero.
    q2 = q**2
    aa2 = mod_fraction(aa, q2)
    enext = (ba * zmod[c * q] - aa2 * bmod[c * q]) % q2
    constant = mod_fraction(Fraction(6, c**3), q2)
    diff = (enext - constant) % q2
    assert diff % q**2 == 0, ("R next endpoint should start q2", q, a, diff)
    return p, m, 0


def actual_digits(
    q: int,
    a: int,
    n: int,
    bsmall: list[int],
    asmall: list[Fraction],
    bmod: list[int],
    zmod: list[int],
) -> tuple[int, int]:
    q5 = q**5
    aa5 = mod_fraction(asmall[a], q5)
    D = (bsmall[a] * zmod[n] - aa5 * bmod[n]) % q5
    assert D % q**3 == 0, ("target projective q3", q, a, n, D)
    return D // q**3 % q, D // q**4 % q


def predicted_direct(
    q: int,
    a: int,
    r: int,
    p: int,
    m: int,
    lam: int,
    bsmall: list[int],
    asmall: list[Fraction],
    F: list[list[Fraction]],
    G: list[list[Fraction]],
) -> tuple[int, int, int, dict[str, int]]:
    beta = bsmall[r] // q
    W, H, J = F[r][1], F[r][2], G[r][4]
    ar = asmall[r]
    omega = Fraction(beta) + a * W
    X = ar + p * omega
    Y = m * omega + a * a * H * p - Fraction(a**3, 6) * ar * lam - 6 * a * J
    k0 = mod_fraction(X, q)
    old_k1 = mod_fraction(Y, q)
    car = carry_digit(X, k0, q)
    k1 = (old_k1 + car) % q
    data = {
        "p": p,
        "m": m,
        "lambda": lam,
        "omega": mod_fraction(omega, q),
        "W": mod_fraction(W, q),
        "H": mod_fraction(H, q),
        "J": mod_fraction(J, q),
        "a_r": mod_fraction(ar, q),
        "carry": car,
    }
    return k0, old_k1, k1, data


def predicted_reflected(
    q: int,
    a: int,
    r: int,
    p: int,
    m: int,
    lam: int,
    bsmall: list[int],
    asmall: list[Fraction],
    F: list[list[Fraction]],
    G: list[list[Fraction]],
) -> tuple[int, int, int, dict[str, int]]:
    beta = bsmall[r] // q
    c = a + 1
    W, H, J = F[r][1], F[r][2], G[r][4]
    ar = asmall[r]
    omega = Fraction(beta) - c * W
    X = ar + p * omega
    # lam is zero through q^5, retained here to make the normalization explicit.
    Y = m * omega + c * c * H * p + Fraction(c**3, 6) * ar * lam + 6 * c * J
    k0 = mod_fraction(X, q)
    old_k1 = mod_fraction(Y, q)
    car = carry_digit(X, k0, q)
    k1 = (old_k1 + car) % q
    data = {
        "p": p,
        "m": m,
        "lambda": lam,
        "omega": mod_fraction(omega, q),
        "W": mod_fraction(W, q),
        "H": mod_fraction(H, q),
        "J": mod_fraction(J, q),
        "a_r": mod_fraction(ar, q),
        "carry": car,
    }
    return k0, old_k1, k1, data


def main() -> None:
    bsmall = exact_apery(Q_LIMIT)
    asmall = exact_companion(Q_LIMIT)
    F, G = transfer_series(Q_LIMIT)
    for r in range(1, Q_LIMIT):
        assert F[r][0] == bsmall[r]
        assert G[r][3] == -asmall[r] / 6

    total_rows = 0
    direct_rows = 0
    reflected_rows = 0
    lower_zero_outer_rows = 0
    q3_failures = []
    old_q4_failures = []
    corrected_q4_failures = []
    hostile = []
    target_primes = []

    for q in primes_upto(Q_LIMIT):
        if q < 7:
            continue
        folded_targets = [
            r for r in range(1, (q - 1) // 2 + 1) if bsmall[r] % q == 0
        ]
        if not folded_targets:
            continue
        target_primes.append((q, tuple(folded_targets)))
        amax = (q - 1) // 2
        limit = (amax + 1) * q
        bmod, zmod = qadic_pair(q, limit, 5)

        for r in folded_targets:
            for a in range(1, amax + 1):
                if bsmall[a] % q == 0:
                    lower_zero_outer_rows += 2

                # Direct row n=aq+r.
                p, m, lam = endpoint_digits_direct(q, a, bsmall, asmall, bmod, zmod)
                n = a * q + r
                actual = actual_digits(q, a, n, bsmall, asmall, bmod, zmod)
                pred0, old1, pred1, data = predicted_direct(
                    q, a, r, p, m, lam, bsmall, asmall, F, G
                )
                rec = (q, a, r, "D", actual, (pred0, old1), (pred0, pred1), data)
                total_rows += 1
                direct_rows += 1
                if actual[0] != pred0:
                    q3_failures.append(rec)
                if actual[1] != old1:
                    old_q4_failures.append(rec)
                if actual != (pred0, pred1):
                    corrected_q4_failures.append(rec)
                if (q, r) in {(11, 5), (17, 3)} and a <= 3:
                    hostile.append(rec)

                # Reflected row n=(a+1)q-1-r.
                p, m, lam = endpoint_digits_reflected(q, a, bsmall, asmall, bmod, zmod)
                n = (a + 1) * q - 1 - r
                actual = actual_digits(q, a, n, bsmall, asmall, bmod, zmod)
                pred0, old1, pred1, data = predicted_reflected(
                    q, a, r, p, m, lam, bsmall, asmall, F, G
                )
                rec = (q, a, r, "R", actual, (pred0, old1), (pred0, pred1), data)
                total_rows += 1
                reflected_rows += 1
                if actual[0] != pred0:
                    q3_failures.append(rec)
                if actual[1] != old1:
                    old_q4_failures.append(rec)
                if actual != (pred0, pred1):
                    corrected_q4_failures.append(rec)
                if (q, r) in {(11, 5), (17, 3)} and a <= 3:
                    hostile.append(rec)

    print("q_limit=", Q_LIMIT)
    print("target_primes=", target_primes)
    print("total_rows=", total_rows)
    print("direct_rows=", direct_rows)
    print("reflected_rows=", reflected_rows)
    print("lower_zero_outer_rows=", lower_zero_outer_rows)
    print("q3_formula_failures=", len(q3_failures))
    print("first_q3_formula_failure=", q3_failures[0] if q3_failures else None)
    print("old_q4_formula_failures=", len(old_q4_failures))
    print("first_old_q4_formula_failure=", old_q4_failures[0] if old_q4_failures else None)
    print("corrected_q4_formula_failures=", len(corrected_q4_failures))
    print("first_corrected_q4_formula_failure=", corrected_q4_failures[0] if corrected_q4_failures else None)
    print("hostile_rows=")
    for rec in hostile:
        print(rec)

    assert not q3_failures
    assert old_q4_failures
    assert old_q4_failures[0][:4] == (11, 1, 5, "D")
    assert not corrected_q4_failures


if __name__ == "__main__":
    main()
