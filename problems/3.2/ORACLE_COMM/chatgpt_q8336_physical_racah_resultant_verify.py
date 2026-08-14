#!/usr/bin/env python3
"""Q8336: exact verifier for the physical-n Racah/resultant obstruction.

Pure Python standard library.  No floating point and no external packages.

The symbolic proofs are in
  problems/3.2/ORACLE_COMM/chatgpt_q8336_physical_racah_resultant_obstruction.md
This program checks every explicit finite identity/counterexample quoted there
and independently stress-tests the exact integral transport formulas.
"""

from __future__ import annotations

from fractions import Fraction
from math import comb, factorial, gcd


def tau(n: int) -> int:
    return n * (n + 1)


def U(n: int, k: int) -> int:
    if k < 0 or k > n:
        return 0
    return comb(n, k) * comb(n + k, k)


def apery(n: int) -> int:
    return sum(U(n, k) ** 2 for k in range(n + 1))


def franel(n: int) -> int:
    return sum(comb(n, k) ** 3 for k in range(n + 1))


def vp(value: int, p: int) -> int:
    assert value != 0
    value = abs(value)
    out = 0
    while value % p == 0:
        value //= p
        out += 1
    return out


def primes_up_to(limit: int) -> list[int]:
    if limit < 2:
        return []
    mark = bytearray(b"\x01") * (limit + 1)
    mark[0:2] = b"\x00\x00"
    d = 2
    while d * d <= limit:
        if mark[d]:
            start = d * d
            mark[start : limit + 1 : d] = b"\x00" * (((limit - start) // d) + 1)
        d += 1
    return [i for i in range(2, limit + 1) if mark[i]]


def icbrt(n: int) -> int:
    lo, hi = 0, n + 1
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if mid**3 <= n:
            lo = mid
        else:
            hi = mid
    return lo


# ---------------------------------------------------------------------------
# Integer polynomials, ascending coefficient order.
# ---------------------------------------------------------------------------


def trim(poly):
    poly = list(poly)
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    return poly


def poly_add(a, b):
    n = max(len(a), len(b))
    out = [0] * n
    for i in range(n):
        out[i] = (a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)
    return trim(out)


def poly_sub(a, b):
    n = max(len(a), len(b))
    out = [0] * n
    for i in range(n):
        out[i] = (a[i] if i < len(a) else 0) - (b[i] if i < len(b) else 0)
    return trim(out)


def poly_scale(a, scalar):
    return trim([scalar * x for x in a])


def poly_mul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return trim(out)


def poly_eval(poly, x):
    out = 0
    for c in reversed(poly):
        out = out * x + c
    return out


def divide_root(poly, root):
    """Divide an integer polynomial known to vanish at root by Y-root."""
    poly = trim(poly)
    assert poly_eval(poly, root) == 0
    degree = len(poly) - 1
    if degree == 0:
        return [0]
    q = [0] * degree
    q[-1] = poly[-1]
    for k in range(degree - 2, -1, -1):
        q[k] = poly[k + 1] + root * q[k + 1]
    assert poly[0] + root * q[0] == 0
    assert poly_mul(q, [-root, 1]) == poly
    return trim(q)


def scaled_racah(s: int) -> list[int]:
    r"""Return \hat R_s(Y)=(s!)^2 R_s(Y) in Z[Y].

    R_s(Y) = sum_{k=0}^s U(s,k) phi_k(Y),
    phi_k(Y) = prod_{a<k}(Y-a(a+1))/(k!)^2.
    """
    sf = factorial(s)
    out = [0]
    node_product = [1]
    for k in range(s + 1):
        coefficient = U(s, k) * (sf // factorial(k)) ** 2
        out = poly_add(out, poly_scale(node_product, coefficient))
        node_product = poly_mul(node_product, [-tau(k), 1])
    return trim(out)


# ---------------------------------------------------------------------------
# Fraction-polynomial helpers and exact Sylvester resultants.
# ---------------------------------------------------------------------------


def ftrim(poly):
    return trim([Fraction(x) for x in poly])


def fpoly_add(a, b):
    return ftrim(poly_add(a, b))


def fpoly_sub(a, b):
    return ftrim(poly_sub(a, b))


def fpoly_scale(a, scalar):
    return ftrim(poly_scale(a, Fraction(scalar)))


def fpoly_mul(a, b):
    return ftrim(poly_mul(a, b))


def monic_racah(s: int) -> list[Fraction]:
    sf2 = factorial(s) ** 2
    poly = [Fraction(c, sf2) for c in scaled_racah(s)]
    lead = poly[-1]
    return [c / lead for c in poly]


def det_fraction(matrix: list[list[Fraction]]) -> Fraction:
    a = [[Fraction(x) for x in row] for row in matrix]
    n = len(a)
    det = Fraction(1)
    for col in range(n):
        pivot = next((r for r in range(col, n) if a[r][col] != 0), None)
        if pivot is None:
            return Fraction(0)
        if pivot != col:
            a[col], a[pivot] = a[pivot], a[col]
            det = -det
        piv = a[col][col]
        det *= piv
        for r in range(col + 1, n):
            if a[r][col] == 0:
                continue
            factor = a[r][col] / piv
            for c in range(col, n):
                a[r][c] -= factor * a[col][c]
    return det


def sylvester_resultant(f, g) -> Fraction:
    f = ftrim(f)
    g = ftrim(g)
    m = len(f) - 1
    n = len(g) - 1
    assert m >= 0 and n >= 0
    if n == 0:
        return g[0] ** m
    if m == 0:
        return f[0] ** n
    fd = list(reversed(f))
    gd = list(reversed(g))
    size = m + n
    rows = []
    for shift in range(n):
        row = [Fraction(0)] * size
        row[shift : shift + m + 1] = fd
        rows.append(row)
    for shift in range(m):
        row = [Fraction(0)] * size
        row[shift : shift + n + 1] = gd
        rows.append(row)
    return det_fraction(rows)


# ---------------------------------------------------------------------------
# Polynomial arithmetic mod p for the finite Hasse--Franel falsification.
# ---------------------------------------------------------------------------


def mtrim(poly, p):
    out = [x % p for x in poly]
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def mmul(a, b, p):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] = (out[i + j] + x * y) % p
    return mtrim(out, p)


def madd(a, b, p):
    n = max(len(a), len(b))
    out = [0] * n
    for i in range(n):
        out[i] = ((a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)) % p
    return mtrim(out, p)


def mscale(a, c, p):
    return mtrim([(c * x) % p for x in a], p)


def mpow(poly, exponent, p):
    out = [1]
    base = mtrim(poly, p)
    e = exponent
    while e:
        if e & 1:
            out = mmul(out, base, p)
        base = mmul(base, base, p)
        e >>= 1
    return out


def meval(poly, x, p):
    out = 0
    for c in reversed(poly):
        out = (out * x + c) % p
    return out


def hasse_franel_check(p: int) -> None:
    K = [franel(m) % p for m in range(p)]
    lhs = mmul(K, K, p)
    rhs = [0]
    one_minus_8x = [1, -8]
    one_plus_x = [1, 1]
    for m in range(p):
        xpow = [0] * m + [1]
        psi = mmul(xpow, mpow(one_minus_8x, m, p), p)
        psi = mmul(psi, mpow(one_plus_x, p - 1 - m, p), p)
        rhs = madd(rhs, mscale(psi, apery(m), p), p)
    assert mtrim(lhs, p) == mtrim(rhs, p)


# ---------------------------------------------------------------------------
# Exact theorem-regression checks.
# ---------------------------------------------------------------------------


def check_transport_identities() -> int:
    checks = 0
    for s in range(0, 11):
        poly = scaled_racah(s)
        sf2 = factorial(s) ** 2
        assert poly_eval(poly, tau(s)) == sf2 * apery(s)
        checks += 1

        # Integral remainder theorem around the diagonal node.
        shifted = list(poly)
        shifted[0] -= sf2 * apery(s)
        q = divide_root(shifted, tau(s))
        reconstructed = poly_add([sf2 * apery(s)], poly_mul(q, [-tau(s), 1]))
        assert reconstructed == poly
        checks += 1

        for n in range(s, 31):
            lhs = poly_eval(poly, tau(n))
            kernel = sum(U(s, k) * U(n, k) for k in range(min(s, n) + 1))
            assert lhs == sf2 * kernel
            assert gcd(abs(lhs), abs(tau(n) - tau(s))) == gcd(
                sf2 * apery(s), abs(tau(n) - tau(s))
            )
            # Integer cross-multiplied self-duality.
            other = poly_eval(scaled_racah(n), tau(s))
            assert factorial(n) ** 2 * lhs == sf2 * other
            checks += 3
    return checks


def check_unique_folded_label() -> int:
    checks = 0
    for p in [q for q in primes_up_to(101) if q % 2]:
        for n in range(0, 151):
            roots = [s for s in range((p - 1) // 2 + 1) if 2 * s < p and (tau(n) - tau(s)) % p == 0]
            r = n % p
            expected = min(r, p - 1 - r)
            assert roots == [expected]
            checks += 1
    return checks


def check_monic_recurrence_and_resultants() -> int:
    checks = 0
    P = [monic_racah(s) for s in range(0, 8)]
    for j in range(1, 7):
        A = Fraction(j**3 + (j + 1) ** 3, 2 * (2 * j + 1))
        beta = Fraction(j**6, 4 * (4 * j * j - 1))
        rhs = fpoly_sub(fpoly_mul(P[j], [A, 1]), fpoly_scale(P[j - 1], beta))
        assert ftrim(rhs) == ftrim(P[j + 1])
        checks += 1

    for s in range(1, 7):
        exact = sylvester_resultant(P[s], P[s - 1])
        predicted = Fraction((-1) ** (s * (s - 1) // 2), 1)
        for j in range(1, s):
            predicted *= Fraction(j**6, 4 * (4 * j * j - 1)) ** j
        assert exact == predicted
        checks += 1
    return checks


def counted_T(n: int, s: int, p: int) -> bool:
    H = icbrt(n)
    return (
        H < s
        and 2 * s <= n - 1
        and p**3 > n**2
        and 2 * s < p <= n
        and apery(s) % p == 0
        and (tau(n) - tau(s)) % p == 0
    )


def check_T_events_and_adjacent_units() -> tuple[int, list[tuple[int, int, int]]]:
    checks = 0
    events = []
    ps = primes_up_to(100)
    for n in range(2, 101):
        for s in range(icbrt(n) + 1, (n - 1) // 2 + 1):
            for p in ps:
                if counted_T(n, s, p):
                    events.append((n, s, p))
                    # Transported value has exactly the original p-divisibility.
                    value = poly_eval(scaled_racah(s), tau(n))
                    assert value % p == 0
                    # The adjacent row is a p-unit; this is the local manifestation
                    # of the exact adjacent-resultant product.
                    if s >= 1:
                        adjacent = poly_eval(scaled_racah(s - 1), tau(n))
                        assert adjacent % p != 0
                    checks += 2
    assert (16, 5, 11) in events
    return checks, events


def check_explicit_counterexamples() -> int:
    # Actual Apéry values modulo 11, using reflection around 5.
    row11 = [apery(m) % 11 for m in range(11)]
    assert row11 == [1, 5, 7, 4, 1, 0, 1, 4, 7, 5, 1]
    assert apery(5) == 819005
    assert vp(apery(5), 11) == 1

    # Coefficient-zero != evaluation-zero in the Apéry Hasse polynomial A_11.
    a = tau(5) % 11
    assert a == 8
    A11_at_8 = sum(row11[m] * pow(a, m, 11) for m in range(11)) % 11
    assert A11_at_8 == 5

    # A genuine T_16 incidence, but physical transport gives no automatic p^2.
    assert counted_T(16, 5, 11)
    transported = poly_eval(scaled_racah(5), tau(16))
    assert transported % 121 == 77
    assert vp(transported, 11) == 1

    # The diagonal value is the original Apéry number; (5!)^2 == 1 mod 11^2.
    assert factorial(5) ** 2 % 121 == 1
    assert poly_eval(scaled_racah(5), tau(5)) % 121 == apery(5) % 121 == 77
    return 11


def check_hasse_franel_finite() -> int:
    for p in (5, 7, 11, 13):
        hasse_franel_check(p)
    return 4


def main() -> None:
    transport = check_transport_identities()
    folded = check_unique_folded_label()
    resultants = check_monic_recurrence_and_resultants()
    event_checks, events = check_T_events_and_adjacent_units()
    explicit = check_explicit_counterexamples()
    hasse = check_hasse_franel_finite()

    print(f"transport_checks={transport}")
    print(f"unique_folded_checks={folded}")
    print(f"recurrence_resultant_checks={resultants}")
    print(f"T_event_checks={event_checks}")
    print(f"T_event_count_n_le_100={len(events)}")
    print(f"first_T_events={events[:12]}")
    print(f"explicit_counterexample_checks={explicit}")
    print(f"hasse_franel_finite_checks={hasse}")
    print("Q8336_PHYSICAL_RACAH_OBSTRUCTION=PASS")


if __name__ == "__main__":
    main()
