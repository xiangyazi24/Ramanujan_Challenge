#!/usr/bin/env python3
from fractions import Fraction
from math import comb, isqrt


def primes_upto(limit):
    mark = bytearray(b"\x01") * (limit + 1)
    mark[0:2] = b"\x00\x00"
    for p in range(2, isqrt(limit) + 1):
        if mark[p]:
            mark[p*p:limit+1:p] = b"\x00" * (((limit-p*p)//p)+1)
    return [p for p in range(2, limit + 1) if mark[p]]


def vp(x, p):
    if x == 0:
        return 10**9
    e = 0
    while x % p == 0:
        x //= p
        e += 1
    return e


def mod_fraction(x, mod):
    x = Fraction(x)
    d = x.denominator % mod
    if __import__('math').gcd(d, mod) != 1:
        raise ValueError("nonunit denominator")
    return (x.numerator % mod) * pow(d, -1, mod) % mod


def apery_small(n):
    return sum(comb(n, k)**2 * comb(n+k, k)**2 for k in range(n+1))


def dot_b(n):
    ans = Fraction(0)
    for k in range(n + 1):
        T = comb(n, k)**2 * comb(n+k, k)**2
        D = sum(Fraction(1, j) for j in range(n-k+1, n+k+1))
        ans += 2 * T * D
    return ans


def companion_small(n):
    a0, a1 = Fraction(0), Fraction(6)
    if n == 0:
        return a0
    if n == 1:
        return a1
    prev, cur = a0, a1
    for m in range(1, n):
        P = 34*m**3 + 51*m**2 + 27*m + 5
        nxt = Fraction(P*cur - m**3*prev, (m+1)**3)
        prev, cur = cur, nxt
    return cur


def stripped_factorials(N, p, mod):
    val = [0] * (N + 1)
    unit = [1] * (N + 1)
    for i in range(1, N + 1):
        x = i
        e = 0
        while x % p == 0:
            x //= p
            e += 1
        val[i] = val[i-1] + e
        unit[i] = unit[i-1] * (x % mod) % mod
    return val, unit


def qbinom(n, k, p, mod, val, unit):
    if k < 0 or k > n:
        return 10**9, 0
    e = val[n] - val[k] - val[n-k]
    u = unit[n]
    u = u * pow(unit[k], -1, mod) % mod
    u = u * pow(unit[n-k], -1, mod) % mod
    return e, u


def qadic_data(n, p, K=5):
    mod = p**K
    val, unit = stripped_factorials(2*n, p, mod)
    w = [0] * (n + 1)
    c_val = [0] * (n + 1)
    c_unit = [0] * (n + 1)
    for k in range(n + 1):
        e1, u1 = qbinom(n, k, p, mod, val, unit)
        e2, u2 = qbinom(n+k, k, p, mod, val, unit)
        e = e1 + e2
        u = u1 * u2 % mod
        c_val[k] = e
        c_unit[k] = u
        if 2*e < K:
            w[k] = pow(p, 2*e, mod) * (u*u % mod) % mod
    suffix = [0] * (n + 2)
    for k in range(n, -1, -1):
        suffix[k] = (suffix[k+1] + w[k]) % mod
    return mod, w, suffix, c_val, c_unit


def apery_mod(n, p, power=3):
    mod, w, _, _, _ = qadic_data(n, p, max(5, power))
    return sum(w) % (p**power)


def companion_scaled_mod(n, p):
    K = 5
    bigmod, w, suffix, cval, cunit = qadic_data(n, p, K)
    mod = p**3
    bn = sum(w) % mod
    a = n // p
    Hblock = 0
    for t in range(1, a + 1):
        Hblock = (Hblock + pow(t, -3, mod)) % mod
    ans = bn * Hblock % mod
    for m in range(1, n + 1):
        vm = vp(m, p)
        e = cval[m]
        S = suffix[m]
        vS = min(vp(S, p), K)
        exponent = 3 + vS - e - 3*vm
        if exponent >= 3:
            continue
        if exponent < 0:
            raise AssertionError((n, p, m, e, vS, exponent))
        need = 3 - exponent
        local_mod = p**need
        Su = (S // (p**vS)) % local_mod
        mu = m // (p**vm)
        denu = 2 * (cunit[m] % local_mod) * pow(mu, 3, local_mod)
        denu %= local_mod
        termu = Su * pow(denu, -1, local_mod) % local_mod
        term = (p**exponent) * termu
        if m % 2 == 0:
            term = -term
        ans = (ans + term) % mod
    return ans


def digits(q, a, r):
    n = a*q + r
    q3 = q**3
    ba = apery_small(a)
    br = apery_small(r)
    db = dot_b(r)
    aa = companion_small(a)
    A, E = aa.numerator, aa.denominator
    bn = apery_mod(n, q, 3)
    first_b = ba * (mod_fraction(br, q3) + a*q*mod_fraction(db, q3))
    delta_b = (bn - first_b) % q3
    assert delta_b % (q*q) == 0
    theta = (delta_b // (q*q)) % q
    zn = companion_scaled_mod(n, q)
    first_a = mod_fraction(aa * (br + a*q*db), q3)
    delta_a = (zn - first_a) % q3
    assert delta_a % (q*q) == 0
    xi_raw = (delta_a // (q*q)) % q
    xi = E * xi_raw % q
    beta = (br // q) % q
    wedge = (E*ba*xi_raw - A*theta) % q
    ratio = None if beta == 0 else theta * pow(beta, -1, q) % q
    return beta, theta, xi, wedge, ratio


def audit(limit=200):
    first_wedge_failure = None
    rows = []
    for q in primes_upto(limit):
        if q < 7:
            continue
        targets = []
        for r in range(1, q-1):
            br = apery_small(r)
            if br % q == 0:
                targets.append(r)
        for r in targets:
            for a in range(1, (q-1)//2 + 1):
                n = a*q + r
                if not (2*a < q and n < q*q):
                    continue
                rec = (q, a, r) + digits(q, a, r)
                rows.append(rec)
                if rec[-2] != 0 and first_wedge_failure is None:
                    first_wedge_failure = rec
                    print("FIRST_WEDGE_FAILURE", first_wedge_failure, flush=True)
    print("records:", len(rows))
    print("first wedge failure:", first_wedge_failure)
    for rec in rows:
        print(rec)


if __name__ == "__main__":
    print("HAND", digits(11, 1, 5), flush=True)
    assert digits(11, 1, 5) == (7, 3, 8, 0, 2)
    audit(200)
