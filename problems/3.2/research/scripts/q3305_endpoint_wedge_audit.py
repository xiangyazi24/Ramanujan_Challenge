#!/usr/bin/env python3
from __future__ import annotations

import argparse
from fractions import Fraction
from functools import lru_cache
from math import comb, gcd, isqrt


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


def vp_int(x: int, p: int) -> int:
    if x == 0:
        return 10**9
    e = 0
    while x % p == 0:
        x //= p
        e += 1
    return e


def mod_fraction(x: Fraction | int, mod: int) -> int:
    x = Fraction(x)
    assert gcd(x.denominator, mod) == 1, (x, mod)
    return x.numerator % mod * pow(x.denominator % mod, -1, mod) % mod


def exact_apery(limit: int) -> list[int]:
    if limit == 0:
        return [1]
    b = [1, 5]
    for n in range(1, limit):
        num = P(n) * b[n] - n**3 * b[n - 1]
        den = (n + 1) ** 3
        assert num % den == 0
        b.append(num // den)
    return b


def exact_companion(limit: int) -> list[Fraction]:
    if limit == 0:
        return [Fraction(0)]
    a = [Fraction(0), Fraction(6)]
    for n in range(1, limit):
        a.append((P(n) * a[n] - n**3 * a[n - 1]) / (n + 1) ** 3)
    return a


def exact_jets(limit: int):
    # F_r(x)=b_r+x W_r+x^2 H_r+O(x^3).
    b = [Fraction(1), Fraction(5)]
    W = [Fraction(0), Fraction(12)]
    H = [Fraction(0), Fraction(0)]
    for n in range(1, limit):
        den = (n + 1) ** 3
        pn = P(n)
        dp = 102 * n**2 + 102 * n + 27
        hp = 102 * n + 51
        nb = (pn * b[n] - n**3 * b[n - 1]) / den
        nW = (
            pn * W[n]
            + dp * b[n]
            - n**3 * W[n - 1]
            - 3 * n**2 * b[n - 1]
            - 3 * (n + 1) ** 2 * nb
        ) / den
        nH = (
            pn * H[n]
            + dp * W[n]
            + hp * b[n]
            - n**3 * H[n - 1]
            - 3 * n**2 * W[n - 1]
            - 3 * n * b[n - 1]
            - 3 * (n + 1) ** 2 * nW
            - 3 * (n + 1) * nb
        ) / den
        b.append(nb)
        W.append(nW)
        H.append(nH)
    return b, W, H


@lru_cache(maxsize=None)
def dot_b_sum(n: int) -> Fraction:
    ans = Fraction(0)
    for k in range(n + 1):
        T = comb(n, k) ** 2 * comb(n + k, k) ** 2
        D = sum(Fraction(1, j) for j in range(n - k + 1, n + k + 1))
        ans += 2 * T * D
    return ans


# ---------- Original Q3228 finite-sum engine ----------

def stripped_factorials(N: int, p: int, mod: int):
    val = [0] * (N + 1)
    unit = [1] * (N + 1)
    for i in range(1, N + 1):
        x = i
        e = 0
        while x % p == 0:
            x //= p
            e += 1
        val[i] = val[i - 1] + e
        unit[i] = unit[i - 1] * (x % mod) % mod
    return val, unit


def qbinom(n: int, k: int, p: int, mod: int, val, unit):
    if k < 0 or k > n:
        return 10**9, 0
    e = val[n] - val[k] - val[n - k]
    u = unit[n]
    u = u * pow(unit[k], -1, mod) % mod
    u = u * pow(unit[n - k], -1, mod) % mod
    return e, u


def qadic_data(n: int, p: int, K: int = 5):
    mod = p**K
    val, unit = stripped_factorials(2 * n, p, mod)
    w = [0] * (n + 1)
    c_val = [0] * (n + 1)
    c_unit = [0] * (n + 1)
    for k in range(n + 1):
        e1, u1 = qbinom(n, k, p, mod, val, unit)
        e2, u2 = qbinom(n + k, k, p, mod, val, unit)
        e = e1 + e2
        u = u1 * u2 % mod
        c_val[k] = e
        c_unit[k] = u
        if 2 * e < K:
            w[k] = pow(p, 2 * e, mod) * (u * u % mod) % mod
    suffix = [0] * (n + 2)
    for k in range(n, -1, -1):
        suffix[k] = (suffix[k + 1] + w[k]) % mod
    return mod, w, suffix, c_val, c_unit


def apery_mod_finite(n: int, p: int, power: int = 3) -> int:
    _, w, _, _, _ = qadic_data(n, p, max(5, power))
    return sum(w) % p**power


def companion_scaled_mod_finite(n: int, p: int) -> int:
    K = 5
    bigmod, w, suffix, cval, cunit = qadic_data(n, p, K)
    mod = p**3
    bn = sum(w) % mod
    high = n // p

    Hblock = 0
    for t in range(1, high + 1):
        Hblock = (Hblock + pow(t, -3, mod)) % mod
    ans = bn * Hblock % mod

    for m in range(1, n + 1):
        vm = vp_int(m, p)
        e = cval[m]
        S = suffix[m]
        vS = min(vp_int(S, p), K)
        exponent = 3 + vS - e - 3 * vm
        if exponent >= 3:
            continue
        if exponent < 0:
            raise AssertionError((n, p, m, e, vS, exponent))
        need = 3 - exponent
        local_mod = p**need
        Su = (S // p**vS) % local_mod
        mu = m // p**vm
        denu = 2 * (cunit[m] % local_mod) * pow(mu, 3, local_mod)
        denu %= local_mod
        termu = Su * pow(denu, -1, local_mod) % local_mod
        term = p**exponent * termu
        if m % 2 == 0:
            term = -term
        ans = (ans + term) % mod
    return ans


def finite_digits(q: int, high: int, r: int, b_small, a_small):
    n = high * q + r
    q3 = q**3
    ba = b_small[high]
    br = b_small[r]
    db = dot_b_sum(r)
    aa = a_small[high]
    A, E = aa.numerator, aa.denominator

    bn = apery_mod_finite(n, q, 3)
    first_b = ba * (mod_fraction(br, q3) + high * q * mod_fraction(db, q3))
    delta_b = (bn - first_b) % q3
    assert delta_b % q**2 == 0, ("finite theta division", q, high, r)
    theta = delta_b // q**2 % q

    zn = companion_scaled_mod_finite(n, q)
    first_a = mod_fraction(aa * (br + high * q * db), q3)
    delta_a = (zn - first_a) % q3
    assert delta_a % q**2 == 0, ("finite xi division", q, high, r)
    xi_raw = delta_a // q**2 % q
    xi = E * xi_raw % q

    beta = br // q % q
    wedge = (ba * xi - A * theta) % q
    return beta, theta, xi, wedge


# ---------- Independent recurrence-only q-adic engine ----------

def qadic_recurrence_pair(q: int, limit: int, final_precision: int = 3):
    seams = limit // q
    precision = final_precision + 3 * seams
    mod = q**precision

    b_prev, b_cur = 1 % mod, 5 % mod
    z_prev, z_cur = 0, (6 * q**3) % mod  # z_n=q^3 a_n
    b = [b_prev, b_cur]
    z = [z_prev, z_cur]

    for n in range(1, limit):
        d = n + 1
        num_b = (P(n) * b_cur - n**3 * b_prev) % mod
        num_z = (P(n) * z_cur - n**3 * z_prev) % mod
        if d % q:
            inv = pow(pow(d, 3, mod), -1, mod)
            nb = num_b * inv % mod
            nz = num_z * inv % mod
        else:
            t = d // q
            assert num_b % q**3 == 0, ("b seam", q, d, precision)
            assert num_z % q**3 == 0, ("z seam", q, d, precision)
            new_precision = precision - 3
            new_mod = q**new_precision
            inv_t3 = pow(pow(t, 3, new_mod), -1, new_mod)
            nb = (num_b // q**3) % new_mod * inv_t3 % new_mod
            nz = (num_z // q**3) % new_mod * inv_t3 % new_mod
            precision = new_precision
            mod = new_mod
            b_prev %= mod
            b_cur %= mod
            z_prev %= mod
            z_cur %= mod
        b_prev, b_cur = b_cur % mod, nb
        z_prev, z_cur = z_cur % mod, nz
        b.append(b_cur % q**3)
        z.append(z_cur % q**3)
    return b, z


def target_map(q_limit: int, b_small):
    out = {}
    for q in primes_upto(q_limit):
        if q < 5:
            continue
        rows = [r for r in range(1, q - 1) if b_small[r] % q == 0]
        if rows:
            out[q] = rows
    return out


def independent_audit(q_limit: int) -> dict:
    b_small = exact_apery(q_limit - 1)
    a_small = exact_companion((q_limit - 1) // 2)
    jet_b, W, H = exact_jets(q_limit - 1)
    targets = target_map(q_limit, b_small)

    checked_endpoints = 0
    checked_targets = 0
    lower_zero_endpoint_rows = 0
    first_b_endpoint_failure = None
    first_scalar_q2_failure = None
    first_full_endpoint_p3_failure = None
    first_b_transport_failure = None
    first_companion_transport_failure = None
    first_theta_failure = None
    first_xi_failure = None
    first_wedge_failure = None
    first_target_projective_failure = None
    lower_zero_details = []

    for q in primes_upto(q_limit):
        if q < 5:
            continue
        high_max = (q - 1) // 2
        limit = high_max * q + (q - 2)
        bmod, zmod = qadic_recurrence_pair(q, limit, 3)
        q2, q3 = q**2, q**3

        for high in range(1, high_max + 1):
            checked_endpoints += 1
            aa = a_small[high]
            aa_mod = mod_fraction(aa, q3)
            ba = b_small[high]
            idx0 = high * q
            if (bmod[idx0] - ba) % q3 and first_b_endpoint_failure is None:
                first_b_endpoint_failure = (q, high, (bmod[idx0] - ba) % q3)
            scalar = (zmod[idx0] - aa_mod) % q3
            if scalar % q2 and first_scalar_q2_failure is None:
                first_scalar_q2_failure = (q, high, scalar)
            delta0 = (zmod[idx0] * (ba % q3) - aa_mod * bmod[idx0]) % q3
            if delta0 and first_full_endpoint_p3_failure is None:
                first_full_endpoint_p3_failure = (q, high, delta0)
            if ba % q == 0:
                lower_zero_endpoint_rows += 1
                lower_zero_details.append((q, high, scalar // q2 % q, delta0))

        for r in targets.get(q, []):
            br = b_small[r]
            Wr = mod_fraction(W[r], q3)
            Hr = mod_fraction(H[r], q3)
            assert mod_fraction(jet_b[r], q3) == br % q3
            for high in range(1, high_max + 1):
                checked_targets += 1
                n = high * q + r
                ba = b_small[high]
                aa = a_small[high]
                A, E = aa.numerator, aa.denominator
                aa_mod = mod_fraction(aa, q3)
                first_b = ba * ((br % q3) + high * q * Wr)
                delta_b = (bmod[n] - first_b) % q3
                if delta_b % q2:
                    if first_b_transport_failure is None:
                        first_b_transport_failure = (q, high, r, delta_b)
                    continue
                theta = delta_b // q2 % q

                endpoint_scalar = (zmod[high * q] - aa_mod) % q3
                if endpoint_scalar % q2:
                    continue
                endpoint_theta = endpoint_scalar // q2 % q
                expected_z = (
                    aa_mod * ((br % q3) + high * q * Wr + high**2 * q2 * Hr)
                    + q2 * endpoint_theta * (br % q3)
                ) % q3
                if (zmod[n] - expected_z) % q3 and first_companion_transport_failure is None:
                    first_companion_transport_failure = (
                        q,
                        high,
                        r,
                        (zmod[n] - expected_z) % q3,
                    )

                first_a = mod_fraction(aa * (br + high * q * W[r]), q3)
                delta_a = (zmod[n] - first_a) % q3
                if delta_a % q2:
                    if first_companion_transport_failure is None:
                        first_companion_transport_failure = (q, high, r, delta_a)
                    continue
                xi = E * (delta_a // q2 % q) % q
                expected_theta = high**2 * (ba % q) * (Hr % q) % q
                expected_xi = high**2 * (A % q) * (Hr % q) % q
                if theta != expected_theta and first_theta_failure is None:
                    first_theta_failure = (q, high, r, theta, expected_theta)
                if xi != expected_xi and first_xi_failure is None:
                    first_xi_failure = (q, high, r, xi, expected_xi)
                wedge = (ba * xi - A * theta) % q
                if wedge and first_wedge_failure is None:
                    first_wedge_failure = (q, high, r, theta, xi, wedge)
                projective = (
                    zmod[n] * (ba % q3) - aa_mod * bmod[n]
                ) % q3
                if projective and first_target_projective_failure is None:
                    first_target_projective_failure = (q, high, r, projective)

    return {
        "q_limit": q_limit,
        "target_primes": targets,
        "checked_endpoints": checked_endpoints,
        "checked_targets": checked_targets,
        "lower_zero_endpoint_rows": lower_zero_endpoint_rows,
        "lower_zero_endpoint_details": lower_zero_details,
        "first_b_endpoint_failure": first_b_endpoint_failure,
        "first_scalar_q2_failure": first_scalar_q2_failure,
        "first_full_endpoint_p3_failure": first_full_endpoint_p3_failure,
        "first_b_transport_failure": first_b_transport_failure,
        "first_companion_transport_failure": first_companion_transport_failure,
        "first_theta_failure": first_theta_failure,
        "first_xi_failure": first_xi_failure,
        "first_wedge_failure": first_wedge_failure,
        "first_target_projective_failure": first_target_projective_failure,
    }


def finite_audit(q_limit: int) -> dict:
    b_small = exact_apery(q_limit - 1)
    a_small = exact_companion((q_limit - 1) // 2)
    _, _, H = exact_jets(q_limit - 1)
    targets = target_map(q_limit, b_small)

    records = 0
    lower_zero_outer_rows = 0
    first_failure = None
    first_wedge_failure = None
    first_theta_formula_failure = None
    first_xi_formula_failure = None

    for q, rs in targets.items():
        for r in rs:
            Hr = mod_fraction(H[r], q)
            for high in range(1, (q - 1) // 2 + 1):
                records += 1
                if b_small[high] % q == 0:
                    lower_zero_outer_rows += 1
                try:
                    beta, theta, xi, wedge = finite_digits(
                        q, high, r, b_small, a_small
                    )
                except Exception as exc:  # report exact first failing tuple
                    first_failure = (q, high, r, repr(exc))
                    return {
                        "q_limit": q_limit,
                        "target_primes": targets,
                        "records": records,
                        "lower_zero_outer_rows": lower_zero_outer_rows,
                        "first_failure": first_failure,
                        "first_wedge_failure": first_wedge_failure,
                        "first_theta_formula_failure": first_theta_formula_failure,
                        "first_xi_formula_failure": first_xi_formula_failure,
                    }
                aa = a_small[high]
                A = aa.numerator
                expected_theta = high**2 * (b_small[high] % q) * Hr % q
                expected_xi = high**2 * (A % q) * Hr % q
                if wedge and first_wedge_failure is None:
                    first_wedge_failure = (q, high, r, beta, theta, xi, wedge)
                if theta != expected_theta and first_theta_formula_failure is None:
                    first_theta_formula_failure = (
                        q,
                        high,
                        r,
                        theta,
                        expected_theta,
                    )
                if xi != expected_xi and first_xi_formula_failure is None:
                    first_xi_formula_failure = (q, high, r, xi, expected_xi)

    return {
        "q_limit": q_limit,
        "target_primes": targets,
        "records": records,
        "lower_zero_outer_rows": lower_zero_outer_rows,
        "first_failure": first_failure,
        "first_wedge_failure": first_wedge_failure,
        "first_theta_formula_failure": first_theta_formula_failure,
        "first_xi_formula_failure": first_xi_formula_failure,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["finite", "projective"], required=True)
    parser.add_argument("--q-limit", type=int, default=300)
    args = parser.parse_args()
    result = (
        finite_audit(args.q_limit)
        if args.mode == "finite"
        else independent_audit(args.q_limit)
    )
    for key, value in result.items():
        print(f"{key}={value}")


if __name__ == "__main__":
    main()
