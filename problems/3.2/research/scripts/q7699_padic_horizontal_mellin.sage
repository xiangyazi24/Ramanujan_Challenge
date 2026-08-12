#!/usr/bin/env sage -python
# Exact Sage experiment for Q7699.  This file is ordinary Python syntax and is
# meant to run as `sage -python ...`; Sage supplies exact arithmetic while no
# preparser writes beside the source file.
#
# For fixed n and n/2 < p <= n, the interior index is r=n-p.  Q7690 identifies
# the defining-characteristic Mellin zero with b_r == 0 mod p.  We therefore
# record the rigorously equivalent Gross--Koblitz/Morita carry profile of
#   b_r = sum_k (binom(r,k) binom(r+k,k))^2.
# The index r is a Mellin exponent/coefficient index, never a geometric fiber.

from sage.all import ZZ, binomial, inverse_mod, prime_range
import argparse
import csv
import io


def vp_int(x, p):
    x = ZZ(x)
    assert x != 0
    v = 0
    while x % p == 0:
        x //= p
        v += 1
    return v


def morita_gamma_integer(a, p):
    assert a >= 1
    prod = ZZ(1)
    for j in range(1, a):
        if j % p:
            prod *= j
    return prod if a % 2 == 0 else -prod


def gamma_unit_mod(r, k, p):
    num = morita_gamma_integer(r + k + 1, p) % p
    den = (
        morita_gamma_integer(k + 1, p) ** 2
        * morita_gamma_integer(r - k + 1, p)
    ) % p
    assert den % p != 0
    return ZZ(num * inverse_mod(den, p) % p)


def apery_term(r, k):
    return ZZ((binomial(r, k) * binomial(r + k, k)) ** 2)


def profile(p, r):
    vals = []
    min_units = []
    gamma_units = []
    bmod = ZZ(0)
    for k in range(r + 1):
        T = apery_term(r, k)
        v = vp_int(T, p)
        vals.append(v)
        bmod = (bmod + T) % p
        if v == 0:
            u = ZZ(T % p)
            min_units.append(u)
            gu = gamma_unit_mod(r, k, p)
            gamma_units.append(ZZ(gu ** 2 % p))
            assert u == gu ** 2 % p
        else:
            assert v == 2

    min_count = min(r, p - 1 - r) + 1
    high_count = max(0, 2 * r + 1 - p)
    assert vals == [0] * min_count + [2] * high_count
    assert len(vals) == r + 1
    assert min_units == gamma_units
    theta = ZZ(sum(min_units) % p)
    assert theta == bmod
    return {
        "min_count": min_count,
        "high_count": high_count,
        "theta": theta,
        "bmod": bmod,
        "zero": int(theta == 0),
        "unit_first": min_units[0],
        "unit_last": min_units[-1],
        "unit_checksum": ZZ(sum((i + 1) * u for i, u in enumerate(min_units)) % p),
    }


def row_for(n, p):
    m = n - 1
    residue = m % (p - 1)
    coeff = n - p
    if coeff == 0:
        return {
            "n": n, "p": p, "m": m, "r": 0, "status": "upper_endpoint",
            "min_count": 1, "high_count": 0, "theta": 1, "bmod": 1,
            "zero": 0, "unit_first": 1, "unit_last": 1, "unit_checksum": 1,
        }
    if coeff == p - 1:
        bmod = sum(apery_term(coeff, k) for k in range(coeff + 1)) % p
        assert bmod == 1
        return {
            "n": n, "p": p, "m": m, "r": 0, "status": "lower_endpoint",
            "min_count": 1, "high_count": p - 1, "theta": 1, "bmod": 1,
            "zero": 0, "unit_first": 1, "unit_last": 1, "unit_checksum": 1,
        }
    assert 1 <= coeff <= p - 2
    assert residue == coeff
    ans = profile(p, coeff)
    ans.update({"n": n, "p": p, "m": m, "r": coeff, "status": "interior"})
    return ans


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=321)
    args = ap.parse_args()
    n = ZZ(args.n)
    rows = [row_for(n, ZZ(p)) for p in prime_range(n // 2 + 1, n + 1)]
    fields = [
        "n", "p", "m", "r", "status", "min_count", "high_count",
        "theta", "bmod", "zero", "unit_first", "unit_last", "unit_checksum",
    ]

    interior = [row for row in rows if row["status"] == "interior"]
    bad = [row for row in rows if row["zero"]]
    print("Q7699_FIXED_N", n)
    print("top_half_primes", len(rows))
    print("interior", len(interior))
    print("bad_count", len(bad))
    print("bad_pairs", [(int(row["p"]), int(row["r"])) for row in bad])
    print("complexity_min", min(int(row["r"]) + 1 for row in interior))
    print("complexity_max", max(int(row["r"]) + 1 for row in interior))
    print("min_slope_terms_min", min(int(row["min_count"]) for row in interior))
    print("min_slope_terms_max", max(int(row["min_count"]) for row in interior))
    print("CSV_BEGIN")
    buf = io.StringIO()
    writer = csv.DictWriter(buf, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    print(buf.getvalue(), end="")
    print("CSV_END")
    print("Q7699_SAGE_VERIFY PASS")


if __name__ == "__main__":
    main()
