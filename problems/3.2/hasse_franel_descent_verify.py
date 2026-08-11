#!/usr/bin/env python3
"""Exact finite-field checks for hasse_franel_descent.tex.

The verifier checks the polynomial (not merely pointwise) Franel pullback,
the two universal-square-root recurrences modulo p, and the CFVZ branch
factorization for every prime 5 <= p <= 199.
"""

from hashlib import sha256
import json


LIMIT = 199
EXPECTED_DIGEST = "90a8aea73b9bda5589e7335aab094bc0b72a82aae8f7a495f35eee78018afa25"


def primes_through(limit):
    out = []
    for n in range(2, limit + 1):
        if all(n % q for q in range(2, int(n**0.5) + 1)):
            out.append(n)
    return out


def trim(poly):
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    return poly


def add(left, right, p):
    out = [0] * max(len(left), len(right))
    for i in range(len(out)):
        out[i] = ((left[i] if i < len(left) else 0)
                  + (right[i] if i < len(right) else 0)) % p
    return trim(out)


def mul(left, right, p):
    out = [0] * (len(left) + len(right) - 1)
    for i, x in enumerate(left):
        for j, y in enumerate(right):
            out[i + j] = (out[i + j] + x * y) % p
    return trim(out)


def apery_truncation(p):
    values = [1]
    for n in range(p - 1):
        previous = values[n - 1] if n else 0
        middle = (2 * n + 1) * (17 * n * n + 17 * n + 5)
        numerator = middle * values[n] - n**3 * previous
        values.append(numerator * pow((n + 1) ** 3, -1, p) % p)
    return values


def franel_truncation(p):
    rows = [[1]]
    for n in range(1, p):
        previous = rows[-1]
        rows.append([1] + [(previous[k - 1] + previous[k]) % p
                           for k in range(1, n)] + [1])
    return [sum(entry**3 for entry in row) % p for row in rows]


def homogenized_pullback(coefficients, p):
    """Return sum b_m q^m l^(p-1-m), q=x(1-8x), l=1+x."""
    q = [0, 1, -8 % p]
    linear_power = [1]
    result = [coefficients[-1]]
    for m in range(p - 2, -1, -1):
        linear_power = mul(linear_power, [1, 1], p)
        result = add(mul(q, result, p),
                     [(coefficients[m] * x) % p for x in linear_power], p)
    return result


def square_root_branch(p, defect):
    degree = (p - 3) // 2 if defect else (p - 1) // 2
    values = [1]
    previous = 0
    for r in range(degree):
        if defect:
            middle = 2 * (68 * r * r + 102 * r + 39)
            trailing = (2 * r + 1) ** 2
        else:
            middle = 2 * (68 * r * r + 34 * r + 5)
            trailing = (2 * r - 1) ** 2
        numerator = middle * values[r] - trailing * previous
        following = numerator * pow(4 * (r + 1) ** 2, -1, p) % p
        previous = values[r]
        values.append(following)
    return values


def legendre(value, p):
    residue = pow(value % p, (p - 1) // 2, p)
    return -1 if residue == p - 1 else residue


def main():
    records = []
    for p in primes_through(LIMIT):
        if p < 5:
            continue
        apery = apery_truncation(p)
        franel = franel_truncation(p)

        left = mul(franel, franel, p)
        right = homogenized_pullback(apery, p)
        assert left == right, (p, "Franel pullback")

        defect = legendre(-6, p) == -1
        root = square_root_branch(p, defect)
        factorized = mul(root, root, p)
        if defect:
            factorized = mul([1, -34 % p, 1], factorized, p)
        assert factorized == trim(apery[:]), (p, "square-root branch")

        zeros = sum(value == 0 for value in apery)
        records.append({
            "p": p,
            "branch": "defect" if defect else "square",
            "degree": len(root) - 1,
            "zeros": zeros,
            "pullback_degree": len(left) - 1,
        })

    payload = json.dumps(records, sort_keys=True, separators=(",", ":"))
    digest = sha256(payload.encode()).hexdigest()
    if EXPECTED_DIGEST is not None:
        assert digest == EXPECTED_DIGEST, (digest, EXPECTED_DIGEST)
    print(f"PASS: {len(records)} primes, 5 <= p <= {LIMIT}")
    print(f"SHA256: {digest}")


if __name__ == "__main__":
    main()
