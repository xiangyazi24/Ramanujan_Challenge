#!/usr/bin/env python3
"""Exact F_181 computation for Q2306.

Computes the Racah/Newton polynomials

    A_n(T) = _4F_3(-n,n+1,-x,x+1;1,1,1;1),  T=x(x+1),

as A_n(T)=sum_k D(n,k) q_k(T), where
q_k(T)=prod_{j<k}(T-j(j+1))/(k!)^2 and
D(n,k)=binom(n,k)binom(n+k,k).

The script independently verifies the resultant by a Sylvester determinant
and by the Euclidean resultant recursion, then decomposes it through the
monic Racah three-term recurrence and the associated Dirichlet block.
Only Python's standard library is used.
"""
from __future__ import annotations

from hashlib import sha256
import json
from math import comb
from typing import Iterable

P = 181
R = 19
S = 47
Poly = list[int]  # coefficients low to high


def inv(a: int) -> int:
    a %= P
    if a == 0:
        raise ZeroDivisionError("inverse of zero")
    return pow(a, P - 2, P)


def trim(a: Iterable[int]) -> Poly:
    out = [x % P for x in a]
    if not out:
        return [0]
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def zero(a: Poly) -> bool:
    return trim(a) == [0]


def degree(a: Poly) -> int:
    a = trim(a)
    return -1 if zero(a) else len(a) - 1


def add(a: Poly, b: Poly) -> Poly:
    n = max(len(a), len(b))
    return trim([
        (a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)
        for i in range(n)
    ])


def sub(a: Poly, b: Poly) -> Poly:
    n = max(len(a), len(b))
    return trim([
        (a[i] if i < len(a) else 0) - (b[i] if i < len(b) else 0)
        for i in range(n)
    ])


def scale(a: Poly, c: int) -> Poly:
    return trim([(c % P) * x for x in a])


def mul(a: Poly, b: Poly) -> Poly:
    if zero(a) or zero(b):
        return [0]
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] = (out[i + j] + x * y) % P
    return trim(out)


def divmod_poly(a: Poly, b: Poly) -> tuple[Poly, Poly]:
    a = trim(a)
    b = trim(b)
    if zero(b):
        raise ZeroDivisionError("polynomial division by zero")
    if degree(a) < degree(b):
        return [0], a
    q = [0] * (degree(a) - degree(b) + 1)
    r = a[:]
    ib = inv(b[-1])
    while not zero(r) and degree(r) >= degree(b):
        d = degree(r) - degree(b)
        c = r[-1] * ib % P
        q[d] = c
        for j, bj in enumerate(b):
            r[j + d] = (r[j + d] - c * bj) % P
        r = trim(r)
    return trim(q), trim(r)


def mod_poly(a: Poly, b: Poly) -> Poly:
    return divmod_poly(a, b)[1]


def monic(a: Poly) -> Poly:
    a = trim(a)
    if zero(a):
        return [0]
    return scale(a, inv(a[-1]))


def gcd_poly(a: Poly, b: Poly) -> Poly:
    a, b = trim(a), trim(b)
    while not zero(b):
        a, b = b, mod_poly(a, b)
    return monic(a)


def deriv(a: Poly) -> Poly:
    if len(a) <= 1:
        return [0]
    return trim([i * a[i] for i in range(1, len(a))])


def eval_poly(a: Poly, x: int) -> int:
    ans = 0
    for c in reversed(a):
        ans = (ans * x + c) % P
    return ans


def powmod_poly(base: Poly, exponent: int, modulus: Poly) -> Poly:
    out = [1]
    cur = mod_poly(base, modulus)
    e = exponent
    while e:
        if e & 1:
            out = mod_poly(mul(out, cur), modulus)
        cur = mod_poly(mul(cur, cur), modulus)
        e >>= 1
    return out


def determinant_mod(matrix: list[list[int]]) -> int:
    a = [[x % P for x in row] for row in matrix]
    n = len(a)
    det = 1
    for col in range(n):
        pivot = next((row for row in range(col, n) if a[row][col]), None)
        if pivot is None:
            return 0
        if pivot != col:
            a[col], a[pivot] = a[pivot], a[col]
            det = -det
        pv = a[col][col] % P
        det = det * pv % P
        ip = inv(pv)
        for row in range(col + 1, n):
            if not a[row][col]:
                continue
            factor = a[row][col] * ip % P
            for j in range(col, n):
                a[row][j] = (a[row][j] - factor * a[col][j]) % P
    return det % P


def sylvester_resultant(f: Poly, g: Poly) -> int:
    f, g = trim(f), trim(g)
    m, n = degree(f), degree(g)
    if m < 0 or n < 0:
        raise ValueError("zero polynomial has no Sylvester resultant here")
    if m == 0:
        return pow(f[0], n, P)
    if n == 0:
        return pow(g[0], m, P)
    fh = list(reversed(f))
    gh = list(reversed(g))
    size = m + n
    mat: list[list[int]] = []
    for shift in range(n):
        row = [0] * size
        row[shift:shift + m + 1] = fh
        mat.append(row)
    for shift in range(m):
        row = [0] * size
        row[shift:shift + n + 1] = gh
        mat.append(row)
    return determinant_mod(mat)


def euclid_resultant(f: Poly, g: Poly, trace: list[dict] | None = None) -> int:
    f, g = trim(f), trim(g)
    m, n = degree(f), degree(g)
    if m < 0 or n < 0:
        raise ValueError("zero polynomial")
    if n == 0:
        return pow(g[0], m, P)
    if m < n:
        sign = P - 1 if (m * n) & 1 else 1
        if trace is not None:
            trace.append({"swap": True, "deg_f": m, "deg_g": n, "sign": sign})
        return sign * euclid_resultant(g, f, trace) % P
    _, rem = divmod_poly(f, g)
    if zero(rem):
        if trace is not None:
            trace.append({"deg_f": m, "deg_g": n, "remainder": "zero"})
        return 0
    k = degree(rem)
    sign = P - 1 if (m * n) & 1 else 1
    factor = sign * pow(g[-1], m - k, P) % P
    if trace is not None:
        trace.append({
            "deg_f": m,
            "deg_g": n,
            "deg_remainder": k,
            "lc_g": g[-1],
            "lc_exponent": m - k,
            "sign": sign,
            "factor": factor,
        })
    return factor * euclid_resultant(g, rem, trace) % P


def lambda_node(n: int) -> int:
    return n * (n + 1) % P


def D(n: int, k: int) -> int:
    if k < 0 or k > n:
        return 0
    return comb(n, k) * comb(n + k, k) % P


def build_q(max_n: int) -> list[Poly]:
    qs: list[Poly] = [[1]]
    for k in range(max_n):
        nxt = mul(qs[-1], [(-lambda_node(k)) % P, 1])
        nxt = scale(nxt, inv((k + 1) ** 2))
        qs.append(nxt)
    return qs


def A(n: int, qs: list[Poly]) -> Poly:
    out = [0]
    for k in range(n + 1):
        out = add(out, scale(qs[k], D(n, k)))
    return trim(out)


def K(n: int, m: int) -> int:
    return sum(D(n, k) * D(m, k) for k in range(min(n, m) + 1)) % P


def lc_formula(n: int) -> int:
    # (2n)!/(n!)^4 modulo p, valid here because 2n<p.
    num = 1
    for j in range(1, 2 * n + 1):
        num = num * j % P
    den = 1
    for j in range(1, n + 1):
        den = den * j % P
    return num * pow(inv(den), 4, P) % P


def alpha(n: int) -> int:
    return (n * n + n + 1) * inv(2) % P


def beta(n: int) -> int:
    if n == 0:
        return 0
    return pow(n, 6, P) * inv(4 * (4 * n * n - 1)) % P


def recurrence_polynomials(max_n: int, As: list[Poly]) -> list[Poly]:
    mons = [monic(a) for a in As]
    assert mons[0] == [1]
    for n in range(max_n):
        rhs = mul([alpha(n), 1], mons[n])
        if n:
            rhs = sub(rhs, scale(mons[n - 1], beta(n)))
        assert trim(rhs) == mons[n + 1], (n, rhs, mons[n + 1])
    return mons


def associated_block(r: int, h: int) -> Poly:
    # S_0=1 and P_{r+j}=P_{r+1} S_{j-1} on P_r=0.
    if h == 1:
        return [1]
    blocks: list[Poly] = [[1], [alpha(r + 1), 1]]
    while len(blocks) <= h - 1:
        j = len(blocks) - 1
        idx = r + j + 1
        nxt = sub(mul([alpha(idx), 1], blocks[j]), scale(blocks[j - 1], beta(idx)))
        blocks.append(trim(nxt))
    return blocks[h - 1]


def consecutive_resultant_formula(n: int) -> int:
    ans = 1
    for j in range(1, n + 1):
        ans = ans * pow((-beta(j)) % P, j, P) % P
    return ans


def roots_in_base_field(f: Poly) -> list[int]:
    return [x for x in range(P) if eval_poly(f, x) == 0]


def canonical_digest(obj: object) -> str:
    raw = json.dumps(obj, sort_keys=True, separators=(",", ":")).encode()
    return sha256(raw).hexdigest()


def discrete_log_table() -> tuple[int, dict[int, int]]:
    factors = [2, 3, 5]
    primitive = next(
        g for g in range(2, P)
        if all(pow(g, (P - 1) // q, P) != 1 for q in factors)
    )
    table: dict[int, int] = {}
    x = 1
    for e in range(P - 1):
        table[x] = e
        x = x * primitive % P
    assert len(table) == P - 1
    return primitive, table


def main() -> None:
    qs = build_q(S)
    As = [A(n, qs) for n in range(S + 1)]
    mons = recurrence_polynomials(S, As)
    ar, ass = As[R], As[S]
    pr, ps = mons[R], mons[S]

    # Normalization and evaluation checks.
    assert degree(ar) == R and degree(ass) == S
    assert ar[-1] == lc_formula(R)
    assert ass[-1] == lc_formula(S)
    for n in (R, S):
        for m in range(S + 1):
            assert eval_poly(As[n], lambda_node(m)) == K(n, m)
    assert K(R, S) == K(S, R)

    br = eval_poly(ar, lambda_node(R))
    bs = eval_poly(ass, lambda_node(S))
    cross1 = eval_poly(ar, lambda_node(S))
    cross2 = eval_poly(ass, lambda_node(R))
    assert br == 0 and bs == 0
    assert cross1 == cross2 == K(R, S)

    # Resultant, independently by two algorithms.
    trace: list[dict] = []
    res_e = euclid_resultant(ar, ass, trace)
    res_s = sylvester_resultant(ar, ass)
    assert res_e == res_s
    g = gcd_poly(ar, ass)
    assert (res_e == 0) == (degree(g) > 0)

    # Monic recurrence/associated-block factorization.
    h = S - R
    block = associated_block(R, h)
    assert degree(block) == h - 1
    congruence = mod_poly(sub(ps, mul(mons[R + 1], block)), pr)
    assert zero(congruence)

    res_monic = euclid_resultant(pr, ps)
    res_adj_direct = euclid_resultant(pr, mons[R + 1])
    res_adj_formula = consecutive_resultant_formula(R)
    assert res_adj_direct == res_adj_formula
    res_block = euclid_resultant(pr, block)
    assert res_monic == res_adj_direct * res_block % P

    scaling = pow(ar[-1], S, P) * pow(ass[-1], R, P) % P
    assert res_e == scaling * res_monic % P

    # More exact structural diagnostics.
    roots_r = roots_in_base_field(ar)
    roots_s = roots_in_base_field(ass)
    common_base = sorted(set(roots_r) & set(roots_s))
    assert not common_base if res_e else True
    sqfree_r = gcd_poly(ar, deriv(ar))
    sqfree_s = gcd_poly(ass, deriv(ass))

    primitive, logs = discrete_log_table()
    nonzero_components = {
        "lc_r": ar[-1],
        "lc_s": ass[-1],
        "scaling": scaling,
        "adjacent_resultant": res_adj_direct,
        "associated_block_resultant": res_block,
        "monic_resultant": res_monic,
        "full_resultant": res_e,
        "cross_K": cross1,
    }
    discrete_logs = {k: logs[v] if v else None for k, v in nonzero_components.items()}

    data = {
        "p": P,
        "r": R,
        "s": S,
        "h": h,
        "lambda_r": lambda_node(R),
        "lambda_s": lambda_node(S),
        "b_r": br,
        "b_s": bs,
        "K_rs": cross1,
        "A_r_degree": degree(ar),
        "A_s_degree": degree(ass),
        "A_r_coefficients_low_to_high": ar,
        "A_s_coefficients_low_to_high": ass,
        "A_r_leading_coefficient": ar[-1],
        "A_s_leading_coefficient": ass[-1],
        "gcd_coefficients_low_to_high": g,
        "gcd_degree": degree(g),
        "resultant_sylvester": res_s,
        "resultant_euclid": res_e,
        "monic_resultant": res_monic,
        "scaling_factor_lc_r_pow_s_lc_s_pow_r": scaling,
        "consecutive_resultant": res_adj_direct,
        "consecutive_resultant_formula": res_adj_formula,
        "associated_block_degree": degree(block),
        "associated_block_coefficients_low_to_high": block,
        "associated_block_resultant": res_block,
        "base_field_roots_A_r": roots_r,
        "base_field_roots_A_s": roots_s,
        "common_base_field_roots": common_base,
        "gcd_A_r_derivative": sqfree_r,
        "gcd_A_s_derivative": sqfree_s,
        "primitive_root": primitive,
        "component_discrete_logs": discrete_logs,
        "euclidean_resultant_trace": trace,
    }
    data["sha256"] = canonical_digest(data)

    print("ANSWER Q2306 9a7673b0")
    print()
    print("# Exact F_181 computation")
    print()
    print("```json")
    print(json.dumps(data, sort_keys=True, separators=(",", ":")))
    print("```")
    print()
    print("# Human-readable summary")
    print()
    for key in (
        "lambda_r", "lambda_s", "b_r", "b_s", "K_rs",
        "resultant_sylvester", "gcd_degree", "A_r_leading_coefficient",
        "A_s_leading_coefficient", "scaling_factor_lc_r_pow_s_lc_s_pow_r",
        "consecutive_resultant", "associated_block_resultant", "monic_resultant",
        "primitive_root",
    ):
        print(f"- `{key}` = `{data[key]}`")
    print(f"- roots of A_{R} in F_{P}: `{roots_r}`")
    print(f"- roots of A_{S} in F_{P}: `{roots_s}`")
    print(f"- gcd(A_{R},A_{S}) = `{g}`")
    print(f"- SHA-256 = `{data['sha256']}`")


if __name__ == "__main__":
    main()
