#!/usr/bin/env python3
"""Exact p=41 Newton/Racah saturated Bezout certificate for Q2273.

Dependency-free: Python standard library only.  Coefficients are always reduced
into {0,...,40}.  A bivariate sparse coefficient triple [i,j,c] means
c*T^i*U^j.
"""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

P = 41
MHALF = 20
UPoly = List[int]                    # low-to-high coefficients
BPoly = Dict[Tuple[int, int], int]   # (T-degree,U-degree) -> coefficient


def inv(a: int) -> int:
    a %= P
    if a == 0:
        raise ZeroDivisionError("inverse of zero")
    return pow(a, P - 2, P)


def utrim(a: Iterable[int]) -> UPoly:
    out = [x % P for x in a]
    if not out:
        return [0]
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def uzero(a: UPoly) -> bool:
    return utrim(a) == [0]


def uadd(a: UPoly, b: UPoly) -> UPoly:
    n = max(len(a), len(b))
    return utrim([
        (a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)
        for i in range(n)
    ])


def usub(a: UPoly, b: UPoly) -> UPoly:
    n = max(len(a), len(b))
    return utrim([
        (a[i] if i < len(a) else 0) - (b[i] if i < len(b) else 0)
        for i in range(n)
    ])


def uscale(a: UPoly, c: int) -> UPoly:
    return utrim([(c % P) * x for x in a])


def umul(a: UPoly, b: UPoly) -> UPoly:
    if uzero(a) or uzero(b):
        return [0]
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] = (out[i + j] + x * y) % P
    return utrim(out)


def udivmod(a: UPoly, b: UPoly) -> Tuple[UPoly, UPoly]:
    a = utrim(a)
    b = utrim(b)
    if uzero(b):
        raise ZeroDivisionError("polynomial division by zero")
    if len(a) < len(b):
        return [0], a
    q = [0] * (len(a) - len(b) + 1)
    r = a[:]
    ib = inv(b[-1])
    while not uzero(r) and len(r) >= len(b):
        shift = len(r) - len(b)
        c = r[-1] * ib % P
        q[shift] = (q[shift] + c) % P
        for i, bi in enumerate(b):
            r[shift + i] = (r[shift + i] - c * bi) % P
        r = utrim(r)
    return utrim(q), utrim(r)


def umonic(a: UPoly) -> UPoly:
    a = utrim(a)
    if uzero(a):
        return [0]
    return uscale(a, inv(a[-1]))


def ugcd(a: UPoly, b: UPoly) -> UPoly:
    a, b = utrim(a), utrim(b)
    while not uzero(b):
        _, r = udivmod(a, b)
        a, b = b, r
    return umonic(a)


def uxgcd(a: UPoly, b: UPoly) -> Tuple[UPoly, UPoly, UPoly]:
    """Return monic d,s,t with d=s*a+t*b."""
    r0, r1 = utrim(a), utrim(b)
    s0, s1 = [1], [0]
    t0, t1 = [0], [1]
    while not uzero(r1):
        q, r2 = udivmod(r0, r1)
        r0, r1 = r1, r2
        s0, s1 = s1, usub(s0, umul(q, s1))
        t0, t1 = t1, usub(t0, umul(q, t1))
    z = inv(r0[-1])
    return uscale(r0, z), uscale(s0, z), uscale(t0, z)


def uderiv(a: UPoly) -> UPoly:
    if len(a) <= 1:
        return [0]
    return utrim([i * a[i] for i in range(1, len(a))])


def ueval(a: UPoly, x: int) -> int:
    x %= P
    y = 0
    for c in reversed(a):
        y = (y * x + c) % P
    return y


def uexpr(a: UPoly, var: str) -> str:
    a = utrim(a)
    terms: List[str] = []
    for i in range(len(a) - 1, -1, -1):
        c = a[i] % P
        if c == 0:
            continue
        if i == 0:
            mon = str(c)
        elif i == 1:
            mon = var if c == 1 else f"{c}*{var}"
        else:
            mon = f"{var}^{i}" if c == 1 else f"{c}*{var}^{i}"
        terms.append(mon)
    return " + ".join(terms) if terms else "0"


def bnorm(f: BPoly) -> BPoly:
    return {ij: c % P for ij, c in f.items() if c % P}


def badd(f: BPoly, g: BPoly) -> BPoly:
    out = dict(f)
    for ij, c in g.items():
        out[ij] = (out.get(ij, 0) + c) % P
    return bnorm(out)


def bsub(f: BPoly, g: BPoly) -> BPoly:
    out = dict(f)
    for ij, c in g.items():
        out[ij] = (out.get(ij, 0) - c) % P
    return bnorm(out)


def bscale(f: BPoly, c: int) -> BPoly:
    return bnorm({ij: c * x for ij, x in f.items()})


def bmul(f: BPoly, g: BPoly) -> BPoly:
    out: BPoly = {}
    for (i, j), x in f.items():
        for (k, ell), y in g.items():
            ij = (i + k, j + ell)
            out[ij] = (out.get(ij, 0) + x * y) % P
    return bnorm(out)


def bembed_t(a: UPoly) -> BPoly:
    return bnorm({(i, 0): c for i, c in enumerate(a)})


def bembed_u(a: UPoly) -> BPoly:
    return bnorm({(0, j): c for j, c in enumerate(a)})


def beval(f: BPoly, t: int, u: int) -> int:
    t %= P
    u %= P
    tp = [1]
    up = [1]
    dt = max((i for i, _ in f), default=0)
    du = max((j for _, j in f), default=0)
    for _ in range(dt):
        tp.append(tp[-1] * t % P)
    for _ in range(du):
        up.append(up[-1] * u % P)
    return sum(c * tp[i] * up[j] for (i, j), c in f.items()) % P


def bdegrees(f: BPoly) -> Tuple[int, int]:
    if not f:
        return (-1, -1)
    return max(i for i, _ in f), max(j for _, j in f)


def bdiv_t(f: BPoly, g: UPoly) -> Tuple[BPoly, BPoly]:
    """Divide by monic g(T), with coefficients in F_41[U]."""
    g = utrim(g)
    if uzero(g):
        raise ZeroDivisionError
    dg = len(g) - 1
    if dg == 0:
        return bscale(f, inv(g[0])), {}
    if g[-1] != 1:
        raise ValueError("bdiv_t expects a monic divisor")
    q: BPoly = {}
    r = dict(f)
    while r:
        top = max(i for i, _ in r)
        if top < dg:
            break
        top_terms = sorted((j, c) for (i, j), c in r.items() if i == top)
        shift = top - dg
        for j, c in top_terms:
            q[(shift, j)] = (q.get((shift, j), 0) + c) % P
            for k, gk in enumerate(g):
                ij = (shift + k, j)
                r[ij] = (r.get(ij, 0) - c * gk) % P
                if r[ij] == 0:
                    r.pop(ij, None)
        r = bnorm(r)
    return bnorm(q), bnorm(r)


def bdiv_u(f: BPoly, g: UPoly) -> Tuple[BPoly, BPoly]:
    """Divide by monic g(U), with coefficients in F_41[T]."""
    g = utrim(g)
    if uzero(g):
        raise ZeroDivisionError
    dg = len(g) - 1
    if dg == 0:
        return bscale(f, inv(g[0])), {}
    if g[-1] != 1:
        raise ValueError("bdiv_u expects a monic divisor")
    q: BPoly = {}
    r = dict(f)
    while r:
        top = max(j for _, j in r)
        if top < dg:
            break
        top_terms = sorted((i, c) for (i, j), c in r.items() if j == top)
        shift = top - dg
        for i, c in top_terms:
            q[(i, shift)] = (q.get((i, shift), 0) + c) % P
            for k, gk in enumerate(g):
                ij = (i, shift + k)
                r[ij] = (r.get(ij, 0) - c * gk) % P
                if r[ij] == 0:
                    r.pop(ij, None)
        r = bnorm(r)
    return bnorm(q), bnorm(r)


def bsparse(f: BPoly) -> List[List[int]]:
    return [[i, j, f[(i, j)] % P] for i, j in sorted(f, key=lambda z: (z[1], z[0]))]


def bmatrix(f: BPoly, dt: int, du: int) -> List[List[int]]:
    """Rows are U^j, columns are T^i."""
    return [[f.get((i, j), 0) % P for i in range(dt + 1)] for j in range(du + 1)]


def canon_u(a: UPoly) -> bytes:
    return json.dumps(utrim(a), separators=(",", ":")).encode()


def canon_b(f: BPoly) -> bytes:
    return json.dumps(bsparse(f), separators=(",", ":")).encode()


def digest_u(a: UPoly) -> str:
    return sha256(canon_u(a)).hexdigest()


def digest_b(f: BPoly) -> str:
    return sha256(canon_b(f)).hexdigest()


def outer(a: UPoly, b: UPoly) -> BPoly:
    return bnorm({(i, j): x * y for i, x in enumerate(a) for j, y in enumerate(b)})


def build_objects() -> Tuple[List[UPoly], List[int], UPoly, BPoly]:
    psi: List[UPoly] = [[1]]
    for k in range(MHALF):
        psi.append(umul(psi[-1], [-(k * (k + 1)), 1]))
    tails: List[int] = []
    weights: List[int] = []
    for k in range(MHALF + 1):
        tail = 1
        for j in range(k + 1, MHALF + 1):
            tail = tail * j % P
        tails.append(tail)
        weights.append(pow(tail, 4, P))
    a: UPoly = [0]
    kernel: BPoly = {}
    for k in range(MHALF + 1):
        a = uadd(a, uscale(umul(psi[k], psi[k]), weights[k]))
        kernel = badd(kernel, bscale(outer(psi[k], psi[k]), weights[k]))
    return psi, weights, a, kernel


def lagrange_basis(g: UPoly, roots: List[int]) -> Dict[int, UPoly]:
    gp = uderiv(g)
    out: Dict[int, UPoly] = {}
    for a in roots:
        q, r = udivmod(g, [-a, 1])
        assert uzero(r)
        denom = ueval(gp, a)
        assert denom != 0
        out[a] = uscale(q, inv(denom))
        for b in roots:
            assert ueval(out[a], b) == (1 if a == b else 0)
    return out


def build_certificate(
    label: str,
    target: BPoly,
    kernel: BPoly,
    a_poly: UPoly,
    node: UPoly,
    g: UPoly,
    roots: List[int],
    bez_s: UPoly,
    bez_t: UPoly,
) -> dict:
    bad = []
    for x in roots:
        for y in roots:
            mv = beval(target, x, y)
            kv = beval(kernel, x, y)
            if mv != 0 and kv == 0:
                bad.append((x, y, mv, kv))
    if bad:
        return {"label": label, "success": False, "bad": bad}

    lag = lagrange_basis(g, roots)
    bco: BPoly = {}
    for x in roots:
        for y in roots:
            mv = beval(target, x, y)
            if mv == 0:
                continue
            kv = beval(kernel, x, y)
            assert kv != 0
            bco = badd(bco, bscale(outer(lag[x], lag[y]), mv * inv(kv)))

    for x in roots:
        for y in roots:
            assert beval(bco, x, y) * beval(kernel, x, y) % P == beval(target, x, y)

    residual = bsub(target, bmul(bco, kernel))
    qt, rem_t = bdiv_t(residual, g)
    qu, rem = bdiv_u(rem_t, g)
    assert not rem
    assert residual == badd(bmul(qt, bembed_t(g)), bmul(qu, bembed_u(g)))

    cco = bmul(qt, bembed_t(bez_s))
    dco = bmul(qu, bembed_u(bez_s))
    eco = bmul(qt, bembed_t(bez_t))
    fco = bmul(qu, bembed_u(bez_t))

    rhs = badd(bmul(bco, kernel), bmul(cco, bembed_t(a_poly)))
    rhs = badd(rhs, bmul(dco, bembed_u(a_poly)))
    rhs = badd(rhs, bmul(eco, bembed_t(node)))
    rhs = badd(rhs, bmul(fco, bembed_u(node)))
    assert bnorm(rhs) == bnorm(target)

    return {
        "label": label,
        "success": True,
        "B": bco,
        "QT": qt,
        "QU": qu,
        "C": cco,
        "D": dco,
        "E": eco,
        "F": fco,
        "residual": residual,
    }


def print_json_block(obj: object) -> None:
    print("```json")
    print(json.dumps(obj, separators=(",", ":")))
    print("```")


def print_poly_data(name: str, f: BPoly) -> None:
    print(f"### `{name}`")
    print(f"- bidegree: `{bdegrees(f)}`")
    print(f"- nonzero terms: `{len(f)}`")
    print(f"- SHA-256 of canonical sparse triples: `{digest_b(f)}`")
    print("- sparse triples `[i,j,c] = c*T^i*U^j`:")
    print_json_block(bsparse(f))


def main() -> None:
    psi, weights, a_poly, kernel = build_objects()
    node = [0] * 42
    node[1] = -1
    node[41] = 1
    node = utrim(node)

    g = ugcd(a_poly, node)
    roots = [x for x in range(P) if ueval(a_poly, x) == 0]
    assert len(roots) == len(g) - 1
    root_product: UPoly = [1]
    for x in roots:
        root_product = umul(root_product, [-x, 1])
    assert umonic(root_product) == g

    d, bez_s, bez_t = uxgcd(a_poly, node)
    assert d == g
    assert uadd(umul(bez_s, a_poly), umul(bez_t, node)) == g

    t_minus_u: BPoly = {(1, 0): 1, (0, 1): -1}
    pole_t: BPoly = {(0, 0): 1, (1, 0): 4}
    pole_u: BPoly = {(0, 0): 1, (0, 1): 4}
    pole = bmul(pole_t, pole_u)
    saturated_m = bmul(t_minus_u, pole)

    cert_diag = build_certificate(
        "strong diagonal/mate divisor M0=T-U",
        t_minus_u, kernel, a_poly, node, g, roots, bez_s, bez_t,
    )
    cert_sat = build_certificate(
        "saturated divisor M=(T-U)(1+4T)(1+4U)",
        saturated_m, kernel, a_poly, node, g, roots, bez_s, bez_t,
    )

    half_nodes = {n * (n + 1) % P: n for n in range(MHALF + 1)}
    central = (-inv(4)) % P

    print("ANSWER Q2273 f6666eaf")
    print()
    print("# Exact p=41 result")
    print()
    print("The computation uses the monic Newton polynomials")
    print()
    print("```text")
    print("Psi_0(T)=1,")
    print("Psi_k(T)=product_{j=0}^{k-1}(T-j(j+1)).")
    print("```")
    print()
    print("Thus `Psi_k(lambda_n)=(k!)^2*binom(n,k)*binom(n+k,k)`.  Since")
    print(f"`(20!)^4 mod 41 = {pow(__import__('math').factorial(20), 4, P)}`, the requested normalization satisfies")
    print("`A_41(lambda_n)=b_n` and `K_41(lambda_r,lambda_s)=K(r,s)`.")
    print()
    print("Q2268 specifies the excluded product only schematically.  In the spectral")
    print("coordinate `lambda_n=n(n+1)`, diagonal and mate have the same divisor because")
    print("`lambda_r-lambda_s=(r-s)(r+s+1)`.  I therefore compute both the stronger")
    print("target `M0=T-U` and the canonical pole-saturated target")
    print("`M=(T-U)(1+4T)(1+4U)`; `1+4T=0` is the central ramification node.")
    print()

    print("# 1. Weights and A_41(T)")
    print()
    print("Weights `w_k=(20!/k!)^4 mod 41`, k=0,...,20:")
    print_json_block(weights)
    print()
    print(f"`deg A_41 = {len(a_poly)-1}`; coefficient list `[T^0,...,T^{len(a_poly)-1}]`:")
    print_json_block(a_poly)
    print()
    print("Expanded polynomial:")
    print("```text")
    print("A_41(T) = " + uexpr(a_poly, "T"))
    print("```")
    print(f"SHA-256: `{digest_u(a_poly)}`")
    print()

    print("# 2. K_41(T,U)")
    print()
    print(f"Bidegree `{bdegrees(kernel)}`; nonzero terms `{len(kernel)}`; SHA-256 `{digest_b(kernel)}`.")
    print("The following 21x21 matrix has row j equal to coefficients of U^j and")
    print("column i equal to coefficients of T^i:")
    print_json_block(bmatrix(kernel, MHALF, MHALF))
    print()

    print("# 3. Field-valued zero set")
    print()
    print("```text")
    print("G(T) = gcd(A_41(T), T^41-T) = " + uexpr(g, "T"))
    if roots:
        print("     = " + "*".join(f"(T-{x})" for x in roots))
    else:
        print("     = 1")
    print("```")
    print(f"Coefficient list: `{json.dumps(g, separators=(',', ':'))}`")
    print(f"Roots in F_41: `{roots}`")
    print(f"Central root value `-1/4`: `{central}`")
    root_info = [
        {
            "T": x,
            "half_grid_index": half_nodes.get(x),
            "central": x == central,
        }
        for x in roots
    ]
    print("Root classification:")
    print_json_block(root_info)
    print()

    print("# 4. Kernel values on zero pairs")
    print()
    pairs = []
    for ix, x in enumerate(roots):
        for y in roots[ix + 1:]:
            pairs.append({
                "T0": x,
                "T1": y,
                "K": beval(kernel, x, y),
                "M0": beval(t_minus_u, x, y),
                "pole_factor": beval(pole, x, y),
                "M_saturated": beval(saturated_m, x, y),
                "admissible_for_M": beval(saturated_m, x, y) != 0,
            })
    print_json_block(pairs)
    admissible = [q for q in pairs if q["admissible_for_M"]]
    failures = [q for q in admissible if q["K"] == 0]
    print(f"Admissible unordered pairs: `{len(admissible)}`; zero K-values: `{len(failures)}`.")
    print()

    print("# 5. Explicit Nullstellensatz certificate")
    print()
    print("Let `X(T)=T^41-T`.  Extended Euclid gives the exact univariate identity")
    print("`G(T)=s(T)A_41(T)+t(T)X(T)` with:")
    print()
    print("`s(T)` coefficient list:")
    print_json_block(bez_s)
    print("`t(T)` coefficient list:")
    print_json_block(bez_t)
    print(f"Checksums: `s={digest_u(bez_s)}`, `t={digest_u(bez_t)}`.")
    print()

    if not cert_sat["success"]:
        print("The requested saturated certificate does not exist for this M because")
        print("the following admissible common zero has K=0:")
        print_json_block(cert_sat["bad"])
    else:
        print("The computation finds exponent `N=1`.  Write the exact sparse polynomials")
        print("`B,Q_T,Q_U` below and define")
        print()
        print("```text")
        print("C(T,U)=Q_T(T,U)*s(T),")
        print("D(T,U)=Q_U(T,U)*s(U),")
        print("E(T,U)=Q_T(T,U)*t(T),")
        print("F(T,U)=Q_U(T,U)*t(U).")
        print("```")
        print()
        print("Then, coefficientwise in F_41[T,U],")
        print()
        print("```text")
        print("(T-U)(1+4T)(1+4U)")
        print(" = B(T,U) K_41(T,U)")
        print(" + C(T,U) A_41(T) + D(T,U) A_41(U)")
        print(" + E(T,U) (T^41-T) + F(T,U) (U^41-U).")
        print("```")
        print()
        print_poly_data("B", cert_sat["B"])
        print_poly_data("Q_T", cert_sat["QT"])
        print_poly_data("Q_U", cert_sat["QU"])
        print("### Expanded-cofactor audit")
        for name in ("C", "D", "E", "F"):
            f = cert_sat[name]
            print(f"- `{name}`: bidegree `{bdegrees(f)}`, terms `{len(f)}`, SHA-256 `{digest_b(f)}`")
        print()
        print("The expanded C,D,E,F triples are reproducibly generated by the code below;")
        print("their factorized definitions above are the exact cofactors and avoid duplicating")
        print("thousands of mechanically multiplied coefficients.")
        print()

    print("## Stronger T-U certificate")
    if cert_diag["success"]:
        print("`K_41` is already nonzero at every distinct pair of F_41-roots of A_41.")
        print("Hence the same construction gives N=1 for `M0=T-U`, before pole saturation.")
        print(f"Its B checksum is `{digest_b(cert_diag['B'])}`.")
        print("Multiplying that identity by `(1+4T)(1+4U)` is an independent derivation")
        print("of the displayed saturated certificate.")
    else:
        print("The stronger unsaturated statement fails exactly at:")
        print_json_block(cert_diag["bad"])
    print()

    print("# Verification assertions")
    checks = {
        "A_degree": len(a_poly) - 1,
        "K_bidegree": bdegrees(kernel),
        "A_lambda_equals_apery_all_half_grid": True,
        "G_degree_equals_number_of_F41_roots": len(g) - 1 == len(roots),
        "G_is_product_of_linear_root_factors": umonic(root_product) == g,
        "extended_euclid_identity": uadd(umul(bez_s, a_poly), umul(bez_t, node)) == g,
        "all_admissible_K_nonzero": not failures,
        "saturated_certificate_verified_coefficientwise": bool(cert_sat["success"]),
        "strong_T_minus_U_certificate_verified": bool(cert_diag["success"]),
    }
    # Directly verify A(lambda_n)=b_n using the binomial definition.
    from math import comb
    for n in range(MHALF + 1):
        lam = n * (n + 1) % P
        apery = sum(comb(n, k) ** 2 * comb(n + k, k) ** 2 for k in range(n + 1)) % P
        assert ueval(a_poly, lam) == apery
    print_json_block(checks)
    print()

    all_payload = (
        canon_u(a_poly) + b"\n" + canon_b(kernel) + b"\n" + canon_u(g)
        + b"\n" + (canon_b(cert_sat["B"]) if cert_sat["success"] else b"FAIL")
    )
    print(f"Combined core-data SHA-256: `{sha256(all_payload).hexdigest()}`")
    print()

    print("# Exact reproducer")
    print()
    print("```python")
    print(Path(__file__).read_text(encoding="utf-8").rstrip())
    print("```")


if __name__ == "__main__":
    main()
