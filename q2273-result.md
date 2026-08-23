ANSWER Q2273 f6666eaf

# Exact p=41 result

The computation uses the monic Newton polynomials

```text
Psi_0(T)=1,
Psi_k(T)=product_{j=0}^{k-1}(T-j(j+1)).
```

Thus `Psi_k(lambda_n)=(k!)^2*binom(n,k)*binom(n+k,k)`.  Since
`(20!)^4 mod 41 = 1`, the requested normalization satisfies
`A_41(lambda_n)=b_n` and `K_41(lambda_r,lambda_s)=K(r,s)`.

Q2268 specifies the excluded product only schematically.  In the spectral
coordinate `lambda_n=n(n+1)`, diagonal and mate have the same divisor because
`lambda_r-lambda_s=(r-s)(r+s+1)`.  I therefore compute both the stronger
target `M0=T-U` and the canonical pole-saturated target
`M=(T-U)(1+4T)(1+4U)`; `1+4T=0` is the central ramification node.

# 1. Weights and A_41(T)

Weights `w_k=(20!/k!)^4 mod 41`, k=0,...,20:
```json
[1,1,18,23,31,40,18,40,31,31,23,16,23,37,4,16,10,23,4,18,1]
```

`deg A_41 = 40`; coefficient list `[T^0,...,T^40]`:
```json
[1,0,0,40,3,7,25,36,33,10,8,22,35,10,24,15,14,11,6,15,3,5,5,25,31,11,16,3,11,13,40,3,15,17,1,32,35,18,36,10,1]
```

Expanded polynomial:
```text
A_41(T) = T^40 + 10*T^39 + 36*T^38 + 18*T^37 + 35*T^36 + 32*T^35 + T^34 + 17*T^33 + 15*T^32 + 3*T^31 + 40*T^30 + 13*T^29 + 11*T^28 + 3*T^27 + 16*T^26 + 11*T^25 + 31*T^24 + 25*T^23 + 5*T^22 + 5*T^21 + 3*T^20 + 15*T^19 + 6*T^18 + 11*T^17 + 14*T^16 + 15*T^15 + 24*T^14 + 10*T^13 + 35*T^12 + 22*T^11 + 8*T^10 + 10*T^9 + 33*T^8 + 36*T^7 + 25*T^6 + 7*T^5 + 3*T^4 + 40*T^3 + 1
```
SHA-256: `c25e2648e656e04d626ee26b591b8e58b991f1af690dad7e65d126a89157b65d`

# 2. K_41(T,U)

Bidegree `(20, 20)`; nonzero terms `389`; SHA-256 `4a40d744ca274c2dc0629629a902b3f999a361ba97da4fe10ebb86ad79fa4dc7`.
The following 21x21 matrix has row j equal to coefficients of U^j and
column i equal to coefficients of T^i:
```json
[[1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,20,1,22,32,36,32,26,3,29,8,14,4,39,4,0,26,5,38,39],[0,20,1,2,1,32,1,23,20,13,27,1,27,32,7,9,1,20,23,32,6],[0,1,2,0,32,40,36,35,6,1,28,40,39,31,1,31,3,16,13,0,21],[0,22,1,32,10,2,21,37,3,34,36,9,6,7,23,36,2,26,19,6,29],[0,32,32,40,2,14,8,1,1,2,11,6,6,27,5,2,40,36,27,38,35],[0,36,1,36,21,8,37,9,19,23,15,28,6,2,5,27,25,38,39,28,22],[0,32,23,35,37,1,9,14,39,15,8,19,0,12,25,25,9,5,31,18,12],[0,26,20,6,3,1,19,39,10,8,7,40,36,1,38,38,4,10,6,20,37],[0,3,13,1,34,2,23,15,8,13,40,39,12,30,21,12,17,20,25,35,6],[0,29,27,28,36,11,15,8,7,40,26,13,10,32,18,1,9,2,39,8,10],[0,8,1,40,9,6,28,19,40,39,13,37,5,24,14,34,18,9,5,6,14],[0,14,27,39,6,6,6,0,36,12,10,5,21,13,24,23,34,34,12,5,1],[0,4,32,31,7,27,2,12,1,30,32,24,13,37,4,27,15,4,23,10,34],[0,39,7,1,23,5,5,25,38,21,18,14,24,4,18,0,9,11,0,39,27],[0,4,9,31,36,2,27,25,38,12,1,34,23,27,0,40,10,3,26,27,35],[0,0,1,3,2,40,25,9,4,17,9,18,34,15,9,10,28,12,4,34,13],[0,26,20,16,26,36,38,5,10,20,2,9,34,4,11,3,12,8,29,14,5],[0,5,23,13,19,27,39,31,6,25,39,5,12,23,0,26,4,29,22,4,17],[0,38,32,0,6,38,28,18,20,35,8,6,5,10,39,27,34,14,4,2,5],[0,39,6,21,29,35,22,12,37,6,10,14,1,34,27,35,13,5,17,5,1]]
```

# 3. Field-valued zero set

```text
G(T) = gcd(A_41(T), T^41-T) = T + 13
     = (T-28)
```
Coefficient list: `[13,1]`
Roots in F_41: `[28]`
Central root value `-1/4`: `10`
Root classification:
```json
[{"T":28,"half_grid_index":10,"central":false}]
```

# 4. Kernel values on zero pairs

```json
[]
```
Admissible unordered pairs: `0`; zero K-values: `0`.

# 5. Explicit Nullstellensatz certificate

Let `X(T)=T^41-T`.  Extended Euclid gives the exact univariate identity
`G(T)=s(T)A_41(T)+t(T)X(T)` with:

`s(T)` coefficient list:
```json
[13,37,21,7,17,30,3,0,9,12,9,21,26,7,7,39,24,8,33,37,28,1,33,37,21,1,16,11,25,15,16,34,9,7,34,5,28,12,19]
```
`t(T)` coefficient list:
```json
[36,21,35,19,6,28,27,28,40,22,40,25,28,9,34,14,33,7,18,9,3,36,12,19,7,2,21,2,37,33,25,12,32,36,7,29,3,22]
```
Checksums: `s=3e893e5fbc9bd3be3d011d8516ace2d7193f8e023b26baf7b09961b3ab014263`, `t=86b92624d8bbcda9da14f6b5dd9736aefa40f51e86b700b40ca48887a944fdc8`.

The computation finds exponent `N=1`.  Write the exact sparse polynomials
`B,Q_T,Q_U` below and define

```text
C(T,U)=Q_T(T,U)*s(T),
D(T,U)=Q_U(T,U)*s(U),
E(T,U)=Q_T(T,U)*t(T),
F(T,U)=Q_U(T,U)*t(U).
```

Then, coefficientwise in F_41[T,U],

```text
(T-U)(1+4T)(1+4U)
 = B(T,U) K_41(T,U)
 + C(T,U) A_41(T) + D(T,U) A_41(U)
 + E(T,U) (T^41-T) + F(T,U) (U^41-U).
```

### `B`
- bidegree: `(-1, -1)`
- nonzero terms: `0`
- SHA-256 of canonical sparse triples: `4f53cda18c2baa0c0354bb5f9a3ecbe5ed12ab4d8e11ba873c2f11161202b945`
- sparse triples `[i,j,c] = c*T^i*U^j`:
```json
[]
```
### `Q_T`
- bidegree: `(1, 2)`
- nonzero terms: `5`
- SHA-256 of canonical sparse triples: `64ee385e6a9a0f1058b85a4917ecdd2d7dea604af87b5e4dc43e348838c09b54`
- sparse triples `[i,j,c] = c*T^i*U^j`:
```json
[[0,0,31],[1,0,4],[0,1,38],[1,1,16],[0,2,25]]
```
### `Q_U`
- bidegree: `(0, 1)`
- nonzero terms: `2`
- SHA-256 of canonical sparse triples: `8b10b2b5e999751eca4c7f2d00dd0ad6f74097201f4b91782a84b681a01f4fcb`
- sparse triples `[i,j,c] = c*T^i*U^j`:
```json
[[0,0,10],[0,1,40]]
```
### Expanded-cofactor audit
- `C`: bidegree `(39, 2)`, terms `117`, SHA-256 `53b1a7cc6386a589bfbcd2e8e83d37aef1ab4148b7a162709afd3a391375157c`
- `D`: bidegree `(0, 39)`, terms `39`, SHA-256 `60b71884d2e6266cba22ea10699d22f6c1a945424a7b1131e73f390f3c0cba00`
- `E`: bidegree `(38, 2)`, terms `114`, SHA-256 `4a0ac33275adcdec789f3fa4e68a2e50a9e603acd4f88729bf8385848a16da99`
- `F`: bidegree `(0, 38)`, terms `37`, SHA-256 `123ffee19a632bd7e2ea81efa20a690909817faf59bf16f95ac33da8785092e2`

The expanded C,D,E,F triples are reproducibly generated by the code below;
their factorized definitions above are the exact cofactors and avoid duplicating
thousands of mechanically multiplied coefficients.

## Stronger T-U certificate
`K_41` is already nonzero at every distinct pair of F_41-roots of A_41.
Hence the same construction gives N=1 for `M0=T-U`, before pole saturation.
Its B checksum is `4f53cda18c2baa0c0354bb5f9a3ecbe5ed12ab4d8e11ba873c2f11161202b945`.
Multiplying that identity by `(1+4T)(1+4U)` is an independent derivation
of the displayed saturated certificate.

# Verification assertions
```json
{"A_degree":40,"K_bidegree":[20,20],"A_lambda_equals_apery_all_half_grid":true,"G_degree_equals_number_of_F41_roots":true,"G_is_product_of_linear_root_factors":true,"extended_euclid_identity":true,"all_admissible_K_nonzero":true,"saturated_certificate_verified_coefficientwise":true,"strong_T_minus_U_certificate_verified":true}
```

Combined core-data SHA-256: `9fac8e29d380445895dd5ecdde1c8def6a2daaa5a668244285d3ad4be326367d`

# Exact reproducer

```python
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
```
