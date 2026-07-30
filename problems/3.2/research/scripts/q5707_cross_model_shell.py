#!/usr/bin/env python3
"""Q5707: exact dependency-free cross-model Newton-carrier audit.

Models
------
Lambda = (1+x)(1+y)(1+z)((1+y)(1+z)+xyz)/(xyz)
Phi    = (x+y)(1+z)(x+y+z)(1+y+z)/(xyz)

For L^M = sum_v c_M(v) X^v, this script evaluates

    C_M^L(d) = sum_{v in d Z^3} c_M(v)

without expanding L^M.  It then forms the primitive Newton carrier

    G_L = sum_{i=0}^ell omega_i C_M^L(d0+i),
    omega_i = (-1)^i C(d0+i,i) C(d0+ell+1,ell-i).

Only the Python standard library is used.  All arithmetic is exact.
"""
from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import lru_cache
from fractions import Fraction
from math import comb, gcd, log, log10
from random import Random
from typing import Dict, Iterable, List, Tuple


# ---------------------------------------------------------------------------
# Small exact infrastructure
# ---------------------------------------------------------------------------

def prime_sieve(n: int) -> List[int]:
    mark = bytearray(b"\x01") * (n + 1)
    if n >= 0:
        mark[0] = 0
    if n >= 1:
        mark[1] = 0
    for p in range(2, int(n**0.5) + 1):
        if mark[p]:
            mark[p * p : n + 1 : p] = b"\x00" * (((n - p * p) // p) + 1)
    return [i for i in range(2, n + 1) if mark[i]]


def apery_numbers(n: int) -> List[int]:
    """b_0,...,b_n from the exact Apéry recurrence."""
    if n == 0:
        return [1]
    out = [1, 5]
    for m in range(1, n):
        num = (34 * m**3 + 51 * m**2 + 27 * m + 5) * out[m] - m**3 * out[m - 1]
        den = (m + 1) ** 3
        q, r = divmod(num, den)
        assert r == 0
        out.append(q)
    return out


def valuation(n: int, p: int) -> int:
    if n == 0:
        return 10**9
    n = abs(n)
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v


def decimal_height(n: int) -> Dict[str, object]:
    n = abs(n)
    if n == 0:
        return {"bits": 0, "digits": 1, "log10": float("-inf"), "ln": float("-inf")}
    s = str(n)
    head_len = min(16, len(s))
    mant = int(s[:head_len]) / 10 ** (head_len - 1)
    lg10 = len(s) - 1 + log10(mant)
    return {"bits": n.bit_length(), "digits": len(s), "log10": lg10, "ln": lg10 * log(10.0)}


def fmt_factorization(factors: Dict[int, int], cofactor: int, status: str) -> str:
    pieces = []
    for p in sorted(factors):
        e = factors[p]
        pieces.append(str(p) if e == 1 else f"{p}^{e}")
    if cofactor != 1:
        pieces.append(f"[{cofactor} ; {status}]")
    return " * ".join(pieces) if pieces else "1"


# ---------------------------------------------------------------------------
# Primality / partial factorization (stdlib only)
# ---------------------------------------------------------------------------

def is_probable_prime(n: int) -> bool:
    if n < 2:
        return False
    small = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47)
    for p in small:
        if n % p == 0:
            return n == p
    d = n - 1
    s = 0
    while d % 2 == 0:
        s += 1
        d //= 2
    # Deterministic for 64-bit; a strong probable-prime audit beyond it.
    bases = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71)
    for a in bases:
        if a >= n:
            continue
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def pollard_brent(n: int, seed: int, iteration_cap: int = 2_000_000) -> int | None:
    if n % 2 == 0:
        return 2
    rng = Random(seed)
    y = rng.randrange(1, n - 1)
    c = rng.randrange(1, n - 1)
    m = 128
    g = r = q = 1
    iterations = 0
    x = ys = 0
    while g == 1 and iterations < iteration_cap:
        x = y
        for _ in range(r):
            y = (y * y + c) % n
        k = 0
        while k < r and g == 1:
            ys = y
            for _ in range(min(m, r - k)):
                y = (y * y + c) % n
                q = q * abs(x - y) % n
                iterations += 1
            g = gcd(q, n)
            k += m
        r <<= 1
    if g == n:
        while iterations < iteration_cap:
            ys = (ys * ys + c) % n
            g = gcd(abs(x - ys), n)
            iterations += 1
            if g > 1:
                break
    return g if 1 < g < n else None


def factor_partial(n: int, trial_limit: int = 1_000_000) -> Tuple[Dict[int, int], int, str]:
    n = abs(n)
    factors: Dict[int, int] = {}
    if n in (0, 1):
        return factors, n, "unit" if n == 1 else "zero"
    for p in prime_sieve(min(trial_limit, int(n**0.5) + 1)):
        if p * p > n:
            break
        while n % p == 0:
            factors[p] = factors.get(p, 0) + 1
            n //= p
    if n == 1:
        return factors, 1, "fully factored"
    if is_probable_prime(n):
        factors[n] = factors.get(n, 0) + 1
        return factors, 1, "fully factored (large factors are strong probable primes)"

    # Pollard-Brent only on moderate cofactors; this is enough for the gcds
    # encountered here and avoids pretending to factor an arbitrary 1000-digit D.
    stack = [n]
    unresolved = 1
    while stack:
        m = stack.pop()
        if m == 1:
            continue
        if is_probable_prime(m):
            factors[m] = factors.get(m, 0) + 1
            continue
        if m.bit_length() > 320:
            unresolved *= m
            continue
        divisor = None
        for seed in range(1, 25):
            divisor = pollard_brent(m, seed)
            if divisor is not None:
                break
        if divisor is None:
            unresolved *= m
        else:
            stack.extend((divisor, m // divisor))
    status = "fully factored (large factors are strong probable primes)" if unresolved == 1 else "unresolved composite cofactor"
    return factors, unresolved, status


# ---------------------------------------------------------------------------
# Exact residue-binomial shell evaluators
# ---------------------------------------------------------------------------

def residue_binom_sum(n: int, residue: int, modulus: int) -> int:
    """sum_{0<=k<=n, k=residue mod modulus} C(n,k)."""
    residue %= modulus
    return sum(comb(n, k) for k in range(residue, n + 1, modulus))


def shell_lambda(M: int, d: int) -> int:
    """Exact C_M^Lambda(d), one outer sum.

    From the one-fold coefficient formula, summing u,v,w over dZ gives

      sum_t C(M,t) R_M(t;d) R_{2M-t}(M;d)^2.
    """
    total = 0
    for t in range(M + 1):
        a = residue_binom_sum(M, t, d)
        b = residue_binom_sum(2 * M - t, M, d)
        total += comb(M, t) * a * b * b
    return total


def shell_phi(M: int, d: int) -> int:
    """Exact C_M^Phi(d) as a finite residue-binomial sum.

    Expand (x+y)^M(1+z)^M(x+y+z)^M(1+y+z)^M/(xyz)^M.
    If R is the non-x degree in the third factor, i is the x degree in
    the first factor, beta is the y degree in that R-block, and delta is
    the y degree in the fourth factor, the three shell congruences are

      i = R                       (mod d),
      beta + delta = i            (mod d),
      k = M - R + beta            (mod d),

    where k is the combined z degree from (1+z)^(2M-delta).
    """
    total = 0
    # Cache the last residue sum for this (M,d); at d>M/2 each entry has
    # at most four binomial terms.
    zcache: Dict[Tuple[int, int], int] = {}

    def zsum(delta: int, residue: int) -> int:
        key = (delta, residue % d)
        val = zcache.get(key)
        if val is None:
            val = residue_binom_sum(2 * M - delta, key[1], d)
            zcache[key] = val
        return val

    binM = [comb(M, k) for k in range(M + 1)]
    for R in range(M + 1):
        cR = binM[R]
        for i in range(R % d, M + 1, d):
            pref = cR * binM[i]
            for beta in range(R + 1):
                inner = 0
                delta0 = (i - beta) % d
                zres = (M - R + beta) % d
                for delta in range(delta0, M + 1, d):
                    inner += binM[delta] * zsum(delta, zres)
                total += pref * comb(R, beta) * inner
    return total


# ---------------------------------------------------------------------------
# Independent tiny polynomial self-test
# ---------------------------------------------------------------------------
Mon = Tuple[int, int, int]
Poly = Dict[Mon, int]


def poly_mul(a: Poly, b: Poly) -> Poly:
    out: Poly = {}
    for (i, j, k), av in a.items():
        for (r, s, t), bv in b.items():
            key = (i + r, j + s, k + t)
            out[key] = out.get(key, 0) + av * bv
    return out


def poly_pow(base: Poly, e: int) -> Poly:
    out: Poly = {(0, 0, 0): 1}
    x = base
    while e:
        if e & 1:
            out = poly_mul(out, x)
        e >>= 1
        if e:
            x = poly_mul(x, x)
    return out


def shifted_product(factors: List[Poly], shift: Mon) -> Poly:
    out: Poly = {(0, 0, 0): 1}
    for f in factors:
        out = poly_mul(out, f)
    return {(a + shift[0], b + shift[1], c + shift[2]): v for (a, b, c), v in out.items()}


def base_models() -> Tuple[Poly, Poly]:
    one_x = {(0, 0, 0): 1, (1, 0, 0): 1}
    one_y = {(0, 0, 0): 1, (0, 1, 0): 1}
    one_z = {(0, 0, 0): 1, (0, 0, 1): 1}
    inner = {(0, 0, 0): 1, (0, 1, 0): 1, (0, 0, 1): 1, (0, 1, 1): 1, (1, 1, 1): 1}
    lam = shifted_product([one_x, one_y, one_z, inner], (-1, -1, -1))

    x_y = {(1, 0, 0): 1, (0, 1, 0): 1}
    xyz = {(1, 0, 0): 1, (0, 1, 0): 1, (0, 0, 1): 1}
    one_y_z = {(0, 0, 0): 1, (0, 1, 0): 1, (0, 0, 1): 1}
    phi = shifted_product([x_y, one_z, xyz, one_y_z], (-1, -1, -1))
    return lam, phi


def shell_from_poly_power(base: Poly, M: int, d: int) -> int:
    return sum(v for exp, v in poly_pow(base, M).items() if all(e % d == 0 for e in exp))


def eval_lambda_fraction(X: Fraction, Y: Fraction, Z: Fraction) -> Fraction:
    return (1 + X) * (1 + Y) * (1 + Z) * ((1 + Y) * (1 + Z) + X * Y * Z) / (X * Y * Z)


def eval_phi_fraction(x: Fraction, y: Fraction, z: Fraction) -> Fraction:
    return (x + y) * (1 + z) * (x + y + z) * (1 + y + z) / (x * y * z)


def self_test() -> None:
    lam, phi = base_models()
    bs = apery_numbers(5)
    for M in range(0, 5):
        assert poly_pow(lam, M).get((0, 0, 0), 0) == bs[M]
        assert poly_pow(phi, M).get((0, 0, 0), 0) == bs[M]
        for d in range(1, 6):
            assert shell_lambda(M, d) == shell_from_poly_power(lam, M, d)
            assert shell_phi(M, d) == shell_from_poly_power(phi, M, d)

    # Correct displayed-model mutation map: the compact Q3051 map followed
    # by the x<->z monomial chart.
    samples = [
        (Fraction(2), Fraction(3), Fraction(5)),
        (Fraction(3, 2), Fraction(5, 3), Fraction(7, 2)),
        (Fraction(-2), Fraction(1, 3), Fraction(4)),
    ]
    for x, y, z in samples:
        assert y + z != 0
        X = (x + y) / z
        Y = 1 / (y + z)
        Z = y / x
        assert eval_lambda_fraction(X, Y, Z) == eval_phi_fraction(x, y, z)


# ---------------------------------------------------------------------------
# Newton carrier and row audit
# ---------------------------------------------------------------------------

def weights(d0: int, ell: int) -> List[int]:
    out = [(-1) ** i * comb(d0 + i, i) * comb(d0 + ell + 1, ell - i) for i in range(ell + 1)]
    assert sum(out) == 1
    g = 0
    for w in out:
        g = gcd(g, abs(w))
    assert g == 1
    return out


def row_audit(n: int, ell: int = 50) -> Dict[str, object]:
    M = n - 1
    d0 = M // 2 + 1
    ell = min(ell, M - d0)
    ds = list(range(d0, d0 + ell + 1))
    ws = weights(d0, ell)
    bs = apery_numbers(M)

    all_primes = prime_sieve(n)
    top_targets = [p for p in all_primes if 2 * p > n and p < n and bs[n - p] % p == 0]
    candidate_primes = [d + 1 for d in ds if d + 1 in set(all_primes)]
    block_targets = [p for p in candidate_primes if bs[n - p] % p == 0]

    CL: List[int] = []
    CP: List[int] = []
    for d in ds:
        CL.append(shell_lambda(M, d))
        CP.append(shell_phi(M, d))

    GL = sum(w * c for w, c in zip(ws, CL))
    GP = sum(w * c for w, c in zip(ws, CP))
    D = GL - GP
    origin = bs[M]
    # Explicit common-origin cancellation check.
    assert D == sum(w * ((a - origin) - (b - origin)) for w, a, b in zip(ws, CL, CP))

    P = 1
    for p in candidate_primes:
        P *= p
    R = 1
    for p in block_targets:
        R *= p

    assert gcd(abs(GL), P) == R
    assert gcd(abs(GP), P) == R
    assert gcd(abs(D), P) == R

    return {
        "n": n,
        "M": M,
        "d0": d0,
        "ell": ell,
        "candidate_primes": candidate_primes,
        "top_targets": top_targets,
        "block_targets": block_targets,
        "candidate_primorial": P,
        "target_product": R,
        "G_lambda": GL,
        "G_phi": GP,
        "D": D,
        "gcd_models": gcd(abs(GL), abs(GP)),
        "gcd_GL_P": gcd(abs(GL), P),
        "gcd_GP_P": gcd(abs(GP), P),
        "gcd_D_P": gcd(abs(D), P),
        "shell_lambda_first_last": (CL[0], CL[-1]),
        "shell_phi_first_last": (CP[0], CP[-1]),
    }


def print_big(label: str, value: int, M: int, H: float) -> None:
    h = decimal_height(value)
    print(f"{label}={value}")
    print(
        f"{label}_bits={h['bits']} {label}_digits={h['digits']} "
        f"{label}_log10={h['log10']:.15f} {label}_ln={h['ln']:.15f} "
        f"{label}_ln_per_M={h['ln']/M:.15f} {label}_ln_per_H={h['ln']/H:.15f}"
    )


def emit_result(r: Dict[str, object]) -> None:
    n = int(r["n"])
    M = int(r["M"])
    H = n ** (1 / 3)
    print("=" * 100)
    print(f"ROW n={n} M={M} d0={r['d0']} ell={r['ell']} H=n^(1/3)={H:.15f}")
    print(f"candidate_primes={r['candidate_primes']}")
    print(f"all_top_half_targets={r['top_targets']}")
    print(f"block_targets={r['block_targets']}")
    print(f"candidate_primorial={r['candidate_primorial']}")
    print(f"target_product={r['target_product']}")
    print(f"shell_lambda_first_last={r['shell_lambda_first_last']}")
    print(f"shell_phi_first_last={r['shell_phi_first_last']}")

    for key in ("G_lambda", "G_phi", "D", "gcd_models"):
        print_big(key, int(r[key]), M, H)

    g = int(r["gcd_models"])
    fac, cof, status = factor_partial(g)
    print(f"gcd_models_factorization={fmt_factorization(fac, cof, status)}")
    print(f"gcd_models_factor_status={status}")

    # A deliberately bounded partial factorization of D.  The exact D is
    # printed above; this line never claims more than it proves.
    facD, cofD, statusD = factor_partial(int(r["D"]), trial_limit=100_000)
    print(f"D_partial_factorization={fmt_factorization(facD, cofD, statusD)}")
    print(f"D_factor_status={statusD}")

    print(f"gcd_G_lambda_candidate_primorial={r['gcd_GL_P']}")
    print(f"gcd_G_phi_candidate_primorial={r['gcd_GP_P']}")
    print(f"gcd_D_candidate_primorial={r['gcd_D_P']}")
    for p in (179, 193, 211):
        print(
            f"p={p} valuations: "
            f"v_p(G_lambda)={valuation(int(r['G_lambda']),p)} "
            f"v_p(G_phi)={valuation(int(r['G_phi']),p)} "
            f"v_p(D)={valuation(int(r['D']),p)} "
            f"v_p(gcd_models)={valuation(int(r['gcd_models']),p)}"
        )


def main() -> None:
    self_test()
    print("SELF_TEST=PASS")
    rows = (200, 272, 300, 321, 755)
    results: Dict[int, Dict[str, object]] = {}
    # Two GitHub-hosted cores: isolate each row so residue caches are released.
    with ProcessPoolExecutor(max_workers=2) as pool:
        jobs = {pool.submit(row_audit, n, 50): n for n in rows}
        for future in as_completed(jobs):
            n = jobs[future]
            results[n] = future.result()
            print(f"ROW_COMPLETE={n}", flush=True)
    for n in rows:
        emit_result(results[n])


if __name__ == "__main__":
    main()
