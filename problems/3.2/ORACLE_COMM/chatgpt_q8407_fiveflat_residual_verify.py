#!/usr/bin/env python3
"""Exact stdlib verifier for the Q8407 five-flat normalization audit.

This verifier deliberately does NOT manufacture the untracked 22-ray Laurent
polynomial Lambda or infer an all-h theorem from finite computation.  It checks
only exact algebraic consequences of supplied rho_{h,s} data and the structural
normalizations used in chatgpt_q8407_fiveflat_residual_report.md:

* half-integer binomial coefficients are dyadic;
* Y_h(x) and the five F_h coefficient gaps agree with the exact triangular map;
* the shifted five roots are -2,-1,0,1,2 and q(z)=z(z^2-1)(z^2-4);
* exact (or mod-P) five-flatness is equivalent to divisibility by q;
* parity and r=h+1 indexing formulas are consistent;
* the abstract rho -> five-gap map has a nonzero integral kernel for h>=3.

Input JSON schema:
    {"h": 7, "rho": ["1", "-3/2", ..., "0"]}
with exactly 2*h+1 entries.  Integers, JSON integers, and strings "a/b" are
accepted.  Floats are rejected to preserve exactness.

Examples:
    python3 chatgpt_q8407_fiveflat_residual_verify.py --self-test
    python3 chatgpt_q8407_fiveflat_residual_verify.py \
        --rho-json exact_rho.json --prime 101 --assert-target

--assert-target is a candidate-counterexample guard only: if P>2h+2 and all
five supplied values vanish mod P, it exits nonzero and reports that the input
would contradict the target lemma.  Passing finite data is not treated as a
proof of the target.
"""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Iterable, Sequence


# ---------------------------------------------------------------------------
# Exact scalar arithmetic
# ---------------------------------------------------------------------------


def parse_fraction(value: object) -> Fraction:
    if isinstance(value, bool):
        raise TypeError("booleans are not valid exact rho entries")
    if isinstance(value, int):
        return Fraction(value, 1)
    if isinstance(value, float):
        raise TypeError("floating-point rho entries are forbidden; use strings a/b")
    if isinstance(value, str):
        text = value.strip()
        if not text:
            raise ValueError("empty rational string")
        return Fraction(text)
    raise TypeError(f"unsupported exact rational value: {value!r}")


def generalized_binom(x: Fraction, n: int) -> Fraction:
    if n < 0:
        return Fraction(0, 1)
    out = Fraction(1, 1)
    for j in range(n):
        out *= x - j
        out /= j + 1
    return out


def is_power_of_two(n: int) -> bool:
    return n > 0 and (n & (n - 1)) == 0


def lambda_h(h: int) -> Fraction:
    if h < 0:
        raise ValueError("h must be nonnegative")
    return Fraction(h // 2, 1) - Fraction(1, 2)


def alpha_h(h: int) -> Fraction:
    return lambda_h(h) - 4


def assert_dyadic_half_binomials(h: int) -> None:
    # Every binomial coefficient actually used in the five-gap/Y transforms.
    lam = lambda_h(h)
    alpha = alpha_h(h)
    for x in [lam - j for j in range(5)] + [alpha]:
        for s in range(2 * h + 1):
            denominator = generalized_binom(x, s).denominator
            assert is_power_of_two(denominator), (h, x, s, denominator)


# ---------------------------------------------------------------------------
# Y_h and F_h coefficient gaps
# ---------------------------------------------------------------------------


def check_rho_shape(h: int, rho: Sequence[Fraction]) -> None:
    if h < 1:
        raise ValueError("Q8407 normalization assumes h>=1")
    if len(rho) != 2 * h + 1:
        raise ValueError(f"expected {2*h+1} rho entries for h={h}, got {len(rho)}")


def y_value(h: int, rho: Sequence[Fraction], x: Fraction) -> Fraction:
    check_rho_shape(h, rho)
    return sum(
        ((-1) ** s) * rho[s] * generalized_binom(x, s)
        for s in range(2 * h + 1)
    )


def gap_c(h: int, rho: Sequence[Fraction], k: int) -> Fraction:
    """C_{h,k} = [u^(h-k)] (1-u)^(lambda_h-4) R_h(u)."""
    check_rho_shape(h, rho)
    if not 0 <= k <= 4:
        raise ValueError("gap index k must be 0..4")
    a = alpha_h(h)
    return sum(
        rho[s] * ((-1) ** (s - k)) * generalized_binom(a, s - k)
        for s in range(k, 2 * h + 1)
    )


def five_y_values(h: int, rho: Sequence[Fraction]) -> list[Fraction]:
    lam = lambda_h(h)
    return [y_value(h, rho, lam - j) for j in range(5)]


def five_gaps(h: int, rho: Sequence[Fraction]) -> list[Fraction]:
    return [gap_c(h, rho, k) for k in range(5)]


def y_from_gaps(gaps: Sequence[Fraction]) -> list[Fraction]:
    if len(gaps) != 5:
        raise ValueError("need exactly five gaps")
    # Output order A_j = Y(lambda-j), j=0..4.
    out: list[Fraction] = []
    for j in range(5):
        m = 4 - j
        out.append(sum(((-1) ** k) * math.comb(m, k) * gaps[k] for k in range(m + 1)))
    return out


def gaps_from_y(values: Sequence[Fraction]) -> list[Fraction]:
    if len(values) != 5:
        raise ValueError("need exactly five Y-values")
    # C_k = sum_{m=0}^k (-1)^m binom(k,m) A_{4-m}.
    out: list[Fraction] = []
    for k in range(5):
        out.append(
            sum(((-1) ** m) * math.comb(k, m) * values[4 - m] for m in range(k + 1))
        )
    return out


def verify_triangular_transform(h: int, rho: Sequence[Fraction]) -> None:
    direct_y = five_y_values(h, rho)
    gaps = five_gaps(h, rho)
    assert y_from_gaps(gaps) == direct_y
    assert gaps_from_y(direct_y) == gaps
    assert (all(x == 0 for x in direct_y)) == (all(x == 0 for x in gaps))


# ---------------------------------------------------------------------------
# Exact polynomials in z, coefficients stored low degree first
# ---------------------------------------------------------------------------


def poly_strip(a: list[Fraction]) -> list[Fraction]:
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a


def poly_add(a: Sequence[Fraction], b: Sequence[Fraction]) -> list[Fraction]:
    n = max(len(a), len(b))
    out = [Fraction(0, 1) for _ in range(n)]
    for i, v in enumerate(a):
        out[i] += v
    for i, v in enumerate(b):
        out[i] += v
    return poly_strip(out)


def poly_scale(a: Sequence[Fraction], c: Fraction) -> list[Fraction]:
    return poly_strip([c * v for v in a])


def poly_mul(a: Sequence[Fraction], b: Sequence[Fraction]) -> list[Fraction]:
    out = [Fraction(0, 1) for _ in range(len(a) + len(b) - 1)]
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return poly_strip(out)


def binom_poly(shift: Fraction, n: int) -> list[Fraction]:
    """Polynomial binom(z+shift,n), low-degree coefficients first."""
    if n < 0:
        return [Fraction(0, 1)]
    out = [Fraction(1, 1)]
    for j in range(n):
        out = poly_mul(out, [shift - j, Fraction(1, 1)])
        out = poly_scale(out, Fraction(1, j + 1))
    return out


def shifted_y_poly(h: int, rho: Sequence[Fraction]) -> list[Fraction]:
    """Y_h(lambda_h-2+z) over Q."""
    check_rho_shape(h, rho)
    shift = lambda_h(h) - 2
    out = [Fraction(0, 1)]
    for s, r in enumerate(rho):
        out = poly_add(out, poly_scale(binom_poly(shift, s), ((-1) ** s) * r))
    return poly_strip(out)


def q_poly() -> list[Fraction]:
    # z^5 - 5 z^3 + 4 z
    return [Fraction(0), Fraction(4), Fraction(0), Fraction(-5), Fraction(0), Fraction(1)]


def q_poly_from_roots() -> list[Fraction]:
    out = [Fraction(1)]
    for root in (-2, -1, 0, 1, 2):
        out = poly_mul(out, [Fraction(-root), Fraction(1)])
    return out


def poly_divmod_monic(
    numerator: Sequence[Fraction], denominator: Sequence[Fraction]
) -> tuple[list[Fraction], list[Fraction]]:
    den = poly_strip(list(denominator))
    if not den or den[-1] != 1:
        raise ValueError("denominator must be monic")
    rem = poly_strip(list(numerator))
    if len(rem) < len(den):
        return [Fraction(0)], rem
    quotient = [Fraction(0) for _ in range(len(rem) - len(den) + 1)]
    while not (len(rem) == 1 and rem[0] == 0) and len(rem) >= len(den):
        shift = len(rem) - len(den)
        coeff = rem[-1]
        quotient[shift] += coeff
        for i, d in enumerate(den):
            rem[i + shift] -= coeff * d
        rem = poly_strip(rem)
    return poly_strip(quotient), poly_strip(rem)


def poly_eval(a: Sequence[Fraction], z: Fraction) -> Fraction:
    out = Fraction(0)
    for coeff in reversed(a):
        out = out * z + coeff
    return out


def verify_centered_q_exact(h: int, rho: Sequence[Fraction]) -> None:
    q = q_poly()
    assert q_poly_from_roots() == q
    # q(z)=120*binom(z+2,5) over Q.
    assert poly_scale(binom_poly(Fraction(2), 5), Fraction(120)) == q

    yshift = shifted_y_poly(h, rho)
    direct = five_y_values(h, rho)
    # z=2,1,0,-1,-2 correspond to j=0,1,2,3,4.
    evals = [poly_eval(yshift, Fraction(z)) for z in (2, 1, 0, -1, -2)]
    assert evals == direct
    _, rem = poly_divmod_monic(yshift, q)
    flat = all(v == 0 for v in direct)
    divisible = all(v == 0 for v in rem)
    assert flat == divisible


# ---------------------------------------------------------------------------
# Modular arithmetic for P>2h+2
# ---------------------------------------------------------------------------


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def fraction_mod(x: Fraction, p: int) -> int:
    den = x.denominator % p
    if den == 0:
        raise ZeroDivisionError(f"denominator {x.denominator} is 0 mod {p}")
    return (x.numerator % p) * pow(den, -1, p) % p


def poly_mod(a: Sequence[Fraction], p: int) -> list[int]:
    out = [fraction_mod(x, p) for x in a]
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def poly_rem_mod(numerator: Sequence[int], denominator: Sequence[int], p: int) -> list[int]:
    den = [x % p for x in denominator]
    while len(den) > 1 and den[-1] == 0:
        den.pop()
    if not den or den[-1] % p != 1:
        raise ValueError("modular denominator must be monic")
    rem = [x % p for x in numerator]
    while len(rem) > 1 and rem[-1] == 0:
        rem.pop()
    while not (len(rem) == 1 and rem[0] == 0) and len(rem) >= len(den):
        shift = len(rem) - len(den)
        coeff = rem[-1] % p
        for i, d in enumerate(den):
            rem[i + shift] = (rem[i + shift] - coeff * d) % p
        while len(rem) > 1 and rem[-1] == 0:
            rem.pop()
    return rem


def verify_mod_prime(h: int, rho: Sequence[Fraction], p: int) -> dict[str, object]:
    if not is_prime(p) or p % 2 == 0:
        raise ValueError(f"P={p} must be an odd prime")
    if p <= 2 * h + 2:
        raise ValueError(f"this verifier's high-prime path requires P>2h+2={2*h+2}")

    ys = [fraction_mod(v, p) for v in five_y_values(h, rho)]
    cs = [fraction_mod(v, p) for v in five_gaps(h, rho)]
    flat_y = all(v == 0 for v in ys)
    flat_c = all(v == 0 for v in cs)
    assert flat_y == flat_c

    ypoly_mod = poly_mod(shifted_y_poly(h, rho), p)
    qmod = [int(x) % p for x in q_poly()]
    rem = poly_rem_mod(ypoly_mod, qmod, p)
    q_divides = all(v % p == 0 for v in rem)
    assert q_divides == flat_y

    # Since P>2h+2, P is also coprime to h; this is the terminal contradiction
    # that a future [TGD-5] Bezout identity would exploit.
    assert h % p != 0

    return {
        "P": p,
        "Y_mod_P": ys,
        "C_mod_P": cs,
        "q_remainder_mod_P": rem,
        "five_flat_mod_P": flat_y,
    }


# ---------------------------------------------------------------------------
# Generic nonzero integral kernel: normalization alone cannot prove target
# ---------------------------------------------------------------------------


def lcm(a: int, b: int) -> int:
    return abs(a * b) // math.gcd(a, b) if a and b else 0


def synthetic_integral_kernel(h: int) -> list[Fraction]:
    if h < 3:
        raise ValueError("explicit rho_5 kernel construction needs h>=3")
    rho = [Fraction(0) for _ in range(2 * h + 1)]
    rho[5] = Fraction(1)
    a = alpha_h(h)

    # Solve C_4=...=C_0=0.  The coefficient of rho[k] in C_k is 1.
    for k in range(4, -1, -1):
        tail = sum(
            rho[s] * ((-1) ** (s - k)) * generalized_binom(a, s - k)
            for s in range(k + 1, 2 * h + 1)
        )
        rho[k] = -tail

    assert any(x != 0 for x in rho)
    assert all(x == 0 for x in five_gaps(h, rho))
    assert all(x == 0 for x in five_y_values(h, rho))

    scale = 1
    for x in rho:
        scale = lcm(scale, x.denominator)
    integer_rho = [x * scale for x in rho]
    assert all(x.denominator == 1 for x in integer_rho)
    assert any(x != 0 for x in integer_rho)
    assert all(x == 0 for x in five_gaps(h, integer_rho))
    assert all(x == 0 for x in five_y_values(h, integer_rho))
    return integer_rho


# ---------------------------------------------------------------------------
# Index/parity audit
# ---------------------------------------------------------------------------


def verify_parity_and_r_index(h: int) -> None:
    r = h + 1
    lam = lambda_h(h)
    a = alpha_h(h)
    assert h == r - 1
    assert 2 * h + 2 == 2 * r
    if h % 2 == 0:
        assert r % 2 == 1
        assert lam == Fraction(r - 2, 2)
        assert a == Fraction(r - 10, 2)
    else:
        assert r % 2 == 0
        assert lam == Fraction(r - 3, 2)
        assert a == Fraction(r - 11, 2)
    assert lambda_h(h + 2) == lam + 1
    assert alpha_h(h + 2) == a + 1


# ---------------------------------------------------------------------------
# Deterministic structural self-test (not theorem evidence)
# ---------------------------------------------------------------------------


def deterministic_rho(h: int) -> list[Fraction]:
    # Non-random exact values chosen only to exercise the identities.
    return [Fraction((s + 1) * (s + h + 2) * ((-1) ** (s % 3)), 1) for s in range(2 * h + 1)]


def run_self_test() -> None:
    assert q_poly_from_roots() == q_poly()
    assert poly_scale(binom_poly(Fraction(2), 5), Fraction(120)) == q_poly()

    for h in range(1, 13):
        verify_parity_and_r_index(h)
        assert_dyadic_half_binomials(h)
        rho = deterministic_rho(h)
        verify_triangular_transform(h, rho)
        verify_centered_q_exact(h, rho)

    for h in range(3, 13):
        kernel = synthetic_integral_kernel(h)
        assert all(x.denominator == 1 for x in kernel)
        verify_triangular_transform(h, kernel)
        verify_centered_q_exact(h, kernel)

    print("PASS: half-integer binomial denominators are dyadic")
    print("PASS: five Y-values <-> five F-gap coefficients are unitriangularly equivalent")
    print("PASS: centered q(z)=z(z^2-1)(z^2-4)=120*binom(z+2,5)")
    print("PASS: exact five-root vanishing <-> q-divisibility")
    print("PASS: parity and r=h+1 indexing identities")
    print("PASS: nonzero integral synthetic five-flat kernel exists for h>=3")
    print("Q8407 STRUCTURAL VERIFIER: PASS (sanity checks only; not an all-h proof)")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def load_rho_json(path: Path) -> tuple[int, list[Fraction]]:
    raw = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(raw, dict):
        raise TypeError("rho JSON must be an object")
    if "h" not in raw or "rho" not in raw:
        raise ValueError("rho JSON must contain keys 'h' and 'rho'")
    h_raw = raw["h"]
    if isinstance(h_raw, bool) or not isinstance(h_raw, int):
        raise TypeError("h must be an integer")
    rho_raw = raw["rho"]
    if not isinstance(rho_raw, list):
        raise TypeError("rho must be a JSON list")
    rho = [parse_fraction(v) for v in rho_raw]
    check_rho_shape(h_raw, rho)
    return h_raw, rho


def fmt_fraction(x: Fraction) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true", help="run exact structural sanity checks")
    parser.add_argument("--rho-json", type=Path, help="JSON file containing exact h and rho entries")
    parser.add_argument("--prime", type=int, help="odd prime P>2h+2 for modular checks")
    parser.add_argument(
        "--assert-target",
        action="store_true",
        help="fail if supplied exact data is five-flat mod the supplied high prime",
    )
    args = parser.parse_args()

    if not args.self_test and args.rho_json is None:
        parser.error("choose --self-test and/or --rho-json")
    if args.assert_target and (args.rho_json is None or args.prime is None):
        parser.error("--assert-target requires --rho-json and --prime")

    if args.self_test:
        run_self_test()

    if args.rho_json is not None:
        h, rho = load_rho_json(args.rho_json)
        verify_parity_and_r_index(h)
        assert_dyadic_half_binomials(h)
        verify_triangular_transform(h, rho)
        verify_centered_q_exact(h, rho)

        ys = five_y_values(h, rho)
        cs = five_gaps(h, rho)
        _, qrem = poly_divmod_monic(shifted_y_poly(h, rho), q_poly())

        print(f"h={h}")
        print(f"r=h+1={h+1}")
        print(f"lambda_h={fmt_fraction(lambda_h(h))}")
        print(f"alpha_h=lambda_h-4={fmt_fraction(alpha_h(h))}")
        print("Y(lambda-j), j=0..4 =", [fmt_fraction(x) for x in ys])
        print("C_[h-k], k=0..4 =", [fmt_fraction(x) for x in cs])
        print("shifted Y remainder mod q over Q =", [fmt_fraction(x) for x in qrem])
        print("exact five-flat =", all(x == 0 for x in ys))

        if args.prime is not None:
            info = verify_mod_prime(h, rho, args.prime)
            print("mod-prime audit =", json.dumps(info, sort_keys=True))
            if args.assert_target and bool(info["five_flat_mod_P"]):
                raise SystemExit(
                    "CANDIDATE COUNTEREXAMPLE: supplied exact rho is five-flat at "
                    f"P={args.prime}>2h+2. Verify the literal Lambda/rho normalization."
                )
            if args.assert_target:
                print("PASS: supplied row is not five-flat at the requested high prime")


if __name__ == "__main__":
    main()
