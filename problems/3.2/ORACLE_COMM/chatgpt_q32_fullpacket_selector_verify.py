#!/usr/bin/env python3
"""Exact verifier for chatgpt_q32_fullpacket_selector.md.

Standard library only.  This checks finite integer identities and locked actual
Apéry regressions.  It is not an asymptotic proof and does not manufacture a
sixteen-node witness when none is supplied.

Optional JSON format:

    {
      "m": 123,
      "g": 2,
      "branch": "direct",
      "h": [ ... sixteen folded rows ... ],
      "p": [ ... optional sixteen primes ... ],
      "require_cross_unit": true
    }

For a Boolean packet one may instead give

    {
      "m": 123,
      "g": 2,
      "branch": "direct",
      "a0": 10,
      "boolean_gaps": [1, 3, 7, 15],
      "require_cross_unit": true
    }

The `h` values are always physical folded rows.  `a0` and `boolean_gaps` are
in the normalized a-coordinate, so h=g*a.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from itertools import product
from math import gcd, prod
from pathlib import Path
from typing import Any, Iterable, Sequence


def P(n: int) -> int:
    return 34 * n**3 + 51 * n**2 + 27 * n + 5


def apery_up_to(N: int) -> list[int]:
    """Return exact Apéry integers b_0,...,b_N."""
    if N < 0:
        raise ValueError("N must be nonnegative")
    if N == 0:
        return [1]
    b = [1, 5]
    for n in range(1, N):
        numerator = P(n) * b[n] - n**3 * b[n - 1]
        denominator = (n + 1) ** 3
        q, r = divmod(numerator, denominator)
        assert r == 0, (n, numerator, denominator)
        b.append(q)
    return b


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


def vp(n: int, p: int) -> int:
    if n == 0:
        raise ValueError("vp(0,p) is not used in this verifier")
    n = abs(n)
    out = 0
    while n % p == 0:
        n //= p
        out += 1
    return out


def branch_prime(m: int, h: int, branch: str) -> int:
    if branch == "direct":
        return m - h
    if branch == "reflected":
        numerator = m + 1 + h
        assert numerator % 2 == 0, (m, h, numerator)
        return numerator // 2
    raise ValueError(f"unknown branch {branch!r}")


def boolean_rows(g: int, a0: int, gaps: Sequence[int]) -> list[int]:
    assert len(gaps) == 4
    assert g > 0 and a0 >= 0 and all(d > 0 for d in gaps)
    out: list[int] = []
    for mask in range(16):
        a = a0 + sum(gaps[r] for r in range(4) if (mask >> r) & 1)
        out.append(g * a)
    assert len(set(out)) == 16, "Boolean subset sums must be distinct"
    return out


def vandermonde_weights(primes: Sequence[int]) -> list[int]:
    return [
        prod(primes[j] - primes[i] for j in range(len(primes)) if j != i)
        for i in range(len(primes))
    ]


def verify_fixed_g_weights(
    g: int, branch: str, hs: Sequence[int], primes: Sequence[int]
) -> list[int]:
    assert all(h % g == 0 for h in hs)
    a = [h // g for h in hs]
    V = vandermonde_weights(primes)
    exponent = len(primes) - 1
    for i in range(len(primes)):
        normalized = prod(a[j] - a[i] for j in range(len(primes)) if j != i)
        if branch == "direct":
            assert V[i] == (-g) ** exponent * normalized
        elif branch == "reflected":
            assert 2**exponent * V[i] == g**exponent * normalized
        else:
            raise ValueError(branch)
    return V


def verify_packet(
    *,
    m: int,
    g: int,
    branch: str,
    hs: Sequence[int],
    expected_primes: Sequence[int] | None = None,
    require_cross_unit: bool = True,
) -> dict[str, Any]:
    """Verify one actual direct/reflected packet and selector/BRJ identities."""
    assert m >= 1 and g >= 1 and (m - 1) % g == 0
    hs = list(hs)
    assert hs and len(set(hs)) == len(hs)
    assert all(h >= 0 and h % g == 0 for h in hs)

    primes = [branch_prime(m, h, branch) for h in hs]
    if expected_primes is not None:
        assert primes == list(expected_primes), (primes, expected_primes)
    assert len(set(primes)) == len(primes)
    assert all(p > 2 and p % 2 == 1 and is_prime(p) for p in primes)

    b = apery_up_to(max(m, max(hs)))
    for i, (h, p) in enumerate(zip(hs, primes)):
        assert b[h] % p == 0, ("target", i, h, p)
        assert b[m] % p == 0, ("common physical row", i, m, p)
        assert gcd(g, p) == 1
        if branch == "direct":
            assert p == m - h
            assert p % g == 1 % g
        else:
            assert 2 * p == m + 1 + h
            assert (2 * p - 2) % g == 0

    if require_cross_unit:
        for i, p in enumerate(primes):
            for j, h in enumerate(hs):
                if i != j:
                    assert b[h] % p != 0, ("cross-unit", i, j, p, h)

    R = prod(primes)
    assert b[m] % R == 0
    B = b[m] // R
    q = [b[h] // p for h, p in zip(hs, primes)]
    assert all(value > 0 for value in q)

    C = sum((R // p) * b[h] for h, p in zip(hs, primes))
    Q = sum(q)
    assert C == R * Q
    assert Q >= len(primes)
    for p in primes:
        assert C % p == 0
        assert vp(C, p) == 1 + vp(Q, p)

    # Exact cancellation criterion for a second p-factor in C.
    for i, p in enumerate(primes):
        criterion = q[i]
        for j in range(len(primes)):
            if i != j:
                criterion += (b[hs[j]] % p) * pow(primes[j], -1, p)
        assert (Q % p == 0) == (criterion % p == 0)

    V = verify_fixed_g_weights(g, branch, hs, primes)
    z = q
    gamma_hat: list[int] = []
    L: list[int] = []
    E: list[int] = []
    for i, (h, p) in enumerate(zip(hs, primes)):
        difference = b[m] - 5 * b[h]
        assert difference % p == 0
        gamma = difference // p
        gamma_hat.append(gamma)
        Li = 5 * z[i] + gamma
        L.append(Li)
        assert Li == b[m] // p == B * (R // p)
        Ei = B * V[i] - Li
        E.append(Ei)
        assert Ei == B * (V[i] - R // p)
        assert Ei % p == 0

    # The first target-divisible nonlinear carrier and first B-eliminating minor.
    for i in range(len(primes)):
        for j in range(i + 1, len(primes)):
            assert (E[i] * E[j]) % (primes[i] * primes[j]) == 0
            minor = V[j] * L[i] - V[i] * L[j]
            eliminated = V[j] * E[i] - V[i] * E[j]
            assert eliminated == -minor
            assert minor == B * (V[j] * (R // primes[i]) - V[i] * (R // primes[j]))

    H = max(hs)
    index_H = hs.index(H)
    p_H = primes[index_H]
    assert Fraction(b[H], max(primes)) <= Q
    assert Q <= Fraction(len(primes) * b[H], min(primes))
    assert q[index_H] == b[H] // p_H

    return {
        "m": m,
        "g": g,
        "branch": branch,
        "h": hs,
        "p": primes,
        "R": R,
        "B": B,
        "C": C,
        "Q": Q,
        "q": q,
        "V": V,
        "L": L,
        "E": E,
        "b": b,
    }


def verify_boolean_packet(
    *,
    m: int,
    g: int,
    branch: str,
    a0: int,
    gaps: Sequence[int],
    expected_primes: Sequence[int] | None = None,
    require_cross_unit: bool = True,
) -> dict[str, Any]:
    hs = boolean_rows(g, a0, gaps)
    data = verify_packet(
        m=m,
        g=g,
        branch=branch,
        hs=hs,
        expected_primes=expected_primes,
        require_cross_unit=require_cross_unit,
    )
    b: list[int] = data["b"]
    primes: list[int] = data["p"]
    q: list[int] = data["q"]
    R: int = data["R"]

    signs = [(-1) ** (4 - mask.bit_count()) for mask in range(16)]
    F = sum(signs[mask] * q[mask] for mask in range(16))
    H = hs[15]
    assert H == max(hs)
    qH = q[15]

    # Independent geometric domination using every intermediate Apéry row.
    if branch == "direct":
        all_q = [Fraction(b[n], m - n) for n in range(H + 1)]
        assert all(all_q[n + 1] >= 5 * all_q[n] for n in range(H))
        assert sum(all_q[:-1], Fraction(0)) <= all_q[-1] / 4
        assert Fraction(F) >= Fraction(3, 4) * qH > 0
    else:
        all_q = [Fraction(2 * b[n], m + 1 + n) for n in range(H + 1)]
        assert all(all_q[n + 1] >= Fraction(5, 2) * all_q[n] for n in range(H))
        assert sum(all_q[:-1], Fraction(0)) <= Fraction(2, 3) * all_q[-1]
        assert Fraction(F) >= Fraction(1, 3) * qH > 0

    # Four physical edge gaps.  In the small-gap shell regime these are R-units.
    physical_edge_denominator = prod(g * d for d in gaps)
    if gcd(R, physical_edge_denominator) == 1:
        numerator_pre_R = abs(R * F) // gcd(abs(R * F), physical_edge_denominator)
        numerator_post_R = abs(F) // gcd(abs(F), physical_edge_denominator)
        assert numerator_pre_R % R == 0
        assert numerator_pre_R == R * numerator_post_R
        assert numerator_pre_R >= R * abs(F) // physical_edge_denominator

    vandermonde = prod(abs(hs[j] - hs[i]) for i in range(16) for j in range(i + 1, 16))
    diameter = max(hs) - min(hs)
    assert vandermonde <= diameter**120
    if gcd(R, vandermonde) == 1:
        numerator_vand = abs(R * F) // gcd(abs(R * F), vandermonde)
        assert numerator_vand % R == 0

    data.update(
        {
            "a0": a0,
            "boolean_gaps": list(gaps),
            "F_square": F,
            "physical_edge_denominator": physical_edge_denominator,
            "vandermonde": vandermonde,
        }
    )
    return data


def universal_selector_holds(primes: Sequence[int], coefficients: Sequence[int]) -> bool:
    """Exhaust the foreign unit residues for a small finite prime tuple."""
    assert len(primes) == len(coefficients)
    for j, p in enumerate(primes):
        foreign = [i for i in range(len(primes)) if i != j]
        for residues in product(range(1, p), repeat=len(foreign)):
            total = 0
            for i, residue in zip(foreign, residues):
                total += coefficients[i] * residue
            if total % p:
                return False
    return True


def verify_linear_saturation() -> None:
    primes = [5, 7, 11]
    R = prod(primes)
    coefficients = [3 * (R // p) for p in primes]
    assert universal_selector_holds(primes, coefficients)
    assert prod(R // p for p in primes) == R ** (len(primes) - 1)

    for i in range(len(primes)):
        for j in range(len(primes)):
            if i != j:
                assert (R // primes[i]) % primes[j] == 0
        assert gcd(R // primes[i], primes[i]) == 1

    bad = coefficients.copy()
    bad[0] //= primes[1]
    assert not universal_selector_holds(primes, bad)


def verify_growth() -> None:
    b = apery_up_to(100)
    for n in range(1, 100):
        assert b[n] > b[n - 1]
        assert b[n + 1] >= 5 * b[n]
        assert b[n] >= (4**n // (2 * n + 1)) ** 2
        assert b[n] <= (n + 1) * 64**n


def verify_actual_regressions() -> tuple[dict[str, Any], dict[str, Any]]:
    # Actual direct defining-characteristic zero from the same-project Q8318 audit.
    direct = verify_packet(
        m=39,
        g=2,
        branch="direct",
        hs=[8],
        expected_primes=[31],
        require_cross_unit=True,
    )
    assert direct["p"] == [31]
    assert direct["Q"] == direct["b"][8] // 31

    # Actual reflected rows from the committed Q8345 regression.
    reflected = verify_packet(
        m=321,
        g=4,
        branch="reflected",
        hs=[36, 64, 100],
        expected_primes=[179, 193, 211],
        # Q8345 locks the needed edge cross-unit; this lower-arity selector
        # regression does not pretend to be the full sixteen-node isolated packet.
        require_cross_unit=False,
    )
    assert reflected["p"] == [179, 193, 211]
    for h, p in zip(reflected["h"], reflected["p"]):
        assert 2 * p == 321 + 1 + h
    return direct, reflected


def packet_from_json(payload: dict[str, Any]) -> dict[str, Any]:
    m = int(payload["m"])
    g = int(payload["g"])
    branch = str(payload["branch"])
    expected = payload.get("p")
    expected_primes = None if expected is None else [int(x) for x in expected]
    require_cross_unit = bool(payload.get("require_cross_unit", True))

    if "boolean_gaps" in payload:
        return verify_boolean_packet(
            m=m,
            g=g,
            branch=branch,
            a0=int(payload["a0"]),
            gaps=[int(x) for x in payload["boolean_gaps"]],
            expected_primes=expected_primes,
            require_cross_unit=require_cross_unit,
        )
    if "h" not in payload:
        raise ValueError("JSON must contain either h or a0+boolean_gaps")
    return verify_packet(
        m=m,
        g=g,
        branch=branch,
        hs=[int(x) for x in payload["h"]],
        expected_primes=expected_primes,
        require_cross_unit=require_cross_unit,
    )


def compact_summary(data: dict[str, Any]) -> dict[str, Any]:
    out = {
        "m": data["m"],
        "g": data["g"],
        "branch": data["branch"],
        "node_count": len(data["h"]),
        "h": data["h"],
        "p": data["p"],
        "R_digits": len(str(abs(data["R"]))),
        "Q_digits": len(str(abs(data["Q"]))),
        "all_targets_divide_C": all(data["C"] % p == 0 for p in data["p"]),
    }
    if "F_square" in data:
        out.update(
            {
                "F_square_sign": (data["F_square"] > 0) - (data["F_square"] < 0),
                "F_square_digits": len(str(abs(data["F_square"]))),
                "boolean_gaps": data["boolean_gaps"],
            }
        )
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--packet-json", type=Path)
    args = parser.parse_args()

    verify_growth()
    verify_linear_saturation()
    direct, reflected = verify_actual_regressions()

    print("LOCKED_DIRECT", json.dumps(compact_summary(direct), sort_keys=True))
    print("LOCKED_REFLECTED", json.dumps(compact_summary(reflected), sort_keys=True))

    if args.packet_json is not None:
        payload = json.loads(args.packet_json.read_text(encoding="utf-8"))
        data = packet_from_json(payload)
        print("INPUT_PACKET", json.dumps(compact_summary(data), sort_keys=True))

    print("Q8373_FULLPACKET_SELECTOR=PASS")


if __name__ == "__main__":
    main()
