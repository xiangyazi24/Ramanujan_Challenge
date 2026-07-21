#!/usr/bin/env python3
"""Fail-closed consistency checks accompanying lemmaBprime_result.tex."""

from __future__ import annotations

import traceback
from fractions import Fraction

from lemmaBprime_explore import (
    EXACT_TRACE_TABLE,
    apery_binomial_mod,
    apery_residues,
    delta_critical_values,
    direct_tlambda_polytope_volume,
    exact_extension_trace_vector,
    exact_probe,
    extension_exponent,
    extension_trace_mod_p,
    factor_coefficients,
    frobenius_slopes,
    hodge_data,
    newton_coefficients_mod_p,
    predicted_power_traces,
    primes_up_to,
    reduce_cyclotomic,
    slope_zero_count,
    small_range_report,
    split_prime_candidate_test_p7,
    stretch_report,
    tlambda_degenerate_witness,
)


class Verifier:
    def __init__(self) -> None:
        self.passed = 0
        self.failed = 0

    def claim(self, name: str, check) -> None:
        try:
            detail = check()
        except Exception:
            self.failed += 1
            print(f"FAIL {name}")
            traceback.print_exc()
        else:
            self.passed += 1
            suffix = f": {detail}" if detail else ""
            print(f"PASS {name}{suffix}")


def check_coordinate_change() -> str:
    samples = ((2, 3, 4), (3, 5, 2), (4, 2, 3), (5, 4, 6))
    checked = 0
    for u0, v0, w0 in samples:
        u, v, w = map(Fraction, (u0, v0, w0))
        if 1 in (u, v, w) or u * v * w == 1:
            continue
        x = -u * v * w
        y = 1 / (u - 1)
        z = 1 / (v - 1)
        original = (
            (1 + x)
            * (1 + y)
            * (1 + z)
            * ((1 + y) * (1 + z) + x * y * z)
            / (x * y * z)
        )
        arrangement = u * v * (w - 1) * (1 - u * v * w) / (
            w * (u - 1) * (v - 1)
        )
        if original != arrangement:
            raise AssertionError((u, v, w, original, arrangement))
        checked += 1
    return f"{checked} exact rational samples of Lambda(u,v,w)"


def check_arrangement_point_count() -> str:
    checked = 0
    for p in primes_up_to(31):
        if p < 5:
            continue
        actual = 0
        for u in range(1, p):
            for v in range(1, p):
                for w in range(1, p):
                    if u != 1 and v != 1 and w != 1 and u * v * w % p != 1:
                        actual += 1
        expected = p**3 - 7 * p**2 + 17 * p - 14
        if actual != expected:
            raise AssertionError((p, actual, expected))

        # Exact cyclotomic check of the exceptional-plane character sum.
        h = p - 1
        factors = []
        work = h
        divisor = 2
        while divisor * divisor <= work:
            if work % divisor == 0:
                factors.append(divisor)
                while work % divisor == 0:
                    work //= divisor
            divisor += 1
        if work > 1:
            factors.append(work)
        generator = next(
            candidate
            for candidate in range(2, p)
            if all(pow(candidate, h // factor, p) != 1 for factor in factors)
        )
        logs = [-1] * p
        value = 1
        for exponent in range(h):
            logs[value] = exponent
            value = value * generator % p
        for j in range(1, h):
            coefficients = [0] * h
            minus_exponent = j * logs[p - 1] % h
            for a in range(1, p):
                for c in range(1, p):
                    middle = (a + 1 + c) % p
                    if middle:
                        exponent = j * (
                            logs[c] + logs[middle] - logs[a]
                        ) + minus_exponent
                        coefficients[exponent % h] += 1
            reduced = reduce_cyclotomic(coefficients)
            if reduced != [p] + [0] * (len(reduced) - 1):
                raise AssertionError((p, j, "exceptional sum", reduced))
        checked += 1
    return f"{checked} primes; U-count polynomial and exceptional Tate sum q"


def check_as_obstruction_and_hodge() -> str:
    if direct_tlambda_polytope_volume() != 38:
        raise AssertionError("wrong direct polytope volume")
    if tlambda_degenerate_witness() != (-1, -1, 1, 1):
        raise AssertionError("wrong degeneracy witness")
    pairs = 0
    for p in primes_up_to(300):
        if p < 5:
            continue
        for j in range(1, p - 1):
            data = hodge_data(p, j)
            if data.residues != tuple(sorted((data.a, data.a, 1 - data.a, 1 - data.a))):
                raise AssertionError((p, j, data))
            if data.rho != (3, 2, 1, 0):
                raise AssertionError((p, j, data.rho))
            if data.generic_slopes != (0, 1, 2, 3):
                raise AssertionError((p, j, data.generic_slopes))
            if data.conifold_slopes != (0, 1, 3):
                raise AssertionError((p, j, data.conifold_slopes))
            if data.slope_zero_multiplicity != 1:
                raise AssertionError((p, j, data.slope_zero_multiplicity))
            pairs += 1
    return f"symbolic fractional-part/LMHS formula gives h0(j)=1 for {pairs} twists"


def check_extension_trace_reduction() -> str:
    triples = 0
    for p in primes_up_to(300):
        if p < 5:
            continue
        values = apery_residues(p)
        for j in range(1, p - 1):
            if values[j] != apery_binomial_mod(j, p):
                raise AssertionError((p, j, values[j], apery_binomial_mod(j, p)))
            for r in (1, 2, 3):
                n = extension_exponent(p, j, r)
                digits = []
                work = n
                for _ in range(r):
                    digits.append(work % p)
                    work //= p
                if digits != [j] * r or work:
                    raise AssertionError((p, j, r, n, digits, work))
                if extension_trace_mod_p(p, j, r) != pow(values[j], r, p):
                    raise AssertionError((p, j, r))
                triples += 1
            e1, e2, e3 = newton_coefficients_mod_p(p, j)
            if (e1, e2, e3) != (values[j], 0, 0):
                raise AssertionError((p, j, (e1, e2, e3), values[j]))
    return f"{triples} Lucas/constant-term congruences; P(T) mod p=1-b_j*T"


def check_signed_unit_sum_and_inclusion() -> str:
    pairs = 0
    zeros = 0
    for p in primes_up_to(300):
        if p < 5:
            continue
        values = apery_residues(p)
        for j in range(1, p - 1):
            b = values[j]
            count = slope_zero_count(p, j)
            effective_unit_sum = b if count else 0
            signed_rgamma_sum = -effective_unit_sum % p
            if signed_rgamma_sum != -b % p:
                raise AssertionError((p, j, b, signed_rgamma_sum))
            if b == 0:
                zeros += 1
                if count == 1:
                    raise AssertionError((p, j, "inclusion failed"))
            elif count != 1:
                raise AssertionError((p, j, "nonzero trace without unit root"))
            pairs += 1
    return (
        f"factor/sign consistency on {pairs} pairs; all {zeros} zeros have m0!=1"
    )


def check_exact_extension_a1() -> str:
    probes = 0
    for p, j, rmax in (
        (5, 1, 3),
        (5, 2, 3),
        (7, 1, 3),
        (7, 2, 3),
        (7, 3, 3),
        (11, 5, 3),
        (13, 1, 1),
    ):
        actual = exact_probe(p, j, rmax)
        expected = list(EXACT_TRACE_TABLE[p, j][:rmax])
        if actual != expected:
            raise AssertionError((p, j, actual, expected))
        probes += 1

    # A non-rational base-field trace, independently counted in Q(zeta_10).
    cyclotomic_s1 = reduce_cyclotomic(exact_extension_trace_vector(11, 1, 1))
    if cyclotomic_s1 != [-27, 0, -8, 8]:
        raise AssertionError(("p=11,j=1 cyclotomic S1", cyclotomic_s1))
    cyclotomic_e2 = [11 * coefficient for coefficient in cyclotomic_s1]
    cyclotomic_e2[0] += 11**3 - 11**2
    if cyclotomic_e2 != [913, 0, -88, 88] or 11**4 != 14641:
        raise AssertionError(("p=11,j=1 factor", cyclotomic_e2))
    # This also catches the forbidden shortcut that lifts only the p-1 base
    # Mellin characters: at (5,1,r=2) the full answer is -125, not -153.
    if EXACT_TRACE_TABLE[5, 1][1] != -125 or EXACT_TRACE_TABLE[5, 1][1] == -153:
        raise AssertionError("Hasse-Davenport shortcut regression")
    return "6 pairs r<=3; p13 r=1; cyclotomic p11,j1 S1"


def check_cubic_factor() -> str:
    rows = 0
    for (p, j), traces in EXACT_TRACE_TABLE.items():
        if list(traces) != predicted_power_traces(p, traces[0], len(traces)):
            raise AssertionError((p, j, traces))
        if len(traces) >= 3:
            s1, s2, s3 = traces[:3]
            numerator2 = s1 * s1 - s2
            numerator3 = s1**3 - 3 * s1 * s2 + 2 * s3
            if numerator2 % 2 or numerator3 % 6:
                raise AssertionError((p, j, "Newton divisibility"))
            actual = (s1, numerator2 // 2, numerator3 // 6)
            expected = factor_coefficients(p, s1)
            if actual != expected:
                raise AssertionError((p, j, actual, expected))
        rows += 1
    independent_cubics = sum(
        (p, j) != (13, 1) and len(traces) >= 3
        for (p, j), traces in EXACT_TRACE_TABLE.items()
    )
    fourth_trace_rows = sum(len(traces) >= 4 for traces in EXACT_TRACE_TABLE.values())
    if (rows, independent_cubics, fourth_trace_rows) != (7, 6, 5):
        raise AssertionError((rows, independent_cubics, fourth_trace_rows))
    return "6 probed cubics; 5 r4 recurrences; p13 predicted"


def check_small_newton_table() -> str:
    report = small_range_report(59)
    expected_exceptional = [
        (5, 1),
        (5, 3),
        (11, 5),
        (17, 3),
        (17, 13),
        (19, 8),
        (19, 10),
        (31, 8),
        (31, 22),
        (37, 17),
        (37, 19),
        (41, 10),
        (41, 30),
        (59, 9),
        (59, 49),
    ]
    expected = (405, 390, 15, expected_exceptional)
    actual = (
        report["pairs"],
        report["ordinary"],
        report["nonordinary"],
        report["exceptional"],
    )
    if actual != expected:
        raise AssertionError((actual, expected))
    for _, _, _, count, slopes in report["records"]:
        if count >= 2:
            raise AssertionError("found two slope-zero roots")
        if count == 1 and slopes != (Fraction(0), Fraction(1), Fraction(3)):
            raise AssertionError((count, slopes))
    return "405 factor classifications; base-A1 v(a)=1 at 15/15 zeros"


def check_stretch_table() -> str:
    report = stretch_report(300, 3_202_026)
    if (report["zeros"], report["random"]) != (63, 1140):
        raise AssertionError((report["zeros"], report["random"]))
    if any(m != 1 or slopes != (Fraction(1), Fraction(1), Fraction(2))
           for _, _, m, slopes in report["zero_records"]):
        raise AssertionError("a zero did not have valuation 1 and slopes (1,1,2)")
    random_zero_hits = sum(m > 0 for _, _, m, _ in report["random_records"])
    random_ordinary = sum(m == 0 for _, _, m, _ in report["random_records"])
    if (random_ordinary, random_zero_hits) != (1127, 13):
        raise AssertionError((random_ordinary, random_zero_hits))
    if any(
        slopes
        != (
            (Fraction(0), Fraction(1), Fraction(3))
            if m == 0
            else (Fraction(1), Fraction(1), Fraction(2))
        )
        for _, _, m, slopes in report["random_records"]
    ):
        raise AssertionError("a random point has the wrong Newton polygon")
    return "base-A1 v(a)=1 at 63/63 zeros; sample 1140, zero hits 13"


def check_delta_and_branches() -> str:
    if delta_critical_values() != ((17, 12), (17, -12)):
        raise AssertionError("wrong Delta critical values")
    result = split_prime_candidate_test_p7()
    if result["unit_residues"] != {1: 5, 2: 3, 3: 3}:
        raise AssertionError(result)
    if result["phi3_remainders"] != {1: (8, 342), 2: (24, 342)}:
        raise AssertionError(result)
    if result["Q3_at_1"] != 320:
        raise AssertionError(result)
    if not result["all_local_eigenvalues_nonsquare"]:
        raise AssertionError(result)
    if result["unit_root_fields"] != (-1291, -83, -199):
        raise AssertionError(result)
    if result["square_obstructions"] != {
        "unit_residue_nonsquare": True,
        "other_valuations_odd": True,
    }:
        raise AssertionError(result)
    return "critical values; p7 exact Teich/Jacobi/Q7-square rejection"


def check_binary_gate() -> str:
    # G2 requires h0>=2.  The derived formula is h0=1 for every nontrivial j.
    for p in (5, 7, 11, 59, 293):
        for j in range(1, p - 1):
            if hodge_data(p, j).slope_zero_multiplicity != 1:
                raise AssertionError((p, j))
    return "G1 (h0=1; the only fixed branch is Tate slope 1)"


def main() -> int:
    verifier = Verifier()
    verifier.claim("marked-coordinate change", check_coordinate_change)
    verifier.claim("toric-arrangement point count", check_arrangement_point_count)
    verifier.claim("symbolic AS/Hodge formula", check_as_obstruction_and_hodge)
    verifier.claim("extension trace congruence", check_extension_trace_reduction)
    verifier.claim("signed slope-zero/inclusion consistency", check_signed_unit_sum_and_inclusion)
    verifier.claim("full A1 probes", check_exact_extension_a1)
    verifier.claim("cubic factor", check_cubic_factor)
    verifier.claim("p<=59 classification", check_small_newton_table)
    verifier.claim("p<=300 sample", check_stretch_table)
    verifier.claim("candidate matches", check_delta_and_branches)
    verifier.claim("binary gate", check_binary_gate)
    print(f"SUMMARY: {verifier.passed} PASS, {verifier.failed} FAIL")
    return int(verifier.failed != 0)


if __name__ == "__main__":
    raise SystemExit(main())
