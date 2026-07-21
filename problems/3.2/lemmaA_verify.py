#!/usr/bin/env python3
"""Fail-closed verifier for every computational claim in lemmaA_result.tex."""

from __future__ import annotations

import sys
import traceback

from lemmaA_explore import (
    apery_binomial_mod,
    apery_residues,
    decision_gate,
    eta_product_coefficients,
    jacobi_mod,
    jacobi_mod_direct,
    primes_up_to,
    skeleton_mod,
    term_valuation,
    term_valuation_simplified,
    torus_coordinate_mod,
    unit_count,
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


def check_jacobi_formula() -> str:
    pairs = 0
    for p in primes_up_to(61):
        if p < 5:
            continue
        for a in range(p - 1):
            for b in range(p - 1):
                actual = jacobi_mod_direct(p, a, b)
                expected = jacobi_mod(p, a, b)
                if actual != expected:
                    raise AssertionError((p, a, b, actual, expected))
                pairs += 1
    return f"{pairs} character pairs, including all degenerate cases"


def check_skeleton() -> str:
    pairs = 0
    for p in primes_up_to(200):
        if p < 5:
            continue
        row = apery_residues(p)
        for j, value in enumerate(row):
            binomial = apery_binomial_mod(j, p)
            if value != binomial:
                raise AssertionError((p, j, "recurrence", value, "binomial", binomial))
            actual = skeleton_mod(p, j)
            if actual != value:
                raise AssertionError((p, j, actual, value))
            pairs += 1
    return (
        f"{pairs} pairs, recurrence=binomial=skeleton for all primes "
        "5<=p<=200 and 0<=j<=p-2"
    )


def check_direct_marked_coordinate() -> str:
    pairs = 0
    for p in primes_up_to(31):
        if p < 5:
            continue
        row = apery_residues(p)
        for j, value in enumerate(row):
            torus = torus_coordinate_mod(p, j)
            skeleton = skeleton_mod(p, j)
            if torus != value or skeleton != value:
                raise AssertionError((p, j, torus, skeleton, value))
            pairs += 1
    return f"{pairs} independent torus sums through p=31"


def check_stickelberger_formula() -> str:
    triples = 0
    for p in primes_up_to(200):
        if p < 5:
            continue
        for j in range(1, p - 1):
            zero_indices = []
            for k in range(p - 1):
                actual = term_valuation(p, j, k)
                expected = term_valuation_simplified(p, j, k)
                if actual != expected:
                    raise AssertionError((p, j, k, actual, expected))
                if actual == 0:
                    zero_indices.append(k)
                triples += 1
            expected_indices = list(range(min(j, p - 1 - j) + 1))
            if zero_indices != expected_indices:
                raise AssertionError((p, j, zero_indices, expected_indices))
            if len(zero_indices) != unit_count(p, j):
                raise AssertionError((p, j, len(zero_indices), unit_count(p, j)))
    return f"{triples} (p,j,k) triples through p=200"


def check_eta_anchor() -> str:
    coefficients = eta_product_coefficients(200)
    if coefficients[11] != -44:
        raise AssertionError(("a_11", coefficients[11], -44))
    primes = 0
    for p in primes_up_to(200):
        if p < 5:
            continue
        midpoint = (p - 1) // 2
        b_midpoint = apery_residues(p)[midpoint]
        if b_midpoint != coefficients[p] % p:
            raise AssertionError((p, b_midpoint, coefficients[p]))
        if skeleton_mod(p, midpoint) != coefficients[p] % p:
            raise AssertionError((p, skeleton_mod(p, midpoint), coefficients[p]))
        primes += 1
    return f"{primes} primes, eta(2z)^4 eta(4z)^4 through p=200"


def check_container_2000() -> str:
    implications = 0
    zeros = 0
    for p in primes_up_to(2000):
        if p < 5:
            continue
        row = apery_residues(p)
        for j in range(1, p - 1):
            u = unit_count(p, j)
            if u == 1:
                implications += 1
                if row[j] == 0:
                    raise AssertionError(f"U=1 counterexample p={p}, j={j}")
            if row[j] == 0:
                zeros += 1
                if u == 1:
                    raise AssertionError(f"container counterexample p={p}, j={j}")
    if implications != 0:
        raise AssertionError("the explicit formula should make U=1 impossible here")
    return f"no counterexample; all {zeros} zeros have U>=2 (U=1 occurs 0 times)"


def check_decision_data_2000() -> str:
    expected_rows = {
        200: (44, 41, 25, 40, [(11, 5)], 4, {0: 25, 1: 1, 2: 16, 4: 2}),
        500: (93, 95, 53, 94, [(11, 5)], 6, {0: 53, 1: 1, 2: 33, 4: 4, 6: 2}),
        1000: (166, 165, 96, 164, [(11, 5)], 6, {0: 96, 1: 1, 2: 58, 4: 9, 6: 2}),
        2000: (
            301,
            283,
            184,
            282,
            [(11, 5)],
            6,
            {0: 184, 1: 1, 2: 95, 4: 17, 6: 4},
        ),
    }
    keys = (
        "prime_count",
        "total_zeros",
        "zero_free",
        "collision_zeros",
        "nonordinary",
        "max_zeros",
        "distribution",
    )
    for bound, expected in expected_rows.items():
        report = decision_gate(bound)
        actual = tuple(report[key] for key in keys)
        if actual != expected:
            raise AssertionError((bound, dict(zip(keys, actual)), dict(zip(keys, expected))))
    return "all four rows through 2000; final row has sum Z=283 and classes 0/282/1"


def check_decision_data_5000() -> str:
    report = decision_gate(5000)
    expected_distribution = {0: 418, 1: 2, 2: 199, 4: 41, 6: 6, 8: 1}
    expected = {
        "prime_count": 667,
        "total_zeros": 608,
        "collision_zeros": 606,
        "nonordinary": [(11, 5), (3137, 1568)],
        "zero_free": 418,
        "max_zeros": 8,
        "distribution": expected_distribution,
    }
    for key, value in expected.items():
        if report[key] != value:
            raise AssertionError((key, report[key], value))
    coefficients = eta_product_coefficients(3137)
    if coefficients[3137] != 207042 or coefficients[3137] != 66 * 3137:
        raise AssertionError((coefficients[3137], 66 * 3137))
    return "667 primes; sum Z=608; classes 0/606/2; a_3137=66*3137"


def main() -> int:
    verifier = Verifier()
    verifier.claim("Jacobi closed formula", check_jacobi_formula)
    verifier.claim("A1 skeleton = b_j", check_skeleton)
    verifier.claim("marked coordinate = skeleton = b_j", check_direct_marked_coordinate)
    verifier.claim("A2 Stickelberger valuation and U(j)", check_stickelberger_formula)
    verifier.claim("gamma_p midpoint anchor", check_eta_anchor)
    verifier.claim("U(j)=1 implies j notin Z_p; container p<=2000", check_container_2000)
    verifier.claim("A3 decision table p<=2000", check_decision_data_2000)
    verifier.claim("A3 stretch table p<=5000", check_decision_data_5000)
    print(f"SUMMARY: {verifier.passed} PASS, {verifier.failed} FAIL")
    if verifier.failed:
        print("LEMMA A VERIFICATION FAILED")
        return 1
    print("ALL LEMMA A CLAIMS PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
