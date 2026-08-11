#!/usr/bin/env python3
"""Q7309: exact p^2 audit of Apéry high-load rows.

This is a standalone, standard-library-only verifier.  It deliberately keeps
three computations separate:

  1. Apéry values below p from the divided three-term recurrence modulo p^2;
  2. the cleared recurrence B_m=(m!)^3 b_m modulo p^2;
  3. b_{qp+r} modulo p^2 from a two-base-p-digit binomial calculation.

The third computation is independent of the Apéry recurrence.  Terms having a
p-adic carry are multiples of p^2 after squaring and may be omitted.  The
remaining no-carry terms are evaluated with the standard two-digit binomial
expansion modulo p^2 and then factorized to O(p) work per hit.

No proof claim is made.  The relation screen is deliberately tiny and fixed by
PREREG_SPEC below before any scan output is examined.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import math
import os
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple


DEFAULT_XS = (128, 256, 512, 1024)
TRAIN_MAX_X = 512
HOLDOUT_X = 1024
COEFF_ALPHABET = (-2, -1, 1, 2)
MAX_SUPPORT = 3
VARIABLES = ("q", "rho", "tau", "eta", "jet_u")

PREREG_SPEC = """Q7309 relation screen v1
scales: X in {128,256,512,1024}
training: all high-load hits with X<=512
holdout: all high-load hits with X=1024
variables: q,rho,tau,eta,jet_u where jet_u=(b_n/p mod p)/u mod p
linear dictionary: 1 plus the five variables
quadratic dictionary: all monomials of total degree at most two
candidate coefficients: primitive, first nonzero positive, each nonzero in {-2,-1,1,2}
support bound: at most three monomials
false positive: exact zero on every training hit and nonzero on at least one holdout hit
matched control: replace q by (q+1) mod X and recompute the Gessel jet, keeping the local zero labels fixed
"""
PREREG_SHA256 = hashlib.sha256(PREREG_SPEC.encode("utf-8")).hexdigest()


def P(n: int) -> int:
    return 34 * n**3 + 51 * n**2 + 27 * n + 5


def Pprime(n: int) -> int:
    return 102 * n**2 + 102 * n + 27


def inv_mod(a: int, modulus: int) -> int:
    a %= modulus
    if math.gcd(a, modulus) != 1:
        raise AssertionError(f"nonunit inversion: {a} mod {modulus}")
    return pow(a, -1, modulus)


def sieve_primes(limit: int) -> List[int]:
    mark = bytearray(b"\x01") * (limit + 1)
    if limit >= 0:
        mark[0] = 0
    if limit >= 1:
        mark[1] = 0
    for a in range(2, math.isqrt(limit) + 1):
        if mark[a]:
            mark[a * a : limit + 1 : a] = b"\x00" * (((limit - a * a) // a) + 1)
    return [i for i in range(2, limit + 1) if mark[i]]


@dataclass
class PrimeData:
    p: int
    p2: int
    b2: List[int]
    deriv: List[int]
    cleared: List[int]
    fact2: List[int]
    invfact2: List[int]
    harmonic: List[int]
    zeros: List[int]

    def choose_small(self, n: int, k: int) -> int:
        """C(n,k) modulo p^2, for 0<=k<=n<p."""
        if k < 0 or k > n or n >= self.p:
            return 0
        return (
            self.fact2[n]
            * self.invfact2[k]
            % self.p2
            * self.invfact2[n - k]
            % self.p2
        )


def build_prime_data(p: int, exhaustive_checks: bool = True) -> PrimeData:
    if p < 5:
        raise ValueError("the independent two-digit formula is used only for p>=5")
    p2 = p * p

    # Divided Apéry recurrence modulo p^2, safe for all indices < p.
    b2 = [0] * p
    b2[0] = 1
    if p > 1:
        b2[1] = 5
    for n in range(1, p - 1):
        den = pow(n + 1, 3, p2)
        num = (P(n) * b2[n] - n**3 * b2[n - 1]) % p2
        b2[n + 1] = num * inv_mod(den, p2) % p2

    # Factorials and inverse factorials below p, modulo p^2.
    fact2 = [1] * p
    for n in range(1, p):
        fact2[n] = fact2[n - 1] * n % p2
    invfact2 = [1] * p
    invfact2[p - 1] = inv_mod(fact2[p - 1], p2)
    for n in range(p - 1, 0, -1):
        invfact2[n - 1] = invfact2[n] * n % p2

    # Cleared recurrence Y_m=(m!)^3 b_m.
    cleared = [0] * p
    cleared[0] = 1
    if p > 1:
        cleared[1] = 5
    for n in range(1, p - 1):
        cleared[n + 1] = (P(n) * cleared[n] - n**6 * cleared[n - 1]) % p2

    # Formal derivative, obtained by differentiating the recurrence.
    bp = [x % p for x in b2]
    deriv = [0] * p
    for n in range(0, p - 1):
        bm1 = bp[n - 1] if n > 0 else 0
        dm1 = deriv[n - 1] if n > 0 else 0
        rhs = (
            P(n) * deriv[n]
            - n**3 * dm1
            + Pprime(n) * bp[n]
            - 3 * (n + 1) ** 2 * bp[n + 1]
            - 3 * n**2 * bm1
        ) % p
        deriv[n + 1] = rhs * inv_mod(pow(n + 1, 3, p), p) % p

    # H_j modulo p, 0<=j<p.
    inverses = [0] * p
    if p > 1:
        inverses[1] = 1
    for n in range(2, p):
        inverses[n] = (-(p // n) * inverses[p % n]) % p
    harmonic = [0] * p
    for n in range(1, p):
        harmonic[n] = (harmonic[n - 1] + inverses[n]) % p

    zeros = [r for r, value in enumerate(bp) if value == 0]

    if exhaustive_checks:
        for n in range(p):
            expected = pow(fact2[n], 3, p2) * b2[n] % p2
            assert cleared[n] == expected, (p, n, cleared[n], expected)
            assert bp[n] == bp[p - 1 - n], ("reflection", p, n)
        for n in range(0, p - 1):
            bm1 = bp[n - 1] if n > 0 else 0
            dm1 = deriv[n - 1] if n > 0 else 0
            lhs = pow(n + 1, 3, p) * deriv[n + 1] % p
            rhs = (
                P(n) * deriv[n]
                - n**3 * dm1
                + Pprime(n) * bp[n]
                - 3 * (n + 1) ** 2 * bp[n + 1]
                - 3 * n**2 * bm1
            ) % p
            assert lhs == rhs, ("formal derivative", p, n)

    return PrimeData(
        p=p,
        p2=p2,
        b2=b2,
        deriv=deriv,
        cleared=cleared,
        fact2=fact2,
        invfact2=invfact2,
        harmonic=harmonic,
        zeros=zeros,
    )


def apery_two_digit_mod_p2(pd: PrimeData, q: int, r: int) -> int:
    """Compute b_{q p+r} mod p^2 independently from the binomial sum.

    Here 0<=q,r<p.  Write k=a p+s.  A summand
        C(qp+r,k)^2 C(qp+r+k,k)^2
    is nonzero modulo p^2 only when there is no carry in either binomial:
        0<=a<=min(q,p-1-q), 0<=s<=min(r,p-1-r).

    For those terms, the two-digit binomial expansion modulo p^2 is used.
    Its correction separates into an a-part and an s-part, reducing the
    double sum to four one-dimensional sums.
    """
    p, p2, H = pd.p, pd.p2, pd.harmonic
    if not (0 <= q < p and 0 <= r < p):
        raise ValueError((p, q, r))
    sigma = min(q, p - 1 - q)
    rho = min(r, p - 1 - r)

    sum_a_mod_p2 = 0
    sum_a_weighted_mod_p = 0
    for a in range(sigma + 1):
        z = pd.choose_small(q, a) * pd.choose_small(q + a, a) % p2
        term = z * z % p2
        sum_a_mod_p2 = (sum_a_mod_p2 + term) % p2
        sum_a_weighted_mod_p = (sum_a_weighted_mod_p + a * (term % p)) % p

    sum_s_mod_p2 = 0
    sum_s_F_mod_p = 0
    sum_s_G_mod_p = 0
    for s in range(rho + 1):
        z = pd.choose_small(r, s) * pd.choose_small(r + s, s) % p2
        term = z * z % p2
        sum_s_mod_p2 = (sum_s_mod_p2 + term) % p2
        F = (H[r + s] - H[r - s]) % p
        G = (H[r - s] + H[r + s] - 2 * H[s]) % p
        term_p = term % p
        sum_s_F_mod_p = (sum_s_F_mod_p + term_p * F) % p
        sum_s_G_mod_p = (sum_s_G_mod_p + term_p * G) % p

    base = sum_a_mod_p2 * sum_s_mod_p2 % p2
    correction_inner = (
        q * (sum_a_mod_p2 % p) * sum_s_F_mod_p
        + sum_a_weighted_mod_p * sum_s_G_mod_p
    ) % p
    return (base + 2 * p * correction_inner) % p2


def apery_exact(n: int) -> int:
    return sum(math.comb(n, k) ** 2 * math.comb(n + k, k) ** 2 for k in range(n + 1))


def independent_formula_selftest() -> Dict[str, int]:
    tested = 0
    for p in (5, 7, 11):
        pd = build_prime_data(p)
        for q in range(p):
            for r in range(p):
                n = q * p + r
                got = apery_two_digit_mod_p2(pd, q, r)
                want = apery_exact(n) % (p * p)
                assert got == want, ("two-digit self-test", p, q, r, got, want)
                tested += 1
    return {"small_primes": 3, "two_digit_values": tested}


def monomial_names(total_degree: int) -> List[str]:
    names = ["1"] + list(VARIABLES)
    if total_degree >= 2:
        for i, a in enumerate(VARIABLES):
            for b in VARIABLES[i:]:
                names.append(f"{a}*{b}")
    return names


def monomial_values(record: Mapping[str, int], total_degree: int) -> List[int]:
    p = int(record["p"])
    vals = [1] + [int(record[name]) % p for name in VARIABLES]
    if total_degree >= 2:
        for i, a in enumerate(VARIABLES):
            av = int(record[a]) % p
            for b in VARIABLES[i:]:
                vals.append(av * (int(record[b]) % p) % p)
    return vals


def control_monomial_values(record: Mapping[str, int], total_degree: int) -> List[int]:
    p = int(record["p"])
    control = dict(record)
    control["q"] = int(record["control_q"])
    control["jet_u"] = int(record["control_jet_u"])
    return monomial_values(control, total_degree)


def sparse_candidates(dimension: int) -> Iterable[Tuple[Tuple[int, ...], Tuple[int, ...]]]:
    for support_size in range(1, min(MAX_SUPPORT, dimension) + 1):
        for support in itertools.combinations(range(dimension), support_size):
            for coeffs in itertools.product(COEFF_ALPHABET, repeat=support_size):
                if coeffs[0] < 0:
                    continue
                if math.gcd(*[abs(c) for c in coeffs]) != 1:
                    continue
                yield support, coeffs


def polynomial_string(names: Sequence[str], support: Sequence[int], coeffs: Sequence[int]) -> str:
    pieces: List[str] = []
    for idx, coefficient in zip(support, coeffs):
        name = names[idx]
        abs_part = name if abs(coefficient) == 1 else f"{abs(coefficient)}*{name}"
        if not pieces:
            pieces.append(abs_part if coefficient > 0 else f"-{abs_part}")
        else:
            pieces.append((" + " if coefficient > 0 else " - ") + abs_part)
    return "".join(pieces) if pieces else "0"


def candidate_failures(
    feature_rows: Sequence[Tuple[int, Sequence[int], Mapping[str, int]]],
    support: Sequence[int],
    coeffs: Sequence[int],
) -> Tuple[int, str]:
    failures = 0
    first = ""
    for p, values, record in feature_rows:
        residue = sum(c * values[j] for j, c in zip(support, coeffs)) % p
        if residue:
            failures += 1
            if not first:
                first = f"X={record['X']};n={record['n']};p={p};residue={residue}"
    return failures, first


def rank_mod(matrix: Sequence[Sequence[int]], p: int) -> int:
    if not matrix:
        return 0
    a = [[x % p for x in row] for row in matrix]
    rows, cols = len(a), len(a[0])
    rank = 0
    for col in range(cols):
        pivot = next((i for i in range(rank, rows) if a[i][col] % p), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        inv = inv_mod(a[rank][col], p)
        a[rank] = [(x * inv) % p for x in a[rank]]
        for i in range(rows):
            if i == rank or a[i][col] == 0:
                continue
            factor = a[i][col]
            a[i] = [(x - factor * y) % p for x, y in zip(a[i], a[rank])]
        rank += 1
        if rank == rows:
            break
    return rank


def write_csv(path: Path, rows: Sequence[Mapping[str, object]], fieldnames: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def scan_scale(X: int, all_primes: Sequence[int]) -> Tuple[dict, List[dict], List[dict], List[dict]]:
    primes = [p for p in all_primes if X < p <= 2 * X]
    pdata: Dict[int, PrimeData] = {p: build_prime_data(p) for p in primes}
    N = X * X

    loads: Dict[int, List[int]] = defaultdict(list)
    total_events = 0
    prime_rows: List[dict] = []
    for p in primes:
        pd = pdata[p]
        prime_rows.append(
            {
                "X": X,
                "p": p,
                "zero_count": len(pd.zeros),
                "zeros": ";".join(map(str, pd.zeros)),
            }
        )
        for r in pd.zeros:
            for n in range(r, N, p):
                loads[n].append(p)
                total_events += 1

    high_rows = [(n, sorted(ps)) for n, ps in loads.items() if len(ps) >= 3]
    high_rows.sort()
    max_load = max((len(ps) for ps in loads.values()), default=0)
    canonical_triples = sum(math.comb(len(ps), 3) for _, ps in high_rows)
    known_triples = {256: 10, 512: 25, 1024: 81}
    if X in known_triples:
        assert canonical_triples == known_triples[X], (
            "census mismatch",
            X,
            canonical_triples,
            known_triples[X],
        )

    row_rows: List[dict] = []
    hit_rows: List[dict] = []
    for n, ps in high_rows:
        certificates: List[str] = []
        for hit_index, p in enumerate(ps):
            pd = pdata[p]
            p2 = p * p
            q, r = divmod(n, p)
            assert q < X < p
            assert r in pd.zeros
            rho = min(r, p - 1 - r)
            ref = p - 1 - rho
            side = "M" if rho == ref else ("L" if r == rho else "R")

            brho2 = pd.b2[rho]
            bref2 = pd.b2[ref]
            br2 = pd.b2[r]
            assert brho2 % p == 0
            assert bref2 % p == 0
            assert br2 % p == 0

            alpha = (brho2 // p) % p
            reflection_delta2 = (bref2 - brho2) % p2
            assert reflection_delta2 % p == 0
            epsilon = (reflection_delta2 // p) % p

            u = pow(rho + 1, 3, p) * (pd.b2[rho + 1] % p) % p
            assert u != 0, ("zero local slope", X, n, p, rho)
            inv_u = inv_mod(u, p)
            tau = alpha * inv_u % p
            eta = epsilon * inv_u % p

            # Independent alpha extraction from the cleared recurrence.
            y = pd.cleared[rho]
            assert y % p == 0
            beta_cleared = (y // p) % p
            fact_rho_cube = pow(pd.fact2[rho] % p, 3, p)
            alpha_cleared = beta_cleared * inv_mod(fact_rho_cube, p) % p
            assert alpha_cleared == alpha

            beta_r = (br2 // p) % p
            beta_expected = (alpha + (epsilon if side == "R" else 0)) % p
            assert beta_r == beta_expected, (
                "oriented reflection lift",
                X,
                n,
                p,
                beta_r,
                beta_expected,
            )

            deriv_r = pd.deriv[r]
            deriv_rho = pd.deriv[rho]
            deriv_ref = pd.deriv[ref]
            delta_r = deriv_r * inv_u % p
            delta_rho = deriv_rho * inv_u % p
            delta_ref = deriv_ref * inv_u % p

            bq2 = pd.b2[q]
            bq = bq2 % p
            full_prediction = (
                br2 * bq2 + p * (q % p) * deriv_r * bq
            ) % p2
            assert full_prediction % p == 0
            jet_prediction = bq * ((beta_r + q * deriv_r) % p) % p
            assert (full_prediction // p) % p == jet_prediction

            direct = apery_two_digit_mod_p2(pd, q, r)
            assert direct % p == 0, ("direct hit not divisible", X, n, p, direct)
            direct_jet = (direct // p) % p
            assert direct == full_prediction, (
                "Gessel p^2 mismatch",
                X,
                n,
                p,
                direct,
                full_prediction,
            )
            assert direct_jet == jet_prediction
            jet_u = direct_jet * inv_u % p

            # q=0 specialization independently recovers the local p^2 value.
            assert apery_two_digit_mod_p2(pd, 0, r) == br2

            control_q = (q + 1) % X
            control_bq = pd.b2[control_q] % p
            control_jet = control_bq * ((beta_r + control_q * deriv_r) % p) % p
            control_jet_u = control_jet * inv_u % p

            hit = {
                "X": X,
                "n": n,
                "K": len(ps),
                "hit_index": hit_index,
                "p": p,
                "q": q,
                "r": r,
                "rho": rho,
                "reflection_endpoint": ref,
                "side": side,
                "b_r_mod_p2": br2,
                "b_rho_mod_p2": brho2,
                "b_ref_mod_p2": bref2,
                "alpha": alpha,
                "alpha_cleared": alpha_cleared,
                "reflection_defect": epsilon,
                "u": u,
                "tau": tau,
                "eta": eta,
                "p2_divides_b_rho": int(alpha == 0),
                "beta_r": beta_r,
                "beta_expected": beta_expected,
                "formal_derivative_r": deriv_r,
                "formal_derivative_rho": deriv_rho,
                "formal_derivative_ref": deriv_ref,
                "delta_r": delta_r,
                "delta_rho": delta_rho,
                "delta_ref": delta_ref,
                "b_q_mod_p": bq,
                "b_q_mod_p2": bq2,
                "gessel_full_prediction_mod_p2": full_prediction,
                "direct_b_n_mod_p2": direct,
                "gessel_jet_prediction": jet_prediction,
                "direct_b_n_over_p_mod_p": direct_jet,
                "jet_u": jet_u,
                "control_q": control_q,
                "control_jet_u": control_jet_u,
                "all_checks": "PASS",
            }
            hit_rows.append(hit)
            certificates.append(
                f"p={p},r={r},rho={rho},q={q},a={alpha},e={epsilon},"
                f"t={tau},h={eta},j={direct_jet}"
            )

        row_rows.append(
            {
                "X": X,
                "n": n,
                "K": len(ps),
                "primes": ";".join(map(str, ps)),
                "certificate": " | ".join(certificates),
            }
        )

    summary = {
        "X": X,
        "N": N,
        "prime_count": len(primes),
        "active_prime_count": sum(bool(pdata[p].zeros) for p in primes),
        "sum_zero_counts": sum(len(pdata[p].zeros) for p in primes),
        "scattered_hit_events": total_events,
        "high_load_row_count": len(high_rows),
        "high_load_hit_count": sum(len(ps) for _, ps in high_rows),
        "K_histogram": dict(sorted(Counter(len(ps) for _, ps in high_rows).items())),
        "max_K": max_load,
        "canonical_triples": canonical_triples,
        "ordered_S3": 6 * canonical_triples,
        "all_checks": "PASS",
    }
    return summary, row_rows, hit_rows, prime_rows


def relation_and_rank_screen(hit_rows: Sequence[dict]) -> Tuple[List[dict], List[dict], dict]:
    rank_rows: List[dict] = []
    relation_rows: List[dict] = []
    relation_summary: dict = {
        "prereg_sha256": PREREG_SHA256,
        "prereg_spec": PREREG_SPEC,
        "dictionaries": {},
    }

    for degree, dictionary_name in ((1, "linear"), (2, "quadratic")):
        names = monomial_names(degree)
        feature_by_record = [
            (int(record["p"]), monomial_values(record, degree), record)
            for record in hit_rows
        ]
        control_by_record = [
            (int(record["p"]), control_monomial_values(record, degree), record)
            for record in hit_rows
        ]

        # Legitimate ranks are per characteristic; there is no single field in
        # which rows belonging to different p can be pooled.
        for split_name, predicate in (
            ("training", lambda rec: int(rec["X"]) <= TRAIN_MAX_X),
            ("holdout", lambda rec: int(rec["X"]) == HOLDOUT_X),
            ("all", lambda rec: True),
        ):
            grouped: Dict[int, List[List[int]]] = defaultdict(list)
            for p, values, record in feature_by_record:
                if predicate(record):
                    grouped[p].append(list(values))
            for p in sorted(grouped):
                matrix = grouped[p]
                rank = rank_mod(matrix, p)
                rank_rows.append(
                    {
                        "dictionary": dictionary_name,
                        "degree": degree,
                        "split": split_name,
                        "p": p,
                        "samples": len(matrix),
                        "columns": len(names),
                        "rank": rank,
                        "nullity": len(names) - rank,
                    }
                )

        subsets: Dict[str, List[Tuple[int, Sequence[int], Mapping[str, int]]]] = {
            "X128": [x for x in feature_by_record if int(x[2]["X"]) == 128],
            "X256": [x for x in feature_by_record if int(x[2]["X"]) == 256],
            "X512": [x for x in feature_by_record if int(x[2]["X"]) == 512],
            "training": [x for x in feature_by_record if int(x[2]["X"]) <= TRAIN_MAX_X],
            "holdout": [x for x in feature_by_record if int(x[2]["X"]) == HOLDOUT_X],
            "all": list(feature_by_record),
            "controls_all": list(control_by_record),
        }
        survivor_counts = Counter()
        candidate_count = 0
        false_positive_count = 0
        all_survivors: List[str] = []
        holdout_false_positives: List[str] = []

        for support, coeffs in sparse_candidates(len(names)):
            candidate_count += 1
            polynomial = polynomial_string(names, support, coeffs)
            failures: Dict[str, int] = {}
            first_counterexamples: Dict[str, str] = {}
            for subset_name, rows in subsets.items():
                failure_count, first = candidate_failures(rows, support, coeffs)
                failures[subset_name] = failure_count
                first_counterexamples[subset_name] = first
                if failure_count == 0:
                    survivor_counts[subset_name] += 1

            is_false_positive = failures["training"] == 0 and failures["holdout"] > 0
            survives_all = failures["all"] == 0
            if is_false_positive:
                false_positive_count += 1
                holdout_false_positives.append(polynomial)
            if survives_all:
                all_survivors.append(polynomial)

            relation_rows.append(
                {
                    "dictionary": dictionary_name,
                    "degree": degree,
                    "polynomial": polynomial,
                    "support": len(support),
                    "coefficients": ";".join(
                        f"{names[j]}:{c}" for j, c in zip(support, coeffs)
                    ),
                    "fail_X128": failures["X128"],
                    "fail_X256": failures["X256"],
                    "fail_X512": failures["X512"],
                    "fail_training": failures["training"],
                    "fail_holdout": failures["holdout"],
                    "fail_all": failures["all"],
                    "fail_controls_all": failures["controls_all"],
                    "false_positive": int(is_false_positive),
                    "survives_all": int(survives_all),
                    "first_holdout_counterexample": first_counterexamples["holdout"],
                }
            )

        relation_summary["dictionaries"][dictionary_name] = {
            "monomials": names,
            "column_count": len(names),
            "candidate_count": candidate_count,
            "survivor_counts": dict(survivor_counts),
            "training_survivors": survivor_counts["training"],
            "holdout_false_positives": false_positive_count,
            "all_survivors": all_survivors,
            "false_positive_polynomials": holdout_false_positives,
        }

    rank_histograms: dict = {}
    for dictionary_name in ("linear", "quadratic"):
        rank_histograms[dictionary_name] = {}
        for split in ("training", "holdout", "all"):
            vals = [
                int(row["rank"])
                for row in rank_rows
                if row["dictionary"] == dictionary_name and row["split"] == split
            ]
            rank_histograms[dictionary_name][split] = {
                "prime_count": len(vals),
                "rank_histogram": dict(sorted(Counter(vals).items())),
                "min_rank": min(vals) if vals else 0,
                "max_rank": max(vals) if vals else 0,
                "sum_rank": sum(vals),
            }
    relation_summary["rank_histograms"] = rank_histograms
    return rank_rows, relation_rows, relation_summary


def markdown_report(
    selftest: Mapping[str, int],
    scale_summaries: Sequence[dict],
    rows: Sequence[dict],
    hits: Sequence[dict],
    relation_summary: Mapping[str, object],
) -> str:
    lines: List[str] = []
    lines.append("# Q7309 exact p-adic high-load scan")
    lines.append("")
    lines.append("Mechanical result only; no proof claim.")
    lines.append("")
    lines.append("## Frozen scope and result")
    lines.append("")
    lines.append(f"Preregistration SHA-256: `{PREREG_SHA256}`.")
    lines.append("")
    lines.append("```text")
    lines.extend(PREREG_SPEC.rstrip().splitlines())
    lines.append("```")
    lines.append("")
    lines.append(
        f"Independent two-digit formula self-test: PASS on {selftest['two_digit_values']} "
        f"values for p=5,7,11."
    )
    lines.append("")
    lines.append("| X | primes | active | sum |Z_p| | K>=3 rows | hits | K histogram | max K | C(K,3) sum |")
    lines.append("|---:|---:|---:|---:|---:|---:|---|---:|---:|")
    for s in scale_summaries:
        lines.append(
            f"| {s['X']} | {s['prime_count']} | {s['active_prime_count']} | "
            f"{s['sum_zero_counts']} | {s['high_load_row_count']} | "
            f"{s['high_load_hit_count']} | `{s['K_histogram']}` | {s['max_K']} | "
            f"{s['canonical_triples']} |"
        )
    lines.append("")
    lines.append(
        f"All {len(hits)} high-load hit records passed the recurrence, cleared-recurrence, "
        "reflection, integrality, nonzero-slope, oriented-lift, independent binomial, "
        "and full Gessel congruence assertions."
    )
    lines.append("")
    lines.append("## Gessel first jet checked")
    lines.append("")
    lines.append("For n=qp+r and p|b_r, the checked prediction is")
    lines.append("")
    lines.append("$$")
    lines.append("b_{qp+r}\\equiv b_r b_q+p q b'_r b_q\\pmod {p^2},")
    lines.append("\\qquad")
    lines.append("\\frac{b_{qp+r}}p\\equiv b_q\\left(\\frac{b_r}p+q b'_r\\right)\\pmod p.")
    lines.append("$$")
    lines.append("")
    lines.append(
        "The left side was computed independently from the defining binomial sum using "
        "the two-base-p-digit no-carry expansion, not from this congruence."
    )
    lines.append("")
    lines.append("## Relation screen")
    lines.append("")
    dictionaries = relation_summary["dictionaries"]  # type: ignore[index]
    lines.append("| dictionary | columns | candidates | training survivors | held-out false positives | all-data survivors |")
    lines.append("|---|---:|---:|---:|---:|---:|")
    for name in ("linear", "quadratic"):
        d = dictionaries[name]
        lines.append(
            f"| {name} | {d['column_count']} | {d['candidate_count']} | "
            f"{d['training_survivors']} | {d['holdout_false_positives']} | "
            f"{len(d['all_survivors'])} |"
        )
    lines.append("")
    lines.append(
        "Ranks are reported per characteristic in `ranks.csv`; pooling rows over different "
        "prime fields into one numerical rank would be meaningless."
    )
    lines.append("")
    for name in ("linear", "quadratic"):
        d = dictionaries[name]
        if d["false_positive_polynomials"]:
            lines.append(f"Held-out rejections for {name}:")
            for poly in d["false_positive_polynomials"]:
                lines.append(f"- `{poly}`")
        if d["all_survivors"]:
            lines.append(f"All-data survivors for {name} (not proof claims):")
            for poly in d["all_survivors"]:
                lines.append(f"- `{poly}`")
    lines.append("")
    lines.append("## Exact row certificates")
    lines.append("")
    for X in DEFAULT_XS:
        subset = [row for row in rows if int(row["X"]) == X]
        lines.append(f"### X={X}")
        lines.append("")
        if not subset:
            lines.append("No K>=3 rows.")
            lines.append("")
            continue
        lines.append("| n | K | primes |")
        lines.append("|---:|---:|---|")
        for row in subset:
            lines.append(f"| {row['n']} | {row['K']} | `{row['primes']}` |")
        lines.append("")
    lines.append("## Files")
    lines.append("")
    lines.append("- `rows.csv`: every exact K>=3 row.")
    lines.append("- `hits.csv`: every p-adic hit certificate and every checked first jet.")
    lines.append("- `primes.csv`: exact zero sets used by the scatter.")
    lines.append("- `ranks.csv`: per-characteristic ranks.")
    lines.append("- `relations.csv`: every preregistered sparse candidate and all failure counts.")
    lines.append("- `summary.json`: machine-readable summary.")
    lines.append("- `MANIFEST.sha256`: hashes of code and exact outputs.")
    lines.append("")
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--xs", default=",".join(map(str, DEFAULT_XS)))
    parser.add_argument("--out", default="problems/3.2/q7309_exact")
    args = parser.parse_args()
    xs = tuple(int(x) for x in args.xs.split(",") if x)
    if xs != DEFAULT_XS:
        raise SystemExit(f"fail closed: preregistered scales are exactly {DEFAULT_XS}, got {xs}")
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    selftest = independent_formula_selftest()
    all_primes = sieve_primes(2 * max(xs))

    scale_summaries: List[dict] = []
    row_rows: List[dict] = []
    hit_rows: List[dict] = []
    prime_rows: List[dict] = []
    for X in xs:
        summary, rows, hits, primes = scan_scale(X, all_primes)
        scale_summaries.append(summary)
        row_rows.extend(rows)
        hit_rows.extend(hits)
        prime_rows.extend(primes)
        print(json.dumps(summary, sort_keys=True), flush=True)

    rank_rows, relation_rows, relation_summary = relation_and_rank_screen(hit_rows)

    row_fields = ["X", "n", "K", "primes", "certificate"]
    hit_fields = [
        "X", "n", "K", "hit_index", "p", "q", "r", "rho",
        "reflection_endpoint", "side", "b_r_mod_p2", "b_rho_mod_p2",
        "b_ref_mod_p2", "alpha", "alpha_cleared", "reflection_defect",
        "u", "tau", "eta", "p2_divides_b_rho", "beta_r", "beta_expected",
        "formal_derivative_r", "formal_derivative_rho", "formal_derivative_ref",
        "delta_r", "delta_rho", "delta_ref", "b_q_mod_p", "b_q_mod_p2",
        "gessel_full_prediction_mod_p2", "direct_b_n_mod_p2",
        "gessel_jet_prediction", "direct_b_n_over_p_mod_p", "jet_u",
        "control_q", "control_jet_u", "all_checks",
    ]
    prime_fields = ["X", "p", "zero_count", "zeros"]
    rank_fields = [
        "dictionary", "degree", "split", "p", "samples", "columns",
        "rank", "nullity",
    ]
    relation_fields = [
        "dictionary", "degree", "polynomial", "support", "coefficients",
        "fail_X128", "fail_X256", "fail_X512", "fail_training",
        "fail_holdout", "fail_all", "fail_controls_all", "false_positive",
        "survives_all", "first_holdout_counterexample",
    ]

    write_csv(out / "rows.csv", row_rows, row_fields)
    write_csv(out / "hits.csv", hit_rows, hit_fields)
    write_csv(out / "primes.csv", prime_rows, prime_fields)
    write_csv(out / "ranks.csv", rank_rows, rank_fields)
    write_csv(out / "relations.csv", relation_rows, relation_fields)

    summary_payload = {
        "status": "PASS",
        "scope": {
            "X": list(xs),
            "training_max_X": TRAIN_MAX_X,
            "holdout_X": HOLDOUT_X,
        },
        "prereg_sha256": PREREG_SHA256,
        "selftest": selftest,
        "scales": scale_summaries,
        "total_high_load_rows": len(row_rows),
        "total_high_load_hits": len(hit_rows),
        "relation_screen": relation_summary,
    }
    (out / "summary.json").write_text(
        json.dumps(summary_payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    (out / "report.md").write_text(
        markdown_report(selftest, scale_summaries, row_rows, hit_rows, relation_summary),
        encoding="utf-8",
    )

    manifest_targets = [
        Path(__file__),
        out / "rows.csv",
        out / "hits.csv",
        out / "primes.csv",
        out / "ranks.csv",
        out / "relations.csv",
        out / "summary.json",
        out / "report.md",
    ]
    manifest_lines = []
    for path in manifest_targets:
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        manifest_lines.append(f"{digest}  {path.as_posix()}")
    (out / "MANIFEST.sha256").write_text("\n".join(manifest_lines) + "\n", encoding="utf-8")

    print(json.dumps({
        "status": "PASS",
        "rows": len(row_rows),
        "hits": len(hit_rows),
        "prereg_sha256": PREREG_SHA256,
        "output": str(out),
    }, sort_keys=True))


if __name__ == "__main__":
    main()
