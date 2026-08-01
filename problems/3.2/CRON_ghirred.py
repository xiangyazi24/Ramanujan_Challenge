#!/usr/bin/env python3
"""Exact absolute-irreducibility certificates for the collision curves G_h.

The default run covers h=2,...,40 and writes the two companion deliverables
``ghirred_certs.json`` and ``CODEX_GHIRRED_report.md`` next to this script.

The script deliberately has two independent exact backends for every h:

1. Singular's characteristic-zero ``absfact.lib`` computes the absolute
   factorization of G_h directly.
2. A finite-field critical-value certificate proves that the geometric
   monodromy of N_h/q_h^3 is S_(3h), hence that the off-diagonal collision
   curve G_h is geometrically irreducible.  Two primes are recorded by
   default as independent certificate instances.

Running this file with ordinary Python automatically re-executes it under
Sage, so both commands below are equivalent:

    python3 CRON_ghirred.py
    sage -python CRON_ghirred.py

All claimed factorizations and polynomial identities use exact arithmetic.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import os
import platform
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


def _enter_sage_if_needed():
    """Re-exec under ``sage -python`` when invoked by system Python."""

    try:
        import sage.all  # noqa: F401
        return
    except ImportError:
        if os.environ.get("GHIRRED_SAGE_BOOTSTRAPPED") == "1":
            raise

    sage = shutil.which("sage")
    if sage is None:
        raise RuntimeError(
            "SageMath is required. Install Sage or run this script in a "
            "Sage-enabled Python environment."
        )
    os.environ["GHIRRED_SAGE_BOOTSTRAPPED"] = "1"
    os.environ["GHIRRED_INVOKER_PYFLINT"] = str(
        importlib.util.find_spec("flint") is not None
    ).lower()
    os.execv(sage, [sage, "-python", os.path.abspath(__file__)] + sys.argv[1:])


_enter_sage_if_needed()

import sympy
from sympy import Poly as SympyPoly
from sympy import ZZ as SympyZZ
from sympy import diff as sympy_diff
from sympy import div as sympy_div
from sympy import expand as sympy_expand
from sympy import factor_list as sympy_factor_list
from sympy import resultant as sympy_resultant
from sympy import symbols as sympy_symbols

from sage.all import GF, QQ, ZZ, PolynomialRing, singular
from sage.env import SAGE_LOCAL
from sage.version import version as sage_version_string


SCRIPT_PATH = Path(__file__).resolve()
DEFAULT_JSON_PATH = SCRIPT_PATH.with_name("ghirred_certs.json")
DEFAULT_REPORT_PATH = SCRIPT_PATH.with_name("CODEX_GHIRRED_report.md")
SCHEMA = "ghirred-certificate-chain-v1"
CANONICAL_JSON = "UTF-8 compact JSON; sort_keys=true; separators=(',', ':')"


def canonical_bytes(value):
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    ).encode("ascii")


def sha256_canonical(value):
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def hashed_certificate(payload):
    return {"payload": payload, "sha256": sha256_canonical(payload)}


def p_polynomial(z):
    return 34 * z**3 + 51 * z**2 + 27 * z + 5


def integer_sparse_payload(poly):
    terms = []
    for exponent, coefficient in poly.dict().items():
        terms.append([int(exponent[0]), int(exponent[1]), int(coefficient)])
    terms.sort(key=lambda row: (-row[0], -row[1]))
    return terms


def integer_sparse_hash(poly):
    return sha256_canonical(integer_sparse_payload(poly))


def rational_sparse_hash(poly):
    terms = []
    for exponent, coefficient in poly.dict().items():
        terms.append([int(exponent[0]), int(exponent[1]), str(coefficient)])
    terms.sort(key=lambda row: (-row[0], -row[1]))
    return sha256_canonical(terms)


def integer_content(poly):
    content = 0
    for coefficient in poly.dict().values():
        content = math.gcd(content, abs(int(coefficient)))
    return content


def is_symmetric_sparse(poly):
    coefficients = {
        (int(exponent[0]), int(exponent[1])): coefficient
        for exponent, coefficient in poly.dict().items()
    }
    return all(
        coefficients.get((exponent[1], exponent[0]), 0) == coefficient
        for exponent, coefficient in coefficients.items()
    )


def sympy_sparse_hash(poly):
    terms = []
    for exponent, coefficient in poly.terms():
        terms.append([int(exponent[0]), int(exponent[1]), int(coefficient)])
    terms.sort(key=lambda row: (-row[0], -row[1]))
    return sha256_canonical(terms)


def sympy_reproduction(max_h=6):
    """Independently reconstruct and factor G_h for the banked h=2,...,6."""

    sx, sy, st = sympy_symbols("x y T")

    def sympy_p(z):
        return 34 * z**3 + 51 * z**2 + 27 * z + 5

    # This is exactly the K_m convention in the specification.
    k_values = {
        0: SympyPoly(1, sx, domain=SympyZZ),
        1: SympyPoly(sympy_p(sx), sx, domain=SympyZZ),
    }
    for m in range(1, max_h):
        k_values[m + 1] = SympyPoly(
            sympy_p(sx + m) * k_values[m].as_expr()
            - (sx + m) ** 6 * k_values[m - 1].as_expr(),
            sx,
            domain=SympyZZ,
        )

    records = {}
    exact_critical = {}
    for h in range(2, max_h + 1):
        started = time.perf_counter()
        numerator = SympyPoly(
            sympy_expand(k_values[h - 1].as_expr().subs(sx, sx + 1)),
            sx,
            domain=SympyZZ,
        )
        qh = SympyPoly(
            sympy_expand(math.prod(sx + j for j in range(1, h + 1))),
            sx,
            domain=SympyZZ,
        )
        r_expr = sympy_expand(
            numerator.as_expr() * qh.as_expr().subs(sx, sy) ** 3
            - numerator.as_expr().subs(sx, sy) * qh.as_expr() ** 3
        )
        r_poly = SympyPoly(r_expr, sx, sy, domain=SympyZZ)
        g_poly, remainder = sympy_div(
            r_poly,
            SympyPoly(sx - sy, sx, sy, domain=SympyZZ),
            domain=SympyZZ,
        )
        quotient_exact = remainder.is_zero
        coefficient, factors = sympy_factor_list(g_poly.as_expr())
        reconstruction = sympy_expand(coefficient)
        for factor, exponent in factors:
            reconstruction *= sympy_expand(factor) ** exponent
        reconstruction_ok = sympy_expand(reconstruction - g_poly.as_expr()) == 0
        single_factor = (
            len(factors) == 1
            and int(factors[0][1]) == 1
            and reconstruction_ok
        )

        record = {
            "h": h,
            "n_coefficients_descending": [
                int(value) for value in numerator.all_coeffs()
            ],
            "g_sha256": sympy_sparse_hash(g_poly),
            "deg_N": int(numerator.degree()),
            "deg_G_total": int(g_poly.total_degree()),
            "exact_division_remainder_zero": bool(quotient_exact),
            "factor_list_unit": str(coefficient),
            "factor_count": len(factors),
            "factor_exponents": [int(exponent) for _, exponent in factors],
            "single_Q_factor": bool(single_factor),
            "wall_seconds": time.perf_counter() - started,
        }
        records[h] = record

        # An exact ZZ baseline checks the modular resultant implementation.
        if h <= 4:
            critical = SympyPoly(
                sympy_resultant(
                    numerator.as_expr() - st * qh.as_expr() ** 3,
                    sympy_diff(numerator.as_expr(), sx) * qh.as_expr()
                    - 3
                    * numerator.as_expr()
                    * sympy_diff(qh.as_expr(), sx),
                    sx,
                ),
                st,
                domain=SympyZZ,
            )
            exact_critical[h] = [
                int(critical.nth(degree))
                for degree in range(critical.degree() + 1)
            ]

    return records, exact_critical


def singular_version_text():
    executable = shutil.which("Singular")
    if executable is None:
        sage_local = os.environ.get("SAGE_LOCAL") or str(SAGE_LOCAL)
        if sage_local:
            candidate = Path(sage_local) / "bin" / "Singular"
            if candidate.exists():
                executable = str(candidate)
    if executable is None:
        return "Singular version unavailable"
    completed = subprocess.run(
        [executable, "--version"],
        check=True,
        capture_output=True,
        text=True,
    )
    return "\n".join(completed.stdout.splitlines()[:4])


def singular_absolute_factorization(
    g_poly, rational_ring, verify_returned_associate=False
):
    """Run and verify Singular's exact characteristic-zero abs factorization."""

    started = time.perf_counter()
    source_ring = singular.ring(0, "(x,y)", "dp")
    source_ring.set_ring()
    singular_poly = singular(str(g_poly))
    output_ring = singular.absFactorizeBCG(singular_poly)
    output_ring.set_ring()

    factor_count = int(str(singular("absolute_factors[4]")))
    multiplicities = []
    returned_hash = None
    associate_unit = None
    associate_verified = None
    if factor_count == 1:
        multiplicities.append(int(str(singular("absolute_factors[2][2]"))))
        if verify_returned_associate:
            returned_text = str(singular("absolute_factors[1][2]"))
            returned_factor = rational_ring(returned_text)
            rational_g = rational_ring(g_poly)
            associate_unit_value = rational_g.lc() / returned_factor.lc()
            associate_verified = (
                rational_g == associate_unit_value * returned_factor
            )
            associate_unit = str(associate_unit_value)
            returned_hash = rational_sparse_hash(returned_factor)
    else:
        for index in range(2, factor_count + 2):
            multiplicities.append(
                int(str(singular("absolute_factors[2][%d]" % index)))
            )

    payload = {
        "method": "Singular_absFactorizeBCG_characteristic_zero",
        "absolute_factor_count": factor_count,
        "multiplicities": multiplicities,
        "single_factor_associate_verified": associate_verified,
        "associate_unit": associate_unit,
        "returned_factor_sha256": returned_hash,
    }
    return payload, time.perf_counter() - started


def singular_sanity_controls():
    """Check that the abs-factor count parser distinguishes a split conic."""

    results = {}
    for label, expression, expected in (
        ("split_over_Qbar", "x2+y2", 2),
        ("smooth_conic", "x2+y2+1", 1),
    ):
        source_ring = singular.ring(0, "(x,y)", "dp")
        source_ring.set_ring()
        output_ring = singular.absFactorizeBCG(singular(expression))
        output_ring.set_ring()
        actual = int(str(singular("absolute_factors[4]")))
        if actual != expected:
            raise RuntimeError(
                "Singular abs-factor sanity control %s returned %d, expected %d"
                % (label, actual, expected)
            )
        results[label] = actual
    return results


def reduce_integer_univariate(poly, target_ring):
    field = target_ring.base_ring()
    return target_ring([field(int(coefficient)) for coefficient in poly.list()])


def critical_value_certificate(h, numerator, qh, prime, exact_baseline=None):
    """Build a finite-field BG(h) certificate at one prime."""

    started = time.perf_counter()
    field = GF(prime)
    x_ring = PolynomialRing(field, "z")
    z = x_ring.gen()
    t_ring = PolynomialRing(field, "T")
    t = t_ring.gen()
    n_mod = reduce_integer_univariate(numerator, x_ring)
    q_mod = reduce_integer_univariate(qh, x_ring)
    critical_derivative = n_mod.derivative() * q_mod - 3 * n_mod * q_mod.derivative()
    denominator = q_mod**3
    expected_degree = 4 * h - 4

    gates = {
        "deg_N": int(n_mod.degree()),
        "expected_deg_N": 3 * h - 3,
        "deg_q": int(q_mod.degree()),
        "expected_deg_q": h,
        "gcd_N_q_degree": int(n_mod.gcd(q_mod).degree()),
        "deg_critical_derivative": int(critical_derivative.degree()),
        "expected_deg_critical_derivative": expected_degree,
        "gcd_critical_derivative_q_degree": int(
            critical_derivative.gcd(q_mod).degree()
        ),
    }

    values = []
    for value in range(1, expected_degree + 2):
        field_value = field(value)
        values.append(
            (
                field_value,
                (n_mod - field_value * denominator).resultant(
                    critical_derivative
                ),
            )
        )
    critical_values = t_ring.lagrange_polynomial(values)

    interpolation_checks = []
    for value in (
        expected_degree + 2,
        expected_degree + 3,
        expected_degree + 7,
    ):
        direct = (n_mod - field(value) * denominator).resultant(
            critical_derivative
        )
        matches = critical_values(field(value)) == direct
        interpolation_checks.append({"T": value, "matches": bool(matches)})

    # At T=0 the degree of N_h-T*q_h^3 drops by three.  The specialization
    # of the fixed-size Sylvester determinant therefore contains the three
    # missing leading rows as lc(A_h)^3.  This is the same correction used by
    # the exact interpolation construction in CODEX_irred_verify.py.
    c0_corrected = (
        critical_derivative.lc() ** 3
        * n_mod.resultant(critical_derivative)
    )
    c0_matches = critical_values(field(0)) == c0_corrected
    gcd_degree = int(critical_values.gcd(critical_values.derivative()).degree())
    coefficients = [
        int(critical_values[degree])
        for degree in range(expected_degree + 1)
    ]

    baseline_matches = None
    if exact_baseline is not None:
        reduced_baseline = [value % prime for value in exact_baseline]
        baseline_matches = reduced_baseline == coefficients

    gates.update(
        {
            "deg_C": int(critical_values.degree()),
            "expected_deg_C": expected_degree,
            "C0_nonzero": bool(critical_values[0] != 0),
            "C0_degree_drop_corrected_resultant_matches": bool(c0_matches),
            "gcd_C_Cprime_degree": gcd_degree,
            "interpolation_checks": interpolation_checks,
            "exact_ZZ_baseline_matches": baseline_matches,
        }
    )
    all_gates = (
        gates["deg_N"] == gates["expected_deg_N"]
        and gates["deg_q"] == gates["expected_deg_q"]
        and gates["gcd_N_q_degree"] == 0
        and gates["deg_critical_derivative"]
        == gates["expected_deg_critical_derivative"]
        and gates["gcd_critical_derivative_q_degree"] == 0
        and gates["deg_C"] == gates["expected_deg_C"]
        and gates["C0_nonzero"]
        and gates["C0_degree_drop_corrected_resultant_matches"]
        and gates["gcd_C_Cprime_degree"] == 0
        and all(check["matches"] for check in interpolation_checks)
        and baseline_matches is not False
    )

    payload = {
        "method": "critical_value_BG_geometric_monodromy",
        "prime": int(prime),
        "gates": gates,
        "C_coefficients_ascending_mod_q": coefficients,
        "criterion_conclusion": (
            "GeomMon(N_h/q_h^3)=S_(3h); G_h is absolutely irreducible"
        ),
        "all_gates_passed": bool(all_gates),
    }
    backend = "%s.%s" % (
        type(z).__module__,
        type(z).__name__,
    )
    return payload, time.perf_counter() - started, backend


def verify_hashed_certificate(certificate):
    return certificate.get("sha256") == sha256_canonical(certificate.get("payload"))


def verify_certificate_tree(data):
    for entry in data["heights"]:
        certificates = [entry["direct_absolute_factorization"]]
        certificates.extend(entry["critical_value_certificates"])
        if entry.get("sympy_Q_reproduction") is not None:
            certificates.append(entry["sympy_Q_reproduction"])
        if not all(verify_hashed_certificate(cert) for cert in certificates):
            raise RuntimeError("certificate hash mismatch at h=%d" % entry["h"])
        bundle_payload = {
            "h": entry["h"],
            "g_sha256": entry["g_sha256"],
            "certificate_sha256": [cert["sha256"] for cert in certificates],
        }
        if entry["bundle_sha256"] != sha256_canonical(bundle_payload):
            raise RuntimeError("bundle hash mismatch at h=%d" % entry["h"])
    expected_chain = sha256_canonical(
        [entry["bundle_sha256"] for entry in data["heights"]]
    )
    if data["chain_sha256"] != expected_chain:
        raise RuntimeError("global chain hash mismatch")


def build_report(data, json_file_sha256):
    complete = data["status"] == "COMPLETE"
    h_min, h_max = data["range"]
    lines = [
        "# Absolute irreducibility certificates for $G_h$, $2\\le h\\le40$",
        "",
        "## Verdict",
        "",
    ]
    if complete:
        lines.extend(
            [
                "The contiguous chain is **complete for every "
                f"$h={h_min},\\ldots,{h_max}$**. For each height, two "
                "independent exact methods certify that $G_h$ is absolutely "
                "irreducible over $\\overline{\\mathbf Q}$:",
                "",
                "1. direct characteristic-zero absolute factorization by "
                "Singular `absfact.lib`;",
                "2. the critical-value/monodromy criterion, certified modulo "
                "two independent good primes.",
                "",
                "Therefore the primitive integral model of $G_h$ remains "
                "absolutely irreducible modulo $p$ for all but finitely many "
                "primes $p$, for every fixed $2\\le h\\le40$.",
            ]
        )
    else:
        lines.append(
            "The run produced a partial chain. See the anomaly section for "
            "the exact failed heights and diagnostics."
        )

    lines.extend(
        [
            "",
            "## Convention and $h=2,\\ldots,6$ reproduction",
            "",
            "The certified convention is exactly the one in the specification "
            "and ledger appendix AS.4:",
            "",
            "```text",
            "N_h(r) = K_(h-1)(r+1)",
            "R_h(x,y) = N_h(x) prod_(j=1)^h (y+j)^3",
            "           - N_h(y) prod_(j=1)^h (x+j)^3",
            "G_h(x,y) = R_h(x,y)/(x-y).",
            "```",
            "",
            "For even $h$, the forced rational factor "
            "$(2X+h+1)$ of $N_h$ is **retained**. No cofactor normalization "
            "is made. Exact SymPy construction from the $K_m$ recurrence was "
            "compared coefficient-for-coefficient with the Sage construction.",
            "",
            "| h | deg $N_h$ | deg $G_h$ | exact remainder | "
            "`factor_list` over $\\mathbf Q$ | convention match | wall (s) |",
            "|---:|---:|---:|:---:|:---:|:---:|---:|",
        ]
    )
    reproduction_by_h = {item["h"]: item for item in data["reproduction"]}
    for h in range(2, 7):
        item = reproduction_by_h[h]
        lines.append(
            "| {h} | {dn} | {dg} | {rem} | {fac} | {match} | {wall:.3f} |".format(
                h=h,
                dn=item["deg_N"],
                dg=item["deg_G_total"],
                rem="0" if item["exact_division_remainder_zero"] else "NONZERO",
                fac=(
                    "one factor"
                    if item["single_Q_factor"]
                    else "%d factors" % item["factor_count"]
                ),
                match="yes" if item["sage_convention_match"] else "NO",
                wall=item["wall_seconds"],
            )
        )

    lines.extend(
        [
            "",
            "## Exact criteria",
            "",
            "### Method AF: direct absolute factorization",
            "",
            "Singular 4.4.1's `absFactorizeBCG` works over characteristic zero "
            "and returns all absolute factors, grouped with their conjugates. "
            "For every $h=2,\\ldots,40$, `absolute_factors[4]` is 1, its "
            "multiplicity is 1. This directly proves absolute irreducibility "
            "and hence $\\mathbf Q$-irreducibility. For $h\\le6$, the returned "
            "factor was additionally parsed back over $\\mathbf Q[x,y]$ and "
            "verified coefficientwise to be an associate of the exact $G_h$.",
            "",
            "Two controls guard the interface/parser: $x^2+y^2$ returns two "
            "absolute factors, while $x^2+y^2+1$ returns one.",
            "",
            "### Method BG: critical values and geometric monodromy",
            "",
            "Put $q_h(X)=\\prod_{j=1}^h(X+j)$, "
            "$\\delta_h=N_h/q_h^3$, and",
            "",
            "$$A_h=N_h'q_h-3N_hq_h',\\qquad "
            "C_h(T)=\\operatorname{Res}_X(N_h-Tq_h^3,A_h).$$",
            "",
            "For a listed prime $q$, the certificate stores every coefficient "
            "of $C_h\\bmod q$ and checks",
            "",
            "$$\\deg C_h=4h-4,\\quad C_h(0)\\ne0,\\quad "
            "\\gcd(C_h,C_h')=1,$$",
            "",
            "together with all construction degrees, pole coprimality, three "
            "out-of-sample resultant evaluations, and (for $h\\le4$) an "
            "independent exact-$\\mathbf Z$ resultant reduction. Expected "
            "degree modulo one prime prevents leading-term loss, so these "
            "checks prove the same nonvanishing statements in characteristic "
            "zero.",
            "",
            "Here is the exact implication used. The checks give the branch "
            "profile of the degree-$3h$ map $\\delta_h$: an index-3 zero at "
            "infinity, $h$ index-3 poles, and $4h-4$ simple ramification "
            "points with pairwise distinct nonzero critical values. If "
            "$\\delta_h=g\\circ u$ with $\\deg u=a>1$, every one of those "
            "$4h-4$ simple points must ramify in $u$, hence "
            "$a\\ge2h-1$ by Riemann--Hurwitz; but the nontrivial outer factor "
            "gives $a\\le3h/2$. This contradicts $h\\ge3$, and the $h=2$ "
            "endpoint is ruled out by its index-3 zero. Thus the geometric "
            "monodromy is primitive. A simple branch point supplies a "
            "transposition, so the monodromy group is $S_{3h}$. Its action "
            "on ordered unequal pairs is transitive; consequently the "
            "off-diagonal factor $G_h$ is geometrically irreducible.",
            "",
            "This is an exact criterion (method (c) in the specification), "
            "not a point-count heuristic.",
            "",
            "## Per-height certificate chain",
            "",
            "`AF` is the direct Singular certificate; `BG(q)` is the "
            "critical-value certificate at prime $q$; `SQ` is the independent "
            "SymPy $\\mathbf Q$-factorization recorded for $h\\le6$. Hashes "
            "are SHA-256 of each certificate's canonical compact JSON payload.",
            "",
            "| h | deg $G_h$ | methods | primes | SHA-256 certificates | wall (s) |",
            "|---:|---:|:---|:---|:---|---:|",
        ]
    )

    for entry in data["heights"]:
        direct = entry["direct_absolute_factorization"]
        bg = entry["critical_value_certificates"]
        methods = "AF + BGx%d" % len(bg)
        hashes = ["AF=`%s`" % direct["sha256"]]
        primes = []
        for cert in bg:
            q = cert["payload"]["prime"]
            primes.append(str(q))
            hashes.append("BG(%d)=`%s`" % (q, cert["sha256"]))
        if entry.get("sympy_Q_reproduction") is not None:
            methods += " + SQ"
            hashes.append("SQ=`%s`" % entry["sympy_Q_reproduction"]["sha256"])
        lines.append(
            "| {h} | {degree} | {methods} | {primes} | {hashes} | {wall:.3f} |".format(
                h=entry["h"],
                degree=entry["degrees"]["G_total"],
                methods=methods,
                primes=", ".join(primes),
                hashes="<br>".join(hashes),
                wall=entry["timing_seconds"]["total"],
            )
        )

    lines.extend(["", "## Anomalies", ""])
    if data["anomalies"]:
        for anomaly in data["anomalies"]:
            lines.append("- `%s`" % json.dumps(anomaly, sort_keys=True))
    else:
        lines.append(
            "None. No factorization, convention mismatch, or failed "
            "certificate gate occurred."
        )

    env = data["environment"]
    pyflint_state = env["invoker_python_flint_importable"]
    if pyflint_state == "false":
        backend_sentence = (
            "The invoking system Python did not have `python-flint`; the run "
            "therefore used Sage's FLINT-backed prime-field univariate "
            "polynomials for modular resultants/interpolation and Singular's "
            "exact `absfact.lib` for bivariate absolute factorization."
        )
    else:
        backend_sentence = (
            "The run used Sage's FLINT-backed prime-field univariate "
            "polynomials for modular resultants/interpolation and Singular's "
            "exact `absfact.lib` for bivariate absolute factorization."
        )
    lines.extend(
        [
            "",
            "## Environment and reproducibility",
            "",
            backend_sentence + " No floating-point arithmetic enters a certificate.",
            "",
            "- Sage: `%s`" % env["sage"],
            "- SymPy: `%s`" % env["sympy"],
            "- modular polynomial backend: `%s`" % env["modular_backend"],
            "- invoking Python `flint` importable: `%s`"
            % env["invoker_python_flint_importable"],
            "- Singular:",
            "",
            "```text",
            env["singular"],
            "```",
            "",
            "Canonical certificate encoding: `%s`." % CANONICAL_JSON,
            "",
            "- chain SHA-256: `%s`" % data["chain_sha256"],
            "- generated JSON file SHA-256: `%s`" % json_file_sha256,
            "",
            "Reproduce from this directory with:",
            "",
            "```console",
            "$ python3 CRON_ghirred.py",
            "```",
            "",
            "The script prints progress for every height, rewrites both "
            "deliverables, rereads the JSON, and verifies every stored "
            "certificate and bundle hash before returning success.",
            "",
        ]
    )
    return "\n".join(lines)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--hmin", type=int, default=2)
    parser.add_argument("--hmax", type=int, default=40)
    parser.add_argument(
        "--primes",
        type=int,
        nargs="+",
        default=[65537, 1000003],
        help="good-prime candidates for independent BG certificates",
    )
    parser.add_argument("--json", type=Path, default=DEFAULT_JSON_PATH)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT_PATH)
    return parser.parse_args()


def main():
    args = parse_args()
    if args.hmin != 2:
        raise ValueError("the certificate chain must start at h=2")
    if args.hmax < 6:
        raise ValueError("hmax must be at least 6 for the reproduction gate")
    if len(args.primes) < 1:
        raise ValueError("at least one BG certificate prime is required")
    if any(not GF(prime).is_prime_field() for prime in args.primes):
        raise ValueError("every --primes entry must be prime")
    if any(prime <= 4 * args.hmax + 7 for prime in args.primes):
        raise ValueError("certificate primes must exceed 4*hmax+7")

    run_started = time.perf_counter()
    print("G_h certificate pipeline started", flush=True)
    print("independent SymPy reproduction h=2..6", flush=True)
    reproduction, exact_critical = sympy_reproduction(6)

    singular.lib("absfact.lib")
    controls = singular_sanity_controls()
    print("Singular absolute-factor controls passed: %s" % controls, flush=True)

    u_ring = PolynomialRing(ZZ, "u")
    u = u_ring.gen()
    xy_ring = PolynomialRing(ZZ, names=("x", "y"), order="degrevlex")
    x, y = xy_ring.gens()
    rational_xy_ring = PolynomialRing(
        QQ, names=("x", "y"), order="degrevlex"
    )
    map_x = u_ring.hom([x], xy_ring)
    map_y = u_ring.hom([y], xy_ring)

    numerator_previous = u_ring(0)
    numerator = u_ring(1)
    qh = u_ring(1)
    heights = []
    anomalies = []
    modular_backend = None

    for h in range(1, args.hmax + 1):
        if h > 1:
            next_numerator = (
                p_polynomial(u + h - 1) * numerator
                - (u + h - 1) ** 6 * numerator_previous
            )
            numerator_previous, numerator = numerator, next_numerator
        qh *= u + h
        if h < args.hmin:
            continue

        height_started = time.perf_counter()
        print("h=%2d: construct G_h" % h, flush=True)
        entry = {"h": h}
        height_anomalies = []

        construction_started = time.perf_counter()
        nx = map_x(numerator)
        ny = map_y(numerator)
        qx_cubed = map_x(qh) ** 3
        qy_cubed = map_y(qh) ** 3
        rh = nx * qy_cubed - ny * qx_cubed
        gh, remainder = rh.quo_rem(x - y)
        exact_division = remainder == 0
        symmetric = is_symmetric_sparse(gh)
        content = integer_content(gh)
        g_sha256 = integer_sparse_hash(gh)
        forced_factor_divides = None
        if h % 2 == 0:
            _, forced_remainder = numerator.quo_rem(2 * u + h + 1)
            forced_factor_divides = forced_remainder == 0
        construction_wall = time.perf_counter() - construction_started

        entry.update(
            {
                "degrees": {
                    "N": int(numerator.degree()),
                    "G_x": int(gh.degree(x)),
                    "G_y": int(gh.degree(y)),
                    "G_total": int(gh.total_degree()),
                },
                "g_monomial_count": len(gh.dict()),
                "g_content": content,
                "g_sha256": g_sha256,
                "exact_division_remainder_zero": bool(exact_division),
                "G_symmetric": bool(symmetric),
                "even_h_forced_N_factor_retained": bool(h % 2 == 0),
                "even_h_forced_factor_divides_N": forced_factor_divides,
            }
        )

        expected_degrees = (
            entry["degrees"]["N"] == 3 * h - 3
            and entry["degrees"]["G_x"] == 3 * h - 1
            and entry["degrees"]["G_y"] == 3 * h - 1
            and entry["degrees"]["G_total"] == 6 * h - 4
        )
        if (
            not exact_division
            or not symmetric
            or not expected_degrees
            or forced_factor_divides is False
        ):
            height_anomalies.append(
                {
                    "stage": "construction",
                    "exact_division": exact_division,
                    "symmetric": symmetric,
                    "degrees": entry["degrees"],
                    "forced_factor_divides": forced_factor_divides,
                }
            )

        sympy_certificate = None
        if h <= 6:
            sympy_record = reproduction[h]
            numerator_match = sympy_record["n_coefficients_descending"] == [
                int(value) for value in reversed(numerator.list())
            ]
            convention_match = (
                numerator_match and sympy_record["g_sha256"] == g_sha256
            )
            sympy_record["sage_convention_match"] = bool(convention_match)
            sympy_payload = {
                "method": "SymPy_factor_list_over_Q",
                "h": h,
                "g_sha256": g_sha256,
                "exact_division_remainder_zero": sympy_record[
                    "exact_division_remainder_zero"
                ],
                "factor_count": sympy_record["factor_count"],
                "factor_exponents": sympy_record["factor_exponents"],
                "single_Q_factor": sympy_record["single_Q_factor"],
                "independent_construction_matches": bool(convention_match),
            }
            sympy_certificate = hashed_certificate(sympy_payload)
            if not (
                convention_match
                and sympy_record["exact_division_remainder_zero"]
                and sympy_record["single_Q_factor"]
            ):
                height_anomalies.append(
                    {"stage": "h2_h6_reproduction", "record": sympy_record}
                )
        entry["sympy_Q_reproduction"] = sympy_certificate

        print("h=%2d: Singular direct absolute factorization" % h, flush=True)
        direct_started = time.perf_counter()
        try:
            direct_payload, direct_wall = singular_absolute_factorization(
                gh, rational_xy_ring, verify_returned_associate=(h <= 6)
            )
        except Exception as exc:
            direct_wall = time.perf_counter() - direct_started
            direct_payload = {
                "method": "Singular_absFactorizeBCG_characteristic_zero",
                "absolute_factor_count": None,
                "multiplicities": [],
                "single_factor_associate_verified": None,
                "associate_unit": None,
                "returned_factor_sha256": None,
                "error": repr(exc),
            }
        direct_payload.update(
            {
                "h": h,
                "g_sha256": g_sha256,
                "g_total_degree": entry["degrees"]["G_total"],
            }
        )
        direct_certificate = hashed_certificate(direct_payload)
        direct_pass = (
            direct_payload["absolute_factor_count"] == 1
            and direct_payload["multiplicities"] == [1]
            and direct_payload["single_factor_associate_verified"] is not False
        )
        if not direct_pass:
            height_anomalies.append(
                {"stage": "direct_absolute_factorization", "result": direct_payload}
            )
        entry["direct_absolute_factorization"] = direct_certificate

        bg_certificates = []
        bg_wall = 0.0
        for prime in args.primes:
            print("h=%2d: BG certificate q=%d" % (h, prime), flush=True)
            baseline = exact_critical.get(h)
            bg_started = time.perf_counter()
            try:
                bg_payload, elapsed, backend = critical_value_certificate(
                    h, numerator, qh, prime, baseline
                )
                modular_backend = backend
            except Exception as exc:
                elapsed = time.perf_counter() - bg_started
                bg_payload = {
                    "method": "critical_value_BG_geometric_monodromy",
                    "prime": int(prime),
                    "all_gates_passed": False,
                    "error": repr(exc),
                }
            bg_payload.update(
                {
                    "h": h,
                    "g_sha256": g_sha256,
                }
            )
            bg_certificate = hashed_certificate(bg_payload)
            bg_certificates.append(bg_certificate)
            bg_wall += elapsed
            if not bg_payload["all_gates_passed"]:
                height_anomalies.append(
                    {
                        "stage": "critical_value_BG",
                        "prime": prime,
                        "gates": bg_payload.get("gates"),
                        "error": bg_payload.get("error"),
                    }
                )
        entry["critical_value_certificates"] = bg_certificates

        certificates = [direct_certificate] + bg_certificates
        if sympy_certificate is not None:
            certificates.append(sympy_certificate)
        bundle_payload = {
            "h": h,
            "g_sha256": g_sha256,
            "certificate_sha256": [cert["sha256"] for cert in certificates],
        }
        entry["bundle_sha256"] = sha256_canonical(bundle_payload)
        entry["timing_seconds"] = {
            "construction": construction_wall,
            "direct_absolute_factorization": direct_wall,
            "critical_value_certificates": bg_wall,
            "total": time.perf_counter() - height_started,
        }
        entry["status"] = "CERTIFIED" if not height_anomalies else "ANOMALY"
        heights.append(entry)
        for anomaly in height_anomalies:
            anomaly_with_h = {"h": h}
            anomaly_with_h.update(anomaly)
            anomalies.append(anomaly_with_h)

        print(
            "h=%2d: %s deg=%d monomials=%d wall=%.3fs bundle=%s"
            % (
                h,
                entry["status"],
                entry["degrees"]["G_total"],
                entry["g_monomial_count"],
                entry["timing_seconds"]["total"],
                entry["bundle_sha256"][:16],
            ),
            flush=True,
        )

    status = (
        "COMPLETE"
        if not anomalies
        and len(heights) == args.hmax - args.hmin + 1
        and all(entry["status"] == "CERTIFIED" for entry in heights)
        else "PARTIAL"
    )
    data = {
        "schema": SCHEMA,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": status,
        "range": [args.hmin, args.hmax],
        "convention": {
            "N_h": "K_(h-1)(r+1)",
            "R_h": "N_h(x)*prod(y+j)^3-N_h(y)*prod(x+j)^3",
            "G_h": "R_h/(x-y)",
            "even_h_forced_factor_of_N_h": "retained",
        },
        "canonical_certificate_encoding": CANONICAL_JSON,
        "environment": {
            "python": platform.python_version(),
            "sage": str(sage_version_string),
            "sympy": sympy.__version__,
            "singular": singular_version_text(),
            "modular_backend": modular_backend,
            "invoker_python_flint_importable": os.environ.get(
                "GHIRRED_INVOKER_PYFLINT", "unknown"
            ),
        },
        "singular_sanity_controls": controls,
        "reproduction": [reproduction[h] for h in range(2, 7)],
        "heights": heights,
        "anomalies": anomalies,
        "total_wall_seconds": time.perf_counter() - run_started,
    }
    data["chain_sha256"] = sha256_canonical(
        [entry["bundle_sha256"] for entry in heights]
    )
    verify_certificate_tree(data)

    args.json.write_text(
        json.dumps(data, indent=2, sort_keys=True, ensure_ascii=True) + "\n",
        encoding="utf-8",
    )
    reread = json.loads(args.json.read_text(encoding="utf-8"))
    verify_certificate_tree(reread)
    json_digest = hashlib.sha256(args.json.read_bytes()).hexdigest()
    args.report.write_text(build_report(reread, json_digest), encoding="utf-8")

    print("wrote %s" % args.json, flush=True)
    print("wrote %s" % args.report, flush=True)
    print("status=%s chain_sha256=%s" % (status, data["chain_sha256"]), flush=True)
    return 0 if status == "COMPLETE" else 1


if __name__ == "__main__":
    sys.exit(main())
