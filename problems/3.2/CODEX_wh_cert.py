#!/usr/bin/env python3
"""Exact and modular certificates for CODEX_SPEC_wh_cert.md.

The h=2,...,6 definition chain is recomputed with SymPy over ZZ.  Exact
primitive value polynomials already computed by CRON_crit2h.py are audited and
factored over QQ through h=12.  For h=13,...,32, resultants are computed
directly over finite fields with Sage/FLINT; factor-degree subset-sum
intersections certify irreducibility when no single irreducible reduction is
found.
"""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import shutil
import sys
import time
import warnings
from typing import Any, Dict, Iterable, List, Sequence, Set, Tuple


ROOT = Path(__file__).resolve().parent
RESULT_JSON = ROOT / "CODEX_WH_CERT_results.json"
REPORT_MD = ROOT / "CODEX_WH_CERT_report.md"
PRIOR_EXACT_JSON = ROOT / "crit2h_results.json"
HEIGHT = 32
EXACT_THROUGH = 12
# Every prime is > 3*32.  The specification requires at least eight tested
# primes in the all-reducible case, so the complete eight-prime panel is used.
PRIMES = (101, 103, 107, 109, 113, 127, 131, 137)
MIN_TESTED_IF_NO_IRREDUCIBLE_REDUCTION = 8
CHECKPOINT_JSON = Path("/private/tmp/CODEX_WH_CERT_modular_checkpoint.json")


def ensure_sage() -> None:
    try:
        from sage.rings.finite_rings.finite_field_constructor import GF  # type: ignore  # noqa: F401
        return
    except ImportError:
        pass
    sage = shutil.which("sage")
    if sage is None:
        raise SystemExit("SageMath is required")
    if os.environ.get("WH_CERT_SAGE_REEXEC") == "1":
        raise SystemExit("sage -python re-exec did not expose sage.all")
    os.environ["WH_CERT_SAGE_REEXEC"] = "1"
    os.execvp(sage, [sage, "-python", str(Path(__file__).resolve()), *sys.argv[1:]])


ensure_sage()

from sage.rings.finite_rings.finite_field_constructor import GF  # type: ignore  # noqa: E402
from sage.rings.integer_ring import ZZ  # type: ignore  # noqa: E402
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing  # type: ignore  # noqa: E402


def log(message: str) -> None:
    print(message, flush=True)


def apery_P(z: Any) -> Any:
    return 34 * z**3 + 51 * z**2 + 27 * z + 5


def product(values: Iterable[Any], one: Any) -> Any:
    result = one
    for value in values:
        result *= value
    return result


def build_numerators(ring: Any, height: int) -> List[Any]:
    x = ring.gen()
    values = [ring.zero(), ring.one()]
    for h in range(1, height):
        values.append(
            apery_P(x + h) * values[h] - (x + h) ** 6 * values[h - 1]
        )
    return values


def q_polynomial(ring: Any, h: int) -> Any:
    x = ring.gen()
    return product((x + j for j in range(1, h + 1)), ring.one())


def digest_coefficients(poly: Any, modulus: int | None = None) -> str:
    if modulus is None:
        payload = ",".join(str(ZZ(poly[i])) for i in range(int(poly.degree()) + 1))
    else:
        payload = ",".join(
            str(int(poly[i]) % modulus) for i in range(int(poly.degree()) + 1)
        )
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def factor_degrees(poly: Any) -> List[int]:
    """Factor a Sage finite-field polynomial with SymPy's galoistools.

    This avoids loading Sage's global ``sage.all`` namespace, whose unrelated
    quaternion extension is broken in the current macOS Sage 10.9 bundle.
    """

    import sympy as sp

    z = sp.symbols("z")
    prime = int(poly.base_ring().characteristic())
    converted = sp.Poly(
        sum((int(poly[i]) % prime) * z**i for i in range(int(poly.degree()) + 1)),
        z,
        modulus=prime,
    )
    _unit, factors = sp.factor_list(converted)
    result: List[int] = []
    for factor, multiplicity in factors:
        result.extend([int(factor.degree())] * int(multiplicity))
    return sorted(result)


def subset_sums(degrees: Sequence[int]) -> Set[int]:
    sums = {0}
    for degree in degrees:
        sums |= {value + degree for value in tuple(sums)}
    return sums


def sympy_definition_chain() -> Tuple[List[Dict[str, Any]], Dict[int, List[int]]]:
    """Recompute every requested definition literally with SymPy for h<=6."""

    import sympy as sp

    x, T, U = sp.symbols("x T U")
    numerators = [sp.Poly(0, x, domain=sp.ZZ), sp.Poly(1, x, domain=sp.ZZ)]
    for h in range(1, 6):
        numerators.append(
            sp.Poly(apery_P(x + h), x, domain=sp.ZZ) * numerators[h]
            - sp.Poly((x + h) ** 6, x, domain=sp.ZZ) * numerators[h - 1]
        )

    rows: List[Dict[str, Any]] = []
    primitive_v_coefficients: Dict[int, List[int]] = {}
    for h in range(2, 7):
        started = time.monotonic()
        log(f"[sympy] h={h} begin exact resultant")
        N = numerators[h]
        q = sp.Poly(sp.prod(x + a for a in range(1, h + 1)), x, domain=sp.ZZ)
        D = q**3
        B = D * N.diff() - D.diff() * N
        C = q * N.diff() - 3 * q.diff() * N
        assert B == q**2 * C
        assert N.degree() == 3 * (h - 1)
        assert D.degree() == 3 * h
        assert C.degree() == 4 * h - 4

        fiber = sp.Poly(N.as_expr() - T * D.as_expr(), x, domain=sp.ZZ[T])
        raw = sp.Poly(sp.resultant(C.as_expr(), fiber.as_expr(), x), T, domain=sp.ZZ)
        content, primitive = raw.primitive()
        if primitive.LC() < 0:
            content, primitive = -content, -primitive
        assert raw == primitive.mul_ground(content)
        assert primitive.degree() == 4 * h - 4
        odd_coefficients_zero = all(
            primitive.nth(exponent) == 0
            for exponent in range(1, primitive.degree() + 1, 2)
        )
        assert odd_coefficients_zero
        W = sp.Poly(
            sum(primitive.nth(2 * k) * U**k for k in range(2 * h - 1)),
            U,
            domain=sp.ZZ,
        )
        assert W.degree() == 2 * h - 2
        assert sp.Poly(W.as_expr().subs(U, T**2), T, domain=sp.ZZ) == primitive
        primitive_v_coefficients[h] = [
            int(primitive.nth(k)) for k in range(primitive.degree() + 1)
        ]
        row = {
            "h": h,
            "N_degree": N.degree(),
            "D_degree": D.degree(),
            "C_degree": C.degree(),
            "quotient_identity": True,
            "V_degree": primitive.degree(),
            "V_even": True,
            "W_degree": W.degree(),
            "raw_content": str(content),
            "W_constant": str(W.nth(0)),
            "W_factor_degrees_Q": sorted(
                int(sp.degree(factor, U))
                for factor, _multiplicity in sp.factor_list(W.as_expr())[1]
            ),
            "W_squarefree_Q": sp.gcd(W, W.diff()).degree() == 0,
            "V_sha256": hashlib.sha256(
                ",".join(map(str, primitive_v_coefficients[h])).encode("ascii")
            ).hexdigest(),
            "elapsed_seconds": round(time.monotonic() - started, 3),
        }
        rows.append(row)
        log(f"[sympy] h={h} chain/evenness PASS ({row['elapsed_seconds']}s)")
    return rows, primitive_v_coefficients


def exact_certificates(
    sympy_coefficients: Dict[int, List[int]],
) -> Dict[int, Dict[str, Any]]:
    """Audit and factor the exact primitive V_h payload through h=12."""

    prior = json.loads(PRIOR_EXACT_JSON.read_text())
    exact_payload = prior["exact_V_coefficients"]
    u_ring = PolynomialRing(ZZ, "U")
    U = u_ring.gen()
    rows: Dict[int, Dict[str, Any]] = {}
    for h in range(2, EXACT_THROUGH + 1):
        entry = exact_payload[str(h)]
        coeffs = [ZZ(value) for value in entry["coefficients_ascending"]]
        assert len(coeffs) == 4 * h - 3
        assert all(coeffs[k] == 0 for k in range(1, len(coeffs), 2))
        if h in sympy_coefficients:
            assert [int(value) for value in coeffs] == sympy_coefficients[h]
        W = u_ring([coeffs[2 * k] for k in range(2 * h - 1)])
        assert int(W.degree()) == 2 * h - 2
        assert abs(ZZ(W.content())) == 1
        import sympy as sp

        z = sp.symbols("z")
        converted = sp.Poly(
            sum(int(W[k]) * z**k for k in range(int(W.degree()) + 1)),
            z,
            domain=sp.ZZ,
        )
        _unit, factors = sp.factor_list(converted)
        degrees = sorted(
            int(factor.degree())
            for factor, multiplicity in factors
            for _ in range(int(multiplicity))
        )
        gcd_degree = int(W.gcd(W.derivative()).degree())
        rows[h] = {
            "h": h,
            "degree": int(W.degree()),
            "constant": str(W[0]),
            "nonzero_at_zero": bool(W[0] != 0),
            "factor_degrees_Q": degrees,
            "irreducible_Q": degrees == [2 * h - 2],
            "squarefree_Q": gcd_degree == 0,
            "gcd_degree_Q": gcd_degree,
            "method": "direct factorization over QQ",
            "W_sha256": digest_coefficients(W),
        }
        log(f"[exact] h={h} factors={degrees} squarefree={gcd_degree == 0}")
    return rows


def modular_prime_rows(prime: int) -> Dict[int, Dict[str, Any]]:
    """Compute all resultants for one prime directly in F_prime."""

    assert ZZ(prime).is_prime() and prime > 3 * HEIGHT
    field = GF(prime)
    t_ring = PolynomialRing(field, "T")
    tx_ring = PolynomialRing(t_ring, "x")
    T = t_ring.gen()
    numerators = build_numerators(tx_ring, HEIGHT)
    rows: Dict[int, Dict[str, Any]] = {}
    log(f"[prime {prime}] begin")
    for h in range(EXACT_THROUGH + 1, HEIGHT + 1):
        started = time.monotonic()
        N = numerators[h]
        q = q_polynomial(tx_ring, h)
        C = q * N.derivative() - 3 * q.derivative() * N
        raw = t_ring(C.resultant(N - T * q**3))
        expected_v_degree = 4 * h - 4
        structural = {
            "ell_gt_3h": prime > 3 * h,
            "N_degree_preserved": int(N.degree()) == 3 * (h - 1),
            "C_degree_preserved": int(C.degree()) == expected_v_degree,
            "gcd_N_q_one": int(N.gcd(q).degree()) == 0,
            "q_squarefree": int(q.gcd(q.derivative()).degree()) == 0,
            "resultant_nonzero": raw != 0,
            "V_degree_preserved": raw != 0 and int(raw.degree()) == expected_v_degree,
        }
        even = raw != 0 and all(raw[k] == 0 for k in range(1, int(raw.degree()) + 1, 2))
        W = t_ring(
            [raw[2 * k] for k in range(2 * h - 1)] if even else [field.zero()]
        )
        structural.update(
            {
                "V_even": bool(even),
                "W_degree_preserved": even and int(W.degree()) == 2 * h - 2,
                "W_constant_nonzero": even and W[0] != 0,
            }
        )
        eligible = all(structural.values())
        degrees: List[int] = []
        gcd_degree: int | None = None
        if eligible:
            degrees = factor_degrees(W)
            assert sum(degrees) == 2 * h - 2
            gcd_degree = int(W.gcd(W.derivative()).degree())
        rows[h] = {
            "h": h,
            "prime": prime,
            "eligible": bool(eligible),
            "structural_checks": structural,
            "factor_degrees": degrees,
            "irreducible": eligible and degrees == [2 * h - 2],
            "gcd_degree": gcd_degree,
            "squarefree": eligible and gcd_degree == 0,
            "W_sha256_mod_ell": digest_coefficients(W.monic(), prime) if eligible else None,
            "elapsed_seconds": round(time.monotonic() - started, 3),
        }
        log(
            f"[mod] h={h} ell={prime} eligible={eligible} "
            f"pattern={degrees} gcddeg={gcd_degree} ({rows[h]['elapsed_seconds']}s)"
        )
    return rows


def combine_modular(all_prime_rows: Dict[int, Dict[int, Dict[str, Any]]]) -> Dict[int, Dict[str, Any]]:
    combined: Dict[int, Dict[str, Any]] = {}
    n_candidates = len(PRIMES)
    for h in range(EXACT_THROUGH + 1, HEIGHT + 1):
        n = 2 * h - 2
        tested = [all_prime_rows[p][h] for p in PRIMES]
        eligible = [row for row in tested if row["eligible"]]
        irreducible_rows = [row for row in eligible if row["irreducible"]]
        squarefree_rows = [row for row in eligible if row["squarefree"]]
        intersection: Set[int] = set(range(n + 1))
        intersection_trace: List[Dict[str, Any]] = []
        for row in eligible:
            if not row["squarefree"]:
                continue
            possible = subset_sums(row["factor_degrees"])
            intersection &= possible
            intersection_trace.append(
                {
                    "prime": row["prime"],
                    "factor_degrees": row["factor_degrees"],
                    "intersection": sorted(intersection),
                }
            )
        if irreducible_rows:
            certificate_method = "irreducible reduction"
            certificate_primes = [irreducible_rows[0]["prime"]]
            irreducible_Q = True
        else:
            # The specification requires at least eight tested primes before
            # using a degree-pattern intersection certificate.
            enough_tests = len(tested) >= MIN_TESTED_IF_NO_IRREDUCIBLE_REDUCTION
            irreducible_Q = enough_tests and intersection == {0, n}
            certificate_method = "degree-pattern intersection"
            certificate_primes = []
            for trace in intersection_trace:
                certificate_primes.append(trace["prime"])
                if trace["intersection"] == [0, n]:
                    break
        combined[h] = {
            "h": h,
            "degree": n,
            "primes_attempted": n_candidates,
            "primes_tested": list(PRIMES),
            "eligible_primes": [row["prime"] for row in eligible],
            "bad_reduction_primes": [row["prime"] for row in tested if not row["eligible"]],
            "patterns": {str(row["prime"]): row["factor_degrees"] for row in eligible},
            "irreducible_reduction_primes": [row["prime"] for row in irreducible_rows],
            "possible_degree_intersection": sorted(intersection),
            "intersection_trace": intersection_trace,
            "certificate_method": certificate_method,
            "certificate_primes": certificate_primes,
            "irreducible_Q": bool(irreducible_Q),
            "nonzero_at_zero": any(
                row["structural_checks"]["W_constant_nonzero"] for row in eligible
            ),
            "squarefree_Q": bool(squarefree_rows),
            "squarefree_certificate_prime": (
                squarefree_rows[0]["prime"] if squarefree_rows else None
            ),
        }
    return combined


def compact_pattern(pattern: Sequence[int]) -> str:
    counts: Dict[int, int] = {}
    for degree in pattern:
        counts[degree] = counts.get(degree, 0) + 1
    return "+".join(
        str(degree) if multiplicity == 1 else f"{degree}^{multiplicity}"
        for degree, multiplicity in sorted(counts.items())
    )


def render_report(data: Dict[str, Any]) -> str:
    symbolic = data["sympy_definition_chain"]
    exact = {int(h): row for h, row in data["exact_certificates"].items()}
    modular = {int(h): row for h, row in data["modular_certificates"].items()}
    all_rows = {**exact, **modular}
    failures = [
        h
        for h, row in sorted(all_rows.items())
        if not (row["nonzero_at_zero"] and row["irreducible_Q"] and row["squarefree_Q"])
    ]
    lines: List[str] = [
        "# W_h irreducibility certificates",
        "",
        "## Verdict",
        "",
    ]
    if failures:
        lines += [f"**FAILURES: h = {failures}.**", ""]
    else:
        lines += [
            "**PASS for every `2 <= h <= 32`.**  In this range `W_h(0) != 0`, "
            "`W_h` is irreducible over `Q`, and `gcd(W_h,W_h')=1`.",
            "",
            "No reducible `W_h` was found.",
            "",
        ]
    lines += [
        "This is an exact finite-range certificate, not a numerical-root test and not an all-h proof.",
        "",
        "## Definitions and symbolic calibration",
        "",
        "The computation uses the recurrence from `campaign3_questions/CTX.txt`:",
        "",
        "```text",
        "P(X)=34X^3+51X^2+27X+5",
        "N_1=1,  N_2=P(X+1)",
        "N_(h+1)=P(X+h)N_h-(X+h)^6 N_(h-1)",
        "q_h=product_(a=1)^h (X+a),  D_h=q_h^3",
        "C_h=q_h N_h' - 3q_h' N_h",
        "V_h(T)=Res_X(C_h,N_h-TD_h)",
        "V_h(T)=content_h W_h(T^2), with W_h primitive",
        "```",
        "",
        "SymPy recomputed the whole chain over `ZZ` for `h=2,...,6`.  It checked "
        "`D_h N_h'-D_h'N_h=q_h^2 C_h`, all expected degrees, exact resultant "
        "formation, vanishing of every odd coefficient of `V_h`, and the exact "
        "descent identity `W_h(T^2)=primitive_part(V_h(T))`.",
        "",
        "| h | deg N | deg D | deg C | deg V | deg W | V even | W(0) != 0 |",
        "|---:|---:|---:|---:|---:|---:|:---:|:---:|",
    ]
    for row in symbolic:
        lines.append(
            f"| {row['h']} | {row['N_degree']} | {row['D_degree']} | {row['C_degree']} | "
            f"{row['V_degree']} | {row['W_degree']} | yes | "
            f"{'yes' if int(row['W_constant']) != 0 else '**NO**'} |"
        )
    lines += [
        "",
        "The h=2,...,6 primitive `V_h` coefficient vectors agree coefficient-for-coefficient "
        "with the independently banked exact payload in `crit2h_results.json`.",
        "",
        "## Certificate logic",
        "",
        "For `h<=12`, the exact primitive integer `W_h` was factored directly over `Q`; "
        "the exact gcd with its derivative was also computed.",
        "",
        "For `h>12`, every resultant was formed directly in `F_ell[T]`, not by reducing "
        "a large characteristic-zero resultant.  A row was admitted only when "
        "`ell>3h`, the degrees of `N_h`, `C_h`, `V_h`, and `W_h` were preserved, "
        "`gcd(N_h,q_h)=1`, `q_h` was squarefree, `V_h` was even, and `W_h(0)` was nonzero. "
        "Thus the modular polynomial is a nonzero scalar multiple of the reduction of "
        "the primitive characteristic-zero `W_h`.",
        "",
        "An irreducible good reduction proves irreducibility over `Q`.  Otherwise, a "
        "factorization with degrees `d_1,...,d_r` restricts the degree of any rational "
        "factor to a subset sum of those degrees.  Intersecting these subset-sum sets "
        "over good squarefree reductions proves irreducibility when only `0` and "
        "`deg W_h` remain.  Eight primes were tested for every `h>12`, as required "
        "in the all-reducible case.  Squarefreeness over `Q` is "
        "independently certified by one squarefree, degree-preserving reduction.",
        "",
        "## Results",
        "",
        "| h | deg W | W(0) != 0 | irreducible over Q | certificate | squarefree |",
        "|---:|---:|:---:|:---:|---|:---:|",
    ]
    for h in range(2, HEIGHT + 1):
        row = all_rows[h]
        if h <= EXACT_THROUGH:
            certificate = "direct Q factorization"
        elif row["certificate_method"] == "irreducible reduction":
            certificate = f"irreducible mod {row['certificate_primes'][0]}"
        else:
            primes = ",".join(map(str, row["certificate_primes"]))
            certificate = f"degree patterns mod {primes}"
        lines.append(
            f"| {h} | {row['degree']} | {'yes' if row['nonzero_at_zero'] else '**NO**'} | "
            f"{'yes' if row['irreducible_Q'] else '**NO**'} | {certificate} | "
            f"{'yes' if row['squarefree_Q'] else '**NO**'} |"
        )
    lines += ["", "## Modular degree-pattern ledger", ""]
    for h in range(EXACT_THROUGH + 1, HEIGHT + 1):
        row = modular[h]
        patterns = "; ".join(
            f"{prime}: {compact_pattern(pattern)}"
            for prime, pattern in row["patterns"].items()
        )
        lines += [
            f"- `h={h}`: {patterns}.  Intersection = "
            f"`{row['possible_degree_intersection']}`; method = "
            f"`{row['certificate_method']}`; squarefree prime = "
            f"`{row['squarefree_certificate_prime']}`.",
        ]
    lines += [
        "",
        "## Reproducibility",
        "",
        "Run:",
        "",
        "```bash",
        "sage -python CODEX_wh_cert.py",
        "```",
        "",
        f"CAS versions for this run: SymPy {data['metadata']['sympy_version']}; "
        f"Sage {data['metadata']['sage_version']}.",
        "",
        "The full machine-readable record, including structural good-reduction gates, "
        "factor patterns, intersection traces, coefficient hashes, and timings, is in "
        "`CODEX_WH_CERT_results.json`.",
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    import sage.version  # type: ignore
    import sympy

    log("[startup] Sage/SymPy imports complete")
    assert all(ZZ(p).is_prime() for p in PRIMES)
    assert min(PRIMES) > 3 * HEIGHT
    symbolic, sympy_coefficients = sympy_definition_chain()
    exact = exact_certificates(sympy_coefficients)
    if CHECKPOINT_JSON.exists():
        checkpoint = json.loads(CHECKPOINT_JSON.read_text())
        prime_rows: Dict[int, Dict[int, Dict[str, Any]]] = {
            int(prime): {int(h): row for h, row in rows.items()}
            for prime, rows in checkpoint.items()
        }
        log(f"[checkpoint] resumed primes={sorted(prime_rows)}")
    else:
        prime_rows = {}
    for prime in PRIMES:
        if prime in prime_rows:
            continue
        prime_rows[prime] = modular_prime_rows(prime)
        CHECKPOINT_JSON.write_text(
            json.dumps(prime_rows, indent=2, sort_keys=True, default=int) + "\n"
        )
    modular = combine_modular(prime_rows)
    data = {
        "metadata": {
            "height": HEIGHT,
            "exact_through": EXACT_THROUGH,
            "primes": list(PRIMES),
            "minimum_tested_if_all_reducible": MIN_TESTED_IF_NO_IRREDUCIBLE_REDUCTION,
            "sympy_version": sympy.__version__,
            "sage_version": sage.version.version,
        },
        "sympy_definition_chain": symbolic,
        "exact_certificates": exact,
        "modular_prime_rows": prime_rows,
        "modular_certificates": modular,
    }
    failures = [
        h
        for h, row in {**exact, **modular}.items()
        if not (row["nonzero_at_zero"] and row["irreducible_Q"] and row["squarefree_Q"])
    ]
    data["failures"] = failures
    RESULT_JSON.write_text(
        json.dumps(data, indent=2, sort_keys=True, default=int) + "\n"
    )
    REPORT_MD.write_text(render_report(data))
    log(f"[done] failures={failures}; wrote {REPORT_MD.name} and {RESULT_JSON.name}")
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
