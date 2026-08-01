#!/usr/bin/env python3
"""Montes/p-adic factor-degree experiments for the polynomials M_h.

This script is deliberately self-contained.  It reproduces the exact SymPy
construction in ``CRON_Mh_galois.py``, runs the required sanity gates, searches
for finite-field Dedekind certificates for h = 15,...,40, obtains exact local
factor degrees from PARI/GP's ``factorpadic`` for p in {2,3,5,17} and
h = 2,...,24, and writes ``CODEX_MONTES_report.md``.

Run from this directory with

    python3 CRON_montes_pipeline.py
"""

from __future__ import annotations

import ast
import hashlib
import json
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from math import gcd
from pathlib import Path
from typing import Iterable, Sequence

import sympy
from sympy import Poly, Rational, ZZ, div, expand, factor_list, isprime, symbols
from sympy.ntheory.generate import primerange


X, T, Y = symbols("X T Y")

H_CERT_MIN = 15
H_CERT_MAX = 40
H_LOCAL_MIN = 2
H_LOCAL_MAX = 24
CERT_PRIME_BOUND = 50_000  # Strict bound: certificate primes satisfy p < 50000.
LOCAL_PRIMES = (2, 3, 5, 17)
PADIC_PRECISIONS = (80, 160)
EXPECTED_SANITY_DEGREES = (1, 3, 4, 6, 7, 9, 10, 12, 13, 15)
REPORT_PATH = Path(__file__).resolve().with_name("CODEX_MONTES_report.md")


class Progress:
    """Emit a timestamped progress line at least approximately every 30 s."""

    def __init__(self, interval: float = 30.0) -> None:
        self.started = time.monotonic()
        self.last = self.started

    def emit(self, message: str, *, force: bool = False) -> None:
        now = time.monotonic()
        if force or now - self.last >= 30.0:
            elapsed = now - self.started
            print(f"[{elapsed:8.1f}s] {message}", flush=True)
            self.last = now


def P_at(e: int):
    t = X + e
    return 34 * t**3 + 51 * t**2 + 27 * t + 5


def build_N(hmax: int, progress: Progress) -> dict[int, Poly]:
    """Construct N_h in ZZ[X] by the recurrence from CRON_Mh_galois.py."""
    N = {1: Poly(1, X, domain=ZZ), 2: Poly(P_at(1), X, domain=ZZ)}
    for d in range(2, hmax):
        N[d + 1] = (
            Poly(P_at(d), X, domain=ZZ) * N[d]
            - Poly((X + d) ** 6, X, domain=ZZ) * N[d - 1]
        )
        progress.emit(f"constructed N_{d + 1}")
    return N


def M_of(N: dict[int, Poly], h: int) -> Poly:
    """Construct primitive M_h in ZZ[Y], exactly matching CRON_Mh_galois.py."""
    f = N[h]
    if h % 2 == 0:
        q, r = div(f, Poly(2 * X + h + 1, X, domain=ZZ), domain="QQ")
        if not r.is_zero:
            raise RuntimeError(f"forced factor missing at h={h}")
        f = Poly(q, X)

    g = expand(f.as_expr().subs(X, (T - h - 1) / Rational(2)))
    gp = Poly(g, T)
    if any(c != 0 for (m,), c in gp.terms() if m % 2 == 1):
        raise RuntimeError(f"h={h}: transformed polynomial is not even in T")

    coeffs = {m // 2: c for (m,), c in gp.terms() if m % 2 == 0}
    cs = [coeffs.get(k, 0) for k in range(max(coeffs), -1, -1)]
    denominator_lcm = 1
    for c in cs:
        if hasattr(c, "q"):
            denominator_lcm = (
                denominator_lcm * c.q // gcd(denominator_lcm, c.q)
            )
    integral_coeffs = [int(c * denominator_lcm) for c in cs]
    content = 0
    for c in integral_coeffs:
        content = gcd(content, abs(c))
    if content == 0:
        raise RuntimeError(f"h={h}: zero polynomial encountered")
    return Poly([c // content for c in integral_coeffs], Y, domain=ZZ)


def polynomial_hash(poly: Poly) -> str:
    """SHA-256 of the canonical compact JSON list of descending coefficients."""
    coeffs = [int(c) for c in poly.all_coeffs()]
    payload = json.dumps(coeffs, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def is_irreducible_over_q(poly: Poly) -> bool:
    _, factors = factor_list(poly.as_expr())
    return len(factors) == 1 and factors[0][1] == 1


def run_sanity_gates(polys: dict[int, Poly], progress: Progress) -> None:
    degrees = tuple(polys[h].degree() for h in range(2, 12))
    if degrees != EXPECTED_SANITY_DEGREES:
        raise RuntimeError(
            "SANITY GATE FAILED: degrees h=2..11 are "
            f"{degrees}, expected {EXPECTED_SANITY_DEGREES}"
        )
    print(f"SANITY degree gate passed: {degrees}", flush=True)

    for h in range(2, 15):
        progress.emit(f"sanity irreducibility check h={h}")
        if not is_irreducible_over_q(polys[h]):
            raise RuntimeError(
                f"SANITY GATE FAILED: M_{h} is not irreducible over Q"
            )
    print("SANITY irreducibility gate passed for h=2..14", flush=True)


def modular_factor_degrees(poly: Poly, p: int) -> tuple[int, ...] | None:
    """Return squarefree-good factor degrees modulo p, or None if p is bad."""
    if int(poly.LC()) % p == 0:
        return None
    _, factors = factor_list(poly.as_expr(), modulus=p)
    if any(multiplicity != 1 for _, multiplicity in factors):
        return None
    degrees = tuple(sorted(Poly(factor, Y, modulus=p).degree() for factor, _ in factors))
    if sum(degrees) != poly.degree():
        raise RuntimeError(
            f"modular factor degrees at p={p} sum to {sum(degrees)}, "
            f"not {poly.degree()}"
        )
    return degrees


@dataclass(frozen=True)
class CertificateResult:
    polynomial_hash: str
    irreducible_prime: int | None
    primitive_cycle: tuple[int, int, tuple[int, ...]] | None
    transposition_cycle: tuple[int, tuple[int, ...]] | None
    primes_tested: int


def find_certificates(poly: Poly, progress: Progress, h: int) -> CertificateResult:
    """Search all p < 50000 until all three requested witnesses are found."""
    n = poly.degree()
    irreducible_prime: int | None = None
    primitive_cycle: tuple[int, int, tuple[int, ...]] | None = None
    transposition_cycle: tuple[int, tuple[int, ...]] | None = None
    tested = 0

    for p in primerange(2, CERT_PRIME_BOUND):
        tested += 1
        degrees = modular_factor_degrees(poly, p)
        if degrees is None:
            continue

        if irreducible_prime is None and degrees == (n,):
            irreducible_prime = p

        if transposition_cycle is None:
            even_parts = [d for d in degrees if d % 2 == 0]
            if even_parts == [2]:
                transposition_cycle = (p, degrees)

        if primitive_cycle is None:
            # This is the precise strong pattern used by CRON_Mh_galois.py:
            # one prime q > n/2 and every remaining cycle is fixed.
            for q in degrees:
                if (
                    isprime(q)
                    and q > n / 2
                    and degrees.count(q) == 1
                    and all(d == 1 or d == q for d in degrees)
                ):
                    primitive_cycle = (p, q, degrees)
                    break

        if irreducible_prime and primitive_cycle and transposition_cycle:
            break
        progress.emit(
            f"h={h}: finite-field search tested {tested} primes (latest p={p})"
        )

    return CertificateResult(
        polynomial_hash=polynomial_hash(poly),
        irreducible_prime=irreducible_prime,
        primitive_cycle=primitive_cycle,
        transposition_cycle=transposition_cycle,
        primes_tested=tested,
    )


def require_gp() -> str:
    gp = shutil.which("gp")
    if gp is None:
        raise RuntimeError(
            "PARI/GP is unavailable. This run refuses to guess local p-adic "
            "factorizations; install gp or implement a certified exact fallback."
        )
    return gp


def pari_version(gp: str) -> str:
    proc = subprocess.run(
        [gp, "-q"],
        input="print(version())\n",
        text=True,
        capture_output=True,
        check=True,
        timeout=30,
    )
    line = proc.stdout.strip()
    return line if line else "unknown"


def gp_factor_degrees(poly: Poly, p: int, precision: int, gp: str) -> tuple[int, ...]:
    """Return exact Q_p irreducible factor degrees using PARI factorpadic."""
    coeffs = ",".join(str(int(c)) for c in poly.all_coeffs())
    program = (
        f"f=Pol([{coeffs}]);"
        f"F=factorpadic(f,{p},{precision});"
        'print("__MONTES__",vector(matsize(F)[1],i,'
        "[poldegree(F[i,1]),F[i,2]]));\n"
    )
    proc = subprocess.run(
        [gp, "-q"],
        input=program,
        text=True,
        capture_output=True,
        timeout=300,
    )
    if proc.returncode != 0 or "***" in proc.stderr:
        raise RuntimeError(
            f"factorpadic failed for p={p}, degree={poly.degree()}, "
            f"precision={precision}:\n{proc.stderr.strip()}"
        )

    marker = "__MONTES__"
    lines = [line.strip() for line in proc.stdout.splitlines() if marker in line]
    if len(lines) != 1:
        raise RuntimeError(
            f"could not parse factorpadic output for p={p}: {proc.stdout!r}"
        )
    rows = ast.literal_eval(lines[0].split(marker, 1)[1])
    degrees: list[int] = []
    for degree, multiplicity in rows:
        degrees.extend([int(degree)] * int(multiplicity))
    degrees.sort()
    if sum(degrees) != poly.degree():
        raise RuntimeError(
            f"Q_{p} factor degrees {degrees} do not sum to {poly.degree()}"
        )
    return tuple(degrees)


def compute_local_degrees(
    polys: dict[int, Poly], progress: Progress, gp: str
) -> dict[tuple[int, int], tuple[int, ...]]:
    """Compute factor degrees twice and require precision stability."""
    result: dict[tuple[int, int], tuple[int, ...]] = {}
    low_precision, high_precision = PADIC_PRECISIONS
    for p in LOCAL_PRIMES:
        for h in range(H_LOCAL_MIN, H_LOCAL_MAX + 1):
            progress.emit(f"factorpadic Q_{p}, h={h}")
            low = gp_factor_degrees(polys[h], p, low_precision, gp)
            high = gp_factor_degrees(polys[h], p, high_precision, gp)
            if low != high:
                raise RuntimeError(
                    f"unstable factorpadic degrees for (p,h)=({p},{h}): "
                    f"precision {low_precision} gives {low}, "
                    f"precision {high_precision} gives {high}"
                )
            result[p, h] = high
        print(f"local factor degrees complete for p={p}", flush=True)
    return result


def subset_sums(degrees: Sequence[int]) -> set[int]:
    sums = {0}
    for degree in degrees:
        sums |= {s + degree for s in tuple(sums)}
    return sums


def local_intersections(
    polys: dict[int, Poly], local: dict[tuple[int, int], tuple[int, ...]]
) -> dict[int, tuple[int, ...]]:
    survivors: dict[int, tuple[int, ...]] = {}
    for h in range(H_LOCAL_MIN, H_LOCAL_MAX + 1):
        common = set.intersection(
            *(subset_sums(local[p, h]) for p in LOCAL_PRIMES)
        )
        n = polys[h].degree()
        if 0 not in common or n not in common:
            raise RuntimeError(f"trivial subset sums missing at h={h}: {common}")
        survivors[h] = tuple(sorted(common - {0, n}))
    return survivors


def fmt_degrees(degrees: Iterable[int]) -> str:
    return "[" + ", ".join(str(d) for d in degrees) + "]"


def fmt_primitive(cert: tuple[int, int, tuple[int, ...]] | None) -> str:
    if cert is None:
        return "NONE"
    p, q, degrees = cert
    return f"{p}; q={q}; {fmt_degrees(degrees)}"


def fmt_transposition(cert: tuple[int, tuple[int, ...]] | None) -> str:
    if cert is None:
        return "NONE"
    p, degrees = cert
    return f"{p}; {fmt_degrees(degrees)}"


def build_report(
    polys: dict[int, Poly],
    certificates: dict[int, CertificateResult],
    local: dict[tuple[int, int], tuple[int, ...]],
    survivors: dict[int, tuple[int, ...]],
    pari: str,
) -> str:
    proved = [h for h, degrees in survivors.items() if not degrees]
    failed = [h for h, degrees in survivors.items() if degrees]

    lines = [
        "# Montes/p-adic factor-degree pipeline for $M_h$",
        "",
        "## Method and sanity gates",
        "",
        "`CRON_montes_pipeline.py` reconstructs every polynomial with exact SymPy "
        "integer/rational arithmetic using the recurrence and `M_of` normalization "
        "from `CRON_Mh_galois.py`. The run reproduced degrees "
        "`1,3,4,6,7,9,10,12,13,15` for $h=2,\\ldots,11$ and exact "
        "factorization over $\\mathbf Q$ verified irreducibility for "
        "$h=2,\\ldots,14$ before the experiment proceeded.",
        "",
        f"Environment: SymPy `{sympy.__version__}`; PARI/GP version vector `{pari}`. "
        "Each local degree list was computed by PARI `factorpadic` from exact "
        "integer coefficients at precisions 80 and 160, with equality required "
        "between the two results.",
        "",
        "The polynomial hash is SHA-256 of the ASCII compact JSON array of all "
        "coefficients in descending degree order (for example, no whitespace: "
        "`[a_n,...,a_0]`). Finite-field searches use only primes $p<50000$ "
        "with nonzero leading coefficient and squarefree reduction.",
        "",
        "## Certificates for $h=15,\\ldots,40$",
        "",
        "An irreducible reduction at `p0` certifies irreducibility over "
        "$\\mathbf Q$. The primitive-cycle column uses the strong pattern from "
        "`CRON_Mh_galois.py`: one prime part $q>n/2$ and all remaining parts 1. "
        "The transposition column has exactly one even part, equal to 2, and all "
        "remaining parts odd. Entries give `prime; cycle data`.",
        "",
        "| h | $m_h$ | SHA-256 | p0 | prime-cycle certificate | transposition certificate |",
        "|---:|---:|:---|---:|:---|:---|",
    ]
    for h in range(H_CERT_MIN, H_CERT_MAX + 1):
        cert = certificates[h]
        p0 = str(cert.irreducible_prime) if cert.irreducible_prime else "NONE"
        primitive = (
            fmt_primitive(cert.primitive_cycle)
            if cert.irreducible_prime is not None
            else "—"
        )
        transposition = (
            fmt_transposition(cert.transposition_cycle)
            if cert.irreducible_prime is not None
            else "—"
        )
        lines.append(
            f"| {h} | {polys[h].degree()} | `{cert.polynomial_hash}` | {p0} | "
            f"{primitive} | {transposition} |"
        )

    lines.extend(
        [
            "",
            "## Exact local factor degrees",
            "",
            "Each list $D_{p,h}$ contains the degrees of the irreducible factors "
            "of $M_h$ over $\\mathbf Q_p$, with multiplicity.",
            "",
            "| h | $m_h$ | $D_{2,h}$ | $D_{3,h}$ | $D_{5,h}$ | $D_{17,h}$ |",
            "|---:|---:|:---|:---|:---|:---|",
        ]
    )
    for h in range(H_LOCAL_MIN, H_LOCAL_MAX + 1):
        cells = " | ".join(fmt_degrees(local[p, h]) for p in LOCAL_PRIMES)
        lines.append(f"| {h} | {polys[h].degree()} | {cells} |")

    lines.extend(
        [
            "",
            "## Mixed-prime subset-sum intersection",
            "",
            "A rational factor of degree $d$ must be a product of local "
            "irreducible factors at every prime, so $d$ must lie in every local "
            "subset-sum set. The table removes the trivial degrees 0 and $m_h$.",
            "",
            "| h | $m_h$ | surviving nontrivial degrees | verdict |",
            "|---:|---:|:---|:---|",
        ]
    )
    for h in range(H_LOCAL_MIN, H_LOCAL_MAX + 1):
        degrees = survivors[h]
        verdict = "irreducible" if not degrees else "not decided"
        lines.append(
            f"| {h} | {polys[h].degree()} | {fmt_degrees(degrees)} | {verdict} |"
        )

    failed_details = "; ".join(
        f"$h={h}$: {fmt_degrees(survivors[h])}" for h in failed
    )
    if not failed_details:
        failed_details = "none"
    answer = "yes" if not failed else "no"
    lines.extend(
        [
            "",
            "## Empirical answer",
            "",
            f"The four-prime intersection proves irreducibility for "
            f"**{len(proved)} of {H_LOCAL_MAX - H_LOCAL_MIN + 1}** values "
            f"$h=2,\\ldots,24$. It fails to decide {len(failed)} values: "
            f"{failed_details}.",
            "",
            f"Thus the empirical answer to whether the fixed mixed-prime set "
            f"$\\{{2,3,5,17\\}}$ is sufficient in practice through $h=24$ is "
            f"**{answer}**.",
            "",
            "## Reproduction",
            "",
            "```console",
            "$ python3 CRON_montes_pipeline.py",
            "```",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> int:
    progress = Progress()
    print("Montes pipeline started", flush=True)
    gp = require_gp()
    pari = pari_version(gp)
    print(f"using gp={gp}, PARI version={pari}", flush=True)

    N = build_N(H_CERT_MAX, progress)
    polys: dict[int, Poly] = {}
    for h in range(H_LOCAL_MIN, H_CERT_MAX + 1):
        polys[h] = M_of(N, h)
        progress.emit(f"constructed M_{h}, degree {polys[h].degree()}")
    print("constructed M_h for h=2..40", flush=True)

    run_sanity_gates(polys, progress)

    certificates: dict[int, CertificateResult] = {}
    for h in range(H_CERT_MIN, H_CERT_MAX + 1):
        print(
            f"certificate search h={h}, degree={polys[h].degree()}", flush=True
        )
        certificates[h] = find_certificates(polys[h], progress, h)
        c = certificates[h]
        print(
            f"  result h={h}: p0={c.irreducible_prime}, "
            f"primitive={c.primitive_cycle}, transposition={c.transposition_cycle}, "
            f"tested={c.primes_tested}",
            flush=True,
        )

    local = compute_local_degrees(polys, progress, gp)
    survivors = local_intersections(polys, local)
    report = build_report(polys, certificates, local, survivors, pari)
    REPORT_PATH.write_text(report, encoding="utf-8")
    print(f"wrote {REPORT_PATH}", flush=True)
    print("Montes pipeline completed", flush=True)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:
        print(f"FATAL: {exc}", file=sys.stderr, flush=True)
        raise
