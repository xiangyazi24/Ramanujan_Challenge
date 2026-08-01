#!/usr/bin/env python3
"""Verification gates for CODEX_MAINTHM_report.md.

The script is deliberately independent of the exploratory scripts in this
directory.  It recomputes every finite datum quoted in the report and checks
the exact algebraic identities on many live primes/indices.  The asymptotic
statements in the report are proved there; the corresponding exponent and
finite-identity bookkeeping is checked here.

Expected runtime on the development machine is below one minute.  The largest
object is the 3001 by 3001 FFT used to extend the vector-Weyl experiment.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
import cmath
import math
import sys
import traceback

import numpy as np


class GateFailure(RuntimeError):
    pass


def check(condition: bool, message: str) -> None:
    if not condition:
        raise GateFailure(message)


def primes_upto(limit: int) -> list[int]:
    sieve = bytearray(b"\x01") * (limit + 1)
    if limit >= 0:
        sieve[0] = 0
    if limit >= 1:
        sieve[1] = 0
    for p in range(2, math.isqrt(limit) + 1):
        if sieve[p]:
            sieve[p * p : limit + 1 : p] = b"\x00" * (
                (limit - p * p) // p + 1
            )
    return [p for p in range(2, limit + 1) if sieve[p]]


def valuation(value: int, prime: int) -> int:
    check(value != 0, "valuation of zero requested")
    value = abs(value)
    answer = 0
    while value % prime == 0:
        value //= prime
        answer += 1
    return answer


def apery_coefficient(n: int) -> int:
    return 34 * n**3 + 51 * n**2 + 27 * n + 5


def exact_apery_sequences(limit: int) -> tuple[list[int], list[Fraction]]:
    b = [1, 5]
    a = [Fraction(0), Fraction(6)]
    for n in range(1, limit):
        den = (n + 1) ** 3
        b_num = apery_coefficient(n) * b[n] - n**3 * b[n - 1]
        check(b_num % den == 0, f"nonintegral b at n={n + 1}")
        b.append(b_num // den)
        a.append(
            (apery_coefficient(n) * a[n] - n**3 * a[n - 1]) / den
        )
    return b, a


def denominator_scales(limit: int) -> list[int]:
    scales = [1] * (limit + 1)
    ell = 1
    for n in range(1, limit + 1):
        ell = math.lcm(ell, n)
        scales[n] = ell**3
    return scales


def fraction_mod(value: Fraction, prime: int) -> int:
    check(value.denominator % prime != 0, "nonunit Fraction denominator")
    return value.numerator * pow(value.denominator, -1, prime) % prime


def exact_g(n: int, b: list[int], a: list[Fraction], d: list[int]) -> int:
    scaled_a = a[n] * d[n]
    check(scaled_a.denominator == 1, f"d_n a_n nonintegral at n={n}")
    return math.gcd(abs(scaled_a.numerator), d[n] * b[n])


def exact_sequence_and_gcd_gate() -> tuple[list[int], list[Fraction], list[int]]:
    max_n = 1200
    b, a = exact_apery_sequences(max_n)
    d = denominator_scales(360)

    # Independent binomial definition and variation-of-parameters formula.
    for n in range(0, 61):
        direct = sum(
            math.comb(n, k) ** 2 * math.comb(n + k, k) ** 2
            for k in range(n + 1)
        )
        check(b[n] == direct, f"binomial Apéry mismatch n={n}")
    for n in range(1, 181):
        ratio = sum(
            (Fraction(6, k**3 * b[k] * b[k - 1]) for k in range(1, n + 1)),
            Fraction(0),
        )
        check(a[n] == b[n] * ratio, f"variation formula n={n}")
        check(a[n].denominator != 0 and d[n] % a[n].denominator == 0,
              f"denominator divisibility n={n}")
        wronskian = a[n] * b[n - 1] - a[n - 1] * b[n]
        check(wronskian == Fraction(6, n**3), f"Wronskian n={n}")

    # Correct target-notation determinant divisibility.  Here d_n already is
    # lcm(1,...,n)^3, so no additional cube belongs on d_n.
    gs = [0] + [exact_g(n, b, a, d) for n in range(1, 181)]
    for n in range(1, 180):
        determinant = Fraction(6 * d[n] * d[n + 1], (n + 1) ** 3)
        check(determinant.denominator == 1, f"integer adjacent determinant n={n}")
        check(determinant.numerator % (gs[n] * gs[n + 1]) == 0,
              f"G_n G_(n+1) divisibility n={n}")
        single = Fraction(6 * d[n] * d[n], n**3)
        check(single.denominator == 1 and single.numerator % gs[n] == 0,
              f"single-index determinant divisibility n={n}")
        for p in primes_upto(n):
            if p >= 5 and p * p > n and gs[n] % p == 0:
                check(valuation(gs[n], p) <= 6, f"large-prime valuation n={n},p={p}")

    # Exact block law and its support consequence at all tested p^2>n.
    block_checks = 0
    for n in range(8, 181):
        g_n = gs[n]
        for p in primes_upto(n):
            if p < 7 or p * p <= n:
                continue
            q, r = divmod(n, p)
            lhs = fraction_mod(a[n] * p**3, p)
            rhs = fraction_mod(a[q], p) * (b[r] % p) % p
            check(lhs == rhs, f"block law n={n},p={p}")
            check((g_n % p == 0) == (rhs == 0), f"block support n={n},p={p}")
            block_checks += 1

    # The companion-height proof must include the denominator.  These values
    # catch the omitted factor in the displayed proof.tex inequality.
    alpha_log = math.log(17 + 12 * math.sqrt(2))
    height_data = []
    for q in (10, 50, 100, 200):
        check(d[q] % a[q].denominator == 0, f"companion denominator q={q}")
        height_data.append(math.log(abs(a[q].numerator)) / q)
    check(height_data[0] > alpha_log and height_data[-1] > 6.4,
          "omitted companion numerator-height factor not detected")

    print(
        "EXACT-SEQUENCES",
        f"wronskian=180 adjacent-G=179 block-checks={block_checks}",
        "numerator-log/q=" + ",".join(f"{x:.6f}" for x in height_data),
    )
    return b, a, d


def modular_orbit(p: int, full: bool = False) -> tuple[np.ndarray, np.ndarray]:
    last = p - 1 if full else p - 2
    b = np.zeros(last + 1, dtype=np.int64)
    c = np.zeros(last + 1, dtype=np.int64)
    b[0], b[1] = 1, 5 % p
    c[0], c[1] = 0, 6 % p
    for n in range(1, last):
        inv = pow(n + 1, -3, p)
        b[n + 1] = (
            (apery_coefficient(n) % p) * int(b[n])
            - n**3 * int(b[n - 1])
        ) * inv % p
        c[n + 1] = (
            (apery_coefficient(n) % p) * int(c[n])
            - n**3 * int(c[n - 1])
        ) * inv % p
    return b, c


def continuant_value(p: int, r: int, h: int) -> int:
    old, current = 0, 1
    for k in range(1, h):
        n = r + k
        old, current = current, (
            apery_coefficient(n) * current - n**6 * old
        ) % p
    return current


def chart_free_algebra_gate() -> None:
    checks = 0
    strict_union_example = None
    for p in (7, 11, 17, 73, 97, 211):
        b, c = modular_orbit(p, full=True)
        for r in range(0, p - 1):
            det = (int(b[r]) * int(c[r + 1]) - int(b[r + 1]) * int(c[r])) % p
            check(det == 6 * pow(r + 1, -3, p) % p,
                  f"Casoratian normalization p={p},r={r}")
            check(det != 0, f"zero adjacent Casoratian p={p},r={r}")
        for h in range(1, min(12, p - 2) + 1):
            for r in range(1, p - 1 - h):
                det = (int(b[r]) * int(c[r + h]) - int(b[r + h]) * int(c[r])) % p
                den = 1
                for j in range(1, h + 1):
                    den = den * (r + j) ** 3 % p
                rhs = 6 * continuant_value(p, r, h) * pow(den, -1, p) % p
                check(det == rhs, f"continuant factor 6 p={p},r={r},h={h}")
                checks += 1

        if p >= 11:
            N = p - 2
            D = 6 if p == 11 else min(18, N // 3)
            rows = []
            incidence = 0
            for gap in range(D // 2 + 1, D + 1):
                z = set()
                for r in range(1, N - gap + 1):
                    det = (int(b[r]) * int(c[r + gap]) - int(b[r + gap]) * int(c[r])) % p
                    if det == 0:
                        z.add(r)
                rows.append(z)
                incidence += len(z)
            union = len(set().union(*rows))
            check(union <= incidence, "union/incidence direction")
            if union < incidence:
                strict_union_example = (p, D, union, incidence)

    check(strict_union_example is not None, "no strict U<S live example found")
    print(
        "CHART-FREE-ALGEBRA",
        f"checks={checks} strict-U<S={strict_union_example} factor=6",
    )


def physical_pair_count(N: int, D: int) -> int:
    return sum(N - gap for gap in range(D // 2 + 1, D + 1))


def physical_histogram(p: int, D: int, rectangular: bool = False) -> tuple[np.ndarray, int]:
    b, c = modular_orbit(p)
    N = p - 2
    hist = np.zeros(p, dtype=np.int64)
    pairs = 0
    for gap in range(D // 2 + 1, D + 1):
        stop = N - D if rectangular else N - gap
        r = np.arange(1, stop + 1, dtype=np.int64)
        u = r + gap
        values = (b[r] * c[u] - b[u] * c[r]) % p
        hist += np.bincount(values, minlength=p)
        pairs += len(r)
    return hist, pairs


def determinant_full_frequency_gate() -> None:
    expected = {
        1009: 141.437464125615,
        2003: 227.222675278363,
        4003: 400.527150763946,
        8009: 676.916905792991,
        16001: 1125.078759647026,
    }
    records = []
    for p, wanted in expected.items():
        N = p - 2
        D = 2 * math.isqrt(N)
        q = D - D // 2
        hist, pairs = physical_histogram(p, D)
        check(pairs == physical_pair_count(N, D) == int(hist.sum()),
              f"physical domain p={p}")
        if D % 2 == 0:
            m = D // 2
            formula = m * N - m * (3 * m + 1) // 2
        else:
            m = D // 2
            formula = (m + 1) * N - (m + 1) * (3 * m + 2) // 2
        check(pairs == formula, f"pair formula p={p}")
        spectrum = np.fft.fft(hist)
        mean = float(np.abs(spectrum[1:]).mean())
        check(abs(mean - wanted) < 2e-8, f"all-frequency mean p={p}")
        check(np.array_equal(hist, hist[(-np.arange(p)) % p]),
              f"determinant histogram symmetry p={p}")
        check(float(np.max(np.abs(spectrum.imag))) < 1e-7 * pairs,
              f"real determinant spectrum p={p}")
        records.append((p, pairs, int(hist[0]), mean, math.sqrt(pairs)))

        if p in (1009, 2003, 4003):
            rect, rect_pairs = physical_histogram(p, D, rectangular=True)
            check(pairs - rect_pairs == q * (q - 1) // 2,
                  f"omitted boundary triangle p={p}")
            old_sample = float(np.abs(np.fft.fft(rect)[1:60]).mean())
            wanted_old = {1009: 148.364859797164, 2003: 195.104129855783,
                          4003: 462.273030476177}[p]
            check(abs(old_sample - wanted_old) < 2e-8,
                  f"old sample59 reproduction p={p}")

    for record in records:
        print(
            "DET-FULLFREQ",
            f"p={record[0]} pairs={record[1]} zeros={record[2]} ",
            f"mean={record[3]:.6f} sqrtM={record[4]:.6f}",
        )


def positive_completion_gate() -> None:
    p, D = 97, 18
    b, c = modular_orbit(p, full=True)
    r = np.arange(p, dtype=np.int64)
    phase_hist = np.zeros((p, p), dtype=np.float64)
    for gap in range(p):
        u = (r + gap) % p
        values = (b[r] * c[u] - b[u] * c[r]) % p
        phase_hist[gap] = np.bincount(values, minlength=p)
    A = np.fft.fft(phase_hist, axis=1)
    F = np.fft.fft(A, axis=0)
    neg = (-np.arange(p)) % p

    lo = D // 2 + 1
    q = D - D // 2
    check(3 * q - 2 < p, "majorant no-wrap condition")
    left = np.zeros(p)
    right = np.zeros(p)
    left[lo - q + 1 : lo + q] = 1
    right[:q] = 1
    convolution = np.zeros(p, dtype=np.int64)
    for x in np.flatnonzero(left):
        for y in np.flatnonzero(right):
            convolution[(x + y) % p] += 1
    W = convolution / q
    check(np.all(W >= 0) and np.all(W[lo : lo + q] == 1),
          "positive plateau majorant")
    check(abs(float(W.sum()) - (2 * q - 1)) < 1e-12,
          "positive majorant mass")
    What = np.fft.fft(W)
    l1 = float(np.abs(What).sum())
    cs_bound = p * math.sqrt((2 * q - 1) / q)
    check(l1 <= cs_bound + 1e-9, "positive majorant Fourier L1")
    direct = W @ A
    completed = (What[neg, None] * F).sum(axis=0) / p
    error = float(np.max(np.abs(direct - completed)))
    check(error < 1e-7, "completion normalization")

    # The exact clock-twisted Weyl identity, independently at a small prime.
    ps, t = 17, 5
    bs, cs = modular_orbit(ps, full=True)
    zeta = cmath.exp(2j * math.pi / ps)
    At = []
    for gap in range(ps):
        At.append(sum(
            zeta ** ((t * (int(bs[x]) * int(cs[(x + gap) % ps])
                           - int(bs[(x + gap) % ps]) * int(cs[x]))) % ps)
            for x in range(ps)
        ))
    Ft = [sum(At[gap] * zeta ** ((xi * gap) % ps) for gap in range(ps))
          for xi in range(ps)]
    for xi in range(ps):
        rhs = 0j
        for x in range(ps):
            alpha, beta = (-t * int(cs[x])) % ps, (t * int(bs[x])) % ps
            twisted = sum(
                zeta ** ((alpha * int(bs[u]) + beta * int(cs[u]) + xi * u) % ps)
                for u in range(ps)
            )
            rhs += zeta ** ((-xi * x) % ps) * twisted
        check(abs(Ft[xi] - rhs) < 1e-8, f"twisted Weyl xi={xi}")

    print(
        "POSITIVE-COMPLETION",
        f"q={q} mass={W.sum():.0f} Fourier-L1={l1:.6f} ",
        f"CS={cs_bound:.6f} inversion-error={error:.2e}",
    )


def cyclic_convolution(a: list[int]) -> list[int]:
    p = len(a)
    return [sum(a[x] * a[(s - x) % p] for x in range(p)) for s in range(p)]


def determinant_moment_gate() -> None:
    records = []
    for p, D in ((31, 8), (43, 10), (61, 12), (101, 18)):
        hist, M = physical_histogram(p, D)
        nu = [int(x) for x in hist]
        e2 = sum(x * x for x in nu)
        conv = cyclic_convolution(nu)
        e4 = sum(x * x for x in conv)
        spectrum = np.fft.fft(hist)
        lhs2 = float(np.sum(np.abs(spectrum[1:]) ** 2))
        lhs4 = float(np.sum(np.abs(spectrum[1:]) ** 4))
        rhs2 = p * e2 - M * M
        rhs4 = p * e4 - M**4
        check(abs(lhs2 - rhs2) < 2e-8 * max(1, rhs2), f"second moment p={p}")
        check(abs(lhs4 - rhs4) < 3e-8 * max(1, rhs4), f"fourth moment p={p}")

        # Exact row/cross-gap centered-variance decomposition.
        b, c = modular_orbit(p)
        N = p - 2
        rows = []
        for gap in range(D // 2 + 1, D + 1):
            row = np.zeros(p, dtype=np.int64)
            rr = np.arange(1, N - gap + 1, dtype=np.int64)
            vals = (b[rr] * c[rr + gap] - b[rr + gap] * c[rr]) % p
            row += np.bincount(vals, minlength=p)
            rows.append(row)
        row_v = 0
        cross_v = 0
        for i, row in enumerate(rows):
            mr = int(row.sum())
            row_v += p * int(row @ row) - mr * mr
            for other in rows[i + 1 :]:
                mo = int(other.sum())
                cross_v += 2 * (p * int(row @ other) - mr * mo)
        check(rhs2 == row_v + cross_v, f"row covariance decomposition p={p}")
        records.append((p, M, rhs2, rhs4))
    print("DET-MOMENTS", records)


def primitive_root(p: int) -> int:
    factors = []
    n = p - 1
    q = 2
    while q * q <= n:
        if n % q == 0:
            factors.append(q)
            while n % q == 0:
                n //= q
        q += 1
    if n > 1:
        factors.append(n)
    for g in range(2, p):
        if all(pow(g, (p - 1) // q, p) != 1 for q in factors):
            return g
    raise GateFailure(f"no primitive root p={p}")


def additive_energy(points: list[tuple[int, int]], p: int) -> int:
    sums = Counter(
        ((x1 + x2) % p, (y1 + y2) % p)
        for x1, y1 in points
        for x2, y2 in points
    )
    return sum(v * v for v in sums.values())


def hyperbola_clock_gate() -> None:
    records = []
    for p, D in ((101, 12), (211, 18), (401, 28)):
        g = primitive_root(p)
        m = (p - 1) // 2
        vectors = [(0, 0)] * (p - 1)
        support = []
        for j in range(1, m + 1):
            v = (pow(g, j, p), pow(g, -j, p))
            support.append(v)
            vectors[j] = v
            vectors[p - 1 - j] = v
        mult = Counter(vectors[1:])
        check(list(mult.values()).count(1) == 1 and all(v in (1, 2) for v in mult.values()),
              f"hyperbola reflection multiplicity p={p}")
        energy = additive_energy(support, p)
        check(energy <= 3 * len(support) ** 2, f"hyperbola additive energy p={p}")
        if p == 101:
            max_line = 0
            for aa in range(p):
                for bb in range(p):
                    if aa == bb == 0:
                        continue
                    counts = Counter((aa * x + bb * y) % p for x, y in support)
                    max_line = max(max_line, max(counts.values()))
            check(max_line <= 2, "hyperbola line intersection")

        N = p - 2
        hist = np.zeros(p, dtype=np.int64)
        for gap in range(D // 2 + 1, D + 1):
            expected = (pow(g, -gap, p) - pow(g, gap, p)) % p
            for r in range(1, m - gap + 1):
                x, y = vectors[r]
                z, w = vectors[r + gap]
                check((x * w - z * y) % p == expected,
                      f"coherent hyperbola layer p={p},gap={gap}")
            rr = range(1, N - gap + 1)
            for r in rr:
                x, y = vectors[r]
                z, w = vectors[r + gap]
                hist[(x * w - z * y) % p] += 1
        M = int(hist.sum())
        spectrum = np.fft.fft(hist)
        mean = float(np.abs(spectrum[1:]).mean())
        e2 = int(hist @ hist)
        second_mean = (p * e2 - M * M) / (p - 1)
        certified_lower = second_mean / M
        check(mean + 1e-8 >= certified_lower and certified_lower > 0.1 * (p - 1),
              f"linear L1 countermodel p={p}")
        records.append((p, len(support), energy, mean, certified_lower))
    print("HYPERBOLA-CLOCK", records)


def vector_parseval_and_full_spectrum_gate() -> None:
    extra = []
    primes = [p for p in primes_upto(1000) if p >= 5]
    for p in primes:
        b, c = modular_orbit(p)
        N = p - 2
        mult = Counter((int(b[r]), int(c[r])) for r in range(1, N + 1))
        square_mass = sum(v * v for v in mult.values())
        forced = 2 * N - 1
        check(square_mass >= forced, f"reflection Parseval floor p={p}")
        if square_mass != forced:
            extra.append((p, square_mass - forced, max(mult.values())))
    check(len(extra) == 17, "extra vector-collision prime count through 1000")
    check(any(p == 73 and mx == 4 for p, _, mx in extra), "p=73 vector collision")
    check(any(p == 997 and mx == 4 for p, _, mx in extra), "p=997 vector collision")

    expected = {
        1009: (2013, 176.033819002445),
        2003: (4001, 244.481930018687),
        3001: (6005, 294.307009318177),
    }
    records = []
    for p, (wanted_mass, wanted_max) in expected.items():
        b, c = modular_orbit(p)
        N = p - 2
        H = np.zeros((p, p), dtype=np.float32)
        np.add.at(H, (b[1:], c[1:]), 1)
        square_mass = int(np.sum(H.astype(np.int64) ** 2))
        check(square_mass == wanted_mass, f"vector multiplicity mass p={p}")
        S = np.fft.fft2(H)
        S[0, 0] = 0
        maximum = float(np.abs(S).max())
        check(abs(maximum - wanted_max) < 2e-7, f"vector full spectrum p={p}")
        records.append((p, square_mass, maximum, maximum / math.sqrt(N)))
        del H, S
    print(f"VECTOR-PARSEVAL extra-primes={len(extra)} examples={extra[:5]}")
    print("VECTOR-FULLFREQ", records)


def poly_trim(poly: list[int], p: int) -> list[int]:
    answer = [x % p for x in poly]
    while len(answer) > 1 and answer[-1] == 0:
        answer.pop()
    return answer


def poly_add(a: list[int], b: list[int], p: int, scale: int = 1) -> list[int]:
    out = [0] * max(len(a), len(b))
    for i, x in enumerate(a):
        out[i] += x
    for i, x in enumerate(b):
        out[i] += scale * x
    return poly_trim(out, p)


def poly_mul(a: list[int], b: list[int], p: int) -> list[int]:
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] = (out[i + j] + x * y) % p
    return poly_trim(out, p)


def shifted_power(shift: int, exponent: int, p: int) -> list[int]:
    out = [1]
    for _ in range(exponent):
        out = poly_mul(out, [shift, 1], p)
    return out


def shifted_P(shift: int, p: int) -> list[int]:
    # P(X+shift), independently expanded by polynomial arithmetic.
    z = [shift, 1]
    return poly_add(
        poly_add(poly_add([5], poly_mul([27], z, p), p),
                 poly_mul([51], poly_mul(z, z, p), p), p),
        poly_mul([34], poly_mul(poly_mul(z, z, p), z, p), p),
    )


def modular_continuants(p: int, height: int) -> list[list[int]]:
    values = [[0], [1]]
    for h in range(1, height):
        values.append(poly_add(
            poly_mul(shifted_P(h, p), values[h], p),
            poly_mul(shifted_power(h, 6, p), values[h - 1], p),
            p,
            scale=-1,
        ))
    return values


def poly_divrem(a: list[int], b: list[int], p: int) -> list[int]:
    a = poly_trim(a, p)
    b = poly_trim(b, p)
    check(b != [0], "polynomial division by zero")
    inv = pow(b[-1], -1, p)
    while len(a) >= len(b) and a != [0]:
        coefficient = a[-1] * inv % p
        shift = len(a) - len(b)
        for i, x in enumerate(b):
            a[i + shift] = (a[i + shift] - coefficient * x) % p
        a = poly_trim(a, p)
    return a


def poly_gcd(a: list[int], b: list[int], p: int) -> list[int]:
    a, b = poly_trim(a, p), poly_trim(b, p)
    while b != [0]:
        a, b = b, poly_divrem(a, b, p)
    inv = pow(a[-1], -1, p)
    return [(x * inv) % p for x in a]


def poly_derivative(a: list[int], p: int) -> list[int]:
    return poly_trim([i * a[i] for i in range(1, len(a))] or [0], p)


def poly_evaluate(a: list[int], x: int, p: int) -> int:
    answer = 0
    for coefficient in reversed(a):
        answer = (answer * x + coefficient) % p
    return answer


def saturation_counterexample_gate() -> None:
    c73 = modular_continuants(73, 4)
    check(poly_evaluate(c73[3], -3, 73) == 0, "N3(-3) mod 73")
    check(poly_evaluate(c73[4], -3, 73) == 0, "N4(-3) mod 73")
    check(poly_gcd(c73[3], c73[4], 73) == [3, 1], "gcd N3,N4 mod73")
    # Integer evaluation follows directly from the same recurrence.
    old, current = 0, 1
    x = -3
    for h in range(1, 3):
        old, current = current, apery_coefficient(x + h) * current - (x + h) ** 6 * old
    check(current == 584 == 8 * 73, "integer N3(-3)")

    c211 = modular_continuants(211, 32)
    repeated = poly_gcd(c211[32], poly_derivative(c211[32], 211), 211)
    check(repeated == [114, 33, 1], "N32 repeated part mod211")
    check(poly_evaluate(repeated, 89, 211) == 0, "physical repeated root 89")
    print("SATURATION p=73 gcd=X+3; p=211 repeated=[114,33,1]")


def original_master_hits(n: int, b: list[int], primes: list[int]) -> list[tuple[int, int]]:
    hits = []
    for r in range((n + 1) // 2):  # exactly 2r<n
        common = math.gcd(b[r], n - r)
        for p in primes:
            if p > common:
                break
            if common % p == 0 and p * p > n and p > r:
                hits.append((p, r))
    return sorted(hits)


def remainder_master_hits(n: int, b: list[int], primes: list[int]) -> list[tuple[int, int]]:
    return sorted(
        (p, n % p)
        for p in primes
        if p <= n and p * p > n and b[n % p] % p == 0
    )


def master_digit_and_average_gate(b: list[int], a: list[Fraction], d: list[int]) -> None:
    primes = primes_upto(1200)
    mapping_checks = 0
    for n in range(8, 361):
        direct = original_master_hits(n, b, primes)
        remainder = remainder_master_hits(n, b, primes)
        check(direct == remainder, f"master expansion n={n}")
        for p, r in remainder:
            check(p > r and 2 * r < n, f"automatic range n={n},p={p}")
        mapping_checks += 1

    # The ledger's R_p incidence has extra representatives and is not M(n).
    n, p = 101, 31
    R = [r for r in range((n + 1) // 2) if r % p == n % p and b[r] % p == 0]
    check(R == [8, 39], "R_p counterexample representatives")
    check((p, 8) in original_master_hits(n, b, primes)
          and (p, 39) not in original_master_hits(n, b, primes),
          "p>r removes the second incidence")

    # q>=2 is a genuine part of the exact master sum and of G_n.
    g37 = exact_g(37, b, a, d)
    check(divmod(37, 17) == (2, 3) and b[3] % 17 == 0,
          "middle quotient witness data")
    check(g37 % 17 == 0, "middle quotient witness in G_37")
    check(not (2 * 17 > 37), "middle witness accidentally in top window")

    # Gessel/Lucas digit criterion, checked against independently generated b_n.
    digit_checks = 0
    for p in [q for q in primes_upto(97) if q >= 2]:
        zero = {r for r in range(p) if b[r] % p == 0}
        for n in range(0, 601):
            digits = []
            x = n
            if x == 0:
                digits = [0]
            while x:
                digits.append(x % p)
                x //= p
            criterion = any(digit in zero for digit in digits)
            check((b[n] % p == 0) == criterion,
                  f"digit criterion p={p},n={n}")
            digit_checks += 1
    top_checks = 0
    for n in range(40, 140):
        for p in primes:
            if p > n or 2 * p <= n or p <= 5:
                continue
            check(b[n] % p == (5 * b[n - p]) % p,
                  f"top Gessel congruence n={n},p={p}")
            check((b[n] % p == 0) == (b[n - p] % p == 0),
                  f"top equivalence n={n},p={p}")
            top_checks += 1

    # Correct truncated averaging identity and its S(X) sandwich (p>5 exact).
    zero_sets = {p: {r for r in range(p) if b[r] % p == 0} for p in primes}

    def A_count(x: int, exclude_small: bool = False) -> int:
        return sum(
            1
            for p in primes
            if p <= x and (not exclude_small or p > 5)
            for r in zero_sets[p]
            if 1 <= r <= x - p
        )

    def S_count(x: int, exclude_small: bool = False) -> int:
        return sum(
            len(zero_sets[p])
            for p in primes
            if p <= x and (not exclude_small or p > 5)
        )

    def actual_top_sum(x: int, exclude_small: bool = False) -> int:
        return sum(
            1
            for n in range(2, x + 1)
            for p in primes
            if p <= n and 2 * p > n and (not exclude_small or p > 5)
            and b[n] % p == 0
        )

    for x in (50, 100, 200, 400, 600):
        check(actual_top_sum(x, True) == A_count(x, True),
              f"corrected averaging identity X={x}")
        check(A_count(x) <= S_count(x) <= A_count(2 * x),
              f"AVG sandwich X={x}")
    check((actual_top_sum(600), actual_top_sum(600, True), A_count(600),
           A_count(600, True), S_count(600), A_count(1200))
          == (85, 80, 82, 80, 109, 137), "X=600 averaging data")

    # A finite reflection-symmetric row model: all selected q=2 primes align
    # at N, while every induced q=1 target has load at most one.
    N = 10000
    selected = [p for p in primes_upto(N) if 5 * p > 2 * N and 20 * p <= 9 * N]
    top_load: Counter[int] = Counter()
    for p in selected:
        q, r = divmod(N, p)
        reflected = p - 1 - r
        row = {r, reflected}
        check(q == 2 and N % p in row and len(row) == 2,
              f"q=2 alignment row p={p}")
        for z in row:
            top_load[p + z] += 1
    check(len(selected) == 60 and max(top_load.values()) == 1,
          "q=2/q=1 separation model")

    print(
        "MASTER-DIGIT-AVG",
        f"mapping={mapping_checks} digit={digit_checks} top={top_checks} ",
        "R31(101)=[8,39] middle=(n,p,q,r)=(37,17,2,3) ",
        "X600=(TOPall85,TOP>5=80,A=82,S=109) q2-model=60/topload1",
    )


def exponent_budget_gate() -> None:
    from fractions import Fraction as Q

    # sqrt(pD) / D^(2-eta), with D=N^(1/2)L.
    eta = Q(2, 5)
    n_exp = -Q(1, 4) + eta / 2
    l_exp = -Q(3, 2) + eta
    check(n_exp < 0 and l_exp < 0, "square-root exponent margin eta=2/5")
    check(-Q(1, 4) + Q(1, 2) / 2 == 0,
          "eta=1/2 square-root boundary")
    # p^(1-kappa) / D^(2-eta).
    kappa = Q(1, 10)
    eta2 = Q(3, 20)
    check(eta2 < 2 * kappa and eta2 < 1,
          "power-saving determinant implication")
    check(eta2 / 2 - kappa < 0, "power-saving N exponent")
    # A p (or p log p) completion error has positive N exponent for eta>0.
    check(eta2 / 2 > 0, "completion obstruction exponent")
    # HM_k closes the worst-case lambda_X budget exactly when k>6.
    check(Q(2, 3) + Q(2, 7) < 1 and Q(2, 3) + Q(2, 6) == 1,
          "high-moment threshold k>6")
    print(
        "EXPONENT-BUDGET",
        "sqrt-ratio=N^(-1/4+eta/2)L^(-3/2+eta); ",
        "power-ratio=N^(eta/2-kappa)L^(eta-2); HM threshold k>6",
    )


def main() -> None:
    b, a, d = exact_sequence_and_gcd_gate()
    chart_free_algebra_gate()
    saturation_counterexample_gate()
    determinant_full_frequency_gate()
    positive_completion_gate()
    determinant_moment_gate()
    vector_parseval_and_full_spectrum_gate()
    hyperbola_clock_gate()
    master_digit_and_average_gate(b, a, d)
    exponent_budget_gate()


if __name__ == "__main__":
    try:
        main()
    except Exception:
        traceback.print_exc()
        print("FAIL")
        sys.exit(1)
    print("PASS")
