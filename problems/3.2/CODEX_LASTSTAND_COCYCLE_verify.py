#!/usr/bin/env python3
"""Exact audit for CODEX_LASTSTAND_COCYCLE_report.md.

The script checks the corrected frame/Borel cocycle, the failure of the
projective action asserted in CODEX_SPEC_laststand_cocycle.md, the exact R1
dispersion sums and private-pair obstruction, the R2 renewal/correlation and
fixed-prime caveats, the R3 no-adiabatic identities, and the displayed Apéry
shell census.  All arithmetic gates are deterministic.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from decimal import Decimal, ROUND_FLOOR, localcontext
from fractions import Fraction
import itertools
import math
import sys

import sympy as sp


def apery_P(n: int | sp.Expr) -> int | sp.Expr:
    return 34 * n**3 + 51 * n**2 + 27 * n + 5


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    return all(n % q for q in range(3, math.isqrt(n) + 1, 2))


def matmul(A: list[list[int]], B: list[list[int]], p: int) -> list[list[int]]:
    return [
        [sum(A[i][k] * B[k][j] for k in range(2)) % p for j in range(2)]
        for i in range(2)
    ]


def companion_matrix(n: int, p: int) -> list[list[int]]:
    inverse_cube = pow(n + 1, -3, p)
    alpha = apery_P(n) * inverse_cube % p
    beta = n**3 * inverse_cube % p
    return [[0, 1], [-beta % p, alpha]]


def apery_orbit(p: int) -> tuple[list[int], list[int]]:
    """Return b,c on 0,...,p-2, with c_0=0 and c_1=6."""

    N = p - 2
    b = [0] * (N + 1)
    c = [0] * (N + 1)
    b[0], b[1] = 1, 5 % p
    c[0], c[1] = 0, 6 % p
    for n in range(1, N):
        inverse_cube = pow(n + 1, -3, p)
        b[n + 1] = (
            (apery_P(n) * b[n] - n**3 * b[n - 1]) * inverse_cube
        ) % p
        c[n + 1] = (
            (apery_P(n) * c[n] - n**3 * c[n - 1]) * inverse_cube
        ) % p
    return b, c


def projective_key(x: int, y: int, p: int) -> tuple[int, int]:
    x %= p
    y %= p
    if y:
        return (x * pow(y, -1, p) % p, 1)
    assert x
    return (1, 0)


def transfer_matrix(r: int, d: int, p: int) -> list[list[int]]:
    T = [[1, 0], [0, 1]]
    for s in range(r, r + d):
        T = matmul(companion_matrix(s, p), T, p)
    return T


def foundational_cocycle_gate() -> None:
    """Refute the false action and certify the corrected frame criterion."""

    alpha_1 = Fraction(117, 8)
    beta_1 = Fraction(1, 8)
    false_image = (Fraction(6), -beta_1 * 5 + alpha_1 * 6)
    true_next = (Fraction(73), Fraction(351, 4))
    false_determinant = false_image[0] * true_next[1] - false_image[1] * true_next[0]
    assert false_image == (Fraction(6), Fraction(697, 8))
    assert false_determinant == Fraction(-46669, 8)
    assert false_determinant

    # Exhaust the corrected identity on a complete small physical block.
    p = 101
    b, c = apery_orbit(p)
    N = p - 2
    for r in range(1, N + 1):
        determinant = (b[r - 1] * c[r] - c[r - 1] * b[r]) % p
        assert determinant == 6 * pow(r, -3, p) % p
        if r < N:
            G_r = [[b[r - 1], c[r - 1]], [b[r], c[r]]]
            G_next = [[b[r], c[r]], [b[r + 1], c[r + 1]]]
            assert matmul(companion_matrix(r, p), G_r, p) == G_next

        T = [[1, 0], [0, 1]]
        N_previous, N_current = 0, 1
        denominator = 1
        for d in range(1, N - r + 1):
            T = matmul(companion_matrix(r + d - 1, p), T, p)
            denominator = denominator * pow(r + d, 3, p) % p
            expected_21 = -r**3 * N_current * pow(denominator, -1, p) % p
            assert T[1][0] == expected_21
            collision = projective_key(b[r], c[r], p) == projective_key(
                b[r + d], c[r + d], p
            )
            assert collision == (T[1][0] == 0)
            N_next = (
                apery_P(r + d) * N_current
                - (r + d) ** 6 * N_previous
            ) % p
            N_previous, N_current = N_current, N_next

    # Live two-run witnesses: both the claimed fixed-point test and v(r,d)
    # fail, while the corrected lower-left-entry test succeeds.
    p, d = 997, 182
    b, c = apery_orbit(p)
    witnesses: list[tuple[object, ...]] = []
    for r in (248, 565):
        assert projective_key(b[r], c[r], p) == projective_key(
            b[r + d], c[r + d], p
        )
        assert projective_key(b[r + 1], c[r + 1], p) == projective_key(
            b[r + d + 1], c[r + d + 1], p
        )
        T = transfer_matrix(r, d, p)
        xi = projective_key(b[r], c[r], p)
        image = projective_key(
            T[0][0] * b[r] + T[0][1] * c[r],
            T[1][0] * b[r] + T[1][1] * c[r],
            p,
        )
        alpha_difference = (
            apery_P(r) * pow(r + 1, -3, p)
            - apery_P(r + d) * pow(r + d + 1, -3, p)
        ) % p
        beta_difference = (
            r**3 * pow(r + 1, -3, p)
            - (r + d) ** 3 * pow(r + d + 1, -3, p)
        ) % p
        proposed_v = projective_key(alpha_difference, beta_difference, p)
        assert T[1][0] == 0
        assert image != xi
        assert proposed_v != xi
        witnesses.append((r, xi, image, proposed_v, T))

    assert witnesses == [
        (248, (758, 1), (344, 1), (798, 1), [[929, 174], [0, 297]]),
        (565, (409, 1), (780, 1), (165, 1), [[355, 66], [0, 323]]),
    ]
    print("FOUNDATIONAL COCYCLE CORRECTION: PASS")
    print("  rational false-action determinant = -46669/8")
    for witness in witnesses:
        print("  p=997 (r,xi,T*xi,v,T) =", witness)


def r1_dispersion_gate() -> None:
    """Check every finite-sum identity used in the R1 calculation."""

    for m in range(1, 80):
        full_sum = sum((m - g) * (g - 1) for g in range(1, m))
        assert full_sum == math.comb(m, 3)
        for G in range(0, m):
            truncated = 3 * sum(
                (m - g) * (g - 1) for g in range(1, G + 1)
            )
            closed = G * (G - 1) * (3 * m - 2 * G - 2) // 2
            assert truncated == closed
        for H in range(1, m + 1):
            sliding = sum(
                (H - g) * (m - g) * (g - 1) for g in range(1, H)
            )
            closed = H * (H - 1) * (H - 2) * (2 * m - H - 1) // 12
            assert sliding == closed

    fibre = (0, 6, 10)
    gaps = {
        x: {y - x for y in fibre if y > x}
        for x in fibre
    }
    assert gaps[0] == {6, 10}
    assert 4 not in gaps[0]
    assert 4 in gaps[6]

    # A literal private-pair incidence/word countermodel to the R1 inputs.
    N = 100_000
    D = math.ceil(math.sqrt(N) * math.log(N))
    t = N // (2 * D)
    assert D * D >= N and 1 <= t <= D // 2
    word = list(range(N))
    cursor = 0
    colour = N
    prescribed: dict[int, set[int]] = {}
    for j in range(t):
        d = D - j
        assert D / 2 < d <= D
        prescribed[d] = set(range(cursor, cursor + d))
        for i in range(d):
            word[cursor + i] = colour
            word[cursor + d + i] = colour
            colour += 1
        cursor += 2 * d
    assert cursor <= N

    fibres: dict[int, list[int]] = defaultdict(list)
    for position, value in enumerate(word):
        fibres[value].append(position)
    observed: dict[int, set[int]] = defaultdict(set)
    for positions in fibres.values():
        if len(positions) == 2:
            r, s = positions
            observed[s - r].add(r)
        else:
            assert len(positions) == 1
    assert dict(observed) == prescribed
    assert all(len(observed[d]) == d <= 3 * (d - 1) for d in observed)
    all_bases = set().union(*observed.values())
    assert sum(map(len, observed.values())) == len(all_bases)
    assert all(observed[d].isdisjoint(observed[e]) for d in observed for e in observed if d < e)
    S = sum(observed[d].__len__() for d in observed)
    assert S == t * D - t * (t - 1) // 2
    assert (N, D, t, S) == (100000, 3641, 13, 47255)
    print("R1 EXACT DISPERSION AND PRIVATE-PAIR OBSTRUCTION: PASS")
    print("  countermodel (N,D,t,S=U) =", (N, D, t, S))


def build_integer_continuants(X: sp.Symbol, height: int) -> list[sp.Expr]:
    values: list[sp.Expr] = [sp.Integer(0), sp.Integer(1)]
    for h in range(1, height):
        values.append(
            sp.expand(
                apery_P(X + h) * values[-1]
                - (X + h) ** 6 * values[-2]
            )
        )
    return values


def poly_trim(a: list[int], p: int) -> list[int]:
    a = [x % p for x in a]
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a


def poly_add(a: list[int], b: list[int], p: int, scale: int = 1) -> list[int]:
    out = [0] * max(len(a), len(b))
    for i, value in enumerate(a):
        out[i] += value
    for i, value in enumerate(b):
        out[i] += scale * value
    return poly_trim(out, p)


def poly_mul(a: list[int], b: list[int], p: int) -> list[int]:
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] = (out[i + j] + x * y) % p
    return poly_trim(out, p)


def shifted_power(h: int, exponent: int, p: int) -> list[int]:
    return [
        math.comb(exponent, i) * pow(h, exponent - i, p) % p
        for i in range(exponent + 1)
    ]


def shifted_P(h: int, p: int) -> list[int]:
    return [
        apery_P(h) % p,
        (102 * h * h + 102 * h + 27) % p,
        (102 * h + 51) % p,
        34 % p,
    ]


def modular_continuants(p: int, height: int) -> list[list[int]]:
    values = [[0], [1]]
    for h in range(1, height):
        values.append(
            poly_add(
                poly_mul(shifted_P(h, p), values[-1], p),
                poly_mul(shifted_power(h, 6, p), values[-2], p),
                p,
                scale=-1,
            )
        )
    return values


def poly_divrem(a: list[int], b: list[int], p: int) -> list[int]:
    a = poly_trim(a[:], p)
    b = poly_trim(b[:], p)
    inverse_lead = pow(b[-1], -1, p)
    while a != [0] and len(a) >= len(b):
        multiplier = a[-1] * inverse_lead % p
        offset = len(a) - len(b)
        for i, value in enumerate(b):
            a[offset + i] = (a[offset + i] - multiplier * value) % p
        a = poly_trim(a, p)
    return a


def poly_gcd(a: list[int], b: list[int], p: int) -> list[int]:
    while b != [0]:
        a, b = b, poly_divrem(a, b, p)
    inverse_lead = pow(a[-1], -1, p)
    return poly_trim([value * inverse_lead for value in a], p)


def poly_derivative(a: list[int], p: int) -> list[int]:
    return poly_trim([i * a[i] for i in range(1, len(a))] or [0], p)


def poly_evaluate(a: list[int], x: int, p: int) -> int:
    result = 0
    for coefficient in reversed(a):
        result = (result * x + coefficient) % p
    return result


def r2_correlation_gate() -> None:
    X = sp.symbols("X")
    continuants = build_integer_continuants(X, 10)
    for a in range(1, 6):
        for g in range(1, 6):
            lhs = continuants[a + g]
            rhs = (
                continuants[g].subs(X, X + a) * continuants[a + 1]
                - (X + a + 1) ** 6
                * continuants[g - 1].subs(X, X + a + 1)
                * continuants[a]
            )
            assert sp.expand(lhs - rhs) == 0

    # Separate the two triangle slices which the report distinguishes.
    # In the abstract fibre {0,6,13}, both consecutive gaps are in
    # B=(5,10], yet both shell bases are singletons and I_B=0.
    abstract_fibre = (0, 6, 13)
    abstract_D = 10
    abstract_B = range(abstract_D // 2 + 1, abstract_D + 1)
    abstract_Z: dict[int, set[int]] = defaultdict(set)
    for i, r in enumerate(abstract_fibre):
        for endpoint in abstract_fibre[i + 1 :]:
            abstract_Z[endpoint - r].add(r)
    abstract_multiplicity = Counter()
    for d in abstract_B:
        for r in abstract_Z[d]:
            abstract_multiplicity[r] += 1
    abstract_I = sum(
        math.comb(value, 2) for value in abstract_multiplicity.values()
    )
    abstract_kappa = sum(
        r in abstract_Z[6] and r + 6 in abstract_Z[7]
        for r in abstract_fibre
    )
    assert abstract_multiplicity == {0: 1, 6: 1}
    assert abstract_I == 0 and abstract_kappa == 1

    # On a live small orbit, check both the same-base I_B slice and the
    # consecutive-edge kappa slice against the full triangle census.
    p, D = 101, 20
    b, c = apery_orbit(p)
    N = p - 2
    fibres: dict[tuple[int, int], list[int]] = defaultdict(list)
    for r in range(1, N + 1):
        fibres[projective_key(b[r], c[r], p)].append(r)
    Z: dict[int, set[int]] = defaultdict(set)
    Q_2D = 0
    for positions in fibres.values():
        for i, r in enumerate(positions):
            for endpoint in positions[i + 1 :]:
                if endpoint - r <= 2 * D:
                    Z[endpoint - r].add(r)
        Q_2D += sum(
            endpoints[-1] - endpoints[0] <= 2 * D
            for endpoints in itertools.combinations(positions, 3)
        )
    B_shell = range(D // 2 + 1, D + 1)
    multiplicities = Counter()
    for d in B_shell:
        for r in Z[d]:
            multiplicities[r] += 1
    I_B = sum(math.comb(value, 2) for value in multiplicities.values())
    same_base_correlations = sum(
        len(Z[h] & Z[k])
        for h in B_shell
        for k in B_shell
        if h < k
    )
    assert I_B == same_base_correlations
    kappa_sum = 0
    for a in B_shell:
        for g in B_shell:
            shifted = {r for r in Z[a] if r + a in Z[g]}
            assert shifted == Z[a] & Z[a + g]
            kappa_sum += len(shifted)
    assert kappa_sum <= Q_2D

    # Nonzero over Z/Q is not nonzero modulo the working prime.
    exact_pole_value = apery_P(-1) * apery_P(-2) - 1
    assert exact_pole_value == 584 == 8 * 73
    continuants_73 = modular_continuants(73, 4)
    assert poly_evaluate(continuants_73[3], -3, 73) == 0
    assert poly_evaluate(continuants_73[4], -3, 73) == 0
    assert poly_gcd(continuants_73[3], continuants_73[4], 73) == [3, 1]

    # Squarefree over Q also need not remain squarefree modulo the live p.
    continuants_211 = modular_continuants(211, 32)
    N_32 = continuants_211[32]
    gcd_32 = poly_gcd(N_32, poly_derivative(N_32, 211), 211)
    assert gcd_32 == [114, 33, 1]  # (X-89)^2 in F_211[X]
    assert poly_evaluate(N_32, 89, 211) == 0
    assert poly_evaluate(poly_derivative(N_32, 211), 89, 211) == 0

    f = sp.Poly(X + 73, X, domain=sp.ZZ)
    g = sp.Poly(X + 146, X, domain=sp.ZZ)
    assert sp.gcd(f, g).degree() == 0
    assert abs(int(sp.resultant(f.as_expr(), g.as_expr(), X))) == 73
    assert sp.gcd(
        sp.Poly(f.as_expr(), X, modulus=73),
        sp.Poly(g.as_expr(), X, modulus=73),
    ).degree() == 1

    # The exact low-multiplicity decomposition used for the surviving target.
    for multiplicity in range(0, 100):
        for K in range(2, 30):
            low = multiplicity if 1 <= multiplicity < K else 0
            high = 2 * math.comb(multiplicity, 2) // (K - 1)
            assert multiplicity <= low + high
    assert 66 * 4 == 264
    assert 264 * 4 == 1056
    assert 2 * 1056 == 2112
    print("R2 RENEWAL, SPECIALIZATION, AND LOW-MULTIPLICITY GATES: PASS")
    print("  p=73: N_3(-3)=584=8*73 and gcd(N_3,N_4)=X+3 mod 73")
    print("  p=211: gcd(N_32,N_32')=(X-89)^2 mod 211")


def r3_no_adiabatic_gate() -> None:
    s, z = sp.symbols("s z")
    alpha = apery_P(s) / (s + 1) ** 3
    beta = s**3 / (s + 1) ** 3
    alpha_next = alpha.subs(s, s + 1)
    beta_next = beta.subs(s, s + 1)
    denominator = (s + 1) ** 3 * (s + 2) ** 3
    A = 51 * s**4 + 252 * s**3 + 435 * s**2 + 306 * s + 77
    B = 3 * s**4 + 12 * s**3 + 15 * s**2 + 6 * s + 1
    assert sp.cancel(alpha_next - alpha - A / denominator) == 0
    assert sp.cancel(beta_next - beta - B / denominator) == 0
    resultant = int(sp.resultant(A, B, s))
    assert resultant == 248832
    assert sp.factorint(resultant) == {2: 10, 3: 5}

    discriminant_numerator = (
        (3 * s**2 + 3 * s + 1)
        * (384 * s**4 + 768 * s**3 + 579 * s**2 + 195 * s + 25)
    )
    assert sp.cancel(alpha**2 - 4 * beta - discriminant_numerator / (s + 1) ** 6) == 0

    R8 = (
        108 * s**8
        + 864 * s**7
        + 2865 * s**6
        + 5094 * s**5
        + 5253 * s**4
        + 3180 * s**3
        + 1095 * s**2
        + 198 * s
        + 16
    )
    chi = z**2 - alpha * z + beta
    chi_next = z**2 - alpha_next * z + beta_next
    symbolic_eigenvector = sp.Matrix([[0, 1], [-beta, alpha]]) * sp.Matrix([1, z])
    assert sp.simplify(symbolic_eigenvector[0] - z) == 0
    assert sp.simplify(symbolic_eigenvector[1] - z**2) == -chi
    assert sp.cancel(
        sp.resultant(chi, chi_next, z)
        + 24 * R8 / ((s + 1) ** 6 * (s + 2) ** 6)
    ) == 0

    a0, a1, b0, b1 = sp.symbols("a0 a1 b0 b1")
    M0 = sp.Matrix([[0, 1], [-b0, a0]])
    M1 = sp.Matrix([[0, 1], [-b1, a1]])
    expected_commutator = sp.Matrix(
        [[b0 - b1, a1 - a0], [a1 * b0 - a0 * b1, b1 - b0]]
    )
    assert M0 * M1 - M1 * M0 == expected_commutator

    A_poly = sp.Poly(A, s)
    B_poly = sp.Poly(B, s)
    disc_poly = sp.Poly(discriminant_numerator, s)
    R8_poly = sp.Poly(R8, s)
    for p in (5, 7, 11, 13, 17, 19, 23, 101, 997):
        regular = [x for x in range(p) if x not in (p - 1, p - 2)]
        assert all(
            not (
                int(A_poly.eval(x)) % p == 0
                and int(B_poly.eval(x)) % p == 0
            )
            for x in regular
        )
        assert sum(int(disc_poly.eval(x)) % p == 0 for x in range(p)) <= 6
        assert sum(int(R8_poly.eval(x)) % p == 0 for x in range(p)) <= 8
        exceptions = {
            x
            for x in regular
            if int(disc_poly.eval(x)) % p == 0
            or int(disc_poly.eval((x + 1) % p)) % p == 0
            or int(R8_poly.eval(x)) % p == 0
        }
        assert len(exceptions) <= 20
    print("R3 NO-P-ADIC-ADIABATIC EIGENFRAME: PASS")
    print("  Res(delta-alpha numerator, delta-beta numerator) = 248832 = 2^10*3^5")
    print("  adjacent characteristic-polynomial resultant numerator = -24*R8")


def decimal_sqrt_log_ceiling(p: int) -> int:
    with localcontext() as context:
        context.prec = 80
        value = Decimal(p).sqrt() * Decimal(p).ln()
        floor = int(value.to_integral_value(rounding=ROUND_FLOOR))
        fractional_part = value - floor
        assert min(fractional_part, 1 - fractional_part) > Decimal("1e-30")
        return floor + 1


def shell_statistics(p: int) -> dict[str, object]:
    assert is_prime(p)
    N = p - 2
    D = decimal_sqrt_log_ceiling(p)
    low = D // 2
    assert 2 * D < N
    b, c = apery_orbit(p)
    fibres: dict[tuple[int, int], list[int]] = defaultdict(list)
    for r in range(1, N + 1):
        fibres[projective_key(b[r], c[r], p)].append(r)

    all_Z: dict[int, set[int]] = defaultdict(set)
    shell_gaps: dict[int, list[int]] = defaultdict(list)
    for positions in fibres.values():
        for index, r in enumerate(positions):
            for endpoint in positions[index + 1 :]:
                d = endpoint - r
                if d > D:
                    break
                all_Z[d].add(r)
                if d > low:
                    shell_gaps[r].append(d)

    assert all(len(roots) <= 3 * (d - 1) for d, roots in all_Z.items())
    for r, gaps in shell_gaps.items():
        for i, d in enumerate(gaps):
            for e in gaps[i + 1 :]:
                assert projective_key(b[r + d], c[r + d], p) == projective_key(
                    b[r + e], c[r + e], p
                )
                assert r + d in all_Z[e - d]

    histogram = Counter(len(gaps) for gaps in shell_gaps.values())
    S = sum(k * number for k, number in histogram.items())
    U = sum(histogram.values())
    singleton = histogram[1]
    pair_mass = sum(math.comb(k, 2) * number for k, number in histogram.items())
    assert S == sum(len(all_Z[d]) for d in range(low + 1, D + 1))
    assert pair_mass == sum(
        math.comb(len(gaps), 2) for gaps in shell_gaps.values()
    )
    return {
        "p": p,
        "N": N,
        "D": D,
        "S": S,
        "U": U,
        "singleton": singleton,
        "max_k": max(histogram, default=0),
        "pair_mass": pair_mass,
        "histogram": dict(sorted(histogram.items())),
    }


def empirical_gate() -> None:
    expected = {
        499: (139, 83, 73, 63, 2, 10, {1: 63, 2: 10}),
        997: (219, 161, 149, 138, 3, 13, {1: 138, 2: 10, 3: 1}),
        1999: (340, 253, 232, 213, 3, 23, {1: 213, 2: 17, 3: 2}),
        4001: (525, 371, 350, 332, 3, 24, {1: 332, 2: 15, 3: 3}),
        7919: (799, 616, 590, 564, 2, 26, {1: 564, 2: 26}),
        16001: (1225, 862, 837, 812, 2, 25, {1: 812, 2: 25}),
        32003: (1856, 1400, 1370, 1340, 2, 30, {1: 1340, 2: 30}),
        65537: (2840, 2090, 2055, 2021, 3, 36, {1: 2021, 2: 33, 3: 1}),
        99991: (3641, 2692, 2647, 2602, 2, 45, {1: 2602, 2: 45}),
    }
    results = []
    for p, target in expected.items():
        result = shell_statistics(p)
        observed = (
            result["D"],
            result["S"],
            result["U"],
            result["singleton"],
            result["max_k"],
            result["pair_mass"],
            result["histogram"],
        )
        assert observed == target
        results.append(result)

    print("APERY DYADIC-SHELL CENSUS: PASS")
    print("  p D S_B U_B singleton max_k I_2 histogram")
    for result in results:
        print(
            " ",
            result["p"],
            result["D"],
            result["S"],
            result["U"],
            result["singleton"],
            result["max_k"],
            result["pair_mass"],
            result["histogram"],
        )


def main() -> int:
    foundational_cocycle_gate()
    r1_dispersion_gate()
    r2_correlation_gate()
    r3_no_adiabatic_gate()
    empirical_gate()
    print("FINAL GATE: PASS -- all programmed exact and numerical gates passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
