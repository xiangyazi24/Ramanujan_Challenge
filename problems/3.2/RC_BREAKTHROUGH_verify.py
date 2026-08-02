#!/usr/bin/env python3
"""Verification gates for RC_BREAKTHROUGH_report.md.

The exact gates use symbolic or rational arithmetic.  Numerical evaluations
are explicitly labelled as regressions and are never used as proofs.  There
are no skipped or placeholder branches.
"""

from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path

import mpmath as mp
import sympy as sp


class GateFailure(RuntimeError):
    pass


def check(condition: bool, message: str) -> None:
    """Raise explicitly on failure, including under ``python -O``."""
    if not condition:
        raise GateFailure(message)


def P(x):
    return 34 * x**3 + 51 * x**2 + 27 * x + 5


def exact_universal_continuant_gate() -> None:
    """Universal induction gates for addition, reflection, and cell split."""

    x, j, s = sp.symbols("x j s")
    check(sp.expand(P(-x - 1) + P(x)) == 0, "P reflection identity")

    # In the addition law
    #
    # N_(j+s)(x) = N_j(x) N_(s+1)(x+j-1)
    #               -(x+j)^6 N_(j-1)(x) N_s(x+j),
    #
    # the two shifted right continuants obey the same recurrence in s.
    nj, njm1 = sp.symbols("nj njm1")
    a_prev, a_now, b_prev, b_now = sp.symbols(
        "a_prev a_now b_prev b_now"
    )
    cut = (x + j) ** 6
    coefficient = P(x + j + s)
    edge = (x + j + s) ** 6
    a_next = coefficient * a_now - edge * a_prev
    b_next = coefficient * b_now - edge * b_prev
    e_prev = nj * a_prev - cut * njm1 * b_prev
    e_now = nj * a_now - cut * njm1 * b_now
    e_next = nj * a_next - cut * njm1 * b_next
    check(
        sp.expand(e_next - coefficient * e_now + edge * e_prev) == 0,
        "universal continuant addition induction",
    )
    check(nj * 1 - cut * njm1 * 0 == nj, "continuant addition base s=0")
    check(
        sp.expand(
            nj * P(x + j)
            - cut * njm1
            - (P(x + j) * nj - cut * njm1)
        )
        == 0,
        "continuant addition base s=1",
    )

    # Expansion from the left endpoint satisfies the same generic induction
    # as expansion from the right endpoint:
    # N_(m+1)(x)=P(x+1)N_m(x+1)-(x+2)^6N_(m-1)(x+2).
    m = sp.symbols("m")
    aa, bb, cc, dd = sp.symbols("aa bb cc dd")
    right_coefficient = P(x + m + 1)
    right_edge = (x + m + 1) ** 6
    left_coefficient = P(x + 1)
    left_edge = (x + 2) ** 6
    via_right = (
        right_coefficient * (left_coefficient * aa - left_edge * cc)
        - right_edge * (left_coefficient * bb - left_edge * dd)
    )
    via_left = (
        left_coefficient * (right_coefficient * aa - right_edge * bb)
        - left_edge * (right_coefficient * cc - right_edge * dd)
    )
    check(sp.expand(via_right - via_left) == 0, "left-end induction step")
    check(
        sp.expand(P(x + 1) - (P(x + 1) * 1 - (x + 2) ** 6 * 0)) == 0,
        "left-end recurrence base",
    )

    # Apply the left-end recurrence at -h-2-x.  The two induction signs are
    # epsilon and -epsilon, yielding the reflection sign -epsilon at h+1.
    h = sp.symbols("h")
    epsilon, nh, nhm1 = sp.symbols("epsilon nh nhm1")
    ordinary_next = P(x + h) * nh - (x + h) ** 6 * nhm1
    reflected_next = (
        P(-h - 1 - x) * (epsilon * nh)
        - (-h - x) ** 6 * (-epsilon * nhm1)
    )
    check(
        sp.expand(P(-h - 1 - x) + P(x + h)) == 0,
        "universal reflected coefficient",
    )
    check(
        sp.expand(reflected_next - (-epsilon) * ordinary_next) == 0,
        "universal continuant reflection induction",
    )
    n1 = sp.Integer(1)
    check(sp.expand(n1 - (-1) ** 0 * n1) == 0, "reflection base h=1")
    check(
        sp.expand(P(-3 - x + 1) + P(x + 1)) == 0,
        "reflection base h=2",
    )

    # Normalize the addition law at x=-j+z.  Reflection gives the minus sign
    # in N_(j-1)/left_den; the two minuses then produce the companion term.
    z = sp.symbols("z", nonzero=True)
    left_den, right_den = sp.symbols(
        "left_den right_den", nonzero=True
    )
    f_left, h_left, f_right, h_right = sp.symbols(
        "f_left h_left f_right h_right"
    )
    n_j = left_den * f_left
    n_jm1 = -left_den * h_left
    n_sp1 = right_den * f_right
    n_s = right_den * h_right
    normalized_cell = sp.cancel(
        (n_j * n_sp1 - z**6 * n_jm1 * n_s)
        / (z**3 * left_den * right_den)
    )
    expected_cell = f_left * f_right / z**3 + z**3 * h_left * h_right
    check(
        sp.cancel(normalized_cell - expected_cell) == 0,
        "universal normalized two-block cell",
    )

    print("UNIVERSAL-CONTINUANT addition_reflection_cell PASS")


def exact_shift_and_casoratian_gate() -> None:
    """Check the finite basis shifts and periodic Casoratian algebra."""

    z = sp.symbols("z")
    check(sp.expand(P(-z - 1) + P(z)) == 0, "P reflection identity")

    def fg_table(u, nmax):
        f = [sp.Integer(1), sp.cancel(P(u) / (u + 1) ** 3)]
        g = [sp.Integer(0), sp.cancel(1 / (u + 1) ** 3)]
        for n in range(1, nmax):
            f.append(
                sp.cancel(
                    (P(u + n) * f[n] - (u + n) ** 3 * f[n - 1])
                    / (u + n + 1) ** 3
                )
            )
            g.append(
                sp.cancel(
                    (P(u + n) * g[n] - (u + n) ** 3 * g[n - 1])
                    / (u + n + 1) ** 3
                )
            )
        return f, g

    checks = 0
    for u in (sp.Rational(2, 7), sp.Rational(-3, 8), sp.Rational(11, 13)):
        f0, g0 = fg_table(u, 9)
        f1, g1 = fg_table(u + 1, 8)
        for n in range(1, 8):
            check(
                sp.cancel(f1[n - 1] - (u + 1) ** 3 * g0[n]) == 0,
                f"sampled F/G shift u={u},n={n}",
            )
            check(
                sp.cancel(
                    g1[n]
                    - (P(u) * g0[n + 1] - f0[n + 1]) / (u + 1) ** 3
                )
                == 0,
                f"sampled companion shift u={u},n={n}",
            )
            checks += 2

    # Universal induction step for F_(n-1)(u+1)=(u+1)^3 G_n(u).
    # The symbols represent two consecutive G-values; no sampled n or u is
    # used in this gate.
    u, nn = sp.symbols("u nn")
    g_nm1, g_n = sp.symbols("g_nm1 g_n")
    shifted_f_nm2 = (u + 1) ** 3 * g_nm1
    shifted_f_nm1 = (u + 1) ** 3 * g_n
    shifted_f_n = (
        P(u + nn) * shifted_f_nm1
        - (u + nn) ** 3 * shifted_f_nm2
    ) / (u + nn + 1) ** 3
    g_np1 = (
        P(u + nn) * g_n - (u + nn) ** 3 * g_nm1
    ) / (u + nn + 1) ** 3
    check(
        sp.cancel(shifted_f_n - (u + 1) ** 3 * g_np1) == 0,
        "universal F/G shift induction",
    )
    check(
        sp.cancel(1 - (u + 1) ** 3 / (u + 1) ** 3) == 0,
        "F/G shift base normalization",
    )

    # Independently reconstruct the gap continuants and check the exact
    # two-block cell formula for every pole split through height seven.
    x = sp.symbols("x")
    numerators = [sp.Integer(0), sp.Integer(1)]
    for h in range(1, 8):
        numerators.append(
            sp.expand(
                P(x + h) * numerators[h]
                - (x + h) ** 6 * numerators[h - 1]
            )
        )
    for h in range(1, 8):
        qh = sp.prod(x + a for a in range(1, h + 1))
        for j in range(1, h + 1):
            r, right = j - 1, h - j
            f_left, g_left = fg_table(-z, max(2, r + 1))
            f_right, g_right = fg_table(z, max(2, right + 1))
            lhs = numerators[h].subs(x, -j + z) / qh.subs(x, -j + z) ** 3
            rhs = (
                f_left[r] * f_right[right] / z**3
                + z**3 * g_left[r] * g_right[right]
            )
            check(
                sp.cancel(lhs - rhs) == 0,
                f"two-block cell identity h={h},j={j}",
            )
            checks += 1

    # Abel identity for any two solutions of the reflected scalar recurrence.
    u0, um, v0, vm = sp.symbols("u0 um v0 vm")
    up = (P(z) * u0 - z**3 * um) / (z + 1) ** 3
    vp = (P(z) * v0 - z**3 * vm) / (z + 1) ** 3
    phi_z = z**3 * (u0 * vm - um * v0)
    phi_z1 = (z + 1) ** 3 * (up * v0 - u0 * vp)
    check(sp.factor(phi_z1 - phi_z) == 0, "Abel identity")

    # Equivalent two-by-two bilinear cocycle.
    lam = sp.symbols("lam", nonzero=True)
    mat = lambda x: lam * sp.Matrix(
        [[0, (x + 1) ** 3], [-(x + 1) ** -3, P(x) * (x + 1) ** -3]]
    )
    dmat = lambda x: sp.diag(1, x**6)
    err = sp.simplify(
        mat(-z - 1).inv().T * dmat(z + 1) * mat(z)
        - ((z + 1) / z) ** 3 * dmat(z)
    )
    check(err == sp.zeros(2), "two-by-two bilinear cocycle")
    print(f"SHIFT-CASORATIAN exact_checks={checks} PASS")


def exact_global_cell_cocycle_gate() -> None:
    """Check the complete two-component BV/Jost cell cocycle exactly."""

    s = sp.symbols("s", nonzero=True)
    d = sp.symbols("d")
    # Formal adjunction sends D to -D and Weyl reordering sends
    # A(D)t^j to t^j A(D+j).  The three displayed identities encode
    # D^3, -tP(D), and t^2(D+1)^3 respectively after reordering.
    check(sp.expand((-d) ** 3 + d**3) == 0, "adjoint D^3 sign")
    middle_adjoint_reordered = -P(-d - 1)
    last_adjoint_reordered = (1 - (d + 2)) ** 3
    check(
        sp.expand(middle_adjoint_reordered - P(d)) == 0,
        "adjoint middle-term reordering",
    )
    check(
        sp.expand(last_adjoint_reordered + (d + 1) ** 3) == 0,
        "adjoint last-term reordering",
    )
    phi_s, psi_s, phi_ms, psi_ms = sp.symbols(
        "phi_s psi_s phi_ms psi_ms"
    )

    # The BV three-term equation is equivalent to these first-order shifts.
    phi_sp1 = (s + 1) ** 3 * psi_s
    psi_sp1 = (P(s) * psi_s - phi_s) / (s + 1) ** 3
    psi_msm1 = -phi_ms / s**3
    phi_msm1 = s**3 * psi_ms - P(s) * psi_msm1

    j_s = phi_ms * phi_s / s**3 + s**3 * psi_ms * psi_s
    j_sp1 = (
        phi_msm1 * phi_sp1 / (s + 1) ** 3
        + (s + 1) ** 3 * psi_msm1 * psi_sp1
    )
    check(sp.factor(j_sp1 - j_s) == 0, "global cell periodicity cocycle")

    # The same calculation starts from Proposition 8 with the exact Apéry
    # differential-operator coefficients and P(-s-1)=-P(s).
    gamma_s, gamma_sp1, gamma_sp2 = sp.symbols(
        "gamma_s gamma_sp1 gamma_sp2"
    )
    proposition8 = (
        (-s) ** 3 * gamma_s
        - P(-s - 1) * gamma_sp1
        + (-s - 1) ** 3 * gamma_sp2
    )
    expected = -(
        s**3 * gamma_s - P(s) * gamma_sp1 + (s + 1) ** 3 * gamma_sp2
    )
    check(
        sp.expand(proposition8 - expected) == 0,
        "Bloch--Vlasenko recurrence sign",
    )

    # The c_0^s normalization cancels from both bilateral products.  Encoding
    # c_0^s by an invertible symbol keeps this an algebraic, branch-free gate.
    gauge, raw_p, raw_m, comp_p, comp_m = sp.symbols(
        "gauge raw_p raw_m comp_p comp_m", nonzero=True
    )
    check(
        sp.cancel((gauge * raw_p) * (raw_m / gauge) - raw_p * raw_m) == 0,
        "principal bilateral gauge cancellation",
    )
    check(
        sp.cancel((gauge * comp_p) * (comp_m / gauge) - comp_p * comp_m)
        == 0,
        "companion bilateral gauge cancellation",
    )

    # Exact local principal part.  The companion is analytic and multiplied by
    # s^3, so it cannot affect the cubic or residue coefficients.
    a1, a2, a3, q0, q1 = sp.symbols("a1 a2 a3 q0 q1")
    kap_p = 1 + a1 * s + a2 * s**2 + a3 * s**3
    kap_m = 1 - a1 * s + a2 * s**2 - a3 * s**3
    companion = s**3 * (q0 - q1 * s) * (q0 + q1 * s)
    local = sp.expand(kap_m * kap_p / s**3 + companion)
    check(sp.expand(local).coeff(s, -3) == 1, "global cell cubic coefficient")
    residue = sp.expand(local).coeff(s, -1)
    check(residue == 2 * a2 - a1**2, "global cell formal residue")
    check(
        residue.subs({a1: 0, a2: -sp.pi**2 / 3})
        == -2 * sp.pi**2 / 3,
        "global cell Apéry residue",
    )

    print("GLOBAL-CELL-COCYCLE recurrence_periodicity_gauge_residue PASS")


def exact_reflection_endpoint_gate() -> None:
    """Check the Apéry reflection exponent needed in BV Proposition 15."""

    t, rho = sp.symbols("t rho")
    c0 = 17 - 12 * sp.sqrt(2)
    # Ordinary-derivative coefficients of y''' and y'' in
    # D^3-tP(D)+t^2(D+1)^3, with D=t*d/dt.
    a3 = t**3 * (t**2 - 34 * t + 1)
    a2 = 3 * t**2 * (2 * t**2 - 51 * t + 1)
    check(sp.simplify(a3.subs(t, c0)) == 0, "reflection point is singular")
    ratio = sp.simplify(a2.subs(t, c0) / sp.diff(a3, t).subs(t, c0))
    check(ratio == sp.Rational(3, 2), "reflection indicial coefficient")
    indicial = sp.factor(rho * (rho - 1) * (rho - 2 + ratio))
    check(
        sp.factor(
            indicial - rho * (rho - 1) * (rho - sp.Rational(1, 2))
        )
        == 0,
        "reflection indicial polynomial",
    )
    print("REFLECTION-ENDPOINT indicial_exponents=0,1/2,1 PASS")


def exact_gz_telescoper_gate() -> None:
    """Check the GZ Proposition 3 certificate through epsilon squared."""

    n, k, e = sp.symbols("n k e")

    def exp2(d1, d2):
        return 1 + d1 * e + (d1**2 + d2) * e**2 / 2

    base_plus = ((n + k + 1) / (n - k + 1)) ** 2
    base_minus = ((n - k) / (n + k)) ** 2
    base_kminus = k**4 / ((n + k) ** 2 * (n - k + 1) ** 2)
    ratio_plus = base_plus * exp2(
        4 / (n + k + 1), -8 / (n + k + 1) ** 2
    )
    ratio_minus = base_minus * exp2(-4 / (n + k), 8 / (n + k) ** 2)
    ratio_kminus = base_kminus * exp2(
        4 * n / (k * (n + k)), -4 / k**2 + 8 / (n + k) ** 2
    )

    def cert(nn, kk):
        return (
            4 * (2 * nn + 1) * (2 * kk**2 + kk - (2 * nn + 1) ** 2)
            + (
                16 * kk**2
                + 8 * (4 * nn + 3) * kk
                - 4 * (2 * nn + 1) * (12 * nn + 5)
            )
            * e
            + 16 * (2 * kk - 5 * nn - 2) * e**2
        )

    lhs = (
        (n + 1 + e) ** 3 * ratio_plus
        - P(n + e)
        + (n + e) ** 3 * ratio_minus
    )
    rhs = cert(n, k) - cert(n, k - 1) * ratio_kminus
    diff = sp.series(lhs - rhs, e, 0, 3).removeO()
    for degree in range(3):
        check(
            sp.factor(sp.together(sp.expand(diff).coeff(e, degree))) == 0,
            f"GZ telescoper epsilon degree {degree}",
        )

    lam = 17 + 12 * sp.sqrt(2)
    check(sp.expand((1 + sp.sqrt(2)) ** 4 - lam) == 0, "Apéry growth root")
    print("GZ-TELESCOPER mod_epsilon^3 PASS")


def exact_hypergeometric_anomaly_gate() -> None:
    """Check the rational Gosper certificate for the balanced 4F3 anomaly."""

    k, x = sp.symbols("k x")
    term_ratio = ((k - x) * (k + x + 1) / (k + 1) ** 2) ** 2
    shifted = (
        (x + 1) ** 3 * ((x + k + 1) / (x + 1 - k)) ** 2
        - P(x)
        + x**3 * ((x - k) / (x + k)) ** 2
    )
    cert = (
        -4
        * k**4
        * (2 * x + 1)
        * (-2 * k**2 + 3 * k + 4 * x**2 + 4 * x)
        / ((k + x) ** 2 * (-k + x + 1) ** 2)
    )
    check(
        sp.factor(term_ratio * cert.subs(k, k + 1) - cert - shifted) == 0,
        "hypergeometric Gosper certificate",
    )
    check(
        sp.factor(sp.limit(cert / k**2, k, sp.oo) - 8 * (2 * x + 1))
        == 0,
        "hypergeometric certificate boundary",
    )

    # Independent high-precision evaluations of the summed identity.
    mp.mp.dps = 60

    def bfun(w):
        return mp.hyper([-w, -w, w + 1, w + 1], [1, 1, 1], 1)

    worst = mp.mpf(0)
    for w in (mp.mpf("0.1"), mp.mpf("0.2"), mp.mpc("0.2", "0.1")):
        got = (w + 1) ** 3 * bfun(w + 1) - P(w) * bfun(w) + w**3 * bfun(w - 1)
        want = 8 * (2 * w + 1) * (mp.sin(mp.pi * w) / mp.pi) ** 2
        worst = max(worst, abs(got - want))
    check(worst < mp.mpf("1e-50"), "hypergeometric numerical regression")
    print(f"HYPERGEOMETRIC-ANOMALY exact_certificate numeric_error={mp.nstr(worst, 5)} PASS")


def exact_positive_remainder_gate() -> None:
    """Check the second WZ identity giving the positive zeta(2) remainder."""

    n, t = sp.symbols("n t")
    r_plus = ((t - n - 1) / (t + n + 1)) ** 2
    r_minus = ((t + n) / (t - n)) ** 2
    r_tplus = t**4 / ((t - n) ** 2 * (t + n + 1) ** 2)
    s = lambda nn, tt: 4 * (2 * nn + 1) * ((2 * nn + 1) ** 2 + tt - 2 * tt**2)
    lhs = (n + 1) ** 3 * r_plus - P(n) + n**3 * r_minus
    rhs = s(n, t + 1) * r_tplus - s(n, t)
    check(sp.factor(sp.together(lhs - rhs)) == 0, "positive-remainder WZ identity")
    check(
        sp.limit(s(n, t) / t**2, t, sp.oo) == -8 * (2 * n + 1),
        "positive-remainder certificate boundary",
    )

    # The explicit positive kernel used in the infinite-tail proof.  These
    # identities are rational in t for every fixed n; checking a range guards
    # the product indexing independently of the Gosper certificate above.
    def kernel(nn):
        out = t**-2
        for k in range(nn + 1):
            out *= ((t - k) / (t + k)) ** 2
        return sp.cancel(out)

    for nn in range(0, 10):
        wn = kernel(nn)
        check(sp.limit(t**2 * wn, t, sp.oo) == 1, f"kernel tail n={nn}")
        for m in range(1, nn + 1):
            check(
                sp.cancel(wn.subs(t, m)) == 0,
                f"kernel zero n={nn},m={m}",
            )
        if nn < 9:
            check(
                sp.factor(
                    kernel(nn + 1) / wn
                    - ((t - nn - 1) / (t + nn + 1)) ** 2
                )
                == 0,
                f"kernel height shift n={nn}",
            )
        check(
            sp.factor(
                wn.subs(t, t + 1) / wn
                - t**4 / ((t - nn) ** 2 * (t + nn + 1) ** 2)
            )
            == 0,
            f"kernel argument shift n={nn}",
        )
        check(
            sp.limit(s(nn, t) * wn, t, sp.oo) == -8 * (2 * nn + 1),
            f"weighted kernel boundary n={nn}",
        )

    # Exact recurrence comparison for the rational part e_n.
    b = [Fraction(1), Fraction(5)]
    e = [Fraction(0), Fraction(8)]
    for j in range(1, 16):
        b.append((P(j) * b[j] - j**3 * b[j - 1]) / (j + 1) ** 3)
        e.append(
            (P(j) * e[j] - j**3 * e[j - 1] + 8 * (2 * j + 1))
            / (j + 1) ** 3
        )
    check(all(x.denominator == 1 for x in b), "Apéry recurrence integrality")
    check(
        e[1:4] == [Fraction(8), Fraction(120), Fraction(21392, 9)],
        "positive-remainder initial data",
    )
    # The positive remainder zeta(2)b_n-e_n is bounded directly in the report;
    # this finite monotonicity check is corroboration, not its proof.  The
    # report denotes this integer sequence e_n by q_n.
    ratios = [e[j] / b[j] for j in range(1, len(b))]
    check(
        all(ratios[j] > ratios[j - 1] for j in range(1, len(ratios))),
        "positive-remainder finite monotonicity",
    )
    print("POSITIVE-REMAINDER exact_WZ_kernel_and_recurrence PASS")


def exact_trigonometric_gate() -> None:
    """Check the Laurent coefficient, factorization, zeros, and critical data."""

    z = sp.symbols("z")
    pi = sp.pi
    phi = pi**3 * sp.cot(pi * z) ** 3 + pi**3 * sp.cot(pi * z) / 3
    factored = (
        pi**3
        * sp.cos(pi * z)
        * (sp.cos(2 * pi * z) + 2)
        / (3 * sp.sin(pi * z) ** 3)
    )
    check(sp.trigsimp(phi - factored) == 0, "trigonometric factorization")
    check(
        sp.trigsimp(phi.subs(z, z + 1) - phi) == 0,
        "trigonometric periodicity",
    )
    check(
        sp.trigsimp(phi.subs(z, -z) + phi) == 0,
        "trigonometric oddness",
    )
    series = sp.series(phi, z, 0, 3).removeO()
    expected = z**-3 - 2 * pi**2 / (3 * z) + 7 * pi**4 * z / 45
    check(sp.simplify(series - expected) == 0, "trigonometric Laurent series")

    derivative = -pi**4 / 3 * sp.csc(pi * z) ** 2 * (
        9 * sp.cot(pi * z) ** 2 + 1
    )
    check(
        sp.trigsimp(sp.diff(phi, z) - derivative) == 0,
        "trigonometric derivative",
    )
    cylinder_cot = -sp.I
    cylinder_value = pi**3 * (
        cylinder_cot**3 + cylinder_cot / 3
    )
    cylinder_derivative = -pi**4 / 3 * (
        1 + cylinder_cot**2
    ) * (9 * cylinder_cot**2 + 1)
    check(
        sp.simplify(cylinder_value - 2 * sp.I * pi**3 / 3) == 0,
        "upper-cylinder limiting value",
    )
    check(
        sp.simplify(cylinder_derivative) == 0,
        "upper-cylinder limiting derivative",
    )
    eta = sp.acosh(2) / (2 * pi)
    ycrit = sp.log(2) / (2 * pi)
    check(sp.simplify(sp.cosh(2 * pi * eta) - 2) == 0, "cell zero height")
    check(
        sp.simplify(sp.tanh(pi * ycrit) - sp.Rational(1, 3)) == 0,
        "critical height equation",
    )
    print("TRIGONOMETRIC-CELL exact_formula_zero_and_critical_heights PASS")


def wall_budget_gate() -> None:
    """Check the two elementary budgets used in the negative wall verdict."""

    # The parity-forced contribution is summable; a linear accidental bound is harmonic.
    for D in (10, 100, 1000):
        forced = sum(Fraction(1, h * h) for h in range(2, D + 1, 2))
        linear = sum(Fraction(3 * h, h * h) for h in range(2, D + 1))
        check(forced < Fraction(1, 2), f"summable parity budget D={D}")
        check(
            linear > 2 * sum(Fraction(1, h) for h in range(2, D + 1)),
            f"harmonic accidental budget D={D}",
        )
    # Sum_{a,g<=D} a*g is quartic scale, explaining why the current resultant
    # height budget cannot remove the Q_D logarithm at D~sqrt(p)L.
    vals = []
    for D in (20, 40, 80):
        vals.append(sum(a * g for a in range(2, D) for g in range(2, D - a + 1)))
    check(
        vals[1] > 12 * vals[0] and vals[2] > 13 * vals[1],
        "quartic resultant-height regression",
    )
    print(f"WALL-BUDGET harmonic_and_quartic={vals} PASS")


def certified_range_manifest_gate() -> None:
    """Audit every archived Arb certificate row for heights 2 through 60.

    This checks the certificate manifest, not the interval arithmetic itself.
    A fresh Arb recomputation uses the command printed in the report.
    """

    path = Path(__file__).with_name("CRON_kinf_results.json")
    data = json.loads(path.read_text())
    scan = data["scan"]
    heights = list(range(2, 61))
    check(scan["min_h"] == 2 and scan["max_h"] == 60, "Arb scan endpoints")
    check(scan["certified_yes"] == heights, "Arb certified-height inventory")
    rows = data["rows"]
    check([row["h"] for row in rows] == heights, "Arb row inventory")
    for row in rows:
        h = row["h"]
        cert = row["certified"]
        roots = 2 * h - 2
        check(cert["status"] == "YES", f"Arb status h={h}")
        check(cert["exhaustive"] is True, f"Arb exhaustive h={h}")
        check(cert["root_balls_disjoint"] is True, f"Arb disjoint roots h={h}")
        check(cert["quotient_roots"] == roots, f"Arb quotient roots h={h}")
        check(cert["critical_zeros"] == 2 * roots, f"Arb critical zeros h={h}")
        check(cert["interval_newton_pass"] == roots, f"Arb Newton count h={h}")
        check(
            cert["denominator_excludes_zero"] == roots,
            f"Arb denominator exclusions h={h}",
        )
        check(
            cert["derivative_excludes_zero"] == roots,
            f"Arb derivative exclusions h={h}",
        )
        check(cert["good_orbits"] == roots, f"Arb good orbits h={h}")
        check(float(cert["certified_margin"]) > 0, f"Arb margin h={h}")
        check(float(cert["relative_margin"]) > 0, f"Arb relative margin h={h}")
    print("ARB-MANIFEST archived_certificates_h=2..60 PASS")


if __name__ == "__main__":
    exact_universal_continuant_gate()
    exact_shift_and_casoratian_gate()
    exact_global_cell_cocycle_gate()
    exact_reflection_endpoint_gate()
    exact_gz_telescoper_gate()
    exact_hypergeometric_anomaly_gate()
    exact_positive_remainder_gate()
    exact_trigonometric_gate()
    wall_budget_gate()
    certified_range_manifest_gate()
    print("PASS")
