#!/usr/bin/env python3
"""Exact standard-library audit for Q7305.

The three independent constructions are:
  (1) triangular elimination of K_p(x)^2 in the Psi basis modulo x^p;
  (2) the Lagrange-inversion binomial formula for [a^j] R(a)^p;
  (3) reduction c_{p,m} == b_m (mod p).

The exploratory identity scan is deliberately finite and printed verbatim: it
checks the natural raw/first-order quantities named in Q7305, rather than
claiming an exhaustive search through arbitrary formulas.
"""

from hashlib import sha256
from math import comb
from pathlib import Path

LIMIT = 199


def primes_upto(n):
    out = []
    for q in range(2, n + 1):
        if all(q % d for d in range(2, int(q ** 0.5) + 1)):
            out.append(q)
    return out


def apery_sequence(nmax):
    """Return b_0,...,b_nmax as exact integers."""
    if nmax == 0:
        return [1]
    b = [1, 5]
    for n in range(1, nmax):
        num = (2 * n + 1) * (17 * n * n + 17 * n + 5) * b[n] - n**3 * b[n - 1]
        den = (n + 1) ** 3
        q, r = divmod(num, den)
        assert r == 0
        b.append(q)
    return b


def franel_sequence(nmax):
    return [sum(comb(n, k) ** 3 for k in range(n + 1)) for n in range(nmax + 1)]


def mul_trunc(a, b, n):
    out = [0] * min(n, len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai:
            stop = min(len(b), n - i)
            for j in range(stop):
                out[i + j] += ai * b[j]
    if len(out) < n:
        out.extend([0] * (n - len(out)))
    return out


def psi_trunc(p, m):
    """Coefficients below x^p of x^m(1+x)^(p-1-m)(1-8x)^m."""
    left = [comb(p - 1 - m, i) for i in range(p - m)]
    right = [comb(m, j) * (-8) ** j for j in range(m + 1)]
    core = mul_trunc(left, right, p - m)
    return [0] * m + core


def basis_elimination(p, franel):
    """Unitriangular elimination in Z[x]/(x^p)."""
    target = mul_trunc(franel[:p], franel[:p], p)
    residual = target[:]
    coords = []
    for m in range(p):
        psi = psi_trunc(p, m)
        assert len(psi) == p and psi[m] == 1 and all(psi[i] == 0 for i in range(m))
        c = residual[m]
        coords.append(c)
        for i in range(m, p):
            residual[i] -= c * psi[i]
    assert residual == [0] * p
    return coords


def r_binomial(p, j):
    """[a^j] R(a)^p from Lagrange inversion; valid here for 0 <= j < p."""
    if j == 0:
        return 1
    inner = 0
    for ell in range(j):
        inner += ((-1) ** (j - 1 - ell) * 8**ell
                  * comb(j + ell - 1, ell)
                  * comb(p - 1 - ell, j - 1 - ell))
    num = -p * inner
    q, rem = divmod(num, j)
    assert rem == 0
    return q


def coefficient_formula(p, apery):
    r = [r_binomial(p, j) for j in range(p)]
    return [sum(apery[m - j] * r[j] for j in range(m + 1)) for m in range(p)]


def first_failure(cases, predicate):
    total = 0
    passed = 0
    first = None
    for item in cases:
        total += 1
        if predicate(*item):
            passed += 1
        elif first is None:
            first = item
    return passed, total, first


def fmt_failure(item):
    return "none" if item is None else repr(item)


def main():
    primes = primes_upto(LIMIT)
    nmax = 2 * LIMIT
    b = apery_sequence(nmax)
    f = franel_sequence(nmax)

    full_hash = sha256()
    zero_hash = sha256()
    all_c = {}
    all_d = {}
    zero_rows = []
    basis_formula_pairs = 0
    congruence_pairs = 0

    for p in primes:
        c_basis = basis_elimination(p, f)
        c_formula = coefficient_formula(p, b)
        assert c_basis == c_formula
        all_c[p] = c_basis
        drow = []
        for m, c in enumerate(c_basis):
            basis_formula_pairs += 1
            assert (c - b[m]) % p == 0
            congruence_pairs += 1
            d = (c - b[m]) // p
            drow.append(d)
            full_hash.update(f"{p},{m},{b[m]},{f[m]},{c},{d}\n".encode())
            if b[m] % p == 0:
                assert c % p == 0
                q0 = b[m] // p
                cp = c // p
                qplus = b[p + m] // p
                qref = b[p - 1 - m] // p
                shift_b = (b[p + m] - 5 * b[m]) // p
                refl_b = (b[p - 1 - m] - b[m]) // p
                shift_f_num = f[p + m] - 2 * f[m]
                assert shift_f_num % p == 0
                shift_f = shift_f_num // p
                row = {
                    "p": p,
                    "m": m,
                    "d": d,
                    "cp": cp,
                    "q0": q0,
                    "qplus": qplus,
                    "qref": qref,
                    "fm": f[m],
                    "fref": f[p - 1 - m],
                    "fplus": f[p + m],
                    "shift_b": shift_b,
                    "refl_b": refl_b,
                    "shift_f": shift_f,
                }
                zero_rows.append(row)
                zero_hash.update(
                    (",".join(str(row[k]) for k in
                              ("p", "m", "d", "cp", "q0", "qplus", "qref",
                               "fm", "fref", "fplus", "shift_b", "refl_b", "shift_f"))
                     + "\n").encode()
                )
        all_d[p] = drow

    print("Q7305 EXACT AUDIT")
    print(f"limit={LIMIT}")
    print(f"prime_count={len(primes)}")
    print(f"coordinate_pairs={basis_formula_pairs}")
    print(f"basis_equals_formula={basis_formula_pairs}/{basis_formula_pairs}")
    print(f"c_congruent_b_mod_p={congruence_pairs}/{congruence_pairs}")
    print(f"zero_pairs={len(zero_rows)}")
    print("zero_positions=" + repr([(r["p"], r["m"]) for r in zero_rows]))
    print(f"full_dataset_sha256={full_hash.hexdigest()}")
    print(f"zero_dataset_sha256={zero_hash.hexdigest()}")
    print(f"script_sha256={sha256(Path(__file__).read_bytes()).hexdigest()}")

    all_pm = [(p, m) for p in primes for m in range(p)]
    all_pm_ge5 = [(p, m) for p in primes if p >= 5 for m in range(p)]
    zero_pm = [(r["p"], r["m"]) for r in zero_rows]

    universal_tests = [
        ("Apery reflection b[p-1-m] == b[m] (mod p)", all_pm,
         lambda p, m: (b[p - 1 - m] - b[m]) % p == 0),
        ("Apery one-digit shift b[p+m] == 5*b[m] (mod p^2), p>=5", all_pm_ge5,
         lambda p, m: (b[p + m] - 5 * b[m]) % (p * p) == 0),
        ("Apery one-digit shift b[p+m] == 5*b[m] (mod p^3), p>=5", all_pm_ge5,
         lambda p, m: (b[p + m] - 5 * b[m]) % (p**3) == 0),
        ("Franel one-digit shift f[p+m] == 2*f[m] (mod p)", all_pm,
         lambda p, m: (f[p + m] - 2 * f[m]) % p == 0),
        ("Franel one-digit shift f[p+m] == 2*f[m] (mod p^2), p>=5", all_pm_ge5,
         lambda p, m: (f[p + m] - 2 * f[m]) % (p * p) == 0),
        ("Apery reflection lifts to mod p^2 on zero positions", zero_pm,
         lambda p, m: (b[p - 1 - m] - b[m]) % (p * p) == 0),
    ]

    print("UNIVERSAL_AND_ZERO_RESTRICTED_TESTS")
    for label, cases, pred in universal_tests:
        passed, total, first = first_failure(cases, pred)
        status = "PASS" if passed == total else "FAIL"
        print(f"{status}: {label}: {passed}/{total}; first_failure={fmt_failure(first)}")

    # Natural first-order quantities on an Apéry zero b_m = 0 (mod p).
    # Every equality below is interpreted modulo the row prime p.
    candidate_tests = [
        ("qplus == 5*q0", lambda r: r["qplus"] - 5 * r["q0"]),
        ("qref == q0", lambda r: r["qref"] - r["q0"]),
        ("d == q0", lambda r: r["d"] - r["q0"]),
        ("d == -q0", lambda r: r["d"] + r["q0"]),
        ("cp == q0", lambda r: r["cp"] - r["q0"]),
        ("cp == qref", lambda r: r["cp"] - r["qref"]),
        ("d == refl_b", lambda r: r["d"] - r["refl_b"]),
        ("d == -refl_b", lambda r: r["d"] + r["refl_b"]),
        ("d == shift_b", lambda r: r["d"] - r["shift_b"]),
        ("d == -shift_b", lambda r: r["d"] + r["shift_b"]),
        ("d == fm", lambda r: r["d"] - r["fm"]),
        ("d == -fm", lambda r: r["d"] + r["fm"]),
        ("cp == fm", lambda r: r["cp"] - r["fm"]),
        ("cp == -fm", lambda r: r["cp"] + r["fm"]),
        ("d == fref", lambda r: r["d"] - r["fref"]),
        ("d == -fref", lambda r: r["d"] + r["fref"]),
        ("d == shift_f", lambda r: r["d"] - r["shift_f"]),
        ("d == -shift_f", lambda r: r["d"] + r["shift_f"]),
        ("cp == 0 (equiv p^2 divides c)", lambda r: r["cp"]),
    ]

    print("ZERO_POSITION_CANDIDATES_MOD_P")
    for label, diff in candidate_tests:
        matches = [r for r in zero_rows if diff(r) % r["p"] == 0]
        first = next((r for r in zero_rows if diff(r) % r["p"] != 0), None)
        status = "PASS" if len(matches) == len(zero_rows) else "CORRELATION_ONLY"
        first_pair = None if first is None else (first["p"], first["m"])
        print(f"{status}: {label}: {len(matches)}/{len(zero_rows)}; first_failure={first_pair}")

    # Constant-multiplier discovery over a fixed small range.  This is a scan,
    # not a theorem generator.  Report only identities valid on every zero row.
    variables = ["d", "cp", "q0", "qplus", "qref", "fm", "fref",
                 "shift_b", "refl_b", "shift_f"]
    multiplier_hits = []
    for lhs in variables:
        for rhs in variables:
            if lhs == rhs:
                continue
            for k in range(-16, 17):
                if all((r[lhs] - k * r[rhs]) % r["p"] == 0 for r in zero_rows):
                    multiplier_hits.append((lhs, k, rhs))
    print("ALL_ZERO_SMALL_MULTIPLIER_IDENTITIES_MOD_P=" + repr(multiplier_hits))

    # Exact equality scan for the named raw quantities.  Trivial definitions
    # (cp=q0+d) are omitted from the pair scan and checked separately.
    exact_vars = ["d", "cp", "q0", "qplus", "qref", "fm", "fref", "fplus"]
    exact_hits = []
    for i, lhs in enumerate(exact_vars):
        for rhs in exact_vars[i + 1:]:
            if all(r[lhs] == r[rhs] for r in zero_rows):
                exact_hits.append((lhs, rhs))
    assert all(r["cp"] == r["q0"] + r["d"] for r in zero_rows)
    print("ALL_ZERO_RAW_EXACT_PAIR_IDENTITIES=" + repr(exact_hits))
    print("DEFINITIONAL_EXACT_IDENTITY=cp=q0+d: PASS")

    # Size criterion.  A nonzero c with |c|<p certifies m is not an Apéry zero.
    small = []
    small_zero = []
    small_interior = []
    for p in primes:
        for m, c in enumerate(all_c[p]):
            if abs(c) < p:
                item = (p, m, c, b[m] % p)
                small.append(item)
                if b[m] % p == 0:
                    small_zero.append(item)
                if 2 <= m <= p - 3 and c != 0:
                    small_interior.append(item)
    grouped = {}
    for p, m, c, residue in small:
        grouped.setdefault(m, []).append((p, c, residue))
    print(f"ABS_C_LT_P_COUNT={len(small)}")
    print("ABS_C_LT_P_BY_M=" + repr(grouped))
    print("ABS_C_LT_P_AT_ACTUAL_ZEROS=" + repr(small_zero))
    print("NONZERO_ABS_C_LT_P_INTERIOR_CERTIFICATES=" + repr(small_interior))

    # A few fixed coefficients catch transcription/sign errors immediately.
    spot = []
    for p in (2, 3, 5, 7, 11, 199):
        upto = min(4, p)
        spot.append((p, all_c[p][:upto], all_d[p][:upto]))
    print("SPOT_CHECKS_(p,c_prefix,d_prefix)=" + repr(spot))
    print("FINAL_STATUS=PASS")


if __name__ == "__main__":
    main()
