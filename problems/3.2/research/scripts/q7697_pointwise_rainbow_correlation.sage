#!/usr/bin/env sage
"""Q7697: exact cross-prime early-CRT correlation experiment for Apéry zeros.

This is deliberately NOT a zero-count census.  The zero sets Z_p are only the
input columns.  The measured objects are the pointwise load

    K_P(m) = #{P < p <= 2P : b_{m mod p} == 0 (mod p)}

on the short parabolic interval 0 <= m < M = 2 L P, its exact centered
variance, and the pair discrepancies

    J_{p,q}(M) - A_p A_q / M,

where A_p is the number of marked m in the interval and J_{p,q}(M) is the
number marked simultaneously by p and q.  Since M << P^2 in the intended
regime L=P^{o(1)}, J_{p,q} is exactly an "early CRT representative" count.

Run, for example:

    sage problems/3.2/research/scripts/q7697_pointwise_rainbow_correlation.sage \
        --P 1000 --L 8 --output /tmp/q7697-P1000-L8.json

All combinatorial counts and moments are exact.  Floating point is used only
for human-readable normalized diagnostics.
"""

from sage.all import Integer, QQ, RR, Integers, prime_range, prod
import argparse
import json
from pathlib import Path


def apery_zero_set(prime):
    """Return Z_p={0<=r<p : b_r=0 mod p} by a division-free recurrence.

    For A_n=(n!)^3 b_n and n<p,

        A_{n+1}=P(n) A_n - n^6 A_{n-1},
        P(n)=34 n^3+51 n^2+27 n+5,

    and multiplication by (n!)^3 is a p-unit.  Hence A_n and b_n have the
    same vanishing set for 0<=n<p.
    """
    p = Integer(prime)
    if p < 5 or not p.is_prime():
        raise ValueError("prime must be a prime >= 5")
    F = Integers(p)
    previous = F(1)   # A_0
    current = F(5)    # A_1
    zeros = []
    if previous == 0:
        zeros.append(0)
    if current == 0:
        zeros.append(1)
    for n0 in range(1, int(p) - 1):
        n = F(n0)
        coefficient = 34*n**3 + 51*n**2 + 27*n + 5
        following = coefficient*current - n**6*previous
        previous, current = current, following
        if current == 0:
            zeros.append(n0 + 1)
    return tuple(zeros)


def validate_zero_set(p, zeros):
    """Check the two exact local structural invariants used in the campaign."""
    S = set(zeros)
    if tuple(sorted(S)) != tuple(zeros):
        raise AssertionError("zero set is not sorted/distinct at p=%s" % p)
    for r in zeros:
        if p - 1 - r not in S:
            raise AssertionError("reflection failed at p=%s, r=%s" % (p, r))
        if (r + 1) % p in S:
            raise AssertionError("consecutive Apéry zeros at p=%s, r=%s" % (p, r))


def representatives_below(M, p, zeros):
    """All m<M whose residue mod p is in Z_p, in sorted order."""
    out = []
    for r in zeros:
        if r >= M:
            continue
        out.extend(range(int(r), int(M), int(p)))
    out.sort()
    return tuple(out)


def build_state(P, L):
    P = Integer(P)
    L = Integer(L)
    if P < 5 or L < 1:
        raise ValueError("require P>=5 and L>=1")
    M = Integer(2) * L * P
    primes = tuple(Integer(p) for p in prime_range(P + 1, 2*P + 1))
    zeros = {}
    positions = {}
    A = {}
    K = [0] * int(M)
    for p in primes:
        zp = apery_zero_set(p)
        validate_zero_set(p, zp)
        pos = representatives_below(M, p, zp)
        zeros[p] = zp
        positions[p] = pos
        A[p] = Integer(len(pos))
        for m in pos:
            K[m] += 1
    return {
        "P": P,
        "L": L,
        "M": M,
        "primes": primes,
        "zeros": zeros,
        "positions": positions,
        "A": A,
        "K": K,
    }


def intersection_size_sorted(left, right):
    """Exact size of the intersection of two sorted integer tuples."""
    i = j = total = 0
    while i < len(left) and j < len(right):
        a = left[i]
        b = right[j]
        if a == b:
            total += 1
            i += 1
            j += 1
        elif a < b:
            i += 1
        else:
            j += 1
    return Integer(total)


def qjson(value):
    value = QQ(value)
    return {"num": str(value.numerator()), "den": str(value.denominator())}


def summarize_state(state, top_pairs=25, top_points=25):
    P = state["P"]
    L = state["L"]
    M = state["M"]
    primes = state["primes"]
    zeros = state["zeros"]
    positions = state["positions"]
    A = state["A"]
    K = state["K"]

    total_marks = Integer(sum(K))
    mu = QQ(total_marks, M)
    lam = sum((QQ(len(zeros[p]), p) for p in primes), QQ(0))
    variance = sum(((QQ(k) - mu)**2 for k in K), QQ(0))
    centered_abs3 = sum((abs(QQ(k) - mu)**3 for k in K), QQ(0))
    f2 = Integer(sum(k*(k-1) for k in K))
    f3 = Integer(sum(k*(k-1)*(k-2) for k in K))

    total_A = Integer(sum(A.values()))
    if total_A != total_marks:
        raise AssertionError("column totals and K totals disagree")
    expected_unordered_pairs = QQ(
        total_A**2 - sum(a*a for a in A.values()), 2*M
    )
    actual_unordered_pairs = QQ(f2, 2)
    aggregate_pair_discrepancy = actual_unordered_pairs - expected_unordered_pairs

    pair_rows = []
    variance_rhs = sum((QQ(A[p]) - QQ(A[p]*A[p], M) for p in primes), QQ(0))
    for i, p in enumerate(primes):
        for q in primes[i+1:]:
            J = intersection_size_sorted(positions[p], positions[q])
            delta = QQ(J) - QQ(A[p]*A[q], M)
            variance_rhs += 2*delta
            pair_rows.append((delta, p, q, J))
    if variance_rhs != variance:
        raise AssertionError("exact variance expansion failed")

    pair_rows.sort(key=lambda row: (row[0], row[1], row[2]), reverse=True)
    pair_preview = []
    for delta, p, q, J in pair_rows[:int(top_pairs)]:
        pair_preview.append({
            "p": int(p),
            "q": int(q),
            "J": int(J),
            "A_p": int(A[p]),
            "A_q": int(A[q]),
            "delta": qjson(delta),
        })

    max_K = max(K) if K else 0
    max_points = [m for m, k in enumerate(K) if k == max_K]
    point_preview = []
    for m in max_points[:int(top_points)]:
        incidences = []
        selected_primes = []
        for p in primes:
            r = Integer(m) % p
            if r in zeros[p]:
                selected_primes.append(p)
                incidences.append({
                    "p": int(p),
                    "quotient": int(Integer(m)//p),
                    "residue": int(r),
                })
        if len(incidences) != K[m]:
            raise AssertionError("incidence reconstruction failed")
        point_preview.append({
            "m": int(m),
            "K": int(K[m]),
            "prime_product": str(prod(selected_primes, z=Integer(1))),
            "incidences": incidences,
        })

    diag = {}
    if lam != 0:
        poisson_scale = QQ(M) * lam
        diag["variance_over_M_lambda"] = float(RR(variance / poisson_scale))
        diag["abs3_over_M_lambda"] = float(RR(centered_abs3 / poisson_scale))
        # The theorem threshold is V << M*lambda*P^(1/3-eta).  We report the
        # eta=0 boundary normalization only as a real diagnostic.
        diag["variance_over_boundary_beta_1_3"] = float(
            RR(variance / poisson_scale) / (RR(P)**(RR(1)/3))
        )
    else:
        diag["variance_over_M_lambda"] = None
        diag["abs3_over_M_lambda"] = None
        diag["variance_over_boundary_beta_1_3"] = None

    return {
        "schema": 1,
        "P": int(P),
        "L": int(L),
        "M": int(M),
        "prime_count": len(primes),
        "active_prime_count": sum(1 for p in primes if zeros[p]),
        "lambda": qjson(lam),
        "finite_mean": qjson(mu),
        "finite_mean_minus_lambda": qjson(mu-lam),
        "total_marks": int(total_marks),
        "variance": qjson(variance),
        "centered_abs3": qjson(centered_abs3),
        "F2_ordered": int(f2),
        "F3_ordered": int(f3),
        "actual_unordered_pairs": qjson(actual_unordered_pairs),
        "finite_margin_expected_unordered_pairs": qjson(expected_unordered_pairs),
        "aggregate_pair_discrepancy": qjson(aggregate_pair_discrepancy),
        "max_K": int(max_K),
        "max_point_count": len(max_points),
        "max_points_preview": point_preview,
        "largest_positive_pair_discrepancies": pair_preview,
        "diagnostics": diag,
        "zero_set_sizes": {str(int(p)): len(zeros[p]) for p in primes},
    }


def run_case(P, L, top_pairs=25, top_points=25):
    return summarize_state(
        build_state(P, L), top_pairs=top_pairs, top_points=top_points
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--P", type=int, required=True,
                        help="dyadic prime scale: P<p<=2P")
    parser.add_argument("--L", type=int, required=True,
                        help="short quotient parameter; M=2LP")
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--top-pairs", type=int, default=25)
    parser.add_argument("--top-points", type=int, default=25)
    args = parser.parse_args()
    result = run_case(args.P, args.L, args.top_pairs, args.top_points)
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(text, end="")
    else:
        args.output.write_text(text, encoding="utf-8")
        print(args.output)


if __name__ == "__main__":
    main()
