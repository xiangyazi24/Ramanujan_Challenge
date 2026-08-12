#!/usr/bin/env python3
"""Q7714: exact audit of quotient-state companion determinants.

This is a finite normalization/counterexample check, not an asymptotic claim.
It uses the first all-distinct-quotient K=3 row in the exact X=256 HM3 census:

    m=12678,
    (p,r,q)=(379,171,33),(443,274,28),(499,203,25).

For the actual Apery pair (a_n,b_n), it verifies the exact quotient Casoratian

  a_i b_j-a_j b_i = -6*N_{j-i}(i)/prod_{t=i+1}^j t^3,

checks p-adic denominator units at every hit prime, and tests whether any hit
prime is forced into any of the three quotient-gap numerators.
"""

from fractions import Fraction
from math import gcd


def P(n: int) -> int:
    return 34*n**3 + 51*n**2 + 27*n + 5


def apery_b(limit: int) -> list[int]:
    b = [1, 5]
    for n in range(1, limit):
        num = P(n)*b[n] - n**3*b[n-1]
        den = (n+1)**3
        assert num % den == 0
        b.append(num // den)
    return b


def apery_a(limit: int) -> list[Fraction]:
    a = [Fraction(0), Fraction(6)]
    for n in range(1, limit):
        a.append((P(n)*a[n] - n**3*a[n-1]) / (n+1)**3)
    return a


def gap_numerator(i: int, h: int) -> int:
    assert h >= 1
    if h == 1:
        return 1
    n0 = 1
    n1 = P(i+1)
    if h == 2:
        return n1
    for k in range(2, h):
        n0, n1 = n1, P(i+k)*n1 - (i+k)**6*n0
    return n1


def cleared_gap_denominator(i: int, j: int) -> int:
    out = 1
    for t in range(i+1, j+1):
        out *= t**3
    return out


def main() -> None:
    m = 12678
    hits = [(379,171,33), (443,274,28), (499,203,25)]
    b = apery_b(max(r for _,r,_ in hits))
    a = apery_a(max(q for _,_,q in hits))
    qs = sorted(q for _,_,q in hits)

    print(f"m={m}")
    print("hits=" + repr(hits))

    for p,r,q in hits:
        assert m == q*p+r
        assert b[r] % p == 0
        assert b[q] % p != 0
        assert gcd(a[q].denominator, p) == 1
        cross_b = [b[t] % p for t in qs]
        cross_a_den = [a[t].denominator % p for t in qs]
        assert all(cross_b)
        assert all(cross_a_den)
        print(
            f"hit p={p} r={r} q={q} "
            f"bq_mod_p={b[q] % p} aq_den_mod_p={a[q].denominator % p} "
            f"all_q_b_mod_p={cross_b} all_q_a_den_mod_p={cross_a_den}"
        )

    pair_data = []
    for x in range(len(qs)):
        for y in range(x+1, len(qs)):
            i,j = qs[x],qs[y]
            h = j-i
            N = gap_numerator(i,h)
            D = cleared_gap_denominator(i,j)
            lhs = a[i]*b[j] - a[j]*b[i]
            rhs = Fraction(-6*N,D)
            assert lhs == rhs
            pair_data.append((i,j,N,lhs))
            print(
                f"pair q=({i},{j}) gap={h} N={N} "
                f"N_digits={len(str(abs(N)))} "
                f"det_num={lhs.numerator} det_den={lhs.denominator}"
            )

    product_N = 1
    for _,_,N,_ in pair_data:
        product_N *= N

    for p,r,q in hits:
        residues = [N % p for _,_,N,_ in pair_data]
        det_den_units = [gcd(det.denominator,p) == 1 for *_,det in pair_data]
        assert all(det_den_units)
        assert all(x != 0 for x in residues), (p, residues)
        assert product_N % p != 0
        print(f"prime p={p}: gap_N_residues={residues}; product_is_unit=True")

    # The quotient-state projective ratios are distinct over Q and remain
    # pairwise distinct modulo every one of the three hit primes.
    slopes = [a[q] / b[q] for q in qs]
    assert len(set(slopes)) == len(slopes)
    for p,_,_ in hits:
        slope_mod = [
            (s.numerator % p) * pow(s.denominator % p, -1, p) % p
            for s in slopes
        ]
        assert len(set(slope_mod)) == len(slope_mod)
        print(f"prime p={p}: quotient_slopes_mod_p={slope_mod}; pairwise_distinct=True")

    for i in range(1,len(qs)):
        q0,q1 = qs[i-1],qs[i]
        exact_increment = slopes[i]-slopes[i-1]
        expected = sum(
            (Fraction(6, k**3 * b[k-1] * b[k])
             for k in range(q0+1,q1+1)),
            Fraction(0),
        )
        assert exact_increment == expected
    print("projective_slope_telescoping=PASS")
    print("quotient_casoratian_identity=PASS")
    print("all_three_hit_primes_are_units_on_all_pair_gap_numerators=PASS")
    print("all_three_projective_quotient_sets_are_separable_mod_each_hit_prime=PASS")
    print("finite_counterexample_only=True")


if __name__ == "__main__":
    main()
