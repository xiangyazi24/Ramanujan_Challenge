#!/usr/bin/env python3
"""E3: center parametrization for the Apéry-number recurrence.

The script computes the fixed rational center basis A_i, B_i, studies

    W(i, j) = A_i B_j - A_j B_i,

and verifies the parametrization modulo selected primes.  Its default output
is /tmp/e3_center_param_results.txt.

Integer factorization note: generic W-numerators in the requested range have
up to 235 decimal digits.  Complete prime factorization of arbitrary integers
of that size is not computationally realistic.  We therefore remove every
prime up to --factor-bound, test the exact remaining cofactor for primality,
and print any unresolved composite cofactor in full.  Thus every printed
decomposition is exact, even when it is only a partial prime factorization.
"""

from __future__ import annotations

import argparse
import os
from fractions import Fraction
from itertools import combinations
from math import gcd, isqrt, prod
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

try:
    from sympy import isprime as sympy_isprime
    from sympy import Poly, QQ, expand, factor_list, factorint, symbols
except ImportError:  # pragma: no cover - SymPy is available in the project env.
    sympy_isprime = None
    Poly = QQ = expand = factor_list = factorint = symbols = None


Q = Fraction


def P(n):
    """P(n) = 34 n^3 + 51 n^2 + 27 n + 5 over Z or Q."""
    return 34 * n**3 + 51 * n**2 + 27 * n + 5


def center_basis_exact(max_i: int) -> Tuple[List[Q], List[Q]]:
    """Return A_i, B_i for 0 <= i <= max_i using exact fractions."""
    if max_i < 0:
        raise ValueError("max_i must be nonnegative")
    if max_i == 0:
        return [Q(1)], [Q(0)]

    A = [Q(1), Q(0)]
    B = [Q(0), Q(1)]
    for i in range(1, max_i):
        denominator = (2 * i + 1) ** 3
        middle = 4 * i * (68 * i * i + 3)
        previous = (2 * i - 1) ** 3
        A.append((middle * A[i] - previous * A[i - 1]) / denominator)
        B.append((middle * B[i] - previous * B[i - 1]) / denominator)
    return A, B


def W_exact(A: Sequence[Q], B: Sequence[Q], i: int, j: int) -> Q:
    """Casorati determinant, with the exact even extension to signed indices."""
    ai, aj = abs(i), abs(j)
    return A[ai] * B[aj] - A[aj] * B[ai]


def p_weight(k: int):
    """p_k=(2k-1)^3, also usable with a symbolic k."""
    return (2 * k - 1) ** 3


def q_weight(k: int):
    """q_k=4k(68k^2+3), also usable with a symbolic k."""
    return 4 * k * (68 * k * k + 3)


def continuant(a: int, b: int) -> int:
    """D[a,b] for diag q_k and adjacent off-diagonal p_{k+1}."""
    if b == a - 2:
        return 0
    if b == a - 1:
        return 1
    if b < a - 2:
        raise ValueError("invalid continuant interval")
    previous_two, previous_one = 0, 1
    for k in range(a, b + 1):
        current = q_weight(k) * previous_one - p_weight(k) ** 2 * previous_two
        previous_two, previous_one = previous_one, current
    return previous_one


def H(n: int) -> int:
    """H_n=prod_{k=2}^n p_k, with H_0=H_1=1."""
    result = 1
    for k in range(2, n + 1):
        result *= p_weight(k)
    return result


def interval_denominator(i: int, j: int) -> int:
    """Natural denominator prod_{k=i+1}^j p_k for 0<=i<j."""
    result = 1
    for k in range(i + 1, j + 1):
        result *= p_weight(k)
    return result


def cleared_W(i: int, j: int) -> int:
    """The raw integer continuant C(i,j)=D[i+1,j-1]."""
    if not 0 <= i < j:
        raise ValueError("cleared_W expects 0<=i<j")
    return continuant(i + 1, j - 1)


def symbolic_gap_scan(max_d: int = 25):
    """Factor fixed-gap continuant polynomials C_d(s) over Q[s]."""
    if symbols is None:
        return None
    s = symbols("s")
    results = []
    polynomials = {}
    for d in range(1, max_d + 1):
        i = (s - d) / 2
        previous_two, previous_one = 0, 1
        for offset in range(1, d):
            k = i + offset
            current = expand(
                q_weight(k) * previous_one - p_weight(k) ** 2 * previous_two
            )
            previous_two, previous_one = previous_one, current
        polynomial = Poly(previous_one, s, domain=QQ)
        coefficient, factors = factor_list(polynomial.as_expr(), s)
        factor_degrees = [
            (Poly(factor, s, domain=QQ).degree(), exponent)
            for factor, exponent in factors
        ]
        parity_sign = -1 if (d - 1) % 2 else 1
        reflected = Poly(polynomial.as_expr().subs(s, -s), s, domain=QQ)
        assert reflected == parity_sign * polynomial
        results.append((d, polynomial.degree(), coefficient, factor_degrees))
        polynomials[d] = polynomial
    return s, polynomials, results


def sieve_primes(limit: int) -> List[int]:
    if limit < 2:
        return []
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[0:2] = b"\x00\x00"
    for p in range(2, isqrt(limit) + 1):
        if sieve[p]:
            start = p * p
            sieve[start : limit + 1 : p] = b"\x00" * (
                (limit - start) // p + 1
            )
    return [p for p in range(2, limit + 1) if sieve[p]]


def is_prime_small(n: int) -> bool:
    if n < 2:
        return False
    for p in range(2, isqrt(n) + 1):
        if n % p == 0:
            return n == p
    return True


class PartialFactorizer:
    """Exact trial factorization with primality testing of the final cofactor."""

    def __init__(self, bound: int, full_factor_digits: int):
        if bound < 2:
            raise ValueError("factor bound must be at least 2")
        if full_factor_digits < 0:
            raise ValueError("full-factor digit threshold must be nonnegative")
        self.bound = bound
        self.full_factor_digits = full_factor_digits
        self.primes = sieve_primes(bound)
        self.cache: Dict[int, Tuple[List[Tuple[int, int]], Optional[int]]] = {}

    def factor(self, n: int) -> Tuple[List[Tuple[int, int]], Optional[int]]:
        """Return (proven/probable prime factors, unresolved composite or None)."""
        n = abs(n)
        if n in self.cache:
            return self.cache[n]
        if n in (0, 1):
            result = ([], None)
            self.cache[n] = result
            return result

        remaining = n
        factors: List[Tuple[int, int]] = []
        for p in self.primes:
            if p * p > remaining:
                factors.append((remaining, 1))
                remaining = 1
                break
            if remaining % p != 0:
                continue
            exponent = 0
            while remaining % p == 0:
                remaining //= p
                exponent += 1
            factors.append((p, exponent))
            if remaining == 1:
                break

        unresolved: Optional[int] = None
        if remaining > 1:
            # SymPy's isprime uses deterministic tests for small inputs and a
            # strong BPSW test for large inputs.  A False result is a proof of
            # compositeness; a True large result is reported below as PRP.
            if sympy_isprime is not None and bool(sympy_isprime(remaining)):
                factors.append((remaining, 1))
            elif remaining <= self.bound * self.bound:
                # All possible prime divisors were exhausted by trial division.
                factors.append((remaining, 1))
            elif factorint is not None and len(str(remaining)) <= self.full_factor_digits:
                try:
                    extra = factorint(remaining)
                    extra_ok = (
                        all(bool(sympy_isprime(prime)) for prime in extra)
                        and all(exponent > 0 for exponent in extra.values())
                        and prod(prime**exponent for prime, exponent in extra.items())
                        == remaining
                    )
                except (ArithmeticError, ValueError):
                    extra_ok = False
                    extra = {}
                if extra_ok:
                    factors.extend(sorted(extra.items()))
                else:
                    unresolved = remaining
            else:
                unresolved = remaining

        reconstructed = 1
        for prime, exponent in factors:
            reconstructed *= prime**exponent
        if unresolved is not None:
            reconstructed *= unresolved
        assert reconstructed == n

        result = (factors, unresolved)
        self.cache[n] = result
        return result

    def format(self, n: int) -> Tuple[str, bool]:
        if n == 0:
            return "0", True
        if abs(n) == 1:
            return str(n), True
        factors, unresolved = self.factor(n)
        pieces = []
        if n < 0:
            pieces.append("-1")
        for p, exponent in factors:
            atom = str(p) if exponent == 1 else f"{p}^{exponent}"
            pieces.append(atom)
        if unresolved is not None:
            pieces.append(f"C{len(str(unresolved))}({unresolved})")
        return " * ".join(pieces), unresolved is None


def apery_mod_p(p: int) -> List[int]:
    """Compute b_0,...,b_{p-1} modulo an odd prime p."""
    if p < 3 or not is_prime_small(p):
        raise ValueError(f"expected an odd prime, got {p}")
    b = [0] * p
    b[0] = 1
    b[1] = 5 % p
    for n in range(1, p - 1):
        numerator = (P(n) * b[n] - n**3 * b[n - 1]) % p
        denominator = pow(n + 1, 3, p)
        # n+1 <= p-1 here, so this inverse always exists.
        b[n + 1] = numerator * pow(denominator, -1, p) % p
    return b


def center_basis_mod_p(p: int) -> Tuple[List[int], List[int]]:
    """Compute A_i, B_i modulo p for the full center range 0 <= i <= c."""
    c = (p - 1) // 2
    if c == 0:
        return [1], [0]
    A = [1, 0]
    B = [0, 1]
    for i in range(1, c):
        denominator = pow(2 * i + 1, 3, p)
        assert denominator != 0  # 2i+1 <= p-2
        inverse = pow(denominator, -1, p)
        middle = 4 * i * (68 * i * i + 3)
        previous = (2 * i - 1) ** 3
        A.append((middle * A[i] - previous * A[i - 1]) * inverse % p)
        B.append((middle * B[i] - previous * B[i - 1]) * inverse % p)
    return A, B


def reduce_fraction_mod_p(value: Q, p: int) -> int:
    assert value.denominator % p != 0
    return value.numerator * pow(value.denominator, -1, p) % p


def fstr(value: Q) -> str:
    return str(value.numerator) if value.denominator == 1 else str(value)


class Report:
    def __init__(self, path: Path):
        self.path = path
        self.tmp_path = path.with_name(path.name + ".tmp")
        self.file = self.tmp_path.open("w", encoding="utf-8")

    def line(self, text: str = "") -> None:
        self.file.write(text + "\n")

    def finish(self) -> None:
        self.file.close()
        os.replace(self.tmp_path, self.path)


def write_derivation_and_patterns(report: Report, A: Sequence[Q], B: Sequence[Q]) -> None:
    half = Q(-1, 2)
    assert P(half) == 0
    # Coefficient check for (2x+1)(17x^2+17x+5).
    assert (2 * 17, 2 * 17 + 17, 2 * 5 + 17, 5) == (34, 51, 27, 5)

    report.line("E3: CENTER PARAMETRIZATION AND CASORATI DETERMINANT")
    report.line("=" * 72)
    report.line()
    report.line("1. Polynomial and centered recurrence")
    report.line("-" * 72)
    report.line("P(n) = 34*n^3 + 51*n^2 + 27*n + 5")
    report.line("P(n) = (2*n+1)*(17*n^2+17*n+5)")
    report.line(f"P(-1/2) = {P(half)}  [verified exactly]")
    report.line("P(-1/2-i) = -i*(68*i^2+3)/2")
    report.line()
    report.line("Writing y_i=b_{c-i}, the fixed recurrence over Q is")
    report.line("  (2*i+1)^3*y_{i+1}")
    report.line("    = 4*i*(68*i^2+3)*y_i - (2*i-1)^3*y_{i-1}.")
    report.line("The two bases have (A_0,A_1)=(1,0), (B_0,B_1)=(0,1).")
    report.line("At i=0 the recurrence gives y_{-1}=y_1; invariance under")
    report.line("i -> -i then gives A_{-i}=A_i and B_{-i}=B_i exactly.")
    report.line()

    # Exact continuant formulas, adjacent determinant, and nonvanishing check.
    for n in range(2, len(A)):
        assert A[n] == Q(-continuant(2, n - 1), H(n))
        assert B[n] == Q(continuant(1, n - 1), H(n))
    continuant_check_max = min(60, len(A) - 1)
    for i in range(continuant_check_max + 1):
        for j in range(i + 1, continuant_check_max + 1):
            raw = cleared_W(i, j)
            natural_denominator = interval_denominator(i, j)
            value = W_exact(A, B, i, j)
            assert value == Q(raw, natural_denominator)
            common = gcd(raw, natural_denominator)
            assert value.numerator == raw // common
            assert value.denominator == natural_denominator // common
    for i in range(len(A) - 1):
        assert W_exact(A, B, i, i + 1) == Q(1, (2 * i + 1) ** 3)
    for i in range(continuant_check_max - 1):
        for j in range(i + 2, continuant_check_max + 1):
            dodgson = (
                W_exact(A, B, i, j) * W_exact(A, B, i + 1, j - 1)
                - W_exact(A, B, i, j - 1) * W_exact(A, B, i + 1, j)
            )
            assert dodgson == Q(-1, p_weight(i + 1) * p_weight(j))
    checked_positive = 0
    for i in range(len(A)):
        for j in range(i + 1, len(A)):
            assert W_exact(A, B, i, j) > 0
            checked_positive += 1

    gap_scan = symbolic_gap_scan(25)
    if gap_scan is not None:
        s, gap_polynomials, gap_results = gap_scan
        assert gap_polynomials[1] == Poly(1, s, domain=QQ)
        assert gap_polynomials[2] == Poly(2 * s * (17 * s**2 + 3), s, domain=QQ)
        assert gap_polynomials[3] == Poly(
            1155 * s**6 - 3060 * s**4 + 3504 * s**2 - 1600,
            s,
            domain=QQ,
        )
        for d, degree, _, factor_degrees in gap_results:
            assert degree == 3 * (d - 1)
            if d == 1:
                assert factor_degrees == []
            elif d % 2 == 0:
                assert sorted(factor_degrees) == [(1, 1), (degree - 1, 1)]
            else:
                assert factor_degrees == [(degree, 1)]

    report.line("2. Exact Casorati structure and nonvanishing")
    report.line("-" * 72)
    report.line("Put p_k=(2*k-1)^3, q_k=4*k*(68*k^2+3), and")
    report.line("H_0=H_1=1, H_n=product(k=2..n) p_k.  Then")
    report.line("  p_(i+1)y_(i+1)=q_i*y_i-p_i*y_(i-1),")
    report.line("  q_i=p_i+p_(i+1)+256*i^3.")
    report.line()
    report.line("Let D[a,b] be the symmetric tridiagonal continuant with diagonal")
    report.line("q_a,...,q_b and adjacent off-diagonal p_(a+1),...,p_b;")
    report.line("D[a,a-1]=1 and D[a,a-2]=0.  Exact formulas are")
    report.line("  A_n=-D[2,n-1]/H_n,  B_n=D[1,n-1]/H_n  (n>=2),")
    report.line("  W(i,j)=H_i/H_j * D[i+1,j-1]")
    report.line("        =D[i+1,j-1]/product(k=i+1..j) p_k  (0<=i<j).")
    report.line("Thus the raw cleared numerator C(i,j)=D[i+1,j-1] is an integer;")
    report.line("the reduced numerator is C/gcd(C,product p_k).")
    report.line("The adjacent case is the exact telescoping product")
    report.line("  W(i,i+1)=1/(2*i+1)^3.")
    report.line()
    report.line("For fixed k, z_j=W(k,j) obeys the same recurrence and starts")
    report.line("z_k=0, z_{k+1}=1/(2*k+1)^3>0.  If")
    report.line("  a_i=4*i*(68*i^2+3)/(2*i+1)^3,")
    report.line("  r_i=(2*i-1)^3/(2*i+1)^3,")
    report.line("then a_i-r_i-1=256*i^3/(2*i+1)^3>0.  Hence")
    report.line("z_{j+1}=(a_j-r_j)z_j+r_j(z_j-z_{j-1})>z_j by induction.")
    report.line("Therefore W(i,j)>0 for every 0<=i<j; in particular it is")
    report.line("nonzero for distinct nonnegative indices.  The diagonal is zero.")
    report.line("Under signed even extension, W(i,j)=0 iff |i|=|j|, so")
    report.line("W(i,-i)=0 is the exact reflection degeneracy.")
    report.line(
        f"Direct exact check: all {checked_positive} pairs 0<=i<j<={len(A)-1} are positive."
    )
    report.line()
    report.line("The continuants obey the exact Dodgson/Pluecker identity")
    report.line("  W(i,j)W(i+1,j-1)-W(i,j-1)W(i+1,j)")
    report.line("    = -1/((2i+1)^3*(2j-1)^3).")
    report.line("When the endpoint weights are invertible modulo p, this says that")
    report.line("W(i,j) cannot vanish together with W(i,j-1), nor with W(i+1,j).")
    report.line()

    # Pure difference/sum separation test.
    pure_lhs = W_exact(A, B, 2, 3) * W_exact(A, B, 2, 5)
    pure_rhs = W_exact(A, B, 3, 4) * W_exact(A, B, 1, 4)
    pure_ratio = pure_lhs / pure_rhs
    assert pure_ratio != 1

    # A stronger test allowing arbitrary left/right endpoint factors.
    broad_lhs = (
        W_exact(A, B, 1, 5)
        * W_exact(A, B, 2, 7)
        * W_exact(A, B, 3, 4)
        * W_exact(A, B, 4, 6)
    )
    broad_rhs = (
        W_exact(A, B, 1, 6)
        * W_exact(A, B, 2, 4)
        * W_exact(A, B, 3, 7)
        * W_exact(A, B, 4, 5)
    )
    broad_ratio = broad_lhs / broad_rhs
    assert broad_ratio != 1
    shifted_broad_ratios = []
    for t in range(54):
        numerator = (
            W_exact(A, B, 1 + t, 5 + t)
            * W_exact(A, B, 2 + t, 7 + t)
            * W_exact(A, B, 3 + t, 4 + t)
            * W_exact(A, B, 4 + t, 6 + t)
        )
        denominator = (
            W_exact(A, B, 1 + t, 6 + t)
            * W_exact(A, B, 2 + t, 4 + t)
            * W_exact(A, B, 3 + t, 7 + t)
            * W_exact(A, B, 4 + t, 5 + t)
        )
        ratio = numerator / denominator
        assert ratio != 1
        shifted_broad_ratios.append(ratio)

    report.line("3. Multiplicative / EDS-like pattern tests")
    report.line("-" * 72)
    report.line("For d=j-i and s=i+j, write C_d(s)=D[i+1,j-1].  It is a")
    report.line("degree 3(d-1) polynomial and center reflection gives")
    report.line("  C_d(-s)=(-1)^(d-1) C_d(s).")
    report.line("The first exact polynomials are")
    report.line("  C_1(s)=1,")
    report.line("  C_2(s)=2*s*(17*s^2+3),")
    report.line("  C_3(s)=1155*s^6-3060*s^4+3504*s^2-1600.")
    if gap_scan is not None:
        report.line("Exact Q[s] factor scan for 2<=d<=25: every odd d>=3 case")
        report.line("is irreducible; every even d has only the forced factor s and")
        report.line("one irreducible degree 3(d-1)-1 cofactor.")
        report.line("Factor degrees by d:")
        for d, degree, _, factor_degrees in gap_results:
            report.line(f"    d={d}: degree={degree}, factors={factor_degrees}")
    else:
        report.line("SymPy unavailable: the optional exact Q[s] factor scan was skipped.")
    report.line()
    report.line("Pure separation W(i,j)=F(j-i)G(i+j) fails.  A necessary 2x2")
    report.line("minor, using (difference,sum)=(1,5),(1,7),(3,5),(3,7), gives")
    report.line("  W(2,3)W(2,5) / (W(3,4)W(1,4))")
    report.line(f"    = {fstr(pure_ratio)} != 1.")
    report.line()
    report.line("Even the most general separated endpoint ansatz")
    report.line("  W(i,j)=L_i*R_j*F(j-i)*G(i+j)")
    report.line("fails.  In the following ratio, the multisets of left endpoints,")
    report.line("right endpoints, differences, and sums all cancel:")
    report.line("  W(1,5)W(2,7)W(3,4)W(4,6)")
    report.line("  -------------------------------- =")
    report.line("  W(1,6)W(2,4)W(3,7)W(4,5)")
    report.line(f"    {fstr(broad_ratio)} != 1.")
    report.line("All 54 translations t=0,...,53 whose indices remain <=60 also fail.")
    report.line()
    report.line("The EDS ansatz psi_i*psi_j/(psi_(j-i)*psi_(i+j)), including")
    report.line("a global scale and arbitrary endpoint normalizations, is a special")
    report.line("case of the rejected L/R/F/G form.  The exact Pluecker and Dodgson")
    report.line("identities above are universal rank-two/continuant identities, not")
    report.line("evidence for an elliptic group law.")
    report.line()
    report.line("Conclusion: there is an endpoint telescoping product, a continuant,")
    report.line("and reflection-forced parity, but no difference/sum product or the")
    report.line("tested EDS law.  This supports switching to p-adic slope methods;")
    report.line("it does not logically exclude every imaginable nonstandard group law.")
    report.line("The large numerators also acquire unrelated composite cofactors;")
    report.line("no stable divisibility pattern was detected in this finite scan.")
    report.line()


def write_sequences(report: Report, A: Sequence[Q], B: Sequence[Q]) -> None:
    report.line("4. Exact center basis sequences")
    report.line("-" * 72)
    report.line(f"Entries i=0,...,{len(A)-1}; b_(c-i)=A_i*b_c+B_i*b_(c-1).")
    for i, (a, b) in enumerate(zip(A, B)):
        report.line(f"i={i}  A_i={fstr(a)}  B_i={fstr(b)}")
    report.line()


def write_w_table(
    report: Report,
    A: Sequence[Q],
    B: Sequence[Q],
    w_max: int,
    factor_bound: int,
    full_factor_digits: int,
) -> None:
    factorizer = PartialFactorizer(factor_bound, full_factor_digits)
    report.line("5. Casorati table and numerator factorizations")
    report.line("-" * 72)
    report.line(f"Range: 0<=i,j<={w_max}.")
    report.line("W(i,i)=0; W(j,i)=-W(i,j).  The table lists i<j, so it")
    report.line("determines every ordered pair in the requested square.")
    report.line(f"All prime factors <= {factor_bound} are removed exactly.")
    report.line(
        f"Residual cofactors with <= {full_factor_digits} digits are fully factored."
    )
    report.line("Ck(N) denotes an exact k-digit residual cofactor N that is known")
    report.line("composite but was not completely split.  Such a line is an exact")
    report.line("partial prime factorization, not a claim that Ck is irreducible.")
    report.line("For factors >2^64, SymPy isprime=True is a BPSW probable-prime")
    report.line("classification; such completed lines are marked complete-or-PRP.")
    report.line()

    complete = 0
    unresolved = 0
    max_digits = 0
    for i in range(w_max + 1):
        for j in range(i + 1, w_max + 1):
            value = W_exact(A, B, i, j)
            assert value > 0
            max_digits = max(max_digits, len(str(abs(value.numerator))))
            factor_text, is_complete = factorizer.format(value.numerator)
            complete += int(is_complete)
            unresolved += int(not is_complete)
            status = "complete-or-PRP" if is_complete else "partial"
            report.line(
                f"i={i} j={j}  W={fstr(value)}  "
                f"num_factors={factor_text}  status={status}"
            )

    report.line()
    report.line(
        f"Factorization summary: {complete} complete/prime-or-PRP residual, "
        f"{unresolved} with an unresolved composite cofactor; "
        f"maximum numerator length {max_digits} digits."
    )
    report.line("Complete prime factorization of generic 100--235 digit composite")
    report.line("residuals is outside a bounded local computation; the residual")
    report.line("integers are printed in full so the data are exact/reproducible.")
    report.line()


def write_prime_verification(
    report: Report, A_exact: Sequence[Q], B_exact: Sequence[Q], primes: Iterable[int]
) -> None:
    report.line("6. Finite-field verification")
    report.line("-" * 72)
    report.line("For signed i=c-j, use the exact even extension A_i=A_|i|,")
    report.line("B_i=B_|i|.  The primary check is projective:")
    report.line("  [b_c:b_(c-1)] = [-B_i:A_i],")
    report.line("equivalently A_i*b_c+B_i*b_(c-1)=0 at a zero.")
    report.line("This remains meaningful when theta=b_c/b_(c-1) and")
    report.line("rho_i=-B_i/A_i are undefined.")
    report.line("For 0<=|i|<=c, each denominator base satisfies 1<=2k-1<=p-2,")
    report.line("so its cube is nonzero mod p; W=0 mod p is exactly")
    report.line("p-divisibility of its reduced integer numerator.")
    report.line()

    for p in primes:
        b = apery_mod_p(p)
        c = (p - 1) // 2
        A, B = center_basis_mod_p(p)
        assert len(A) == c + 1 and len(B) == c + 1
        zeros = [j for j, value in enumerate(b) if value == 0]

        # Check the center basis at every index, not only at zeros.
        all_param_ok = all(
            b[j] == (A[abs(c - j)] * b[c] + B[abs(c - j)] * b[c - 1]) % p
            for j in range(p)
        )
        reflection_ok = all(b[j] == b[p - 1 - j] for j in range(p))
        assert all_param_ok and reflection_ok

        # The modular recurrence agrees with all saved exact fractions in range.
        exact_overlap = min(c, len(A_exact) - 1)
        for i in range(exact_overlap + 1):
            assert reduce_fraction_mod_p(A_exact[i], p) == A[i]
            assert reduce_fraction_mod_p(B_exact[i], p) == B[i]

        theta = None
        if b[c - 1] != 0:
            theta = b[c] * pow(b[c - 1], -1, p) % p

        report.line(f"p={p}, c={c}")
        theta_text = theta if theta is not None else "undefined"
        report.line(f"  b_c={b[c]}, b_(c-1)={b[c-1]}, theta={theta_text}")
        report.line(f"  zeros={zeros}")
        report.line(f"  signed_offsets={[c-j for j in zeros]}")
        report.line(
            f"  all-{p}-positions parametrization={all_param_ok}; reflection={reflection_ok}"
        )

        for j in zeros:
            signed_i = c - j
            i = abs(signed_i)
            linear = (A[i] * b[c] + B[i] * b[c - 1]) % p
            projective_ok = linear == 0
            assert projective_ok
            rho = None if A[i] == 0 else (-B[i] * pow(A[i], -1, p)) % p
            if theta is not None and rho is not None:
                assert theta == rho
                slope_status = f"rho={rho}=theta"
            else:
                slope_status = "rho/theta undefined; projective match"
            report.line(
                f"    zero j={j}: i={signed_i}, |i|={i}, "
                f"A_i={A[i]}, B_i={B[i]}, {slope_status}, linear={linear}"
            )

        pair_count = 0
        nontrivial_count = 0
        for j1, j2 in combinations(zeros, 2):
            i1, i2 = c - j1, c - j2
            determinant = (
                A[abs(i1)] * B[abs(i2)] - A[abs(i2)] * B[abs(i1)]
            ) % p
            assert determinant == 0
            exact_reflection_zero = abs(i1) == abs(i2)
            pair_count += 1
            nontrivial_count += int(not exact_reflection_zero)
            detail = (
                "exact |i|=|i'| degeneracy"
                if exact_reflection_zero
                else "nonzero-Q W divisible mod p"
            )

            # If both indices are in the saved exact range, directly check that
            # p divides the reduced numerator and not its denominator.
            exact_detail = ""
            if max(abs(i1), abs(i2)) < len(A_exact):
                exact_w = W_exact(A_exact, B_exact, i1, i2)
                assert exact_w.denominator % p != 0
                assert exact_w.numerator % p == 0
                exact_detail = (
                    f", numerator_mod_p={exact_w.numerator % p}, "
                    f"denominator_mod_p={exact_w.denominator % p}"
                )
            report.line(
                f"    pair offsets ({i1},{i2}): W_mod_p={determinant}; "
                f"{detail}{exact_detail}"
            )

        report.line(
            f"  pair summary: {pair_count} pairs, {nontrivial_count} with distinct |i|; all pass"
        )
        report.line()

    report.line("Interpretation of the selected-prime test:")
    report.line("The five requested primes have no pair with two distinct absolute")
    report.line("center distances; their two-zero cases are only reflection pairs.")
    report.line("The added p=181 case has distances 43 and 71 and verifies the")
    report.line("nontrivial assertion p | numerator(W(43,71)).")
    report.line("The added p=19 case exercises the projective zero-denominator case.")
    report.line()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-index", type=int, default=200)
    parser.add_argument("--w-max", type=int, default=60)
    parser.add_argument("--factor-bound", type=int, default=100_000)
    parser.add_argument("--full-factor-digits", type=int, default=30)
    parser.add_argument(
        "--output", type=Path, default=Path("/tmp/e3_center_param_results.txt")
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.max_index < 60:
        raise ValueError("--max-index must be at least 60 for the pattern tests")
    if args.w_max < 0:
        raise ValueError("--w-max must be nonnegative")
    if args.max_index < args.w_max:
        raise ValueError("--max-index must be at least --w-max")

    A, B = center_basis_exact(args.max_index)
    report = Report(args.output)
    try:
        write_derivation_and_patterns(report, A, B)
        write_sequences(report, A, B)
        write_w_table(
            report,
            A,
            B,
            args.w_max,
            args.factor_bound,
            args.full_factor_digits,
        )
        # Requested primes, followed by a genuinely nontrivial determinant test
        # and a projective zero-denominator test.
        write_prime_verification(
            report, A, B, [11, 101, 1009, 3137, 10007, 181, 19]
        )
        report.finish()
    except BaseException:
        report.file.close()
        report.tmp_path.unlink(missing_ok=True)
        raise

    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
