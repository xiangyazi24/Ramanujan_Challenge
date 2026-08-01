#!/usr/bin/env python3
"""Exact stress test for the proposed two-layer Mellin inverse theorem.

The exact column compares algebraic integers in Q(zeta_(p-1)).  It never
evaluates a root of unity numerically: a cyclotomic polynomial remainder is
the canonical key.  The mod-p column independently reduces the same sums
through zeta_(p-1) -> g in F_p and prints centered representatives.

The Apéry compatible-system trace is not explicitly available in the archive.
Accordingly ``AperyLift`` is the documented centered integral lift of the
Apéry Hasse--Witt polynomial.  It is an exact integer-valued test vector, but
its characteristic-zero matches must not be presented as Frobenius-trace
matches of the missing compatible system.  The other core families are
literal elliptic or symmetric-square Frobenius trace functions.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import ceil, comb, gcd, isqrt
from random import Random


PRIMES = (29, 37, 41, 53, 61, 73, 89, 101)
RANDOM_FAMILY_COUNT = 8
POSITIVE_PROPORTION = 1 / 8


def legendre(value: int, prime: int) -> int:
    value %= prime
    if value == 0:
        return 0
    return 1 if pow(value, (prime - 1) // 2, prime) == 1 else -1


def centered(value: int, prime: int) -> int:
    value %= prime
    return value if value <= prime // 2 else value - prime


def prime_factors(value: int) -> list[int]:
    factors: list[int] = []
    divisor = 2
    while divisor * divisor <= value:
        if value % divisor == 0:
            factors.append(divisor)
            while value % divisor == 0:
                value //= divisor
        divisor += 1
    if value > 1:
        factors.append(value)
    return factors


def primitive_root(prime: int) -> int:
    order = prime - 1
    factors = prime_factors(order)
    for candidate in range(2, prime):
        if all(pow(candidate, order // factor, prime) != 1 for factor in factors):
            return candidate
    raise AssertionError(f"no primitive root modulo {prime}")


def divisors(value: int) -> list[int]:
    return [divisor for divisor in range(1, value + 1) if value % divisor == 0]


def poly_div_exact(dividend: list[int], divisor: list[int]) -> list[int]:
    """Divide ascending-coefficient integer polynomials by a monic divisor."""
    assert divisor and divisor[-1] == 1
    remainder = dividend[:]
    quotient = [0] * max(1, len(dividend) - len(divisor) + 1)
    for top in range(len(remainder) - 1, len(divisor) - 2, -1):
        coefficient = remainder[top]
        shift = top - len(divisor) + 1
        quotient[shift] = coefficient
        if coefficient:
            for index, value in enumerate(divisor):
                remainder[index + shift] -= coefficient * value
    assert all(value == 0 for value in remainder[: len(divisor) - 1])
    while len(quotient) > 1 and quotient[-1] == 0:
        quotient.pop()
    return quotient


_CYCLOTOMIC_CACHE: dict[int, list[int]] = {1: [-1, 1]}


def cyclotomic(order: int) -> list[int]:
    if order not in _CYCLOTOMIC_CACHE:
        polynomial = [-1] + [0] * (order - 1) + [1]
        for divisor in divisors(order)[:-1]:
            polynomial = poly_div_exact(polynomial, cyclotomic(divisor))
        _CYCLOTOMIC_CACHE[order] = polynomial
    return _CYCLOTOMIC_CACHE[order]


def reduce_cyclotomic(polynomial: list[int], modulus: list[int]) -> tuple[int, ...]:
    """Canonical remainder modulo a monic cyclotomic polynomial."""
    remainder = polynomial[:] + [0] * max(0, len(modulus) - len(polynomial))
    for top in range(len(remainder) - 1, len(modulus) - 2, -1):
        coefficient = remainder[top]
        if coefficient:
            shift = top - len(modulus) + 1
            for index, value in enumerate(modulus):
                remainder[index + shift] -= coefficient * value
    return tuple(remainder[: len(modulus) - 1])


def apery_coefficients(prime: int) -> list[int]:
    return [
        sum(
            comb(index, k) ** 2 * comb(index + k, k) ** 2
            for k in range(index + 1)
        )
        % prime
        for index in range(prime)
    ]


def evaluate_polynomial(coefficients: list[int], argument: int, prime: int) -> int:
    answer = 0
    for coefficient in reversed(coefficients):
        answer = (answer * argument + coefficient) % prime
    return answer


def branch_polynomial(value: int, prime: int) -> int:
    return (value * value - 34 * value + 1) % prime


def elliptic_trace_general(a: int, b: int, prime: int) -> int:
    """Trace of y^2=x^3+a*x+b; caller excludes discriminant-zero fibres."""
    return -sum(legendre(x**3 + a * x + b, prime) for x in range(prime))


def legendre_trace(parameter: int, prime: int) -> int:
    if parameter == 1:
        return 0  # extension by zero at the singular parameter
    return -sum(
        legendre(x * (x - 1) * (x - parameter), prime) for x in range(prime)
    )


def franel_elliptic_trace(parameter: int, prime: int) -> int:
    singular = {0, prime - 1, pow(8, -1, prime)}
    if parameter in singular:
        return 0
    a1 = (1 - 2 * parameter) % prime
    a3 = parameter * parameter % prime
    points = 1
    for x in range(prime):
        linear_y = (a1 * x + a3) % prime
        points += 1 + legendre(linear_y * linear_y + 4 * x**3, prime)
    return prime + 1 - points


def hesse_trace(parameter: int, prime: int) -> int:
    if pow(parameter, 3, prime) == 1:
        return 0
    affine = 0
    for x in range(prime):
        for y in range(prime):
            affine += (x**3 + y**3 + 1 - 3 * parameter * x * y) % prime == 0
    infinity = sum((x**3 + 1) % prime == 0 for x in range(prime))
    return prime + 1 - affine - infinity


@dataclass(frozen=True)
class ZooObject:
    name: str
    family: str
    rank: int
    values: tuple[int, ...]  # indexed by t=0,...,p-1; only G_m is transformed
    shift: int = 0  # Kummer twist: Mellin index k is sent to k+shift
    source: str = ""


def character_label(shift: int, order: int) -> str:
    character_order = order // gcd(shift, order)
    if character_order == 1:
        return "1"
    phase = shift // (order // character_order)
    return f"K{character_order}.{phase}"


def small_kummer_shifts(order: int) -> list[int]:
    return [
        shift
        for shift in range(order)
        if order // gcd(shift, order) <= 6
    ]


def build_core_zoo(prime: int) -> list[ZooObject]:
    order = prime - 1
    coefficients = apery_coefficients(prime)
    apery = [0] + [
        centered(evaluate_polynomial(coefficients, value, prime), prime)
        for value in range(1, prime)
    ]
    apery_q = [0] + [
        legendre(branch_polynomial(value, prime), prime) * apery[value]
        for value in range(1, prime)
    ]

    objects: list[ZooObject] = []
    for family, prefix, values in (
        ("apery", "A", apery),
        ("apery_graph", "Aq", apery_q),
    ):
        for shift in small_kummer_shifts(order):
            objects.append(
                ZooObject(
                    name=f"{prefix}*{character_label(shift, order)}",
                    family=family,
                    rank=3,
                    values=tuple(values),
                    shift=shift,
                    source="centered Apéry Hasse--Witt lift",
                )
            )

    legendre_values = tuple(
        [0] + [legendre_trace(value, prime) for value in range(1, prime)]
    )
    legendre_sym2 = tuple(
        0 if value == 1 else legendre_values[value] ** 2 - prime
        for value in range(prime)
    )
    franel_values = tuple(
        [0] + [franel_elliptic_trace(value, prime) for value in range(1, prime)]
    )
    franel_sym2 = tuple(
        0
        if value in {0, prime - 1, pow(8, -1, prime)}
        else franel_values[value] ** 2 - prime
        for value in range(prime)
    )
    franel_sym2_q = tuple(
        legendre(branch_polynomial(value, prime), prime) * franel_sym2[value]
        for value in range(prime)
    )
    hesse_values = tuple(
        [0] + [hesse_trace(value, prime) for value in range(1, prime)]
    )
    objects.extend(
        [
            ZooObject("Leg2", "legendre", 2, legendre_values, source="2F1(1/2,1/2)"),
            ZooObject("LegS2", "legendre_sym2", 3, legendre_sym2, source="Sym^2 Legendre"),
            ZooObject("Fra2", "franel", 2, franel_values, source="Franel/Beauville IV"),
            ZooObject("FraS2", "franel_sym2", 3, franel_sym2, source="Sym^2 Franel"),
            ZooObject("FraS2q", "franel_sym2_graph", 3, franel_sym2_q, source="q(t)-graph twist of Sym^2 Franel"),
            ZooObject("Hes2", "hesse", 2, hesse_values, source="Hesse 2F1(1/3,2/3)"),
        ]
    )
    return objects


def random_coefficients() -> list[tuple[int, int, int, int]]:
    """Fixed small linear a(t),b(t) families used only in the counterexample hunt."""
    rng = Random(3206339)
    answer: list[tuple[int, int, int, int]] = []
    while len(answer) < RANDOM_FAMILY_COUNT:
        candidate = tuple(rng.randint(-5, 5) for _ in range(4))
        if candidate not in answer and candidate[:2] != (0, 0) and candidate[2:] != (0, 0):
            answer.append(candidate)  # type: ignore[arg-type]
    return answer


RANDOM_COEFFICIENTS = random_coefficients()


def build_random_zoo(prime: int) -> list[ZooObject]:
    answer: list[ZooObject] = []
    for index, (a0, a1, b0, b1) in enumerate(RANDOM_COEFFICIENTS):
        values = [0]
        for parameter in range(1, prime):
            a = (a0 + a1 * parameter) % prime
            b = (b0 + b1 * parameter) % prime
            discriminant = (-16 * (4 * a**3 + 27 * b**2)) % prime
            values.append(0 if discriminant == 0 else elliptic_trace_general(a, b, prime))
        answer.append(
            ZooObject(
                f"Rnd{index}",
                f"random_{index}",
                2,
                tuple(values),
                source=f"y^2=x^3+({a0}+{a1}t)x+({b0}+{b1}t)",
            )
        )
    return answer


def exact_mellin_keys(objects: list[ZooObject], prime: int, generator: int) -> list[list[tuple[int, ...]]]:
    order = prime - 1
    modulus = cyclotomic(order)
    powers = [pow(generator, exponent, prime) for exponent in range(order)]
    all_keys: list[list[tuple[int, ...]]] = []
    for obj in objects:
        keys: list[tuple[int, ...]] = []
        for character in range(order):
            exponent = (character + obj.shift) % order
            polynomial = [0] * order
            for log_t, parameter in enumerate(powers):
                polynomial[(-exponent * log_t) % order] += obj.values[parameter]
            keys.append(reduce_cyclotomic(polynomial, modulus))
        all_keys.append(keys)
    return all_keys


def mod_mellin_keys(objects: list[ZooObject], prime: int, generator: int) -> list[list[int]]:
    order = prime - 1
    powers = [pow(generator, exponent, prime) for exponent in range(order)]
    answer: list[list[int]] = []
    for obj in objects:
        row: list[int] = []
        for character in range(order):
            exponent = (character + obj.shift) % order
            value = sum(
                obj.values[parameter] * pow(parameter, -exponent, prime)
                for parameter in powers
            )
            row.append(centered(value, prime))
        answer.append(row)
    return answer


def match_matrix(keys: list[list[object]]) -> list[list[int]]:
    size = len(keys)
    return [
        [sum(left == right for left, right in zip(keys[i], keys[j])) for j in range(size)]
        for i in range(size)
    ]


@dataclass
class PrimeResult:
    prime: int
    objects: list[ZooObject]
    exact: list[list[int]]
    mod_p: list[list[int]]
    random_objects: list[ZooObject]
    random_exact: list[list[int]]
    random_mod_p: list[list[int]]


def run_prime(prime: int) -> PrimeResult:
    generator = primitive_root(prime)
    objects = build_core_zoo(prime)
    random_objects = build_random_zoo(prime)
    expanded = objects + random_objects
    expanded_exact = match_matrix(exact_mellin_keys(expanded, prime, generator))
    expanded_mod = match_matrix(mod_mellin_keys(expanded, prime, generator))
    core_size = len(objects)
    exact = [row[:core_size] for row in expanded_exact[:core_size]]
    mod_p = [row[:core_size] for row in expanded_mod[:core_size]]
    return PrimeResult(
        prime,
        objects,
        exact,
        mod_p,
        random_objects,
        expanded_exact,
        expanded_mod,
    )


def print_matrix(result: PrimeResult) -> None:
    print(f"\n=== p={result.prime}; entries=exact/mod-p ===")
    print("LEGEND " + " ".join(f"{index:02d}={obj.name}" for index, obj in enumerate(result.objects)))
    header = "    " + " ".join(f"{index:>5}" for index in range(len(result.objects)))
    print(header)
    for i in range(len(result.objects)):
        cells = " ".join(
            f"{result.exact[i][j]}/{result.mod_p[i][j]}".rjust(5)
            for j in range(len(result.objects))
        )
        print(f"{i:02d}  {cells}")


def base_index(result: PrimeResult, name: str) -> int:
    return next(index for index, obj in enumerate(result.objects) if obj.name == name)


def pair_sequence(results: list[PrimeResult], left: str, right: str) -> tuple[list[int], list[int]]:
    exact: list[int] = []
    mod_p: list[int] = []
    for result in results:
        i = base_index(result, left)
        j = base_index(result, right)
        exact.append(result.exact[i][j])
        mod_p.append(result.mod_p[i][j])
    return exact, mod_p


def counterexample_hunt(
    results: list[PrimeResult],
) -> tuple[
    list[tuple[str, str, list[int], list[int]]],
    list[tuple[str, str, list[int], list[int]]],
    list[tuple[str, str, list[int], list[int]]],
]:
    """Return pre-registered positive-proportion flags and the top ten pairs."""
    expanded_by_prime = [result.objects + result.random_objects for result in results]
    common_names = set(obj.name for obj in expanded_by_prime[0])
    for expanded in expanded_by_prime[1:]:
        common_names &= {obj.name for obj in expanded}
    names = [obj.name for obj in expanded_by_prime[0] if obj.name in common_names]
    records: list[tuple[str, str, list[int], list[int]]] = []
    for left in range(len(names)):
        for right in range(left + 1, len(names)):
            exact_counts: list[int] = []
            mod_counts: list[int] = []
            unrelated = True
            for result, expanded in zip(results, expanded_by_prime):
                index = {obj.name: position for position, obj in enumerate(expanded)}
                i, j = index[names[left]], index[names[right]]
                expanded = result.objects + result.random_objects
                obj_left, obj_right = expanded[i], expanded[j]
                if obj_left.family == obj_right.family:
                    unrelated = False
                if {obj_left.family, obj_right.family} == {"apery", "apery_graph"}:
                    unrelated = False
                if {obj_left.family, obj_right.family} == {"franel_sym2", "franel_sym2_graph"}:
                    unrelated = False
                if {obj_left.family, obj_right.family} <= {"legendre", "legendre_sym2"}:
                    unrelated = False
                if {obj_left.family, obj_right.family} <= {
                    "franel",
                    "franel_sym2",
                    "franel_sym2_graph",
                }:
                    unrelated = False
                exact_counts.append(result.random_exact[i][j])
                mod_counts.append(result.random_mod_p[i][j])
            if unrelated:
                records.append((names[left], names[right], exact_counts, mod_counts))

    exact_flags = [
        record
        for record in records
        if all(
            record[2][index] >= ceil(POSITIVE_PROPORTION * (results[index].prime - 1))
            for index in range(len(results) - 3, len(results))
        )
    ]
    mod_flags = [
        record
        for record in records
        if all(
            record[3][index] >= ceil(POSITIVE_PROPORTION * (results[index].prime - 1))
            for index in range(len(results) - 3, len(results))
        )
    ]
    top = sorted(records, key=lambda record: (max(record[2]), sum(record[2])), reverse=True)[:10]
    return exact_flags, mod_flags, top


def print_summary(results: list[PrimeResult]) -> None:
    exact_graph, mod_graph = pair_sequence(results, "A*1", "Aq*1")
    exact_literal_graph, mod_literal_graph = pair_sequence(results, "FraS2", "FraS2q")
    primes = [result.prime for result in results]
    print("\nSUMMARY")
    print(f"primes={primes}")
    print(f"core_zoo_sizes={[len(result.objects) for result in results]}")
    print(f"AperyLift/AperyQ exact={exact_graph}")
    print(f"AperyLift/AperyQ mod-p={mod_graph}")
    print(f"literal-rank3 graph exact={exact_literal_graph}")
    print(f"literal-rank3 graph mod-p={mod_literal_graph}")
    graph_positive = all(
        count >= ceil(POSITIVE_PROPORTION * (result.prime - 1))
        for count, result in zip(exact_graph[-3:], results[-3:])
    )
    print(f"Apery graph positive-proportion(>=1/8 on last 3 primes)={graph_positive}")

    distinct_kummer: list[tuple[int, int, int]] = []
    for result in results:
        counts = []
        for i, left in enumerate(result.objects):
            for j in range(i + 1, len(result.objects)):
                right = result.objects[j]
                if left.family == right.family and left.family in {"apery", "apery_graph"}:
                    counts.append(result.exact[i][j])
        distinct_kummer.append((result.prime, min(counts), max(counts)))
    print(f"distinct Kummer-related exact match ranges(p,min,max)={distinct_kummer}")
    print("same-index Kummer-twist implication=REFUTED" if any(high < (p - 1) / 8 for p, _, high in distinct_kummer) else "same-index Kummer-twist implication=NOT REFUTED")

    exact_flags, mod_flags, top = counterexample_hunt(results)
    print(
        "counterexample criterion=unrelated exact matches >= ceil((p-1)/8) "
        "for each of p=73,89,101"
    )
    print(f"exact_counterexample_flags={len(exact_flags)}")
    for left, right, exact, mod_p in exact_flags:
        print(f"EXACT FLAG {left} vs {right}: exact={exact}; mod-p={mod_p}")
    print(f"mod-p_positive-proportion_flags={len(mod_flags)}")
    for left, right, exact, mod_p in mod_flags:
        print(f"MOD-P FLAG {left} vs {right}: exact={exact}; mod-p={mod_p}")
    print("top_unrelated_pairs_by_exact_matches:")
    for left, right, exact, mod_p in top:
        print(f"  {left} vs {right}: exact={exact}; mod-p={mod_p}")
    print(
        f"random_rank2_families={RANDOM_FAMILY_COUNT}; coefficients={RANDOM_COEFFICIENTS}"
    )

    # Internal arithmetic and matrix sanity gates.
    for result in results:
        order = result.prime - 1
        assert all(result.exact[i][i] == order for i in range(len(result.objects)))
        assert all(result.mod_p[i][i] == order for i in range(len(result.objects)))
        assert result.exact == [list(row) for row in zip(*result.exact)]
        assert result.mod_p == [list(row) for row in zip(*result.mod_p)]
        assert len(cyclotomic(order)) - 1 > 0
        assert reduce_cyclotomic(cyclotomic(order), cyclotomic(order)) == (
            0,
        ) * (len(cyclotomic(order)) - 1)
    print("SANITY=PASS (cyclotomic exactness, centered mod-p separation, symmetric matrices)")


def main() -> None:
    results = [run_prime(prime) for prime in PRIMES]
    for result in results:
        print_matrix(result)
    print_summary(results)


if __name__ == "__main__":
    main()
