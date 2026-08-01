#!/usr/bin/env python3
"""Check the Franel tensor-square pushforward and its small self-twists.

The computation is deliberately exact.  Integer Frobenius traces come from
point counts on

    E_u : y^2 + (1 - 2u)xy + u^2 y = x^3,

and every finite-field identity is checked with integer modular arithmetic.
For nonsplit fibres of the quadratic cover, a tiny F_{p^2} implementation
checks the indispensable factor

    A_p(phi(x)) = H_p(x)^2 / (1 + x)^(p-1).

No floating-point character values are used.  A character is represented by
its exponent relative to the least primitive root modulo p, so the self-twist
test is the exact kernel-support criterion.
"""

from __future__ import annotations

import argparse
import hashlib
from dataclasses import dataclass
from math import comb, gcd
from typing import Iterable


PRIMES = (29, 37, 41, 53, 61, 73, 89, 101)
Pair = tuple[int, int]


def legendre(value: int, prime: int) -> int:
    value %= prime
    if value == 0:
        return 0
    return 1 if pow(value, (prime - 1) // 2, prime) == 1 else -1


def evaluate(coefficients: list[int], argument: int, prime: int) -> int:
    value = 0
    for coefficient in reversed(coefficients):
        value = (value * argument + coefficient) % prime
    return value


def franel_coefficients(prime: int) -> list[int]:
    return [
        sum(comb(index, k) ** 3 for k in range(index + 1)) % prime
        for index in range(prime)
    ]


def apery_coefficients(prime: int) -> list[int]:
    return [
        sum(
            comb(index, k) ** 2 * comb(index + k, k) ** 2
            for k in range(index + 1)
        )
        % prime
        for index in range(prime)
    ]


def centered_lift(value: int, prime: int) -> int:
    value %= prime
    return value if value <= prime // 2 else value - prime


def cover_map(parameter: int, prime: int) -> int:
    assert (1 + parameter) % prime != 0
    return (
        parameter
        * (1 - 8 * parameter)
        * pow(1 + parameter, -1, prime)
    ) % prime


def branch_polynomial(value: int, prime: int) -> int:
    return (value * value - 34 * value + 1) % prime


def elliptic_trace(parameter: int, prime: int) -> int:
    """Return p+1-#E_parameter(F_p), including the nodal fibre convention."""
    a1 = (1 - 2 * parameter) % prime
    a3 = parameter * parameter % prime
    points = 1
    for x_value in range(prime):
        linear_y = (a1 * x_value + a3) % prime
        discriminant = (linear_y * linear_y + 4 * x_value**3) % prime
        points += 1 + legendre(discriminant, prime)
    return prime + 1 - points


def qadd(left: Pair, right: Pair, prime: int) -> Pair:
    return ((left[0] + right[0]) % prime, (left[1] + right[1]) % prime)


def qscale(scalar: int, value: Pair, prime: int) -> Pair:
    return (scalar * value[0] % prime, scalar * value[1] % prime)


def qmul(left: Pair, right: Pair, nonsquare: int, prime: int) -> Pair:
    return (
        (left[0] * right[0] + nonsquare * left[1] * right[1]) % prime,
        (left[0] * right[1] + left[1] * right[0]) % prime,
    )


def qinv(value: Pair, nonsquare: int, prime: int) -> Pair:
    norm = (value[0] * value[0] - nonsquare * value[1] * value[1]) % prime
    assert norm != 0
    inverse_norm = pow(norm, -1, prime)
    return (value[0] * inverse_norm % prime, -value[1] * inverse_norm % prime)


def qpow(value: Pair, exponent: int, nonsquare: int, prime: int) -> Pair:
    assert exponent >= 0
    answer = (1, 0)
    base = value
    while exponent:
        if exponent & 1:
            answer = qmul(answer, base, nonsquare, prime)
        base = qmul(base, base, nonsquare, prime)
        exponent >>= 1
    return answer


def qevaluate(
    coefficients: list[int], argument: Pair, nonsquare: int, prime: int
) -> Pair:
    value = (0, 0)
    for coefficient in reversed(coefficients):
        value = qadd(
            qmul(value, argument, nonsquare, prime),
            (coefficient, 0),
            prime,
        )
    return value


def normalized_franel_square(
    value: int,
    fibres: list[list[int]],
    franel: list[int],
    apery: list[int],
    prime: int,
) -> tuple[int, bool]:
    """Recover A_p(value) from Franel-square data on the quadratic cover.

    The Boolean records whether the naive H_p(x)^2 expression fails before
    division by (1+x)^(p-1).  Such failures can only occur on nonsplit fibres.
    """
    direct = evaluate(apery, value, prime)
    if fibres[value]:
        parameter = fibres[value][0]
        h_value = evaluate(franel, parameter, prime)
        assert pow(1 + parameter, prime - 1, prime) == 1
        recovered = h_value * h_value % prime
        assert recovered == direct
        return recovered, False

    discriminant = branch_polynomial(value, prime)
    assert legendre(discriminant, prime) == -1
    inverse_sixteen = pow(16, -1, prime)
    # x=(1-t+sqrt(q(t)))/16 in F_p[sqrt(q(t))].
    parameter = (
        (1 - value) * inverse_sixteen % prime,
        inverse_sixteen,
    )
    parameter_square = qmul(parameter, parameter, discriminant, prime)
    cover_equation = qadd(
        qadd(
            qscale(8, parameter_square, prime),
            qscale(value - 1, parameter, prime),
            prime,
        ),
        (value, 0),
        prime,
    )
    assert cover_equation == (0, 0)

    h_value = qevaluate(franel, parameter, discriminant, prime)
    naive_square = qmul(h_value, h_value, discriminant, prime)
    one_plus_parameter = qadd((1, 0), parameter, prime)
    denominator = qpow(
        one_plus_parameter, prime - 1, discriminant, prime
    )
    recovered_pair = qmul(
        naive_square, qinv(denominator, discriminant, prime), discriminant, prime
    )
    assert recovered_pair[1] == 0
    assert recovered_pair[0] == direct
    return recovered_pair[0], naive_square != (direct, 0)


def prime_factors(value: int) -> list[int]:
    factors = []
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
    raise AssertionError(f"no primitive root found modulo {prime}")


def character_order(index: int, group_order: int) -> int:
    return group_order // gcd(index, group_order)


def character_label(index: int, group_order: int) -> str:
    order = character_order(index, group_order)
    if order == 1:
        return "trivial"
    step = group_order // order
    assert index % step == 0
    phase = index // step
    assert gcd(phase, order) == 1
    return f"ord{order}:phase{phase}"


def tested_characters(group_order: int) -> list[int]:
    return [
        index
        for index in range(group_order)
        if character_order(index, group_order) <= 30
        or 24 % character_order(index, group_order) == 0
    ]


def passing_self_twists(
    values: list[int], prime: int, logarithm: dict[int, int]
) -> list[int]:
    group_order = prime - 1
    support = [value for value in range(1, prime) if values[value] != 0]
    passing = []
    for character in tested_characters(group_order):
        if all(character * logarithm[value] % group_order == 0 for value in support):
            passing.append(character)
    return passing


def mellin_mod(values: list[int], prime: int) -> list[int]:
    group_order = prime - 1
    return [
        sum(
            (values[value] % prime)
            * pow(value, (-character) % group_order, prime)
            for value in range(1, prime)
        )
        % prime
        for character in range(group_order)
    ]


def inversion_shifts(transform: list[int]) -> list[int]:
    order = len(transform)
    return [
        shift
        for shift in range(order)
        if all(
            transform[index] == transform[(-index - shift) % order]
            for index in range(order)
        )
    ]


def digest(values: Iterable[int]) -> str:
    payload = ",".join(str(value) for value in values).encode()
    return hashlib.sha256(payload).hexdigest()[:16]


@dataclass
class PrimeResult:
    prime: int
    primitive_root: int
    fibre_histogram: tuple[int, int, int]
    smooth_trace_checks: int
    nonsplit_normalizations: int
    naive_normalization_failures: int
    pushforward_support: int
    pushforward_sum: int
    pushforward_energy: int
    pushforward_digest: str
    mellin_digest: str
    raw_mellin_failures: tuple[int, ...]
    direct_q6375_matches: int
    tested_character_count: int
    pushforward_self_twists: tuple[str, ...]
    companion_self_twists: tuple[str, ...]
    q_twist_is_multiplicative: bool
    q_partner_difference_support: int
    inversion_shifts_g: tuple[int, ...]
    inversion_shifts_companion: tuple[int, ...]
    inversion_shifts_virtual: tuple[int, ...]
    fibres: list[list[int]]
    pushforward: list[int]
    q_companion: list[int]
    normalized_apery: list[int]
    virtual_mod: list[int]


def check_tasks_one_and_two(prime: int) -> tuple[dict[str, object], list[int]]:
    assert prime not in (2, 3)
    franel = franel_coefficients(prime)
    apery = apery_coefficients(prime)
    pole = prime - 1
    singular = {0, pole, pow(8, -1, prime)}

    traces = [0] * prime
    trace_checks = 0
    for parameter in range(prime):
        if parameter == pole:
            continue
        trace = elliptic_trace(parameter, prime)
        h_value = evaluate(franel, parameter, prime)
        assert trace % prime == h_value
        assert centered_lift(h_value, prime) == trace
        if parameter not in singular:
            assert trace * trace < 4 * prime
            trace_checks += 1
        traces[parameter] = trace

    fibres = [[] for _ in range(prime)]
    for parameter in range(prime):
        if parameter != pole:
            fibres[cover_map(parameter, prime)].append(parameter)

    histogram = [0, 0, 0]
    for value in range(prime):
        actual = len(fibres[value])
        predicted = 1 + legendre(branch_polynomial(value, prime), prime)
        assert actual == predicted
        histogram[actual] += 1
    assert sum(len(fibre) for fibre in fibres) == prime - 1

    pushforward = [
        sum(traces[parameter] ** 2 for parameter in fibre) for fibre in fibres
    ]
    normalized_apery = [0] * prime
    naive_failures = 0
    nonsplit = 0
    for value in range(prime):
        recovered, naive_failed = normalized_franel_square(
            value, fibres, franel, apery, prime
        )
        normalized_apery[value] = recovered
        if not fibres[value]:
            nonsplit += 1
            naive_failures += int(naive_failed)
        assert recovered == evaluate(apery, value, prime)
        predicted_pushforward = (
            (1 + legendre(branch_polynomial(value, prime), prime)) * recovered
        ) % prime
        assert pushforward[value] % prime == predicted_pushforward

    virtual_mod = [
        (
            -pushforward[value]
            + legendre(branch_polynomial(value, prime), prime)
            * normalized_apery[value]
        )
        % prime
        for value in range(prime)
    ]
    assert all(
        virtual_mod[value] == -evaluate(apery, value, prime) % prime
        for value in range(prime)
    )

    raw_failures = []
    for exponent in range(prime):
        raw_mellin = sum(
            virtual_mod[value]
            * pow(value, (-exponent) % (prime - 1), prime)
            for value in range(1, prime)
        ) % prime
        if raw_mellin != apery[exponent]:
            raw_failures.append(exponent)
        endpoint_alias = 0
        if exponent == 0:
            endpoint_alias = apery[prime - 1]
        elif exponent == prime - 1:
            endpoint_alias = apery[0]
        corrected_mellin = (raw_mellin - endpoint_alias) % prime
        assert corrected_mellin == apery[exponent]
    assert raw_failures == [0, prime - 1]
    assert all(
        apery[index] == apery[prime - 1 - index]
        for index in range(prime)
    )

    q_companion = [
        legendre(branch_polynomial(value, prime), prime) * pushforward[value]
        for value in range(prime)
    ]
    direct_q6375_matches = 0
    for exponent in range(prime):
        direct_difference = sum(
            (pushforward[value] - q_companion[value])
            * pow(value, (-exponent) % (prime - 1), prime)
            for value in range(1, prime)
        ) % prime
        direct_q6375_matches += direct_difference == apery[exponent]

    data: dict[str, object] = {
        "franel": franel,
        "apery": apery,
        "traces": traces,
        "fibres": fibres,
        "histogram": tuple(histogram),
        "trace_checks": trace_checks,
        "nonsplit": nonsplit,
        "naive_failures": naive_failures,
        "pushforward": pushforward,
        "normalized_apery": normalized_apery,
        "virtual_mod": virtual_mod,
        "q_companion": q_companion,
        "raw_failures": tuple(raw_failures),
        "direct_q6375_matches": direct_q6375_matches,
    }
    return data, apery


def finish_tasks_three_and_four(prime: int, data: dict[str, object]) -> PrimeResult:
    group_order = prime - 1
    generator = primitive_root(prime)
    logarithm = {
        pow(generator, exponent, prime): exponent for exponent in range(group_order)
    }
    assert len(logarithm) == group_order

    fibres = data["fibres"]
    pushforward = data["pushforward"]
    normalized_apery = data["normalized_apery"]
    virtual_mod = data["virtual_mod"]
    q_companion = data["q_companion"]
    assert isinstance(fibres, list)
    assert isinstance(pushforward, list)
    assert isinstance(normalized_apery, list)
    assert isinstance(virtual_mod, list)
    assert isinstance(q_companion, list)

    passing_g = passing_self_twists(pushforward, prime, logarithm)
    passing_companion = passing_self_twists(q_companion, prime, logarithm)
    labels_g = tuple(character_label(index, group_order) for index in passing_g)
    labels_companion = tuple(
        character_label(index, group_order) for index in passing_companion
    )

    q_factors = [
        legendre(branch_polynomial(value, prime), prime)
        for value in range(prime)
    ]
    q_twist_is_multiplicative = all(
        q_factors[left * right % prime] == q_factors[left] * q_factors[right]
        for left in range(1, prime)
        for right in range(1, prime)
    )
    for value in range(1, prime):
        inverse = pow(value, -1, prime)
        assert (
            branch_polynomial(inverse, prime) * value * value
            - branch_polynomial(value, prime)
        ) % prime == 0
        assert q_factors[inverse] == q_factors[value]
        assert pushforward[inverse] == pushforward[value]
        assert q_companion[inverse] == q_companion[value]
        assert normalized_apery[inverse] == normalized_apery[value]
        assert virtual_mod[inverse] == virtual_mod[value]

    transform_g = mellin_mod(pushforward, prime)
    transform_companion = mellin_mod(q_companion, prime)
    transform_virtual = mellin_mod(virtual_mod, prime)
    shifts_g = tuple(inversion_shifts(transform_g))
    shifts_companion = tuple(inversion_shifts(transform_companion))
    shifts_virtual = tuple(inversion_shifts(transform_virtual))
    assert shifts_g == (0,)
    assert shifts_companion == (0,)
    assert shifts_virtual == (0,)

    histogram = data["histogram"]
    raw_failures = data["raw_failures"]
    assert isinstance(histogram, tuple)
    assert isinstance(raw_failures, tuple)
    return PrimeResult(
        prime=prime,
        primitive_root=generator,
        fibre_histogram=histogram,
        smooth_trace_checks=int(data["trace_checks"]),
        nonsplit_normalizations=int(data["nonsplit"]),
        naive_normalization_failures=int(data["naive_failures"]),
        pushforward_support=sum(pushforward[value] != 0 for value in range(1, prime)),
        pushforward_sum=sum(pushforward[1:]),
        pushforward_energy=sum(value * value for value in pushforward[1:]),
        pushforward_digest=digest(pushforward),
        mellin_digest=digest(transform_g),
        raw_mellin_failures=raw_failures,
        direct_q6375_matches=int(data["direct_q6375_matches"]),
        tested_character_count=len(tested_characters(group_order)),
        pushforward_self_twists=labels_g,
        companion_self_twists=labels_companion,
        q_twist_is_multiplicative=q_twist_is_multiplicative,
        q_partner_difference_support=sum(
            pushforward[value] != q_companion[value] for value in range(1, prime)
        ),
        inversion_shifts_g=shifts_g,
        inversion_shifts_companion=shifts_companion,
        inversion_shifts_virtual=shifts_virtual,
        fibres=fibres,
        pushforward=pushforward,
        q_companion=q_companion,
        normalized_apery=normalized_apery,
        virtual_mod=virtual_mod,
    )


def run_all() -> list[PrimeResult]:
    # The Task-2 sanity gate is completed for every prime before any
    # self-twist or duality test is allowed to run.
    partial: list[tuple[int, dict[str, object]]] = []
    for prime in PRIMES:
        data, _ = check_tasks_one_and_two(prime)
        partial.append((prime, data))
    assert len(partial) >= 3
    return [finish_tasks_three_and_four(prime, data) for prime, data in partial]


def print_summary(results: list[PrimeResult]) -> None:
    print(
        "p  fibres(N=0/1/2)  smooth  Fp2 naive-fail  supp(G)  "
        "chars  self(G)  self(Gq)  c"
    )
    for result in results:
        fibre_text = "/".join(str(value) for value in result.fibre_histogram)
        normalization = (
            f"{result.nonsplit_normalizations}/"
            f"{result.naive_normalization_failures}"
        )
        print(
            f"{result.prime:<3}{fibre_text:>17}"
            f"{result.smooth_trace_checks:>8}{normalization:>17}"
            f"{result.pushforward_support:>9}{result.tested_character_count:>7}  "
            f"{','.join(result.pushforward_self_twists):<9}"
            f"{','.join(result.companion_self_twists):<10}0"
        )
        print(
            f"    TG sha={result.pushforward_digest}, Mellin sha={result.mellin_digest}, "
            f"sum={result.pushforward_sum}, energy={result.pushforward_energy}, "
            f"raw endpoint failures={list(result.raw_mellin_failures)}, "
            f"Q6375-direct matches={result.direct_q6375_matches}/{result.prime}"
        )
    assert all(
        result.pushforward_self_twists == ("trivial",)
        and result.companion_self_twists == ("trivial",)
        for result in results
    )
    assert all(not result.q_twist_is_multiplicative for result in results)
    print("SANITY TASK 1: PASS for every t and every requested prime")
    print(
        "SANITY TASK 2: PASS for every 0<=r<p after the two explicit endpoint "
        "alias corrections; raw formula passes every 1<=r<=p-2"
    )
    print("TASK 3: only the trivial tested t-Kummer self-twist survives")
    print("TASK 4: unique inversion shift c=0 for G, Gq, and the virtual trace")
    print(
        "TASK 5: m=1 Mellin vectors computed; extension-field traces and "
        "L-function degree probe skipped"
    )
    print(
        "RANK: Q6375 rank 18 is refuted: rank(F)=2 gives rank(phi_*(F tensor F))=8; "
        "the Sym^2 reduction has rank 6"
    )


def dump_traces(results: list[PrimeResult]) -> None:
    print("p,t,N,T_G,T_G_chi_q,A_from_Franel,virtual_mod")
    for result in results:
        for value in range(result.prime):
            print(
                result.prime,
                value,
                len(result.fibres[value]),
                result.pushforward[value],
                result.q_companion[value],
                result.normalized_apery[value],
                result.virtual_mod[value],
                sep=",",
            )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--dump-traces",
        action="store_true",
        help="print the full reproducible trace substrate as CSV",
    )
    arguments = parser.parse_args()
    results = run_all()
    if arguments.dump_traces:
        dump_traces(results)
    else:
        print_summary(results)


if __name__ == "__main__":
    main()
