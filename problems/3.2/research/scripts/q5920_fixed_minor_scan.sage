#!/usr/bin/env sage
"""Exact modular scan of det(R_{-1,m}, R_{1,m}) beyond n=500.

This uses exactly the operator and shell definitions at commit 55a92ba.
For each central core it computes all shell data modulo Q_I^2 once, then
extracts the p^2 jet separately for every core prime p.  No factorization or
floating-point arithmetic is used.
"""

from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
load(str(HERE / "q32_doubled_period_gauge_audit.sage"))


def primitive_integral_coefficients(operator):
    coefficients = [K(operator[i]) for i in range(operator.order() + 1)]
    polynomial_denominator = lcm(value.denominator() for value in coefficients)
    polynomials = [R(polynomial_denominator * value) for value in coefficients]
    rational_denominator = lcm(
        coefficient.denominator()
        for polynomial in polynomials
        for coefficient in polynomial.list()
    )
    polynomials = [R(rational_denominator * value) for value in polynomials]
    content = gcd([
        ZZ(coefficient)
        for polynomial in polynomials
        for coefficient in polynomial.list()
    ])
    polynomials = [R(value / content) for value in polynomials]
    assert gcd(polynomials) == 1
    return tuple(polynomials)


common_operator = P.lclm(apery_operator)
assert common_operator.order() == 5
ghost_coefficients = primitive_integral_coefficients(common_operator)
assert len(ghost_coefficients) == 6


def arg_value(name, default):
    prefix = "--" + name + "="
    for argument in sys.argv:
        if argument.startswith(prefix):
            return ZZ(argument[len(prefix):])
    return ZZ(default)


START = arg_value("start", 520)
STOP = arg_value("stop", 1000)
STEP = arg_value("step", 20)


def carrier_mod(values, lower, length, modulus):
    return ZZ(sum(
        (-1)^offset
        * binomial(lower + offset, offset)
        * binomial(lower + length + 1, length - offset)
        * ZZ(values[lower + offset])
        for offset in range(length + 1)
    ) % modulus)


def quotient_coordinate(value_mod_p2, p, cofactor):
    value_mod_p2 = ZZ(value_mod_p2 % (p^2))
    assert value_mod_p2 % p == 0
    return ZZ((value_mod_p2 // p) * inverse_mod(cofactor % p, p) % p)


def scan_index(index):
    D = (13 * index) // 20
    N = index // 5
    moment = index - 1
    maximum_margin = min(
        D - moment // 2,
        moment - D - N + 2,
    )
    assert maximum_margin >= 2
    margins = (maximum_margin - 1, maximum_margin)

    candidates = tuple(prime_range(D + 1, D + N + 1))
    if not candidates:
        return []
    candidate_product = ZZ(prod(candidates))
    modulus = candidate_product^2

    lower = D - maximum_margin + 1
    upper = D + N + maximum_margin - 2
    nodes = tuple(range(lower, upper + 1))
    neighbours = tuple(range(lower - 1, upper + 2))

    Y = shell_batch(moment, nodes, modulus=modulus)
    time_shell_rows = {
        time: shell_batch(time, neighbours, modulus=modulus)
        for time in range(index, index + 7)
    }

    def shell_pair(time, node):
        return ZZ(
            time_shell_rows[time][node - 1]
            + time_shell_rows[time][node + 1]
        ) % modulus

    origin_lift = {
        node: ZZ(sum(
            ZZ(integer_coefficients[shift](n=index))
            * shell_pair(index + shift, node)
            for shift in range(origin_operator.order() + 1)
        ) % modulus)
        for node in nodes
    }
    multiplier_value = ZZ(multiplier(n=index))
    W = {
        node: ZZ((origin_lift[node] - multiplier_value * Y[node]) % modulus)
        for node in nodes
    }
    X1 = {
        node: ZZ(sum(
            ZZ(ghost_coefficients[shift](n=index + 1))
            * shell_pair(index + 1 + shift, node)
            for shift in range(len(ghost_coefficients))
        ) % modulus)
        for node in nodes
    }

    failures = []
    for margin in margins:
        d = D - margin + 1
        length = N + margin - 2
        W_left_all = carrier_mod(W, d, length, modulus)
        W_right_all = carrier_mod(W, D, length, modulus)
        X_left_all = carrier_mod(X1, d, length, modulus)
        X_right_all = carrier_mod(X1, D, length, modulus)

        for p in candidates:
            p2 = p^2
            cofactor = candidate_product // p
            W_left = quotient_coordinate(W_left_all, p, cofactor)
            W_right = quotient_coordinate(W_right_all, p, cofactor)
            X_left = quotient_coordinate(X_left_all, p, cofactor)
            X_right = quotient_coordinate(X_right_all, p, cofactor)
            determinant = ZZ((W_left * X_right - W_right * X_left) % p)
            if determinant == 0:
                failures.append((index, D, N, margin, p, index - p,
                                 W_left, W_right, X_left, X_right))
                print("FIXED_MINOR_FAILURE", failures[-1])
                return failures

    print("FIXED_MINOR_PASS", index, (D, N), margins, len(candidates))
    return failures


all_failures = []
incidences = 0
for index in range(START, STOP + 1, STEP):
    failures = scan_index(ZZ(index))
    all_failures.extend(failures)
    if failures:
        break

print("Q5920_FIXED_MINOR_SCAN_DONE", {
    "range": (START, STOP, STEP),
    "failures": tuple(all_failures),
})
assert not all_failures
