#!/usr/bin/env sage
"""Modular guessing for the complete 14-ray first-cell correction.

This is an exploratory companion to ``q32_first_cell_ray_telescopers.sage``.
It does not turn a modular guess into a theorem.  Its purpose is to locate
the generic minimal order/degree before attempting exact reconstruction.

For fixed ``M`` it computes, modulo a large prime,

    F_M(r) = sum_{kappa in P(Z), kappa != 0}
             c_M((M-r) kappa)

on the full range ``0 <= r <= 2*M`` and searches for

    sum_{j=0}^R a_j(r) F_M(r+j) = 0,
    deg_r(a_j) <= D.

The full range is legitimate for the ray sum itself; it coincides with
``C_M(M-r)-b_M`` on the first cell ``0 <= r < M/2``.
"""

import argparse
from math import comb


RAY_CLASSES = (
    ((-1, -1, -1), 1),
    ((-1, -1, 0), 2),
    ((-1, -1, 1), 2),
    ((-1, 0, 0), 1),
    ((-1, 0, 1), 2),
    ((-1, 1, 1), 1),
    ((0, -1, -1), 1),
    ((0, -1, 0), 2),
    ((0, -1, 1), 2),
    ((0, 0, 1), 2),
    ((0, 1, 1), 1),
    ((1, 0, 0), 1),
    ((1, 0, 1), 2),
    ((1, 1, 1), 1),
)


def binomial_rows(moment, prime):
    rows = {}
    for upper in range(moment, 2 * moment + 1):
        rows[upper] = [
            comb(upper, lower) % prime
            for lower in range(upper + 1)
        ]
    return rows


def lookup(rows, upper, lower):
    return rows[upper][lower] if 0 <= lower <= upper else 0


def ray_sum_values(moment, prime):
    rows = binomial_rows(moment, prime)
    moment_row = rows[moment]
    values = []
    for residue in range(2 * moment + 1):
        node = moment - residue
        total = 0
        for (u, v, w), multiplicity in RAY_CLASSES:
            ray = 0
            for index in range(moment + 1):
                upper = 2 * moment - index
                ray += (
                    moment_row[index]
                    * lookup(rows, moment, index - node * u)
                    * lookup(rows, upper, moment - node * v)
                    * lookup(rows, upper, moment - node * w)
                )
            total += multiplicity * ray
        values.append(total % prime)
    return values


def guess(values, prime, order, degree, holdout):
    field = GF(prime)
    columns = (order + 1) * (degree + 1)
    equation_count = len(values) - order
    training_count = equation_count - holdout
    if training_count <= columns:
        return None

    rows = []
    for residue in range(training_count):
        powers = [1]
        for _ in range(degree):
            powers.append(powers[-1] * residue % prime)
        rows.append(
            [
                values[residue + shift] * powers[power] % prime
                for shift in range(order + 1)
                for power in range(degree + 1)
            ]
        )
    kernel = Matrix(field, rows).right_kernel()
    if kernel.dimension() == 0:
        return None

    for vector in kernel.basis():
        valid = True
        for residue in range(training_count, equation_count):
            total = 0
            for shift in range(order + 1):
                polynomial = 0
                for power in reversed(range(degree + 1)):
                    polynomial = (
                        polynomial * residue
                        + ZZ(vector[shift * (degree + 1) + power])
                    ) % prime
                total += polynomial * values[residue + shift]
            if total % prime:
                valid = False
                break
        if valid:
            return vector
    return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--moment", type=int, default=199)
    parser.add_argument("--prime", type=int, default=1000000007)
    parser.add_argument("--max-order", type=int, default=38)
    parser.add_argument("--max-degree", type=int, default=40)
    parser.add_argument("--holdout", type=int, default=24)
    args = parser.parse_args()

    values = ray_sum_values(args.moment, args.prime)
    print("VALUES", args.moment, len(values), args.prime)
    for order in range(1, args.max_order + 1):
        for degree in range(args.max_degree + 1):
            vector = guess(
                values,
                args.prime,
                order,
                degree,
                args.holdout,
            )
            if vector is not None:
                print(
                    "GUESS",
                    "order",
                    order,
                    "degree",
                    degree,
                    "kernel_vector",
                    list(vector),
                )
                return
        print("NO_GUESS_THROUGH", order, args.max_degree)
    print("NO_GUESS")


if __name__ == "__main__":
    main()
