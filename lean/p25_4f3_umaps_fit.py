#!/usr/bin/env python3
"""Numerically recover the free initial row in the Catalan 4F3 gauge.

The public limit rows S satisfy S U(0) = A_challenge.  This fixes the last
two rows of U(0) and leaves its first row free.  We propagate the four affine
basis matrices by

    C_4F3(n) U(n+1) = U(n) T_25(n)

and choose the free row so that all projective entries of U(n) are rational
functions of n.  This is the rank-three analogue of UMAPS' empirical step.
"""

from __future__ import annotations

import argparse
from fractions import Fraction

import numpy as np
from scipy.optimize import least_squares
import sympy as sp

import p25_4f3_interpolate as data


def affine_samples(count: int) -> np.ndarray:
    constant = sp.Matrix([[0, 0, 0], *data.lower])
    bases = [constant]
    for column in range(3):
        value = sp.zeros(3)
        value[0, column] = 1
        bases.append(value)

    samples = []
    for index in range(count):
        coefficients = np.empty((9, 4), dtype=float)
        for basis_index, value in enumerate(bases):
            coefficients[:, basis_index] = [float(entry) for entry in value]
        # A common scale at each n is projectively irrelevant and prevents
        # overflow.  Keep the four affine directions on the same scale.
        coefficients /= np.max(np.abs(coefficients))
        samples.append(coefficients)
        if index + 1 < count:
            left_inverse = data.source(index).inv()
            right = data.challenge(index)
            bases = [left_inverse * value * right for value in bases]
    return np.asarray(samples)


def rational_holdout(values: np.ndarray, degree: int, holdout: int) -> np.ndarray:
    train = 2 * degree + 1
    x = np.linspace(-1.0, 1.0, train + holdout)
    matrix = np.empty((train, 2 * degree + 1))
    # p_0+...+p_d*x^d = y*(1+q_1*x+...+q_d*x^d)
    for row in range(train):
        powers = x[row] ** np.arange(degree + 1)
        matrix[row, : degree + 1] = powers
        matrix[row, degree + 1 :] = -values[row] * powers[1:]
    try:
        coefficients = np.linalg.solve(matrix, values[:train])
    except np.linalg.LinAlgError:
        return np.full(holdout, 1e6)
    answer = []
    for row in range(train, train + holdout):
        powers = x[row] ** np.arange(degree + 1)
        numerator = coefficients[: degree + 1] @ powers
        denominator = 1 + coefficients[degree + 1 :] @ powers[1:]
        prediction = numerator / denominator
        answer.append((prediction - values[row]) / (1 + abs(values[row])))
    return np.asarray(answer)


def fit(samples: np.ndarray, degree: int, holdout: int, seed: int,
        starts_count: int, reference: int):
    rng = np.random.default_rng(seed)
    affine = np.array([1.0, 0.0, 0.0, 0.0])

    def residual(parameters: np.ndarray) -> np.ndarray:
        affine[1:] = parameters
        entries = samples @ affine
        # A fixed reference keeps the residual smooth for least squares.
        candidates = []
        for current_reference in ([reference] if reference >= 0 else range(9)):
            denominator = entries[:, current_reference]
            if np.min(np.abs(denominator)) < 1e-12 * np.max(np.abs(denominator)):
                continue
            ratios = entries / denominator[:, None]
            pieces = []
            for entry in range(9):
                if entry != current_reference:
                    pieces.append(rational_holdout(ratios[:, entry], degree, holdout))
            candidates.append(np.concatenate(pieces))
        if not candidates:
            return np.full(8 * holdout, 1e6)
        return min(candidates, key=lambda value: np.linalg.norm(value))

    starts = [np.zeros(3), np.ones(3)]
    starts += [rng.uniform(-80, 80, 3) for _ in range(starts_count)]
    best = None
    for start in starts:
        result = least_squares(
            residual,
            start,
            max_nfev=800,
            xtol=1e-14,
            ftol=1e-14,
            gtol=1e-14,
        )
        score = np.linalg.norm(result.fun)
        if best is None or score < best[0]:
            best = score, result.x
            print("degree", degree, "best", score, result.x, flush=True)
    return best


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-degree", type=int, default=7)
    parser.add_argument("--min-degree", type=int, default=1)
    parser.add_argument("--holdout", type=int, default=5)
    parser.add_argument("--seed", type=int, default=2505)
    parser.add_argument("--starts", type=int, default=6)
    parser.add_argument("--reference", type=int, default=3)
    args = parser.parse_args()
    count = 2 * args.max_degree + 1 + args.holdout
    samples = affine_samples(count)
    for degree in range(args.min_degree, args.max_degree + 1):
        best = fit(samples[: 2 * degree + 1 + args.holdout], degree,
                   args.holdout, args.seed + degree, args.starts,
                   args.reference)
        print("RESULT", degree, best[0], best[1],
              [Fraction(float(value)).limit_denominator(1000000)
               for value in best[1]], flush=True)


if __name__ == "__main__":
    main()
