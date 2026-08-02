#!/usr/bin/env python3
"""Search the Catalan 4F3 CMF for the Problem 2.5 projective spectrum.

The public CMF_SCANNER experiment fixes the Catalan start point.  A trajectory
is identified here only through the eigenvalue ratios of its matrix at
infinity, which are invariant under scalar normalization and change of basis.
"""

from __future__ import annotations

import argparse
import itertools
import math
import random

import numpy as np
import sympy as sp

from ramanujantools.cmf import pFq


TARGET = np.array([1.0, 17.0 - 12.0 * math.sqrt(2.0),
                   (17.0 - 12.0 * math.sqrt(2.0)) ** 2])


def primitive(values: tuple[int, ...]) -> bool:
    g = 0
    for value in values:
        g = math.gcd(g, abs(value))
    return g == 1


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bound", type=int, default=3)
    parser.add_argument("--random", type=int, default=0,
                        help="sample this many unsorted trajectories")
    parser.add_argument("--seed", type=int, default=25042026)
    parser.add_argument("--n", type=float, default=1.0e6)
    parser.add_argument("--keep", type=int, default=30)
    args = parser.parse_args()

    cmf = pFq(4, 3, 1)
    axes = sorted(cmf.axes(), key=str)
    symbols = sp.symbols("x:4") + sp.symbols("y:3")
    axis_functions = {}
    for axis in axes:
        for sign in (False, True):
            matrix = sp.Matrix(cmf.M(axis, sign).tolist())
            axis_functions[(axis, sign)] = sp.lambdify(symbols, matrix, "numpy")

    n = args.n
    balance = np.fromfunction(lambda i, j: n ** (i - j), (3, 3), dtype=float)

    def score(values: tuple[int, ...]):
        coords = tuple(n * value for value in values)
        product = np.eye(3)
        try:
            for axis, amount in zip(axes, values):
                if amount == 0:
                    continue
                step = np.asarray(
                    axis_functions[(axis, amount > 0)](*coords), dtype=float
                ) * balance
                if not np.all(np.isfinite(step)):
                    return None
                for _ in range(abs(amount)):
                    product = product @ step
                    norm = np.max(np.abs(product))
                    if not math.isfinite(norm) or norm == 0:
                        return None
                    product /= norm
            eig = np.linalg.eigvals(product)
        except (FloatingPointError, ZeroDivisionError, ValueError):
            return None
        if np.max(np.abs(eig.imag)) > 1.0e-5 * max(1.0, np.max(np.abs(eig.real))):
            return None
        ratios = np.sort(np.abs(eig.real))[::-1]
        if ratios[0] == 0 or ratios[-1] == 0:
            return None
        ratios /= ratios[0]
        distance = float(np.max(np.abs(np.log(ratios / TARGET))))
        return distance, ratios, eig

    if args.random:
        rng = random.Random(args.seed)

        def candidates():
            for _ in range(args.random):
                values = tuple(rng.randint(-args.bound, args.bound) for _ in axes)
                if values != (0,) * len(axes):
                    yield values
    else:
        domain = range(-args.bound, args.bound + 1)

        def candidates():
            # Numerator and denominator parameters are separately symmetric.
            for xs in itertools.combinations_with_replacement(domain, 4):
                for ys in itertools.combinations_with_replacement(domain, 3):
                    yield xs + ys

    best: list[tuple[float, tuple[int, ...], np.ndarray, np.ndarray]] = []
    checked = 0
    for values in candidates():
        if not primitive(values):
            continue
        # Generic rays are exactly those used by CMF_SCANNER's sampler.
        if set(values[:4]) & set(values[4:]):
            continue
        if sum(values[:4]) == sum(values[4:]):
            continue
        result = score(values)
        if result is None:
            continue
        checked += 1
        distance, ratios, eig = result
        item = (distance, values, ratios, eig)
        if len(best) < args.keep or distance < best[-1][0]:
            best.append(item)
            best.sort(key=lambda item: item[0])
            del best[args.keep:]

    print("checked", checked)
    print("target", TARGET)
    for distance, values, ratios, eig in best:
        print(f"{distance:.12g}", values, ratios, eig)


if __name__ == "__main__":
    main()
