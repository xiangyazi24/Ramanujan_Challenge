#!/usr/bin/env python3
"""Scan short trajectories in the derived 3F2 CMF for the P2.5 spectrum."""

import itertools
import math
import sys

import numpy as np
import sympy as sp

from ramanujantools.cmf import known_cmfs


def compositions(total, parts, prefix=()):
    if parts == 1:
        yield prefix + (total,)
        return
    for value in range(total + 1):
        yield from compositions(total - value, parts - 1, prefix + (value,))


def main():
    max_length = int(sys.argv[1])
    cmf = known_cmfs.hypergeometric_derived_3F2().subs({sp.Symbol("z"): -1})
    axes = sp.symbols("x0 x1 x2 y0 y1")
    forms = {
        (index, sign): sp.lambdify(axes, cmf.M(axis, sign > 0), modules="numpy")
        for index, axis in enumerate(axes)
        for sign in (1, -1)
    }
    # G = 3F2(1/2,1/2,1;3/2,3/2;-1), with y_i = b_i - 1.
    base = np.array([0.5, 0.5, 1.0, 0.5, 0.5])

    def step_matrix(trajectory, depth=100000.0):
        position = base + depth * np.array(trajectory, dtype=float)
        result = np.eye(3)
        with np.errstate(all="ignore"):
            for index in reversed(range(5)):
                sign = 1 if trajectory[index] >= 0 else -1
                for _ in range(abs(trajectory[index])):
                    result = result @ np.asarray(
                        forms[index, sign](*position), dtype=float
                    )
                    position[index] += sign
        return result

    target1 = 1225.0
    target2 = 42875.0
    hits = []
    best = []
    for length in range(1, max_length + 1):
        for magnitudes in compositions(length, 5):
            nonzero = [index for index, value in enumerate(magnitudes) if value]
            for signs in itertools.product((-1, 1), repeat=len(nonzero)):
                trajectory = list(magnitudes)
                for index, sign in zip(nonzero, signs):
                    trajectory[index] *= sign
                try:
                    eigenvalues = np.linalg.eigvals(step_matrix(trajectory))
                    e1 = sum(eigenvalues)
                    e2 = (
                        eigenvalues[0] * eigenvalues[1]
                        + eigenvalues[0] * eigenvalues[2]
                        + eigenvalues[1] * eigenvalues[2]
                    )
                    e3 = np.prod(eigenvalues)
                    if not np.isfinite(eigenvalues).all() or abs(e3) < 1e-100:
                        continue
                    invariant1 = float(np.real_if_close(e1 * e2 / e3))
                    invariant2 = float(np.real_if_close(e1**3 / e3))
                    if invariant1 == 0 or invariant2 == 0:
                        continue
                    score = abs(math.log(abs(invariant1 / target1))) + abs(
                        math.log(abs(invariant2 / target2))
                    )
                    item = (score, tuple(trajectory), invariant1, invariant2, eigenvalues)
                    best.append(item)
                    if score < 1e-4:
                        hits.append(item)
                except (TypeError, ValueError, np.linalg.LinAlgError):
                    pass
        best.sort(key=lambda item: item[0])
        best = best[:20]
        print("length", length, "hits", hits[-20:], "best", best[:5], flush=True)


if __name__ == "__main__":
    main()
