#!/usr/bin/env python3
"""Temporary spectral scan for an ascended 4F3(1) source of P2.5."""

import itertools
import math
import sys

import numpy as np
import sympy as sp

from ramanujantools.cmf import pFq


def canonical_vectors(count, max_length):
    values = range(-max_length, max_length + 1)
    for vector in itertools.combinations_with_replacement(values, count):
        length = sum(abs(value) for value in vector)
        if 0 < length <= max_length:
            yield vector


def main():
    max_length = int(sys.argv[1])
    axes = sp.symbols("x0 x1 x2 x3 y0 y1 y2")
    cmf = pFq(4, 3, 1)
    forms = {
        (index, sign): sp.lambdify(axes, cmf.M(axis, sign > 0), modules="numpy")
        for index, axis in enumerate(axes)
        for sign in (1, -1)
    }
    # Ascend 3F2(1/2,1/2,1/2;1,3/2;1); the last pair is stationary.
    base = np.array([0.5, 0.5, 0.5, 2.0, 1.0, 1.5, 2.0])

    def step_matrix(trajectory, n=100000.0):
        position = base + n * np.array(trajectory, dtype=float)
        result = np.eye(3)
        with np.errstate(all="ignore"):
            for index in reversed(range(7)):
                sign = 1 if trajectory[index] >= 0 else -1
                for _ in range(abs(trajectory[index])):
                    current = np.asarray(forms[index, sign](*position), dtype=float)
                    result = result @ current
                    position[index] += sign
        return result

    x_vectors = list(canonical_vectors(4, max_length))
    y_vectors = list(canonical_vectors(3, max_length)) + [(0, 0, 0)]
    hits = []
    tested = 0
    for xs in x_vectors:
        x_length = sum(abs(value) for value in xs)
        for ys in y_vectors:
            if x_length + sum(abs(value) for value in ys) > max_length:
                continue
            trajectory = xs + ys
            tested += 1
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
                score = abs(math.log(abs(invariant1 / 1225))) + abs(
                    math.log(abs(invariant2 / 42875))
                )
                if score < 0.02:
                    hits.append((score, trajectory, invariant1, invariant2, eigenvalues))
            except (TypeError, ValueError, np.linalg.LinAlgError):
                pass
    print("tested", tested)
    for hit in sorted(hits, key=lambda item: item[0])[:200]:
        print(hit)


if __name__ == "__main__":
    main()
