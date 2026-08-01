#!/usr/bin/env python3
"""Temporary spectral scan for the direct Catalan 3F2 at z = -1."""

import itertools
import math
import sys

import numpy as np
import sympy as sp


def compositions(total, parts, prefix=()):
    if parts == 1:
        yield prefix + (total,)
        return
    for value in range(total + 1):
        yield from compositions(total - value, parts - 1, prefix + (value,))


def main():
    max_length = int(sys.argv[1])
    theta = sp.symbols("theta")
    axes = sp.symbols("x0 x1 x2 y0 y1")
    x0, x1, x2, y0, y1 = axes
    # theta(theta+y0-1)(theta+y1-1) + prod(theta+xi), since z=-1.
    polynomial = sp.Poly(
        sp.expand(theta * (theta + y0 - 1) * (theta + y1 - 1)
                  + (theta + x0) * (theta + x1) * (theta + x2)),
        theta,
    ).monic()
    coefficients = polynomial.all_coeffs()
    companion = sp.Matrix([
        [0, 0, -coefficients[3]],
        [1, 0, -coefficients[2]],
        [0, 1, -coefficients[1]],
    ])
    eye = sp.eye(3)
    native = [eye + companion / x0, eye + companion / x1,
              eye + companion / x2, eye + companion / (y0 - 1),
              eye + companion / (y1 - 1)]
    positive = native[:3] + [
        native[3].subs(y0, y0 + 1).inv(),
        native[4].subs(y1, y1 + 1).inv(),
    ]
    negative = [
        native[0].subs(x0, x0 - 1).inv(),
        native[1].subs(x1, x1 - 1).inv(),
        native[2].subs(x2, x2 - 1).inv(),
        native[3], native[4],
    ]
    forms = {
        (index, sign): sp.lambdify(axes,
            positive[index] if sign > 0 else negative[index], modules="numpy")
        for index in range(5) for sign in (-1, 1)
    }
    base = np.array([0.5, 0.5, 1.0, 1.5, 1.5])

    def step_matrix(trajectory, n=100000.0):
        position = base + n * np.array(trajectory, dtype=float)
        result = np.eye(3)
        with np.errstate(all="ignore"):
            for index in reversed(range(5)):
                direction = 1 if trajectory[index] >= 0 else -1
                for _ in range(abs(trajectory[index])):
                    result = result @ np.asarray(
                        forms[index, direction](*position), dtype=float)
                    position[index] += direction
        return result

    hits = []
    tested = 0
    for length in range(1, max_length + 1):
        for magnitudes in compositions(length, 5):
            nonzero = [i for i, value in enumerate(magnitudes) if value]
            for signs in itertools.product((-1, 1), repeat=len(nonzero)):
                trajectory = list(magnitudes)
                for index, sign in zip(nonzero, signs):
                    trajectory[index] *= sign
                tested += 1
                try:
                    eigenvalues = np.linalg.eigvals(step_matrix(trajectory))
                    e1 = sum(eigenvalues)
                    e2 = (eigenvalues[0] * eigenvalues[1]
                          + eigenvalues[0] * eigenvalues[2]
                          + eigenvalues[1] * eigenvalues[2])
                    e3 = np.prod(eigenvalues)
                    if not np.isfinite(eigenvalues).all() or abs(e3) < 1e-100:
                        continue
                    invariant1 = float(np.real_if_close(e1 * e2 / e3))
                    invariant2 = float(np.real_if_close(e1**3 / e3))
                    if invariant1 == 0 or invariant2 == 0:
                        continue
                    score = abs(math.log(abs(invariant1 / 1225))) + abs(
                        math.log(abs(invariant2 / 42875)))
                    if score < 0.01:
                        hits.append((score, tuple(trajectory), invariant1,
                                     invariant2, eigenvalues))
                except (TypeError, ValueError, np.linalg.LinAlgError):
                    pass
        print("length", length, "tested", tested, "hits", len(hits), flush=True)
    for hit in sorted(hits, key=lambda item: item[0])[:200]:
        print(hit, flush=True)


if __name__ == "__main__":
    main()
