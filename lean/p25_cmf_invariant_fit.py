#!/usr/bin/env python3
"""Temporary projective spectral fit for a P2.5 CMF trajectory."""

import itertools
import sys

import numpy as np
from scipy.optimize import least_squares

import p25_meijer_direct_fit as direct


def invariants(matrix):
    e1 = np.trace(matrix)
    e2 = (e1*e1-np.trace(matrix@matrix))/2
    e3 = np.linalg.det(matrix)
    return np.array([e1*e2/e3, e1**3/e3])


def fit(name, forms, axes, directions):
    rng = np.random.default_rng(2525)
    targets = {
        "forward": [invariants(direct.target(n)) for n in range(5)],
        "inverse": [invariants(direct.inverse3(direct.target(n))) for n in range(5)],
    }
    for direction in directions:
        for orientation, wanted in targets.items():
            def residual(start):
                try:
                    values = []
                    for n in range(5):
                        current = invariants(
                            direct.trajectory_value(forms, start, direction, n))
                        values.extend(current/wanted[n]-1)
                    values = np.asarray(values)
                    if not np.isfinite(values).all():
                        raise ValueError
                    return values
                except (ValueError, ZeroDivisionError, np.linalg.LinAlgError):
                    return np.full(10, 1e6)

            seeds = [np.zeros(len(axes)), np.full(len(axes), 0.5)]
            seeds += [rng.integers(-8, 9, len(axes))/2 for _ in range(8)]
            best = None
            for seed in seeds:
                result = least_squares(
                    residual, seed, max_nfev=350,
                    xtol=1e-13, ftol=1e-13, gtol=1e-13,
                )
                score = np.linalg.norm(result.fun)
                if best is None or score < best[0]:
                    best = score, result.x
                if score < 1e-9:
                    break
            print(name, direction, orientation, best, flush=True)


def main():
    family = sys.argv[1]
    if family == "derived":
        axes, forms = direct.derived_pfq_elementary(-1)
        directions = [(-2, 2, 0, 0, 0), (0, 0, 0, -2, 2)]
    elif family == "pfq":
        axes, forms = direct.pfq_elementary(-1)
        directions = [(-2, 2, 0, 0, 0), (0, 0, 0, -2, 2)]
    elif family == "meijer":
        axes, forms = direct.meijer_elementary(-1)
        directions = [(-2, 0, 2, 0, 0, 0), (0, 0, 0, -2, 0, 2)]
    else:
        raise SystemExit(f"unknown family {family}")
    fit(family, forms, axes, directions)


if __name__ == "__main__":
    main()
