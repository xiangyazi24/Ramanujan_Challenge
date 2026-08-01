#!/usr/bin/env python3
"""Temporary spectral scan for a rank-three Meijer-G trajectory behind P2.5."""

import itertools
import math
import sys

import numpy as np
import sympy as sp

from ramanujantools.cmf import MeijerG


def canonical_vectors(count, max_length):
    values = range(-max_length, max_length + 1)
    for vector in itertools.combinations_with_replacement(values, count):
        length = sum(abs(value) for value in vector)
        if length <= max_length:
            yield vector


def scan(p, q, sign_class, max_length):
    # m=0; choose n parity so (-1)^(p-m-n) equals sign_class.
    n_param = next(value for value in range(p + 1) if (-1) ** (p - value) == sign_class)
    cmf = MeijerG(0, n_param, p, q, 1)
    axes = tuple(MeijerG.a_axes(p)) + tuple(MeijerG.b_axes(q))
    forms = {
        (index, sign): sp.lambdify(axes, cmf.M(axis, sign > 0), modules="numpy")
        for index, axis in enumerate(axes)
        for sign in (1, -1)
    }
    base = np.arange(1, len(axes) + 1, dtype=float) / 7.0 + 0.13

    def step_matrix(trajectory, sample_n=100000.0):
        position = base + sample_n * np.array(trajectory, dtype=float)
        result = np.eye(3)
        with np.errstate(all="ignore"):
            for index in reversed(range(len(axes))):
                direction = 1 if trajectory[index] >= 0 else -1
                for _ in range(abs(trajectory[index])):
                    current = np.asarray(forms[index, direction](*position), dtype=float)
                    result = result @ current
                    position[index] += direction
        scale = np.max(np.abs(result))
        return result / scale if scale else result

    avectors = list(canonical_vectors(p, max_length))
    bvectors = list(canonical_vectors(q, max_length))
    hits = []
    tested = 0
    for aa in avectors:
        alength = sum(abs(value) for value in aa)
        for bb in bvectors:
            length = alength + sum(abs(value) for value in bb)
            if length == 0 or length > max_length:
                continue
            trajectory = aa + bb
            tested += 1
            try:
                eigenvalues = np.linalg.eigvals(step_matrix(trajectory))
                e1 = sum(eigenvalues)
                e2 = sum(
                    eigenvalues[i] * eigenvalues[j]
                    for i in range(3)
                    for j in range(i + 1, 3)
                )
                e3 = np.prod(eigenvalues)
                if not np.isfinite(eigenvalues).all() or abs(e3) < 1e-100:
                    continue
                invariant1 = float(np.real_if_close(e1 * e2 / e3))
                invariant2 = float(np.real_if_close(e1**3 / e3))
                score = abs(math.log(abs(invariant1 / 1225))) + abs(
                    math.log(abs(invariant2 / 42875))
                )
                if score < 0.02:
                    hits.append((score, trajectory, invariant1, invariant2, eigenvalues))
            except (TypeError, ValueError, np.linalg.LinAlgError):
                pass
    print("CASE", p, q, sign_class, "tested", tested)
    for hit in sorted(hits, key=lambda item: item[0])[:100]:
        print(hit)


def main():
    max_length = int(sys.argv[1])
    for p, q in ((3, 1), (3, 2), (3, 3), (1, 3), (2, 3)):
        for sign_class in (-1, 1):
            scan(p, q, sign_class, max_length)


if __name__ == "__main__":
    main()
