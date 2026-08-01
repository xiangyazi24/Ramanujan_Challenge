#!/usr/bin/env python3
"""Temporary pure-Numpy spectral scan for a rank-three Catalan Meijer G CMF."""

import itertools
import math
import sys

import numpy as np
import sympy as sp


def canonical_vectors(count, max_length):
    values = range(-max_length, max_length + 1)
    for vector in itertools.combinations_with_replacement(values, count):
        if sum(abs(value) for value in vector) <= max_length:
            yield vector


def axis_forms():
    theta = sp.symbols("theta")
    axes = sp.symbols("a0 a1 a2 b0 b1 b2")
    aa, bb = axes[:3], axes[3:]
    polynomial = sp.Poly(
        -sp.prod(theta-a+1 for a in aa)-sp.prod(theta-b for b in bb),
        theta,
    ).monic()
    coefficients = polynomial.all_coeffs()
    companion = sp.Matrix([
        [0, 0, -coefficients[3]],
        [1, 0, -coefficients[2]],
        [0, 1, -coefficients[1]],
    ])
    eye = sp.eye(3)
    native = [companion-(a-1)*eye for a in aa]
    native += [companion-b*eye for b in bb]
    positive = []
    negative = []
    for index, axis in enumerate(axes):
        if index < 3:
            negative.append(native[index])
            positive.append(native[index].subs(axis, axis+1).inv())
        else:
            positive.append(native[index])
            negative.append(native[index].inv().subs(axis, axis-1))
    return axes, {
        (index, sign): sp.lambdify(
            axes, positive[index] if sign > 0 else negative[index],
            modules="numpy",
        )
        for index in range(6) for sign in (-1, 1)
    }


def main():
    max_length = int(sys.argv[1])
    axes, forms = axis_forms()
    base = np.arange(1, len(axes)+1, dtype=float)/7.0+0.13

    def step_matrix(trajectory, sample_n=100000.0):
        position = base+sample_n*np.array(trajectory, dtype=float)
        result = np.eye(3)
        with np.errstate(all="ignore"):
            for index in reversed(range(6)):
                direction = 1 if trajectory[index] >= 0 else -1
                for _ in range(abs(trajectory[index])):
                    result = result@np.asarray(
                        forms[index, direction](*position), dtype=float)
                    position[index] += direction
        scale = np.max(np.abs(result))
        return result/scale if scale else result

    avectors = list(canonical_vectors(3, max_length))
    bvectors = list(canonical_vectors(3, max_length))
    hits = []
    tested = 0
    for aa in avectors:
        alength = sum(abs(value) for value in aa)
        for bb in bvectors:
            length = alength+sum(abs(value) for value in bb)
            if length == 0 or length > max_length:
                continue
            trajectory = aa+bb
            tested += 1
            try:
                eigenvalues = np.linalg.eigvals(step_matrix(trajectory))
                e1 = sum(eigenvalues)
                e2 = sum(eigenvalues[i]*eigenvalues[j]
                         for i in range(3) for j in range(i+1, 3))
                e3 = np.prod(eigenvalues)
                if not np.isfinite(eigenvalues).all() or abs(e3) < 1e-100:
                    continue
                invariant1 = float(np.real_if_close(e1*e2/e3))
                invariant2 = float(np.real_if_close(e1**3/e3))
                score = abs(math.log(abs(invariant1/1225)))+abs(
                    math.log(abs(invariant2/42875)))
                if score < 0.02:
                    hits.append((score, trajectory, invariant1, invariant2,
                                 eigenvalues))
            except (TypeError, ValueError, np.linalg.LinAlgError):
                pass
    print("tested", tested)
    for hit in sorted(hits, key=lambda item: item[0])[:200]:
        print(hit)


if __name__ == "__main__":
    main()
