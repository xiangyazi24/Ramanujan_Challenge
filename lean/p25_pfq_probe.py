#!/usr/bin/env python3
"""Search asymptotic Catalan 3F2 CMF directions for the P2.5 spectrum."""

import itertools

import numpy as np


target = 17 + 12 * np.sqrt(2)


def roots_for_direction(values):
    a = values[:3]
    b = values[3:]
    # z=-1: t(t+b0)(t+b1) + (t+a0)(t+a1)(t+a2) = 0.
    first = np.polymul([1, 0], np.polymul([1, b[0]], [1, b[1]]))
    second = np.polymul(np.polymul([1, a[0]], [1, a[1]]), [1, a[2]])
    theta_roots = np.roots(first + second)
    eigenvalues = []
    for theta in theta_roots:
        value = 1.0 + 0j
        for step in a:
            value *= (1 + theta / step) ** step
        for step in b:
            value *= (1 + theta / step) ** (-step)
        eigenvalues.append(value)
    return sorted(eigenvalues, key=abs, reverse=True)


candidates = []
domain = [value for value in range(-6, 7) if value]
for xs in itertools.combinations_with_replacement(domain, 3):
    for ys in itertools.combinations_with_replacement(domain, 2):
        values = xs + ys
        if sum(abs(value) for value in values) > 20:
            continue
        try:
            eigenvalues = roots_for_direction(values)
            ratios = [abs(eigenvalues[0] / eigenvalues[1]),
                      abs(eigenvalues[1] / eigenvalues[2])]
            score = max(abs(np.log(ratio / target)) for ratio in ratios)
            score += max(abs(value.imag) / max(1, abs(value.real))
                         for value in eigenvalues)
        except (FloatingPointError, OverflowError, ZeroDivisionError):
            continue
        if not np.isfinite(score):
            continue
        candidates.append((score, values, eigenvalues, ratios))
        candidates.sort(key=lambda item: item[0])
        del candidates[40:]

for score, values, eigenvalues, ratios in candidates:
    print(score, values, eigenvalues, ratios)
