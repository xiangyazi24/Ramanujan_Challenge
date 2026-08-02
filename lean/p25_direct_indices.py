#!/usr/bin/env python3
"""Exact formal indices for the short direct Catalan 3F2 trajectories."""

import sympy as sp


n = sp.symbols("n", positive=True)
parameters = sp.symbols("x0 x1 x2 y0 y1")
x = parameters[:3]
y = parameters[3:]

c2 = (sum(x) + sum(y) - 2) / 2
c1 = (sum(x[i] * x[j] for i in range(3) for j in range(i + 1, 3))
      + (y[0] - 1) * (y[1] - 1)) / 2
c0 = x[0] * x[1] * x[2] / 2
theta = sp.Matrix([[0, 0, -c0], [1, 0, -c1], [0, 1, -c2]])
identity = sp.eye(3)
native = [identity + theta / value for value in x]
native += [identity + theta / (value - 1) for value in y]
positive = native[:3]
negative = [native[index].subs(x[index], x[index] - 1).inv()
            for index in range(3)]
for index in range(2):
    positive.append(native[index + 3].subs(y[index], y[index] + 1).inv())
    negative.append(native[index + 3])


def trajectory_matrix(trajectory):
    start = [sp.Rational(1, 2), sp.Rational(1, 2), 1,
             sp.Rational(3, 2), sp.Rational(3, 2)]
    position = [start[index] + n * trajectory[index] for index in range(5)]
    result = identity
    for index in reversed(range(5)):
        direction = 1 if trajectory[index] >= 0 else -1
        source = positive[index] if direction > 0 else negative[index]
        for _ in range(abs(trajectory[index])):
            substitution = dict(zip(parameters, position))
            result = result * source.subs(substitution)
            position[index] += direction
    return result.applyfunc(sp.cancel)


def coefficient_at_infinity(value, order):
    t = sp.symbols("t")
    transformed = sp.cancel(value.subs(n, 1 / t))
    return sp.cancel(sp.diff(transformed, t, order).subs(t, 0)
                     / sp.factorial(order))


def formal_data(trajectory):
    matrix = trajectory_matrix(trajectory)
    balance = sp.diag(1, n, n**2)
    balance_next = sp.diag(1, n + 1, (n + 1)**2)
    balanced = (balance * matrix * balance_next.inv()).applyfunc(sp.cancel)
    constant = balanced.applyfunc(lambda value: coefficient_at_infinity(value, 0))
    first = balanced.applyfunc(lambda value: coefficient_at_infinity(value, 1))
    roots = [1, 17 + 12 * sp.sqrt(2), 17 - 12 * sp.sqrt(2)]
    indices = []
    for root in roots:
        right = (constant - root * identity).nullspace()[0]
        left = (constant.T - root * identity).nullspace()[0]
        index = sp.factor((left.T * first * right)[0]
                          / (root * (left.T * right)[0]))
        indices.append(index)
    return matrix, constant, indices


for candidate in [
    (-2, 2, 0, 0, 0),
    (0, -2, 2, 0, 0),
    (-2, 0, 2, 0, 0),
    (0, 2, -2, 0, 0),
    (2, 0, -2, 0, 0),
    (0, 0, 0, -2, 2),
]:
    print("trajectory", candidate, flush=True)
    source, constant, indices = formal_data(candidate)
    print("det", sp.factor(source.det()), flush=True)
    print("constant", constant, flush=True)
    print("indices", indices, flush=True)
