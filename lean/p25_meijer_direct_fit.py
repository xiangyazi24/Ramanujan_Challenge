#!/usr/bin/env python3
"""Temporary numerical recovery of a raw CMF trajectory behind P2.5."""

import itertools
import sys

import numpy as np
import sympy as sp
from scipy.optimize import least_squares


def inverse3(matrix):
    """Explicit 3-by-3 inverse; much faster than Accelerate for tiny matrices."""
    a, b, c = matrix[0]
    d, e, f = matrix[1]
    g, h, i = matrix[2]
    determinant = a*(e*i-f*h)-b*(d*i-f*g)+c*(d*h-e*g)
    return np.array([
        [e*i-f*h, c*h-b*i, b*f-c*e],
        [f*g-d*i, a*i-c*g, c*d-a*f],
        [d*h-e*g, b*g-a*h, a*e-b*d],
    ], dtype=float)/determinant


def companion_from_poly(poly, theta):
    coefficients = sp.Poly(sp.expand(poly), theta).monic().all_coeffs()
    return sp.Matrix([
        [0, 0, -coefficients[3]],
        [1, 0, -coefficients[2]],
        [0, 1, -coefficients[1]],
    ])


def meijer_elementary(sign_class=-1):
    theta = sp.symbols("theta")
    axes = sp.symbols("a0 a1 a2 b0 b1 b2")
    aa, bb = axes[:3], axes[3:]
    companion = companion_from_poly(
        sign_class * sp.prod(theta-a+1 for a in aa)
        - sp.prod(theta-b for b in bb), theta)
    companion_value = sp.lambdify(axes, companion, modules="numpy")

    def raw(index, values):
        matrix = np.asarray(companion_value(*values), dtype=float)
        scalar = values[index]-1 if index < 3 else values[index]
        return matrix-scalar*np.eye(3)

    forms = {}
    for index in range(6):
        def elementary(*values, index=index, step=1):
            shifted = list(values)
            if index < 3 and step > 0:
                shifted[index] += 1
                return inverse3(raw(index, shifted))
            if index >= 3 and step < 0:
                shifted[index] -= 1
                return inverse3(raw(index, shifted))
            return raw(index, shifted)
        forms[index, 1] = lambda *values, elementary=elementary: elementary(*values, step=1)
        forms[index, -1] = lambda *values, elementary=elementary: elementary(*values, step=-1)
    return axes, forms


def pfq_elementary(z_value=-1):
    theta = sp.symbols("theta")
    axes = sp.symbols("x0 x1 x2 y0 y1")
    xx, yy = axes[:3], axes[3:]
    companion = companion_from_poly(
        theta*sp.prod(theta+y-1 for y in yy)
        - z_value*sp.prod(theta+x for x in xx), theta)
    companion_value = sp.lambdify(axes, companion, modules="numpy")

    def raw(index, values):
        matrix = np.asarray(companion_value(*values), dtype=float)
        denominator = values[index] if index < 3 else values[index]-1
        return np.eye(3)+matrix/denominator

    forms = {}
    for index in range(5):
        def elementary(*values, index=index, step=1):
            shifted = list(values)
            if index < 3 and step < 0:
                shifted[index] -= 1
                return inverse3(raw(index, shifted))
            if index >= 3 and step > 0:
                shifted[index] += 1
                return inverse3(raw(index, shifted))
            return raw(index, shifted)
        forms[index, 1] = lambda *values, elementary=elementary: elementary(*values, step=1)
        forms[index, -1] = lambda *values, elementary=elementary: elementary(*values, step=-1)
    return axes, forms


def derived_pfq_elementary(z_value=-1):
    """The older hypergeometric-derived 3F2 CMF basis in ResearchTools."""
    axes = sp.symbols("x0 x1 x2 y0 y1")

    def base(values):
        xx = values[:3]
        yy = values[3:]
        sx = sum(xx)
        sy = sum(yy)
        tx = xx[0]*xx[1]+xx[0]*xx[2]+xx[1]*xx[2]
        ty = yy[0]*yy[1]
        px = xx[0]*xx[1]*xx[2]
        z = z_value
        return np.array([
            [0, 0, px/((1-z)*z)],
            [z, 1, ((tx+sx+1)*z-ty)/((1-z)*z)],
            [0, z, ((sx+1)*z+sy+1)/(1-z)],
        ], dtype=float)

    def raw(index, values):
        matrix = base(values)
        if index < 3:
            return np.eye(3)+matrix/values[index]
        return np.eye(3)-matrix/(values[index]+1)

    forms = {}
    for index in range(5):
        forms[index, 1] = lambda *values, index=index: raw(index, values)

        def negative(*values, index=index):
            shifted = list(values)
            shifted[index] -= 1
            return inverse3(raw(index, shifted))

        forms[index, -1] = negative
    return axes, forms


def target(n):
    return np.array([
        [(-2*n-5)*(n+3)**2*(136*n**4+1424*n**3+5548*n**2+9551*n+6141),
         384*n**6+6384*n**5+44168*n**4+162698*n**3+336377*n**2+369933*n+169011,
         -(480*n**4+4980*n**3+19210*n**2+32690*n+20730)],
        [(n+2)**2*(n+3)**2*(4*n+10)*(48*n**3+386*n**2+1017*n+879),
         (n+2)**2*(-272*n**5-3848*n**4-21732*n**3-61184*n**2-85761*n-47808),
         (n+2)**2*(320*n**3+2540*n**2+6610*n+5640)],
        [(-4*n-10)*(n+2)**2*(n+3)**2*(32*n**4+302*n**3+1037*n**2+1530*n+813),
         (n+2)**2*(192*n**6+2984*n**5+19116*n**4+64452*n**3+120256*n**2+117279*n+46476),
         (n+2)**2*(-16*n**5-408*n**4-2912*n**3-8884*n**2-12254*n-6240)],
    ], dtype=float)


def normalized(matrix):
    scale = np.linalg.norm(matrix)
    return matrix/scale if np.isfinite(scale) and scale > 1e-200 else matrix


def trajectory_value(forms, start, direction, n):
    position = np.asarray(start, dtype=float) + n*np.asarray(direction, dtype=float)
    result = np.eye(3)
    with np.errstate(all="ignore"):
        for index in reversed(range(len(direction))):
            step = 1 if direction[index] >= 0 else -1
            for _ in range(abs(direction[index])):
                result = result @ np.asarray(forms[index, step](*position), dtype=float)
                position[index] += step
    return result


def variants(matrix):
    yield "plain", matrix
    yield "transpose", matrix.T
    try:
        yield "dual", inverse3(matrix).T
        yield "inverse", inverse3(matrix)
    except np.linalg.LinAlgError:
        pass


def fit_family(name, axes, forms, directions, attempts=16):
    rng = np.random.default_rng(2505)
    targets = [normalized(target(k)) for k in range(4)]
    signs = np.diag([1, -1, 1])
    target_sets = {
        "raw": targets,
        "positive": [normalized(-signs@target(k)@signs) for k in range(4)],
    }
    for direction in directions:
        for transform_name in ("plain", "transpose", "dual", "inverse"):
            for target_name, desired in target_sets.items():
                def residual(values):
                    output = []
                    try:
                        for k in range(4):
                            raw = trajectory_value(forms, values, direction, k)
                            transformed = dict(variants(raw))[transform_name]
                            current = normalized(transformed)
                            if not np.isfinite(current).all():
                                raise ValueError
                            scalar_sign = 1 if np.vdot(current, desired[k]) >= 0 else -1
                            output.extend((current-scalar_sign*desired[k]).ravel())
                        return np.asarray(output)
                    except (ValueError, ZeroDivisionError, np.linalg.LinAlgError):
                        return np.full(36, 1e3)
                best = None
                seeds = [np.zeros(len(axes))]
                seeds += [rng.uniform(-6, 6, size=len(axes)) for _ in range(attempts)]
                for initial in seeds:
                    fit = least_squares(residual, initial, max_nfev=700,
                                        xtol=1e-12, ftol=1e-12, gtol=1e-12)
                    score = np.linalg.norm(fit.fun)
                    if best is None or score < best[0]:
                        best = (score, fit.x)
                print(name, direction, transform_name, target_name, best, flush=True)
                if best[0] < 1e-8:
                    return


def main():
    family = sys.argv[1] if len(sys.argv) > 1 else "meijer-minus"
    if family.startswith("meijer"):
        sign_class = -1 if family.endswith("minus") else 1
        axes, forms = meijer_elementary(sign_class)
        directions = [(-2, 0, 2, 0, 0, 0), (0, 0, 0, -2, 0, 2)]
        fit_family(family, axes, forms, directions, attempts=8)
    elif family == "pfq-minus-one":
        axes, forms = pfq_elementary(-1)
        seed_directions = [(-2, 0, 2, 0, 0), (0, 0, 0, -2, 2)]
        directions = sorted(set(itertools.chain.from_iterable(
            set(itertools.permutations(direction[:3]))
            for direction in seed_directions[:1])))
        directions = [direction+(0, 0) for direction in directions]
        directions += [(0, 0, 0, -2, 2), (0, 0, 0, 2, -2)]
        fit_family(family, axes, forms, directions, attempts=8)
    elif family == "derived-pfq-minus-one":
        axes, forms = derived_pfq_elementary(-1)
        directions = []
        for pair in itertools.permutations(range(3), 2):
            direction = [0]*5
            direction[pair[0]] = -2
            direction[pair[1]] = 2
            directions.append(tuple(direction))
        directions += [(0, 0, 0, -2, 2), (0, 0, 0, 2, -2)]
        fit_family(family, axes, forms, directions, attempts=3)
    else:
        raise SystemExit(f"unknown family: {family}")


if __name__ == "__main__":
    main()
