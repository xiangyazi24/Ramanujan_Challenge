#!/usr/bin/env sage -python
"""Temporary exact rational-gauge recovery for Problem 2.5."""

import sys

from sage.all import QQ, PolynomialRing, matrix, vector

sys.path.insert(0, "/Users/huangx/Library/SageMath-10-9/lib/python3.14/site-packages")
from ore_algebra import OreAlgebra


R = QQ["n"]
n = R.gen()
K = R.fraction_field()
OAK = OreAlgebra(K, names=("Sn",))


def sigma(value, amount=1):
    return value(n=n + amount)


def sigma_matrix(value, amount=1):
    return value.apply_map(lambda entry: sigma(entry, amount))


def flatten(value):
    return vector(K, [value[row, column] for row in range(3) for column in range(3)])


def challenge_matrix():
    entries = [
        (-2*n-5)*(n+3)**2*(136*n**4+1424*n**3+5548*n**2+9551*n+6141),
        384*n**6+6384*n**5+44168*n**4+162698*n**3+336377*n**2+369933*n+169011,
        -(480*n**4+4980*n**3+19210*n**2+32690*n+20730),
        (n+2)**2*(n+3)**2*(4*n+10)*(48*n**3+386*n**2+1017*n+879),
        (n+2)**2*(-272*n**5-3848*n**4-21732*n**3-61184*n**2-85761*n-47808),
        (n+2)**2*(320*n**3+2540*n**2+6610*n+5640),
        (-4*n-10)*(n+2)**2*(n+3)**2*(32*n**4+302*n**3+1037*n**2+1530*n+813),
        (n+2)**2*(192*n**6+2984*n**5+19116*n**4+64452*n**3+120256*n**2+117279*n+46476),
        (n+2)**2*(-16*n**5-408*n**4-2912*n**3-8884*n**2-12254*n-6240),
    ]
    delta = -2*(n+2)**2*(n+3)**2*(2*n+5)*(2*n+7)**2
    return matrix(K, 3, 3, entries) / delta


def pfq_matrix():
    d0 = (n+1)**3*(2*n+3)**2
    d1 = (n+1)**3*(2*n+1)*(2*n+3)**2
    entries = [
        (256*n**5+7936*n**4+24608*n**3+28688*n**2+14561*n+2726)/(64*d0),
        5*(4*n+3)*(8*n+5)*(128*n**3+316*n**2+237*n+59)/(128*d0),
        5*(4*n+3)*(8*n+5)*(384*n**4+1184*n**3+1356*n**2+731*n+150)/(256*d0),
        (8*n+5)*(384*n**4+1824*n**3+2908*n**2+1831*n+378)/(32*d1),
        (4*n+3)*(2176*n**5+11296*n**4+20944*n**3+17592*n**2+6875*n+1017)/(64*d1),
        (4*n+3)*(8*n+5)*(768*n**5+4032*n**4+7664*n**3+6896*n**2+3093*n+522)/(128*d1),
        (8*n+5)*(128*n**3+332*n**2+257*n+58)/(16*d1),
        (4*n+3)*(8*n+5)*(96*n**3+224*n**2+149*n+31)/(32*d1),
        (4*n+3)*(2176*n**5+7712*n**4+10768*n**3+7736*n**2+2831*n+402)/(64*d1),
    ]
    return matrix(K, 3, 3, entries)


def direct_pfq_axis_matrices():
    """Exact 3F2 CMF matrices at z=-1 in the theta basis."""
    parameter_ring = PolynomialRing(K, names=("x0,x1,x2,y0,y1"))
    parameter_field = parameter_ring.fraction_field()
    params = parameter_field.gens()
    xx = params[:3]
    yy = params[3:]
    c2 = (sum(xx) + sum(yy) - 2) / 2
    c1 = (sum(xx[i] * xx[j] for i in range(3) for j in range(i + 1, 3))
          + (yy[0] - 1) * (yy[1] - 1)) / 2
    c0 = xx[0] * xx[1] * xx[2] / 2
    companion = matrix(parameter_field, [
        [0, 0, -c0],
        [1, 0, -c1],
        [0, 1, -c2],
    ])
    identity = matrix.identity(parameter_field, 3)
    native = [identity + companion / value for value in xx]
    native += [identity + companion / (value - 1) for value in yy]
    positive = native[:3]
    negative = [native[index](**{str(xx[index]): xx[index] - 1}).inverse()
                for index in range(3)]
    for index in range(2):
        axis = yy[index]
        positive.append(native[index + 3](**{str(axis): axis + 1}).inverse())
        negative.append(native[index + 3])
    return params, positive, negative


DIRECT_PFQ_PARAMS, DIRECT_PFQ_POSITIVE, DIRECT_PFQ_NEGATIVE = direct_pfq_axis_matrices()


def direct_pfq_matrix(trajectory, offset=0):
    """One trajectory step based at the direct Catalan 3F2 parameters."""
    start = [K(1)/2, K(1)/2, 1, K(3)/2, K(3)/2]
    position = [start[index] + (n + offset) * trajectory[index]
                for index in range(5)]
    result = matrix.identity(K, 3)
    for index in reversed(range(5)):
        direction = 1 if trajectory[index] >= 0 else -1
        source = DIRECT_PFQ_POSITIVE[index] if direction > 0 else DIRECT_PFQ_NEGATIVE[index]
        for _ in range(abs(trajectory[index])):
            substitution = {DIRECT_PFQ_PARAMS[j]: position[j] for j in range(5)}
            result *= source.subs(substitution)
            position[index] += direction
    return result.apply_map(K)


def meijer_axis_matrices():
    """Exact CMF matrices for G^1,3_3,3 at z=1.

    Parameters are ordered a0,a1,a2,b0,b1,b2.  This duplicates the small
    DFinite construction in ramanujantools, but stays inside Sage so that the
    rational-solution code uses one compatible CAS runtime.
    """
    parameter_ring = PolynomialRing(K, names=("a0,a1,a2,b0,b1,b2"))
    parameter_field = parameter_ring.fraction_field()
    params = parameter_field.gens()
    aa = params[:3]
    bb = params[3:]

    # p-m-n = -1, hence the differential polynomial is
    # -prod(theta-a_i+1)-prod(theta-b_j).  Its monic coefficients are
    # written explicitly to avoid a nested polynomial-ring coercion.
    ar = [value - 1 for value in aa]
    c2 = -(sum(ar) + sum(bb)) / 2
    c1 = (sum(ar[i] * ar[j] for i in range(3) for j in range(i + 1, 3))
          + sum(bb[i] * bb[j] for i in range(3) for j in range(i + 1, 3))) / 2
    c0 = -(ar[0] * ar[1] * ar[2] + bb[0] * bb[1] * bb[2]) / 2
    companion = matrix(parameter_field, [
        [0, 0, -c0],
        [1, 0, -c1],
        [0, 1, -c2],
    ])
    identity = matrix.identity(parameter_field, 3)

    raw = []
    for index, axis in enumerate(params):
        is_a = index < 3
        local_index = index if is_a else index - 3
        sign = 1 if is_a else (-1 if local_index < 1 else 1)
        multiplier = axis - 1 if is_a else axis
        raw.append(sign * (companion - multiplier * identity))

    positive = []
    negative = []
    for index, axis in enumerate(params):
        if index < 3:  # a axes are natively decrementing
            negative.append(raw[index])
            positive.append(raw[index](**{str(axis): axis + 1}).inverse())
        else:          # b axes are natively incrementing
            positive.append(raw[index])
            negative.append(raw[index].inverse()(**{str(axis): axis - 1}))
    return params, positive, negative


MEIJER_PARAMS, MEIJER_POSITIVE, MEIJER_NEGATIVE = meijer_axis_matrices()


def meijer_matrix(trajectory, offset=0):
    """One exact trajectory step based at Catalan's Meijer-G parameters."""
    start = [K(1)/2, K(1)/2, 0, 0, -K(1)/2, -K(1)/2]
    position = [start[index] + (n + offset) * trajectory[index]
                for index in range(6)]
    result = matrix.identity(K, 3)
    # Any path is valid by conservation.  This fixed path mirrors the numeric
    # spectral scan and avoids singular intermediate points where possible.
    for index in reversed(range(6)):
        direction = 1 if trajectory[index] >= 0 else -1
        source = MEIJER_POSITIVE[index] if direction > 0 else MEIJER_NEGATIVE[index]
        for _ in range(abs(trajectory[index])):
            substitution = {MEIJER_PARAMS[j]: position[j] for j in range(6)}
            result *= source.subs(substitution)
            position[index] += direction
    return result.apply_map(K)


def recover(name, left, right):
    print(name, flush=True)
    basis = []
    for index in range(9):
        value = matrix(K, 3, 3, 0)
        value[index // 3, index % 3] = 1
        basis.append(value)
    left_inverse = left.inverse()
    transition = matrix(
        K,
        9,
        9,
        lambda row, column: flatten(left_inverse * basis[column] * right)[row],
    )
    print("transition", flush=True)
    for seed_index in range(9):
        print("seed", seed_index, flush=True)
        seed = [K.zero()] * 9
        seed[seed_index] = K.one()
        rows = [vector(K, seed)]
        for index in range(9):
            rows.append(vector(K, [sigma(entry) for entry in rows[-1]]) * transition)
        krylov = matrix(K, rows[:9])
        if krylov.det() == 0:
            print("noncyclic", flush=True)
            continue
        relation = krylov.transpose().solve_right(-rows[9])
        operator = OAK(list(relation) + [1]).normalize()
        solutions = operator.rational_solutions()
        print("solutions", len(solutions), flush=True)
        for solution in solutions:
            scalar = solution[0]
            rhs = vector(K, [sigma(scalar, shift) for shift in range(9)])
            gauge_vector = krylov.solve_right(rhs)
            gauge = matrix(K, 3, 3, list(gauge_vector))
            assert sigma_matrix(gauge) == left_inverse * gauge * right
            print("gauge_rank", gauge.rank(), flush=True)
            if gauge.det() == 0:
                continue
            print("gauge_det", gauge.det().factor(), flush=True)
            for row in gauge.rows():
                print([entry.factor() for entry in row], flush=True)
            return gauge
    return None


if __name__ == "__main__":
    A = challenge_matrix()
    print("challenge det", A.det().factor(), flush=True)
    for trajectory in [
        (-2, 2, 0, 0, 0), (2, -2, 0, 0, 0),
        (0, 0, 0, -2, 2), (0, 0, 0, 2, -2),
    ]:
        B = direct_pfq_matrix(trajectory, 0)
        print("pfq", trajectory, "det", B.det().factor(), flush=True)
        gauge = recover("challenge_to_pfq", A, B)
        if gauge is not None:
            raise SystemExit(0)
