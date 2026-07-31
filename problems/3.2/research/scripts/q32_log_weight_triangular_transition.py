#!/usr/bin/env python3
"""Full transition graph in the zero-block triangular logarithmic gauge.

The y and x source rows are gauged so that every zero-shift ideal output
vanishes.  The z source has its diagonal zero-shift removed; its two
remaining zero-shift outputs point only to y,x.  This script verifies
those claims by exact Laurent division, emits every transition, and
tests strict affine-linear ranking feasibility.
"""

from pathlib import Path

import sympy as sp


HELPER = Path(__file__).with_name("q32_log_weight_zero_block_gauge.py")
prefix = HELPER.read_text().split("original =")[0]
namespace = {}
namespace["__file__"] = str(HELPER)
exec(prefix, namespace)

x, y, z = namespace["x"], namespace["y"], namespace["z"]
a, b, c = namespace["a"], namespace["b"], namespace["c"]
fields = namespace["fields"]
pair_fields = namespace["pair_fields"]
theta = namespace["theta"]
canonical_division = namespace["canonical_division"]
laurent_ledger = namespace["laurent_ledger"]

H_y, H_x, H_z = fields
K_yx, K_yz, K_xz = pair_fields
g_y = y - z
g_x = x - 2 * z - 1
g_z = z**2 - sp.Rational(1, 2)


def scale_field(scalar, field):
    return tuple(sp.expand(scalar * entry) for entry in field)


def add_fields(*summands):
    return tuple(
        sp.expand(sum(field[index] for field in summands))
        for index in range(3)
    )


K_y_row = add_fields(
    scale_field(1 + y + sp.Rational(1, 2) * z**-3 - z**-1, K_yx),
    scale_field(2 - 2 * x, K_yz),
    scale_field(-1 + x + sp.Rational(1, 2) * z**-4, K_xz),
)
K_x_row = add_fields(
    scale_field(
        sp.Rational(7, 2)
        + sp.Rational(3, 2) * y
        + sp.Rational(1, 2) * z**-3
        + 2 * z**-2
        - sp.Rational(7, 2) * z**-1,
        K_yx,
    ),
    scale_field(5 - 6 * x, K_yz),
    scale_field(
        -3 + 3 * x
        + sp.Rational(7, 2) * z**-4
        - 2 * z**-3,
        K_xz,
    ),
)
K_z_diagonal = add_fields(
    scale_field(-sp.Rational(1, 2), K_yz),
    scale_field(
        -sp.Rational(19, 4) * z**-3
        + 3 * z**-2
        - sp.Rational(1, 2) * z**-1,
        K_xz,
    ),
)

triangular_fields = (
    add_fields(H_y, K_y_row),
    add_fields(H_x, K_x_row),
    add_fields(H_z, K_z_diagonal),
)

# Exact weight verification from the defining q-vector.
F = (
    (1 + x) * (1 + y) * (1 + z)
    * ((1 + y) * (1 + z) + x * y * z)
)
q = tuple(theta(F, variable) - F for variable in (x, y, z))
weights = (g_y, g_x, g_z)
for field, weight in zip(triangular_fields, weights):
    assert sp.cancel(
        sum(field[index] * q[index] for index in range(3))
        - weight * F
    ) == 0


names = ("y", "x", "z")
edges = []
zero_block = {}
for source, field in enumerate(triangular_fields):
    beta = sum(
        theta(field[index], (x, y, z)[index])
        for index in range(3)
    )
    coefficient = sp.expand(
        a * field[0] + b * field[1] + c * field[2] + beta
    )
    quotients, remainder = canonical_division(coefficient)
    print("SOURCE", names[source])
    for target, quotient in enumerate(quotients):
        ledger = laurent_ledger(quotient)
        constant = ledger.pop((0, 0, 0), 0)
        zero_block[(source, target)] = sp.factor(constant)
        print(
            " TO", names[target],
            "ZERO", sp.factor(constant),
            "COUNT", len(ledger),
        )
        for shift, value in sorted(ledger.items()):
            value = sp.factor(value)
            print("  EDGE", shift, value)
            edges.append((source, target, shift, value))
    print(" REMAINDER", laurent_ledger(remainder))

for source in (0, 1):
    assert all(zero_block[(source, target)] == 0 for target in range(3))
assert zero_block[(2, 2)] == 0
for target in (0, 1):
    if zero_block[(2, target)] != 0:
        edges.append((2, target, (0, 0, 0), zero_block[(2, target)]))

print("ZERO_BLOCK")
for source in range(3):
    print(
        names[source],
        tuple(zero_block[(source, target)] for target in range(3)),
    )
print("EDGE_COUNT", len(edges))


from scipy.optimize import linprog

# Variables wx,wy,wz,t_x,t_z; t_y=0.
rows = []
for source, target, shift, _ in edges:
    row = [
        float(shift[0]),
        float(shift[1]),
        float(shift[2]),
        0.0,
        0.0,
    ]
    if target == 1:
        row[3] += 1
    elif target == 2:
        row[4] += 1
    if source == 1:
        row[3] -= 1
    elif source == 2:
        row[4] -= 1
    rows.append(row)

result = linprog(
    c=[0.0] * 5,
    A_ub=rows,
    b_ub=[-1.0] * len(rows),
    bounds=[(None, None)] * 5,
    method="highs",
)
print("STRICT_LINEAR_RANKING", result.success)
if result.success:
    print(" RANKING_FLOAT", tuple(result.x))


def short_zero_cycle():
    # Length one.
    for edge in edges:
        if edge[0] == edge[1] and edge[2] == (0, 0, 0):
            return [edge]

    # Length two.
    reverse = {}
    for edge in edges:
        reverse.setdefault((edge[0], edge[1], edge[2]), edge)
    for edge in edges:
        source, target, shift, _ = edge
        opposite = tuple(-entry for entry in shift)
        partner = reverse.get((target, source, opposite))
        if partner is not None:
            return [edge, partner]

    # Length three by hashing the closing edge.
    for first in edges:
        for second in edges:
            if first[1] != second[0]:
                continue
            needed = tuple(
                -first[2][index] - second[2][index]
                for index in range(3)
            )
            third = reverse.get((second[1], first[0], needed))
            if third is not None:
                return [first, second, third]

    # Length four: match two directed paths with opposite total shifts.
    paths = {}
    for first in edges:
        for second in edges:
            if first[1] != second[0]:
                continue
            total = tuple(
                first[2][index] + second[2][index]
                for index in range(3)
            )
            paths.setdefault((first[0], second[1], total), (first, second))
    for (source, target, total), path in paths.items():
        opposite = tuple(-entry for entry in total)
        partner = paths.get((target, source, opposite))
        if partner is not None:
            return list(path + partner)
    return None


cycle = short_zero_cycle()
print("ZERO_TOTAL_SHIFT_CYCLE", cycle is not None)
if cycle is not None:
    for edge in cycle:
        print(
            " CYCLE",
            names[edge[0]], "->", names[edge[1]],
            edge[2], edge[3],
        )

print("Q32_LOG_WEIGHT_TRIANGULAR_TRANSITION=PASS")
