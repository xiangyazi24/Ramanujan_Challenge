#!/usr/bin/env sage
"""Probe the gauge freedom in the same-moment logarithmic reduction.

For

    q_i = theta_i(F) - F,

the logarithmic weight map is H |-> (sum H_i q_i)/F.  Two fields with
the same weight differ by an element of syz(q_x,q_y,q_z).  This script
prints a polynomial syzygy basis and verifies it exactly.
"""

R.<x, y, z> = PolynomialRing(QQ, order="degrevlex")
F = (
    (1 + x) * (1 + y) * (1 + z)
    * ((1 + y) * (1 + z) + x * y * z)
)


def theta(polynomial, variable):
    return variable * polynomial.derivative(variable)


q = tuple(theta(F, variable) - F for variable in (x, y, z))
kernel = list(R.ideal(q).syzygy_module())

print("Q_FACTORS")
for entry in q:
    print(factor(entry))

print("KERNEL_SIZE", len(kernel))
for index, row in enumerate(kernel):
    row = tuple(R(entry) for entry in row)
    assert sum(row[coordinate] * q[coordinate] for coordinate in range(3)) == 0
    print("KERNEL", index, row)
    print("DEGREES", tuple(entry.degree() for entry in row))

print("Q32_LOG_WEIGHT_KERNEL_PROBE=PASS")
