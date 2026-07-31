#!/usr/bin/env sage
"""Exact signed-binomial transforms of the first-cell ray recurrences.

For a nonzero Newton ray ``kappa`` put

    R_kappa(M,r) = [X^((M-r)kappa)] Lambda^M,
    B_kappa(M,L) = sum_{r=0}^L (-1)^r binom(L,r) R_kappa(M,r).

The script derives a certified telescoper for ``R_kappa`` and transports
it through the exact ordinary-generating-function identity

    sum_L B_kappa(M,L) z^L
      = (1-z)^(-1) A_kappa(M,-z/(1-z)),

where ``A_kappa(M,u)=sum_r R_kappa(M,r)u^r``.  The resulting operators
are genuine recurrences for the constant-term ray packets, not
recurrences in a Jacobian quotient.

The fourteen y/z-orbits form an exact finite first-order system for the
summed terminal packet.  We print and audit every endpoint factor; this
is the denominator ledger needed at ``M=p+r-1, L=r``.
"""

from math import comb

from ore_algebra import *
from ore_algebra import nullspace
from ore_algebra.ore_algebra import OreAlgebra_generic


def _associated_commutative_algebra(self):
    try:
        return self._commutative_ring
    except AttributeError:
        self._commutative_ring = PolynomialRing(
            self.base_ring(), self.variable_names()
        )
        return self._commutative_ring


OreAlgebra_generic.associated_commutative_algebra = (
    _associated_commutative_algebra
)

_original_kronecker = nullspace.kronecker


def _qq_safe_kronecker(subsolver, presolver=None):
    return nullspace.clear(_original_kronecker(subsolver, presolver))


nullspace.kronecker = _qq_safe_kronecker

R.<M, r, t> = QQ[]
OA.<Sr, St> = OreAlgebra(R)
d = M - r

RAY_CLASSES = (
    ((-1, -1, -1), 1),
    ((-1, -1, 0), 2),
    ((-1, -1, 1), 2),
    ((-1, 0, 0), 1),
    ((-1, 0, 1), 2),
    ((-1, 1, 1), 1),
    ((0, -1, -1), 1),
    ((0, -1, 0), 2),
    ((0, -1, 1), 2),
    ((0, 0, 1), 2),
    ((0, 1, 1), 1),
    ((1, 0, 0), 1),
    ((1, 0, 1), 2),
    ((1, 1, 1), 1),
)
assert sum(multiplicity for _, multiplicity in RAY_CLASSES) == 21


def shifted_binomial_ratio(upper, lower, shift):
    if shift == 1:
        return (upper - lower) / (lower + 1)
    if shift == -1:
        return lower / (upper - lower + 1)
    assert shift == 0
    return R.one()


def ray_operator(point):
    """Return a certified r-telescoper for one coefficient ray."""

    u, v, w = point
    kx = t - u * d
    upper = 2 * M - t
    ky = M - v * d
    kz = M - w * d
    t_ratio = (
        (M - t)
        / (t + 1)
        * (M - kx)
        / (kx + 1)
        * (upper - ky)
        / upper
        * (upper - kz)
        / upper
    )
    r_ratio = (
        shifted_binomial_ratio(M, kx, u)
        * shifted_binomial_ratio(upper, ky, v)
        * shifted_binomial_ratio(upper, kz, w)
    )
    ideal = OA.ideal(
        [
            t_ratio.denominator() * St - t_ratio.numerator(),
            r_ratio.denominator() * Sr - r_ratio.numerator(),
        ]
    )
    telescopers, certificates = ideal.ct(
        St - 1,
        certificates=True,
        early_termination=True,
        iteration_limit=16,
    )
    assert len(telescopers) == len(certificates) == 1
    return telescopers[0], certificates[0]


# Use M as a constant parameter while converting recurrence -> ODE ->
# rational pullback -> recurrence.
PM.<MM> = QQ[]
Rr.<rr> = PM[]
Ar.<Srr> = OreAlgebra(Rr)
Rz.<zz> = PM[]
Az.<Dzz> = OreAlgebra(Rz)
Kz = Rz.fraction_field()
AzK = OreAlgebra(
    Kz,
    ("DzzK", lambda value: value,
     lambda value: value.derivative(zz)),
)
DzzK = AzK.gen()
RL.<LL> = PM[]
AL.<SLL> = OreAlgebra(RL)
to_parameter_recurrence = R.hom([MM, rr, 0], Rr)


def signed_binomial_transform(operator):
    """Transport an r-recurrence to its signed binomial transform."""

    recurrence = Ar(
        sum(
            to_parameter_recurrence(operator[shift]) * Srr**shift
            for shift in range(operator.order() + 1)
        )
    )
    generating_ode = recurrence.to_D(Az)
    # Direct chain-rule pullback.  The generic composition routine uses
    # a costly nullspace solver and is fragile over the parameter ring
    # QQ[M].  If phi=-z/(1-z), then D_u acts on A(phi(z)) as
    # (phi'(z))^(-1) D_z.  Ore multiplication automatically differentiates
    # the variable coefficients in successive powers.
    phi = Kz(-zz / (1 - zz))
    pulled_derivative = Kz.one() / phi.derivative(zz) * DzzK
    pulled_ode_fraction = sum(
        Kz(generating_ode[index](zz=phi))
        * pulled_derivative**index
        for index in range(generating_ode.order() + 1)
    )
    # B(z)=C(z)/(1-z), so C(z)=(1-z)B(z).  Right multiplication
    # by 1-z is the exact gauge transform L(C)=L((1-z)B).
    transformed_fraction = pulled_ode_fraction * Kz(1 - zz)
    common_denominator = lcm(
        Kz(transformed_fraction[index]).denominator()
        for index in range(transformed_fraction.order() + 1)
    )
    transformed_ode = Az(
        sum(
            Rz(
                common_denominator
                * Kz(transformed_fraction[index])
            )
            * Dzz**index
            for index in range(transformed_fraction.order() + 1)
        )
    )
    transformed = transformed_ode.to_S(AL)
    if transformed.valuation():
        transformed = AL(
            sum(
                transformed[index + transformed.valuation()] * SLL**index
                for index in range(
                    transformed.order() - transformed.valuation() + 1
                )
            )
        )
    return transformed.primitive_part()


def C(n, k):
    return comb(n, k) if 0 <= k <= n else 0


def ray_value(moment, residue, point):
    u, v, w = point
    node = moment - residue
    return sum(
        C(moment, index)
        * C(moment, index - node * u)
        * C(2 * moment - index, moment - node * v)
        * C(2 * moment - index, moment - node * w)
        for index in range(moment + 1)
    )


def ray_packet(moment, order, point):
    return sum(
        (-1) ** residue
        * C(order, residue)
        * ray_value(moment, residue, point)
        for residue in range(order + 1)
    )


def evaluate_coefficient(polynomial, moment, order):
    return QQ(polynomial(MM=moment, LL=order))


def audit_transformed(operator, point):
    for moment in tuple(range(10, 23)) + (31, 37):
        values = [
            ray_packet(moment, order, point)
            for order in range(moment + 1)
        ]
        for order in range(moment - operator.order()):
            total = sum(
                evaluate_coefficient(operator[shift], moment, order)
                * values[order + shift]
                for shift in range(operator.order() + 1)
            )
            assert total == 0, (point, moment, order, total)


transformed_operators = []
trailing_product = RL.one()
leading_product = RL.one()
for point, multiplicity in RAY_CLASSES:
    recurrence, certificate = ray_operator(point)
    transformed = signed_binomial_transform(recurrence)
    audit_transformed(transformed, point)
    transformed_operators.append(transformed)
    trailing = factor(transformed[0])
    leading = factor(transformed[transformed.order()])
    trailing_product *= transformed[0]
    leading_product *= transformed[transformed.order()]
    print(
        "BINOMIAL_RAY",
        point,
        "MULTIPLICITY",
        multiplicity,
        "R_ORDER",
        recurrence.order(),
        "L_ORDER",
        transformed.order(),
    )
    print("  L_TRAILING", trailing)
    print("  L_LEADING", leading)

# The direct sum of the fourteen companion systems is an exact system
# for the all-rays-summed packet.  Its two transition determinants have
# these products.  We keep them factored rather than expanding a large
# scalar LCLM.
print("SUMMED_SYSTEM_DIMENSION", sum(op.order() for op in transformed_operators))
print("SUMMED_SYSTEM_TRAILING_PRODUCT", factor(trailing_product))
print("SUMMED_SYSTEM_LEADING_PRODUCT", factor(leading_product))

# At a target, M=p+r-1 and L=r.  Reducing the endpoint products modulo p
# is the exact singular-factor test requested by the arithmetic argument.
Rp.<p, ell> = QQ[]
target_substitution = {
    MM: p + ell - 1,
    LL: ell,
}
target_trailing = Rp(trailing_product.subs(target_substitution))
target_leading = Rp(leading_product.subs(target_substitution))
print("TARGET_TRAILING_PRODUCT", factor(target_trailing))
print("TARGET_LEADING_PRODUCT", factor(target_leading))
print("Q32_TERMINAL_BINOMIAL_TELESCOPERS=PASS")
