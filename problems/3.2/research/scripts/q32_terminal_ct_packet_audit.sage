#!/usr/bin/env sage
"""Exact terminal constant-term and origin-coordinate audit.

For a first-cell shell and endpoint ``Q`` put

    T(N,Q,L) = sum_{kappa in P(Z)} sum_{r=0}^L
                 (-1)^r binom(L,r)c_N((Q-r)kappa).

This is the coefficient form of

    sum_kappa CT[Lambda^N X^(-Q*kappa)(1-X^kappa)^L].

The script checks the terminal formula, its two-dimensional Pascal law,
and the compact formula for the origin-cancelled W coordinate.
"""

from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

load(str(HERE / "q32_doubled_period_gauge_audit.sage"))
from q32_cartier_packet_audit import coefficient, polytope_points

# The unweighted lattice-point polynomial of P has no hidden product
# decomposition.  Clearing its Laurent monomial gives a primitive
# irreducible polynomial in QQ[x,y,z].
RP.<x, y, z> = PolynomialRing(QQ)
polytope_numerator = sum(
    x ** (kappa[0] + 1)
    * y ** (kappa[1] + 1)
    * z ** (kappa[2] + 1)
    for kappa in polytope_points(1)
)
polytope_factorization = list(polytope_numerator.factor())
assert len(polytope_factorization) == 1
assert polytope_factorization[0][1] == 1


def terminal_packet(moment, endpoint, order):
    """Evaluate T(moment, endpoint, order) coefficient by coefficient."""

    return ZZ(
        sum(
            (-1) ** residue
            * binomial(order, residue)
            * coefficient(
                moment,
                (endpoint - residue) * kappa[0],
                (endpoint - residue) * kappa[1],
                (endpoint - residue) * kappa[2],
            )
            for kappa in polytope_points(1)
            for residue in range(order + 1)
        )
    )


def direct_difference(moment, endpoint, order):
    start = endpoint - order
    values = shell_batch(moment, range(start, endpoint + 1))
    return ZZ(
        sum(
            (-1) ** residue
            * binomial(order, residue)
            * values[endpoint - residue]
            for residue in range(order + 1)
        )
    )


def terminal_origin_coordinate(global_index, endpoint, order):
    """The compact terminal-packet formula for Delta^L W."""

    moment = global_index - 1
    return ZZ(
        sum(
            ZZ(integer_coefficients[shift](n=global_index))
            * (
                terminal_packet(
                    global_index + shift, endpoint - 1, order
                )
                + terminal_packet(
                    global_index + shift, endpoint + 1, order
                )
            )
            for shift in range(origin_operator.order() + 1)
        )
        - ZZ(multiplier(n=global_index))
        * terminal_packet(moment, endpoint, order)
    )


def direct_origin_coordinate(global_index, endpoint, order):
    moment = global_index - 1
    start = endpoint - order
    nodes = list(range(start, endpoint + 1))
    required = list(range(start - 1, endpoint + 2))
    old = shell_batch(moment, nodes)
    shifted = {
        shift: shell_batch(global_index + shift, required)
        for shift in range(origin_operator.order() + 1)
    }
    values = {
        node: sum(
            ZZ(integer_coefficients[shift](n=global_index))
            * (
                shifted[shift][node - 1]
                + shifted[shift][node + 1]
            )
            for shift in range(origin_operator.order() + 1)
        )
        - ZZ(multiplier(n=global_index)) * old[node]
        for node in nodes
    }
    return ZZ(
        sum(
            (-1) ** residue
            * binomial(order, residue)
            * values[endpoint - residue]
            for residue in range(order + 1)
        )
    )


# The minimum node is in the first cell for every shifted moment below.
for global_index, order in ((21, 4), (31, 8)):
    endpoint = global_index - 1
    assert endpoint - order - 1 > (
        global_index + origin_operator.order()
    ) // 2
    assert terminal_packet(endpoint, endpoint, order) == (
        direct_difference(endpoint, endpoint, order)
    )
    assert terminal_origin_coordinate(
        global_index, endpoint, order
    ) == direct_origin_coordinate(global_index, endpoint, order)

    # Phi_{Q,L}=Phi_{Q,L-1}-Phi_{Q-1,L-1}.
    for surplus in range(3):
        shell_moment = endpoint + surplus
        assert terminal_packet(
            shell_moment, endpoint, order
        ) == terminal_packet(
            shell_moment, endpoint, order - 1
        ) - terminal_packet(
            shell_moment, endpoint - 1, order - 1
        )


def ray_value_formula(moment, residue, point):
    """Closed one-fold formula for c_M((M-r)point)."""

    u, v, w = point
    if u == -1:
        total = 0
        for index in range(residue + 1):
            packets = []
            for coordinate in (v, w):
                if coordinate == -1:
                    packets.append(
                        binomial(2 * moment - index, residue - index)
                    )
                elif coordinate == 0:
                    packets.append(
                        binomial(2 * moment - index, moment)
                    )
                else:
                    packets.append(
                        binomial(2 * moment - index, residue)
                    )
            total += (
                binomial(moment, index)
                * binomial(moment, residue - index)
                * packets[0]
                * packets[1]
            )
        return ZZ(total)

    if u == 1:
        total = 0
        for index in range(residue + 1):
            packets = [
                (
                    binomial(
                        moment + residue - index,
                        residue - index,
                    )
                    if coordinate == 0
                    else binomial(
                        moment + residue - index, residue
                    )
                )
                for coordinate in (v, w)
            ]
            total += (
                binomial(moment, index)
                * binomial(moment, residue - index)
                * packets[0]
                * packets[1]
            )
        return ZZ(total)

    assert u == 0
    upper_limit = residue if -1 in (v, w) else moment
    total = 0
    for index in range(upper_limit + 1):
        packets = []
        for coordinate in (v, w):
            if coordinate == -1:
                packets.append(
                    binomial(2 * moment - index, residue - index)
                )
            elif coordinate == 0:
                packets.append(
                    binomial(2 * moment - index, moment)
                )
            else:
                packets.append(
                    binomial(2 * moment - index, residue)
                )
        total += (
            binomial(moment, index) ** 2
            * packets[0]
            * packets[1]
        )
    return ZZ(total)


for moment in range(3, 12):
    for residue in range((moment - 1) // 2 + 1):
        node = moment - residue
        for point in polytope_points(1):
            if point == (0, 0, 0):
                continue
            assert ray_value_formula(
                moment, residue, point
            ) == coefficient(
                moment,
                node * point[0],
                node * point[1],
                node * point[2],
            )

# The (0,0,1) terminal ray begins with the zeta(2) Apery sequence.
for moment in range(12):
    ray_origin = ray_value_formula(moment, 0, (0, 0, 1))
    apery_zeta2 = sum(
        binomial(moment, index) ** 2
        * binomial(moment + index, index)
        for index in range(moment + 1)
    )
    assert ray_origin == apery_zeta2

LQ.<q> = LaurentPolynomialRing(QQ)


def singular_curve_restriction(endpoint, order):
    """Restrict Phi_{endpoint,order} to y=z=-1, x=q."""

    out = LQ.zero()
    for u, v, w in polytope_points(1):
        monomial = (-1) ** (v + w) * q**u
        out += monomial ** (-endpoint) * (1 - monomial) ** order
    return out


def singular_curve_closed_form(endpoint, order):
    """The grouped three-u-level formula for the same restriction."""

    parity = (-1) ** endpoint
    return (
        4 * parity * 2**order
        + q ** (endpoint - order)
        * (
            5 * (q - 1) ** order
            + 4 * parity * (q + 1) ** order
        )
        + 2
        * q ** (-endpoint)
        * (
            (1 - q) ** order
            + parity * (1 + q) ** order
        )
    )


# On the torus curve y=z=-1, Lambda and all its logarithmic
# derivatives vanish.  Hence restriction to this curve descends through
# every toric integration-by-parts quotient.  The formula below is thus a
# certificate for the summed packet, not for an artificial labelled lift.
for endpoint in range(4, 13):
    for order in range(1, endpoint):
        restriction = singular_curve_restriction(endpoint, order)
        assert restriction == singular_curve_closed_form(endpoint, order)
        # The lowest positive exponent is endpoint-order and its
        # coefficient is one of +/-1,+/-9, hence is never zero.
        assert restriction[endpoint - order] == (
            5 * (-1) ** order + 4 * (-1) ** endpoint
        )

print("Q32_TERMINAL_CT_PACKET_AUDIT=PASS")
print("SUMMED_SINGULAR_CURVE_OBSTRUCTION=PASS")

if "--operator-polynomial-factor" in sys.argv:
    RT.<T> = PolynomialRing(R)
    ordinary_symbol = sum(
        RT(integer_coefficients[shift]) * T**shift
        for shift in range(origin_operator.order() + 1)
    )
    assert ordinary_symbol.is_irreducible()
    print(
        "ORIGIN_ORDINARY_SYMBOL",
        "degree_T",
        ordinary_symbol.degree(),
        "degree_n",
        max(coefficient.degree() for coefficient in ordinary_symbol),
        "irreducible",
        True,
    )


def terminal_difference_table(values):
    """Return backward terminal differences B_L from F(M-L)."""

    row = [values[node] for node in sorted(values)]
    out = [row[-1]]
    while len(row) > 1:
        row = [
            row[index + 1] - row[index]
            for index in range(len(row) - 1)
        ]
        out.append(row[-1])
    return out


def joint_modular_data(prime=1000000007):
    """Build (M,L,B_Y,B_W,b_M) rows for recurrence discovery."""

    rows = []
    for moment in range(16, 81):
        global_index = moment + 1
        maximum_order = (moment - 8) // 2
        nodes = list(
            range(moment - maximum_order, moment + 1)
        )
        required = list(
            range(nodes[0] - 1, nodes[-1] + 2)
        )
        old = shell_batch(moment, nodes, modulus=prime)
        shifted = {
            shift: shell_batch(
                global_index + shift,
                required,
                modulus=prime,
            )
            for shift in range(origin_operator.order() + 1)
        }
        mu = ZZ(multiplier(n=global_index)) % prime
        w_values = {
            node: (
                sum(
                    (
                        ZZ(
                            integer_coefficients[shift](
                                n=global_index
                            )
                        )
                        % prime
                    )
                    * (
                        shifted[shift][node - 1]
                        + shifted[shift][node + 1]
                    )
                    for shift in range(
                        origin_operator.order() + 1
                    )
                )
                - mu * old[node]
            )
            % prime
            for node in nodes
        }
        by = terminal_difference_table(old)
        bw = terminal_difference_table(w_values)
        b0 = apery(moment) % prime
        for order in range(1, maximum_order + 1):
            rows.append((moment, order, by, bw, b0))
    return rows


def total_degree_monomials(field, moment, order, degree):
    return [
        field(moment) ** left * field(order) ** right
        for left in range(degree + 1)
        for right in range(degree + 1 - left)
    ]


def guess_joint_relation():
    """Search for a low-order Y/W relation, optionally sourced by b_M."""

    prime = 1000000007
    field = GF(prime)
    rows = joint_modular_data(prime)
    holdout = 48
    for recurrence_order in range(0, 9):
        eligible = [
            row
            for row in rows
            if row[1] + recurrence_order < len(row[2])
        ]
        for degree in range(7):
            monomial_count = (degree + 1) * (degree + 2) // 2
            block_count = 2 * (recurrence_order + 1) + 1
            column_count = block_count * monomial_count
            if len(eligible) <= column_count + holdout:
                continue
            matrix_rows = []
            for moment, order, by, bw, b0 in eligible[:-holdout]:
                monomials = total_degree_monomials(
                    field, moment, order, degree
                )
                matrix_rows.append(
                    [
                        value * monomial
                        for values in (
                            by[
                                order : order
                                + recurrence_order
                                + 1
                            ],
                            bw[
                                order : order
                                + recurrence_order
                                + 1
                            ],
                            (b0,),
                        )
                        for value in values
                        for monomial in monomials
                    ]
                )
            kernel = matrix(field, matrix_rows).right_kernel()
            if kernel.dimension() == 0:
                continue
            for candidate in kernel.basis():
                valid = True
                for moment, order, by, bw, b0 in eligible[-holdout:]:
                    monomials = total_degree_monomials(
                        field, moment, order, degree
                    )
                    vector_values = (
                        by[
                            order : order
                            + recurrence_order
                            + 1
                        ]
                        + bw[
                            order : order
                            + recurrence_order
                            + 1
                        ]
                        + [b0]
                    )
                    total = sum(
                        candidate[
                            block * monomial_count + column
                        ]
                        * value
                        * monomial
                        for block, value in enumerate(vector_values)
                        for column, monomial in enumerate(monomials)
                    )
                    if total:
                        valid = False
                        break
                if valid:
                    print(
                        "JOINT_RELATION_GUESS",
                        recurrence_order,
                        degree,
                        list(candidate),
                    )
                    return
        print("NO_JOINT_RELATION", recurrence_order, 6)
    print("NO_JOINT_RELATION_THROUGH_SEARCH")


if "--guess-joint" in sys.argv:
    guess_joint_relation()


def guess_homogeneous_z_relation():
    """Search a scalar recurrence for Z=W+mu*Y, the only curve-compatible direction."""

    prime = 1000000007
    field = GF(prime)
    raw_rows = joint_modular_data(prime)
    rows = []
    for moment, order, by, bw, b0 in raw_rows:
        mu = field(multiplier(n=moment + 1))
        bz = [
            field(bw[index]) + mu * field(by[index])
            for index in range(len(by))
        ]
        rows.append((moment, order, bz))
    holdout = 64
    for recurrence_order in range(1, 13):
        eligible = [
            row
            for row in rows
            if row[1] + recurrence_order < len(row[2])
        ]
        for degree in range(9):
            monomial_count = (degree + 1) * (degree + 2) // 2
            column_count = (
                recurrence_order + 1
            ) * monomial_count
            if len(eligible) <= column_count + holdout:
                continue
            matrix_rows = []
            for moment, order, bz in eligible[:-holdout]:
                monomials = total_degree_monomials(
                    field, moment, order, degree
                )
                matrix_rows.append(
                    [
                        bz[order + shift] * monomial
                        for shift in range(recurrence_order + 1)
                        for monomial in monomials
                    ]
                )
            kernel = matrix(field, matrix_rows).right_kernel()
            for candidate in kernel.basis():
                valid = True
                for moment, order, bz in eligible[-holdout:]:
                    monomials = total_degree_monomials(
                        field, moment, order, degree
                    )
                    total = sum(
                        candidate[
                            shift * monomial_count + column
                        ]
                        * bz[order + shift]
                        * monomial
                        for shift in range(recurrence_order + 1)
                        for column, monomial in enumerate(monomials)
                    )
                    if total:
                        valid = False
                        break
                if valid:
                    print(
                        "Z_RELATION_GUESS",
                        recurrence_order,
                        degree,
                        list(candidate),
                    )
                    return
        print("NO_Z_RELATION", recurrence_order, 8)
    print("NO_Z_RELATION_THROUGH_SEARCH")


if "--guess-z" in sys.argv:
    guess_homogeneous_z_relation()


def terminal_yw_table(moment, maximum_order):
    """Return the exact terminal high-difference vectors (B_Y,B_W)."""

    global_index = moment + 1
    nodes = list(range(moment - maximum_order, moment + 1))
    required = list(range(nodes[0] - 1, nodes[-1] + 2))
    old = shell_batch(moment, nodes)
    shifted = {
        shift: shell_batch(global_index + shift, required)
        for shift in range(origin_operator.order() + 1)
    }
    mu = ZZ(multiplier(n=global_index))
    w_values = {
        node: sum(
            ZZ(integer_coefficients[shift](n=global_index))
            * (
                shifted[shift][node - 1]
                + shifted[shift][node + 1]
            )
            for shift in range(origin_operator.order() + 1)
        )
        - mu * old[node]
        for node in nodes
    }
    return (
        terminal_difference_table(old),
        terminal_difference_table(w_values),
    )


def audit_terminal_wedges():
    """Measure exact adjacent wedges after primitive vector saturation."""

    for moment in (20, 30, 40, 50, 60, 80):
        maximum_order = (moment - 8) // 2
        by, bw = terminal_yw_table(moment, maximum_order)
        primitive_vectors = []
        for order in range(1, maximum_order + 1):
            content = gcd(by[order], bw[order])
            primitive_vectors.append(
                (
                    by[order] // content,
                    bw[order] // content,
                )
            )
        wedges = [
            left[0] * right[1] - left[1] * right[0]
            for left, right in zip(
                primitive_vectors, primitive_vectors[1:]
            )
        ]
        last = wedges[-1]
        all_wedge_gcd = gcd(wedges)
        component_bits = max(
            abs(entry).nbits()
            for vector in primitive_vectors[-2:]
            for entry in vector
        )
        print(
            "TERMINAL_WEDGE",
            "M",
            moment,
            "L",
            maximum_order - 1,
            "component_bits",
            component_bits,
            "last_wedge_bits",
            abs(last).nbits(),
            "all_wedge_gcd_bits",
            abs(all_wedge_gcd).nbits(),
            "all_wedge_gcd",
            all_wedge_gcd,
        )


if "--wedge-height" in sys.argv:
    audit_terminal_wedges()


def audit_terminal_exterior_sequence():
    """Measure the diagonal once-divided Y/W exterior sequence.

    If ``V_L=(f_L(Y),f_L(W))`` and

        V_L-V_{L-1}=(-1)^L binom(M+1,L) B_L,

    then the adjacent Newton determinant divided by its universal
    Pascal factor is ``det(B_L,V_{L-1})``.  This is the terminal
    specialization of the normalized exterior carrier.
    """

    for moment in (20, 30, 40, 50, 60, 80, 100):
        maximum_order = (moment - 8) // 2
        by, bw = terminal_yw_table(moment, maximum_order)
        fy = [terminal_packet(moment, moment, 0)]
        fw = [
            terminal_origin_coordinate(
                moment + 1, moment, 0
            )
        ]
        for order in range(1, maximum_order + 1):
            pascal = (-1) ** order * binomial(
                moment + 1, order
            )
            fy.append(fy[-1] + pascal * by[order])
            fw.append(fw[-1] + pascal * bw[order])
        exterior = [
            by[order] * fw[order - 1]
            - bw[order] * fy[order - 1]
            for order in range(1, maximum_order + 1)
        ]
        exterior_gcd = gcd(exterior)
        prefix_gcds = []
        running_gcd = ZZ.zero()
        for value in exterior[:10]:
            running_gcd = gcd(running_gcd, value)
            prefix_gcds.append(abs(running_gcd).nbits())
        operator_content = gcd(
            [
                ZZ(multiplier(n=moment + 1)),
                *[
                    ZZ(
                        integer_coefficients[shift](
                            n=moment + 1
                        )
                    )
                    for shift in range(
                        origin_operator.order() + 1
                    )
                ],
            ]
        )
        print(
            "TERMINAL_EXTERIOR",
            "M",
            moment,
            "Lmax",
            maximum_order,
            "gcd_bits",
            abs(exterior_gcd).nbits(),
            "gcd",
            factor(exterior_gcd),
            "operator_content",
            factor(operator_content),
            "quotient",
            factor(exterior_gcd // gcd(
                exterior_gcd, operator_content
            )),
            "first_bits",
            [abs(value).nbits() for value in exterior[:5]],
            "prefix_gcd_bits",
            prefix_gcds,
        )


if "--terminal-exterior" in sys.argv:
    audit_terminal_exterior_sequence()


def audit_hostile_terminal_exterior_prefix():
    """Check the first seven exterior values at hostile target indices."""

    for moment in (199, 271, 299, 320):
        maximum_order = 7
        by, bw = terminal_yw_table(moment, maximum_order)
        fy = [terminal_packet(moment, moment, 0)]
        fw = [
            terminal_origin_coordinate(
                moment + 1, moment, 0
            )
        ]
        exterior = []
        for order in range(1, maximum_order + 1):
            exterior.append(
                by[order] * fw[-1]
                - bw[order] * fy[-1]
            )
            pascal = (-1) ** order * binomial(
                moment + 1, order
            )
            fy.append(fy[-1] + pascal * by[order])
            fw.append(fw[-1] + pascal * bw[order])
        exterior_gcd = gcd(exterior)
        print(
            "HOSTILE_TERMINAL_EXTERIOR_PREFIX",
            "M",
            moment,
            "gcd_bits",
            abs(exterior_gcd).nbits(),
            "gcd_factorization",
            factor(exterior_gcd),
        )


if "--hostile-terminal-exterior" in sys.argv:
    audit_hostile_terminal_exterior_prefix()


def terminal_yw_table_mod(moment, maximum_order, modulus):
    """Return terminal Y/W differences modulo ``modulus``."""

    global_index = moment + 1
    nodes = list(range(moment - maximum_order, moment + 1))
    required = list(range(nodes[0] - 1, nodes[-1] + 2))
    old = shell_batch(moment, nodes, modulus=modulus)
    shifted = {
        shift: shell_batch(
            global_index + shift,
            required,
            modulus=modulus,
        )
        for shift in range(origin_operator.order() + 1)
    }
    mu = ZZ(multiplier(n=global_index)) % modulus
    w_values = {
        node: (
            sum(
                (
                    ZZ(
                        integer_coefficients[shift](
                            n=global_index
                        )
                    )
                    % modulus
                )
                * (
                    shifted[shift][node - 1]
                    + shifted[shift][node + 1]
                )
                for shift in range(
                    origin_operator.order() + 1
                )
            )
            - mu * old[node]
        )
        % modulus
        for node in nodes
    }
    return (
        [
            value % modulus
            for value in terminal_difference_table(old)
        ],
        [
            value % modulus
            for value in terminal_difference_table(w_values)
        ],
        old[moment] % modulus,
        w_values[moment] % modulus,
    )


def audit_target_boundary_jets():
    """Print the first divided terminal exterior digit at exact targets."""

    targets = (
        (199, 139),
        (199, 181),
        (271, 191),
        (271, 233),
        (299, 191),
        (299, 227),
        (320, 179),
        (320, 193),
        (320, 211),
    )
    for moment, prime in targets:
        global_index = moment + 1
        residue = global_index - prime
        modulus = prime**2
        by, bw, fy, fw = terminal_yw_table_mod(
            moment, residue, modulus
        )
        for order in range(1, residue + 1):
            pascal = (
                (-1) ** order
                * binomial(moment + 1, order)
            ) % modulus
            previous_y, previous_w = fy, fw
            fy = (fy + pascal * by[order]) % modulus
            fw = (fw + pascal * bw[order]) % modulus
        exterior = (
            by[residue] * previous_w
            - bw[residue] * previous_y
        ) % modulus
        assert fy % prime == 0
        assert fw % prime == 0
        assert exterior % prime == 0
        print(
            "TARGET_BOUNDARY_JET",
            "M",
            moment,
            "p",
            prime,
            "r",
            residue,
            "Y_div",
            (fy // prime) % prime,
            "W_div",
            (fw // prime) % prime,
            "E_div",
            (exterior // prime) % prime,
            "B",
            (by[residue] % prime, bw[residue] % prime),
        )


if "--target-boundary-jets" in sys.argv:
    audit_target_boundary_jets()


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


def search_polynomial_ray_syzygies():
    """Look for low-degree relations among the 14 summed-label ray classes."""

    prime = 1000000007
    field = GF(prime)
    data = []
    for moment in range(8, 51):
        for residue in range((moment - 1) // 2 + 1):
            node = moment - residue
            values = [
                (
                    multiplicity
                    * coefficient(
                        moment,
                        node * point[0],
                        node * point[1],
                        node * point[2],
                        modulus=prime,
                    )
                )
                % prime
                for point, multiplicity in RAY_CLASSES
            ]
            data.append((moment, residue, values))
    holdout = 80
    for degree in range(7):
        monomial_count = (degree + 1) * (degree + 2) // 2
        column_count = len(RAY_CLASSES) * monomial_count
        if len(data) <= column_count + holdout:
            continue
        rows = []
        for moment, residue, values in data[:-holdout]:
            monomials = total_degree_monomials(
                field, moment, residue, degree
            )
            rows.append(
                [
                    value * monomial
                    for value in values
                    for monomial in monomials
                ]
            )
        kernel = matrix(field, rows).right_kernel()
        valid_dimension = 0
        for candidate in kernel.basis():
            if all(
                sum(
                    candidate[
                        block * monomial_count + column
                    ]
                    * value
                    * monomial
                    for block, value in enumerate(values)
                    for column, monomial in enumerate(
                        total_degree_monomials(
                            field, moment, residue, degree
                        )
                    )
                )
                == 0
                for moment, residue, values in data[-holdout:]
            ):
                valid_dimension += 1
        print(
            "RAY_SYZYGY_SEARCH",
            "degree",
            degree,
            "columns",
            column_count,
            "training_kernel",
            kernel.dimension(),
            "holdout_valid",
            valid_dimension,
        )


if "--ray-syzygy" in sys.argv:
    search_polynomial_ray_syzygies()
