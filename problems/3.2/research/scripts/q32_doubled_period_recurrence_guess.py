#!/usr/bin/env python3
"""Reproduce the exact data behind the doubled-period recurrence guess.

Define

    h_s = CT(Lambda(X)^s Lambda(X^2)),
    J_s = h_s + 40*b_s.

An exact rational-nullspace calculation on the first 103 recurrence
positions found the order-three, degree-21 candidate recorded below.
The first 100 positions were used to reconstruct the one-dimensional
nullspace; positions 100--122 were held out.  This dependency-free script
does not redo the rational nullspace calculation.  It independently
generates the actual coefficients from formula (49.1) and verifies the
candidate on both the training and held-out positions.

This script alone is a finite exact audit, not a telescoping certificate.
The companion ``q32_doubled_period_telescoper.sage`` now derives an exact
Zeilberger certificate and proves that its operator is precisely ``-4``
times the primitive operator stored here.
"""

from __future__ import annotations

from math import gcd

from q32_cartier_packet_audit import LAMBDA, coefficient


RECURRENCE = (
    (
        9585039744000,
        68753480659200,
        210146974558560,
        334475579541528,
        222682732634392,
        -172898826462224,
        -554675454152281,
        -585156141287484,
        -276537292394282,
        74111519480120,
        238959524223092,
        219098839094610,
        130005921803994,
        56606796649202,
        18847904358725,
        4860889665898,
        968247842904,
        146595361934,
        16340444704,
        1265923648,
        60936192,
        1373568,
    ),
    (
        -9585039744000,
        -27218308435200,
        -260648986415040,
        -1378737502607520,
        -4654090268946900,
        -12300680777950668,
        -25619361686795425,
        -40747787754762200,
        -49243180538883564,
        -45766980047585484,
        -33231105984106752,
        -19104872754366546,
        -8777026645183450,
        -3236610761159462,
        -957449588327967,
        -225810562663662,
        -41917224674694,
        -5992228655946,
        -636667895776,
        -47350565312,
        -2200071552,
        -48074880,
    ),
    (
        10504629504000,
        37252155571200,
        111314783295360,
        535298131943568,
        2053939023458892,
        5941318973833770,
        13124660500594643,
        21780535169045784,
        27195869550502460,
        25989849537062710,
        19404481058967302,
        11519261125173136,
        5502206777359876,
        2125781298488182,
        663242737472047,
        165705503144168,
        32638509241804,
        4945939617098,
        555365357088,
        43458515264,
        2113536768,
        48074880,
    ),
    (
        -476824320000,
        -1637953920000,
        -4287069511200,
        -20923661696520,
        -81479915366144,
        -236142522921518,
        -519324884985025,
        -852051982005388,
        -1046246697712302,
        -980483217652986,
        -717194046914066,
        -417131089829528,
        -195254706803116,
        -73925664253738,
        -22594128598005,
        -5526099615540,
        -1064671483358,
        -157662660110,
        -17280084064,
        -1317744064,
        -62309760,
        -1373568,
    ),
)


def apery_values(limit: int) -> list[int]:
    values = [1, 5]
    for n in range(1, limit):
        numerator = (
            (34 * n**3 + 51 * n**2 + 27 * n + 5) * values[n]
            - n**3 * values[n - 1]
        )
        value, remainder = divmod(numerator, (n + 1) ** 3)
        assert remainder == 0
        values.append(value)
    return values[: limit + 1]


def doubled_period(index: int) -> int:
    return sum(
        weight * coefficient(index, -2 * u, -2 * v, -2 * w)
        for (u, v, w), weight in LAMBDA.items()
    )


def polynomial_value(coefficients: tuple[int, ...], index: int) -> int:
    out = 0
    for coefficient_value in reversed(coefficients):
        out = out * index + coefficient_value
    return out


def main() -> None:
    final_index = 125
    apery = apery_values(final_index)
    doubled = [doubled_period(index) for index in range(final_index + 1)]
    values = [
        doubled[index] + 40 * apery[index]
        for index in range(final_index + 1)
    ]

    assert doubled[:8] == [
        5,
        25,
        545,
        14917,
        429029,
        12570545,
        372777785,
        11164475165,
    ]

    recurrence_checks = 0
    for index in range(len(values) - 3):
        total = sum(
            polynomial_value(RECURRENCE[shift], index)
            * values[index + shift]
            for shift in range(4)
        )
        assert total == 0, index
        recurrence_checks += 1

    # These are the first two large-prime common factors of
    # gcd(b_{s+1}, J_s).  They are exactly the first two exceptional
    # shortest-jet targets in q32_newton_second_layer_audit.py.
    exceptional_checks = []
    for index, prime in ((55, 61), (76, 139)):
        common = gcd(apery[index + 1], values[index])
        assert common % prime == 0
        exceptional_checks.append((index, prime, common))

    print("Q32_DOUBLED_PERIOD_RECURRENCE_GUESS=PASS")
    print("EXACT_TERMS", len(values))
    print("TRAINING_POSITIONS", 100)
    print("HELD_OUT_POSITIONS", recurrence_checks - 100)
    print("ORDER", 3)
    print("DEGREE", 21)
    print("EXCEPTIONAL_COMMON_FACTORS", exceptional_checks)
    print("STATUS", "CERTIFIED_BY_Q32_DOUBLED_PERIOD_TELESCOPER")


if __name__ == "__main__":
    main()
