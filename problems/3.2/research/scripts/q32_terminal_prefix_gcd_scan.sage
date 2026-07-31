#!/usr/bin/env sage
"""Exact factor scan for the fixed terminal exterior prefix.

This is diagnostic only.  It reuses the certified Y/W terminal packet
definitions and records

    gcd_{1 <= L <= K} det(B_L, V_{L-1})

for a fixed small K.  The scan is intended to reveal the polynomial
boundary factors that a later symbolic Bezout certificate must explain.
"""

from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
load(str(HERE / "q32_terminal_ct_packet_audit.sage"))


def exterior_prefix_gcd(moment, prefix_length=7):
    by, bw = terminal_yw_table(moment, prefix_length)
    fy = [terminal_packet(moment, moment, 0)]
    fw = [terminal_origin_coordinate(moment + 1, moment, 0)]
    exterior = []
    for order in range(1, prefix_length + 1):
        exterior.append(
            by[order] * fw[order - 1]
            - bw[order] * fy[order - 1]
        )
        pascal = (-1) ** order * binomial(moment + 1, order)
        fy.append(fy[-1] + pascal * by[order])
        fw.append(fw[-1] + pascal * bw[order])
    raw = []
    for residue in range(prefix_length + 1):
        raw.append(
            (
                sum(
                    (-1) ** order
                    * binomial(residue, order)
                    * by[order]
                    for order in range(residue + 1)
                ),
                sum(
                    (-1) ** order
                    * binomial(residue, order)
                    * bw[order]
                    for order in range(residue + 1)
                ),
            )
        )
    raw_minors = {
        (left, right):
        raw[left][0] * raw[right][1]
        - raw[left][1] * raw[right][0]
        for left in range(prefix_length + 1)
        for right in range(left + 1, prefix_length + 1)
    }
    return gcd(exterior), exterior, gcd(raw_minors.values()), raw_minors


moments = (
    23, 29, 37, 43, 53, 61, 71, 83, 97, 109,
    127, 149, 173, 199, 223, 251, 271, 299, 320,
)
if "--short" in sys.argv:
    moments = (199, 320)
for moment in moments:
    value, exterior, raw_gcd, raw_minors = exterior_prefix_gcd(moment)
    assert value % raw_gcd == 0
    smallest_pairs = sorted(
        raw_minors,
        key=lambda pair: abs(raw_minors[pair]).nbits(),
    )[:3]
    print(
        "PREFIX_GCD",
        "M", moment,
        "N", moment + 1,
        "BITS", abs(value).nbits(),
        "FACTOR", factor(value),
        "RAW_GCD_FACTOR", factor(raw_gcd),
        "RATIO_FACTOR", factor(value // raw_gcd),
        "SMALLEST_RAW_PAIRS", tuple(
            (pair, abs(raw_minors[pair]).nbits())
            for pair in smallest_pairs
        ),
        "EXTERIOR_BITS", tuple(abs(entry).nbits() for entry in exterior),
    )
    if moment in (199, 320):
        running = ZZ.zero()
        ledger = []
        for residue in range(7):
            minor = raw_minors[(residue, residue + 1)]
            running = gcd(running, minor)
            ledger.append(
                (
                    residue,
                    abs(minor).nbits(),
                    factor(running),
                )
            )
        print("ADJACENT_RAW_GCD_LEDGER", "M", moment, ledger)

print("Q32_TERMINAL_PREFIX_GCD_SCAN=PASS")
