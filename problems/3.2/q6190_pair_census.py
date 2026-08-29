#!/usr/bin/env python3
"""Q6190 scratch-only structural census of the saved exact Apéry zero pairs.

This does not change theorem sources.  It reads the repository's canonical
headerless little-endian (p,r) file, checks its SHA256, and reports only exact
finite combinatorics of the stored actual zero predicates b_r == 0 (mod p).
"""
from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import struct

HERE = Path(__file__).resolve().parent
DATA = HERE / "data_zp_pairs.bin"
EXPECTED_SHA256 = "8746d0b400c1b669b001eae955c602908a10c9ee4cb3cac62c6676ea2ddd874d"


def main() -> None:
    raw = DATA.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    assert digest == EXPECTED_SHA256
    assert len(raw) % 8 == 0

    grouped: dict[int, list[int]] = defaultdict(list)
    previous = (-1, -1)
    for pair in struct.iter_unpack("<II", raw):
        p, r = pair
        assert pair > previous
        assert 0 <= r < p
        grouped[p].append(r)
        previous = pair
    zero_sets = {p: tuple(rs) for p, rs in grouped.items()}

    zero_count_dist = Counter(map(len, zero_sets.values()))
    gap_lengths: Counter[int] = Counter()
    gap2_words: Counter[tuple[int, int]] = Counter()
    gap3_words: Counter[tuple[int, int, int]] = Counter()
    triples = 0
    triple_adj_centered = 0
    triple_outer_reflection_pair = 0
    triple_reflection_fixed = 0
    triple_orbit_representatives = 0
    quadruples = 0
    quad_centered_adjacent = 0
    quad_off_center = 0
    quad_short_off_center = 0
    closest_quad = None
    closest_ratio = None
    midpoint_primes = []
    symmetry_failures = []
    consecutive_failures = []

    for p, zeros in zero_sets.items():
        assert tuple(sorted(set(zeros))) == zeros
        zset = set(zeros)
        mid = (p - 1) // 2
        if mid in zset:
            midpoint_primes.append(p)
        for r in zeros:
            if p - 1 - r not in zset:
                symmetry_failures.append((p, r))
            if r + 1 in zset:
                consecutive_failures.append((p, r))

        for x, y in zip(zeros, zeros[1:]):
            gap_lengths[y - x] += 1
        for i in range(len(zeros) - 2):
            x, y, z = zeros[i:i+3]
            gaps = (y - x, z - y)
            gap2_words[gaps] += 1
            triples += 1
            if x + y == p - 1 or y + z == p - 1:
                triple_adj_centered += 1
            if x + z == p - 1:
                triple_outer_reflection_pair += 1
            if x + z == p - 1 and 2 * y == p - 1:
                triple_reflection_fixed += 1
            mirror = (p - 1 - z, p - 1 - y, p - 1 - x)
            if (x, y, z) <= mirror:
                triple_orbit_representatives += 1

        for i in range(len(zeros) - 3):
            chain = zeros[i:i+4]
            gaps = tuple(chain[j+1] - chain[j] for j in range(3))
            gap3_words[gaps] += 1
            quadruples += 1
            centered = [j for j in range(3) if chain[j] + chain[j+1] == p - 1]
            if centered:
                quad_centered_adjacent += 1
            else:
                quad_off_center += 1
                span = chain[-1] - chain[0]
                ratio = Fraction(span * span, p)
                if closest_ratio is None or ratio < closest_ratio:
                    closest_ratio = ratio
                    closest_quad = (p, chain, gaps)
                if span * span <= p:
                    quad_short_off_center += 1

    assert not symmetry_failures
    assert not consecutive_failures
    assert quad_short_off_center == 0

    max_item = max(zero_sets.items(), key=lambda kv: (len(kv[1]), -kv[0]))
    payload = {
        "data_sha256": digest,
        "records": len(raw) // 8,
        "active_primes": len(zero_sets),
        "minimum_active_prime": min(zero_sets),
        "maximum_active_prime": max(zero_sets),
        "zero_count_distribution": sorted(zero_count_dist.items()),
        "maximum_zero_set": [max_item[0], list(max_item[1])],
        "midpoint_zero_primes": midpoint_primes,
        "gap_occurrences": sum(gap_lengths.values()),
        "distinct_gap_lengths": len(gap_lengths),
        "top_gap_lengths": gap_lengths.most_common(20),
        "consecutive_triples": triples,
        "distinct_gap2_words": len(gap2_words),
        "top_gap2_words": [[list(k), v] for k, v in gap2_words.most_common(20)],
        "triple_adjacent_reflection_pair": triple_adj_centered,
        "triple_outer_reflection_pair": triple_outer_reflection_pair,
        "triple_reflection_fixed": triple_reflection_fixed,
        "triple_reflection_orbit_representatives": triple_orbit_representatives,
        "consecutive_quadruples": quadruples,
        "distinct_gap3_words": len(gap3_words),
        "top_gap3_words": [[list(k), v] for k, v in gap3_words.most_common(20)],
        "quad_centered_adjacent": quad_centered_adjacent,
        "quad_off_center": quad_off_center,
        "quad_short_off_center": quad_short_off_center,
        "closest_off_center_quad": None if closest_quad is None else [closest_quad[0], list(closest_quad[1]), list(closest_quad[2])],
        "closest_span_squared_over_p": None if closest_ratio is None else [closest_ratio.numerator, closest_ratio.denominator],
        "symmetry_failures": len(symmetry_failures),
        "consecutive_zero_failures": len(consecutive_failures),
    }
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    payload_sha = hashlib.sha256(canonical).hexdigest()

    print("Q6190_PAIR_CENSUS")
    print(json.dumps(payload, sort_keys=True, indent=2))
    print(f"payload_sha256={payload_sha}")


if __name__ == "__main__":
    main()
