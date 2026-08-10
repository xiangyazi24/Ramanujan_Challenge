#!/usr/bin/env python3
"""Verify the external-second-window condition on every known selected chain.

``primitive_projective_prime_verify.py`` records the complete raw selected
chain census through ``p <= 500000``.  This script takes those fixed records,
recomputes the projective orbit independently, removes each chain's own
reflection support, and asks whether any remaining later occurrence begins
another four-consecutive-occurrence window of span at most ``sqrt(p)``.
"""

from __future__ import annotations

from bisect import bisect_left
import hashlib
import json
import math

from long_bridge_incidence_verify import projective_states, solution_pair
from primitive_projective_prime_verify import RECORDS


EXPECTED_SHA256 = "13f2f869de772d2f30de9d78a1818ad149e9cef4f92adec9ec2eef7c4e25ea1b"


def main() -> None:
    assert len(RECORDS) == 20
    for left, right in zip(RECORDS[::2], RECORDS[1::2]):
        prime, start, gaps = left
        reflected_prime, reflected_start, reflected_gaps = right
        assert reflected_prime == prime
        assert reflected_gaps == tuple(reversed(gaps))
        assert reflected_start == prime - 1 - start - sum(gaps)

    state_cache: dict[int, list[int]] = {}
    rows = []
    for prime, start, gaps in RECORDS:
        if prime not in state_cache:
            apery, companion = solution_pair(prime)
            state_cache[prime] = projective_states(prime, apery, companion)
        states = state_cache[prime]
        state = states[start]
        occurrences = [
            index for index, value in enumerate(states) if value == state
        ]

        first, second, third = gaps
        chain = (
            start,
            start + first,
            start + first + second,
            start + first + second + third,
        )
        position = bisect_left(occurrences, start)
        assert tuple(occurrences[position : position + 4]) == chain
        support = set(chain)
        support.update(prime - 1 - index for index in chain)

        external = [
            index
            for index in occurrences
            if index >= chain[-1] + 2 and index not in support
        ]
        external_short_four = []
        for index in external:
            occurrence_index = bisect_left(occurrences, index)
            assert occurrences[occurrence_index] == index
            if (
                occurrence_index + 3 < len(occurrences)
                and occurrences[occurrence_index + 3] - index
                <= math.isqrt(prime)
            ):
                external_short_four.append(index)

        rows.append(
            {
                "prime": prime,
                "start": start,
                "gaps": gaps,
                "fiber_size": len(occurrences),
                "external": external,
                "external_short_four": external_short_four,
            }
        )

    assert not any(row["external_short_four"] for row in rows)
    encoded = json.dumps(rows, sort_keys=True, separators=(",", ":")).encode()
    digest = hashlib.sha256(encoded).hexdigest()
    assert digest == EXPECTED_SHA256
    print(
        "LONG_BRIDGE_SELECTED_EXTENSION_VERIFY"
        f" records={len(rows)} primes={len(state_cache)}"
        f" external={sum(len(row['external']) for row in rows)}"
        " external_short_four=0"
    )
    print(f"sha256={digest}")
    for row in rows:
        print(json.dumps(row, sort_keys=True))


if __name__ == "__main__":
    main()
