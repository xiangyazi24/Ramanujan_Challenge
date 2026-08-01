#!/usr/bin/env python3
"""Gate-only audit for the proposed joint Tannakian moment computation.

The finite-field Mellin inversion gate is exact and passes.  The moment
computation is intentionally not run: the eigen-trace formula in
CODEX_SPEC_joint_tannaka.md conflicts with the deck decomposition already
verified by CRON_pushforward_check.py.  See CODEX_JOINT_TANNAKA.md.
"""

from __future__ import annotations

from pathlib import Path
import runpy


ROOT = Path(__file__).resolve().parents[2]
SOURCE = runpy.run_path(str(ROOT / "CRON_pushforward_check.py"))
PRIMES = (29, 37, 41, 53, 61, 73, 89, 101)


def check_prime(prime: int) -> tuple[int, int, int, int, int]:
    """Run the exact residual gate and exhibit a split-fibre contradiction."""
    data, apery = SOURCE["check_tasks_one_and_two"](prime)
    legendre = SOURCE["legendre"]
    branch_polynomial = SOURCE["branch_polynomial"]
    centered_lift = SOURCE["centered_lift"]

    fibres = data["fibres"]
    traces = data["traces"]
    tensor_pushforward = data["pushforward"]
    apery_values = data["normalized_apery"]
    assert isinstance(fibres, list)
    assert isinstance(traces, list)
    assert isinstance(tensor_pushforward, list)
    assert isinstance(apery_values, list)

    # Exact weight-two convention: remove the determinant/Tate trace p once
    # for every rational source fibre.
    sym2_pushforward = [
        tensor_pushforward[t] - prime * len(fibres[t])
        for t in range(prime)
    ]
    chi_q = [
        legendre(branch_polynomial(t, prime), prime)
        for t in range(prime)
    ]

    # This is the exact pointwise congruence established by the source script.
    assert all(
        sym2_pushforward[t] % prime
        == (1 + chi_q[t]) * apery_values[t] % prime
        for t in range(prime)
    )

    # Formal residues obtained from the two equations printed in the joint
    # spec.  They exist modulo p, but the gate below cannot lift them to
    # characteristic-zero rank-three Frobenius traces.
    inverse_two = pow(2, -1, prime)
    specified_difference = [
        chi_q[t] * apery_values[t] % prime for t in range(prime)
    ]
    plus_residue = [
        (sym2_pushforward[t] + specified_difference[t]) * inverse_two % prime
        for t in range(prime)
    ]
    minus_residue = [
        (sym2_pushforward[t] - specified_difference[t]) * inverse_two % prime
        for t in range(prime)
    ]
    assert all(
        (plus_residue[t] + minus_residue[t]) % prime
        == sym2_pushforward[t] % prime
        and (plus_residue[t] - minus_residue[t]) % prime
        == specified_difference[t]
        for t in range(prime)
    )

    # The requested Mellin gate is the already verified virtual identity
    # -P + chi(q) A = -A modulo p.  Check every exponent, including the two
    # endpoint aliases exactly as in CRON_pushforward_check.py.
    virtual = [
        (-sym2_pushforward[t] + specified_difference[t]) % prime
        for t in range(prime)
    ]
    assert all(virtual[t] == -apery_values[t] % prime for t in range(prime))
    raw_failures: list[int] = []
    for exponent in range(prime):
        raw_mellin = sum(
            virtual[t]
            * pow(t, (-exponent) % (prime - 1), prime)
            for t in range(1, prime)
        ) % prime
        if raw_mellin != apery[exponent]:
            raw_failures.append(exponent)
        endpoint_alias = 0
        if exponent == 0:
            endpoint_alias = apery[prime - 1]
        elif exponent == prime - 1:
            endpoint_alias = apery[0]
        assert (raw_mellin - endpoint_alias) % prime == apery[exponent]
    assert raw_failures == [0, prime - 1]
    print(
        f"GATE VERIFIED p={prime}: all r=0..{prime - 1}; "
        f"raw endpoint aliases={raw_failures}"
    )

    # At a split fibre, the two exact source Sym^2 traces agree (the explicit
    # isogeny relation checked elsewhere in the repository explains this).
    # Hence the two deck descents have equal trace there and difference zero.
    # The spec instead assigns the generally nonzero residue chi(q) A_p.
    witness = None
    for t in range(1, prime):
        if chi_q[t] != 1 or apery_values[t] == 0 or len(fibres[t]) != 2:
            continue
        source_values = [traces[u] ** 2 - prime for u in fibres[t]]
        if source_values[0] == source_values[1]:
            witness = (t, apery_values[t], source_values[0])
            break
    assert witness is not None
    witness_t, witness_apery, source_trace = witness
    assert source_trace % prime == witness_apery
    assert specified_difference[witness_t] == witness_apery
    assert witness_apery != 0

    # A centered residue lift is the only canonical integer lift present in
    # the supplied arithmetic, and it often makes the proposed halves
    # nonintegral.  Other lifts differ by p and are not selected by the
    # residual Mellin gate.
    parity_failures = sum(
        (
            sym2_pushforward[t]
            + chi_q[t] * centered_lift(apery_values[t], prime)
        )
        % 2
        != 0
        for t in range(1, prime)
    )
    assert parity_failures > 0
    print(
        f"EIGEN-TRACE RECONSTRUCTION REFUTED p={prime}: split witness "
        f"t={witness_t}, source Sym2 traces=({source_trace},{source_trace}), "
        f"true deck difference=0, specified chi(q)A residue={witness_apery}; "
        f"centered-lift parity failures={parity_failures}"
    )
    return (
        prime,
        witness_t,
        witness_apery,
        source_trace,
        parity_failures,
    )


def main() -> None:
    results = [check_prime(prime) for prime in PRIMES]
    assert len(results) == len(PRIMES)
    print(
        "NORMALIZATION REFUTED: a weight-2 local trace has Mellin weight 3; "
        "the Weil normalization is p^(3/2), not p"
    )
    print(
        "STALL: no moments computed because neither an exact integral lift "
        "nor the stated deck-eigen decomposition is valid"
    )


if __name__ == "__main__":
    main()
