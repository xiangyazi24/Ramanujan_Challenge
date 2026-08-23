#!/usr/bin/env python3
"""Verify and print the minimal Q2273 certificate after G(T)=T+13 is known."""

from pathlib import Path
import json

import q2273_certificate as q


def block(obj: object) -> None:
    print("```json")
    print(json.dumps(obj, separators=(",", ":")))
    print("```")


def main() -> None:
    _, _, a_poly, kernel = q.build_objects()
    node = [0] * 42
    node[1] = -1
    node[41] = 1
    node = q.utrim(node)
    g = q.ugcd(a_poly, node)
    d, s, t = q.uxgcd(a_poly, node)
    assert d == g == [13, 1]

    h = q.bmul(
        {(0, 0): 1, (1, 0): 4},
        {(0, 0): 1, (0, 1): 4},
    )
    m = q.bmul({(1, 0): 1, (0, 1): -1}, h)

    # Minimal saturated cofactors.
    bco = {}
    cco = q.bmul(h, q.bembed_t(s))
    dco = q.bscale(q.bmul(h, q.bembed_u(s)), -1)
    eco = q.bmul(h, q.bembed_t(t))
    fco = q.bscale(q.bmul(h, q.bembed_u(t)), -1)

    rhs = q.badd(q.bmul(bco, kernel), q.bmul(cco, q.bembed_t(a_poly)))
    rhs = q.badd(rhs, q.bmul(dco, q.bembed_u(a_poly)))
    rhs = q.badd(rhs, q.bmul(eco, q.bembed_t(node)))
    rhs = q.badd(rhs, q.bmul(fco, q.bembed_u(node)))
    assert q.bnorm(rhs) == q.bnorm(m)

    # Strong unsaturated certificate.
    strong_rhs = q.badd(
        q.bmul(q.bembed_t(s), q.bembed_t(a_poly)),
        q.bscale(q.bmul(q.bembed_u(s), q.bembed_u(a_poly)), -1),
    )
    strong_rhs = q.badd(strong_rhs, q.bmul(q.bembed_t(t), q.bembed_t(node)))
    strong_rhs = q.badd(
        strong_rhs,
        q.bscale(q.bmul(q.bembed_u(t), q.bembed_u(node)), -1),
    )
    assert strong_rhs == q.bnorm({(1, 0): 1, (0, 1): -1})

    print()
    print("# Minimal certificate (the actual p=41 outcome)")
    print()
    print("The root computation gives `G(T)=T+13`, so the full two-node problem")
    print("collapses before the Racah kernel is needed.  The extended Euclid output")
    print("is the exact identity")
    print()
    print("```text")
    print("T+13 = s(T) A_41(T) + t(T) (T^41-T).")
    print("```")
    print()
    print("Subtract its U-copy.  This immediately gives the stronger certificate")
    print()
    print("```text")
    print("T-U = 0*K_41(T,U)")
    print("    + s(T)A_41(T) - s(U)A_41(U)")
    print("    + t(T)(T^41-T) - t(U)(U^41-U).")
    print("```")
    print()
    print("Hence `N=1` and, for `H(T,U)=(1+4T)(1+4U)`, a minimal saturated")
    print("choice is")
    print()
    print("```text")
    print("B(T,U) = 0,")
    print("C(T,U) = H(T,U) s(T),")
    print("D(T,U) = -H(T,U) s(U),")
    print("E(T,U) = H(T,U) t(T),")
    print("F(T,U) = -H(T,U) t(U).")
    print("```")
    print()
    print("Thus")
    print()
    print("```text")
    print("(T-U)(1+4T)(1+4U)")
    print(" = B*K_41 + C*A_41(T) + D*A_41(U)")
    print(" + E*(T^41-T) + F*(U^41-U)")
    print("```")
    print()
    print("with coefficient lists (low to high)")
    print("`s(T) =`")
    block(s)
    print("`t(T) =`")
    block(t)
    print()
    print("Machine checks:")
    block({
        "G": g,
        "minimal_saturated_identity": True,
        "strong_T_minus_U_identity": True,
        "B_is_zero": not bco,
        "C_sha256": q.digest_b(cco),
        "D_sha256": q.digest_b(dco),
        "E_sha256": q.digest_b(eco),
        "F_sha256": q.digest_b(fco),
    })
    print()
    print("## Minimal-certificate checker")
    print()
    print("```python")
    print(Path(__file__).read_text(encoding="utf-8").rstrip())
    print("```")


if __name__ == "__main__":
    main()
