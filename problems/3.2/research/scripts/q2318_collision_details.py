#!/usr/bin/env python3
"""Exact local details for the Q2318 obstruction examples."""
from q2318_associated_resultant import (
    S, apery_half, beta, build_P, eval_poly, gcd_poly, lam,
    resultant, transfer_block,
)


def details(p: int, r: int, s: int) -> dict[str, object]:
    ps = build_P(s, p)
    h = s - r
    c = transfer_block(r, s, p)
    shifted = S(r + 1, h - 2, p)
    lower = S(r, h - 2, p)
    roots = [x for x in range(p) if eval_poly(ps[r], x, p) == 0 and eval_poly(c, x, p) == 0]
    vals = apery_half(p)
    out: dict[str, object] = {
        "p": p, "r": r, "s": s, "h": h,
        "b_r": vals[r], "b_s": vals[s],
        "lambda_r": lam(r, p), "lambda_s": lam(s, p),
        "associated_resultant": resultant(ps[r], c, p),
        "gcd": gcd_poly(ps[r], c, p),
        "common_base_field_roots": roots,
        "beta_r_plus_1": beta(r + 1, p),
    }
    if roots:
        a = roots[0]
        out.update({
            "alpha": a,
            "P_r_alpha": eval_poly(ps[r], a, p),
            "C_alpha": eval_poly(c, a, p),
            "C_shift_alpha": eval_poly(shifted, a, p),
            "C_lower_alpha": eval_poly(lower, a, p),
            "P_r_plus_1_alpha": eval_poly(ps[r + 1], a, p),
            "P_s_minus_1_alpha": eval_poly(ps[s - 1], a, p),
            "P_s_alpha": eval_poly(ps[s], a, p),
        })
    return out


print(details(17, 4, 6))
print(details(139, 31, 61))
