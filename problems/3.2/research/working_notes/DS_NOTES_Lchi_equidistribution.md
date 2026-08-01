# DS note: L(chi) mod p — Claude's reformulation verified + quantified (2026-07-31)

## Verification (p=11, p=101)
- **b_m ≡ −L(χ_m) mod p** (p=11: m=1..9 all; p=101: sample m=1,2,3,5,17,50 all) — the Mellin
  trace connection (matches Q6125's moment formula). L(χ_m)=Σ_a χ_m(a)N_Λ(a).
- **Vertical Weil bound holds**: p=101, |C_p|=|Σ_{a=0}^{p-1}e_p(b_a)| = 8.09 ≈ √p = 10.
- **Phases e_p(L(χ_m)) are equidistributed**: |mean|=0.076 (uniform-circle pred ~0.10),
  variance=0.994 (pred 1.0). This IS the "equidistribution of L(χ) mod p" equivalent to the
  Weil bound.
- **NUANCE**: L(χ_m) mod p takes only 41 distinct values over the 100 characters (max
  multiplicity 4) — NOT uniform over F_p. But the phases still cancel (the distinct values
  are spread + multiplicities balanced). So "equidistribution" = the PHASES e_p(L(χ)) are
  spread on the circle, not L(χ) uniform over F_p.

## Connection to the mod-p-of-Frobenius-trace problem
The vertical Weil bound is a statement about the Mellin moments L(χ) of the point-count
function N_Λ. Its mod-p behavior (41 values, phases spread) is the "mod-p Frobenius trace"
structure. The Dwork/Cartier module (DS_NOTES_dwork_setup) controls exactly this: the
Cartier matrix on periods is the identity (Frobenius-simple), and the L(χ) mod p spread is
what the character sum must explain.

## Files
DS_NOTES_dwork_setup.md (Cartier matrix), CLAUDE_NOTES_two_prime_weil.md (Claude's Weil route),
this note.
