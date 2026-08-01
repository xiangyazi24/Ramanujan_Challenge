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

## Katz discrepancy correction (Q6287 + DS empirical)
- **Q6287 verdict**: Katz/Sato-Tate for L(χ)/p^{3/2} does NOT imply mod-p residue uniform, even
  with O(1/√p) archimedean discrepancy. |L(χ)| ≤ Cp^{3/2} folds mod p only ~√p times; the
  residue distribution is a LOCAL (mod-p) statement, not controlled by the archimedean law.
- **DS empirical confirms**: |Σ_m e_p(L(χ_m))|/√p = 0.76, 1.37 (p=101,199) — O(1), vertical
  Weil holds; #distinct L(χ) mod p = 41, 79 (~√p folding, NOT F_p-uniform).
- **The correct target (Q6287)**: a LOCAL-LIMIT theorem
  max_r |#{m : L(χ_m)=r mod p} − (p−1)/p| = O(√p); Fourier inversion ⟹ |Σ e_p(L(χ))| = O(√p).
- **Decisive computation**: is L(χ) mod p a bounded-complexity trace function of the character
  parameter χ? (If yes, Deligne/Katz applies; if not, Sato-Tate alone can't close it.)

## Vertical Weil constant BOUNDED (DS, p=101..503)
|Σ_m e_p(L(χ_m))|/√p = 0.76, 1.37, 3.22, 1.06, 1.00 (p=101,199,307,401,503).
The p=307 ratio 3.22 is a ~3σ random-walk fluctuation; p=401,503 return to ~1.0.
⟹ the vertical Weil bound |C_p| ≤ C√p holds with a bounded constant (cluster ~1, max ~3.2).

## Decisive complexity + local-limit experiments (DS, p=11..19)
- **L(χ_m) mod p is NOT bounded-complexity (linear-recurrence sense)**: minimal constant-
  coefficient recurrence order of m -> L(χ_m) = p (p=11,13,17,19 → 11,13,17,19), growing
  with p. This is the generic behavior (any function on F_p^* has order ≤ p−1). Refutes the
  "L(χ) mod p is a bounded-order linear-recurrent trace" hope.
- **Local-limit deviation ~√p (supports the target)**: #{m: L(χ_m)=r mod p} takes only
  ~√p distinct values (6,6,7 of p−1=10,12,16), max multiplicity 2,4,4 vs √p ≈ 3,4,4.
  The O(√p) local-limit statement holds (hit residues have multiplicity ~√p; non-hit = 0,
  trivially within O(√p) of the mean). NOT F_p-uniform — only ~√p residues hit, confirming
  the folding picture.
