# DS note: L(chi) mod p — Claude's reformulation verified + quantified (2026-07-31)

**[CORRECTION BLOCK 2026-07-31 Fable — supersedes the folding story below;
full data in FABLE_NOTES_vertical_value_law.md, script fable_vertical_value_law.py]**
1. The "41 distinct values ≈ √p folding" reading is WRONG. The law is LINEAR:
   D(p) = (1−e^{−1/2})·p = 0.3935p — Poisson occupancy conditioned on the
   reflection FE b_r ≡ b_{p−1−r} (Q3.2_density_theorem.md §5.9). Measured over
   164 primes ≤ 1000: D/p = 0.3963 ± 0.018; p=997: 392 observed vs 391.8 predicted.
   41 = 0.406·101, 79 = 0.397·199. No archimedean p^{3/2}-folding content at all.
2. The value multiset {b_m mod p} is statistically "iid uniform mod FE" four
   constants deep: E(p) fluct std 2.67 vs pred 2√2; mean|C_p(1)|/√p = 1.261 vs
   Rayleigh 1.2533 (this explains the measured 1.26); max_h = Gumbel
   √(2p ln(p/2)) (so uniform-h C√p is false; fixed-h / E(p)=O(p) is the target).
3. The "decisive computation" (bounded-complexity trace function?) is answered
   NO by proof: e_p(L(χ_m)) is the trace function of the rank-1 Artin–Schreier
   sheaf of f(X) = Σ_s N_Λ(g^s)X^s with Swan_∞ = deg f ≈ p (the measured
   recurrence order = p IS this degree). Katz/Deligne in the χ-aspect returns
   exactly the archimedean p^{3/2} — same bound, two languages.

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

## ⚠️ CORRECTION (Fable, 8b7bbf5) — my "√p folding" story was WRONG
Fable (strongest model) refuted the folding mechanism: D(p) = #distinct L(χ) mod p = 
0.3963p ± 0.018 (measured 164 primes ≤ 1000), EXACTLY p(1−e^{−1/2}) = 0.3935p — the Poisson
occupancy law under the reflection FE (b_r ≡ b_{p−1−r}), NOT archimedean folding. My p=11..19
was too small to discriminate (the two laws coincide at p~11). Verified: D/p = 41/101=0.406,
79/199=0.397, 121/307=0.394 — confirms 0.3963 (41 ≫ √101=10, so NOT √p). The "|L|≤Cp^{3/2}
folds only √p times" story is REMOVED. Distinct-value deficit (0.39 vs Poisson 0.63) = 100%
reflection FE, 0% folding.
Also per Fable: mean|C_p(1)|/√p = Rayleigh 1.2533 (my measured 1.26 matches); max_h|C_p(h)|
= Gumbel √(2p·ln(p/2)) (3.27 explained); the uniform-in-h C√p statement is FALSE (fixed-h /
E(p)=O(p) is the correct target); e_p(hL(χ_m)) IS a rank-1 Artin-Schreier trace with
Swan_∞=deg f ~ p (my recurrence-order-p confirmed), so (iv) fails at maximal strength and
GOS gives p^{3/2} — the ℓ-adic and archimedean bounds are the same object in two languages.
**E(p) = O(p) is the single cleanest open vertical statement** (Parseval: Σ_{h≠0}|C_p(h)|²
= p·E(p) − p²).
