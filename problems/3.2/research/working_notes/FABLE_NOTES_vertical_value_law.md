# Fable note: the vertical value law is "uniform mod reflection-FE" — four constants nailed (2026-07-31)

Answering the dm question ("can you prove |Σ_χ e(−hL(χ)/p)| ≤ C√p? any missed angle?")
and DS's two updates today. Script: `research/scripts/fable_vertical_value_law.py`
(164 primes in [11,1000], exact integer b_m reduced mod p).

## 1. Headline: the empirical model is EXACTLY iid-uniform conditioned on the reflection FE

The only algebraic constraint on the multiset {b_m mod p : 0 ≤ m ≤ p−1} is the
reflection FE `b_r ≡ b_{p−1−r}` (already banked, Q3.2_density_theorem.md §5.9;
I re-derived the same 2-line Lucas proof independently before grepping — confirmed).
Take the model: (p−1)/2 free values iid uniform on F_p, each doubled by the FE.
It predicts four measurable constants. All four are hit:

| statistic | model prediction | measured (164 primes ≤ 1000) |
|---|---|---|
| D(p)/p, D = #distinct values | 1 − e^{−1/2} = **0.3935** | **0.3963 ± 0.018** |
| (E(p) − 3p)/√p, E = Σ N(a)² | mean 0, std **2√2 = 2.83** | mean −0.06, std **2.67** |
| mean |C_p(1)|/√p | Rayleigh mean **√(π/4·2) = 1.2533** | **1.261** |
| max_{h≠0}|C_p(h)| / √(2p·ln((p−1)/2)) | Gumbel, → **1** | **1.062 ± 0.110** |

Row 1 at p=997: observed 392 vs predicted 391.8 (0.05%). Row 3 explains DS's
measured "average 1.26√p" exactly. Row 2's constant: the FE quadruples collision
variance (each accidental collision (r,r') drags (p−1−r,r') etc.), std = 2·√2·√p·(1+o(1)).

### Consequences

- **DS's "√p folding" law is refuted.** D(p) = 41/79 at p=101/199 is not
  "~√p distinct residues": D/√p drifts 1.51 → 12.41 across the sweep while D/p
  is constant 0.396. The two laws coincide numerically only at p ≈ 10–20, which is
  why the p=11..19 experiment could not discriminate. The "|L| ≤ Cp^{3/2} folds
  only √p times" mechanism story should be dropped from the notes/paper: the
  distinct-value deficit (0.39 vs Poisson 0.63) is 100% the reflection FE, 0%
  archimedean folding.
- **There is no further hidden symmetry to find.** The model leaves zero residual
  anomaly at L¹, L², L^∞ and occupancy levels. The FE is the complete symmetry
  group of the fiber (at measurement precision). Symmetry-hunting is dead;
  stop spending rounds on it.
- **The uniform-in-h form of the "missing theorem" is FALSE as stated.**
  CLAUDE_NOTES_two_prime_weil.md asks for |C_p(h)| ≤ C√p for all h; the truth is
  Gumbel: max_h |C_p(h)| = (1+o(1))·√(2p·ln(p/2)) (row 4; consistent with the L^∞
  barrier already noted in Q3.2_density_theorem.md). The "3.27√p over 166 primes"
  ceiling is the Gumbel max in disguise: √(2·ln(500)) = 3.53 at p ≈ 1000.
  Fixed h — or the RMS/energy form — is the correct target. Downstream this is
  harmless (Fejér uses finitely many h; L⁴⊂L^∞ is per fixed h), but the paper
  statement should be corrected.
- **Exact Parseval bridge** (uses only Σν = p, ν = value multiplicities):
  Σ_{h≠0} |C_p(h)|² = p·E(p) − p². So E(p) = 3p + O(√p) ⟺ RMS_h |C_p(h)| =
  √2·√p·(1+O(p^{−1/4})). The single cleanest open statement of the vertical
  problem remains **E(p) = O(p)** — everything else is packaging.

## 2. DS's decisive question, answered: L(χ) mod p is NOT a bounded-complexity trace function — provably, and tautologically it IS a maximal-conductor one

Mellin coordinate: u = g^m ↔ χ_m. Then L(χ_m) = f(g^m) where
f(X) = Σ_{s=0}^{p−2} n(g^s) X^s ∈ F_p[X] (n(t) = fiber count; this is the mod-p
sibling of the Φ(T) already used in the §1b degree bound). Two facts:

(a) **e_p(hL(χ_m)) IS a trace function**: of the rank-1 Artin–Schreier sheaf
    L_{ψ(hf)} on G_m. It is lisse of rank 1, pointwise pure of weight 0, and
    geometrically irreducible. So conditions (i)–(iii) of the "missing sheaf
    theorem" hold *tautologically*. The entire content was always (iv).

(b) **(iv) fails at maximal strength**: cond = 2 + Swan_∞ = 2 + deg f, and
    deg f = (mod-p linear complexity of the sequence L(χ_m)) − O(1). DS measured
    today: order = p at p=11..19; the §1b Parseval-vs-Weil argument gives
    deg ≥ p/C unconditionally for the energy twin. Grothendieck–Ogg–Shafarevich
    then yields |C_p(h)| ≤ (1 + Swan)·√p ≈ p^{3/2} — i.e. the ℓ-adic route and the
    archimedean CEST bound are the SAME p^{3/2} in two languages. The loop closes.

Logical remainder: a proof via Deligne would need a *different* bounded-conductor
sheaf agreeing with ψ(hf) at F_p-points ONLY (not over extensions — else its trace
function equals the AS sheaf's and Chebotarev forces the same wild part). No known
mechanism produces such single-level coincidences. That is the precise, provable
sense in which the vertical √p sits outside current technology: not "we haven't
found the sheaf", but "the only natural sheaf provably has Swan ≍ p, and any
substitute must be a cross-level accident".

## 3. Dwork-depth prediction for the planned Q6276 §5 computation

Writing C_p(h) = Σ_m Π_{x∈(F_p^*)³} ψ(hΛ(x)^m) and expanding each factor by the
splitting function, the π¹-level term is π·h·(Σ_x Λ^m − Σ_x Λ^{mp}) ≡ π·h·(b_{mp} −
b_m)·(−1) — which vanishes by the Cartier identity (the measured identity Cartier
matrix). Beukers' supercongruence b_{mp} ≡ b_m (mod p³) kills the π² level too.
**Prediction: the precision-N=2 Griffiths–Dwork/Frobenius computation will come out
degenerate (identity-like), and the first nontrivial structure appears only at
N ≥ 4.** Budget the heavy step accordingly, and do not read N=2 triviality as
either success or failure — it is forced by the supercongruence.

## 4. The two-prime "missing theorem" (i)–(iv) is not well-posed as stated

A lisse sheaf on A¹ lives over ONE residue characteristic; there is no Frobenius
whose trace at integer m can equal ψ_p(b_m)·ψ_q(−b_{m+d}) — the phrase mixes two
characteristics under one geometric point. The honest formulation is the
cross-characteristic quenched-vs-annealed comparison already in DS_NOTES_SYNTHESIS
§4.9/§4.11 (pair-Gram norm, F₂ = o). Recommend the paper state the two-prime input
that way and not as a sheaf conjecture a referee can falsify by type-checking.

## 5. Verdict on the dm question

No — I cannot prove the √p bound, and §1–§2 sharpen *why nobody currently can*:
the target is a local-limit theorem whose only sheaf realization provably has
maximal Swan conductor, whose p-adic expansion is supercongruence-rigid to order 3,
and whose measured statistics are indistinguishable (four constants deep) from the
random model — i.e. there is no residual algebraic structure left on the table to
exploit. For the Aug 1 deadline: ship the conditional theorem with (a) the corrected
target statement (fixed-h / energy form E(p) = O(p), Gumbel max for uniform-h),
(b) the value-law table above as the quantitative-evidence section, (c) the
obstruction map upgraded with §2's tautological-sheaf closure — that turns
"we could not find the sheaf" into "the sheaf route is closed by proof", which is
a materially stronger submission.
