# DS note: additive energy + 2D anti-concentration — Q6124 verdict + verification

Date 2026-07-31, DS window. Source: ChatGPT Q6124 (ds2), archived
`chatgpt-answers/Q6124-68cc7066.md`.

## Setup
I = {(r,p) : p prime, r<p, p|b_r}, H(n) = #{r+p=n}. Want max_n H(n) = o(n/log n).
Energy E = Σ_{N<n≤2N} H(n)². Trivial bound H ≤ (1/2+o(1))n/log n.

## Verdict (ChatGPT)
1. **No existing theorem** controls the diagonal energy from the stated fibre
   hypotheses — and NO such theorem can exist: there is an explicit reflected
   "star" (a_p = 2N−p, a_p' = p−1−a_p, all hitting level n₀=2N) with
   H(n₀) ≍ N/log N satisfying bounded fibre, reflection, no-consecutive.
   Even 1D-uniform marginals can be preserved by compensating the star with
   symmetric fibres (density h(u) = (1−cg(u))/(1−c)) while keeping a spike of
   height c·N/log N. **So 1D uniformity of r/p is NOT the missing theorem;
   a JOINT statement is required.**
2. Equivalence refinement: E = o(T²) ⟺ pointwise is NOT unconditional from the
   current p^{2/3}. The right object is the factorial second moment
   F₂ = Σ(H(n))₂ (= E − M_N), and the energy condition is equivalent to the
   pointwise condition only after a bounded-average z_p theorem.
3. Best route: prove the **no-line-mass** statement — the normalized incidence
   measure μ_N (in (p/N, r/N)) has no line-supported singular component
   (sup_t μ_N{|x+y−t|≤η} → 0). This is "mesoscopic anti-concentration",
   weaker than full Poisson but sufficient for pointwise. Fourier core:
   Σ_p (1/p) Σ_{a≠0} Σ_r Φ e_p(a b_r) = o(N/log N).
4. All large-sieve / inverse-sieve / BFI-Zhang-Maynard toolkits have the wrong
   quantifiers (fixed set vs moving Z_p; prime is modulus AND family parameter).

## DS verification (2026-07-31)
From the (p,z) data bank (149,112 pairs, p ≤ 1,999,993):
1. **Strip/line-mass test PASSES**: sup_t μ_N{|x+y−t|≤η} = (1.2–1.35)·2η·log2,
   i.e. O(η), stable across N = 10⁵, 2·10⁵, 4·10⁵. No line atoms — a genuine
   star (H(n₀)≍N/log N) would give normalized strip mass ≍1, NOT O(η). So the
   star mechanism is empirically ruled out. [Caveat: O(η) is also just the
   first-moment/average; the pointwise content is carried by max H(n).]
2. **max_n H(n) = 3 for all n ≤ 8·10⁵** (extends Claude4.6's n ≤ 4·10⁵);
   ratio to n/log n ~ 1e-4. #levels with H≥2 grows ≈ N/log²N (183,262,553,981
   at N=10⁵..8·10⁵).
3. **Exponential sums are random**: |Σ_{r<p} e_p(a b_r)| ~ √p for a≠0
   (measured p=101,211,503,997, a=1,2,3,(p−1)/2). No algebraic cancellation →
   the additive-character core of the Fourier route is essentially tautological;
   a saving must come from joint/cross-p structure.

## Key structural insight (DS derivation)
The no-line-mass statement is EXACTLY equivalent to the pointwise bound (not a
weakening): a point mass in the projected measure ν = (π)_*μ (π = x+y) is
precisely a star. The localized exceptional-set theorem (O_ε(1) exceptional
levels per window of length c_ε N/log N) makes exceptional levels SPARSE but
does NOT rule out a single heavy level: each exceptional level can carry up to
(3/2)N/logN mass (normalized O(1)), so sparse-exceptionals alone give strip
mass O(η log N), not O(η). The gap to pointwise is exactly the missing
"no single level is heavy" statement.

## Follow-up dispatched
Q6129 (ds2, 2nd round): strongest PROVABLE pointwise bound today, easiest
constant < 1/2, and any theorem on prime-divisor-count of a single holonomic
integer in a short window.

## Relevance to the campaign
- Re-routes the attack: not energy/higher moments (no Markov shortcut), but
  joint equidistribution / no-line-mass, which is the pointwise statement.
- Confirms the "compensated star" is the ONLY structural way the conjecture
  could fail pointwise — and it's empirically absent.
- Consistent with Q6125: a genuinely new defining-characteristic
  equidistribution/local-limit input is needed (Dwork–Gross–Koblitz complexity,
  or vertical Apéry large sieve).
