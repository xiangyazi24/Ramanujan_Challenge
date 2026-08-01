# Two-prime Weil correlation — the potential proof path (2026-07-31)

## Discovery
Numerical computation reveals that both the VERTICAL complete exponential sum
and the TWO-PRIME shifted correlation of the Apéry sequence satisfy Weil-type bounds.

### Data 1: Vertical complete sum
```
C_p(h) = Σ_{a=0}^{p-1} e(h·b_a/p)
|C_p(1)| ≤ 3.27·√p   for ALL 166 primes p ≤ 1000
Average |C_p(1)|/√p = 1.26
```

**[CORRECTED 2026-07-31 Fable — see FABLE_NOTES_vertical_value_law.md]** The right
reading of Data 1: fixed-h |C_p(h)| is Rayleigh with mean 1.2533√p (measured 1.261 —
the 1.26 above IS this constant); max over h is Gumbel √(2p·ln(p/2)), NOT C·√p
(the 3.27 ceiling is the Gumbel max at p ≤ 1000: √(2 ln 500) = 3.53). So the
uniform-in-h form "|C_p(h)| ≤ C√p for all h" is FALSE asymptotically; the correct
targets are the fixed-h bound and the energy form E(p) = Σ_a N_p(a)² = O(p)
(exact Parseval bridge: Σ_{h≠0}|C_p(h)|² = p·E(p) − p²).

### Data 2: Two-prime shifted correlation
```
Corr(p,q,d) = Σ_{m=0}^{M-1} e(b_m/p - b_{m+d}/q),   d = |p-q|, M = min(p,q)-|d|-1
|Corr(p,q,d)| ≤ 2.09·√M   for ALL 127 tested prime pairs (p,q ≤ 709)
102 pairs in [200,600] tested, max ratio 2.09, 1 pair over 2.0
```

## The proof chain (if the two-prime bound is proved)

1. **Two-prime Weil bound** (DATA 2): |Σ_m ψ_p(b_m)·ψ̄_q(b_{m+d})| ≤ C√M
2. → **4th-moment bound**: Σ_{N<n≤2N} |S_h(n)|^4 ≤ C'·N^3/log^2 N
   (cross terms in the 4th-moment expansion are O(N^{3/2}) per pair)
3. → **Pointwise bound**: max |S_h(n)| = O(N^{3/4}/log^{1/2}N) = o(N/log N)
4. → **Fejér kernel**: H(n) ≤ P_n/K + (1/K)Σ |S_j| = o(P_n)
5. → **The GCD conjecture**: log G_n = o(n)

## The precise missing theorem (from Q6261 + Q6262)

**[CORRECTED 2026-07-31 Fable]** As stated below this is NOT well-posed: a lisse
sheaf on A¹ lives over ONE residue characteristic, and no Frobenius at a point m
can produce the two-characteristic value ψ_p(b_m)·ψ_q(−b_{m+d}) (m is an integer
in an interval, reduced mod p AND mod q simultaneously). The honest form of the
two-prime input is the cross-characteristic quenched-vs-annealed comparison
(DS_NOTES_SYNTHESIS §4.9/§4.11: pair-Gram norm bounded, F₂ = o(N²/log²N)).
Keep the original text below as the record of the (type-incoherent) sheaf dream:

**Theorem (to prove — NOT well-posed as stated, see correction above).** For every
pair of distinct large primes p, q and every fixed
shift d, there exists a lisse sheaf G_{p,q,d} on A^1 such that:

(i) Tr(Frob_m | G_{p,q,d}) = ψ_p(b_m)·ψ_q(-b_{m+d})

(ii) G is pointwise pure of weight 0

(iii) G is geometrically irreducible (monodromy acts irreducibly after base change)

(iv) cond(G) = rank + Σ_x(1 + Swan_x) = O(1) independent of p, q

Then Deligne's Riemann Hypothesis gives the needed bound.

**[ADDENDUM 2026-07-31 Fable — one-prime version resolved tautologically]** For the
VERTICAL (one-prime) analogue the sheaf exists and is unique-in-effect: with
u = g^m and f(X) = Σ_s N_Λ(g^s)X^s ∈ F_p[X], the phase e_p(hL(χ_m)) IS the trace
function of the rank-1 Artin–Schreier sheaf L_{ψ(hf)} on G_m — pure of weight 0,
geometrically irreducible, i.e. (i)–(iii) hold for free. The content was always
(iv), and (iv) FAILS provably: Swan_∞ = deg f = mod-p linear complexity − O(1) ≈ p
(DS measured order = p; §1b degree bound gives deg ≥ p/C unconditionally). GOS then
returns exactly the archimedean p^{3/2} — the ℓ-adic and archimedean routes are the
same bound in two languages. Any Deligne-style proof would need a DIFFERENT
bounded-conductor sheaf agreeing with ψ(hf) on F_p-points only (not over
extensions); no known mechanism produces such coincidences.

## The construction route

The Apéry numbers b_n = Σ C(n,k)^2 C(n+k,k)^2 are the diagonal of a rational function
(Furstenberg-Deligne theory). The associated Picard-Fuchs/Gauss-Manin sheaf provides the
ℓ-adic object. The two-prime version would be a tensor product of two independent
Apéry period sheaves (one mod p, one mod q), with the shift d handled by a Tate twist or
translation functor.

The decisive step: prove the tensor product has no trivial component
(geometric irreducibility = "two-prime shifted independence of the Apéry motive").

Relevant technology: Katz "Exponential Sums and Differential Equations" (1990),
Katz "Gauss Sums, Kloosterman Sums, and Monodromy Groups" (1988),
Adolphson-Sperber toric exponential sums, Fu-Wan polynomial recursion sums.

## Extended attack session (Q6261-Q6271): routes explored

| Q# | Route | Verdict |
|----|-------|---------|
| Q6261 | Two-prime correlation as Deligne trace-function theorem | ✅ Framework correct, missing sheaf construction |
| Q6262 | ℓ-adic sheaf construction for Apéry | ✅ Precise theorem identified; Apéry is Calabi-Yau, not rigid |
| Q6263 | Katz rigid local systems | Apéry is NOT rigid; NOT hypergeometric |
| Q6264 | Elementary (Stepanov, collision counting) | All reduce to the same sheaf problem |
| Q6265 | Katz rigidity / middle convolution | NOT applicable; non-rigid |
| Q6266 | Hadamard reduction to ζ(2) | b_n ≠ u_n²; different motivic structure |
| Q6267 | Finite-dimensional matrix realization | b_a ≠ Tr(M^a) for bounded-rank M |
| Q6269 | Coefficient extraction via ℓ-adic Fourier transform | ✅ Standard operation (Deligne-Laumon) |
| Q6270 | Fourier constraint from recurrence | Gives Picard-Fuchs ODE, not bound |
| Q6271 | Dwork unit-root = period coefficients? | **FALSE**: α_p(t)=F(t)/F(t^p)≠F(t) |

## The precise gap (after 20+ rounds)

The ℓ-adic Fourier transform CAN extract coefficients from a trace function (Q6269).
The Dwork theory gives a sheaf for the UNIT ROOT α_p(t) = F(t)/F(t^p).
But b_a mod p = [t^a]F(t), NOT [t^a]α_p(t).

**Missing theorem**: realize the PERIOD COEFFICIENTS (not unit root) as trace functions
of a geometrically controlled ℓ-adic sheaf. This is a problem in p-adic Hodge theory:
comparison between de Rham period sections and ℓ-adic Frobenius traces at the
coefficient level.

## Moment bridge reformulation (Q6278-Q6282, verified p=5,7,11)

**Key identity**: b_m ≡ -L(χ_m) (mod p) for m = 0,...,p-2, where:
- L(χ) = Σ_{a ∈ F_p^*} χ(a) N_Λ(a) is the twisted point count
- χ_m(a) = a^m is the m-th multiplicative character
- N_Λ(a) = #{x ∈ (F_p^*)^3 : Λ(x,y,z) = a}
- Λ(x,y,z) = (1+x)(1+y)(1+z)((1+y)(1+z)+xyz)/(xyz)

**Consequence**: C_p(h) = Σ_χ ψ(-h·L(χ)) + O(1)

The vertical Weil bound (fixed h; uniform-in-h is false, see Data 1 correction)
|C_p(h)| ≤ C_h√p is EQUIVALENT to: |Σ_χ e(-h·L(χ)/p)| ≤ C_h√p

This is a sum of additive characters of L-function values — a "mod p of Frobenius
trace" equidistribution problem. No existing theorem (Katz, Deligne, Sato-Tate)
gives this: archimedean equidistribution controls |L(χ)|/p^{3/2}, not L(χ) mod p.

**Status**: Numerically verified (|C_p| ≤ 3.27√p for 166 primes), but unproved.
The reformulation is new and connects the Apéry coefficient problem to families
of L-functions in a precise way.

## Additional routes explored (Q6285-Q6293)

| Q# | Route | Verdict |
|----|-------|---------|
| Q6285 | Bezout / bounded-degree orbit curve | NO: orbit degree grows with p |
| Q6288 | Katz discrepancy → mod-p residues | NO: test function frequency ~p^{1/2} too high |
| Q6289 | Modularity / Hecke structure | Identifies object but doesn't close additive twist |
| Q6291 | RS/MDS coding theory | NO: MDS controls Hamming distance, not Fourier |
| Q6293 | Weyl differencing on L(χ_m) | NO: "degree" = p-1, all frequencies present |

## Collision count data (46 primes ≤ 199)
Σ N_t² / (2p) ∈ [0.90, 2.01], mean ≈ 1.4. Max multiplicity ≤ 10 (mostly ≤ 6).
Random prediction: ratio → 1.0. Observed: close to random with mild excess.

## DS data (p=11,13,17,101,...,503)
- |C_p|/√p bounded: max 3.22 at p=307, cluster ~1.0
- ~~L(χ_m) takes ~√p distinct residues mod p (folding from range Cp^{3/2})~~
  **[CORRECTED 2026-07-31 Fable]** The distinct-value law is LINEAR, not √p:
  D(p) = (1−e^{−1/2})p = 0.3935p (measured 0.3963p ± 0.018 over 164 primes;
  p=997: 392 obs vs 391.8 pred). It is the Poisson occupancy law conditioned on
  the reflection FE b_r ≡ b_{p−1−r}, with ZERO archimedean-folding content.
  The p=11..19 data could not discriminate (the two laws cross at p ≈ 11).
- ~~Local-limit: max multiplicity ~√p, confirmed for p=11,13,17~~
  **[CORRECTED]** max multiplicity grows like the Poisson-max ~log p/log log p
  (measured 4→8 from p=101→997), far below √p.
- Linear-recurrence order of L(χ_m) = p (growing, not bounded-complexity) —
  consistent with Swan_∞ = deg f ≈ p in the Artin–Schreier tautology above.
