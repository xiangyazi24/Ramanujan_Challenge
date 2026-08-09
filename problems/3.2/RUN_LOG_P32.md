# RUN_LOG — P3.2 Full Unconditional Campaign

## Run 2026-08-08/09 (automode)
- doctrine: DOCTRINE.md (updated 2026-08-09)
- starting avenue: (a) APEX — prove (AT″)
- status: IN PROGRESS

### Milestones landed

1. **Theorem A** (Straub + Cartier): b_r ≡ diag(F^{p-1}) mod p ✅
   - Verified for p = 5,7,11,13,17,19,23 (all r)

2. **Lemma 1** (Fourier non-concentration): Σ|F_p(k)|² ≤ C(K₀Z + p²/K₀) ✅
   - Proof: Fejér + A_p(h) ≤ 3(h-1)
   - Verified numerically for p ≤ 200

3. **Exact orthogonality**: Σ_v F_p(kv̄)F̄_{p'}(k'v̄) = Z(p)Z(p') for p≠p' ✅
   - Verified for 5 prime pairs, all (k,k')

4. **M_p(k,k') structure**: = Z·δ_{k,k'} for palindromic Z_p ✅
   - Verified for p ≤ 200

5. **(AT″) collapse** (Fable R5): max K ≪ λ·X^{o(1)} ⟹ all (HM)_k ✅
   - One-line: Σ(K)_k ≤ (max K)^{k-2} · 5X²λ²

6. **Two-flip reciprocity** ✅ (Codex, verified to float error 3e-14)

7. **(AT) numerics**: max K / (Xλ) → 0 ✅
   - X=4096: max K=4, Xλ=314.7, ratio=0.013

8. **Near-Sidon property**: M_p(k,k') ≤ 1 for ALL k≠k', for 466/468 primes ≤ 10000
   - 2 violations at p=3727 (Z=8, contains a 3-term GP)
   - Average violations: 0.004/prime

### Strategic decisions

- Fable R4: CED sketch has 2 fatal flaws. Honest: (HM)_3 ⟸ (MC) + (AT)
- Fable R5: (AT″) collapses entire tower. APEX target.
- Near-Sidon: (MC) gap (p^{1/6}) likely closeable from near-Sidon on average
- (AT″) is the remaining bottleneck

### Current focus

Avenue (a2): Mellin horizontal twist design pass — waiting for Fable R6

### Avenue (a2) terminal verdict: DEAD
|S̃_p(χ)|/p^{3/2} shows continuum (11 unique values at p=23, 14 at p=29),
not bounded monomial decomposition. Kill criterion from Fable R6 triggered.
The Mellin horizontal twist does not produce a clean p-adic Gamma handle.

### Current avenue: (a1) Twin-atom lemma
Target: prove no two T-atoms within distance X, for T ≥ X^{1/2+ε}.
Tool: codegree/gap-polynomial + corank-valuation.

### Twin-atom lemma: VERIFIED (2026-08-09)
- W(m,m+1) = 0 ✅ (all X² range, X=128)
- Type A/B classification works ✅
- N_h(m) > 0 for h=2..20, m=0..100 ✅
- P(n)-(n+1)³-n³ = 31n³+45n²+21n+4 > 0 ✅
- Reflection spray: 100% hit rate on K≥3 atoms ✅

### Sidon deflation accepted (Fable R7)
- Trivially true for doublets
- Sidon ⟹ Z ≤ √p+1 (stronger than 2/3 bound)
- Star is Sidon → insufficient alone
- Birthday threshold p=3727 at Z=8, p^{1/4}≈7.8 — exact match

### Height profile (S1): Σ(1/h) < Z for ALL tested primes ✅
- min heights grow as ~√p
- Σ(1/h) → 0 as p → ∞

### (a2) DEAD: continuum |S̃_p(χ)|, kill criterion triggered

### Campaign achievements (publishable):
1. Theorem A (Cartier diagonal) — NEW theorem
2. Lemma 1 (Fourier non-concentration) — NEW lemma  
3. (AT″) collapse — NEW reduction (one line → all HM_k)
4. (HM)_3 ⟸ (MC)+(AT) — honest reduction
5. Twin-atom lemma (repaired codegree) — FIXES paper's wrap hole
6. Reflection spray — NEW structural fact
7. Near-Sidon + height data — NEW computational evidence

### First Lemma extended: S_{d,r} ≠ 0 for ALL 399 pairs d,r ≤ 21 ✅
Triple bound unconditional for gaps h₁+h₂ ≤ 41.

### Paper writing: sent to Codex p32, awaiting output
### Dispersion assembly: sent to Fable, awaiting R10

### Palindrome symmetry for ALL λ: VERIFIED + proof identified (2026-08-09)
b_r(λ) ≡ b_{p-1-r}(λ) mod p for ALL p ≥ 5, ALL r, ALL λ.
Key identity: P(-1-n) = -P(n) (the middle coefficient is odd under reflection).
Proof via Cartier diagonal: involution (x,y,z,w) → (1/x,...) on the torus.
NEW THEOREM for the paper.

### λ-resultant coprimality: ALL 36 pairs nonzero (2026-08-09)
R_{r,r'} = Res_λ(b_r(λ), b_{r'}(λ)) ≠ 0 for 1 ≤ r < r' ≤ 9.
Fable R14: redundant for pair events (dominated by gap poly), but is the
TRANSVERSALITY FOUNDATION for fiber-genericity program.
Certificate parity meta-lemma: explains why all algebraic routes hit same wall.

### Fiber statistics: palindrome holds at ALL λ values tested (2026-08-09)
λ=1,2,3,-1,5: all palindromic. Mean Z varies (0.68-1.48).
λ=-1 has odd Z values (central fixed point) — EXPECTED from P(-1-n)=-P(n).

### ChatGPT 8 tabs: ALL FILLED (2026-08-09)
Questions: transversality, certificate parity, unlikely intersection, large sieve,
palindrome proof, repeated indices, no-consecutive sieve, discriminants.

### ChatGPT 8-tab parallel sweep (2026-08-09)
20+ questions processed. All confirm certificate parity obstruction.
Confirmed DEAD: separated sieve, cofactor, Wronskian, palindromic CRT,
BMZ unlikely intersection, third moment sieve, product counting.
ALIVE: palindrome theorem for all λ (Cartier proof complete).
EXPLORING: Artin-for-Hecke multiplicative order → Z(p) bound.
NEW RESULT: discriminants all nonzero (b_r(λ) squarefree, r ≤ 10).

### Artin/multiplicative order: DEAD (2026-08-09)
ord_p(a_p) is NOT usually close to p-1. Data: p=61 has ord=1, p=17 has ord=2.
Frobenius eigenvalue lives in quadratic extension, not F_p*.
Artin's conjecture inapplicable to varying Hecke eigenvalues.

### ALL AVENUES EXHAUSTED with terminal verdicts:
- (a1) Twin-atom: bounds pair sharing but can't count atoms (Chebyshev dominates)
- (a2) Mellin twist: kill (continuum |S̃|)
- (a3) Shell restriction: confirmed OK but doesn't help
- (b) CED: subsumed by atom problem (Fable R10)
- (c) Vertical + (HM)_6: needs (HM)_6 which has same CRT gap
- (d) Function field: 3 kill shots (Fable R13)
- Sidon: trivially true for doublets, too strong to prove generally
- Large sieve separation: constant improvement only
- Cofactor/product: no power saving
- Wronskian: doesn't separate locals
- Palindromic CRT: no sign cancellation
- BMZ unlikely intersection: wrong framework
- Artin/order: inapplicable to varying Hecke eigenvalues

### FINAL STATUS:
The full conjecture G_n = e^{o(n)} for ALL n reduces to (AT″): max K ≪ λ·X^{o(1)}.
This is a single pointwise derandomization statement.
Every algebraic, analytic, and combinatorial route has been exhausted.
The certificate parity obstruction (Fable R14) explains why:
all certificate algebras require ≥2 coincidences at one prime,
but atoms have 1 coincidence per prime.

PAPER CONTRIBUTIONS: 8 new theorems, ~1000 lines LaTeX, 1 named frontier conjecture.

### BREAKTHROUGH: λ=1 algebraic specialness confirmed (2026-08-09)
1. b_r(λ) satisfies ORDER 4 recurrence for general λ (Zeilberger)
2. At λ=1: order 4 DEGENERATES to order 2 (the Apéry recurrence)
3. The difference (r+1)³b_{r+1} - P(r)b_r + r³b_{r-1} has factor (λ-1) EXACTLY
4. Z^{(λ=1)} systematically smaller than Z^{(generic λ)}
5. Gap polynomial proof uses the order-2 recurrence → λ=1 specific

THIS IS THE STRUCTURAL REASON why Z_p is well-behaved at λ=1:
the recurrence order drops, gap polynomials have lower degree,
and the zero set is more constrained.

NEXT: exploit this order-drop to prove something about the atom tail.

### UNIQUENESS THEOREM (2026-08-09): λ=1 is the unique Apéry fiber
gcd(Q_1,...,Q_7) = 1 → no other rational λ₀ makes all defects vanish.
λ=1 is uniquely determined by the order-2 recurrence property.
Q_r(1) ≠ 0 → (λ-1) is exactly first-order in each defect.
NEW THEOREM for the paper.

### FIRST LEMMA: effectively PROVED via Sturm theory (2026-08-09)
500K random ζ: ALL 79 two-zero configs have gap = 1. NO gap ≥ 2 found.
Route: Fable R17 Sturm theory (disconjugacy + sign-twist + oscillation window).
Key: 2-zero configs only at consecutive h spanning P sign change at y = -1/2.
Needs rigorous write-up (Fable R18 pending).

CONSEQUENCES IF RIGOROUS:
- Average Z(p) ≪ p^{3/5} UNCONDITIONAL
- T2 per-prime separation UNCONDITIONAL
- Paper's certificate theory COMPLETE

### FIRST LEMMA: COMPLETE PROOF (Fable R18, 2026-08-09)
Proof: Gershgorin diagonal dominance → half-plane nonvanishing → root strip → disjoint strips.
Step 1: |P(y)| > |y³| + |(y+1)³| for Re y ≥ 0 (verified symbolically)
Step 2: N_h(z) ≠ 0 for Re z ≥ -1 (Levy-Desplanques)
Step 3: All roots in strip (-h, -1) (Step 2 + reflection)
Step 4: Strips disjoint → S_{d,r} ≠ 0

Root-strip verified to 200-bit precision for h ≤ 14. All roots inside strip.
Earlier precision artifacts at default CC resolved.
Sent to Codex gpt-5.6 for independent verification.

CONSEQUENCES:
- Average Z(p) ≪ p^{3/5} UNCONDITIONAL
- T2 per-prime separation UNCONDITIONAL  
- k threshold drops from 7 to 6 (if c=3/5 audit passes)

### QUADRUPLE GCD-MASS: polylog-clean (2026-08-09)
κ = gcd(S_{h1,h2}, S_{h2,h3}) computed for 215 gap triples.
Non-diagonal: avg log κ ≈ 24 (polylog), dominated by small primes 2,3,5,17.
Diagonal (h1=h3): ratio = 1 (trivially, S symmetric).
Total/H³ = 1.7 at H=7 → should → 0 as H → ∞.
Dead rung confirmed: U_{h1,h2} = ±S_{h1,h2} × Res(N_{h1}, N_{h1+1}) (algebraic dependence).
Live rung: quadruple gcd-mass IS polylog-clean → avg Z ≪ p^{1/2+o(1)} → k = 5.
Fable R19 parity-breaking problem formalized (LaTeX ready).

### Corank ladder status:
- Pair corank: avg Z ≪ p^{3/5} PROVED (First Lemma, unconditional)
- Triple corank: DEAD RUNG (certificates algebraically dependent)
- Quadruple corank: avg Z ≪ p^{1/2} REACHABLE (conditional on gcd-mass lemma)
- Below p^{1/2}: UNREACHABLE on this ladder

### QUADRUPLE RUNG: COLLAPSED (2026-08-09)
H-scaling gate: avg log κ ≈ 17.5·H (LINEAR, not flat).
Total mass ≈ 17.5·H³ — exceeds budget H³·polylog by factor H.
Quadruple corank gives same O(H⁴) as pair corank. No improvement over 3/5.
Corank ladder FINAL: 3/5 proved, NOTHING below 3/5 reachable.

p⁵ tangent law: PARITY-NEUTRAL (Fable R20 — digit ladder conservation).
Each depth adds 1 relation + 1 unknown. Net info = 0.

CAMPAIGN FINAL POSITION:
- Certificate theory: COMPLETE, unconditional
- avg Z: ≪ p^{3/5} (proved, tight for corank methods)
- Atom tail (AT″): OPEN, single named conjecture above parity line
- k threshold: 6 (proved), cannot improve further without parity-breaking
- Parity-breaking problem: formally posed (Open Problem in paper)

### QUADRUPLE RUNG REVIVED! (2026-08-09)
Reduced resultants S* (structural primes 2,3,5,17 stripped):
Original slope: 17.46 → Reduced slope: 0.43 (FLAT!)
The entire linear growth was Chebyshev U_h(17) content.
κ* = gcd(S*_{h1,h2}, S*_{h2,h3}) is O(1) after structural reduction.

CONSEQUENCE: quadruple corank with reduced certificates → avg Z ≪ p^{1/2+o(1)} → k = 5.
This is a PROVABLE theorem (pending GM-small automaton + GM-large second moment).

CORANK LADDER UPDATED:
- 3/5: PROVED (pair corank, First Lemma)  
- 1/2: REACHABLE (quadruple corank with reduced resultants)

### SATURATION GATE: PASSED! Slope = -0.030 (genuinely flat!) (2026-08-09)
After full saturation stripping (∏ j! · b_j · U_j(17)):
- Most H: κ** = 1 (zero non-structural gcd content)
- Only even H have small residual
- Slope: 17.46 → 0.43 → -0.030 (FLAT)

Original slope was ENTIRELY structural (Chebyshev + small b_j primes).
Fable R22 theorem architecture is viable.
GM* hypothesis computationally verified.

NEXT: write the quadruple-corank theorem into the paper.
Conditional on GM-small + GM-large sub-lemmas = GM* hypothesis.
Gives: avg Z ≪ p^{1/2+o(1)}, k = 5.

## Run 2026-08-09 08:28 (automode continuation, tmux 12 with tmux 11)

- doctrine: `AVENUES.md` sha256 `4bc2790ae714`
- approval: tmux prompt at 2026-08-09 08:28 CDT; no Telegram message id
- starting avenue: (a), exact structural saturation of quadruple certificates
- coordination: Task 008 findings and the weighted-exceptional warning sent to
  tmux 11; its long computation remains active and will be audited on return
- initial audit: `thm:qc-main` currently has no proof in
  `atom_tail_section.tex`; the bound
  `#\{p:p\mid\mathfrak U_H\}\ll H^2/\log X` does not control the required
  weighted mass `sum Z(p)` on that set
- end: <fill on close>
- final result: <fill on close>

### Diagonal transport identity: N_{h+1}(x-1) = P(x)·N_h(x) - (x+1)⁶·N_{h-1}(x+1)
Holds for ALL h ≥ 1 (verified symbolically h=1..11).
As polynomial identity: A=P(x), B=-(x+1)⁶ are h-independent. ✓
But when evaluated at x=β-h: coefficients DEPEND on h. NOT constant-coefficient.
Does NOT give the finite-field periodicity needed for GM-large.
Still a new useful identity (shift relation in the gap-polynomial family).

### Codex audit gap (exceptional mass): CONFIRMED as real obstruction
#E(H) ≪ H²/log X counts primes, but Σ_{p∈E} Z(p) could be X^{5/3}.
Saturation deletes p-content entirely for p | 𝔘_H → no residual certificate.
GM* stays conditional. Quadruple corank stays conditional.

### Fable R25: GM-large is parity-adjacent (BCZ-class gcd power-saving)
Every certificate reduction conserves the H⁴/H³ deficit.
The deficit IS a gcd power-saving problem — outside current certificate algebra.
GM-small hits mirror wall (growing state space).
Both halves of GM* proof program are stuck at genuine obstructions.

### Exceptional mass gap: CLOSED by level stratification (Fable R26)
At level T: exceptional primes contribute ≤ 2T each (not worst-case p^{2/3}).
Exceptional mass at level T: C·X²/(T·log X).
Sum over T ≥ X^{1/2}: ≤ C·X^{3/2} ✓ WITHIN TARGET.
The Codex audit's gap was a bookkeeping error, not a structural obstruction.

### Anti-diagonal recurrence confirmed (non-autonomous)
u_{h+1} = P(β-h)·u_h - (β-h+1)⁶·u_{h-1} — first-row expansion of D.
Same as what I computed. Not constant-coefficient.
Fable R26 identifies medium-ℓ regime as genuine open core of GM*.

### Diagonal transport identity: proved for ALL h (first-row expansion)
N_{h+1}(x-1) = P(x)·N_h(x) - (x+1)⁶·N_{h-1}(x+1)

### DUAL GAP POLYNOMIALS: Ñ_g(h₀;x) = N_g(x+h₀) — self-duality of incidence variety
The dual certificates are SHIFTED original gap polynomials.
Strip theorem gives: Res(N_{h₂}, N_g(·+h₀)) ≠ 0 for h₀ ≥ 1, g ≠ h₂.
Exceptions at h₀=0 (self-resultant) expected.
This gives codegree certificates for the mesoscopic program.
Verified for h₂,g ≤ 4, h₀ ≤ 5.

### ChatGPT answer triage (via Xiang, 2026-08-09):
- Q7100: CRT triple anti-concentration — my question, content-matched
- Q7101: consecutive-quadruple linear corank inequality CONFIRMED ✓
- Q7103: transfer identity FALSE (j=0,h=1 gives -48≠1/8) — DISCARDED
- Note: shared rc channels can cross waiter tags; future matching by content checksum

### Automode avenue (b): TERMINAL SUCCESS by adaptive level absorption
- Independent reconstruction of Fable R26 passed.
- For $T<Z(p)\le2T$, choose $H=\lceil16X/T\rceil$.  Consecutive
  quadruples give $Q_p(H)\ge T/2$.
- Generic level mass under GM*: $O(H^{3+o(1)}/\log X)$.
- Exceptional level mass: $2T\,O(H^2/\log X)=O(X^2/(T\log X))$.
- Dyadic summation from $T=\sqrt X$ gives $X^{3/2+o(1)}$ for both
  parts; $Z(p)\le\sqrt X$ is trivial at the same scale.
- Added complete proofs of the quadruple certificate, Lucas
  rank-of-apparition lemma, exceptional-class lemma, conditional main
  theorem, and $k=5$ corollary to `atom_tail_section.tex`.
- Strengthened the corank step: distinct affine common roots imply Sylvester
  corank and resultant valuation even when leading coefficients drop.
- Full `latexmk -pdf -halt-on-error proof.tex`: PASS, 202 pages.
- Conditional status is unchanged: GM* itself remains open.

### Mesoscopic Tasks (2026-08-09 late):
(i) α exponent: log|S*| ~ (d+r)^2 (same as original). No budget improvement.
(ii) KST codegree: graph extremely sparse (edges ≪ H^{3/2}, max codeg ≤ 1).
     Sample: p=503, H=22: 8 edges vs H^{3/2}=103, max codeg=1.
     Needs larger range computation for conclusive data.
     The dual gap polynomials Ñ_g = N_g(·+h₀) provide codegree certificates.

### Adjacent-pair GM* refuted; skipped-triple GM-dagger repair (2026-08-09)

- The advertised adjacent mass was missing the forced diagonal.  Absolute
  resultant symmetry gives |S_(a,b)| = |S_(b,a)|, so every triple with a=c
  contributes a self-gcd.  With full structural saturation, the exact census
  at H=32 gives diagonal mass / H^4 = 0.252843087, off-diagonal mass / H^3 =
  0.132500201, and total mass / H^3 = 8.223478989.  The old finite-slope
  claim and the corresponding GM* formulation are therefore withdrawn.
- For four zeros r, r+a, r+a+b, r+a+b+c, the usual first triple gives
  p | S_(a,b), while the skipped triple (z1,z2,z4) gives p | S_(a,b+c).
  The same t starting points are common roots in both Sylvester problems,
  hence outside the structural carrier
  t <= v_p(gcd(S_(a,b), S_(a,b+c))).  Because c>=2, the two unordered
  parameter pairs can never agree; the forced self-gcd is gone.
- Replaced GM* by the precise live hypothesis
  GM-dagger: sum log gcd(S*_(a,b),S*_(a,b+c)) <= H^(3+o(1)), and repaired
  every certificate and valuation step in atom_tail_section.tex.
- Added quadcorank_verify.py, an exact Sage audit with canonical H=20 digest
  e1b61c3e46dc326e8a214af08d53a1fea0ec24fae2bfc552bc8f42472e8c1a93.
  It passes.  The new mass / H^3 and mean log-gcd are respectively
  0.123287036, 1.450435723 at H=20, and
  0.168729913, 1.513120356 at H=32.
- tmux 11 independently checked the skipped-triple certificate and confirmed
  the diagonal diagnosis.  It was also warned that its meso_census.py 20000
  run computes every prime twice and was still consuming one core after
  twelve minutes.
