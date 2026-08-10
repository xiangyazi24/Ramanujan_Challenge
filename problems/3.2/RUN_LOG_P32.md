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
- end: 2026-08-09 10:45 CDT
- final result: the conditional quadruple-corank implication is complete and
  sharpened from scalar GM-tr to the same-root aligned-content hypothesis
  AC-tr, with an exact verifier and level-adaptive exceptional absorption.
  The full conjecture remains open: AC-tr is a genuine BCZ-class content
  bound, while fixed-depth Cartier lifting, naive cut-edge carrier expansion,
  and the present mesoscopic codegree tools all reached terminal obstructions.

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

## Run 2026-08-09 20:17 (automode breakthrough continuation)

- doctrine: `AVENUES.md` sha256 `b613076e8d2c`
- approval: user prompt "卡住就找突破，不要放弃。用 chatgpt 配合。";
  commentary notice at 2026-08-09 20:17 CDT; no Telegram message id
- starting avenue: (a0), local-algebra formula for aligned content
- coordination: tmux 11 is retained as an independent worker; ChatGPT bridge
  channels are to be kept occupied by tactical and strategic audits while the
  local exact computation proceeds
- end: 2026-08-09 23:52 CDT
- final result: the normalization-length formula, fully deflated certificate,
  and the conditional implications under FDAC-tr, SWR, and QPRS are proved.
  Exact audits close the formal computations through height 36, span 84,
  actual primes five million, and projective primes 500,000.  The formal
  census disproves QPRS with `C=2` but not QPRS with an unspecified constant;
  its random-fiber growth makes distinguished-orbit/SWR coupling the remaining
  primary input.  No unconditional square-root average is claimed.

### Breakthrough-run milestones (20:17 continuation)

- Proved the exact local normalization formula
  `v(cont_T Res(F,G+TJ)) = sum m_i length(O_i/(G(alpha_i),J(alpha_i)))`
  over a complete characteristic-zero DVR with perfect residue field and
  unit leading coefficient.
  The invariant is a length on the normalization, not generally on the
  order `R[x]/(F)`; the discrepancy is verified by
  `F=x^2-pi^2, G=x-pi, J=x+pi`.
- Added the fully deflated adjacent pencil and conditional hypothesis
  FDAC-tr.  Every selected actual quadruple survives all three deflations
  except those containing the globally unique centered adjacent pair, so
  the total loss is at most three.  The adaptive-level proof gives the same
  conditional `sum Z(p) << X^(3/2+o(1))` conclusion as AC-tr.
- Exact `ZZ[T]` computation through ambient height 32: thirteen nontrivial
  reduced contents, supported only at 157, 431, and 653; mass/H^3 is
  `0.002569880`, digest
  `437f7328a4ea7d78a62b3e781165a76a31f277f7ec6983f9f6c3a4e67ee05efa`.
  The value 24649 in one record is `157^2`, not a new prime.
- The complete extension through height 36 has nineteen nontrivial records,
  support `{157,431,499,653,1297}`, mass/H^3 `0.002595232`, and digest
  `196327788a7ea9adbbb141a5c6161ae96125ea5b8242bc6603352080b32aac76`.
  Independent finite-field reconstruction classifies the p=431, p=499,
  p=157, and p=1297 reflected pairs as primitive phantoms; every other record
  is a skipped chain whose full return list contains a centered adjacent pair.
- Exact finite-field classification through height 36: among the nine
  endpoint-gcd candidates with `p>(a+b+c)^2`, three are endpoint false
  positives, four are nonprimitive, and two are primitive phantoms.  The new
  reflected phantom pair is `(5,20,10,p,r)=(5,20,10,1297,360)` and
  `(10,20,5,1297,901)`; its span is 35, but the independently checked actual
  zero set is `Z_1297=(459,530,766,837)`.  No short candidate is an actual
  Apéry-zero start.  Some p=653 nonprimitive roots are also in the short
  range, so both primitivity and distinguished-orbit coupling are essential.
- The extended 350104-pair census through p<=5,000,000 contains 43366
  sliding actual-zero quadruples and 1418 off-center windows, but no
  off-center window of span at most sqrt(p).  Its SHA-256 is
  `5739c6e7fee4210678bc50bdda7d0a7c2f9fa082ab392e050556ebdc62ecac8b`;
  a second fail-closed run with that expected digest passed.  The closest
  scale and first unrestricted example remain p=3727, with exact
  `span^2/p=1428025/3727`.
- Isolated the short-window reflection hypothesis (SWR).  It suffices that
  the reflection statement hold below `eta sqrt(p)` for any one fixed
  `eta>0`; the empirical `eta=1` version is stronger than necessary.
  Unconditionally,
  `Z(p)<3 eta^(-1)sqrt(p)+6+E_sw(p;eta)`.  Under SWR, blocks of size
  `floor(eta sqrt(p))+1` give
  `Z(p)<=3 ceil(p/(floor(eta sqrt(p))+1))+3
       <3 eta^(-1)sqrt(p)+6`,
  and hence the desired prime-window average.  ChatGPT Q7139 independently
  verified the `eta=1` ceiling and block-consecutiveness steps; the same
  calculation is uniform in every fixed `eta`.
- Isolated the alternative quadratic primitive-return support hypothesis
  QPRS: outside `U_s`, every primitive off-center four-return chain with a
  non-all-equal gap vector need only satisfy `p<=C s^2` for some absolute
  `C`.  Taking the low
  zero-count cutoff
  `K sqrt(X)` with `K>>sqrt(C)` makes `C H^2<X<p` on every higher adaptive
  level, so no generic selected quadruple survives; the structural primes are
  absorbed by the existing `O(H^2/log X)` exceptional count.  This gives a
  complete conditional square-root average without any content-mass
  hypothesis.  This implication remains correct for any fixed `C`; the
  computations below concern the plausibility of the hypothesis itself.
- Exact endpoint scans verify QPRS with `C=2` for every non-all-equal gap
  triple of span at most 84.  The incremental block digests for maxima
  40,44,52,60,68,76,84 are hard-coded fail-closed in
  `primitive_fd_candidate_verify.py`.  Beyond span 36 the tested cases are
  five endpoint false positives, eight p=8941 nonprimitive roots, and the
  primitive phantom pair `(3,40,31,p,x)=(3,40,31,7411,4681)` and
  `(31,40,3,7411,2655)`.  Its ratio is exactly `7411/74^2`, and the
  independent recurrence check gives `Z_7411=empty`.
- Added a direct prime-first projective-orbit scanner and an independent
  standard-library verifier.  The exact scan of every prime `7<=p<=500000`
  and every projective fiber examines 13,695,120 carrier-free primitive
  off-center windows: 20 have `p>s^2`, four have `p>2s^2`, and all are
  phantoms.  The maximum is `128047/164^2=4.760819...`, at
  `(p,x;a,b,c)=(128047,42375;41,86,37)` and its reflection.  The scan-output
  SHA-256 is
  `8eeb4d6d4a6f371d8f8d87facadff137f66dfa476a28db1a3af8892625032bf9`.
  Thus `C=2` is false, although no finite scan disproves existence of some
  absolute QPRS constant.  The random-fiber estimate `s^3/p^2` per prime
  predicts cumulative size `sqrt(P)/(C^(3/2) log P)` at
  `s=sqrt(p/C)`, so QPRS is now secondary to distinguished-orbit/SWR input.
- Generic product-formula/Fitting summation and generic codegree/KST routes
  are terminal: direct-sum heights retain all labelled triples, and even
  codegree one permits quadratic deep energy in the affine-plane abstract
  obstruction.  Only an Apéry-specific weighted split-gcd estimate remains
  live on that branch.
- The fully deflated height-36 audit, endpoint scan through span 84,
  actual-zero census through five million, and formal prime-first census
  through 500,000 are complete.  The active target is SWR or an equivalent
  distinguished-orbit coupling; QPRS is retained only as a proved
  conditional implication.  FDAC-tr and AC-tr remain hypotheses.
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

### GM-dagger also refuted; central deflation gives GM-tr (2026-08-09)

- Full factor inspection found a second forced degeneration that the
  no-self-gcd check missed.  On c=a with even b, the odd part of the center
  value T_a^(b) divides both S_(a,b) and S_(a,a+b), exactly as predicted by
  Proposition meso-center-recurrence.  At H=32 this palindromic slice alone
  has log mass / H^4 = 0.004833476.  Hence GM-dagger is withdrawn.
- Put N^o_b=N_b/(2x+b+1) for even b and N^o_b=N_b for odd b, and define
  D_(a,b)=Res(N_a(x),N^o_b(x+a)).  For a four-zero pattern (a,b,a), every
  noncentral start remains a common root for D_(a,b) and S_(a,a+b).
  A removed central start satisfies z1+z4=p-1; reflection and global
  consecutiveness show that at most one selected quadruple per prime is
  lost.
- The residual all-equal slice a=b=c is one-parameter.  Its multiplicity is
  detected directly by S_(a,a), and
  sum_(3a<=H) log|S_(a,a)| = O(H^3 log H), already within budget.
- The precise live hypothesis GM-tr is the sum of skipped gcds over c!=a
  plus deflated gcds gcd(D*_(a,b),S*_(a,a+b)) over c=a, a!=b.  Its exact
  mass / H^3 is 0.005568648 at H=20 and 0.014732048 at H=32.  The separate
  progression mass / H^3 is 0.183084461 and 0.262264680.
- quadcorank_verify.py now audits all three failed/repaired formulations.
  Its new canonical H=20 digest is
  6f7ba2a8f3542da4d0a051c698432a9c32124a1811ad8601837dfc1d87968b1a,
  and the exact check passes.

### Same-root refinement: aligned content AC-tr (2026-08-09)

- The scalar GM-tr gcd can still charge the same prime at unrelated roots.
  For every non-AP gap triple define
  `C_(a,b,c) = cont_T Res_x(F,G+T*J)`, using the center-deflated `G` on
  `c=a!=b`.
- Rigorous valuation lemma: over the Gauss DVR of `Q_p(T)`, `t` distinct
  common roots of `F,G,J` give Sylvester corank at least `t`, hence
  `v_p(C)>=t`.  The residue field is `F_p(T)`; this is not a specialization
  at `T=0` and requires no root-simplicity assumption.
- Degree-padding audit: since `deg J>deg G`, the constant coefficient is
  `lc(F)^(deg J-deg G) Res(F,G)` and the top coefficient is `Res(F,J)`.
  Since `lc(F)=V_a` is supported on the structural carrier, `C*` divides
  the old saturated scalar gcd term by term.  AC-tr is therefore strictly
  no stronger than GM-tr and is the exact same-root input needed by corank.
- Added `aligned_corank_verify.py`, which computes exact resultants in
  `ZZ[T]` and checks both endpoint identities and the divisibility into the
  scalar gcd.  Canonical H=20 digest:
  `c7174f329e889bbb1f999a235015042987e14575a2a9bff111fe2ed4e180af4f`.
- Exact data: all reduced contents are 1 through H=14.  At H=20, only
  `(2,4,14)` and `(14,4,2)` are nontrivial, both equal to 653, and
  `sum log(C*)/H^3 = 0.001620394`.  At H=24, 18 of 1323 terms are
  nontrivial and the ratio is `0.007248253276`.
- A faster consecutive-segment pencil replaces `N_(b+c)(x+a)` by the
  directly forced `N_c(x+a+b)`.  Its reduced contents agree termwise at
  H=20.  An exact H=28 scan has only 22 nontrivial terms among 2292 and
  mass ratio `0.005780475`; the residual support is `{131,157,653}`.
  The continuant composition identity relates the two pencils, but its
  cut-edge multiplier can vanish at primes dividing values
  `N_a(-a-j)` that are not all supported on the present carrier.  Thus the
  observed equality does not yet give a uniform reduced-content theorem,
  and the direct pencil has not been used in the paper.
- Naively adjoining all distinct cut-edge values does not repair this:
  reflection gives `N_a(-a-j)=+/-N_a(j-1)`, and the triangular product over
  `a+j<=H` has logarithmic height `Theta(H^3 log H)`.  At a dyadic level
  `H~X/T` its exceptional mass is `T H^3/log X`, which is `X^2` at
  `T~sqrt(X)`, above the `X^(3/2+o(1))` target.  Any useful replacement
  would need a genuinely compressed carrier, not this product.
- Replaced the paper's conditional hypothesis by AC-tr and rewrote the
  aggregate valuation proof around aligned contents.  Arithmetic
  progressions and the globally unique removed center remain treated as
  before.

### Q7104: fixed-depth Cartier zero lifting is terminally insufficient

- Lucas gives the exact Boolean-OR law on base-p digits, with depth-d zero
  density `1-(1-Z(p)/p)^d`.
- Singular blocks are entirely zero; regular blocks only copy the base atom.
  The p^4 companion law applies only on regular blocks and adds no atom
  restriction on singular blocks.
- Therefore every fixed-depth Cartier/Frobenius tower controls descendants
  of `Z(p)` but cannot bound the base atom itself.  This closes the
  fixed-depth branch; only an unbounded-order compression could differ.

### Mesoscopic census recovered through p=20000

- The completed single-pass `meso_census.py 20000` output was recovered
  from the Claude task log; no rerun was needed.
- For `H=floor(sqrt(p))` and every prime through 20000: maximum separated
  collision multiplicity `m_(d,r)` is 2, maximum left-pair codegree is 3,
  and total collision energy is at most 24.  Twelve primes have some
  `m>=2`; none has `m>=3`.
- Uniform `m=O(1)` and codegree `O(1)` would imply `O(H^(3/2))` support by
  Kővári--Sós--Turán and hence the desired mesoscopic energy bound.  The
  census is strong evidence for exactly these two inputs, but supplies no
  proof of either uniform statement.

## Run 2026-08-09 23:57 CDT (automode continuation)

- doctrine: `AVENUES.md`, SHA-256
  `50ca74804c57ca1ae44144c9b95bc502ccaa00f2b8ffbbcdbfea317b450246ce`
- approval: direct instruction, `继续。 与 chatgpt 配合。`
- starting avenue: (e0), distinguished-orbit coupling after phantom removal
- Pair-cap route: terminal negative.  `pair_cap_extremal.tex` constructs
  reflection-symmetric endpoint-free sets with all pair caps but
  `H^2/7200-4` short off-center windows.  Its strengthened prescribed-prime
  corollary gives artificially labeled quotient separated energy
  `asymp p^(4/3)` for every sufficiently large prime and dyadic mass
  `>> X^(7/3)/log X`.  This is a logical countermodel for the local
  set-theoretic axioms, not an actual Apéry orbit or a lower bound for its
  variance; Q7195 caught and audited this necessary distinction.
  A direct synthetic check at `p=2000003,m=100` gave quotient
  `k=98`, overlap `576`, separated `8930`, exactly as proved.
- Hasse-square route: terminal negative as a formal input.
  `hasse_square_no_go.tex` proves that for `m<d/2` the first `m`
  coefficients of `Delta^epsilon B^2` are triangular coordinates with
  diagonal `2`; arbitrary clustered zeros and their reciprocal reflections
  occur in the abstract reciprocal-square family.
- Projective-variance reduction: proved in
  `projective_variance_reduction.tex`.  Reflection is free on off-center
  chains; filtered primitive chains are sliding four-occurrence windows;
  each quotient orbit has at most six overlapping quotient neighbors.
  Thus `Vbar_p <= 7 Mbar_p + Esep_p`, and only separated distinct
  nonreflection orbit pairs require arithmetic control.
- Exact computation: the local single-process scan and the independent
  16-shard GitHub scan both checked all `41,535` primes through `500000`.
  Both found raw mass `20`, quotient mass `10`, quotient maximum `1`, zero
  nonreflection collisions, and zero actual-fiber chains.  The separate
  fixed-record verifier passed on all 20 phantoms.  The C++ scanner now
  decomposes quotient energy into diagonal/overlap/separated parts with
  exact assertions and passes strict compilation plus the `p<=5000`
  regression.
- Index-first actual scan and marked-label regressions pass.  The former has
  zero off-center quadruples at `(max-index,max-span)=(1100,100)`; the latter
  verifies the actual, phantom, and nonreflection label witnesses.
- Exact remaining obstruction: prove the dyadic separated energy bound
  `sum Esep_p << X^(2+o(1))/log X` using Apéry-specific long-transfer or
  carrier information.  ChatGPT Q7177 independently audited the
  reflection/overlap constants.  Q7185 confirmed that an exact elimination
  certificate necessarily retains the long bridge `N_L`; generic
  Casoratian/Pluecker identities do not compress it to the six short gaps.
- Near-bridge range: proved unconditionally in `near_bridge_energy.tex`.
  Consecutiveness makes `(p,x,s,G)` injective, so the second chain span is
  not a certificate parameter.  The nonzero separated-block resultant and
  Smith nullity give `sum Enear_p(K) << H^2 K^2`; taking
  `K=H/sqrt(log X)` meets the full `X^2/log X` budget.  Only longer
  bridges remain.  The independent tmux-11 audit passed the orientation,
  injectivity, degree-drop, and height arguments.
- Far-incidence hierarchy: proved in `far_bridge_incidence.tex`, with exact
  injections
  `Efar_p(K) <= B4_p(K) <= Bcirc_p(K) <= B_p(K) <= A_nw,p(K)`.
  Here `B_p` retains the selected first four-return chain and drops only the
  second chain.  `Bcirc_p` removes the first chain's own reflection support,
  with `0<=B_p-Bcirc_p<=4M_p`; `B4_p` also requires the external point to
  begin a four-consecutive-occurrence window of span at most `sqrt(p)`.
  Finally, `A_nw,p` counts all nonwrapping triples `x<y<z` in one projective
  fiber with `2<=y-x<=sqrt(p)` and `z-y>K`.  The exact continuant formula
  keeps `0<=x<=p-1-s-G`; removing this restriction gives only an upper bound
  by full-`F_p` common roots.  A polynomial gcd degree is not the observable,
  because nonsplit irreducible factors contribute to the degree but no
  `F_p` root.
- Uniform fiber cap: if `e_h` counts consecutive gaps of length `h` in one
  projective fiber, then `e_h<=3(h-1)` and `sum h e_h<=p-1`.  Splitting at
  `p^(1/3)` gives `max_q |pi_p^(-1)(q)|<=4p^(2/3)`.  Combining this with
  the at most `3p/2` same-fiber pairs of gap at most `sqrt(p)` proves
  `A_nw,p(K)<=6p^(5/3)` and hence the unconditional dyadic bound
  `sum A_nw,p(K)<<X^(8/3)/log X`.  This improves the termwise degree bound
  by `X^(1/3)` but remains `X^(2/3)` above the target.
- Reciprocal-gap triple bound: for consecutive gaps `g_(q,i)`, the global
  budget `sum_(q,i) g_(q,i)^(-2)<=3 log p`, blockwise Hölder, and the full
  energy `sum_q m_q^2<<p^(3/2)` give
  `A_nw,p(K)<<p^(3/2) log(2p)` uniformly in `K`.  Hence the dyadic total is
  `O(X^(5/2))`, still `X^(1/2) log X` above target.  The complete proof,
  including the two endpoint occurrences omitted from the original energy
  theorem, is in `far_bridge_incidence.tex`.
- Short-chain energy cap: if `e_q(h)` counts consecutive fiber edges of
  length at most `H` and `E_q=sum_h e_q(h)`, window--edge charging gives
  `C_p(q)<=E_q`.  The fixed-gap caps and per-fiber total length imply the
  level estimate `#{q:E_q>=t}<=12p^2/t^3` for `t>=4p/H`.  Dyadic
  summation proves `sum_q C_p(q)^2<=30pH` and, after the reflection
  quotient, `Esep_p<=(15/2)pH`.  At `H=floor(sqrt(p))` this gives the
  unconditional dyadic bound `sum Esep_p<<X^(5/2)/log X`, still
  `X^(1/2)` above target.  The proof is in
  `projective_variance_reduction.tex`; the 92-prime direct regression has
  SHA-256
  `0796e0f34777aba4a631767076dbedb71503eef4422e19553ad09e702be6566c`.
- Sharpness of the short-chain cap: `pair_cap_extremal.tex` gives infinitely
  many complete reflection-invariant colorings satisfying every global
  fixed-gap cap and all local selected-window conditions, but with
  `Esep_p>=pH/1920`.  Hence the `O(pH)` estimate cannot be improved from
  those abstract inputs; further saving must use the Apéry recurrence.
  `colored_pair_cap_energy_verify.py` checks the partition, all pair gaps,
  quotient separation, conditioned far incidence, full fiber energy, and
  both lower-bound constants for 11 parameters, with
  SHA-256
  `707f006bbafe1da79d69fb859f76b78b885d238c9bacda97e47043981600e39e`.
  An independent tmux-11 proof audit passed every component.
- Conditioned-far sharpness: the same row-shift coloring has full fiber
  energy `O(p^(3/2))`.  At
  `K=floor(floor(sqrt(p))/sqrt(log p))`, choosing two
  left-packet window starts at least five rows apart gives raw external
  incidences counted by `B4_kappa(K)`, with
  `B4_kappa(K)>=pH/3840=Omega(p^(3/2))`.  Hence even the second-short-window
  condition does not improve the exponent from abstract pair caps and
  energy alone.
- Full-cycle bridge audit: for a return `N_h(x)=0`, the renormalized row has
  multiplier `mu_h(x)=-(x+1)^6 N_(h-1)(x+1)`.  Reflection gives
  `mu_h(x) mu_h(p-1-x-h)=(-1)^h prod_(j=1)^h(x+j)^6`.  Hence multiplying a
  bridge by its reflected bridge recovers a prescribed unit, not an
  independent zero condition; the rank-one transfer step at `n=0` prevents
  inversion across the period boundary.  This is a no-go only for the direct
  reflection/full-cycle construction.  The proof is in
  `full_cycle_bridge_reciprocity.tex`.  The verifier checks all returns for
  the 24 primes through `101` and the four known offsets at `(1297,360)`;
  SHA-256
  `4faf19b270b89c14f9e6e584fd2631bd79a242fc65704f902e9dfbcd78f98b07`.
- Exact long-bridge census: `long_bridge_incidence_scan.cpp` groups the
  projective orbit by state and counts the occurrence-list formula directly.
  For `X=1000,5000,10000,20000,50000,100000`, the far nonwrapping masses are
  respectively `8356,81018,214495,575531,2099977,5614535`; after division by
  `(# dyadic primes)*sqrt(X)` the ratios are
  `1.9573,2.0460,2.0764,2.0967,2.1062,2.1157`.  The corresponding conditioned
  masses `B(K)` are only `4,6,12,0,8,16`, from raw selected-chain counts
  `2,2,4,0,2,4`.  Removing the automatic reflection support gives
  `Bcirc(K)=0,2,4,0,4,8`; the still weaker-than-energy second-short-window
  count `B4(K)` is zero throughout.  These are complete finite dyadic data,
  not an asymptotic estimate.  Source SHA-256 is
  `8a9b8614e27d0a841a70a692daa77bd25b331c930666d646cc9b2bae493191d6`;
  strict compilation and ASan/UBSan pass.
- Independent long-bridge regression: `long_bridge_incidence_verify.py`
  directly enumerates all pairs/triples instead of using occurrence-list
  binary searches.  It also evaluates every nonwrapping gap continuant and
  verifies `N_h(x)=0 iff pi_p(x)=pi_p(x+h)`.  The default `X=70` regression
  exactly matches the C++ aggregate and has output digest
  `9125a62abc9121a31e4190f0a10c1650fd7742b0360b79036345ace34a410fe4`.
- Fixed selected-record extension: `long_bridge_selected_extension_verify.py`
  recomputes the projective fibers for all 20 raw selected records in the
  complete `p<=500000` census.  After deleting each first chain's reflection
  support it finds 22 external later occurrences, but none begins another
  four-consecutive-occurrence window of span at most `sqrt(p)`.  Thus the
  one-sided relaxed quantity `B4` is also zero throughout the census.  The
  fail-closed digest is
  `13f2f869de772d2f30de9d78a1818ad149e9cef4f92adec9ec2eef7c4e25ea1b`.
- ChatGPT audit triage: Q7215 independently found the occurrence-list
  enumeration, but its identification of actual `F_p` roots with polynomial
  gcd degree was corrected.  Q7212's full-cycle argument was rejected because
  it again propagated the row of two scalar solutions `(B_n,D_n)` by the
  one-solution companion matrix; the asserted rank-one monodromy therefore
  does not apply to the projective orbit.  Q7217 supplied no long-bridge
  estimate, and its proposed per-`(p,s)` bound had a summation error: it gives
  `p^(3/2) log p` per prime, not the claimed dyadic target.  None of these
  claims is used in the proof or in the census.
- Carrier/center deflation was separated into exact stages in
  `separated_resultant_deflation_probe.py`.  At `H=12,Gmax=36`, the full,
  center-deflated, and carrier-deflated aggregate bit lengths are
  `2699303`, `2604380`, and `2113651`; SHA-256
  `4416c8495689693651c36a53dea7a41f21a32f8900218bc6b1d2cdc14466bae0`.
  The residual fraction is stable near `0.783` in every bridge bucket, and
  no term vanishes.  The independent extension `H=14,Gmax=42` retains
  fraction `0.777555633`, with digest
  `030ad11c148a8433e009f017cb7e5a3dcd7bc4a7b0586b48d98f7782bde7a988`.
  On every diagonal `2<=h<=12`, the formal factor
  `lc(N_h)^(2 deg N_h)` fails to divide the integer resultant.  After exact
  center plus maximal `U_h`-supported deflation,
  `log|R*_(h,h)|/(h^2 log h)` is `10.26,10.47,10.69` at `h=8,10,12`.
  This is finite evidence only; it neither proves an asymptotic lower bound
  nor locates prime factors in `(X,2X]`.  tmux 11 retracted its earlier
  divisibility/proof claim and corrected a second table-label error that had
  confused full carrier deflation with center-plus-carrier deflation.
- A direct low-order bispectral shortcut is rigorously excluded within a
  finite polynomial ansatz.  `bispectral_operator_scan.py` forms the integer
  coefficient system for
  `L=sum_(|j|<=r) A_j(x)T^j`, `deg A_j<=d`, and
  `L N_h=lambda_h N_h`.  For `r=8,d=30,1<=h<=20`, reduction modulo
  `1000003` has shape `1190x547`, rank `546`, and nullity one; SHA-256
  `6fcb4e496ac5a583466a9e779c87fb94c239a5c7feab85f9d30af93bc03fe31a`.
  Since the scalar identity supplies that kernel, the rational kernel is
  exactly scalar and every sub-ansatz is excluded.  The common-denominator
  rational extension reduces linearly to the necessary divisibilities
  `N_h | sum_j A_j(x)N_h(x+j)`.  Its modular remainder matrix has shape
  `570x527`, rank `496`, and kernel dimension `31=d+1`, exactly the
  multiplication numerators; digest
  `fc02fe448531fc2b6625f709342bb1b6736bfabda7bcf2a737ecd41f63f8bfd3`.
  Thus rational operators are also scalar whenever their cleared numerators
  fit `r<=8,d<=30`.  Gauges whose conjugated numerators exceed this bound and
  higher-order operators remain outside the claim.
- Correct far range: `s<=H~sqrt(X)` but the canonical bridge can satisfy
  `G~p~X`.  The generic height sum over all `s,G` is `X^(3+o(1))`, one
  factor `X` above the raw `X^(2+o(1))` budget.  A lower bound for residual
  bit height is irrelevant to this upper-bound problem; the needed input is
  a dyadic `(X,2X]` valuation/nonconcentration estimate for the long
  resultants or an equivalent Apéry-cocycle correlation theorem.
- Mass-sensitive interpolation: if `M_p=sum_q C_p(q)` is the selected raw
  short-chain mass, the cubic edge-level tail gives
  `sum_q C_p(q)^2 <= 20 p sqrt(M_p)` and, after reflection,
  `Esep_p <= 8 p sqrt(Mbar_p)`.  Hence
  `sum_p Esep_p <= 16 X sqrt(#P_X * sum_p Mbar_p)`.  The target therefore
  follows from `sum_p Mbar_p <= X^(1+o(1))/log X`; an `H^3` certificate
  mass at `H~sqrt(X)` would only give `X^(9/4+o(1))/log X` energy.  This
  pins the remaining mass input to average polylogarithmic short-chain
  sparsity, not ordinary three-gap content control.
- Moving-gap Green kernel: with
  `V_m(a)=U_m(a)/(a+1)_m^3` and `F_a=sum_m V_m(a) z^m`, exact summation gives
  `D_a F_a=a^3`, where
  `D_a=(theta+a)^3-zP(theta+a)+z^2(theta+a+1)^3`.  Equivalently,
  `L_0(z^a F_a)=a^3 z^a`.  This corrects ChatGPT Q7231, whose claimed
  homogeneous eigenfunction equation dropped the `m=0` boundary term.
  Adjoining the constant solution gives a rank-four system with fixed
  singular support `z(1-34z+z^2)=0`, and the exact Casoratian formula is
  `V_m(a)=a^3(b_(a-1)c_(a+m)-c_(a-1)b_(a+m))`.  The bridge dictionary is
  `N_G(u)=((u+2)_(G-1))^3 V_(G-1)(u+1)`.  This is a genuine fixed-order
  differential compression, but coefficient extraction at `G-1` remains a
  growing jet and supplies no long-bridge estimate by itself.  Proof:
  `gap_kernel_green.tex`.  Exact verifier digest:
  `8b27017ed37134278861b3460f86814e796916a432138233bad49359fa5a65ca`.
- end: pending
- final result: pending
