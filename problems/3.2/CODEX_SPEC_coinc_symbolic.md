# CODEX SPEC — difference-curve symbolic factorization (campaign 3, avenue a2)

## Mission
Determine the COMPONENT STRUCTURE of the difference curves X_{h,h'} over Q.
Write report to `CODEX_COINC_SYMBOLIC_report.md` in this directory (overwrite ok).
Background: read `campaign3_questions/CTX.txt` (authoritative, contains (D1)-(D6)).
Sympy is available (import sympy). This is a HARD symbolic task: be patient,
use resultants/gcds cleverly, never brute-force degree-60 factorizations blindly.

## Objects
N_h(r): N_1=1, N_2 = P(r+1), N_{h+1}(r) = P(r+h) N_h(r) - (r+h)^6 N_{h-1}(r),
P(u) = 34u^3+51u^2+27u+5. deg N_h = 3(h-1).
F_{h,h'}(r,r') = N_h(r) * prod_{j=1}^{h'}(r'+j)^3 - N_{h'}(r') * prod_{j=1}^{h}(r+j)^3.
X_{h,h'} = zero locus of F_{h,h'} in (r,r').

## Tasks

### S1 — same-gap curves F_{h,h}
For h = 2..8: factor F_{h,h}(r,r') over Q (sympy factor; it is antisymmetric-ish in
(r,r') so (r-r') should divide — verify). List ALL irreducible factors with degrees.
Detect the mirror component: is there a factor of the form (r + r' + c_h) (linear)?
If yes give c_h as a function of h; if a nonlinear extra factor appears, display it.
Empirical component count c(h,h) = number of irreducible factors (excluding constants).

### S2 — distinct-gap curves F_{h,h'}
For all 2 <= h < h' <= 8 (and the h=1 column: N_1 = 1 so F_{1,h'} =
prod(r+j)^3... wait h=1 means Delta_{r,1} = 1/(r+1)^3, never 0; still include:
F_{1,h'}(r,r') = prod_{j=1}^{h'}(r'+j)^3 - N_{h'}(r') (r+1)^3): factor over Q.
Report: irreducible? If reducible, list factors + degrees. Look for patterns:
parity-induced factors from N_{2k} = (2r+2k+1) M_{2k}: does the linear factor of
N_h interact with the linear factor of N_{h'} to produce a rational component
(2r+h+1)*something - (2r'+h'+1)*something?

### S3 — genus/singularity budget
For the pairs in S2 that are irreducible: compute the degree and (if cheap via
sympy) the singular points count of the projective closure — we need the constant
in #X(F_p) = p + O(g sqrt p): report deg and a crude genus upper bound
(deg-1)(deg-2)/2 per pair. No need for exact genus.

### S4 — update-identity leverage (theory task, no compute)
Using G_{m+g}(s) = G_g(s+m) G_m(s) (transfer product), derive the algebraic
identity expressing F_{h,h'} in terms of lower-index N's when h' = h + d:
substitute N_{h+d} via the composition (write N_{h+d}(r) in terms of
N_{h}(r), N_{d}(r+h), N_{d-1}(r+h+1)-type products — derive the EXACT identity
from the matrix product, verify symbolically for h=3,d=2). State whether this
exhibits X_{h,h+d} as related to X-objects at gap d — a self-similarity that
could power induction on h. Write the cleanest identity you can prove.

### S5 — verdict
One paragraph: for generic (h,h') is X absolutely irreducible (evidence)?
Complete list of discovered systematic components with formulas, and what each
contributes to N_coinc main term (p points per component per pair).

## Rules
- sympy over QQ; keep h,h' <= 8 unless trivial to go higher.
- If factor() stalls > a few min on one pair, use modular factorization mod a few
  primes ~ 10^4 to DETECT reducibility pattern instead, and mark as "mod-p evidence".
- No early stop; all five sections in the report. ASCII only.
