# HANDOFF — Ramanujan Challenge submission prep

```
project: Ramanujan_Challenge
workdir: ~/repos/Ramanujan_Challenge
updated:  2026-08-01 America/Chicago
deadline: 2026-08-01 23:59 UTC  ( = 2026-08-01 18:59 America/Chicago )
```

If your current project is not `Ramanujan_Challenge`, **stop** and hand this to
the right window — do not execute it here.

A separate, narrower handoff for the 3.1 Lean development lives at
`SUBMIT/3.1/HANDOFF_LEAN.md`. This file covers the submission as a whole.

---

## 0. The one thing to internalise before touching anything

**A formal statement counts only if it mentions the challenge's own recurrence,
coefficients, or sequences.**

Until today this repo reported "0 sorry ✅" for six problems on the strength of
theorems like

```lean
theorem problem25_identity : ∃ (p q : ℕ → ℝ), Tendsto (fun n => p n / q n) atTop (𝓝 L) :=
  ⟨fun _ => L, fun _ => 1, by simp⟩
```

That is witnessed by constant sequences and says nothing. P2.2's version
asserted the challenge's initial values and then converged to `179/306`, not γ.
There was also `sign_flip_P … : True := trivial` with four unused hypotheses,
and `u1_value : x = x := rfl`. All of these are now **deleted** (commit
`109a033`). P2.2 and P2.3 have since been rebuilt as unconditional proofs.

Do not reintroduce this pattern in any form. A `sorry` is fine and legible; a
vacuous existential shipped in a submission reads as a faked formalization. If a
statement cannot be proved, either state it as an explicit **hypothesis of the
theorem** (the pattern used in 2.1, 2.3 and 2.8 — see §2) or leave a `sorry` and
say so in the README.

Same rule for prose: `STATUS.md` used to mark 2.1–2.8 and 3.1 all "unconditional".
Reading the `.tex` showed 2.2's decisive step was written as *"can be checked by
comparing sufficiently many evaluations"* (it was not done) and 2.6 says of
itself *"verified numerically to 39 digits"*. `STATUS.md` now separates
**audited this session** from **inherited claim**. Keep that distinction; do not
promote a row without doing the work.

---

## 1. Current state: what is shippable right now

Packaged submission set = **2.1, 2.3, 2.8, 3.1**, under `SUBMIT/`.
P2.2 is now unconditionally proved in the main Lean project and has a rewritten
standalone PDF, but has not yet been copied into a `SUBMIT/2.2` package.

```bash
cd ~/repos/Ramanujan_Challenge/SUBMIT
bash make_zip.sh 2.1 2.3 2.8 3.1     # -> dist/ramanujan-huang.zip, 2.0 MB
```

The script refuses to build an archive for a problem with no `solution.pdf` or
`solution.tex`, and fails loudly over 50 MB (the challenge's limit).

Each package has: `solution.pdf` + `solution.tex`, a self-contained `lean/`,
`verify.py`, and a `README.md` stating the axiom audit, the build scope, and
exactly which classical result is cited.

| # | topic | cited input (the only one) | Lean |
|---|---|---|---|
| 2.1 | PCF → 6/(3−π) | Cohen Entry 5.3.22 (arXiv:2607.06581) | 0 sorry, 3 std axioms |
| 2.3 | π+e | Lambert CF value π/4 | 0 sorry, 3 std axioms |
| 2.8 | √10005/π | Chudnovsky CM period-derivative evaluation | 0 sorry, +5 `native_decide` on one theorem |
| 3.1 | knot π² | see `SUBMIT/3.1/HANDOFF_LEAN.md` | 0 sorry |

**Regression gate — run these before and after any change:**

```bash
# main Lean project (v4.30): 3427 jobs, 0 errors, warnings = exactly 5 real sorry
cd ~/repos/Ramanujan_Challenge/lean && ~/.elan/bin/lake build

# the submitted standalone packages
cd ~/repos/Ramanujan_Challenge/SUBMIT/2.1/lean && ~/.elan/bin/lake build && ~/.elan/bin/lake env lean AxiomCheck.lean
cd ~/repos/Ramanujan_Challenge/SUBMIT/2.3/lean && ~/.elan/bin/lake build && ~/.elan/bin/lake env lean AxiomCheck.lean

# numerics
python3 ~/repos/Ramanujan_Challenge/SUBMIT/2.1/verify.py
python3 ~/repos/Ramanujan_Challenge/SUBMIT/2.3/verify.py
```

The 5 legitimate `sorry`s in the main project are 2.4 (1), 2.8 (1),
3.2/Main (3). They are declared-hard statements, not evasions.

### Build traps

- **Toolchain pins differ across packages.** 2.1 and 2.3 are `v4.30.0`; 2.8 and
  3.1 are `v4.29.0`. Building 2.8 against v4.30 fails —`ModularForm.delta` was
  deprecated in favour of `ModularForm.discriminant` (2026-03-23) and the
  level-one group presentation changed. Do not "helpfully" bump a pin.
- `SUBMIT/*/lean/.lake/packages` is a **symlink** to the main project's package
  tree, so the standalone packages build without re-downloading Mathlib. That is
  why the READMEs say the build was verified as a standalone Lake project *but
  not* from a cold cache on a fresh machine. Keep that wording accurate. The
  symlink and build products are gitignored (`SUBMIT/*/lean/.lake/`); if you
  recreate a package, recreate the symlink:
  ```bash
  ln -s ~/repos/Ramanujan_Challenge/lean/.lake/packages SUBMIT/<p>/lean/.lake/packages
  ```
- `make_zip.sh` excludes `.lake`, `*.aux/log/out`, `__pycache__`, `*.olean`.

---

## 2. The pattern that closed 2.1 and 2.3 — reuse it

Both were closed the same way, and it is the template for the rest:

1. **Find the structure**, not an asymptotic story.
   - 2.3: order 4 = 2×2. The operator annihilates **every** product
     `X_{n+2}Y_{n+3}` of a Lambert-recurrence solution with a
     derangement-recurrence solution. `m!` satisfies the *same* recurrence as
     `D_m`, which is why π and e appear together and additively. Result: the
     ratio splits **exactly**,
     `p_n/q_n = 4·B_{n+2}/A_{n+2} + (n+3)!/D_{n+3}` — no Poincaré–Perron.
   - 2.1: the challenge PCF is the sign-flip of the tail of a published CF, and
     the sign-flip lemma is proved **at the level of convergents**
     (`P̃_n = (−1)^{n+1}P_n`, `Q̃_n = (−1)^n Q_n`), so no tail-convergence
     question ever arises.
2. **Prove the identity generically**, with free initial values. `tensor_rec`
   takes arbitrary `X, Y`; `cf_neg_convergent` takes arbitrary `c, d`. An
   identity that holds only for the specific sequences is much weaker evidence
   and usually means you have not found the structure.
3. **Isolate the one classical input as an explicit hypothesis** of the Lean
   theorem, so the dependency is in the statement:
   ```lean
   theorem problem23_pi_add_e
       (hLambert : Tendsto (fun m => (lambertB m : ℝ) / (lambertA m : ℝ))
                     atTop (𝓝 (Real.pi / 4))) : …
   ```
4. **Do not trust the citation — fetch it.** For 2.1 I pulled
   arXiv:2607.06581 and confirmed Entry 5.3.22 verbatim, including its
   machine-readable form and displayed quotients 42, 396, 1047, 38400, 4340.
   This matters: see §3, where the inherited citation for 2.2 turned out to be
   simply wrong.
5. **`verify.py` must forward-solve the challenge's recurrence from the problem
   statement's own initial values** and compare against the closed form. That is
   the transcription/faithfulness gate, and it is what caught nothing this time
   only because it was run.

---

## 3. Per-problem state and next actions

### 2.2 (γ) — DONE unconditionally; package next

The old direct Aptekarev identification was false. The replacement proof uses
an exact first-order Ore transform of Rivoal's recurrence,

```text
C_n = (T_{n+1} + (n+1)(3n+4)T_n)/(8n+11),
```

whose transformed initial triples are exactly the challenge's. The fixed
Rivoal solutions have WZ-verified finite sums with positive weights
`(2n+k+1) choose(n,k)^2/k!`.

For the normalized weights, the adjacent ratio yields a finite birth–death
Stein identity. With `G=k³-(n-k)²`, it proves
`E[G²] ≤ 81 n³√n`. A good/bad saddle decomposition then gives the explicit
harmonic error bound

```text
5/M + 8δ + 81(6+5 log n)/(δ²√n),
```

which tends to zero. The rational correction is `≤1/(2n+1)`, and the Ore
quotient is a positive weighted average of adjacent Rivoal quotients.

Lean modules: `Problem22.lean`, `Problem22Concentration.lean`,
`Problem22Moment.lean`, `Problem22Harmonic.lean`. Public theorem:
`problem22_solved : Problem22Claim`; 0 `sorry`; standard three-axiom audit in
`lean/P22Audit.lean`. The matching proof is `problems/2.2/proof.tex` and builds
to a five-page PDF. Remaining submission work is packaging only.

### 2.6 (ζ(2)+ζ(3) series) — real gap, honestly self-declared

`problems/2.6/proof.tex` says of itself "verified numerically to 39 digits" and
its "complete algebraic proof" ends at "reduces to a sum that decomposes as
V + W … reducible via Abel summation, the substitution z = y/(1+y)², and
cyclotomic distribution identities at ω = e^{2πi/3}". That last step is not done.

The solid part *is* solid and is worth building on: the recurrence completely
factors into two first-order Ore factors over ℚ(n), the recessive solution is
hypergeometric with `v_{n+1}/v_n = (n+3)²/(2(n+4)(2n+7))`, normalised
`v_n = 1120((n+2)!)²/(2n+6)!`, and reduction of order gives an explicit closed
form for `u_n`. What is missing is the evaluation `Σ u_n = ζ(2)+ζ(3) − 2077/720`.
Note the write-up's own claim that the evaluation has **total weight three, no
Li₄(1/2)** — that is a useful constraint on the answer.

### 2.4, 2.5, 2.7 — not re-audited this session

Treat the ✅ in any older doc as unverified. Before doing new work on one of
these, audit it the way 2.2 was audited:

1. Read `problems/<p>/proof.tex` looking specifically for *"can be shown"*,
   *"can be verified"*, *"verified numerically"*, *"it follows that"* at the load-bearing step.
2. Forward-solve the challenge's own statement numerically and confirm the claim
   is even true.
3. Check every citation against the actual source.

Known from a quick pass: **2.5** (10 pp) reduces, by its own text, to "the single
scalar identity" `(p₀ − G·q₀)·w₊(0) = 0`, with a lot of "verified numerically"
around it. **2.7** (9 pp) claims unconditional via a rational gauge transfer
`R(n) ∈ GL₃(ℚ(n))` from Zudilin's recurrence, verified symbolically in Sage —
that one may well be genuine, and is the best candidate of the three; the thing
to check is whether the Birkhoff/subdominance step is proved or asserted.
**2.4** hinges on an "algebraic identity that can be verified by symbolic
computation" — i.e. asserted; if the identity is true, a Zeilberger certificate
would make it real.

### 3.2 — see `problems/3.2/`, unchanged this session (partial by design: density-1 unconditional, all-n conditional).

---

## 4. Environment notes

- Lean: `~/.elan/bin/lake`, main project `~/repos/Ramanujan_Challenge/lean`,
  Mathlib already built locally (`.lake` ≈ 7.4 GB), rev `c5ea00351c28` (v4.30.0).
  Single-file check: `~/.elan/bin/lake env lean RamanujanChallenge/ProblemNN.lean`.
  Acceptance is `lake build` + `#print axioms`, not `lake env lean` alone.
- LaTeX: `PATH="/Library/TeX/texbin:$PATH" pdflatex`. Run twice for refs; then
  grep the `.log` for `^!`, `Undefined`, `LaTeX Warning: Reference`,
  `LaTeX Warning: Citation`.
- **ChatGPT bridge was unavailable** in this window all session
  (`[BRIDGE: no channels for window 'ccdex']`). Also note `ask-gpt.py` takes the
  question as **argv**, not stdin — `ask-gpt.py "$(cat q.txt)"`. If the bridge
  is up in your window, the Lambert-remainder question in §5 is a good one to
  send.
- `pdftotext -layout` works well on the arXiv PDFs and is how Cohen's Entry
  5.3.22 and Aptekarev's recurrence were extracted. `WebFetch` on an arXiv PDF
  returns binary; download with `curl -sL -o x.pdf https://arxiv.org/pdf/XXXX`
  then `pdftotext`.

---

## 5. Open improvements to the two finished problems (optional, not blocking)

Both are already shippable; these would only strengthen them.

- **2.3**: the cited input is the *value* π/4 of Lambert's CF. We prove
  unconditionally that the limit **exists** (Casorati
  `A_mB_{m−1}−A_{m−1}B_m = (−1)^{m+1}(m!)²` ⟹ alternating series with term ratio
  ≤ 1/5). What would remove the citation entirely is an explicit integral
  representation for the remainder `R_m = A_m·π/4 − B_m`, which satisfies the
  same recurrence with `R_{−1} = π/4`, `R_0 = π/4 − 1`; numerically
  `R_m ≈ 1.28·m!·(1−√2)^m`. I tried four candidate kernels (see
  §"remainder hunt" reasoning; none matched: `∫₀¹ xᵐ(1−x)ᵐ/(1+x²)^{m+1}`,
  `(m!)²` times it, `∫₀^{π/4} tan^{2m+2}`, and `m!` times that). The `x = tanθ`
  substitution is where π/4 enters and is probably the right frame:
  `R_0 = −∫₀^{π/4} tan²θ dθ`. Worth one focused attempt.
- **2.1**: nothing missing beyond Cohen's entry itself. Formalizing Cohen's
  derivation (arcsin seed → Bauer–Muir → even contraction) would be a large
  independent project; not worth it before the deadline.

---

## 6. Commit discipline

- Commit per milestone with a message that says what changed *and what was
  wrong before*; the git history is the audit trail for exactly the kind of
  overclaiming described in §0.
- `git diff --check` before committing (trailing whitespace).
- `.gitignore` now covers `SUBMIT/*/lean/.lake/`; if you add a package, confirm
  build products are not staged — `git status --short` should not show
  `.lake/build/...`.
- **Nothing has been pushed.** Four commits sit locally on `main`:
  `437f258` (2.3), `109a033` (purge + status), `f1f45b5` (2.1), `b7280a8`
  (STATUS/RUN_LOG). Xiang decides when to push.
- Append to `RUN_LOG.md` at the end of a working session.
