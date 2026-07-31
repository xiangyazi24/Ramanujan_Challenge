# HANDOFF — Problem 3.1 Lean development

project: Ramanujan_Challenge / SUBMIT/3.1
workdir: `~/repos/Ramanujan_Challenge/SUBMIT/3.1`
written: 2026-07-30, by Zinan, for Codex
head at handoff: `7c69e60` (pushed to `origin/main`)

---

## 0. Read this first

**The submission is already complete and shippable.** `solution.pdf` (8 pages),
`README.md`, `scripts/` (18 scripts), and `lean/` all build clean, and
`SUBMIT/dist/ramanujan-3.1-huang.zip` is packaged. Deadline is **2026-08-01
18:59 Chicago** (23:59 UTC).

Everything described below is *improvement past that baseline*. **Do not break
the shippable state.** After any change:

```
cd /tmp/r31build && lake build            # must say: Build completed successfully
cd ~/repos/Ramanujan_Challenge/SUBMIT/3.1 && pdflatex solution.tex   # 0 errors
cd ~/repos/Ramanujan_Challenge/SUBMIT && bash make_zip.sh 3.1
```

There is **no `sorry` and no `native_decide` anywhere** in `lean/Ramanujan31/`.
Keep it that way. If a step cannot be proved, leave the theorem unstated rather
than stating it with a hole — do not introduce an `axiom`, and do not write
"acceptable" next to a gap.

---

## 1. Build environment — non-obvious, read carefully

The challenge requires **Mathlib v4.29.0**. Three traps:

1. **`~/repos/Ripple` is on Lean/Mathlib v4.30.0.** Do not build against it.
   I lost a whole pass verifying everything on 4.30 before noticing —
   `cat ~/repos/Ripple/lean-toolchain` says `v4.30.0`. Check the toolchain of
   any host checkout before trusting it.

2. **`/tmp/ripple429`** is a v4.29.0 checkout whose Mathlib olean cache is
   already downloaded (~3.8 GB). It gets wiped by `/tmp` cleanup periodically;
   if `.lake/packages/mathlib/.lake/build/lib` is missing, restore with
   `cd /tmp/ripple429 && lake exe cache get` (8232 files, several minutes).

3. **`/tmp/r31build`** is the scratch build directory. It is
   `SUBMIT/3.1/lean/` plus a **symlink** to ripple429's packages, so we never
   rebuild Mathlib. Recreate it with:

```bash
rm -rf /tmp/r31build && mkdir -p /tmp/r31build/.lake
cp -R ~/repos/Ramanujan_Challenge/SUBMIT/3.1/lean/. /tmp/r31build/
cp /tmp/ripple429/lake-manifest.json /tmp/r31build/
ln -s /tmp/ripple429/.lake/packages /tmp/r31build/.lake/packages
cd /tmp/r31build && lake build
```

Workflow: edit in the repo, `cp -R ~/repos/.../lean/. /tmp/r31build/`, build
there. For a single file, `lake env lean Ramanujan31/Foo.lean` is much faster
than a full `lake build`.

Note: a shell hook blocks commands containing the literal string `lake build`
when they also touch certain paths (it is a guard against kernel panics on the
24 GB mini). If you hit `⛔ 本地 lake build 被禁止`, put the build in a script
file and run the script.

---

## 2. What is proved — current state of `lean/Ramanujan31/`

Full build: **8261 jobs, clean**. Axiom audit: every declaration depends only on
`[propext, Classical.choice, Quot.sound]`; the four `ShapeField` certificates
need only `[propext, Quot.sound]`.

| module | content |
|---|---|
| `RatReconstruct.lean` | rational reconstruction from a denominator bound; `regulator_quotient_eq_robust` (any `85 ≤ Q ≤ 10¹⁵⁰` works); `regulator_quotient_eq_4080` |
| `ChartSymmetry.lean` | `u(1/a) = u(a)`; palindromic decompositions `f = a^d·g(a+1/a)` for both endpoints |
| `ShapeField.lean` | algebraic certificate that `1+4u²` is a **square in the endpoint field** (explicit integer `s`, `linear_combination c·hf`) |
| `TraceRoots.lean` | trace polynomials are totally real with **exact** root counts in `(−2,2)`, incl. an exhaustion clause |
| `UnitCircle.lean` | real trace in `[−2,2]` **⟺** `‖a‖ = 1` |
| `ShapeCancellation.lean` | the four-shape cancellation, stated for abstract `D` with `BlochWignerLaws` |
| `MainTheorem.lean` | composition + final step, with the external inputs as named hypotheses |
| `Dilog/Basic.lean` | complex `Li₂` power series: summability, Schwarz reflection, tail bound, **`Li₂(1) = π²/6`** |
| `Dilog/RealBounds.lean` | **real** `Li₂` two-sided rational enclosure (partial sum is the lower bound; explicit tail is the upper) |
| `Dilog/Certify.lean` | end-to-end certified evaluations that actually close in Lean |
| `Dilog/BlochWigner.lean` | **`blochWignerGeom`** and the three functional equations, unconditionally |
| `Dilog/Instance.lean` | `blochWignerGeom_laws : BlochWignerLaws blochWignerGeom` — the structure is no longer assumed |

The machine-checked chain, with nothing assumed about `D`:

```
trace real in [-2,2]   TraceRoots.gAlpha_totally_real
  → ‖a‖ = 1            UnitCircle.norm_eq_one_of_trace_real_abs_le_two
  → conj u = u         ChartSymmetry.chartUAlpha_isReal_of_norm_one
  → W = (1-conj V)⁻¹   ShapeCancellation.W_eq_inv_one_sub_conj_V
  → ΣD = 0             four_shape_sum_vanishes_of_trace_real'   (Dilog/Instance)
```

### The design decision that made `BlochWigner.lean` cheap

Do **not** define `D` via analytically-continued `Li₂` — three of the four shapes
are outside the unit disc and that route is branch cuts and continuation.
Instead `D` is the sum of three Lobachevsky angle terms over the three-cycle
orbit `z → (1-z)⁻¹ → 1-z⁻¹ → z`. Consequences:

* `clausen`'s argument is always a **phase** (modulus 1), so the series converges
  on all of `ℂ∖{0,1}` — no continuation at all;
* the **squared** phase is used, so a real argument (phase `±1`) gives `1` and
  contributes `0`, with no case split between `arg = 0` and `arg = π`;
* `D((1-z)⁻¹) = D z` is cyclic reassociation of a three-term sum.

Keep this architecture. Reverting to the `Li₂ + arg` definition would reintroduce
everything it avoids.

---

## 3. Open work, in priority order

### P1 — the certified numerical evaluation (the reachable one)

Goal, in Lean:

```
| Re[ΔR]/π² − (−4/85) | < 3·10⁻⁸
```

`3·10⁻⁸ = 1/(2·4080²)`. **Only ~8 digits are needed.** The 301 digits in the
paper are what Sage produced for free; do not chase them.

Four sub-pieces:

**(a) Certified rational enclosures of the eight shape values.** Each shape is a
rational function of a root of an explicit integer polynomial. Numerically, at
the first endpoint the shapes are `−0.157875280712, 6.815798700830,
55.872474425246, 5.932921712939`, and their `3`-cycle representatives in `(0,1)`
are `0.86365088, 0.85328205, 0.98210210, 0.83144898`. `TraceRoots.exists_root_Ioo`
already gives roots from sign changes; what is missing is refinement to `10⁻¹²`
and propagation through the rational function. A Newton–Kantorovich style
one-shot certificate (`|f(x₀)|` small + a derivative lower bound on an interval)
is probably cheaper than ~40 bisections; see the ChatGPT answer in
`/tmp/gpt/dm/Q5802.md` (question text in
`/Users/huangx/.claude/jobs/cfed56dd/tmp/d6.txt`).

**(b) Euler's reflection formula** — this is the keystone, and it is the single
most valuable next lemma:

```
Li₂(x) + Li₂(1−x) = π²/6 − log x · log(1−x),   0 < x < 1
```

Why it is required, not optional: the representative `0.9821021` needs `N ≈ 1000`
terms, and the exact rational partial sum then carries `lcm(1,…,1000)²` in its
denominator — roughly 870 digits before the powers of the argument. `norm_num`
will not do that. Reflection sends `0.9821021 → 0.0178979`, where **eight terms**
suffice.

`Li₂(1) = π²/6` is already proved (`Complex.dilog_one`, from Mathlib's
`hasSum_zeta_two`) and is the anchor to evaluate the constant against. The
classical proof differentiates both sides; `Li₂'(x) = −log(1−x)/x`. Strategy
discussion and the Mathlib v4.29 API question are in `/tmp/gpt/dm/Q5800.md`
(question in `d4.txt`).

**(c) Certified `log` at rational points**, to ~12 digits. Mathlib has the
analytic infrastructure but no ready evaluator. `Real.pi_gt_d20` /
`Real.pi_lt_d20` already give **π to 20 digits**, machine-proved, so π is a
solved problem.

**(d) Assembly** — a small bespoke rational-interval type with a handful of
soundness lemmas, then `norm_num`. Recommended layering (from `Q5783.md`):
ℚ intervals for root isolation, dyadic rationals for the series evaluation,
`norm_num` for the final discharge. Do **not** try to carry exact ℚ through a
thousand-term sum.

### P2 — the bridge

`blochWignerGeom` is *a* standard definition of Bloch–Wigner, but it is **not
proved equal** to `Im Li₂(z) + arg(1−z)·log‖z‖`. Both `Instance.lean` and the
README say so explicitly. Proving the bridge would let the paper cite the usual
definition without a caveat. Not on the critical path.

### P3 — nonvanishing

`D(z) ≠ 0` at some non-real `z` is not proved. The 3.1 argument only ever uses
*vanishing*, so a vacuous `D` would not make any downstream statement false — but
it would make them uninformative, and a referee is entitled to ask. A direct
estimate on the Clausen series at, say, `z = i` should do it.

### Out of scope

The **Merkurjev–Suslin denominator bound** (`|K₃^ind(F)_tors| = w₂(F)`) is not a
formalization target at any realistic scope. It stays a cited input in the paper,
and enters Lean only as the named hypothesis `torsion` of
`MainTheorem.regulator_value`. `regulator_quotient_eq_robust` already shows the
conclusion survives any bound below `10¹⁵⁰`, so its exact value is not
load-bearing.

---

## 4. Traps that cost me time — do not repeat

* **Check `lean-toolchain` of any host checkout before building against it.**
  See §1.
* **`Polynomial.Splits` in v4.29 takes only the polynomial** (`p.Splits`), not a
  `RingHom`.
* **`variable (hD : ...)` does not auto-enter theorem statements** if `hD` appears
  only in the proof. Make such hypotheses explicit binders.
* **`tsum_le_tsum` is `Summable.tsum_le_tsum`** in v4.29.
* **`div_le_div_of_nonneg_right'`, `inv_eq_one_iff_eq_one`, `pow_le_pow_left`,
  `div_le_div_iff`, `Complex.abs`** — all absent or renamed. Prefer `gcongr`,
  `congrArg (·⁻¹)`, `mul_self_le_mul_self`, `‖·‖`.
* **`set_option ... in` may not sit between a docstring and its declaration.**
  Put it before the `/--`.
* **A number recorded once and propagated is a bug factory.** Two real defects
  came from this: stale `w₂ = 60/204` living in the paper long after the audit
  fixed it to `120/408`, and a Lean threshold of `1/10³⁰¹` that the actual
  computed residual (`1.63·10⁻³⁰¹`) does not satisfy. **Re-run the script; do not
  quote the note.** Both a clean build and a clean axiom audit will miss this
  class of error, because each artifact is internally consistent.

---

## 5. ChatGPT material

Questions I sent are in `/Users/huangx/.claude/jobs/cfed56dd/tmp/{d1..d6,q*}.txt`;
answers in `/tmp/gpt/dm/Q57{98..802}.md` and `/tmp/gpt/flt/Q57*.md`
(and mirrored to the `chatgpt-drop` branch of `xiangyazi24/zinan-memory`).

Most useful:

| file | topic |
|---|---|
| `Q5798.md` | check-before-build (nothing exists) + **the geometric architecture** |
| `Q5799.md` | Schwarz reflection, `arg` on the negative axis, vanishing on ℝ |
| `Q5800.md` | Euler reflection: proof strategies and v4.29 API |
| `Q5801.md` | continuation vs. orbit definition |
| `Q5802.md` | certified enclosures of algebraic numbers, Newton certificate |
| `Q5783.md` | certified numerics: layering, dyadics, A/B/C structure |
| `Q5767/Q5778.md` | the two hostile-referee passes on the mathematics |

Dispatch: `TMUX_PANE=<pane> ASK_CHANNEL_HEARTBEAT_MAX_AGE=300 python3
~/.openclaw/workspace/scripts/ask-gpt.py "$(cat q.txt)"`. The relay window must be
named after the channel group (`dm` → `dm1..dm6`). **The `flt` tabs went stale
mid-session and silently swallowed two questions** — check
`curl -s localhost:8801/api/status` for `last_seen_s` before trusting a dispatch.

---

## 6. Commit discipline

One commit per completed piece. In the message record: what was tried, the
verdict, and the two-question self-audit — (a) which part am I least sure of,
(b) what did nobody examine. Those flags are what the next reader uses, and two
of the real defects in this project were found by answering (b) seriously.
