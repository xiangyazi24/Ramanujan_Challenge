# UNDERSTANDING — Ramanujan Challenge

The one maintained doc for this repo. Update it in place; do not spawn parallel
`DOCTRINE_*` / `HANDOFF_*` / dated-snapshot files (they rot and then mislead).

## SUBMIT/ is the deliverable, not a workspace

**`SUBMIT/<problem>/` holds the cleanest, self-contained submission package for
that problem. Never work in it.**

- Work happens in `lean/`, `problems/<n>/`, `notes/`, `working_notes/`.
- A problem's package is *copied into* `SUBMIT/` only once it is clean:
  builds green on the pinned toolchain, no `sorry`, `#print axioms` showing only
  `propext / Classical.choice / Quot.sound`, README stating exactly what is and
  is not proved.
- Nothing half-finished, no scratch files, no experiment leftovers under
  `SUBMIT/`. If you find yourself editing a file under `SUBMIT/` to make it
  compile, you are in the wrong directory — fix it in the working tree and
  re-copy.

## Toolchain is pinned to Lean v4.29.0

The competition requires it, and `SUBMIT/2.8/README.md` and
`SUBMIT/3.1/README.md` both pin `leanprover/lean4:v4.29.0`.

`lean/lean-toolchain` is the source of truth. Two hazards seen in practice:

- An agent working in a **clone on a different toolchain** (a v4.30.0 sandbox
  copy) verifies against the wrong Mathlib; its "0 errors" does not transfer.
  Check `lean-toolchain` before trusting any verification result.
- **Never rsync `.olean` files between clones on different revisions.** Lean
  version-tags oleans, so the failure is loud (missing-module / unknown-prefix
  errors) rather than silent — but it wrecks the build tree. Recovery is
  `lake exe cache get`, which refetches the oleans for the pinned revision.

## Repo layout and who writes where

Several agents work this repo concurrently. Content-wise they stay out of each
other's way; mechanically they do not, so:

- **Commit with path scoping** (`git commit -- <paths>`), never `git add -A`.
  A full-tree commit sweeps up another agent's uncommitted work under your
  message. This has already happened.
- Expect your index to be cleared by a concurrent commit. Re-stage and retry
  rather than assuming your edit was lost.
- A sandboxed agent that can only write under its own clone has to be merged in
  explicitly. Wire its clone as a git remote and cherry-pick by path.

## Verification gate

A result counts as done only when all of these hold:

1. `lake build <Module>` green on v4.29.0 — not merely `lake env lean` on one file.
2. `#print axioms <theorem>` shows only `propext, Classical.choice, Quot.sound`.
   **A green build with a right-looking statement is not a proof**: this repo has
   produced theorems that compiled clean, matched their target value to 25 digits
   numerically, and still depended on `sorryAx` through three layers.
3. Carried hypotheses read and judged satisfiable. `#print axioms` is silent on
   them, so a conditional theorem can be axiom-clean and vacuous.

`scripts/truth.sh <leanDir> [modulePrefix]` reports 1 and 2 in plain language;
`scripts/bank.sh` is the gate it runs; `scripts/deps.sh` + `scripts/impact.sh`
give the real dependency graph ("I changed X, what must be re-verified").

## Before proving anything: grep the repo

This repo is large and several problems share machinery. Two false gaps were
found in a single session — a continuity obligation that followed in five lines
from a lemma proved fifteen minutes earlier, and a Landen identity that already
existed in `Problem26WeightThree` while a hundred lines were spent re-deriving
it. Assume more exist. Grep by statement shape, not just by name.

## `lake env lean` is not a gate on a large file — use `lake build`

On `Problem24QuadraticAlt.lean` this bit four separate times in one session:
`lake env lean <file>` reported 0 errors while `lake build <Module>` reported
real ones — forward references (a lemma used ~500 lines above its definition),
an `Unknown identifier` after a block was moved, and `open`-scope differences
when a block written in a scratch file was pasted in. It appears to be reusing a
stale `.olean` rather than elaborating from scratch.

Use `lake env lean` for fast iteration on a *scratch* file only. The gate is

```
/Users/huangx/.elan/toolchains/leanprover--lean4---v4.29.0/bin/lake build \
  RamanujanChallenge.<Module>
```

Corollary: a block that compiled in a scratch file is **not** verified for the
main file. Scratch files usually have different `open`s (the main file opens only
`Filter Set Topology`), and their declaration order differs.

## Endpoint singularities: reach for the tendsto-FTC, not continuity

`intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le` wants `ContinuousOn f` on
the CLOSED interval. When the antiderivative only has *limits* at the endpoints —
which is the usual situation here, since Lean's junk conventions (`log 0 = 0`,
`0/0 = 0`) tend to make the endpoint value equal the limit by accident — use

```
intervalIntegral.integral_eq_sub_of_hasDerivAt_of_tendsto
```

instead. It takes the two one-sided limits in place of closed-interval
continuity. Two `sorry`s in this file were deleted outright this way rather than
discharged: proving the limit needs only first-order information, while proving
continuity at the endpoint needs second-order.

For integrability with an endpoint singularity, the corresponding move is
`intervalIntegrable_iff_integrableOn_Ioo_of_le` (drops both endpoints) plus a
*local* majorant. A single global majorant usually does not exist: near `t = 1`,
`W0 = o(1-t)` does not bound `W0 · H1/(1-t)` by `C(1 + log² t)`, but it does give
`|W0| ≤ 1-t` and hence domination by `H1`, which is integrable.

## Numerical verification: substitute before trusting PSLQ

Direct quadrature of `∫₀¹ W0(t) H(t)/t dt` stalls near 14 digits at the `log²`
endpoint. At that precision PSLQ reports "no integer relation" for `I10` — which
reads as "this constant is not in the weight-4 basis" and is simply false. The
substitution `t = e^{-u}` (resp. `1-t = e^{-u}` for the `1/(1-t)` rows) restores
full precision and PSLQ converges immediately.

Rule: before concluding that a constant is *not* expressible in a basis, check
that the quadrature is actually converged — substitute away the endpoint
singularity and re-run.

## Talking to the ChatGPT bridge from a background job

`ask-gpt.py` routes by the tmux window name it reads from `$TMUX_PANE`. A
background job has no stable pane, so the window resolves inconsistently (`rc`,
`dm`, `ds` on three consecutive calls) and questions land in whatever channel
group that resolves to — one of ours ended up in `life2`. Pin it:

```
TMUX_PANE=%7 python3 scripts/ask-gpt.py "<question>"
```

where `%7` is the pane of the `ccdex` window (`tmux list-panes -a -F '#{pane_id} #W'`).

`⚡ ALL CONNECTORS FAILED` is a delivery/polling timeout, not a failure — the tab
is still working. Do not resend.
