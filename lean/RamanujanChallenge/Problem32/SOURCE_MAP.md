# Problem 3.2 source map

Lean declaration names are stable semantic identifiers. Paper theorem numbers
are snapshot metadata only: when the manuscript is reordered, update this map
rather than renaming declarations. The LaTeX `\label` is the primary link to
the manuscript claim.

## Manuscript snapshot

- Repository revision: `6ba7ed097a2f6823e873f7a62904189a526b24ec`
- `proof.tex` blob: `b0db45bb6f9383cf1646437ba8e6279e60abed3f`
- Display numbers: read from `proof.aux` dated 2026-08-09

## Zero-count chain

| Claim ID | Current paper claim | LaTeX label | Lean declaration(s) | Status |
|---|---|---|---|---|
| `P32-ZERO-NO-CONSEC` | Lemma 6, “No consecutive zeros” | `lem:no-consec` | `no_consecutive_zeros` | `exact` |
| `P32-GAP-RECURRENCE` | Lemma 7, “Gap polynomial” (canonical numerator recurrence) | `lem:gap-poly` | `gapPolynomial`, `gapPolynomial_succ`, `eval_gapPolynomial_succ` | `exact-component` |
| `P32-GAP-DEGREE-UPPER` | Lemma 7, degree at most `3(h-1)` | `lem:gap-poly` | `gapPolynomial_natDegree_le` | `exact-component` |
| `P32-GAP-DEGREE-EXACT` | Lemma 7, exact degree and positive leading coefficient | `lem:gap-poly` | — | `open` |
| `P32-GAP-NONVANISH` | Lemma 8, “Nonvanishing over F_p” | `lem:nonvanish` | — | `open` |
| `P32-ZERO-EFFECTIVE` | Proposition 9, effective `Z(p)` bound | `prop:zp-bound` | — | `open` |
| `P32-ZERO-SUBINTERVAL` | Corollary 10, “Sub-interval zero count” | `cor:subinterval` | — | `open` |
| `P32-ZERO-SUBLINEAR` | explicit real-rpow corollary used by the formal capstone | `prop:zp-bound` | `zero_count_sublinear` | `sorry` |

Status meanings:

- `exact`: the Lean theorem matches the complete paper claim.
- `exact-component`: a separately stated algebraic component matches exactly,
  but the surrounding paper lemma has further unformalized conclusions.
- `open`: no corresponding Lean theorem has yet been accepted.
- `sorry`: the declaration exists but is not yet kernel-checked.

## Source-maintenance rule

For a paper edit, first classify the change:

1. Renumbering only: update “Current paper claim” and the snapshot metadata.
2. Wording or proof change with the same mathematical statement: update the
   snapshot metadata; Lean names remain unchanged.
3. Mathematical statement change: mark the row `stale`, add the replacement
   claim as a new row, and re-check every dependent declaration before changing
   the status back to `exact` or `exact-component`.
