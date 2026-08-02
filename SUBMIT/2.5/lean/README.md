# Lean source snapshot

These files mirror the canonical Problem 2.5 development in
`../../../lean/RamanujanChallenge/`.

They are provided as a readable source bundle. The checked project, Lake
configuration, and pinned Lean/Mathlib versions live in the repository's
top-level `lean/` directory. Verify the final theorem there with:

```bash
lake env lean RamanujanChallenge/Problem25Moment.lean
```

The final theorem is:

```lean
theorem problem25_solved : Problem25Claim
```

The snapshot has no `sorry`, added axiom, or `admit`. The final theorem's
axiom audit is `{propext, Classical.choice, Quot.sound}`.
