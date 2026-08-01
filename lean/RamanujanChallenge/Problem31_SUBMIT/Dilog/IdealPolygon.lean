import Ramanujan31.Dilog.CyclicRogers

/-!
# Rogers sums for ideal polygons

The main lemma in this file is the structural vertex-removal identity for an
ideal polygon.  Two adjacent angular gaps `A` and `B` are merged.  The Rogers
terms involving the removed vertex form a telescoping sum of five-term
relations, and their total is exactly `R(1)`.

This is the induction step behind the finite ideal-polygon orthospectrum
identity.  It replaces endpoint-specific linear combinations of five-term
relations by one reusable geometric argument.
-/

open scoped BigOperators

namespace Real

private theorem sum_range_succ_sub (F : ℕ → ℝ) (n : ℕ) :
    (∑ i ∈ Finset.range n, (F (i + 1) - F i)) = F n - F 0 := by
  induction n with
  | zero => simp
  | succ n ih =>
      rw [Finset.sum_range_succ, ih]
      ring

/--
Removing a vertex of an ideal polygon contributes one copy of `R(1)`.

`C i` is the angular gap accumulated to the left of the moving side and
`D i` is the remaining angular gap to its right.  Thus

```
C 0 = B,       D N = B,
C (i+1) = C i + B,
D i = B + D (i+1).
```

For every `i`, the four gaps `A, B, C i, D i` fill a semicircle.  Each
summand is the difference between the old and new polygon after merging the
first two gaps.  Abel's five-term relation turns it into `F (i+1) - F i`;
the sum telescopes, and the two boundary terms are Euler complements.
-/
theorem rogers_idealPolygon_remove_vertex
    {N : ℕ} {A B : ℝ} {C D : ℕ → ℝ}
    (hA : 0 < A) (hB : 0 < B)
    (hC : ∀ i ≤ N, 0 < C i)
    (hD : ∀ i ≤ N, 0 < D i)
    (hC0 : C 0 = B)
    (hDN : D N = B)
    (hCstep : ∀ i < N, C (i + 1) = C i + B)
    (hDstep : ∀ i < N, D i = B + D (i + 1))
    (hsum : ∀ i ≤ N, A + B + C i + D i = Real.pi) :
    rogers (sineCross A B (C 0) (D 0))
        + (∑ i ∈ Finset.range N,
            (rogers (sineCross B (C i) B (D (i + 1) + A))
              + rogers (sineCross A (B + C i) B (D (i + 1)))
              - rogers (sineCross (A + B) (C i) B (D (i + 1)))))
        + rogers (sineCross B (C N) B A)
      = rogers 1 := by
  let F : ℕ → ℝ := fun i => rogers (sineCross A B (C i) (D i))

  have hlocal :
      ∀ i < N,
        rogers (sineCross B (C i) B (D (i + 1) + A))
            + rogers (sineCross A (B + C i) B (D (i + 1)))
            - rogers (sineCross (A + B) (C i) B (D (i + 1)))
          = F (i + 1) - F i := by
    intro i hi
    have hCi : 0 < C i := hC i (Nat.le_of_lt hi)
    have hDi1 : 0 < D (i + 1) := hD (i + 1) (by omega)
    have hsumi : A + B + C i + B + D (i + 1) = Real.pi := by
      have hi' := hsum i (Nat.le_of_lt hi)
      have hd' := hDstep i hi
      linarith
    have hfive := rogers_sineCross_five
      (A := A) (B := B) (C := C i) (D := B) (E := D (i + 1))
      hA hB hCi hB hDi1 hsumi
    have hFnext :
        F (i + 1) =
          rogers (sineCross A B (C i + B) (D (i + 1))) := by
      simp only [F, hCstep i hi]
    have hFi :
        F i = rogers (sineCross A B (C i) (B + D (i + 1))) := by
      simp only [F, hDstep i hi]
    rw [hFnext, hFi]
    linarith

  have htel :
      (∑ i ∈ Finset.range N,
          (rogers (sineCross B (C i) B (D (i + 1) + A))
            + rogers (sineCross A (B + C i) B (D (i + 1)))
            - rogers (sineCross (A + B) (C i) B (D (i + 1)))))
        = F N - F 0 := by
    calc
      _ = ∑ i ∈ Finset.range N, (F (i + 1) - F i) := by
        apply Finset.sum_congr rfl
        intro i hi
        exact hlocal i (Finset.mem_range.mp hi)
      _ = F N - F 0 := sum_range_succ_sub F N

  have hboundary :
      rogers (sineCross B (C N) B A) + F N = rogers 1 := by
    have hsumN : A + B + C N + B = Real.pi := by
      have hn := hsum N le_rfl
      rw [hDN] at hn
      exact hn
    have hrot :
        sineCross B (C N) B A =
          1 - sineCross A B (C N) (D N) := by
      rw [hDN]
      exact sineCross_rot_one hA hB (hC N le_rfl) hB hsumN
    have hmem :=
      sineCross_mem_Ioo hA hB (hC N le_rfl) (hD N le_rfl) (hsum N le_rfl)
    simp only [F]
    rw [hrot]
    calc
      rogers (1 - sineCross A B (C N) (D N))
            + rogers (sineCross A B (C N) (D N)) =
          rogers (sineCross A B (C N) (D N))
            + rogers (1 - sineCross A B (C N) (D N)) := add_comm _ _
      _ = Real.pi ^ 2 / 6 := rogers_add_rogers_one_sub hmem.1 hmem.2
      _ = rogers 1 := rogers_one.symm

  rw [htel]
  simp only [F]
  rw [hC0]
  linarith

/--
The vertex-removal identity for a polygon with one angular gap `A` and all
remaining gaps equal to `B`.  The index `N` counts the interior moving
positions after the two boundary terms have been removed.
-/
theorem rogers_oneLargePolygon_step_terms
    {N : ℕ} {A B : ℝ}
    (hA : 0 < A) (hB : 0 < B)
    (hsum : A + (N + 4 : ℕ) * B = Real.pi) :
    rogers (sineCross A B B ((N + 1 + 1 : ℕ) * B))
        + (∑ i ∈ Finset.range (N + 1),
            (rogers
                (sineCross B ((i + 1 : ℕ) * B) B
                  (((N + 1 - i : ℕ) * B) + A))
              + rogers
                  (sineCross A (B + ((i + 1 : ℕ) * B)) B
                    ((N + 1 - i : ℕ) * B))
              - rogers
                  (sineCross (A + B) ((i + 1 : ℕ) * B) B
                    ((N + 1 - i : ℕ) * B))))
        + rogers (sineCross B ((N + 1 + 1 : ℕ) * B) B A)
      = rogers 1 := by
  let C : ℕ → ℝ := fun i => (i + 1 : ℕ) * B
  let D : ℕ → ℝ := fun i => (N + 2 - i : ℕ) * B
  have hC : ∀ i ≤ N + 1, 0 < C i := by
    intro i hi
    exact mul_pos (Nat.cast_pos.mpr (by omega)) hB
  have hD : ∀ i ≤ N + 1, 0 < D i := by
    intro i hi
    exact mul_pos (Nat.cast_pos.mpr (by omega)) hB
  have hC0 : C 0 = B := by
    simp [C]
  have hDN : D (N + 1) = B := by
    simp [D]
  have hCstep : ∀ i < N + 1, C (i + 1) = C i + B := by
    intro i hi
    simp only [C]
    push_cast
    ring
  have hDstep : ∀ i < N + 1, D i = B + D (i + 1) := by
    intro i hi
    have hn : N + 2 - i = 1 + (N + 2 - (i + 1)) := by omega
    simp only [D, hn]
    push_cast
    ring
  have hfill : ∀ i ≤ N + 1, A + B + C i + D i = Real.pi := by
    intro i hi
    have hn : (i + 1) + (N + 2 - i) = N + 3 := by omega
    have hnR :
        ((i + 1 : ℕ) : ℝ) + ((N + 2 - i : ℕ) : ℝ) =
          ((N + 3 : ℕ) : ℝ) := by
      exact_mod_cast hn
    simp only [C, D]
    push_cast at hnR hsum ⊢
    nlinarith [hnR]
  have h := rogers_idealPolygon_remove_vertex
    hA hB hC hD hC0 hDN hCstep hDstep hfill
  simpa [C, D] using h

/--
The Rogers sum over all unordered pairs of non-adjacent sides of a polygon
whose cyclic gaps are

```
A, B, ..., B
```

with `N + 3` copies of `B`.  The first sum contains the pairs involving the
distinguished `A`-side.  In the second sum, sides at a fixed cyclic
separation have been grouped; `N + 1 - i` is their multiplicity.
-/
noncomputable def oneLargePolygonSum (A B : ℝ) (N : ℕ) : ℝ :=
  (∑ i ∈ Finset.range (N + 1),
      rogers
        (sineCross A ((i + 1 : ℕ) * B) B
          ((N + 1 - i : ℕ) * B)))
    + ∑ i ∈ Finset.range (N + 1),
        (N + 1 - i : ℕ) *
          rogers
            (sineCross B ((i + 1 : ℕ) * B) B
              (A + (N - i : ℕ) * B))

/-- The quadrilateral case of `oneLargePolygonSum`. -/
theorem oneLargePolygonSum_zero
    {A B : ℝ} (hA : 0 < A) (hB : 0 < B)
    (hsum : A + 3 * B = Real.pi) :
    oneLargePolygonSum A B 0 = rogers 1 := by
  have hrot :
      sineCross B B B A = 1 - sineCross A B B B :=
    sineCross_rot_one hA hB hB hB (by linarith)
  have hmem :=
    sineCross_mem_Ioo hA hB hB hB (by linarith)
  simp only [oneLargePolygonSum, Finset.sum_range_one, Nat.zero_add,
    Nat.sub_zero, Nat.sub_self, Nat.cast_one, Nat.cast_zero, one_mul, zero_mul,
    add_zero]
  rw [hrot]
  calc
    rogers (sineCross A B B B) + rogers (1 - sineCross A B B B) =
        Real.pi ^ 2 / 6 := rogers_add_rogers_one_sub hmem.1 hmem.2
    _ = rogers 1 := rogers_one.symm

private theorem weighted_sum_range_succ (f : ℕ → ℝ) (n : ℕ) :
    (∑ i ∈ Finset.range (n + 1), (n + 1 - i : ℕ) * f i)
      = (∑ i ∈ Finset.range n, (n - i : ℕ) * f i)
          + ∑ i ∈ Finset.range (n + 1), f i := by
  calc
    (∑ i ∈ Finset.range (n + 1), (n + 1 - i : ℕ) * f i) =
        (∑ i ∈ Finset.range n, (((n - i : ℕ) : ℝ) + 1) * f i) + f n := by
      rw [Finset.sum_range_succ]
      congr 1
      · apply Finset.sum_congr rfl
        intro i hi
        have hin : i < n := Finset.mem_range.mp hi
        have hn : n + 1 - i = (n - i) + 1 := by omega
        rw [hn]
        push_cast
        rfl
      · simp
    _ = (∑ i ∈ Finset.range n, (n - i : ℕ) * f i)
          + (∑ i ∈ Finset.range n, f i) + f n := by
      simp_rw [add_mul, one_mul, Finset.sum_add_distrib]
    _ = (∑ i ∈ Finset.range n, (n - i : ℕ) * f i)
          + ∑ i ∈ Finset.range (n + 1), f i := by
      rw [Finset.sum_range_succ]
      ring

/-- Merging the distinguished gap with its first neighbouring `B`-gap
decreases `oneLargePolygonSum` by exactly one copy of `R(1)`. -/
theorem oneLargePolygonSum_succ
    {N : ℕ} {A B : ℝ} (hA : 0 < A) (hB : 0 < B)
    (hsum : A + (N + 4 : ℕ) * B = Real.pi) :
    oneLargePolygonSum A B (N + 1) =
      oneLargePolygonSum (A + B) B N + rogers 1 := by
  let aOld : ℕ → ℝ := fun i =>
    rogers
      (sineCross A ((i + 1 : ℕ) * B) B
        ((N + 2 - i : ℕ) * B))
  let aNew : ℕ → ℝ := fun i =>
    rogers
      (sineCross (A + B) ((i + 1 : ℕ) * B) B
        ((N + 1 - i : ℕ) * B))
  let bTerm : ℕ → ℝ := fun i =>
    rogers
      (sineCross B ((i + 1 : ℕ) * B) B
        (A + (N + 1 - i : ℕ) * B))

  have ha :
      (∑ i ∈ Finset.range (N + 2), aOld i) =
        rogers (sineCross A B B ((N + 1 + 1 : ℕ) * B))
          + ∑ i ∈ Finset.range (N + 1),
              rogers
                (sineCross A (B + ((i + 1 : ℕ) * B)) B
                  ((N + 1 - i : ℕ) * B)) := by
    have hhead :
        aOld 0 =
          rogers (sineCross A B B ((N + 1 + 1 : ℕ) * B)) := by
      simp only [aOld, Nat.zero_add, Nat.cast_one, one_mul, Nat.sub_zero]
    have htail :
        (∑ i ∈ Finset.range (N + 1), aOld (i + 1)) =
          ∑ i ∈ Finset.range (N + 1),
            rogers
              (sineCross A (B + ((i + 1 : ℕ) * B)) B
                ((N + 1 - i : ℕ) * B)) := by
      apply Finset.sum_congr rfl
      intro i hi
      have hin : i < N + 1 := Finset.mem_range.mp hi
      have hsub : N + 2 - (i + 1) = N + 1 - i := by omega
      simp only [aOld, hsub]
      congr 2
      push_cast
      ring
    rw [Finset.sum_range_succ', hhead, htail]
    ring

  have hbNew :
      (∑ i ∈ Finset.range (N + 1), (N + 1 - i : ℕ) *
          rogers
            (sineCross B ((i + 1 : ℕ) * B) B
              ((A + B) + (N - i : ℕ) * B))) =
        ∑ i ∈ Finset.range (N + 1), (N + 1 - i : ℕ) * bTerm i := by
    apply Finset.sum_congr rfl
    intro i hi
    have hin : i < N + 1 := Finset.mem_range.mp hi
    have hsub : N + 1 - i = (N - i) + 1 := by omega
    simp only [bTerm]
    congr 2

  have hb :
      (∑ i ∈ Finset.range (N + 2), (N + 2 - i : ℕ) * bTerm i) =
        (∑ i ∈ Finset.range (N + 1), (N + 1 - i : ℕ) * bTerm i)
          + ∑ i ∈ Finset.range (N + 2), bTerm i := by
    simpa [Nat.add_assoc] using weighted_sum_range_succ bTerm (N + 1)

  have hbStep :
      (∑ i ∈ Finset.range (N + 1),
          rogers
            (sineCross B ((i + 1 : ℕ) * B) B
              (((N + 1 - i : ℕ) * B) + A))) =
        ∑ i ∈ Finset.range (N + 1), bTerm i := by
    apply Finset.sum_congr rfl
    intro i hi
    simp only [bTerm]
    congr 2
    ring
  have hnew :
      (∑ i ∈ Finset.range (N + 1),
          rogers
            (sineCross (A + B) ((i + 1 : ℕ) * B) B
              ((N + 1 - i : ℕ) * B))) =
        ∑ i ∈ Finset.range (N + 1), aNew i := by
    rfl
  have hfinal :
      rogers (sineCross B ((N + 1 + 1 : ℕ) * B) B A) =
        bTerm (N + 1) := by
    simp [bTerm]
  have hbSplit :
      (∑ i ∈ Finset.range (N + 2), bTerm i) =
        (∑ i ∈ Finset.range (N + 1), bTerm i) + bTerm (N + 1) := by
    simpa [Nat.add_assoc] using Finset.sum_range_succ bTerm (N + 1)

  have hstep := rogers_oneLargePolygon_step_terms
    (N := N) hA hB hsum
  rw [Finset.sum_sub_distrib, Finset.sum_add_distrib] at hstep
  rw [hbStep, hnew, hfinal] at hstep

  simp only [oneLargePolygonSum, Nat.add_assoc]
  change
    (∑ i ∈ Finset.range (N + 2), aOld i)
        + (∑ i ∈ Finset.range (N + 2), (N + 2 - i : ℕ) * bTerm i) =
      (∑ i ∈ Finset.range (N + 1), aNew i)
        + (∑ i ∈ Finset.range (N + 1), (N + 1 - i : ℕ) *
            rogers
              (sineCross B ((i + 1 : ℕ) * B) B
                ((A + B) + (N - i : ℕ) * B)))
        + rogers 1
  rw [ha, hb, hbNew, hbSplit]
  linarith

/-- The ideal-polygon Rogers identity for one distinguished gap. -/
theorem oneLargePolygonSum_eq
    {N : ℕ} {A B : ℝ} (hA : 0 < A) (hB : 0 < B)
    (hsum : A + (N + 3 : ℕ) * B = Real.pi) :
    oneLargePolygonSum A B N = (N + 1 : ℕ) * rogers 1 := by
  induction N generalizing A with
  | zero =>
      simpa using oneLargePolygonSum_zero hA hB hsum
  | succ N ih =>
      have hsumStep : A + (N + 4 : ℕ) * B = Real.pi := by
        convert hsum using 1 <;> omega
      have hsumMerged : (A + B) + (N + 3 : ℕ) * B = Real.pi := by
        push_cast at hsumStep ⊢
        nlinarith
      rw [oneLargePolygonSum_succ hA hB hsumStep]
      rw [ih (by linarith) hsumMerged]
      push_cast
      ring

/-- The ideal-polygon identity for a regular `(N + 4)`-gon. -/
theorem regularPolygonSum_eq
    {N : ℕ} {B : ℝ} (hB : 0 < B)
    (hsum : (N + 4 : ℕ) * B = Real.pi) :
    oneLargePolygonSum B B N = (N + 1 : ℕ) * rogers 1 := by
  apply oneLargePolygonSum_eq hB hB
  push_cast at hsum ⊢
  nlinarith

end Real
