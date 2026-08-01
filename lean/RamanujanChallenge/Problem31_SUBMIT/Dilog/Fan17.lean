import Ramanujan31.Dilog.Regular17

/-!
# A structural fan exchange in the regular seventeen-gon

The beta endpoint is governed by two fan decompositions of the regular ideal
seventeen-gon, with coarse gaps `6θ` and `11θ`.  Their difference is changed
to five standard terms by six pentagon flips.  This file records that short
geometric path; no endpoint-specific search data occur in the theorem.
-/

open scoped BigOperators

namespace Real

/-- A Rogers sine cross-ratio with all four angular gaps measured in units of
`θ`. -/
noncomputable def scaledCyclicRogers
    (θ a b c d : ℝ) : ℝ :=
  rogers (sineCross (a * θ) (b * θ) (c * θ) (d * θ))

private theorem scaledCyclicRogers_reflect
    {θ a b c d : ℝ}
    (hsum : (a + b + c + d) * θ = Real.pi) :
    scaledCyclicRogers θ a d c b =
      scaledCyclicRogers θ a b c d := by
  unfold scaledCyclicRogers
  apply congrArg rogers
  apply sineCross_reflect
  calc
    a * θ + b * θ + c * θ + d * θ =
        (a + b + c + d) * θ := by ring
    _ = Real.pi := hsum

private theorem scaledCyclicRogers_rot_two
    {θ a b c d : ℝ}
    (hsum : (a + b + c + d) * θ = Real.pi) :
    scaledCyclicRogers θ c d a b =
      scaledCyclicRogers θ a b c d := by
  unfold scaledCyclicRogers
  apply congrArg rogers
  apply sineCross_rot_two
  calc
    a * θ + b * θ + c * θ + d * θ =
        (a + b + c + d) * θ := by ring
    _ = Real.pi := hsum

private theorem scaledCyclicRogers_swap_outer
    (θ a b c d : ℝ) :
    scaledCyclicRogers θ a b c d =
      scaledCyclicRogers θ c b a d := by
  unfold scaledCyclicRogers sineCross
  ring

private theorem scaledCyclicRogers_euler
    {θ a b c d : ℝ}
    (hθ : 0 < θ) (ha : 0 < a) (hb : 0 < b)
    (hc : 0 < c) (hd : 0 < d)
    (hsum : (a + b + c + d) * θ = Real.pi) :
    scaledCyclicRogers θ a b c d
        + scaledCyclicRogers θ b a d c = rogers 1 := by
  have hsum' :
      a * θ + b * θ + c * θ + d * θ = Real.pi := by
    calc
      a * θ + b * θ + c * θ + d * θ =
          (a + b + c + d) * θ := by ring
      _ = Real.pi := hsum
  have hrot := sineCross_rot_one
    (mul_pos ha hθ) (mul_pos hb hθ) (mul_pos hc hθ) (mul_pos hd hθ) hsum'
  have href := sineCross_reflect
    (A := b * θ) (B := c * θ) (C := d * θ) (D := a * θ)
    (by linarith)
  have hmem := sineCross_mem_Ioo
    (mul_pos ha hθ) (mul_pos hb hθ) (mul_pos hc hθ) (mul_pos hd hθ) hsum'
  unfold scaledCyclicRogers
  rw [href, hrot]
  calc
    rogers (sineCross (a * θ) (b * θ) (c * θ) (d * θ))
          + rogers (1 - sineCross (a * θ) (b * θ) (c * θ) (d * θ)) =
        Real.pi ^ 2 / 6 := rogers_add_rogers_one_sub hmem.1 hmem.2
    _ = rogers 1 := rogers_one.symm

private theorem scaledCyclicRogers_five
    {θ a b c d e : ℝ}
    (hθ : 0 < θ) (ha : 0 < a) (hb : 0 < b) (hc : 0 < c)
    (hd : 0 < d) (he : 0 < e)
    (hsum : (a + b + c + d + e) * θ = Real.pi) :
    scaledCyclicRogers θ b c d (e + a)
      - scaledCyclicRogers θ (a + b) c d e
      + scaledCyclicRogers θ a (b + c) d e
      - scaledCyclicRogers θ a b (c + d) e
      + scaledCyclicRogers θ a b c (d + e) = 0 := by
  have hsum' :
      a * θ + b * θ + c * θ + d * θ + e * θ = Real.pi := by
    calc
      a * θ + b * θ + c * θ + d * θ + e * θ =
          (a + b + c + d + e) * θ := by ring
      _ = Real.pi := hsum
  have h := rogers_sineCross_five
    (mul_pos ha hθ) (mul_pos hb hθ) (mul_pos hc hθ)
    (mul_pos hd hθ) (mul_pos he hθ) hsum'
  unfold scaledCyclicRogers
  simpa only [add_mul] using h

/-- The five positive Rogers terms used at the beta endpoint. -/
noncomputable def betaStandard17Sum (θ : ℝ) : ℝ :=
  (∑ i ∈ Finset.range 4,
      scaledCyclicRogers θ 1 (i + 1) 1 (14 - i))
    + scaledCyclicRogers θ 1 5 6 5

/-- The `6θ` fan based at a unit side of the regular seventeen-gon. -/
noncomputable def fanSix17 (θ : ℝ) : ℝ :=
  ∑ i ∈ Finset.range 9,
    scaledCyclicRogers θ 1 (i + 1) 6 (9 - i)

/-- The complementary `11θ` fan. -/
noncomputable def fanEleven17 (θ : ℝ) : ℝ :=
  ∑ i ∈ Finset.range 4,
    scaledCyclicRogers θ 1 (i + 1) 11 (4 - i)

/--
Six pentagon flips change the beta five-term expression into the difference
of the `6θ` and `11θ` fans.
-/
theorem betaStandard17Sum_eq_fanDifference
    {θ : ℝ} (hθ : 0 < θ) (hsum : 17 * θ = Real.pi) :
    betaStandard17Sum θ = fanSix17 θ - fanEleven17 θ := by
  have hE1 :
      scaledCyclicRogers θ 1 3 2 11
          + scaledCyclicRogers θ 3 1 11 2 = rogers 1 :=
    scaledCyclicRogers_euler hθ (by norm_num) (by norm_num)
      (by norm_num) (by norm_num) (by norm_num at hsum ⊢; linarith)
  have hE2 :
      scaledCyclicRogers θ 1 4 6 6
          + scaledCyclicRogers θ 4 1 6 6 = rogers 1 :=
    scaledCyclicRogers_euler hθ (by norm_num) (by norm_num)
      (by norm_num) (by norm_num) (by norm_num at hsum ⊢; linarith)

  have hF1 := scaledCyclicRogers_five
    hθ (by norm_num) (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    (show (1 + 1 + 1 + 11 + 3 : ℝ) * θ = Real.pi by norm_num at hsum ⊢; linarith)
  have hF2 := scaledCyclicRogers_five
    hθ (by norm_num) (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    (show (1 + 1 + 6 + 6 + 3 : ℝ) * θ = Real.pi by norm_num at hsum ⊢; linarith)
  have hF3 := scaledCyclicRogers_five
    hθ (by norm_num) (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    (show (1 + 1 + 11 + 1 + 3 : ℝ) * θ = Real.pi by norm_num at hsum ⊢; linarith)
  have hF4 := scaledCyclicRogers_five
    hθ (by norm_num) (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    (show (1 + 2 + 1 + 11 + 2 : ℝ) * θ = Real.pi by norm_num at hsum ⊢; linarith)
  have hF5 := scaledCyclicRogers_five
    hθ (by norm_num) (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    (show (1 + 2 + 6 + 6 + 2 : ℝ) * θ = Real.pi by norm_num at hsum ⊢; linarith)
  have hF6 := scaledCyclicRogers_five
    hθ (by norm_num) (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    (show (1 + 3 + 6 + 6 + 1 : ℝ) * θ = Real.pi by norm_num at hsum ⊢; linarith)
  norm_num at hF1 hF2 hF3 hF4 hF5 hF6

  have hr1664 :
      scaledCyclicRogers θ 1 6 6 4 =
        scaledCyclicRogers θ 1 4 6 6 :=
    scaledCyclicRogers_reflect (by norm_num at hsum ⊢; linarith)
  have hr2663 :
      scaledCyclicRogers θ 2 6 6 3 =
        scaledCyclicRogers θ 2 3 6 6 :=
    scaledCyclicRogers_reflect (by norm_num at hsum ⊢; linarith)
  have hr1763 :
      scaledCyclicRogers θ 1 7 6 3 =
        scaledCyclicRogers θ 1 3 6 7 :=
    scaledCyclicRogers_reflect (by norm_num at hsum ⊢; linarith)
  have hr1362 :
      scaledCyclicRogers θ 1 3 11 2 =
        scaledCyclicRogers θ 1 2 11 3 :=
    scaledCyclicRogers_reflect (by norm_num at hsum ⊢; linarith)
  have hr3662 :
      scaledCyclicRogers θ 3 6 6 2 =
        scaledCyclicRogers θ 3 2 6 6 :=
    scaledCyclicRogers_reflect (by norm_num at hsum ⊢; linarith)
  have hr1862 :
      scaledCyclicRogers θ 1 8 6 2 =
        scaledCyclicRogers θ 1 2 6 8 :=
    scaledCyclicRogers_reflect (by norm_num at hsum ⊢; linarith)
  have hr4661 :
      scaledCyclicRogers θ 4 6 6 1 =
        scaledCyclicRogers θ 4 1 6 6 :=
    scaledCyclicRogers_reflect (by norm_num at hsum ⊢; linarith)
  have hr1961 :
      scaledCyclicRogers θ 1 9 6 1 =
        scaledCyclicRogers θ 1 1 6 9 :=
    scaledCyclicRogers_reflect (by norm_num at hsum ⊢; linarith)
  have hr13121 :
      scaledCyclicRogers θ 1 3 12 1 =
        scaledCyclicRogers θ 1 1 12 3 :=
    scaledCyclicRogers_reflect (by norm_num at hsum ⊢; linarith)
  have hr14111 :
      scaledCyclicRogers θ 1 4 11 1 =
        scaledCyclicRogers θ 1 1 11 4 :=
    scaledCyclicRogers_reflect (by norm_num at hsum ⊢; linarith)
  have ht1114 :
      scaledCyclicRogers θ 1 11 1 4 =
        scaledCyclicRogers θ 1 4 1 11 :=
    scaledCyclicRogers_rot_two (by norm_num at hsum ⊢; linarith)
  have ht2113 :
      scaledCyclicRogers θ 2 11 1 3 =
        scaledCyclicRogers θ 1 3 2 11 :=
    scaledCyclicRogers_rot_two (by norm_num at hsum ⊢; linarith)
  have ht1123 :
      scaledCyclicRogers θ 1 12 1 3 =
        scaledCyclicRogers θ 1 3 1 12 :=
    scaledCyclicRogers_rot_two (by norm_num at hsum ⊢; linarith)

  rw [hr1664, hr2663, hr1763] at hF2
  rw [ht1114, ht2113, ht1123] at hF3
  rw [hr1362] at hF4
  rw [hr2663, hr3662, hr1862] at hF5
  rw [hr3662, hr4661, hr1961, hr13121] at hF6

  norm_num [betaStandard17Sum, fanSix17, fanEleven17,
    Finset.sum_range_succ]
  rw [hr1664, hr1763, hr1862, hr1961, hr1362, hr14111]
  linarith

private theorem oneLargePolygonSum_six_decompose (θ : ℝ) :
    oneLargePolygonSum (6 * θ) θ 8 =
      fanSix17 θ
        + ∑ i ∈ Finset.range 9,
            (9 - i : ℕ) * rogers (regular17Argument θ i) := by
  unfold oneLargePolygonSum fanSix17
  apply congrArg₂ (· + ·)
  · apply Finset.sum_congr rfl
    intro i hi
    change
      scaledCyclicRogers θ 6 (i + 1) 1 (9 - i) =
        scaledCyclicRogers θ 1 (i + 1) 6 (9 - i)
    exact scaledCyclicRogers_swap_outer _ _ _ _ _
  · apply Finset.sum_congr rfl
    intro i hi
    simp only [regular17Argument, sineCross]

private theorem oneLargePolygonSum_eleven_decompose (θ : ℝ) :
    oneLargePolygonSum (11 * θ) θ 3 =
      fanEleven17 θ
        + ∑ i ∈ Finset.range 4,
            (4 - i : ℕ) * rogers (regular17Argument θ i) := by
  unfold oneLargePolygonSum fanEleven17
  apply congrArg₂ (· + ·)
  · apply Finset.sum_congr rfl
    intro i hi
    change
      scaledCyclicRogers θ 11 (i + 1) 1 (4 - i) =
        scaledCyclicRogers θ 1 (i + 1) 11 (4 - i)
    exact scaledCyclicRogers_swap_outer _ _ _ _ _
  · apply Finset.sum_congr rfl
    intro i hi
    simp only [regular17Argument, sineCross]

private theorem regular17_weight_difference
    {θ : ℝ} (hsum : 17 * θ = Real.pi) :
    (∑ i ∈ Finset.range 9,
        (9 - i : ℕ) * rogers (regular17Argument θ i))
      - (∑ i ∈ Finset.range 4,
          (4 - i : ℕ) * rogers (regular17Argument θ i))
        =
      5 * ∑ i ∈ Finset.range 7, rogers (regular17Argument θ i) := by
  have h85 :
      rogers (regular17Argument θ 8) =
        rogers (regular17Argument θ 5) := by
    have h := scaledCyclicRogers_reflect
      (θ := θ) (a := 1) (b := 6) (c := 1) (d := 9)
      (by norm_num at hsum ⊢; linarith)
    simpa [scaledCyclicRogers, regular17Argument] using h.symm
  have h76 :
      rogers (regular17Argument θ 7) =
        rogers (regular17Argument θ 6) := by
    have h := scaledCyclicRogers_reflect
      (θ := θ) (a := 1) (b := 7) (c := 1) (d := 8)
      (by norm_num at hsum ⊢; linarith)
    simpa [scaledCyclicRogers, regular17Argument] using h.symm
  norm_num [Finset.sum_range_succ]
  rw [h85, h76]
  ring

/-- The two complementary fans have regulator `15/17` of `R(1)`. -/
theorem rogers_fanDifference17
    {θ : ℝ} (hθ : 0 < θ) (hsum : 17 * θ = Real.pi) :
    17 * (fanSix17 θ - fanEleven17 θ) = 15 * rogers 1 := by
  have hsix :
      oneLargePolygonSum (6 * θ) θ 8 = 9 * rogers 1 := by
    apply oneLargePolygonSum_eq (by positivity) hθ
    norm_num at hsum ⊢
    linarith
  have heleven :
      oneLargePolygonSum (11 * θ) θ 3 = 4 * rogers 1 := by
    apply oneLargePolygonSum_eq (by positivity) hθ
    norm_num at hsum ⊢
    linarith
  rw [oneLargePolygonSum_six_decompose] at hsix
  rw [oneLargePolygonSum_eleven_decompose] at heleven
  have hweight := regular17_weight_difference hsum
  have hregular := rogers_regular17 hθ hsum
  linarith

/-- The five beta terms have the exact regular-seventeen-gon value. -/
theorem rogers_betaStandard17
    {θ : ℝ} (hθ : 0 < θ) (hsum : 17 * θ = Real.pi) :
    17 * betaStandard17Sum θ = 15 * rogers 1 := by
  rw [betaStandard17Sum_eq_fanDifference hθ hsum]
  exact rogers_fanDifference17 hθ hsum

end Real
