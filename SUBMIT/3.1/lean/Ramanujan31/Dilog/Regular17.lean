import Ramanujan31.Dilog.IdealPolygon

/-!
# The regular ideal seventeen-gon identity

The general vertex-removal theorem specializes to a regular ideal
seventeen-gon.  Reflection pairs its fourteen possible side separations,
leaving seven sine-ratio arguments, each with multiplicity seventeen.
-/

open scoped BigOperators

namespace Real

private theorem weighted_reflection_fourteen
    (f : ℕ → ℝ)
    (hsym : ∀ i < 14, f (13 - i) = f i) :
    (∑ i ∈ Finset.range 14, ((15 - i : ℕ) : ℝ) * f i) =
      17 * ∑ i ∈ Finset.range 7, f i := by
  let W : ℕ → ℝ := fun i => ((15 - i : ℕ) : ℝ) * f i
  have hmirror :
      (∑ i ∈ Finset.range 14, W i) =
        ∑ i ∈ Finset.range 14, ((i + 2 : ℕ) : ℝ) * f i := by
    calc
      (∑ i ∈ Finset.range 14, W i) =
          ∑ i ∈ Finset.range 14, W (14 - 1 - i) :=
        (Finset.sum_range_reflect W 14).symm
      _ = ∑ i ∈ Finset.range 14, ((i + 2 : ℕ) : ℝ) * f i := by
        apply Finset.sum_congr rfl
        intro i hi
        have hi' : i < 14 := Finset.mem_range.mp hi
        have hindex : 14 - 1 - i = 13 - i := by omega
        have hweight : 15 - (13 - i) = i + 2 := by omega
        simp only [W, hindex, hweight, hsym i hi']

  have htotal :
      (∑ i ∈ Finset.range 14, f i) =
        2 * ∑ i ∈ Finset.range 7, f i := by
    have hsecond :
        (∑ i ∈ Finset.range 7, f (i + 7)) =
          ∑ i ∈ Finset.range 7, f i := by
      calc
        (∑ i ∈ Finset.range 7, f (i + 7)) =
            ∑ i ∈ Finset.range 7, f ((7 - 1 - i) + 7) :=
          (Finset.sum_range_reflect (fun i => f (i + 7)) 7).symm
        _ = ∑ i ∈ Finset.range 7, f i := by
          apply Finset.sum_congr rfl
          intro i hi
          have hi' : i < 7 := Finset.mem_range.mp hi
          have hindex : (7 - 1 - i) + 7 = 13 - i := by omega
          rw [hindex, hsym i (by omega)]
    calc
      (∑ i ∈ Finset.range 14, f i) =
          (∑ i ∈ Finset.range 7, f i)
            + ∑ i ∈ Finset.range 7, f (i + 7) := by
        simpa using Finset.sum_range_add f 7 7
      _ = 2 * ∑ i ∈ Finset.range 7, f i := by rw [hsecond]; ring

  have hdouble :
      2 * (∑ i ∈ Finset.range 14, W i) =
        17 * ∑ i ∈ Finset.range 14, f i := by
    calc
      2 * (∑ i ∈ Finset.range 14, W i) =
          (∑ i ∈ Finset.range 14, W i)
            + ∑ i ∈ Finset.range 14, W i := by ring
      _ = (∑ i ∈ Finset.range 14, W i)
            + ∑ i ∈ Finset.range 14, ((i + 2 : ℕ) : ℝ) * f i := by
        rw [hmirror]
      _ = ∑ i ∈ Finset.range 14,
            (W i + ((i + 2 : ℕ) : ℝ) * f i) := by
        rw [Finset.sum_add_distrib]
      _ = ∑ i ∈ Finset.range 14, 17 * f i := by
        apply Finset.sum_congr rfl
        intro i hi
        have hi' : i < 14 := Finset.mem_range.mp hi
        have hweight : 15 - i + (i + 2) = 17 := by omega
        have hweightR :
            (((15 - i : ℕ) : ℝ) + ((i + 2 : ℕ) : ℝ)) = 17 := by
          exact_mod_cast hweight
        simp only [W]
        calc
          ((15 - i : ℕ) : ℝ) * f i + ((i + 2 : ℕ) : ℝ) * f i =
              (((15 - i : ℕ) : ℝ) + ((i + 2 : ℕ) : ℝ)) * f i := by ring
          _ = 17 * f i := by rw [hweightR]
      _ = 17 * ∑ i ∈ Finset.range 14, f i := by
        rw [Finset.mul_sum]

  simp only [W] at hdouble ⊢
  rw [htotal] at hdouble
  linarith

/-- The cross-ratio attached to separation `i + 2` in the regular
seventeen-gon. -/
noncomputable def regular17Argument (θ : ℝ) (i : ℕ) : ℝ :=
  sineCross θ ((i + 1 : ℕ) * θ) θ ((14 - i : ℕ) * θ)

theorem regular17Argument_eq_sineRatio (θ : ℝ) (i : ℕ) :
    regular17Argument θ i =
      (sin θ / sin ((i + 2 : ℕ) * θ)) ^ 2 := by
  have hangle :
      θ + (i + 1 : ℕ) * θ = (i + 2 : ℕ) * θ := by
    push_cast
    ring
  have hangle' :
      (i + 1 : ℕ) * θ + θ = (i + 2 : ℕ) * θ := by
    push_cast
    ring
  unfold regular17Argument sineCross
  rw [hangle, hangle']
  ring

private theorem oneLargePolygonSum_self_thirteen
    (θ : ℝ) :
    oneLargePolygonSum θ θ 13 =
      ∑ i ∈ Finset.range 14, ((15 - i : ℕ) : ℝ) *
        rogers (regular17Argument θ i) := by
  simp only [oneLargePolygonSum]
  rw [← Finset.sum_add_distrib]
  apply Finset.sum_congr rfl
  intro i hi
  have hi' : i < 14 := Finset.mem_range.mp hi
  have hrem : 14 - i = (13 - i) + 1 := by omega
  have hweight : 15 - i = (14 - i) + 1 := by omega
  have hremR :
      (((14 - i : ℕ) : ℝ)) = ((13 - i : ℕ) : ℝ) + 1 := by
    exact_mod_cast hrem
  have hweightR :
      (((15 - i : ℕ) : ℝ)) = ((14 - i : ℕ) : ℝ) + 1 := by
    exact_mod_cast hweight
  simp only [regular17Argument]
  rw [hrem, hweight]
  push_cast
  rw [hremR]
  ring

/-- The regular ideal seventeen-gon Rogers identity. -/
theorem rogers_regular17
    {θ : ℝ} (hθ : 0 < θ) (hsum : 17 * θ = Real.pi) :
    17 * (∑ i ∈ Finset.range 7, rogers (regular17Argument θ i)) =
      14 * rogers 1 := by
  let f : ℕ → ℝ := fun i => rogers (regular17Argument θ i)
  have hsym : ∀ i < 14, f (13 - i) = f i := by
    intro i hi
    have hfill :
        θ + (i + 1 : ℕ) * θ + θ + (14 - i : ℕ) * θ = Real.pi := by
      push_cast at hsum ⊢
      have hiR : ((i : ℕ) : ℝ) ≤ 13 := by exact_mod_cast (show i ≤ 13 by omega)
      have hsubR : (((14 - i : ℕ) : ℝ)) = 14 - i := by
        rw [Nat.cast_sub (by omega)]
        norm_num
      rw [hsubR]
      nlinarith
    have href := sineCross_reflect hfill
    simp only [f, regular17Argument]
    have hleft : 13 - i + 1 = 14 - i := by omega
    have hright : 14 - (13 - i) = i + 1 := by omega
    rw [hleft, hright]
    exact congrArg rogers href
  have hpair := weighted_reflection_fourteen f hsym
  have hpolygon :=
    regularPolygonSum_eq (N := 13) hθ (by simpa using hsum)
  rw [oneLargePolygonSum_self_thirteen θ] at hpolygon
  have hpolygon' :
      (∑ i ∈ Finset.range 14, ((15 - i : ℕ) : ℝ) *
          rogers (regular17Argument θ i)) =
        14 * rogers 1 := by
    convert hpolygon using 1
  simp only [f] at hpair
  linarith [hpolygon']

/-- Sine-ratio form of the regular seventeen-gon identity. -/
theorem rogers_regular17_sineRatios :
    17 * (∑ i ∈ Finset.range 7,
      rogers
        ((sin (Real.pi / 17) /
          sin ((i + 2 : ℕ) * (Real.pi / 17))) ^ 2)) =
      14 * rogers 1 := by
  have hpi : 0 < Real.pi := Real.pi_pos
  have h :=
    rogers_regular17 (θ := Real.pi / 17) (by positivity) (by field_simp)
  simpa only [regular17Argument_eq_sineRatio] using h

end Real
