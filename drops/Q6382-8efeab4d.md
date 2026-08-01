ANSWER Q6382 8efeab4d

# Verdict

The one-strip contour shift is now implemented as a self-contained Lean 4.29 scratch module and compile-tested against the repository's pinned toolchain and current `Problem27Barnes.lean`.

The implementation proves, without a residue theorem or hidden residue hypothesis:

1. `verticalIntegral_eq_of_horizontal_tendsto` — rectangle Cauchy plus vanishing horizontal edges;
2. a `dslope` extension of `sin (π z)/(z-m)`;
3. nonvanishing of that extension on the closed half-integer strip;
4. equality of the raw Barnes quotient with the extension away from the crossed integer;
5. `barnes_one_strip_shift` — equality of the two whole vertical-line integrals.

The raw totalized quotient is **not** asserted equal to the extension at `z=m`; in Lean it generally has the junk value caused by division by zero. Instead, the extension is holomorphic at `m`, and the two boundary lines do not contain `m`. This is the correct formulation and is exactly what avoids smuggling in a residue assumption.

The current project kernel is

```lean
def zudilinBarnesSquaredSineKernel27 (s : ℂ) : ℂ :=
  ((Real.pi : ℂ) / Complex.sin ((Real.pi : ℂ) * s)) ^ 2
```

so its scalar in the generic code below is `κ = (Real.pi : ℂ) ^ 2`. For the normalization written in the question, `π / sin(πt)^2`, use `κ = (Real.pi : ℂ)`.

# Exact compile test

Audit source commit:

```text
1efb9c29ac20ddc80fea4e9acf8b5e28e6b6c2fb
```

Final audit-workflow commit:

```text
c3bf67a3b581321a43ce1d75467a4213d8dc9706
```

GitHub Actions run:

```text
30698114952
```

The workflow used the repository's `lean/lean-toolchain`:

```text
leanprover/lean4:v4.29.0
```

and Mathlib revision from the current manifest:

```text
8a178386ffc0f5fef0b77738bb5449d50efeea95
```

Exact commands:

```bash
cd lean
lake exe cache get
lake build RamanujanChallenge.Q6382Scratch
```

Relevant output:

```text
✔ [2719/2721] Built RamanujanChallenge.Problem27 (36s)
⚠ [2720/2721] Built RamanujanChallenge.Problem27Barnes (12s)
⚠ [2721/2721] Built RamanujanChallenge.Q6382Scratch (9.1s)
...
Build completed successfully (2721 jobs).
```

The `⚠` markers are linter warnings, not build failures. The scratch warnings are only an unused `push_cast`, an unnecessary `simpa`, and two sequencing-style suggestions. There are no errors, `sorry`, or user axioms.

The four `#print axioms` results were:

```text
'RamanujanChallenge.P27.Q6382.verticalIntegral_eq_of_horizontal_tendsto'
depends on axioms: [propext, Classical.choice, Quot.sound]

'...sineSlope_ne_zero_on_strip'
depends on axioms: [propext, Classical.choice, Quot.sound]

'...barnesRaw_eq_extension_of_mem_strip_of_ne'
depends on axioms: [propext, Classical.choice, Quot.sound]

'RamanujanChallenge.P27.Q6382.barnes_one_strip_shift'
depends on axioms: [propext, Classical.choice, Quot.sound]
```

In particular, none depends on `sorryAx`.

# Complete compiled scratch file

This is the exact source that built in the audit. The helper declarations are `private` only to avoid polluting the project namespace during the audit; for production, remove `private` from the reusable interface and move it to `Problem27BarnesShift.lean`.

```lean
import RamanujanChallenge.Problem27Barnes
import Mathlib.Analysis.Complex.CauchyIntegral
import Mathlib.Analysis.Complex.RemovableSingularity
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Complex
import Mathlib.MeasureTheory.Integral.IntegralEqImproper

open Filter Set Topology
open scoped Interval Real

noncomputable section

namespace RamanujanChallenge.P27.Q6382

private def closedVerticalStrip (a b : ℝ) : Set ℂ :=
  {z | a ≤ z.re ∧ z.re ≤ b}

private def verticalPoint (x y : ℝ) : ℂ :=
  (x : ℂ) + (y : ℂ) * Complex.I

/-- Rectangle Cauchy plus vanishing horizontal sides. These are the
parameterized vertical integrals; multiplying both sides by `I` gives the
usual oriented line integrals. -/
theorem verticalIntegral_eq_of_horizontal_tendsto
    {F : ℂ → ℂ} {a b : ℝ} (hab : a ≤ b)
    (hF : DifferentiableOn ℂ F (closedVerticalStrip a b))
    (hleft : MeasureTheory.Integrable
      (fun y : ℝ => F (verticalPoint a y)))
    (hright : MeasureTheory.Integrable
      (fun y : ℝ => F (verticalPoint b y)))
    (htop : Tendsto
      (fun T : ℝ => ∫ x in a..b,
        F ((x : ℂ) + (T : ℂ) * Complex.I))
      atTop (𝓝 0))
    (hbottom : Tendsto
      (fun T : ℝ => ∫ x in a..b,
        F ((x : ℂ) - (T : ℂ) * Complex.I))
      atTop (𝓝 0)) :
    (∫ y : ℝ, F (verticalPoint a y)) =
      ∫ y : ℝ, F (verticalPoint b y) := by
  have hfinite : ∀ T : ℝ, 0 ≤ T →
      (∫ y in -T..T, F (verticalPoint a y)) -
          (∫ y in -T..T, F (verticalPoint b y)) =
        Complex.I *
          ((∫ x in a..b,
              F ((x : ℂ) + (T : ℂ) * Complex.I)) -
           (∫ x in a..b,
              F ((x : ℂ) - (T : ℂ) * Complex.I))) := by
    intro T hT
    let z : ℂ := (a : ℂ) - (T : ℂ) * Complex.I
    let w : ℂ := (b : ℂ) + (T : ℂ) * Complex.I
    have hrect :=
      Complex.integral_boundary_rect_eq_zero_of_differentiableOn
        F z w (hF.mono (by
          intro u hu
          have hre : u.re ∈ [[a, b]] := by
            simpa [z, w] using hu.1
          simpa [closedVerticalStrip, uIcc_of_le hab] using hre))
    have hrect' :
        (∫ x in a..b,
            F ((x : ℂ) - (T : ℂ) * Complex.I)) -
          (∫ x in a..b,
            F ((x : ℂ) + (T : ℂ) * Complex.I)) +
          Complex.I * (∫ y in -T..T, F (verticalPoint b y)) -
          Complex.I * (∫ y in -T..T, F (verticalPoint a y)) = 0 := by
      simpa [z, w, verticalPoint, smul_eq_mul] using hrect
    have hIR :
        Complex.I *
            ((∫ y in -T..T, F (verticalPoint b y)) -
              (∫ y in -T..T, F (verticalPoint a y))) =
          (∫ x in a..b,
              F ((x : ℂ) + (T : ℂ) * Complex.I)) -
            (∫ x in a..b,
              F ((x : ℂ) - (T : ℂ) * Complex.I)) := by
      linear_combination hrect'
    calc
      (∫ y in -T..T, F (verticalPoint a y)) -
          (∫ y in -T..T, F (verticalPoint b y)) =
          -((∫ y in -T..T, F (verticalPoint b y)) -
            (∫ y in -T..T, F (verticalPoint a y))) := by ring
      _ = Complex.I ^ 2 *
          ((∫ y in -T..T, F (verticalPoint b y)) -
            (∫ y in -T..T, F (verticalPoint a y))) := by
        rw [Complex.I_sq]
        ring
      _ = Complex.I *
          (Complex.I *
            ((∫ y in -T..T, F (verticalPoint b y)) -
              (∫ y in -T..T, F (verticalPoint a y)))) := by ring
      _ = Complex.I *
          ((∫ x in a..b,
              F ((x : ℂ) + (T : ℂ) * Complex.I)) -
            (∫ x in a..b,
              F ((x : ℂ) - (T : ℂ) * Complex.I))) := by rw [hIR]

  have hleft_lim : Tendsto
      (fun T : ℝ => ∫ y in -T..T, F (verticalPoint a y))
      atTop (𝓝 (∫ y : ℝ, F (verticalPoint a y))) :=
    MeasureTheory.intervalIntegral_tendsto_integral
      hleft tendsto_neg_atTop_atBot tendsto_id

  have hright_lim : Tendsto
      (fun T : ℝ => ∫ y in -T..T, F (verticalPoint b y))
      atTop (𝓝 (∫ y : ℝ, F (verticalPoint b y))) :=
    MeasureTheory.intervalIntegral_tendsto_integral
      hright tendsto_neg_atTop_atBot tendsto_id

  have hfinite_eventually : ∀ᶠ T : ℝ in atTop,
      (∫ y in -T..T, F (verticalPoint a y)) -
          (∫ y in -T..T, F (verticalPoint b y)) =
        Complex.I *
          ((∫ x in a..b,
              F ((x : ℂ) + (T : ℂ) * Complex.I)) -
           (∫ x in a..b,
              F ((x : ℂ) - (T : ℂ) * Complex.I))) :=
    (eventually_ge_atTop (0 : ℝ)).mono hfinite

  have hlhs := hleft_lim.sub hright_lim
  have hrhs : Tendsto
      (fun T : ℝ => Complex.I *
        ((∫ x in a..b,
            F ((x : ℂ) + (T : ℂ) * Complex.I)) -
         (∫ x in a..b,
            F ((x : ℂ) - (T : ℂ) * Complex.I))))
      atTop (𝓝 0) := by
    simpa using tendsto_const_nhds.mul (htop.sub hbottom)
  have hlhs_zero := hrhs.congr'
    (hfinite_eventually.mono fun _ h => h.symm)
  have hsub :
      (∫ y : ℝ, F (verticalPoint a y)) -
          (∫ y : ℝ, F (verticalPoint b y)) = 0 :=
    tendsto_nhds_unique hlhs hlhs_zero
  exact sub_eq_zero.mp hsub

private def sinePi (z : ℂ) : ℂ :=
  Complex.sin ((Real.pi : ℂ) * z)

/-- Holomorphic extension of `sin (π z) / (z-m)` at the integer `m`. -/
private def sineSlope (m : ℤ) : ℂ → ℂ :=
  dslope sinePi (m : ℂ)

@[simp] private theorem sinePi_int (m : ℤ) :
    sinePi (m : ℂ) = 0 := by
  rw [sinePi, Complex.sin_eq_zero_iff]
  exact ⟨m, by push_cast; ring⟩

private theorem sinePi_eq_sub_mul_sineSlope (m : ℤ) (z : ℂ) :
    sinePi z = (z - (m : ℂ)) * sineSlope m z := by
  have h := sub_smul_dslope sinePi (m : ℂ) z
  simpa [sineSlope, smul_eq_mul] using h.symm

private theorem sinePi_differentiable : Differentiable ℂ sinePi := by
  intro z
  simpa [sinePi] using
    (((hasDerivAt_id z).const_mul (Real.pi : ℂ)).csin.differentiableAt)

private theorem sineSlope_differentiable (m : ℤ) :
    Differentiable ℂ (sineSlope m) := by
  rw [← differentiableOn_univ]
  exact (Complex.differentiableOn_dslope
      (f := sinePi) (s := Set.univ) (c := (m : ℂ)) univ_mem).2
    sinePi_differentiable.differentiableOn

private theorem sineSlope_at_int_ne_zero (m : ℤ) :
    sineSlope m (m : ℂ) ≠ 0 := by
  have hderiv :
      deriv sinePi (m : ℂ) =
        (Real.pi : ℂ) *
          Complex.cos ((Real.pi : ℂ) * (m : ℂ)) := by
    have hinner :
        HasDerivAt (fun z : ℂ => (Real.pi : ℂ) * z)
          (Real.pi : ℂ) (m : ℂ) := by
      simpa only [id_eq, mul_comm, mul_one] using
        ((hasDerivAt_id (m : ℂ)).const_mul (Real.pi : ℂ))
    have h := hinner.csin.deriv
    rw [show sinePi =
      (fun z : ℂ => Complex.sin ((Real.pi : ℂ) * z)) by rfl]
    simpa only [mul_comm] using h
  have harg :
      (Real.pi : ℂ) * (m : ℂ) =
        (((m : ℝ) * Real.pi : ℝ) : ℂ) := by
    push_cast
    ring
  have hcos :
      Complex.cos ((Real.pi : ℂ) * (m : ℂ)) =
        ((((-1 : ℝ) ^ m : ℝ)) : ℂ) := by
    rw [harg, ← Complex.ofReal_cos, Real.cos_int_mul_pi]
  have hpow : ((-1 : ℝ) ^ m) ≠ 0 := by
    intro hzero
    have habs := Real.abs_cos_int_mul_pi m
    rw [Real.cos_int_mul_pi, hzero, abs_zero] at habs
    norm_num at habs
  have hcos0 : Complex.cos ((Real.pi : ℂ) * (m : ℂ)) ≠ 0 := by
    rw [hcos]
    exact Complex.ofReal_ne_zero.mpr hpow
  rw [sineSlope, dslope_same, hderiv]
  exact mul_ne_zero (Complex.ofReal_ne_zero.mpr Real.pi_ne_zero) hcos0

private def halfIntegerStrip (m : ℤ) : Set ℂ :=
  {z | (m : ℝ) - 1 / 2 ≤ z.re ∧ z.re ≤ (m : ℝ) + 1 / 2}

private theorem sineSlope_ne_zero_on_strip
    (m : ℤ) {z : ℂ} (hz : z ∈ halfIntegerStrip m) :
    sineSlope m z ≠ 0 := by
  by_cases hzm : z = (m : ℂ)
  · subst z
    exact sineSlope_at_int_ne_zero m
  intro hslope
  have hsin : sinePi z = 0 := by
    rw [sinePi_eq_sub_mul_sineSlope m z, hslope, mul_zero]
  rcases Complex.sin_eq_zero_iff.mp (by simpa [sinePi] using hsin) with ⟨k, hk⟩
  have hpi0 : (Real.pi : ℂ) ≠ 0 :=
    Complex.ofReal_ne_zero.mpr Real.pi_ne_zero
  have hzk : z = (k : ℂ) := by
    calc
      z = (Real.pi : ℂ)⁻¹ * ((Real.pi : ℂ) * z) := by
        field_simp
      _ = (Real.pi : ℂ)⁻¹ * ((k : ℂ) * (Real.pi : ℂ)) := by
        rw [hk]
      _ = (k : ℂ) := by
        field_simp
  subst z
  have hl : (m : ℝ) - 1 / 2 ≤ (k : ℝ) := by
    simpa [halfIntegerStrip] using hz.1
  have hu : (k : ℝ) ≤ (m : ℝ) + 1 / 2 := by
    simpa [halfIntegerStrip] using hz.2
  have hlowR : (-1 : ℝ) < ((k - m : ℤ) : ℝ) := by
    push_cast
    linarith
  have huppR : ((k - m : ℤ) : ℝ) < (1 : ℝ) := by
    push_cast
    linarith
  have hlowZ : (-1 : ℤ) < k - m := by
    exact_mod_cast hlowR
  have huppZ : k - m < (1 : ℤ) := by
    exact_mod_cast huppR
  have hkm : k = m := by omega
  apply hzm
  simpa [hkm]

private def barnesRaw (κ : ℂ) (A : ℂ → ℂ) (z : ℂ) : ℂ :=
  κ * A z / sinePi z ^ 2

private def barnesExtension
    (κ : ℂ) (m : ℤ) (P : ℂ → ℂ) (z : ℂ) : ℂ :=
  κ * P z / sineSlope m z ^ 2

private theorem barnesRaw_eq_extension_of_mem_strip_of_ne
    {κ : ℂ} {A P : ℂ → ℂ} (m : ℤ) {z : ℂ}
    (hfactor : ∀ w : ℂ,
      A w = (w - (m : ℂ)) ^ 2 * P w)
    (hz : z ∈ halfIntegerStrip m) (hzm : z ≠ (m : ℂ)) :
    barnesRaw κ A z = barnesExtension κ m P z := by
  have hslope := sineSlope_ne_zero_on_strip m hz
  rw [barnesRaw, barnesExtension, hfactor,
    sinePi_eq_sub_mul_sineSlope m z]
  field_simp [sub_ne_zero.mpr hzm, hslope]

private theorem barnesExtension_differentiableOn
    {κ : ℂ} {P : ℂ → ℂ} (m : ℤ)
    (hP : DifferentiableOn ℂ P (halfIntegerStrip m)) :
    DifferentiableOn ℂ (barnesExtension κ m P) (halfIntegerStrip m) := by
  have hconst : DifferentiableOn ℂ (fun _ : ℂ => κ) (halfIntegerStrip m) :=
    differentiableOn_const κ
  have hnum : DifferentiableOn ℂ (fun z => κ * P z) (halfIntegerStrip m) :=
    hconst.mul hP
  have hden : DifferentiableOn ℂ (fun z => sineSlope m z ^ 2)
      (halfIntegerStrip m) :=
    (sineSlope_differentiable m).differentiableOn.pow 2
  simpa [barnesExtension] using
    hnum.div hden (fun z hz => pow_ne_zero 2 (sineSlope_ne_zero_on_strip m hz))

/-- One-strip Barnes contour shift across an integer whose apparent pole is
removed by the supplied quadratic factorization. -/
theorem barnes_one_strip_shift
    {κ : ℂ} {A P : ℂ → ℂ} (m : ℤ)
    (hfactor : ∀ z : ℂ,
      A z = (z - (m : ℂ)) ^ 2 * P z)
    (hP : DifferentiableOn ℂ P (halfIntegerStrip m))
    (hleft : MeasureTheory.Integrable
      (fun y : ℝ => barnesRaw κ A
        (verticalPoint ((m : ℝ) - 1 / 2) y)))
    (hright : MeasureTheory.Integrable
      (fun y : ℝ => barnesRaw κ A
        (verticalPoint ((m : ℝ) + 1 / 2) y)))
    (htop : Tendsto
      (fun T : ℝ =>
        ∫ x in ((m : ℝ) - 1 / 2)..((m : ℝ) + 1 / 2),
          barnesExtension κ m P
            ((x : ℂ) + (T : ℂ) * Complex.I))
      atTop (𝓝 0))
    (hbottom : Tendsto
      (fun T : ℝ =>
        ∫ x in ((m : ℝ) - 1 / 2)..((m : ℝ) + 1 / 2),
          barnesExtension κ m P
            ((x : ℂ) - (T : ℂ) * Complex.I))
      atTop (𝓝 0)) :
    (∫ y : ℝ, barnesRaw κ A
      (verticalPoint ((m : ℝ) - 1 / 2) y)) =
    ∫ y : ℝ, barnesRaw κ A
      (verticalPoint ((m : ℝ) + 1 / 2) y) := by
  let a : ℝ := (m : ℝ) - 1 / 2
  let b : ℝ := (m : ℝ) + 1 / 2
  have hab : a ≤ b := by dsimp [a, b]; linarith
  have hstrip : closedVerticalStrip a b = halfIntegerStrip m := by
    ext z
    simp [closedVerticalStrip, halfIntegerStrip, a, b]
  have hleftEq :
      (fun y : ℝ => barnesExtension κ m P (verticalPoint a y)) =
        fun y : ℝ => barnesRaw κ A (verticalPoint a y) := by
    funext y
    symm
    apply barnesRaw_eq_extension_of_mem_strip_of_ne m hfactor
    · simp [halfIntegerStrip, verticalPoint, a] <;> linarith
    · intro h
      have hre := congrArg Complex.re h
      simp [verticalPoint, a] at hre
  have hrightEq :
      (fun y : ℝ => barnesExtension κ m P (verticalPoint b y)) =
        fun y : ℝ => barnesRaw κ A (verticalPoint b y) := by
    funext y
    symm
    apply barnesRaw_eq_extension_of_mem_strip_of_ne m hfactor
    · simp [halfIntegerStrip, verticalPoint, b] <;> linarith
    · intro h
      have hre := congrArg Complex.re h
      simp [verticalPoint, b] at hre
  have hleftExt : MeasureTheory.Integrable
      (fun y : ℝ => barnesExtension κ m P (verticalPoint a y)) := by
    rw [hleftEq]
    simpa [a] using hleft
  have hrightExt : MeasureTheory.Integrable
      (fun y : ℝ => barnesExtension κ m P (verticalPoint b y)) := by
    rw [hrightEq]
    simpa [b] using hright
  have hshift := verticalIntegral_eq_of_horizontal_tendsto
    (F := barnesExtension κ m P) hab
    (by rw [hstrip]; exact barnesExtension_differentiableOn m hP)
    hleftExt hrightExt
    (by simpa [a, b] using htop)
    (by simpa [a, b] using hbottom)
  rw [hleftEq, hrightEq] at hshift
  simpa [a, b] using hshift

#print axioms verticalIntegral_eq_of_horizontal_tendsto
#print axioms sineSlope_ne_zero_on_strip
#print axioms barnesRaw_eq_extension_of_mem_strip_of_ne
#print axioms barnes_one_strip_shift

end RamanujanChallenge.P27.Q6382
```

# Exact Mathlib signatures and orientation

## Rectangle Cauchy theorem

The code uses:

```lean
Complex.integral_boundary_rect_eq_zero_of_differentiableOn
```

with opposite corners

```lean
z = a - T * I
w = b + T * I.
```

After simplification, Mathlib's orientation is exactly

```text
bottom - top + I * right - I * left = 0.
```

Thus

```text
left - right = I * (top - bottom).
```

The code does this in two explicit steps: first derives `hIR`, then uses `I^2=-1`. There is no informal orientation convention left to infer.

The theorem compares parameterized vertical integrals

```text
∫ y, F(a+i y),    ∫ y, F(b+i y).
```

If the project's contour notation includes `dt`, then each side is multiplied by the same `I`, since `dt=I dy`; the equality is unchanged.

## Improper-limit bridge

The exact theorem used is:

```lean
MeasureTheory.intervalIntegral_tendsto_integral
```

instantiated as

```lean
MeasureTheory.intervalIntegral_tendsto_integral
  hleft tendsto_neg_atTop_atBot tendsto_id
```

and similarly on the right. This converts `∫_{-T}^T` to the whole real-line Bochner integral.

## Removable sine quotient

The extension is

```lean
sineSlope m = dslope sinePi (m : ℂ),
sinePi z = sin (π z).
```

The identity

```lean
sinePi z = (z-m) * sineSlope m z
```

comes directly from

```lean
sub_smul_dslope.
```

Holomorphy is obtained from the exact equivalence

```lean
Complex.differentiableOn_dslope
```

on `Set.univ`. No piecewise continuity or hand-built `Function.update` proof is needed.

## Nonvanishing and casts

At the center, the extension equals the derivative:

```lean
dslopesame:
  sineSlope m m = deriv sinePi m
              = π * cos(πm).
```

The proof uses:

```lean
Real.cos_int_mul_pi
Real.abs_cos_int_mul_pi
Complex.ofReal_cos
Complex.ofReal_ne_zero
```

so no unproved claim about `(-1)^m` is introduced.

Away from the center, if `sineSlope m z=0`, then `sin(πz)=0`. The exact zero classification

```lean
Complex.sin_eq_zero_iff
```

gives `z=k` for an integer `k`. Strip membership yields

```text
-1 < k-m < 1.
```

The real inequalities are transported to `ℤ` with `exact_mod_cast`, and `omega` proves `k=m`, contradicting `z≠m`. This resolves the full cast/parity block without an appeal to geometric intuition.

# Why a triple zero is exactly enough

The generic theorem asks for

```lean
A z = (z-m)^2 * P z.
```

In the Barnes application,

```text
A = S - S'/2.
```

If the Zeilberger certificate has a triple zero

```text
S(z) = (z-m)^3 Q(z),
```

then differentiating gives

```text
S'(z) = 3(z-m)^2 Q(z) + (z-m)^3 Q'(z),
```

hence

```text
S(z) - S'(z)/2
  = (z-m)^2 *
      ((z-m) Q(z) - 3 Q(z)/2 - (z-m) Q'(z)/2).
```

Therefore a canonical specialization is

```lean
P z :=
  (z - (m : ℂ)) * Q z - 3 * Q z / 2
    - (z - (m : ℂ)) * deriv Q z / 2
```

provided the derivative identity for `S` and `Q` is established on the strip. For a finite rational expression, it is usually shorter to prove the denominator-cleared identity for `A=(z-m)^2P` directly with `field_simp` and `ring` than to formalize an abstract order-of-zero API.

# Remaining specialized obligations

The generic contour machinery is complete. To apply it to a concrete Zudilin/Zeilberger term, the remaining work is exactly the following.

Let

```text
a = m - 1/2,
b = m + 1/2,
A(z) = S(z) - S'(z)/2,
κ = π²                    -- current project kernel
```

and choose `P` with `A=(z-m)^2P`.

## 1. Algebraic factor certificate

Prove globally, or at least on the strip:

```lean
hfactor : ∀ z : ℂ,
  A z = (z - (m : ℂ)) ^ 2 * P z
```

For the current rational functions, this should be a denominator-cleared finite-product/rational identity. Do not define `S'` via a logarithmic derivative at zeros; use the actual `deriv`, as `Problem27Barnes.lean` already does for `zudilinBarnesFPrime27`.

## 2. Holomorphy of `P` on the strip

Prove:

```lean
hP : DifferentiableOn ℂ P (halfIntegerStrip m)
```

This reduces to showing every rational denominator factor is nonzero throughout

```text
m-1/2 ≤ Re z ≤ m+1/2.
```

For the positive-integer strips crossed by the Apéry/Zudilin numerator zeros, the poles of the Barnes rational functions lie to the left, so these are elementary real-part inequalities of the same form already used in `Problem27Barnes.lean`.

## 3. Vertical integrability on both boundary lines

Prove exactly:

```lean
MeasureTheory.Integrable
  (fun y : ℝ => barnesRaw κ A
    (((m : ℝ) - 1/2 : ℝ) + y * Complex.I))

MeasureTheory.Integrable
  (fun y : ℝ => barnesRaw κ A
    (((m : ℝ) + 1/2 : ℝ) + y * Complex.I))
```

On either half-integer line,

```text
|sin(π(m±1/2+i y))| = cosh(π y),
```

so it suffices to dominate the rational-polynomial numerator by a fixed polynomial in `|y|` and use the exponential decay of `cosh(πy)^{-2}`. This is simpler than the uniform-in-`n` midpoint bound already formalized in `Problem27Barnes.lean`, because here `n` and the strip are fixed during one contour shift.

## 4. Top horizontal edge

Prove precisely:

```lean
Tendsto
  (fun T : ℝ =>
    ∫ x in ((m : ℝ) - 1/2)..((m : ℝ) + 1/2),
      barnesExtension κ m P
        ((x : ℂ) + (T : ℂ) * Complex.I))
  atTop (𝓝 0)
```

## 5. Bottom horizontal edge

Prove precisely:

```lean
Tendsto
  (fun T : ℝ =>
    ∫ x in ((m : ℝ) - 1/2)..((m : ℝ) + 1/2),
      barnesExtension κ m P
        ((x : ℂ) - (T : ℂ) * Complex.I))
  atTop (𝓝 0)
```

If all rational coefficients are real, conjugation identifies the norms of the two edges, so one growth lemma can discharge both.

# A reusable horizontal-decay lemma to add next

The most useful next generic theorem is not another contour theorem, but the following bounded-edge wrapper.

```lean
theorem horizontalIntegral_tendsto_zero_of_bound
    {F : ℂ → ℂ} {a b C : ℝ} {d : ℕ}
    (hab : a ≤ b) (hC : 0 ≤ C)
    (hbound : ∀ᶠ T : ℝ in atTop,
      ∀ x ∈ Set.Icc a b,
        ‖F ((x : ℂ) + (T : ℂ) * Complex.I)‖ ≤
          C * (1 + T) ^ d * Real.exp (-2 * Real.pi * T)) :
    Tendsto
      (fun T : ℝ => ∫ x in a..b,
        F ((x : ℂ) + (T : ℂ) * Complex.I))
      atTop (𝓝 0) := by
  -- Use `intervalIntegral.norm_integral_le_of_norm_le_const`,
  -- then squeeze by `C * |b-a| * (1+T)^d * exp(-2πT)`.
  -- Finish with `Real.tendsto_pow_mul_exp_neg_atTop_nhds_zero`
  -- after the positive linear change of variable `2πT`.
  ...
```

For a one-strip interval, `|b-a|=1`. The pointwise estimate needed for a rational-polynomial Barnes factor is of the form

```text
‖A(x±iT)‖ ≤ C₁ (1+T)^d,
```

uniformly for `x∈[m-1/2,m+1/2]`, together with

```text
1 / ‖sin(π(x±iT))‖² ≤ 16 exp(-2πT)
```

eventually in `T`.

The sine estimate follows from

```text
sin(π(x+iT))
 = sin(πx) cosh(πT) + i cos(πx) sinh(πT),
```

hence

```text
‖sin(π(x+iT))‖²
 = sin²(πx) + sinh²(πT)
 ≥ sinh²(πT).
```

For sufficiently large `T`,

```text
sinh(πT) ≥ exp(πT)/4,
```

which gives the displayed factor `16`. The exact constant is irrelevant.

On a horizontal edge `T>0`, no point equals `m`, so `barnesRaw=barnesExtension`; it is usually easier to prove the growth estimate using the raw `sin^{-2}` expression and then rewrite to the extension.

# Reuse in the four-term Zeilberger recurrence

For a four-term telescoper, keep the contour theorem completely separate from the recurrence algebra.

For each recurrence term `j=0,1,2,3`, package:

```text
κ_j,
A_j(t) = S_j(t) - S_j'(t)/2,
P_j,
m_j,
A_j(t) = (t-m_j)^2 P_j(t),
holomorphy of P_j on the crossed strip,
vertical integrability,
horizontal decay.
```

Then:

1. invoke `barnes_one_strip_shift` once for each crossed strip;
2. iterate it if a term must move across several consecutive integers;
3. align all four terms on one common half-integer line;
4. use linearity of the real-line integral to insert the four rational recurrence coefficients;
5. apply the pointwise Zeilberger identity on the common line.

For the Zudilin rational functions, the numerator contains a product of cubic factors. Every crossed positive integer is therefore a triple zero of the certificate, so each one-strip invocation is residue-free. The theorem never computes, assumes, or cancels a residue; it replaces the apparent pole with an explicitly constructed holomorphic function before Cauchy is applied.

A practical production API is:

```lean
structure BarnesStripCertificate where
  κ : ℂ
  A P : ℂ → ℂ
  m : ℤ
  factor : ∀ z, A z = (z - (m : ℂ)) ^ 2 * P z
  diffP : DifferentiableOn ℂ P (halfIntegerStrip m)
  integrable_left : ...
  integrable_right : ...
  top_zero : ...
  bottom_zero : ...
```

with

```lean
theorem BarnesStripCertificate.shift
    (c : BarnesStripCertificate) : ... :=
  barnes_one_strip_shift c.m c.factor c.diffP
    c.integrable_left c.integrable_right c.top_zero c.bottom_zero
```

This keeps the 4-term recurrence proof readable and makes the same machinery reusable for every shifted `R_k` contour.

# Recommended repository placement

Add the generic material to:

```text
lean/RamanujanChallenge/Problem27BarnesShift.lean
```

with imports:

```lean
import RamanujanChallenge.Problem27Barnes
import Mathlib.Analysis.Complex.CauchyIntegral
import Mathlib.Analysis.Complex.RemovableSingularity
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Complex
import Mathlib.MeasureTheory.Integral.IntegralEqImproper
```

Do not place it in `SUBMIT/`; the repository's maintained policy says work belongs under `lean/` until the full normalization chain is green, sorry-free, and axiom-audited.

The temporary audit PR was used only to run the pinned-toolchain compile. It is not intended for merge.