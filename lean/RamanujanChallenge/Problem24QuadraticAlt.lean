/-
  Problem 2.4, Q⁻: the outer-alternating quadratic level-two Euler sum.

  Target theorem (the missing hypothesis `hAlternatingQuadratic` of
  `Problem24.of_standard_euler_and_wz`):

      HasSum alternatingQuadraticEulerTerm24 alternatingQuadraticEulerValue24

  with
      alternatingQuadraticEulerTerm24 n = (-1)^(n+1) * (P(n+1)^2 - H(n+1)^(2)) / (n+1)^2
      P_n = H_n + 2 * signedHarmonic24 n,   H_n^(2) = harmonicSquare24 n

  and
      alternatingQuadraticEulerValue24 =
        -22*Li4(1/2) - (11/12)*log^4 2 - (13/2)*log^2 2 * zeta2
        - (7/4)*log 2 * zeta3 + (67/10)*zeta2^2.

  Certificate (harvested from chatgpt-answers/Q6047-7668a207.md,
  "Alternating quadratic level-two sum: a six-integral real certificate"):

  * closed generating function Q(x) = sum_n S_n x^n,  S_n = P_n^2 - H_n^(2):
        Q(x) = 2 J(x)/(1-x),
        J(x) = (1/2) a^2 + b^2 + 2ab + 2L b + Z2 - L^2 - 2 Li2((1+x)/2) + Li2(x^2),
        a = log(1-x), b = log(1+x), L = log 2.
  * coefficient integration:  Aquad = integral_0^1 (-log x)/x * Q(-x) dx.
  * one integration by parts with W(x) = Z2 - 2 Li2(u) - log(u)^2, u = x/(1+x);
  * Mobius map t = 2x/(1+x) reduces to six integrals
        I10..I22 = integral_0^1 W0(t) * H1(t)/H2(t) over t, (1-t), (2-t),
        W0(t) = Z2 - 2 Li2(t/2) - log(t/2)^2,  H1(t) = -log(1-t), H2(t) = -log(1-t/2);
  * the bridge:  Aquad = -2*I10 - 2*I11 + 2*I12 + 4*I20 + 6*I21 - 5*I22;
  * six endpoint evaluations in terms of the already-proved linear sums
    Tplus = cubicLinearEulerValue24,  Tminus = alternatingCubicLinearEulerValue24;
  * final bridge:  Aquad = 11*Tminus + 2*Tplus - 12 L^2 Z2 - 14 L Z3 + (41/10) Z2^2.

  All names below are public definitions from Problem24 / Problem24Euler;
  nothing here claims the endpoint evaluations yet — they are the open work.
-/
import RamanujanChallenge.Problem24
import RamanujanChallenge.Problem24Euler
import RamanujanChallenge.Dilogarithm
import RamanujanChallenge.Problem26WeightThree
import Mathlib.Analysis.Calculus.LHopital

noncomputable section

open Filter Set Topology
open scoped Interval

set_option maxHeartbeats 800000

namespace RamanujanChallenge.P24QuadAlt

/-! ## Layer A: the quadratic remainder sequence `S` (Q6047 §1-2)

`S n = P(n+1)^2 - H(n+1)^(2)` is the unnormalized quadratic Euler summand
(`quadraticEulerTerm24 n = S n / (n+1)^2`).  The increment identity
(Q6047 (2.2)) is the algebraic core of the generating-function proof. -/

/-- `S n = P(n+1)² − H(n+1)²`, with `P = parityRemainder24`,
`H² = harmonicSquare24`. -/
def quadAltS (n : ℕ) : ℝ :=
  parityRemainder24 (n + 1) ^ 2 - harmonicSquare24 (n + 1)

/-- `quadraticEulerTerm24 n = S n / (n+1)²`. -/
theorem quadAltS_eq_quadraticEulerTerm (n : ℕ) :
    quadAltS n / (n + 1 : ℝ) ^ 2 = quadraticEulerTerm24 n := by
  unfold quadAltS quadraticEulerTerm24
  ring

/-- `((-1)^m)^2 = 1`. -/
theorem quadAlt_neg_one_pow_sq (m : ℕ) :
    ((-1 : ℝ) ^ m) ^ 2 = 1 := by
  rw [sq, ← pow_add]
  rw [show m + m = 2 * m by omega]
  rw [pow_mul]
  norm_num

/-- `(−1)^(m·2) = 1` (even power). -/
theorem quadAlt_neg_one_pow_even (m : ℕ) : (-1 : ℝ) ^ (m * 2) = 1 := by
  rw [mul_comm, pow_mul]
  norm_num

/-- `S(n) − S(n−1) = 2c·P(n)/(n+1) + (c²−1)/(n+1)²` with
`c = 1 + 2(−1)^{n+1}` (Q6047 (2.2), index shifted by one). -/
theorem quadAltS_succ_sub (n : ℕ) :
    quadAltS (n + 1) - quadAltS n =
      2 * (1 + 2 * (-1 : ℝ) ^ (n + 2)) * parityRemainder24 (n + 1) /
          (n + 2 : ℝ) +
        ((1 + 2 * (-1 : ℝ) ^ (n + 2)) ^ 2 - 1) /
          (n + 2 : ℝ) ^ 2 := by
  let t : ℝ := parityRemainder24 (n + 1)
  let s : ℝ := (-1 : ℝ) ^ (n + 2)
  let c : ℝ := 1 + 2 * s
  let d : ℝ := (n + 2 : ℝ)
  have hsq : s ^ 2 = 1 := by
    simpa [s] using quadAlt_neg_one_pow_sq (n + 2)
  have hstepP : parityRemainder24 ((n + 1) + 1) = t + c / d := by
    dsimp [t, c, d, s]
    rw [parityRemainder24_succ (n + 1)]
    have hpow : (-1 : ℝ) ^ ((n + 1) + 1) = (-1 : ℝ) ^ (n + 2) := by rfl
    rw [hpow]
    push_cast
    ring
  have hstepH : harmonicSquare24 ((n + 1) + 1) =
      harmonicSquare24 (n + 1) + 1 / d ^ 2 := by
    dsimp [d]
    rw [harmonicSquare24_succ (n + 1)]
    congr 1
    push_cast
    ring
  calc
    quadAltS (n + 1) - quadAltS n
        = (t + c / d) ^ 2 - (harmonicSquare24 (n + 1) + 1 / d ^ 2) -
            (t ^ 2 - harmonicSquare24 (n + 1)) := by
          unfold quadAltS
          rw [hstepP, hstepH]
    _ = 2 * t * c / d + (c ^ 2 - 1) / d ^ 2 := by ring
    _ = 2 * c * t / d + (c ^ 2 - 1) / d ^ 2 := by ring
    _ = 2 * (1 + 2 * (-1 : ℝ) ^ (n + 2)) * parityRemainder24 (n + 1) /
            (n + 2 : ℝ) +
          ((1 + 2 * (-1 : ℝ) ^ (n + 2)) ^ 2 - 1) / (n + 2 : ℝ) ^ 2 := by
          dsimp [t, c, s, d]

/-- `(c²−1) = 4(1+(−1)^{n+2})` (Q6047 (2.3)). -/
theorem quadAltC_sq_sub_one (n : ℕ) :
    (1 + 2 * (-1 : ℝ) ^ (n + 2)) ^ 2 - 1 =
      4 * (1 + (-1 : ℝ) ^ (n + 2)) := by
  have hsq : ((-1 : ℝ) ^ (n + 2)) ^ 2 = 1 :=
    quadAlt_neg_one_pow_sq (n + 2)
  calc
    (1 + 2 * (-1 : ℝ) ^ (n + 2)) ^ 2 - 1
        = 4 * (-1 : ℝ) ^ (n + 2) + 4 * ((-1 : ℝ) ^ (n + 2)) ^ 2 := by
          ring
    _ = 4 * (-1 : ℝ) ^ (n + 2) + 4 := by
      rw [hsq]
      ring
    _ = 4 * (1 + (-1 : ℝ) ^ (n + 2)) := by ring

/-! ## Layer B: closed generating functions (Q6047 §2)

With `a(x) = log(1-x)`, `b(x) = log(1+x)`, `L = log 2`, `Z2 = π²/6`:

    F(x) = Σ P_{n+1} x^n          = -(a + 2b)/(x(1-x))            (2.1)
    M(x) = Σ c_n P_{n-1} x^n / n  = a²/2 + b² + 2ab + 2Lb + Z2 - L²
                                     - 2 Li2((1+x)/2)              (2.7)
    J(x) = M(x) + Li2(x²)                                          (2.11)
    Q(x) = Σ S_{n+1} x^{n+1}      = 2 J(x)/(1-x)                   (2.10)

`M'(x) = F(x) - 2F(-x)` (2.6) is the derivative certificate. -/

/-- `Mclosed` (2.7): the normalized primitive of `F(x) − 2F(−x)`. -/
def quadAltMclosed (x : ℝ) : ℝ :=
  Real.log (1 - x) ^ 2 / 2 + Real.log (1 + x) ^ 2 +
    2 * Real.log (1 - x) * Real.log (1 + x) +
    2 * Real.log 2 * Real.log (1 + x) +
    Real.pi ^ 2 / 6 - Real.log 2 ^ 2 -
    2 * dilog ((1 + x) / 2)

/-- `Jclosed` (2.11): `M + Li2(x²)`. -/
def quadAltJclosed (x : ℝ) : ℝ :=
  quadAltMclosed x + dilog (x ^ 2)

/-- `Qclosed` (2.10): `2 J(x)/(1-x)`. -/
def quadAltQclosed (x : ℝ) : ℝ :=
  2 * quadAltJclosed x / (1 - x)

/-! ## Layer B: derivative certificate for `M` (Q6047 §2.2)

`F(x) = Σ_{n≥1} P_n xⁿ = −(a + 2b)/(1−x)` with `a = log(1−x)`,
`b = log(1+x)` (Q6047 (2.1)); note this is `x` times the value of
`parityRemainder24_generating_hasSum`.  `M'(x) = F(x) − 2F(−x)` (2.6). -/

/-- `F(x) = −(a + 2b)/(1−x)`, `a = log(1−x)`, `b = log(1+x)` (Q6047 (2.1)). -/
def quadAltFclosed (x : ℝ) : ℝ :=
  -(Real.log (1 - x) + 2 * Real.log (1 + x)) / (1 - x)

/-- `Mclosed` is the normalized primitive: `M'(x) = F(x) − 2F(−x)`
(Q6047 (2.6)–(2.8)). -/
theorem quadAltMclosed_hasDerivAt {x : ℝ} (hx0 : 0 < x) (hx1 : x < 1) :
    HasDerivAt quadAltMclosed (quadAltFclosed x - 2 * quadAltFclosed (-x)) x := by
  have h1mx : 0 < 1 - x := by nlinarith
  have h1px : 0 < 1 + x := by nlinarith
  have hmid0 : 0 < (1 + x) / 2 := by nlinarith
  have hmid1 : (1 + x) / 2 < 1 := by nlinarith
  have ha : HasDerivAt (fun y : ℝ => Real.log (1 - y)) (-1 / (1 - x)) x := by
    have hc : HasDerivAt (fun y : ℝ => 1 - y) (-1) x := by
      convert (hasDerivAt_const x 1).sub (hasDerivAt_id x) using 1
      simp
    have hlog := Real.hasDerivAt_log (ne_of_gt h1mx)
    have hD : HasDerivAt (fun y : ℝ => Real.log (1 - y)) ((1 - x)⁻¹ * (-1)) x :=
      HasDerivAt.comp (h := fun y : ℝ => 1 - y) x hlog hc
    convert hD using 1
    field_simp
  have hb : HasDerivAt (fun y : ℝ => Real.log (1 + y)) (1 / (1 + x)) x := by
    have hc : HasDerivAt (fun y : ℝ => 1 + y) 1 x := by
      convert (hasDerivAt_const x 1).add (hasDerivAt_id x) using 1
      simp
    have hlog := Real.hasDerivAt_log (ne_of_gt h1px)
    have hD : HasDerivAt (fun y : ℝ => Real.log (1 + y)) ((1 + x)⁻¹ * 1) x :=
      HasDerivAt.comp (h := fun y : ℝ => 1 + y) x hlog hc
    convert hD using 1
    field_simp
  have hd : HasDerivAt (fun y : ℝ => (1 + y) / 2) (1 / 2) x := by
    have hc : HasDerivAt (fun y : ℝ => 1 + y) 1 x := by
      convert (hasDerivAt_const x 1).add (hasDerivAt_id x) using 1
      simp
    convert hc.div_const 2 using 1
  have hD : HasDerivAt (fun y : ℝ => dilog ((1 + y) / 2))
      (-Real.log ((1 - x) / 2) / (1 + x)) x := by
    have hmid := dilog_hasDerivAt hmid0 hmid1
    have hD2 : HasDerivAt (fun y : ℝ => dilog ((1 + y) / 2))
        (-(Real.log (1 - (1 + x) / 2)) / ((1 + x) / 2) * (1 / 2)) x :=
      HasDerivAt.comp (h := fun y : ℝ => (1 + y) / 2) x hmid hd
    convert hD2 using 1
    field_simp
    congr 1
    ring
  -- assemble the derivative of Mclosed
  have hconst : HasDerivAt (fun y : ℝ => Real.pi ^ 2 / 6 - Real.log 2 ^ 2) 0 x := by
    exact hasDerivAt_const x (Real.pi ^ 2 / 6 - Real.log 2 ^ 2)
  -- a²/2
  have h1 : HasDerivAt (fun y : ℝ => Real.log (1 - y) ^ 2 / 2)
      (-Real.log (1 - x) / (1 - x)) x := by
    convert (ha.pow 2).div_const 2 using 1
    all_goals norm_num
    all_goals field_simp [h1mx.ne', h1px.ne']
    all_goals ring_nf
  -- b²
  have h2 : HasDerivAt (fun y : ℝ => Real.log (1 + y) ^ 2)
      (2 * Real.log (1 + x) / (1 + x)) x := by
    convert hb.pow 2 using 1
    all_goals norm_num
    all_goals field_simp [h1mx.ne', h1px.ne']
    all_goals ring_nf
  -- 2ab
  have h3 : HasDerivAt (fun y : ℝ => 2 * Real.log (1 - y) * Real.log (1 + y))
      (2 * (-Real.log (1 + x) / (1 - x) + Real.log (1 - x) / (1 + x))) x := by
    have hm := ha.mul hb
    convert hm.const_mul 2 using 1
    all_goals norm_num
    all_goals field_simp [h1mx.ne', h1px.ne']
    all_goals ring_nf
  -- 2Lb
  have h4 : HasDerivAt (fun y : ℝ => 2 * Real.log 2 * Real.log (1 + y))
      (2 * Real.log 2 / (1 + x)) x := by
    convert hb.const_mul (2 * Real.log 2) using 1
    all_goals norm_num
    all_goals field_simp [h1mx.ne', h1px.ne']
    all_goals ring_nf
  -- −2·dilog((1+y)/2)
  have h5 : HasDerivAt (fun y : ℝ => -2 * dilog ((1 + y) / 2))
      (2 * Real.log ((1 - x) / 2) / (1 + x)) x := by
    convert hD.const_mul (-2) using 1
    all_goals norm_num
    all_goals field_simp [h1mx.ne', h1px.ne']
    all_goals ring_nf
  have hM : HasDerivAt (fun y : ℝ => quadAltMclosed y)
      (-Real.log (1 - x) / (1 - x) +
        2 * Real.log (1 + x) / (1 + x) +
        2 * (-Real.log (1 + x) / (1 - x) + Real.log (1 - x) / (1 + x)) +
        2 * Real.log 2 / (1 + x) +
        0 +
        2 * Real.log ((1 - x) / 2) / (1 + x)) x := by
    unfold quadAltMclosed
    convert ((((h1.add h2).add h3).add h4).add hconst).add h5 using 2
    all_goals simp
    all_goals field_simp [h1mx.ne', h1px.ne']
    all_goals ring_nf
  convert hM using 1
  unfold quadAltFclosed
  norm_num
  rw [show Real.log ((1 - x) / 2) = Real.log (1 - x) - Real.log 2 by
    rw [Real.log_div]
    · exact h1mx.ne'
    · norm_num]
  ring_nf

/-! ## Layer D: integration by parts + Möbius normalization (Q6047 §4-5)

`W(x) = Z2 − 2 Li2(u) − log(u)²` with `u = x/(1+x)` (Landen form, (4.7));
after the Möbius map `t = 2x/(1+x)` this becomes `W0(t)`. -/

/-- `Dminus(x) = J'(−x)` (Q6047 (4.8)): the derivative of `J` at `−x`,
with `a = log(1−x)`, `b = log(1+x)`. -/
def quadAltDminus (x : ℝ) : ℝ :=
  -(Real.log (1 + x) + 2 * Real.log (1 - x)) / (1 + x) +
    2 * (Real.log (1 - x) + 2 * Real.log (1 + x)) / (1 - x) +
    2 * (Real.log (1 - x) + Real.log (1 + x)) / x


/-- `W0(t) = Z2 − 2 Li2(t/2) − log(t/2)²` (Q6047 (5.2)). -/
def W0 (t : ℝ) : ℝ :=
  Real.pi ^ 2 / 6 - 2 * dilog (t / 2) - Real.log (t / 2) ^ 2

/-- `H1(t) = -log(1-t)`. -/
def H1 (t : ℝ) : ℝ := -Real.log (1 - t)

/-- `H2(t) = -log(1-t/2)`. -/
def H2 (t : ℝ) : ℝ := -Real.log (1 - t / 2)

/-- The Möbius identity (Q6047 (5.4)): `Dminus(t/(2−t))·2/(2−t)²`
equals the `H1`/`H2` combination. -/
theorem quadAltMobius_identity {t : ℝ} (ht0 : 0 < t) (ht1 : t < 1) :
    quadAltDminus (t / (2 - t)) * 2 / (2 - t) ^ 2 =
      H1 t * (-2 / t - 2 / (1 - t) + 2 / (2 - t)) +
      H2 t * (4 / t + 6 / (1 - t) - 5 / (2 - t)) := by
  have h2mt : 2 - t ≠ 0 := by linarith
  have h1mt : 1 - t ≠ 0 := by linarith
  have htm1 : t - 1 ≠ 0 := by linarith
  unfold quadAltDminus H1 H2
  have h1mx : 1 - t / (2 - t) = 2 * (1 - t) / (2 - t) := by
    field_simp [h2mt]
    ring
  have h1px : 1 + t / (2 - t) = 2 / (2 - t) := by
    field_simp [h2mt]
    ring
  have hlog1mx : Real.log (2 * (1 - t) / (2 - t)) = Real.log 2 + Real.log (1 - t) - Real.log (2 - t) := by
    have h21 : 2 * (1 - t) ≠ 0 := by
      exact mul_ne_zero (by norm_num) (ne_of_gt (by linarith : 0 < 1 - t))
    have hdiv := Real.log_div h21 h2mt
    have hmul := Real.log_mul (by norm_num : (2 : ℝ) ≠ 0) h1mt
    rw [hdiv, hmul]
  have hlog1px : Real.log (2 / (2 - t)) = Real.log 2 - Real.log (2 - t) := by
    have hdiv := Real.log_div (by norm_num : (2 : ℝ) ≠ 0) h2mt
    rw [hdiv]
  have hlog2mt : Real.log (2 - t) = Real.log 2 + Real.log (1 - t / 2) := by
    have hdiv' : 2 - t = 2 * (1 - t / 2) := by ring
    rw [hdiv']
    rw [Real.log_mul]
    · norm_num
    · exact ne_of_gt (by linarith : 0 < 1 - t / 2)
  rw [h1mx, h1px, hlog1mx, hlog1px, hlog2mt]
  field_simp [h2mt, h1mt, ne_of_gt ht0, htm1]
  ring

/-- `I10 = ∫₀¹ W0(t)·H1(t)/t dt`. -/
def I10 : ℝ := ∫ t in (0 : ℝ)..1, W0 t * H1 t / t

/-- `I11 = ∫₀¹ W0(t)·H1(t)/(1-t) dt`. -/
def I11 : ℝ := ∫ t in (0 : ℝ)..1, W0 t * H1 t / (1 - t)

/-- `I12 = ∫₀¹ W0(t)·H1(t)/(2-t) dt`. -/
def I12 : ℝ := ∫ t in (0 : ℝ)..1, W0 t * H1 t / (2 - t)

/-- `I20 = ∫₀¹ W0(t)·H2(t)/t dt`. -/
def I20 : ℝ := ∫ t in (0 : ℝ)..1, W0 t * H2 t / t

/-- `I21 = ∫₀¹ W0(t)·H2(t)/(1-t) dt`. -/
def I21 : ℝ := ∫ t in (0 : ℝ)..1, W0 t * H2 t / (1 - t)

/-- `I22 = ∫₀¹ W0(t)·H2(t)/(2-t) dt`. -/
def I22 : ℝ := ∫ t in (0 : ℝ)..1, W0 t * H2 t / (2 - t)

/-! ## The bridge combination (Q6047 (5.6)) -/

/-- The six-integral combination that equals Aquad (Q6047 (5.6)). -/
def sixIntegralCombination : ℝ :=
  -2 * I10 - 2 * I11 + 2 * I12 + 4 * I20 + 6 * I21 - 5 * I22

/-! ## The endpoint evaluation bridge (Q6047 §7, (7.1)-(7.6)) -/

/-- The bridge value `11*Tminus + 2*Tplus - 12 L² Z2 - 14 L Z3 + (41/10) Z2²`
(Q6047 (V)).  `Tplus = cubicLinearEulerValue24`,
`Tminus = alternatingCubicLinearEulerValue24`. -/
def bridgeValue : ℝ :=
  11 * alternatingCubicLinearEulerValue24 +
    2 * cubicLinearEulerValue24 -
    12 * Real.log 2 ^ 2 * (Real.pi ^ 2 / 6) -
    14 * Real.log 2 * zeta3_24 +
    (41 / 10) * (Real.pi ^ 2 / 6) ^ 2

/-- Sanity: the bridge value is exactly the declared target value. -/
theorem bridgeValue_eq_alternatingQuadraticEulerValue24 :
    bridgeValue = alternatingQuadraticEulerValue24 := by
  unfold bridgeValue alternatingQuadraticEulerValue24
    cubicLinearEulerValue24 alternatingCubicLinearEulerValue24
  ring

/-! ## Layer E: the six endpoint values collapse to ONE weight-four integral -/

/-- `K = ∫₀¹ log²x · log(1+x)/(1+x) dx`.

Modulo the lower-product space `span_ℚ{L²ζ₂, Lζ₃, ζ₂²}` the six endpoint values
`I10 … I22` span a ONE-dimensional space, and this is its generator.  `L⁴`
cancels out of every row once they are written through `K`, which is the
structural reason `I10, I20` share the leading pair `(-2, -1/12)` and
`I12, I21, I22` share `(-6, -1/4)` in the `Li₄(1/2)` basis. -/
noncomputable def quadAltK : ℝ :=
  ∫ x in (0 : ℝ)..1, Real.log x ^ 2 * Real.log (1 + x) / (1 + x)

/-- Layer E, algebraic core: the six endpoint evaluations in `K` normal form,
plus the reduction of `K` to the repo's already-proved cubic-linear constants,
give the bridge value.

This step contains no analysis at all — it is rational arithmetic over the
atoms `polylog4 (1/2)`, `Real.log 2`, `zeta3_24` and `π`.  Its purpose is to
pin down exactly what Layer E still owes: seven integral evaluations and
nothing else.  Every hypothesis below is a true statement, each confirmed
numerically at 40 digits against its integral, and `I10, I11, I20, I22`
additionally have independent analytic derivations that use neither quadrature
nor integer-relation detection. -/
theorem quadAltSixIntegral_eq_bridgeValue
    (hK : quadAltK = (1 / 5) * (Real.pi ^ 2 / 6) ^ 2
        - 2 * alternatingCubicLinearEulerValue24 - 2 * cubicLinearEulerValue24)
    (h10 : I10 = -(1 / 2) * quadAltK
        - (3 / 2) * Real.log 2 ^ 2 * (Real.pi ^ 2 / 6)
        - (13 / 20) * (Real.pi ^ 2 / 6) ^ 2)
    (h11 : I11 = -(7 / 2) * Real.log 2 * zeta3_24
        + (3 / 4) * (Real.pi ^ 2 / 6) ^ 2)
    (h12 : I12 = -(3 / 2) * quadAltK
        + (3 / 2) * Real.log 2 ^ 2 * (Real.pi ^ 2 / 6)
        - (9 / 20) * (Real.pi ^ 2 / 6) ^ 2)
    (h20 : I20 = -(1 / 2) * quadAltK - (1 / 2) * (Real.pi ^ 2 / 6) ^ 2)
    (h21 : I21 = -(3 / 2) * quadAltK
        - 3 * Real.log 2 ^ 2 * (Real.pi ^ 2 / 6)
        + (7 / 4) * Real.log 2 * zeta3_24
        + (3 / 10) * (Real.pi ^ 2 / 6) ^ 2)
    (h22 : I22 = -(3 / 2) * quadAltK + (1 / 20) * (Real.pi ^ 2 / 6) ^ 2) :
    sixIntegralCombination = bridgeValue := by
  unfold sixIntegralCombination bridgeValue
  rw [h10, h11, h12, h20, h21, h22, hK]
  unfold cubicLinearEulerValue24 alternatingCubicLinearEulerValue24
  ring

/-! ## The target theorem -/

/-- Summability of the outer-alternating quadratic terms (already proved in
`Problem24`). -/
theorem summable_quadAlt : Summable alternatingQuadraticEulerTerm24 :=
  summable_alternatingQuadraticEulerTerm24

/-! ## Layer B2: normalization `M(0) = 0` (Q6047 (2.9)) -/

/-- `Li2(1/2) = π²/12 − L²/2` (from the dilog reflection identity at z=1/2). -/
theorem quadAlt_dilog_half :
    dilog (1 / 2) = Real.pi ^ 2 / 12 - Real.log 2 ^ 2 / 2 := by
  have hrefl := dilog_reflection
    (by norm_num : 0 < (1 : ℝ) / 2) (by norm_num : (1 : ℝ) / 2 < 1)
  rw [show (1 : ℝ) - (1 : ℝ) / 2 = (1 : ℝ) / 2 by norm_num] at hrefl
  have hlog : Real.log (1 / 2) = -Real.log 2 := by
    rw [Real.log_div]
    · simp
    · norm_num
    · norm_num
  rw [hlog] at hrefl
  simp at hrefl
  nlinarith

/-- `M(0) = 0`: the normalization constant of the primitive (Q6047 (2.9)). -/
theorem quadAltMclosed_zero : quadAltMclosed 0 = 0 := by
  unfold quadAltMclosed
  norm_num
  rw [quadAlt_dilog_half]
  ring

/-! ## Layer B3: the generating function of `S` (Q6047 §2.3) -/

/-- `S_0 = P(0)² − H(0)² = 0`, so the telescoping prefix has zero start. -/
theorem quadAltS_zero : quadAltS 0 = 0 := by
  unfold quadAltS
  norm_num [parityRemainder24, harmonicSquare24, signedHarmonic24,
    harmonicNumber_one]

/-- The increment `S_{n+1} − S_n` at index `n` (Q6047 (2.2) in `quadAltS` indexing). -/
def quadAltIncrement (n : ℕ) : ℝ :=
  quadAltS n - if n = 0 then 0 else quadAltS (n - 1)

/-- Telescoping: `S_{n+1} = Σ_{j≤n} increment j` (Q6047 §2.3 prefix). -/
theorem quadAltS_eq_sum_increment (n : ℕ) :
    quadAltS n = ∑ j ∈ Finset.range (n + 1), quadAltIncrement j := by
  induction n with
  | zero =>
      simp [quadAltIncrement, quadAltS_zero]
  | succ n ih =>
      rw [Finset.sum_range_succ]
      rw [← ih]
      unfold quadAltIncrement
      have hn : n + 1 ≠ 0 := by omega
      simp [hn]


/-! ## Layer B3b: the level-two `P` series and `M` series (Q6047 §2) -/

/-- `P(0) = 0`. -/
theorem quadAltP_zero : parityRemainder24 0 = 0 := by
  unfold parityRemainder24
  simp [harmonicNumber_zero, signedHarmonic24]

/-- `Σ_{k≥0} P(k) t^k = F(t)` for `|t| < 1`, `t ≠ 0` (Q6047 (2.1) shifted). -/
theorem quadAltP_hasSum {t : ℝ} (ht : |t| < 1) (htne : t ≠ 0) :
    HasSum (fun k : ℕ => parityRemainder24 k * t ^ k) (quadAltFclosed t) := by
  have hgen := parityRemainder24_generating_hasSum ht htne
  -- hgen : HasSum (fun n => P(n+1) t^n) (F(t)/t)
  have hmain : HasSum (fun n : ℕ => parityRemainder24 (n + 1) * t ^ (n + 1))
      (quadAltFclosed t) := by
    convert hgen.mul_left t using 1
    · funext n
      rw [pow_succ]
      ring
    · have h1mt : 1 - t ≠ 0 := by
        nlinarith [abs_lt.mp ht]
      unfold quadAltFclosed
      field_simp [htne, h1mt]
      ring
  -- reindex: Σ_{k≥0} P(k) t^k = Σ_{n≥0} P(n+1) t^(n+1) + P(0) t^0, P(0) = 0
  have htail : HasSum (fun n : ℕ => parityRemainder24 (n + 1) * t ^ (n + 1))
      (quadAltFclosed t - ∑ i ∈ Finset.range 1, parityRemainder24 i * t ^ i) := by
    simpa [quadAltP_zero] using hmain
  exact (hasSum_nat_add_iff' (f := fun k : ℕ => parityRemainder24 k * t ^ k)
    (g := quadAltFclosed t) (1 : ℕ)).mp htail


/-- `Σ_{k≥0} c_{k+1} P(k) t^k = F(t) − 2F(−t)` with `c_{k+1} = 1 + 2(−1)^{k+1}`
(Q6047 (2.6)). -/
theorem quadAltC_hasSum {t : ℝ} (ht : |t| < 1) (htne : t ≠ 0) :
    HasSum (fun k : ℕ => (1 + 2 * (-1 : ℝ) ^ (k + 1)) * parityRemainder24 k * t ^ k)
      (quadAltFclosed t - 2 * quadAltFclosed (-t)) := by
  have hP := quadAltP_hasSum ht htne
  have hPn : HasSum (fun k : ℕ => parityRemainder24 k * (-t) ^ k)
      (quadAltFclosed (-t)) :=
    quadAltP_hasSum (by simpa [abs_neg] using ht) (by simpa [neg_ne_zero] using htne)
  -- Σ (−1)^{k+1} P(k) t^k = −Σ P(k) (−t)^k = −F(−t)
  have hneg : HasSum (fun k : ℕ => (-1 : ℝ) ^ (k + 1) * parityRemainder24 k * t ^ k)
      (-quadAltFclosed (-t)) := by
    convert hPn.neg using 1
    · funext k
      rw [pow_succ, neg_pow]
      ring
  -- 组合
  convert (hP.add (hneg.mul_left 2)) using 1
  · funext k
    ring
  · ring


/-- `Li2(x) + Li2(−x) = Li2(x²)/2` for `0 < x < 1` (even-term identity). -/
theorem quadAlt_dilog_add_neg {x : ℝ} (hx0 : 0 < x) (hx1 : x < 1) :
    dilog x + dilog (-x) = dilog (x ^ 2) / 2 := by
  let f : ℝ → ℝ := fun t => dilog t + dilog (-t) - dilog (t ^ 2) / 2
  have hf0 : f 0 = 0 := by
    simp [f, dilog_zero]
  have hderiv : ∀ t ∈ Set.Ioo (0 : ℝ) x, HasDerivAt f 0 t := by
    intro t ht
    have ht0 : 0 < t := ht.1
    have ht1 : t < 1 := lt_trans ht.2 hx1
    have habs : |t| < 1 := by
      rw [abs_of_pos ht0]
      exact ht1
    have hd1 : HasDerivAt dilog (-Real.log (1 - t) / t) t :=
      dilog_hasDerivAt ht0 ht1
    have hd2 : HasDerivAt (fun s : ℝ => dilog (-s)) (-Real.log (1 + t) / t) t := by
      have hd2a : HasDerivAt dilog (-Real.log (1 + t) / (-t)) (-t) := by
        convert dilog_hasDerivAt_of_abs_lt_one (by simpa using habs)
          (neg_ne_zero.mpr (ne_of_gt ht0)) using 1
        congr 1
        ring
      have hneg : HasDerivAt (fun s : ℝ => -s) (-1) t := by
        simpa using (hasDerivAt_id t).neg
      have hc := HasDerivAt.comp (h := fun s : ℝ => -s) t hd2a hneg
      convert hc using 1
      field_simp
    have hd3 : HasDerivAt (fun s : ℝ => dilog (s ^ 2))
        (-2 * Real.log (1 - t ^ 2) / t) t := by
      have ht2 : t ^ 2 < 1 := (sq_lt_one_iff_abs_lt_one t).mpr habs
      have hd3a : HasDerivAt dilog (-Real.log (1 - t ^ 2) / t ^ 2) (t ^ 2) :=
        dilog_hasDerivAt (by positivity) ht2
      have hpow2 : HasDerivAt (fun s : ℝ => s ^ 2) (2 * t) t := by
        simpa using (hasDerivAt_id t).pow 2
      have hc := HasDerivAt.comp (h := fun s : ℝ => s ^ 2) t hd3a hpow2
      convert hc using 1
      field_simp
    -- f' = hd1 + hd2 − hd3/2 = 0
    have hcomb : HasDerivAt f
        (-Real.log (1 - t) / t + (-Real.log (1 + t) / t) -
          (-2 * Real.log (1 - t ^ 2) / t) / 2) t := by
      unfold f
      convert (((hd1.add hd2).sub (hd3.div_const 2))) using 1
    have hlog : Real.log (1 - t ^ 2) = Real.log (1 - t) + Real.log (1 + t) := by
      rw [show 1 - t ^ 2 = (1 - t) * (1 + t) by ring]
      rw [Real.log_mul]
      · exact ne_of_gt (by nlinarith : 0 < 1 - t)
      · exact ne_of_gt (by nlinarith : 0 < 1 + t)
    have hD0 : -Real.log (1 - t) / t + (-Real.log (1 + t) / t) -
          (-2 * Real.log (1 - t ^ 2) / t) / 2 = 0 := by
      rw [hlog]
      calc
        -Real.log (1 - t) / t + (-Real.log (1 + t) / t) -
              (-2 * (Real.log (1 - t) + Real.log (1 + t)) / t) / 2
            = (-Real.log (1 - t) - Real.log (1 + t) +
                Real.log (1 - t) + Real.log (1 + t)) / t := by
              ring_nf
        _ = 0 := by
          ring
    convert hcomb using 1
    exact hD0.symm
  have hcont : ContinuousOn f (Set.Icc (0 : ℝ) x) := by
    unfold f
    have hc1 : ContinuousOn (fun t : ℝ => dilog t) (Set.Icc (0 : ℝ) x) :=
      dilog_continuousOn_unit.mono (by
        rintro y ⟨hy0, hyx⟩
        constructor
        · linarith
        · exact le_trans hyx (le_of_lt hx1))
    have hc2 : ContinuousOn (fun t : ℝ => dilog (-t)) (Set.Icc (0 : ℝ) x) :=
      dilog_continuousOn_unit.comp (by fun_prop) (by
        rintro y ⟨hy0, hyx⟩
        constructor
        · linarith [hyx, le_of_lt hx1]
        · linarith)
    have hc3 : ContinuousOn (fun t : ℝ => dilog (t ^ 2) / 2) (Set.Icc (0 : ℝ) x) :=
      (dilog_continuousOn_unit.comp (by fun_prop) (by
        intro y hy
        rcases (Set.mem_Icc.mp hy) with ⟨hy0, hyx⟩
        have hy_abs : |y| < 1 := by
          rw [abs_of_nonneg hy0]
          exact lt_of_le_of_lt hyx hx1
        have hsq : y ^ 2 < 1 := (sq_lt_one_iff_abs_lt_one y).mpr hy_abs
        constructor
        · nlinarith [sq_nonneg y]
        · exact le_of_lt hsq)).div_const 2
    refine ((hc1.add (hc2.sub hc3)).congr ?_)
    intro t ht
    simp
    ring
  have hint : IntervalIntegrable (fun _ : ℝ => (0 : ℝ)) MeasureTheory.volume 0 x := by
    exact intervalIntegral.intervalIntegrable_const
  have hfund := intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le
    (a := 0) (b := x) (f := f) (f' := fun _ : ℝ => (0 : ℝ))
    (le_of_lt hx0) hcont hderiv hint
  have hzero : (∫ _t : ℝ in (0 : ℝ)..x, (0 : ℝ)) = 0 := by
    simp
  rw [hzero] at hfund
  have hfx : f x = 0 := by
    linarith
  unfold f at hfx
  linarith


/-- Moment term for the `M` series: `c_{k+1} P(k) t^k` on `[0,x]`. -/
def quadAltMMoment (k : ℕ) (t : ℝ) : ℝ :=
  (1 + 2 * (-1 : ℝ) ^ (k + 1)) * parityRemainder24 k * t ^ k

/-- Pointwise HasSum for the `M` moment: `Σ c_{k+1}P(k) t^k = F(t) − 2F(−t)`
for `0 < t < 1` (Q6047 (2.6)). -/
theorem quadAltMMoment_hasSum_pointwise {x t : ℝ} (ht0 : 0 < t) (htx : t ≤ x) (hx1 : x < 1) :
    HasSum (fun k : ℕ => quadAltMMoment k t) (quadAltFclosed t - 2 * quadAltFclosed (-t)) := by
  have ht1 : t < 1 := lt_of_le_of_lt htx hx1
  exact quadAltC_hasSum
    (by rw [abs_of_pos ht0]; exact ht1) (ne_of_gt ht0)

/-- `Mclosed` is the FTC integral of `F − 2F(−t)`: `∫₀ˣ (F − 2F(−t)) dt = Mclosed x`. -/
theorem quadAltMclosed_eq_integral {x : ℝ} (hx0 : 0 < x) (hx1 : x < 1) :
    (∫ t : ℝ in (0 : ℝ)..x, (quadAltFclosed t - 2 * quadAltFclosed (-t))) =
      quadAltMclosed x := by
  have hd : ContinuousOn (fun t : ℝ => dilog ((1 + t) / 2)) (Set.Icc (0 : ℝ) x) := by
    refine (dilog_continuousOn_unit.comp (by fun_prop) ?_)
    intro t ht
    rcases (Set.mem_Icc.mp ht) with ⟨ht0', htx⟩
    constructor
    · nlinarith
    · nlinarith [htx, hx1]
  have hlog1 : ContinuousOn (fun t : ℝ => Real.log (1 - t)) (Set.Icc (0 : ℝ) x) := by
    refine ContinuousOn.log (by fun_prop) ?_
    intro t ht
    rcases (Set.mem_Icc.mp ht) with ⟨_, htx⟩
    exact ne_of_gt (by nlinarith [htx, hx1])
  have hlog2 : ContinuousOn (fun t : ℝ => Real.log (1 + t)) (Set.Icc (0 : ℝ) x) := by
    refine ContinuousOn.log (by fun_prop) ?_
    intro t ht
    rcases (Set.mem_Icc.mp ht) with ⟨ht0', _⟩
    exact ne_of_gt (by nlinarith [ht0'])
  have hconst : ContinuousOn (fun _ : ℝ => Real.pi ^ 2 / 6 - Real.log 2 ^ 2)
      (Set.Icc (0 : ℝ) x) := continuousOn_const
  -- a²/2 + b² + 2ab + 2Lb + c − 2·dilog((1+t)/2)
  have hcomp := ((((hlog1.pow 2).div_const 2).add (hlog2.pow 2)).add
    ((hlog1.mul hlog2).const_mul 2)).add ((hlog2.const_mul (2 * Real.log 2)))
  have hcomp2 := (hcomp.add hconst).sub (hd.const_mul 2)
  have hcontM : ContinuousOn quadAltMclosed (Set.Icc (0 : ℝ) x) := by
    unfold quadAltMclosed
    refine (hcomp2.congr ?_)
    intro t ht
    simp
    ring
  have hderivM : ∀ t ∈ Set.Ioo (0 : ℝ) x,
      HasDerivAt quadAltMclosed (quadAltFclosed t - 2 * quadAltFclosed (-t)) t := by
    intro t ht
    exact quadAltMclosed_hasDerivAt ht.1 (lt_trans ht.2 hx1)
  have hcontF : ContinuousOn
      (fun t : ℝ => quadAltFclosed t - 2 * quadAltFclosed (-t)) (Set.Icc (0 : ℝ) x) := by
    have hden1 : ContinuousOn (fun t : ℝ => 1 - t) (Set.Icc (0 : ℝ) x) := by fun_prop
    have hF1 : ContinuousOn (fun t : ℝ => quadAltFclosed t) (Set.Icc (0 : ℝ) x) := by
      unfold quadAltFclosed
      have hnum : ContinuousOn (fun t : ℝ => -(Real.log (1 - t) + 2 * Real.log (1 + t)))
          (Set.Icc (0 : ℝ) x) := (hlog1.add (hlog2.const_mul 2)).neg
      exact hnum.div hden1 (by
        intro t ht
        rcases (Set.mem_Icc.mp ht) with ⟨_, htx⟩
        exact ne_of_gt (by nlinarith [htx, hx1]))
    have hden2 : ContinuousOn (fun t : ℝ => 1 + t) (Set.Icc (0 : ℝ) x) := by fun_prop
    have hF2 : ContinuousOn (fun t : ℝ => quadAltFclosed (-t)) (Set.Icc (0 : ℝ) x) := by
      unfold quadAltFclosed
      have hnum : ContinuousOn (fun t : ℝ => -(Real.log (1 + t) + 2 * Real.log (1 - t)))
          (Set.Icc (0 : ℝ) x) := (hlog2.add (hlog1.const_mul 2)).neg
      have hdiv := hnum.div hden2 (by
        intro t ht
        rcases (Set.mem_Icc.mp ht) with ⟨ht0', _⟩
        exact ne_of_gt (by nlinarith [ht0']))
      simpa [sub_neg_eq_add] using hdiv
    refine (hF1.sub (hF2.const_mul 2)).congr ?_
    intro t ht
    rfl
  have hcontF' : ContinuousOn
      (fun t : ℝ => quadAltFclosed t - 2 * quadAltFclosed (-t))
      (Set.uIcc (0 : ℝ) x) := by
    simpa [uIcc_of_le (le_of_lt hx0)] using hcontF
  have hintM : IntervalIntegrable
      (fun t : ℝ => quadAltFclosed t - 2 * quadAltFclosed (-t)) MeasureTheory.volume 0 x :=
    hcontF'.intervalIntegrable
  have hfund := intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le
    (a := 0) (b := x) (f := quadAltMclosed)
    (f' := fun t : ℝ => quadAltFclosed t - 2 * quadAltFclosed (-t))
    (le_of_lt hx0) hcontM hderivM hintM
  rw [quadAltMclosed_zero] at hfund
  linarith


/-- Integral of the `M` moment over `[0,x]` = the `M` series term. -/
theorem quadAltMMoment_integral {x : ℝ} (hx0 : 0 < x) (k : ℕ) :
    (∫ t : ℝ in (0 : ℝ)..x, quadAltMMoment k t) =
      (1 + 2 * (-1 : ℝ) ^ (k + 1)) * parityRemainder24 k *
        x ^ (k + 1) / (k + 1 : ℝ) := by
  unfold quadAltMMoment
  rw [show (fun t : ℝ => (1 + 2 * (-1 : ℝ) ^ (k + 1)) * parityRemainder24 k * t ^ k) =
      fun t => ((1 + 2 * (-1 : ℝ) ^ (k + 1)) * parityRemainder24 k) * t ^ k by
        funext t
        ring,
    intervalIntegral.integral_const_mul, integral_pow]
  ring

/-- Bound `|c_{k+1}| ≤ 3`. -/
theorem quadAltC_abs_le (k : ℕ) : |1 + 2 * (-1 : ℝ) ^ (k + 1)| ≤ 3 := by
  have h2 : |2 * (-1 : ℝ) ^ (k + 1)| = 2 := by
    rw [abs_mul, abs_pow]
    norm_num
  calc
    |1 + 2 * (-1 : ℝ) ^ (k + 1)| ≤ |(1 : ℝ)| + |2 * (-1 : ℝ) ^ (k + 1)| := abs_add_le _ _
    _ = 1 + 2 := by
      rw [h2]
      norm_num
    _ = 3 := by norm_num


/-- Bound `|moment k t| ≤ 9·H(k+1)·x^k` on `0 ≤ t ≤ x`. -/
theorem quadAltMMoment_norm_le {x t : ℝ} (hx0 : 0 ≤ x) (ht0 : 0 ≤ t) (htx : t ≤ x) (k : ℕ) :
    ‖quadAltMMoment k t‖ ≤ 9 * (harmonicNumber (k + 1) : ℝ) * x ^ k := by
  unfold quadAltMMoment
  rw [Real.norm_eq_abs, abs_mul, abs_mul]
  calc
    |1 + 2 * (-1 : ℝ) ^ (k + 1)| * |parityRemainder24 k| * |t ^ k|
        ≤ 3 * (3 * harmonicNumber k) * |t ^ k| := by
          gcongr
          · exact quadAltC_abs_le k
          · exact abs_parityRemainder24_le k
    _ = 9 * (harmonicNumber k : ℝ) * |t ^ k| := by ring
    _ ≤ 9 * (harmonicNumber (k + 1) : ℝ) * |t ^ k| := by
      have h1 : 9 * harmonicNumber k ≤ 9 * harmonicNumber (k + 1) := by
        exact mul_le_mul_of_nonneg_left (harmonicNumber_mono (Nat.le_succ k)) (by positivity)
      exact mul_le_mul_of_nonneg_right h1 (by positivity)
    _ ≤ 9 * (harmonicNumber (k + 1) : ℝ) * x ^ k := by
      have h2 : |t ^ k| ≤ x ^ k := by
        rw [abs_of_nonneg (pow_nonneg ht0 k)]
        exact pow_le_pow_left₀ ht0 htx k
      have hc : 0 ≤ 9 * harmonicNumber (k + 1) := by
        exact mul_nonneg (by norm_num) (harmonicNumber_nonneg (k + 1))
      exact mul_le_mul_of_nonneg_left h2 hc

/-- Summability of the `M`-moment norm integrals (dominance by `Σ 9H(k+1)x^k`). -/
theorem quadAltMMoment_integral_norm_summable {x : ℝ} (hx0 : 0 ≤ x) (hx1 : x < 1) :
    Summable (fun k : ℕ => ∫ t : ℝ in Set.Ioc (0 : ℝ) x, ‖quadAltMMoment k t‖) := by
  have hH (k : ℕ) : (harmonicNumber (k + 1) : ℝ) ≤ (k + 1 : ℝ) := by
    have hsum : (∑ j ∈ Finset.range (k + 1), (1 : ℝ) / (j + 1 : ℝ)) ≤ (k + 1 : ℝ) := by
      calc
        (∑ j ∈ Finset.range (k + 1), (1 : ℝ) / (j + 1 : ℝ))
            ≤ ∑ _j ∈ Finset.range (k + 1), (1 : ℝ) := by
              exact Finset.sum_le_sum (fun j hj => by
                rw [one_div]
                have hj1 : (1 : ℝ) ≤ (j + 1 : ℝ) := by
                  exact_mod_cast (by omega : 1 ≤ j + 1)
                have hjpos : (0 : ℝ) < (j + 1 : ℝ) := by
                  exact_mod_cast (by omega : 0 < j + 1)
                exact (inv_le_one₀ hjpos).mpr hj1)
        _ = (k + 1 : ℝ) := by simp
    simpa [harmonicNumber] using hsum
  have hgeom : Summable (fun k : ℕ => 9 * (harmonicNumber (k + 1) : ℝ) * x ^ k) := by
    have hxnorm : ‖x‖ < 1 := by
      rw [Real.norm_eq_abs, abs_of_nonneg hx0]
      exact hx1
    have hk : Summable (fun k : ℕ => (k : ℝ) * x ^ k) := by
      simpa using (summable_pow_mul_geometric_of_norm_lt_one 1 (r := x) hxnorm)
    have h1 : Summable (fun k : ℕ => x ^ k) :=
      summable_geometric_of_lt_one hx0 hx1
    have hsum : Summable (fun k : ℕ => (k + 1 : ℝ) * x ^ k) := by
      -- (k+1)x^k = k·x^k + x^k
      simpa [add_mul, one_mul, Nat.cast_add, Nat.cast_one] using
        (hk.add h1).congr (fun k => by ring)
    refine (hsum.mul_left 9).of_nonneg_of_le ?_ ?_
    · intro k
      exact mul_nonneg (mul_nonneg (by norm_num) (harmonicNumber_nonneg (k + 1)))
        (pow_nonneg hx0 k)
    · intro k
      have hxnonneg : 0 ≤ x ^ k := pow_nonneg hx0 k
      calc
        9 * harmonicNumber (k + 1) * x ^ k = 9 * (harmonicNumber (k + 1) * x ^ k) := by ring
        _ ≤ 9 * ((k + 1 : ℝ) * x ^ k) := by
          exact mul_le_mul_of_nonneg_left
            (mul_le_mul_of_nonneg_right (hH k) hxnonneg) (by norm_num)
  refine hgeom.of_nonneg_of_le ?_ ?_
  · intro k
    exact MeasureTheory.integral_nonneg (fun t => norm_nonneg _)
  · intro k
    -- 目标是集合积分 ∫ t in Ioc 0 x, ‖moment k t‖ ≤ 常数, 不是两个区间积分比大小,
    -- 所以走 setIntegral_mono_on 压成常数积分, 再用 |Ioc 0 x| = x ≤ 1 收口。
    have hC : (0 : ℝ) ≤ 9 * (harmonicNumber (k + 1) : ℝ) * x ^ k :=
      mul_nonneg (mul_nonneg (by norm_num) (harmonicNumber_nonneg (k + 1))) (pow_nonneg hx0 k)
    have hcont : ContinuousOn (fun t : ℝ => ‖quadAltMMoment k t‖) (Set.uIcc (0 : ℝ) x) := by
      unfold quadAltMMoment
      exact ((continuousOn_pow k).const_mul
        ((1 + 2 * (-1 : ℝ) ^ (k + 1)) * parityRemainder24 k)).norm
    have hIon : MeasureTheory.IntegrableOn (fun t : ℝ => ‖quadAltMMoment k t‖)
        (Set.Ioc (0 : ℝ) x) MeasureTheory.volume :=
      (hcont.integrableOn_compact isCompact_uIcc).mono_set (by
        intro t ht
        rcases (Set.mem_Ioc.mp ht) with ⟨ht0, htx⟩
        simp [Set.uIcc_of_le hx0, ht0.le, htx])
    calc (∫ t : ℝ in Set.Ioc (0 : ℝ) x, ‖quadAltMMoment k t‖)
        ≤ ∫ _t : ℝ in Set.Ioc (0 : ℝ) x, 9 * (harmonicNumber (k + 1) : ℝ) * x ^ k :=
          MeasureTheory.setIntegral_mono_on hIon
            (MeasureTheory.integrableOn_const
              (hs := by simp [Real.volume_Ioc]))
            measurableSet_Ioc
            (fun t ht => quadAltMMoment_norm_le hx0 ht.1.le ht.2 k)
      _ = 9 * (harmonicNumber (k + 1) : ℝ) * x ^ k * x := by
          rw [MeasureTheory.setIntegral_const]
          simp [MeasureTheory.measureReal_def, Real.volume_Ioc, hx0, smul_eq_mul,
            mul_comm]
      _ ≤ 9 * (harmonicNumber (k + 1) : ℝ) * x ^ k := by nlinarith [hC, hx0, hx1.le]

/-- The `M` series: `Σ_{k≥0} c_{k+1}P(k) x^(k+1)/(k+1) = Mclosed x` for `0 < x < 1`
(Q6047 (2.5)–(2.7), via termwise integration). -/
theorem quadAltM_hasSum {x : ℝ} (hx0 : 0 < x) (hx1 : x < 1) :
    HasSum (fun k : ℕ => (1 + 2 * (-1 : ℝ) ^ (k + 1)) * parityRemainder24 k *
        x ^ (k + 1) / (k + 1 : ℝ))
      (quadAltMclosed x) := by
  have hint : ∀ k : ℕ, MeasureTheory.Integrable (quadAltMMoment k)
      (MeasureTheory.volume.restrict (Set.Ioc (0 : ℝ) x)) := by
    intro k
    have hcont : ContinuousOn (quadAltMMoment k) (Set.uIcc (0 : ℝ) x) := by
      unfold quadAltMMoment
      exact (continuousOn_pow k).const_mul
        ((1 + 2 * (-1 : ℝ) ^ (k + 1)) * parityRemainder24 k)
    have hI : MeasureTheory.IntegrableOn (quadAltMMoment k)
        (Set.uIcc (0 : ℝ) x) MeasureTheory.volume :=
      hcont.integrableOn_compact isCompact_uIcc
    exact hI.mono_set (by
      intro t ht
      rcases (Set.mem_Ioc.mp ht) with ⟨ht0, htx⟩
      simp [Set.uIcc_of_le (le_of_lt hx0), ht0.le, htx])
  have hnorm_sum : Summable (fun k : ℕ =>
      ∫ t : ℝ in Set.Ioc (0 : ℝ) x, ‖quadAltMMoment k t‖) :=
    quadAltMMoment_integral_norm_summable hx0.le hx1
  have hsum := MeasureTheory.hasSum_integral_of_summable_integral_norm
    (μ := MeasureTheory.volume.restrict (Set.Ioc (0 : ℝ) x)) hint hnorm_sum
  -- hsum : HasSum (fun k => ∫ t in Ioc 0 x, moment k t) (∫ t in Ioc 0 x, ∑' k, moment k t)
  -- 全程留在集合积分上, 只在需要闭式时才换成区间积分 —— 混用两种记号正是之前卡住的原因。
  have hterm : ∀ k : ℕ, (∫ t : ℝ in Set.Ioc (0 : ℝ) x, quadAltMMoment k t)
      = (1 + 2 * (-1 : ℝ) ^ (k + 1)) * parityRemainder24 k * x ^ (k + 1) / (k + 1 : ℝ) := by
    intro k
    rw [← intervalIntegral.integral_of_le (le_of_lt hx0), quadAltMMoment_integral hx0]
  have hlim : (∫ t : ℝ in Set.Ioc (0 : ℝ) x, ∑' k : ℕ, quadAltMMoment k t)
      = quadAltMclosed x := by
    rw [← quadAltMclosed_eq_integral hx0 hx1,
      intervalIntegral.integral_of_le (le_of_lt hx0)]
    refine MeasureTheory.setIntegral_congr_ae measurableSet_Ioc ?_
    filter_upwards with t ht
    exact (quadAltMMoment_hasSum_pointwise ht.1 ht.2 hx1).tsum_eq
  simpa only [hterm, hlim] using hsum


/-! ## Layer C: coefficient integration (Q6047 §3) -/

/-- `Σ_{k≥1} c_{k+1}P(k)y^k/(k+1) = Mclosed(y)/y`: the `P`-part of the
increment series, reindexed from `quadAltM_hasSum` (Q6047 (2.5)). -/
theorem quadAltMPart_hasSum {y : ℝ} (hy0 : 0 < y) (hy1 : y < 1) :
    HasSum (fun m : ℕ => (1 + 2 * (-1 : ℝ) ^ (m + 2)) * parityRemainder24 (m + 1) *
        y ^ (m + 1) / (m + 2 : ℝ))
      (quadAltMclosed y / y) := by
  have H1 := quadAltM_hasSum hy0 hy1
  -- H1 : HasSum (fun k => c_{k+1} P(k) y^(k+1)/(k+1)) (Mclosed y)
  have hf0 : (1 + 2 * (-1 : ℝ) ^ (0 + 1)) * parityRemainder24 0 * y ^ (0 + 1) / (0 + 1 : ℝ) = 0 := by
    simp [quadAltP_zero]
  have H1a : HasSum (fun n : ℕ => (1 + 2 * (-1 : ℝ) ^ (n + 1 + 1)) * parityRemainder24 (n + 1) *
        y ^ (n + 1 + 1) / ((n + 1 : ℝ) + 1))
      (quadAltMclosed y - (1 + 2 * (-1 : ℝ) ^ (0 + 1)) * parityRemainder24 0 *
        y ^ (0 + 1) / (0 + 1 : ℝ)) := by
    simpa [hf0] using
      (hasSum_nat_add_iff' (f := fun k : ℕ =>
        (1 + 2 * (-1 : ℝ) ^ (k + 1)) * parityRemainder24 k * y ^ (k + 1) / (k + 1 : ℝ))
        (1 : ℕ)).mpr H1
  -- multiply by 1/y: y^(n+2)/((n+2)·y) = y^(n+1)/(n+2)
  have H1b := H1a.mul_left (1 / y)
  convert H1b using 1
  · funext n
    field_simp [ne_of_gt hy0]
    ring
  · field_simp [ne_of_gt hy0]
    simp [quadAltP_zero]


/-- `Σ_{j≥0} y^(j+1)/(j+2)² = (dilog y − y)/y` for `|y| < 1`, `y ≠ 0`.
(NOTE: `(j + 2 : ℝ)` must be written as `↑(j+2)` — the single-cast form —
because `(j + 2 : ℝ)` elaborates to `↑j + 2` and forces a whnf blowup when
unified with the `f (n+1)` terms from `hasSum_nat_add_iff'`.) -/
theorem quadAlt_dilog_tail_hasSum {y : ℝ} (hy : |y| < 1) (hyne : y ≠ 0) :
    HasSum (fun j : ℕ => y ^ (j + 1) / (↑(j + 2) : ℝ) ^ 2) ((dilog y - y) / y) := by
  have H : HasSum (fun m : ℕ => y ^ (m + 1) / (↑(m + 1) : ℝ) ^ 2) (dilog y) := by
    unfold dilog
    exact (dilog_summable hy.le).hasSum
  have H1 : HasSum (fun n : ℕ => y ^ (n + 1 + 1) / (↑n + 1 + 1) ^ 2)
      (dilog y - y) := by
    simpa [show y ^ (0 + 1) / (↑(0 + 1) : ℝ) ^ 2 = y by norm_num] using
      (hasSum_nat_add_iff' (f := fun m : ℕ => y ^ (m + 1) / (↑(m + 1) : ℝ) ^ 2)
        (1 : ℕ)).mpr H
  have H2 := H1.mul_left (1 / y)
  have hterm : ∀ n : ℕ, y ^ (n + 1) / (↑(n + 2) : ℝ) ^ 2 =
      1 / y * (y ^ (n + 1 + 1) / (↑n + 1 + 1) ^ 2) := by
    intro n
    field_simp [hyne]
    push_cast
    ring
  have H3 : HasSum (fun n : ℕ => y ^ (n + 1) / (↑(n + 2) : ℝ) ^ 2)
      (1 / y * (dilog y - y)) :=
    H2.congr_fun (fun n => hterm n)
  convert H3 using 1
  ring


/-- `Σ_{j≥0} (−1)^j y^(j+1)/(j+2)² = (dilog(−y) + y)/y` for `|y| < 1`, `y ≠ 0`. -/
theorem quadAlt_dilog_neg_tail_hasSum {y : ℝ} (hy : |y| < 1) (hyne : y ≠ 0) :
    HasSum (fun j : ℕ => (-1 : ℝ) ^ j * y ^ (j + 1) / (↑(j + 2) : ℝ) ^ 2)
      ((dilog (-y) + y) / y) := by
  have H : HasSum (fun m : ℕ => (-y) ^ (m + 1) / (↑(m + 1) : ℝ) ^ 2) (dilog (-y)) := by
    unfold dilog
    exact (dilog_summable (by simpa using hy.le)).hasSum
  have H1 : HasSum (fun n : ℕ => (-y) ^ (n + 1 + 1) / (↑n + 1 + 1) ^ 2)
      (dilog (-y) - (-y)) := by
    simpa [show (-y) ^ (0 + 1) / (↑(0 + 1) : ℝ) ^ 2 = -y by norm_num] using
      (hasSum_nat_add_iff' (f := fun m : ℕ => (-y) ^ (m + 1) / (↑(m + 1) : ℝ) ^ 2)
        (1 : ℕ)).mpr H
  have H2 := H1.mul_left (1 / y)
  have hterm : ∀ n : ℕ, (-1 : ℝ) ^ n * y ^ (n + 1) / (↑(n + 2) : ℝ) ^ 2 =
      1 / y * ((-y) ^ (n + 1 + 1) / (↑n + 1 + 1) ^ 2) := by
    intro n
    field_simp [hyne]
    push_cast
    rw [neg_pow]
    ring
  have H3 : HasSum (fun n : ℕ => (-1 : ℝ) ^ n * y ^ (n + 1) / (↑(n + 2) : ℝ) ^ 2)
      (1 / y * (dilog (-y) - -y)) :=
    H2.congr_fun (fun n => hterm n)
  convert H3 using 1
  ring


/-- Diagonal part of the increment series:
`Σ 4(1+(−1)^{j+2})y^(j+1)/(j+2)² = 4(dilog y + dilog(−y))/y` (Q6047 (2.4)). -/
theorem quadAltDiag_hasSum {y : ℝ} (hy : |y| < 1) (hyne : y ≠ 0) :
    HasSum (fun j : ℕ => 4 * (1 + (-1 : ℝ) ^ (j + 2)) * y ^ (j + 1) / (↑(j + 2) : ℝ) ^ 2)
      (4 * (dilog y + dilog (-y)) / y) := by
  have hp := quadAlt_dilog_tail_hasSum hy hyne
  have hm := quadAlt_dilog_neg_tail_hasSum hy hyne
  have hcomb := hp.add hm
  have hpow (j : ℕ) : (-1 : ℝ) ^ (j + 2) = (-1 : ℝ) ^ j := by
    rw [pow_add]
    norm_num
  convert hcomb.mul_left 4 using 1
  · funext j
    rw [hpow]
    ring
  · ring


/-- `Li2(z) + Li2(−z) = Li2(z²)/2` for `|z| < 1` (general version; the
positive version `quadAlt_dilog_add_neg` covers only `0 < z < 1`). -/
theorem quadAlt_dilog_add_neg' {z : ℝ} (hz : |z| < 1) :
    dilog z + dilog (-z) = dilog (z ^ 2) / 2 := by
  rcases lt_or_ge z 0 with hzneg | hzpos
  · -- z < 0：令 t = −z > 0
    have ht : |(-z)| < 1 := by simpa using hz
    have hzlt1 : z < 1 := lt_of_lt_of_le hzneg (by linarith)
    have htpos : 0 < -z := by linarith
    have h := quadAlt_dilog_add_neg (x := -z) htpos (by
      rw [abs_lt] at hz
      linarith)
    -- h : dilog (−z) + dilog (−(−z)) = dilog ((−z)²)/2
    convert h using 1
    all_goals ring
  · -- z ≥ 0
    rcases lt_or_eq_of_le hzpos with hzpos' | rfl
    · have hz1 : z < 1 := (abs_lt.mp hz).2
      exact quadAlt_dilog_add_neg hzpos' hz1
    · simp [dilog_zero]


/-- `M' = F − 2F(−x)` for all `|x| < 1` (general version; the positive-only
`quadAltMclosed_hasDerivAt` is kept for the `[0,x]` FTC in `quadAltM_hasSum`). -/
theorem quadAltMclosed_hasDerivAt' {x : ℝ} (hx : |x| < 1) :
    HasDerivAt quadAltMclosed (quadAltFclosed x - 2 * quadAltFclosed (-x)) x := by
  have h1mx : 0 < 1 - x := by
    rw [abs_lt] at hx
    nlinarith
  have h1px : 0 < 1 + x := by
    rw [abs_lt] at hx
    nlinarith
  have hmid0 : 0 < (1 + x) / 2 := by
    rw [abs_lt] at hx
    nlinarith
  have hmid1 : (1 + x) / 2 < 1 := by
    rw [abs_lt] at hx
    nlinarith
  have ha : HasDerivAt (fun y : ℝ => Real.log (1 - y)) (-1 / (1 - x)) x := by
    have hc : HasDerivAt (fun y : ℝ => 1 - y) (-1) x := by
      convert (hasDerivAt_const x 1).sub (hasDerivAt_id x) using 1
      simp
    have hlog := Real.hasDerivAt_log (ne_of_gt h1mx)
    have hD : HasDerivAt (fun y : ℝ => Real.log (1 - y)) ((1 - x)⁻¹ * (-1)) x :=
      HasDerivAt.comp (h := fun y : ℝ => 1 - y) x hlog hc
    convert hD using 1
    field_simp
  have hb : HasDerivAt (fun y : ℝ => Real.log (1 + y)) (1 / (1 + x)) x := by
    have hc : HasDerivAt (fun y : ℝ => 1 + y) 1 x := by
      convert (hasDerivAt_const x 1).add (hasDerivAt_id x) using 1
      simp
    have hlog := Real.hasDerivAt_log (ne_of_gt h1px)
    have hD : HasDerivAt (fun y : ℝ => Real.log (1 + y)) ((1 + x)⁻¹ * 1) x :=
      HasDerivAt.comp (h := fun y : ℝ => 1 + y) x hlog hc
    convert hD using 1
    field_simp
  have hd : HasDerivAt (fun y : ℝ => (1 + y) / 2) (1 / 2) x := by
    have hc : HasDerivAt (fun y : ℝ => 1 + y) 1 x := by
      convert (hasDerivAt_const x 1).add (hasDerivAt_id x) using 1
      simp
    convert hc.div_const 2 using 1
  have hD : HasDerivAt (fun y : ℝ => dilog ((1 + y) / 2))
      (-Real.log ((1 - x) / 2) / (1 + x)) x := by
    have hmid := dilog_hasDerivAt hmid0 hmid1
    have hD2 : HasDerivAt (fun y : ℝ => dilog ((1 + y) / 2))
        (-(Real.log (1 - (1 + x) / 2)) / ((1 + x) / 2) * (1 / 2)) x :=
      HasDerivAt.comp (h := fun y : ℝ => (1 + y) / 2) x hmid hd
    convert hD2 using 1
    field_simp
    congr 1
    ring
  have hconst : HasDerivAt (fun y : ℝ => Real.pi ^ 2 / 6 - Real.log 2 ^ 2) 0 x := by
    exact hasDerivAt_const x (Real.pi ^ 2 / 6 - Real.log 2 ^ 2)
  have h1 : HasDerivAt (fun y : ℝ => Real.log (1 - y) ^ 2 / 2)
      (-Real.log (1 - x) / (1 - x)) x := by
    convert (ha.pow 2).div_const 2 using 1
    all_goals norm_num
    all_goals field_simp [h1mx.ne', h1px.ne']
    all_goals ring_nf
  have h2 : HasDerivAt (fun y : ℝ => Real.log (1 + y) ^ 2)
      (2 * Real.log (1 + x) / (1 + x)) x := by
    convert hb.pow 2 using 1
    all_goals norm_num
    all_goals field_simp [h1mx.ne', h1px.ne']
    all_goals ring_nf
  have h3 : HasDerivAt (fun y : ℝ => 2 * Real.log (1 - y) * Real.log (1 + y))
      (2 * (-Real.log (1 + x) / (1 - x) + Real.log (1 - x) / (1 + x))) x := by
    have hm := ha.mul hb
    convert hm.const_mul 2 using 1
    all_goals norm_num
    all_goals field_simp [h1mx.ne', h1px.ne']
    all_goals ring_nf
  have h4 : HasDerivAt (fun y : ℝ => 2 * Real.log 2 * Real.log (1 + y))
      (2 * Real.log 2 / (1 + x)) x := by
    convert hb.const_mul (2 * Real.log 2) using 1
    all_goals norm_num
    all_goals field_simp [h1mx.ne', h1px.ne']
    all_goals ring_nf
  have h5 : HasDerivAt (fun y : ℝ => -2 * dilog ((1 + y) / 2))
      (2 * Real.log ((1 - x) / 2) / (1 + x)) x := by
    convert hD.const_mul (-2) using 1
    all_goals norm_num
    all_goals field_simp [h1mx.ne', h1px.ne']
    all_goals ring_nf
  have hM : HasDerivAt (fun y : ℝ => quadAltMclosed y)
      (-Real.log (1 - x) / (1 - x) +
        2 * Real.log (1 + x) / (1 + x) +
        2 * (-Real.log (1 + x) / (1 - x) + Real.log (1 - x) / (1 + x)) +
        2 * Real.log 2 / (1 + x) +
        0 +
        2 * Real.log ((1 - x) / 2) / (1 + x)) x := by
    unfold quadAltMclosed
    convert ((((h1.add h2).add h3).add h4).add hconst).add h5 using 2
    all_goals simp
    all_goals field_simp [h1mx.ne', h1px.ne']
    all_goals ring_nf
  convert hM using 1
  unfold quadAltFclosed
  norm_num
  rw [show Real.log ((1 - x) / 2) = Real.log (1 - x) - Real.log 2 by
    rw [Real.log_div]
    · exact h1mx.ne'
    · norm_num]
  ring_nf


/-- `M` series for `−1 < y < 0` (termwise integration on `[y,0]`; the
positive-only `quadAltM_hasSum` covers `0 < y < 1`). -/
theorem quadAltM_hasSum_neg {y : ℝ} (hy0 : -1 < y) (hy1 : y < 0) :
    HasSum (fun k : ℕ => (1 + 2 * (-1 : ℝ) ^ (k + 1)) * parityRemainder24 k *
        y ^ (k + 1) / (k + 1 : ℝ))
      (quadAltMclosed y) := by
  have hyabs : |y| < 1 := by
    rw [abs_of_neg hy1]
    linarith
  have hle : y ≤ (0 : ℝ) := le_of_lt hy1
  -- 1. integrability on Ioc y 0
  have hint : ∀ k : ℕ, MeasureTheory.Integrable (quadAltMMoment k)
      (MeasureTheory.volume.restrict (Set.Ioc y (0 : ℝ))) := by
    intro k
    have hcont : ContinuousOn (quadAltMMoment k) (Set.uIcc y (0 : ℝ)) := by
      unfold quadAltMMoment
      exact (continuousOn_pow k).const_mul
        ((1 + 2 * (-1 : ℝ) ^ (k + 1)) * parityRemainder24 k)
    have hI : MeasureTheory.IntegrableOn (quadAltMMoment k)
        (Set.uIcc y (0 : ℝ)) MeasureTheory.volume :=
      hcont.integrableOn_compact isCompact_uIcc
    exact hI.mono_set (by
      intro t ht
      rcases (Set.mem_Ioc.mp ht) with ⟨hyt, ht0⟩
      simp [Set.uIcc, min_eq_left hle, max_eq_right hle, hyt.le, ht0])
  -- 2. norm summability: |t| ≤ |y| on (y, 0]
  have hnorm_sum : Summable (fun k : ℕ =>
      ∫ t : ℝ in Set.Ioc y (0 : ℝ), ‖quadAltMMoment k t‖) := by
    have hH (k : ℕ) : (harmonicNumber (k + 1) : ℝ) ≤ (k + 1 : ℝ) := by
      have hsum : (∑ j ∈ Finset.range (k + 1), (1 : ℝ) / (j + 1 : ℝ)) ≤ (k + 1 : ℝ) := by
        calc
          (∑ j ∈ Finset.range (k + 1), (1 : ℝ) / (j + 1 : ℝ))
              ≤ ∑ _j ∈ Finset.range (k + 1), (1 : ℝ) := by
                exact Finset.sum_le_sum (fun j hj => by
                  rw [one_div]
                  have hj1 : (1 : ℝ) ≤ (j + 1 : ℝ) := by exact_mod_cast (by omega : 1 ≤ j + 1)
                  have hjpos : (0 : ℝ) < (j + 1 : ℝ) := by exact_mod_cast (by omega : 0 < j + 1)
                  exact (inv_le_one₀ hjpos).mpr hj1)
          _ = (k + 1 : ℝ) := by simp
      simpa [harmonicNumber] using hsum
    have hyabs0 : 0 ≤ |y| := abs_nonneg y
    have hgeom : Summable (fun k : ℕ => 9 * (harmonicNumber (k + 1) : ℝ) * |y| ^ k) := by
      have hxnorm : ‖|y|‖ < 1 := by simpa using hyabs
      have hk : Summable (fun k : ℕ => (k : ℝ) * |y| ^ k) := by
        simpa using (summable_pow_mul_geometric_of_norm_lt_one 1 (r := |y|) hxnorm)
      have h1 : Summable (fun k : ℕ => |y| ^ k) :=
        summable_geometric_of_lt_one hyabs0 hyabs
      have hsum : Summable (fun k : ℕ => (k + 1 : ℝ) * |y| ^ k) := by
        simpa [add_mul, one_mul, Nat.cast_add, Nat.cast_one] using
          (hk.add h1).congr (fun k => by ring)
      refine (hsum.mul_left 9).of_nonneg_of_le ?_ ?_
      · intro k
        exact mul_nonneg (mul_nonneg (by norm_num) (harmonicNumber_nonneg (k + 1)))
          (pow_nonneg hyabs0 k)
      · intro k
        have hyk : 0 ≤ |y| ^ k := pow_nonneg hyabs0 k
        have hHk : harmonicNumber (k + 1) ≤ (k + 1 : ℝ) := hH k
        have hmid : harmonicNumber (k + 1) * |y| ^ k ≤ (k + 1 : ℝ) * |y| ^ k :=
          mul_le_mul_of_nonneg_right hHk hyk
        have hres : 9 * (harmonicNumber (k + 1) * |y| ^ k) ≤ 9 * ((k + 1 : ℝ) * |y| ^ k) :=
          mul_le_mul_of_nonneg_left hmid (by norm_num)
        simpa [mul_assoc] using hres
    refine hgeom.of_nonneg_of_le ?_ ?_
    · intro k
      exact MeasureTheory.setIntegral_nonneg measurableSet_Ioc (fun t ht => norm_nonneg _)
    · intro k
      -- ∫_{Ioc y 0} ‖moment‖ ≤ 9H(k+1)|y|^k·(−y) ≤ 9H(k+1)|y|^k（−y ≤ 1）
      have hb (t : ℝ) (ht : t ∈ Set.Ioc y (0 : ℝ)) :
          ‖quadAltMMoment k t‖ ≤ 9 * (harmonicNumber (k + 1) : ℝ) * |y| ^ k := by
        rcases (Set.mem_Ioc.mp ht) with ⟨hyt, ht0⟩
        have htabs : |t| ≤ |y| := by
          rw [abs_of_nonpos ht0, abs_of_neg hy1]
          linarith
        have htnonneg : 0 ≤ |t| := abs_nonneg t
        unfold quadAltMMoment
        rw [Real.norm_eq_abs, abs_mul, abs_mul, abs_pow]
        calc
          |1 + 2 * (-1 : ℝ) ^ (k + 1)| * |parityRemainder24 k| * |t| ^ k
              ≤ 3 * (3 * harmonicNumber k) * |t| ^ k := by
                gcongr
                · exact quadAltC_abs_le k
                · exact abs_parityRemainder24_le k
          _ = 9 * (harmonicNumber k : ℝ) * |t| ^ k := by ring
          _ ≤ 9 * (harmonicNumber (k + 1) : ℝ) * |y| ^ k := by
            have h1 : 9 * harmonicNumber k ≤ 9 * harmonicNumber (k + 1) := by
              exact mul_le_mul_of_nonneg_left (harmonicNumber_mono (Nat.le_succ k)) (by positivity)
            have h2 : |t| ^ k ≤ |y| ^ k := pow_le_pow_left₀ htnonneg htabs k
            have h3 : 0 ≤ 9 * harmonicNumber (k + 1) := by
              exact mul_nonneg (by norm_num) (harmonicNumber_nonneg (k + 1))
            exact (mul_le_mul_of_nonneg_right h1 (pow_nonneg htnonneg k)).trans
              (mul_le_mul_of_nonneg_left h2 h3)
      -- ∫_{Ioc y 0} 常数 9H|y|^k = 9H|y|^k·(−y)
      calc
        (∫ t : ℝ in Set.Ioc y (0 : ℝ), ‖quadAltMMoment k t‖)
            ≤ ∫ _t : ℝ in Set.Ioc y (0 : ℝ), 9 * (harmonicNumber (k + 1) : ℝ) * |y| ^ k := by
              have hI_abs : MeasureTheory.IntegrableOn (fun t : ℝ => ‖quadAltMMoment k t‖)
                  (Set.Ioc y (0 : ℝ)) MeasureTheory.volume := by
                have hcont_abs : ContinuousOn (fun t : ℝ => ‖quadAltMMoment k t‖)
                    (Set.uIcc y (0 : ℝ)) := by
                  unfold quadAltMMoment
                  exact ((continuousOn_pow k).const_mul
                    ((1 + 2 * (-1 : ℝ) ^ (k + 1)) * parityRemainder24 k)).norm
                have hI : MeasureTheory.IntegrableOn (fun t : ℝ => ‖quadAltMMoment k t‖)
                    (Set.uIcc y (0 : ℝ)) MeasureTheory.volume :=
                  hcont_abs.integrableOn_compact isCompact_uIcc
                exact hI.mono_set (by
                  intro t ht
                  rcases (Set.mem_Ioc.mp ht) with ⟨hyt, ht0⟩
                  simp [Set.uIcc, min_eq_left hle, max_eq_right hle, hyt.le, ht0])
              have hgI : MeasureTheory.IntegrableOn
                  (fun _t : ℝ => 9 * (harmonicNumber (k + 1) : ℝ) * |y| ^ k)
                  (Set.Ioc y (0 : ℝ)) MeasureTheory.volume := by
                have hc : ContinuousOn
                    (fun _t : ℝ => 9 * (harmonicNumber (k + 1) : ℝ) * |y| ^ k)
                    (Set.uIcc y (0 : ℝ)) := continuousOn_const
                have hI : MeasureTheory.IntegrableOn
                    (fun _t : ℝ => 9 * (harmonicNumber (k + 1) : ℝ) * |y| ^ k)
                    (Set.uIcc y (0 : ℝ)) MeasureTheory.volume :=
                  hc.integrableOn_compact isCompact_uIcc
                exact hI.mono_set (by
                  intro t ht
                  rcases (Set.mem_Ioc.mp ht) with ⟨hyt, ht0⟩
                  simp [Set.uIcc, min_eq_left hle, max_eq_right hle, hyt.le, ht0])
              exact MeasureTheory.setIntegral_mono_on hI_abs hgI measurableSet_Ioc (by
                intro t ht
                exact hb t ht)
        _ = 9 * (harmonicNumber (k + 1) : ℝ) * |y| ^ k * (-y) := by
          rw [MeasureTheory.setIntegral_const]
          have hvol : MeasureTheory.volume.real (Set.Ioc y (0 : ℝ)) = -y := by
            change (MeasureTheory.volume (Set.Ioc y (0 : ℝ))).toReal = -y
            rw [Real.volume_Ioc]
            simp [ENNReal.toReal_ofReal, hy1.le]
          rw [hvol]
          simp only [smul_eq_mul]
          ring
        _ ≤ 9 * (harmonicNumber (k + 1) : ℝ) * |y| ^ k := by
          -- −y ≤ 1（y > −1）
          have hneg : -y ≤ 1 := by linarith
          have hnn : 0 ≤ 9 * (harmonicNumber (k + 1) : ℝ) * |y| ^ k := by
            exact mul_nonneg (mul_nonneg (by norm_num) (harmonicNumber_nonneg (k + 1)))
              (pow_nonneg hyabs0 k)
          calc
            9 * harmonicNumber (k + 1) * |y| ^ k * -y = (-y) * (9 * harmonicNumber (k + 1) * |y| ^ k) := by ring
            _ ≤ 1 * (9 * harmonicNumber (k + 1) * |y| ^ k) := by
              exact mul_le_mul_of_nonneg_right hneg hnn
            _ = 9 * harmonicNumber (k + 1) * |y| ^ k := by ring
  -- 3. hasSum_integral
  have hsum := MeasureTheory.hasSum_integral_of_summable_integral_norm
    (μ := MeasureTheory.volume.restrict (Set.Ioc y (0 : ℝ))) hint hnorm_sum
  -- 4. moment integral = − M-term
  have hIoc_pow (k : ℕ) : (∫ t : ℝ in Set.Ioc y (0 : ℝ), t ^ k) = -y ^ (k + 1) / (k + 1 : ℝ) := by
    calc
      (∫ t : ℝ in Set.Ioc y (0 : ℝ), t ^ k) = -∫ t : ℝ in (0 : ℝ)..y, t ^ k := by
        rw [intervalIntegral.integral_of_ge hle]
        ring
      _ = -y ^ (k + 1) / (k + 1 : ℝ) := by
        rw [integral_pow]
        ring
  have hterm : ∀ k : ℕ, (∫ t : ℝ in Set.Ioc y (0 : ℝ), quadAltMMoment k t)
      = -((1 + 2 * (-1 : ℝ) ^ (k + 1)) * parityRemainder24 k * y ^ (k + 1) / (k + 1 : ℝ)) := by
    intro k
    unfold quadAltMMoment
    rw [show (fun t : ℝ => (1 + 2 * (-1 : ℝ) ^ (k + 1)) * parityRemainder24 k * t ^ k) =
        fun t => ((1 + 2 * (-1 : ℝ) ^ (k + 1)) * parityRemainder24 k) * t ^ k by
          funext t; ring,
      MeasureTheory.integral_const_mul, hIoc_pow k]
    ring
  -- 5. limit: ∫_{Ioc y 0} tsum = −Mclosed y
  have hlim : (∫ t : ℝ in Set.Ioc y (0 : ℝ), ∑' k : ℕ, quadAltMMoment k t)
      = -quadAltMclosed y := by
    -- 桥：∫_{Ioc y 0} g = −∫₀₍y₎ g = ∫₍y₎₀ g
    have hftc : (∫ t : ℝ in (0 : ℝ)..y, (quadAltFclosed t - 2 * quadAltFclosed (-t))) = quadAltMclosed y := by
      -- FTC on [y,0]：∫₍y₎₀ f' = f 0 − f y = −Mclosed y；∫₀ʸ = −∫₍y₎₀ = Mclosed y
      have hcontM : ContinuousOn quadAltMclosed (Set.Icc y (0 : ℝ)) := by
        intro t ht
        have htabs : |t| < 1 := by
          rw [abs_lt]
          rcases Set.mem_Icc.mp ht with ⟨hty, ht0⟩
          constructor <;> linarith
        exact ((quadAltMclosed_hasDerivAt' htabs).continuousAt).continuousWithinAt
      have hderivM : ∀ t ∈ Set.Ioo y (0 : ℝ),
          HasDerivAt quadAltMclosed (quadAltFclosed t - 2 * quadAltFclosed (-t)) t := by
        intro t ht
        rcases (Set.mem_Ioo.mp ht) with ⟨hyt, ht1⟩
        exact quadAltMclosed_hasDerivAt' (by
          rw [abs_lt]
          constructor <;> linarith)
      have hintM : IntervalIntegrable
          (fun t : ℝ => quadAltFclosed t - 2 * quadAltFclosed (-t)) MeasureTheory.volume y 0 := by
        apply ContinuousOn.intervalIntegrable
        rw [Set.uIcc_of_le (le_of_lt hy1)]
        intro t ht
        rcases Set.mem_Icc.mp ht with ⟨hty, ht0⟩
        have htabs : |t| < 1 := by rw [abs_lt]; constructor <;> linarith
        have hF : ContinuousAt quadAltFclosed t := by
          unfold quadAltFclosed
          have h1 : (1 : ℝ) - t ≠ 0 := by linarith
          have h2 : (0 : ℝ) < 1 + t := by linarith
          fun_prop (disch := first | linarith | positivity)
        have hFneg : ContinuousAt (fun s : ℝ => quadAltFclosed (-s)) t := by
          unfold quadAltFclosed
          have h1 : (1 : ℝ) + t ≠ 0 := by linarith
          have h2 : (0 : ℝ) < 1 - t := by linarith
          fun_prop (disch := first | linarith | positivity)
        exact ((hF.sub (hFneg.const_smul (2:ℝ))).continuousWithinAt).congr
          (fun s _ => by simp [smul_eq_mul]) (by simp [smul_eq_mul])
      have hfund := intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le
        (a := y) (b := 0) (f := quadAltMclosed)
        (f' := fun t : ℝ => quadAltFclosed t - 2 * quadAltFclosed (-t))
        hle hcontM hderivM hintM
      rw [quadAltMclosed_zero] at hfund
      -- hfund : ∫₍y₎₀ f' = −Mclosed y；∫₀₍y₎ = −∫₍y₎₀ = Mclosed y
      rw [intervalIntegral.integral_symm y 0, hfund]
      ring
    calc
      (∫ t : ℝ in Set.Ioc y (0 : ℝ), ∑' k : ℕ, quadAltMMoment k t)
          = -∫ t : ℝ in (0 : ℝ)..y, (∑' k : ℕ, quadAltMMoment k t) := by
            rw [intervalIntegral.integral_of_ge hle]
            ring
      _ = -quadAltMclosed y := by
        have hcon : (∫ t : ℝ in (0 : ℝ)..y, ∑' k : ℕ, quadAltMMoment k t)
            = (∫ t : ℝ in (0 : ℝ)..y, (quadAltFclosed t - 2 * quadAltFclosed (-t))) := by
          apply intervalIntegral.integral_congr_ae
          filter_upwards with t ht
          have htabs : |t| < 1 := by
            rw [abs_lt]
            rcases (Set.mem_uIoc.mp ht) with ⟨hyt, ht0⟩ | ⟨hyt, ht0⟩
            · exfalso
              linarith
            · constructor <;> linarith
          by_cases ht0 : t = 0
          · subst t
            have hterm0 : ∀ k : ℕ, quadAltMMoment k 0 = 0 := by
              intro k
              unfold quadAltMMoment
              by_cases hk : k = 0
              · subst hk
                simp [quadAltP_zero]
              · have hk' : (0 : ℝ) ^ k = 0 := zero_pow (by omega : k ≠ 0)
                rw [hk']
                ring
            have htsum : (∑' k : ℕ, quadAltMMoment k 0) = 0 := by
              simp [hterm0]
            rw [htsum]
            simp [quadAltFclosed]
          · exact (quadAltC_hasSum htabs ht0).tsum_eq
        rw [hcon, hftc]
  -- 6. assemble
  have hmain : HasSum (fun k : ℕ => -((1 + 2 * (-1 : ℝ) ^ (k + 1)) * parityRemainder24 k *
        y ^ (k + 1) / (k + 1 : ℝ))) (-quadAltMclosed y) := by
    simpa only [hterm, hlim] using hsum
  convert hmain.neg using 1
  · funext k
    ring
  · ring


/-- `M` series for all `|y| < 1`, `y ≠ 0` (general version, by sign split). -/
theorem quadAltM_hasSum' {y : ℝ} (hy : |y| < 1) (hyne : y ≠ 0) :
    HasSum (fun k : ℕ => (1 + 2 * (-1 : ℝ) ^ (k + 1)) * parityRemainder24 k *
        y ^ (k + 1) / (k + 1 : ℝ))
      (quadAltMclosed y) := by
  rcases lt_or_gt_of_ne hyne with hyneg | hypos
  · have hy0 : -1 < y := (abs_lt.mp hy).1
    exact quadAltM_hasSum_neg hy0 hyneg
  · have hy1 : y < 1 := (abs_lt.mp hy).2
    exact quadAltM_hasSum hypos hy1

/-- `Σ_{k≥1} c_{k+1}P(k)y^k/(k+1) = Mclosed(y)/y` for `|y| < 1`, `y ≠ 0`
(general version of `quadAltMPart_hasSum`). -/
theorem quadAltMPart_hasSum' {y : ℝ} (hy : |y| < 1) (hyne : y ≠ 0) :
    HasSum (fun m : ℕ => (1 + 2 * (-1 : ℝ) ^ (m + 2)) * parityRemainder24 (m + 1) *
        y ^ (m + 1) / (m + 2 : ℝ))
      (quadAltMclosed y / y) := by
  have H1 := quadAltM_hasSum' hy hyne
  have hf0 : (1 + 2 * (-1 : ℝ) ^ (0 + 1)) * parityRemainder24 0 * y ^ (0 + 1) / (0 + 1 : ℝ) = 0 := by
    simp [quadAltP_zero]
  have H1a : HasSum (fun n : ℕ => (1 + 2 * (-1 : ℝ) ^ (n + 1 + 1)) * parityRemainder24 (n + 1) *
        y ^ (n + 1 + 1) / ((n + 1 : ℝ) + 1))
      (quadAltMclosed y - (1 + 2 * (-1 : ℝ) ^ (0 + 1)) * parityRemainder24 0 *
        y ^ (0 + 1) / (0 + 1 : ℝ)) := by
    simpa [hf0] using
      (hasSum_nat_add_iff' (f := fun k : ℕ =>
        (1 + 2 * (-1 : ℝ) ^ (k + 1)) * parityRemainder24 k * y ^ (k + 1) / (k + 1 : ℝ))
        (1 : ℕ)).mpr H1
  have H2 := H1a.mul_left (1 / y)
  have hterm : ∀ n : ℕ, (1 + 2 * (-1 : ℝ) ^ (n + 2)) * parityRemainder24 (n + 1) *
        y ^ (n + 1) / (n + 2 : ℝ) =
      1 / y * ((1 + 2 * (-1 : ℝ) ^ (n + 1 + 1)) * parityRemainder24 (n + 1) *
        y ^ (n + 1 + 1) / ((n + 1 : ℝ) + 1)) := by
    intro n
    field_simp [hyne]
    ring
  have H3 : HasSum (fun n : ℕ => (1 + 2 * (-1 : ℝ) ^ (n + 2)) * parityRemainder24 (n + 1) *
        y ^ (n + 1) / (n + 2 : ℝ))
      (1 / y * (quadAltMclosed y - (1 + 2 * (-1 : ℝ) ^ (0 + 1)) * parityRemainder24 0 *
        y ^ (0 + 1) / (0 + 1 : ℝ))) :=
    H2.congr_fun (fun n => hterm n)
  convert H3 using 1
  field_simp [hyne]
  simp [quadAltP_zero]


/-- Increment series: `Σ_{k≥0} (S_{k+1}−S_k)y^k = 2J(y)/y` (Q6047 (2.2)–(2.5)). -/
theorem quadAltIncrement_hasSum {y : ℝ} (hy : |y| < 1) (hyne : y ≠ 0) :
    HasSum (fun k : ℕ => quadAltIncrement k * y ^ k) (2 * quadAltJclosed y / y) := by
  have h0 : quadAltIncrement 0 = 0 := by
    unfold quadAltIncrement
    simp [quadAltS_zero]
  have hP := quadAltMPart_hasSum' hy hyne
  have hD := quadAltDiag_hasSum hy hyne
  have hmain : HasSum (fun j : ℕ =>
        2 * ((1 + 2 * (-1 : ℝ) ^ (j + 2)) * parityRemainder24 (j + 1) * y ^ (j + 1) / (j + 2 : ℝ)) +
        4 * (1 + (-1 : ℝ) ^ (j + 2)) * y ^ (j + 1) / (↑(j + 2) : ℝ) ^ 2)
      (2 * (quadAltMclosed y / y) + 4 * (dilog y + dilog (-y)) / y) :=
    (hP.mul_left 2).add hD
  have hterm (j : ℕ) :
      2 * ((1 + 2 * (-1 : ℝ) ^ (j + 2)) * parityRemainder24 (j + 1) * y ^ (j + 1) / (j + 2 : ℝ)) +
        4 * (1 + (-1 : ℝ) ^ (j + 2)) * y ^ (j + 1) / (↑(j + 2) : ℝ) ^ 2
      = quadAltIncrement (j + 1) * y ^ (j + 1) := by
    unfold quadAltIncrement
    have hj : j + 1 ≠ 0 := by omega
    simp [hj]
    rw [quadAltS_succ_sub]
    rw [quadAltC_sq_sub_one]
    ring
  have hval : 2 * (quadAltMclosed y / y) + 4 * (dilog y + dilog (-y)) / y = 2 * quadAltJclosed y / y := by
    unfold quadAltJclosed
    have hd := quadAlt_dilog_add_neg' hy
    rw [hd]
    ring
  have hshift : HasSum (fun k : ℕ => quadAltIncrement (k + 1) * y ^ (k + 1))
      (2 * quadAltJclosed y / y) := by
    convert hmain using 1
    · funext j
      exact (hterm j).symm
    · exact hval.symm
  refine (hasSum_nat_add_iff' (f := fun k : ℕ => quadAltIncrement k * y ^ k)
    (g := 2 * quadAltJclosed y / y) (1 : ℕ)).mp ?_
  simpa [h0] using hshift

/-- `S`-series: `Σ_{n≥0} S_{n+1} y^n = Qclosed(y)/y` (Q6047 (2.10), Cauchy product
of the increment series with the geometric series). -/
theorem quadAltS_generating_hasSum {y : ℝ} (hy : |y| < 1) (hyne : y ≠ 0) :
    HasSum (fun n : ℕ => quadAltS n * y ^ n) (quadAltQclosed y / y) := by
  have hInc := quadAltIncrement_hasSum hy hyne
  have hf : Summable (fun k : ℕ => ‖quadAltIncrement k * y ^ k‖) :=
    hInc.summable.norm
  have hg : Summable (fun m : ℕ => ‖y ^ m‖) := by
    simpa [Real.norm_eq_abs, abs_pow, abs_of_nonneg (abs_nonneg y)] using
      (summable_geometric_of_lt_one (abs_nonneg y) hy)
  have hproduct := hasSum_sum_range_mul_of_summable_norm
    (f := fun k : ℕ => quadAltIncrement k * y ^ k) (g := fun m : ℕ => y ^ m) hf hg
  have hcoeff (n : ℕ) : (∑ k ∈ Finset.range (n + 1), quadAltIncrement k * y ^ k * y ^ (n - k))
      = quadAltS n * y ^ n := by
    calc
      (∑ k ∈ Finset.range (n + 1), quadAltIncrement k * y ^ k * y ^ (n - k))
          = (∑ k ∈ Finset.range (n + 1), quadAltIncrement k * y ^ n) := by
            apply Finset.sum_congr rfl
            intro k hk
            have hkn : k < n + 1 := Finset.mem_range.mp hk
            have hpow : y ^ k * y ^ (n - k) = y ^ n := by
              rw [← pow_add]
              congr 1
              omega
            rw [mul_assoc, hpow]
      _ = (∑ k ∈ Finset.range (n + 1), quadAltIncrement k) * y ^ n := by
            rw [Finset.sum_mul]
      _ = quadAltS n * y ^ n := by
        rw [quadAltS_eq_sum_increment]
  have hval : (∑' k : ℕ, quadAltIncrement k * y ^ k) * (∑' m : ℕ, y ^ m) = quadAltQclosed y / y := by
    have hgeom : (∑' m : ℕ, y ^ m) = 1 / (1 - y) := by
      have hs := hasSum_geometric_of_norm_lt_one (show ‖y‖ < 1 by simpa using hy)
      simpa using hs.tsum_eq
    rw [hInc.tsum_eq, hgeom]
    unfold quadAltQclosed
    field_simp [hyne]
  convert hproduct using 1
  · funext n
    exact (hcoeff n).symm
  · exact hval.symm


/-! ## Layer C main: coefficient integration (Q6047 §3.3) -/

/-- Coefficient moment: `(−log x)/x · S_{n+1}(−x)^{n+1}` on `[0,1]`. -/
def quadAltCoeffMoment (n : ℕ) (x : ℝ) : ℝ :=
  (-Real.log x) / x * quadAltS n * (-x) ^ (n + 1)

/-- Pointwise HasSum: `Σ moment n x = (−log x)/x · Q(−x)` for `0 < x < 1`. -/
theorem quadAltCoeffMoment_hasSum_pointwise {x : ℝ} (hx0 : 0 < x) (hx1 : x < 1) :
    HasSum (fun n : ℕ => quadAltCoeffMoment n x)
      ((-Real.log x) / x * quadAltQclosed (-x)) := by
  have hS := quadAltS_generating_hasSum (y := -x)
    (by simpa [abs_of_pos hx0] using hx1) (by linarith)
  -- hS : HasSum (fun n => S_{n+1}·(−x)^n) (Qclosed(−x)/(−x))
  -- moment = (−log x)/x·(−x)·S·(−x)^n
  have hmul := hS.mul_left ((-Real.log x) / x * (-x))
  convert hmul using 1
  · funext n
    unfold quadAltCoeffMoment
    rw [pow_succ]
    field_simp [ne_of_gt hx0]
  · field_simp [ne_of_gt hx0]

/-- Integral of the coefficient moment over `[0,1]` = the alternating quadratic term. -/
theorem quadAltCoeffMoment_integral (n : ℕ) :
    (∫ x : ℝ in (0 : ℝ)..1, quadAltCoeffMoment n x) =
      alternatingQuadraticEulerTerm24 n := by
  unfold quadAltCoeffMoment alternatingQuadraticEulerTerm24 quadraticEulerTerm24
  -- ∫₀¹ (−log x)/x·S·(−x)^{n+1} = S·(−1)^{n+1}·∫₀¹ (−log x)x^n
  rw [show (fun x : ℝ => (-Real.log x) / x * quadAltS n * (-x) ^ (n + 1)) =
      fun x => quadAltS n * (-1 : ℝ) ^ (n + 1) * (x ^ n * (-Real.log x)) by
        funext x
        rw [neg_pow]
        by_cases hx : x = 0
        · subst x
          simp
        · field_simp [hx]
          ring]
  rw [intervalIntegral.integral_const_mul]
  -- ∫₀¹ x^n·(−log x) = 1/(n+1)²
  have hlog : (∫ x : ℝ in (0 : ℝ)..1, x ^ n * (-Real.log x)) = 1 / (n + 1 : ℝ) ^ 2 := by
    have h := RamanujanChallenge.P26.integral_pow_mul_log26 n
    -- h : ∫₀¹ x^n·log x = −1/(n+1)²
    have hneg : (∫ x : ℝ in (0 : ℝ)..1, x ^ n * (-Real.log x)) =
        -∫ x : ℝ in (0 : ℝ)..1, x ^ n * Real.log x := by
      rw [← intervalIntegral.integral_neg]
      congr 1
      funext x
      ring
    rw [hneg, h]
    ring
  rw [hlog]
  unfold quadAltS
  ring


/-- Layer C main: `tsum(alternatingQuadraticEulerTerm24) = ∫₀¹ (−log x)/x · Q(−x) dx`
(the coefficient integration, Q6047 §3.3). -/
theorem quadAlt_tsum_eq_coeff_integral :
    (∑' n : ℕ, alternatingQuadraticEulerTerm24 n) =
      ∫ x : ℝ in (0 : ℝ)..1, (-Real.log x) / x * quadAltQclosed (-x) := by
  have hint : ∀ n : ℕ, MeasureTheory.Integrable (quadAltCoeffMoment n)
      (MeasureTheory.volume.restrict (Set.Ioc (0 : ℝ) 1)) := by
    intro n
    have hlogI : IntervalIntegrable (fun x : ℝ => x ^ n * Real.log x)
        MeasureTheory.volume 0 1 :=
      (intervalIntegral.intervalIntegrable_log').continuousOn_mul (continuousOn_pow n)
    have hlogI' : IntervalIntegrable (fun x : ℝ => (-Real.log x) * x ^ n)
        MeasureTheory.volume 0 1 := by
      simpa [mul_comm] using hlogI.neg
    have hEq (x : ℝ) : quadAltCoeffMoment n x =
        quadAltS n * (-1 : ℝ) ^ (n + 1) * ((-Real.log x) * x ^ n) := by
      unfold quadAltCoeffMoment
      rw [neg_pow]
      by_cases hx : x = 0
      · subst x
        simp
      · field_simp [hx]
        ring
    have hconst := (hlogI'.const_mul (quadAltS n * (-1 : ℝ) ^ (n + 1)))
    have hI : IntervalIntegrable (quadAltCoeffMoment n) MeasureTheory.volume 0 1 := by
      convert hconst using 1
      funext x
      rw [hEq]
    -- IntervalIntegrable 0 1 = IntegrableOn (Ioc 0 1) ∧ …；第一分量即所需
    exact hI.1.integrable
  have hnorm_sum : Summable (fun n : ℕ =>
      ∫ t : ℝ in Set.Ioc (0 : ℝ) 1, ‖quadAltCoeffMoment n t‖) := by
    -- ∫|moment| = |alternatingQuadraticEulerTerm24|（用 ∫₀¹ (−log x)x^n = 1/(n+1)² 的精确值）
    have hnorm (n : ℕ) : (∫ t : ℝ in Set.Ioc (0 : ℝ) 1, ‖quadAltCoeffMoment n t‖)
        = |alternatingQuadraticEulerTerm24 n| := by
      have hEq (t : ℝ) (ht0 : 0 < t) : ‖quadAltCoeffMoment n t‖ =
          |quadAltS n| * ‖(-Real.log t) * t ^ n‖ := by
        unfold quadAltCoeffMoment
        rw [Real.norm_eq_abs, Real.norm_eq_abs]
        simp only [abs_mul, abs_div, abs_pow, abs_neg]
        rw [abs_of_pos ht0]
        field_simp [ne_of_gt ht0]
        ring
      have hlog (n : ℕ) : (∫ x : ℝ in (0 : ℝ)..1, (-Real.log x) * x ^ n) = 1 / (n + 1 : ℝ) ^ 2 := by
        have h := RamanujanChallenge.P26.integral_pow_mul_log26 n
        have hneg : (∫ x : ℝ in (0 : ℝ)..1, x ^ n * (-Real.log x)) =
            -∫ x : ℝ in (0 : ℝ)..1, x ^ n * Real.log x := by
          rw [← intervalIntegral.integral_neg]
          congr 1
          funext x
          ring
        have hswap : (∫ x : ℝ in (0 : ℝ)..1, (-Real.log x) * x ^ n)
            = ∫ x : ℝ in (0 : ℝ)..1, x ^ n * (-Real.log x) := by
          apply intervalIntegral.integral_congr_ae
          filter_upwards with x hx
          ring
        rw [hswap, hneg, h]
        ring
      -- ∫_{Ioc} |moment| = ∫₀¹ |S|·|(−log t)·t^n| = |S|·∫₀¹ (−log t)·t^n
      calc
        (∫ t : ℝ in Set.Ioc (0 : ℝ) 1, ‖quadAltCoeffMoment n t‖)
            = ∫ t : ℝ in (0 : ℝ)..1, ‖quadAltCoeffMoment n t‖ := by
              rw [intervalIntegral.integral_of_le (by norm_num : (0 : ℝ) ≤ 1)]
        _ = ∫ t : ℝ in (0 : ℝ)..1, |quadAltS n| * ‖(-Real.log t) * t ^ n‖ := by
              apply intervalIntegral.integral_congr_ae
              filter_upwards [MeasureTheory.Measure.ae_ne MeasureTheory.volume (1 : ℝ)] with t htne ht
              exact hEq t (by simpa using ht.1)
        _ = |quadAltS n| * (∫ t : ℝ in (0 : ℝ)..1, ‖(-Real.log t) * t ^ n‖) := by
              rw [intervalIntegral.integral_const_mul]
        _ = |quadAltS n| * (1 / (n + 1 : ℝ) ^ 2) := by
              congr 1
              -- ∫₀¹ |(−log t)·t^n| = ∫₀¹ (−log t)·t^n（t ∈ (0,1] 时 −log t ≥ 0）
              have hnonneg : ∀ t ∈ Set.Ioc (0 : ℝ) 1, 0 ≤ (-Real.log t) * t ^ n := by
                intro t ht
                rcases (Set.mem_Ioc.mp ht) with ⟨ht0, ht1⟩
                have hlogt : Real.log t ≤ 0 := Real.log_nonpos ht0.le ht1
                exact mul_nonneg (by linarith) (pow_nonneg ht0.le n)
              rw [← hlog]
              apply intervalIntegral.integral_congr_ae
              filter_upwards [MeasureTheory.Measure.ae_ne MeasureTheory.volume (1 : ℝ)] with t htne ht
              rw [Real.norm_eq_abs]
              congr 1
              rw [abs_of_nonneg (hnonneg t ⟨by simpa using ht.1, by simpa using ht.2⟩)]
        _ = |alternatingQuadraticEulerTerm24 n| := by
              unfold alternatingQuadraticEulerTerm24 quadraticEulerTerm24 quadAltS
              rw [abs_mul, abs_div, abs_pow]
              norm_num
              ring
    exact summable_quadAlt.norm.congr (fun n => (hnorm n).symm)
  have hsum := MeasureTheory.hasSum_integral_of_summable_integral_norm
    (μ := MeasureTheory.volume.restrict (Set.Ioc (0 : ℝ) 1)) hint hnorm_sum
  have hterm (n : ℕ) : (∫ t : ℝ in Set.Ioc (0 : ℝ) 1, quadAltCoeffMoment n t)
      = alternatingQuadraticEulerTerm24 n := by
    rw [← intervalIntegral.integral_of_le (by norm_num : (0 : ℝ) ≤ 1),
      quadAltCoeffMoment_integral n]
  have hlim : (∫ t : ℝ in Set.Ioc (0 : ℝ) 1, ∑' n : ℕ, quadAltCoeffMoment n t)
      = ∫ x : ℝ in (0 : ℝ)..1, (-Real.log x) / x * quadAltQclosed (-x) := by
    rw [← intervalIntegral.integral_of_le (by norm_num : (0 : ℝ) ≤ 1)]
    apply intervalIntegral.integral_congr_ae
    filter_upwards [MeasureTheory.Measure.ae_ne MeasureTheory.volume (1 : ℝ)] with x hxne hx
    have hx0 : 0 < x := by simpa using hx.1
    have hx1 : x < 1 := lt_of_le_of_ne (by simpa using hx.2) hxne
    exact (quadAltCoeffMoment_hasSum_pointwise hx0 hx1).tsum_eq
  have hsum' : HasSum (fun n : ℕ => ∫ t : ℝ in Set.Ioc (0 : ℝ) 1, quadAltCoeffMoment n t)
      (∫ x : ℝ in (0 : ℝ)..1, (-Real.log x) / x * quadAltQclosed (-x)) := by
    simpa [hlim] using hsum
  have hfinal : HasSum alternatingQuadraticEulerTerm24
      (∫ x : ℝ in (0 : ℝ)..1, (-Real.log x) / x * quadAltQclosed (-x)) := by
    convert hsum' using 1
    funext n
    exact (hterm n).symm
  rw [← hfinal.tsum_eq]

/-! ## Layer D bricks grafted from the sandbox agent (verified numerically) -/

/-- `W0'(t) = −2·log(t/(2−t))/t` for `0 < t < 1` (Q6047 (5.3)). -/
theorem quadAltW0_hasDerivAt {t : ℝ} (ht0 : 0 < t) (ht1 : t < 2) :
    HasDerivAt W0 (-2 * Real.log (t / (2 - t)) / t) t := by
  have hdlog : HasDerivAt (fun s : ℝ => dilog (s / 2)) (-Real.log (1 - t / 2) / (t / 2) * (1 / 2)) t := by
    have hmid0 : 0 < t / 2 := by positivity
    have hmid1 : t / 2 < 1 := by linarith
    have hd := dilog_hasDerivAt hmid0 hmid1
    have hhalf : HasDerivAt (fun s : ℝ => s / 2) (1 / 2) t := by
      simpa using (hasDerivAt_id t).div_const 2
    have hcomp := HasDerivAt.comp (h := fun s : ℝ => s / 2) t hd hhalf
    convert hcomp using 1
  have hlog2 : HasDerivAt (fun s : ℝ => Real.log (s / 2) ^ 2) (2 * Real.log (t / 2) * (1 / t)) t := by
    have hlogd : HasDerivAt (fun s : ℝ => Real.log (s / 2)) (1 / t) t := by
      have hd1 : HasDerivAt (fun s : ℝ => s / 2) (1 / 2) t := by
        simpa using (hasDerivAt_id t).div_const 2
      have hl := Real.hasDerivAt_log (ne_of_gt (by positivity : 0 < t / 2))
      have hcomp := HasDerivAt.comp (h := fun s : ℝ => s / 2) t hl hd1
      convert hcomp using 1
      field_simp
    have hp := hlogd.pow 2
    convert hp using 1
    all_goals norm_num
    all_goals field_simp
    all_goals ring
  unfold W0
  have hconst : HasDerivAt (fun _ : ℝ => Real.pi ^ 2 / 6) 0 t := hasDerivAt_const _ _
  have hmain : HasDerivAt (fun s : ℝ => Real.pi ^ 2 / 6 - 2 * dilog (s / 2) - Real.log (s / 2) ^ 2)
      (0 - 2 * (-Real.log (1 - t / 2) / (t / 2) * (1 / 2)) - 2 * Real.log (t / 2) * (1 / t)) t := by
    convert (hconst.sub (hdlog.const_mul 2)).sub hlog2 using 1
  have hlogcomb : Real.log (1 - t / 2) - Real.log (t / 2) = -Real.log (t / (2 - t)) := by
    rw [← Real.log_div (by linarith : 1 - t / 2 ≠ 0) (ne_of_gt (by positivity : 0 < t / 2))]
    have hdiv' : (1 - t / 2) / (t / 2) = (t / (2 - t))⁻¹ := by
      rw [← inv_div]
      field_simp [ne_of_gt ht0, ne_of_gt (by positivity : 0 < t / 2),
        ne_of_gt (by linarith : 0 < 2 - t)]
    rw [hdiv']
    exact Real.log_inv (t / (2 - t))
  have hD : 0 - 2 * (-Real.log (1 - t / 2) / (t / 2) * (1 / 2)) - 2 * Real.log (t / 2) * (1 / t)
      = -2 * Real.log (t / (2 - t)) / t := by
    have hA' : (t / 2)⁻¹ * (1 / 2) = t⁻¹ := by
      ring_nf
    have hlogA : (-Real.log (1 - t / 2)) / (t / 2) * (1 / 2) = (-Real.log (1 - t / 2)) / t := by
      calc
        (-Real.log (1 - t / 2)) / (t / 2) * (1 / 2) = (-Real.log (1 - t / 2)) * ((t / 2)⁻¹ * (1 / 2)) := by
          rw [div_eq_mul_inv]
          rw [mul_assoc]
        _ = (-Real.log (1 - t / 2)) * t⁻¹ := by rw [hA']
        _ = (-Real.log (1 - t / 2)) / t := by
          simp [div_eq_mul_inv]
    have hlogB : Real.log (t / 2) * (1 / t) = Real.log (t / 2) / t := by
      simp [div_eq_mul_inv]
    rw [hlogA]
    rw [mul_assoc]
    rw [hlogB]
    field_simp [ne_of_gt ht0]
    rw [show (2 - t) / 2 = 1 - t / 2 by ring]
    calc
      0 * t - -(2 * Real.log (1 - t / 2)) - 2 * Real.log (t / 2)
          = 2 * (Real.log (1 - t / 2) - Real.log (t / 2)) := by
            ring
      _ = -(2 * Real.log (t / (2 - t))) := by
            rw [hlogcomb]
            ring
  exact hmain.congr_deriv hD


/-- `V(x) = log(x)²/2 − log(x)·log(1+x) − Li2(−x) − Z2/2` (Q6047 (4.1)). -/
def quadAltV (x : ℝ) : ℝ :=
  Real.log x ^ 2 / 2 - Real.log x * Real.log (1 + x) - dilog (-x) - Real.pi ^ 2 / 12

/-- `V(1) = 0` (Q6047 (4.4), since `Li2(−1) = −Z2/2`). -/
theorem quadAltV_one : quadAltV 1 = 0 := by
  unfold quadAltV
  rw [show Real.log (1 : ℝ) = 0 by norm_num]
  have hneg : dilog (-1) = -Real.pi ^ 2 / 12 := by
    exact RamanujanChallenge.P26.dilog_neg_one26
  rw [hneg]
  ring


/-- `V'(x) = log(x)/(x(1+x))` for `0 < x < 1` (Q6047 (4.2)). -/
theorem quadAltV_hasDerivAt {x : ℝ} (hx0 : 0 < x) (hx1 : x < 1) :
    HasDerivAt quadAltV (Real.log x / (x * (1 + x))) x := by
  have hxne : x ≠ 0 := ne_of_gt hx0
  have h1px : 0 < 1 + x := by positivity
  -- log x 的导数
  have hl := Real.hasDerivAt_log hxne
  -- log(1+x) 的导数
  have hl1 : HasDerivAt (fun y : ℝ => Real.log (1 + y)) (1 / (1 + x)) x := by
    have hc : HasDerivAt (fun y : ℝ => 1 + y) 1 x := by
      convert (hasDerivAt_const x 1).add (hasDerivAt_id x) using 1
      simp
    have hlog := Real.hasDerivAt_log (ne_of_gt h1px)
    have hD : HasDerivAt (fun y : ℝ => Real.log (1 + y)) ((1 + x)⁻¹ * 1) x :=
      HasDerivAt.comp (h := fun y : ℝ => 1 + y) x hlog hc
    convert hD using 1
    field_simp
  -- dilog(−x) 的导数：dilog'(−x)·(−1) = −log(1+x)/x
  have hdneg : HasDerivAt (fun y : ℝ => dilog (-y)) (-Real.log (1 + x) / x) x := by
    have habs : |x| < 1 := by
      rw [abs_of_pos hx0]
      exact hx1
    have hd2a : HasDerivAt dilog (Real.log (1 + x) / x) (-x) := by
      simpa using dilog_hasDerivAt_of_abs_lt_one (by simpa using habs) (neg_ne_zero.mpr hxne)
    have hneg : HasDerivAt (fun y : ℝ => -y) (-1) x := by
      simpa using (hasDerivAt_id x).neg
    have hc := HasDerivAt.comp (h := fun y : ℝ => -y) x hd2a hneg
    convert hc using 1
    field_simp
  -- 组合：V = log²/2 − log·log(1+x) − dilog(−x) − π²/12
  have h1 : HasDerivAt (fun y : ℝ => Real.log y ^ 2 / 2) (Real.log x / x) x := by
    convert (hl.pow 2).div_const 2 using 1
    all_goals norm_num
    all_goals field_simp [hxne]
    all_goals ring
  have h2 : HasDerivAt (fun y : ℝ => Real.log y * Real.log (1 + y))
      (Real.log (1 + x) / x + Real.log x / (1 + x)) x := by
    convert (hl.mul hl1) using 1
    field_simp [hxne]
  have h3 : HasDerivAt (fun y : ℝ => dilog (-y)) (-Real.log (1 + x) / x) x := hdneg
  have hconst : HasDerivAt (fun _ : ℝ => Real.pi ^ 2 / 12) 0 x := hasDerivAt_const _ _
  unfold quadAltV
  have hmain : HasDerivAt (fun y : ℝ => Real.log y ^ 2 / 2 - Real.log y * Real.log (1 + y) -
        dilog (-y) - Real.pi ^ 2 / 12)
      (Real.log x / x - (Real.log (1 + x) / x + Real.log x / (1 + x)) -
        (-Real.log (1 + x) / x) - 0) x := by
    convert (((h1.sub h2).sub h3).sub hconst) using 1
  convert hmain using 1
  field_simp [hxne, ne_of_gt h1px]
  ring


/-- `d/dx J(−x) = −Dminus(x)` for `0 < x < 1` (Q6047 (4.8); note the
chain-rule sign: `J'(−x) = Dminus`, so `(J∘neg)'(x) = −Dminus`). -/
theorem quadAltJclosed_neg_hasDerivAt {x : ℝ} (hx0 : 0 < x) (hx1 : x < 1) :
    HasDerivAt (fun y : ℝ => quadAltJclosed (-y)) (-quadAltDminus x) x := by
  have hM := quadAltMclosed_hasDerivAt' (x := -x) (by simpa [abs_neg, abs_of_pos hx0] using hx1)
  -- hM : HasDerivAt (fun y => Mclosed y) (F(-x) − 2F(x)) (−x)
  have hMneg : HasDerivAt (fun y : ℝ => quadAltMclosed (-y))
      (-(quadAltFclosed (-x) - 2 * quadAltFclosed x)) x := by
    have hneg : HasDerivAt (fun y : ℝ => -y) (-1) x := by
      simpa using (hasDerivAt_id x).neg
    have hcomp := HasDerivAt.comp (h := fun y : ℝ => -y) x hM hneg
    convert hcomp using 1
    ring
  have ht2 : x ^ 2 < 1 := (sq_lt_one_iff_abs_lt_one x).mpr (by rw [abs_of_pos hx0]; exact hx1)
  have hd : HasDerivAt (fun y : ℝ => dilog (y ^ 2)) (-2 * Real.log (1 - x ^ 2) / x) x := by
    have hd2 := dilog_hasDerivAt (by positivity : 0 < x ^ 2) ht2
    have hpow2 : HasDerivAt (fun y : ℝ => y ^ 2) (2 * x) x := by
      simpa using (hasDerivAt_id x).pow 2
    have hcomp := HasDerivAt.comp (h := fun y : ℝ => y ^ 2) x hd2 hpow2
    convert hcomp using 1
    field_simp
  have hJ : HasDerivAt (fun y : ℝ => quadAltJclosed (-y))
      (-(quadAltFclosed (-x) - 2 * quadAltFclosed x) + -2 * Real.log (1 - x ^ 2) / x) x := by
    unfold quadAltJclosed
    convert (hMneg.add hd) using 1
    · funext y
      simp
  have hD' : -(quadAltFclosed (-x) - 2 * quadAltFclosed x) + -2 * Real.log (1 - x ^ 2) / x
      = -quadAltDminus x := by
    unfold quadAltDminus quadAltFclosed
    have hlog : Real.log (1 - x ^ 2) = Real.log (1 - x) + Real.log (1 + x) := by
      have hfac : 1 - x ^ 2 = (1 - x) * (1 + x) := by ring
      rw [hfac]
      rw [Real.log_mul]
      · exact ne_of_gt (by linarith : 0 < 1 - x)
      · exact ne_of_gt (by positivity : 0 < 1 + x)
    rw [hlog]
    rw [show 1 + -x = 1 - x by ring]
    simp [sub_neg_eq_add]
    rw [div_eq_mul_inv]
    ring
  exact hJ.congr_deriv hD'


/-- Pointwise bridge for (4.9): `(−log x)/x·Q(−x) = −2·V'(x)·J(−x)`
(`Q(−x) = 2J(−x)/(1+x)` and `V' = log x/(x(1+x))`). -/
theorem quadAltCoeffKernel_eq_VJ {x : ℝ} (hx0 : 0 < x) (hx1 : x < 1) :
    (-Real.log x) / x * quadAltQclosed (-x) =
      -2 * (Real.log x / (x * (1 + x))) * quadAltJclosed (-x) := by
  unfold quadAltQclosed
  field_simp [ne_of_gt hx0, ne_of_gt (by positivity : 0 < 1 + x)]
  ring


/-! ## Endpoint limit atoms for the IBP boundary conditions

The IBP on `[0,1]` has removable singularities at both ends: `V ~ (log x)²/2`
blows up at `0` but `J(-x) → 0` fast enough to kill it, and at `1` the blow-up of
`Dminus` is killed by `V 1 = 0`. These atoms make that precise. -/

theorem logSq_mul_self_tendsto :
    Tendsto (fun x : ℝ => Real.log x ^ 2 * x) (𝓝[>] (0:ℝ)) (𝓝 0) := by
  have h := tendsto_log_mul_rpow_nhdsGT_zero (r := (1:ℝ)/2) (by norm_num)
  have hsq := h.mul h
  simp only [mul_zero] at hsq
  refine hsq.congr' ?_
  filter_upwards [self_mem_nhdsWithin] with x hx
  have hx0 : (0:ℝ) < x := hx
  have hh : x ^ ((1:ℝ)/2) * x ^ ((1:ℝ)/2) = x := by
    rw [← Real.rpow_add hx0]; norm_num
  calc Real.log x * x ^ ((1:ℝ)/2) * (Real.log x * x ^ ((1:ℝ)/2))
      = Real.log x ^ 2 * (x ^ ((1:ℝ)/2) * x ^ ((1:ℝ)/2)) := by ring
    _ = Real.log x ^ 2 * x := by rw [hh]

theorem log_mul_self_tendsto :
    Tendsto (fun x : ℝ => Real.log x * x) (𝓝[>] (0:ℝ)) (𝓝 0) := by
  have h := tendsto_log_mul_rpow_nhdsGT_zero (r := (1:ℝ)) one_pos
  refine h.congr' ?_
  filter_upwards [self_mem_nhdsWithin] with x hx
  rw [Real.rpow_one]

/-- If `f` is differentiable at `a` with `f a = 0`, the difference quotient
`f x / (x - a)` converges to the derivative. -/
theorem slope_tendsto_of_hasDerivAt_eq_zero (f : ℝ → ℝ) (a d : ℝ)
    (h : HasDerivAt f d a) (h0 : f a = 0) :
    Tendsto (fun x : ℝ => f x / (x - a)) (𝓝[≠] a) (𝓝 d) := by
  refine (hasDerivAt_iff_tendsto_slope.mp h).congr ?_
  intro x
  rw [slope_def_field, h0]
  ring

/-- If `f` is differentiable at `0` with `f 0 = 0`, the difference quotient
`f x / x` converges to the derivative. -/
theorem slope_tendsto_of_hasDerivAt_zero (f : ℝ → ℝ) (d : ℝ)
    (h : HasDerivAt f d 0) (h0 : f 0 = 0) :
    Tendsto (fun x : ℝ => f x / x) (𝓝[≠] (0:ℝ)) (𝓝 d) := by
  have := slope_tendsto_of_hasDerivAt_eq_zero f 0 d h h0
  simpa using this

/-- `V x * x → 0` as `x → 0⁺`. -/
theorem quadAltV_mul_self_tendsto :
    Tendsto (fun x : ℝ => quadAltV x * x) (𝓝[>] (0:ℝ)) (𝓝 0) := by
  have t1 : Tendsto (fun x : ℝ => Real.log x ^ 2 / 2 * x) (𝓝[>] (0:ℝ)) (𝓝 0) := by
    have := logSq_mul_self_tendsto.div_const 2
    simpa [mul_comm, mul_div_assoc, mul_assoc] using this
  have t2 : Tendsto (fun x : ℝ => Real.log x * Real.log (1+x) * x) (𝓝[>] (0:ℝ)) (𝓝 0) := by
    have hlog1 : Tendsto (fun x : ℝ => Real.log (1+x)) (𝓝[>] (0:ℝ)) (𝓝 0) := by
      have hc : ContinuousAt (fun x : ℝ => Real.log (1+x)) 0 := by
        apply ContinuousAt.log (by fun_prop); norm_num
      simpa using (hc.tendsto.mono_left nhdsWithin_le_nhds)
    have := log_mul_self_tendsto.mul hlog1
    simpa [mul_zero, mul_comm, mul_assoc, mul_left_comm] using this
  have t3 : Tendsto (fun x : ℝ => dilog (-x) * x) (𝓝[>] (0:ℝ)) (𝓝 0) := by
    have hd : Tendsto (fun x : ℝ => dilog (-x)) (𝓝[>] (0:ℝ)) (𝓝 0) := by
      have hca : ContinuousAt dilog 0 :=
        (RamanujanChallenge.P26.dilog_hasDerivAt_zero26).continuousAt
      have hneg : ContinuousAt (fun x : ℝ => -x) 0 := by fun_prop
      have hca' : ContinuousAt dilog (-(0:ℝ)) := by rwa [neg_zero]
      have hcomp := hca'.comp hneg
      have := hcomp.tendsto.mono_left (nhdsWithin_le_nhds (a := (0:ℝ)) (s := Set.Ioi 0))
      simpa [dilog_zero] using this
    have := hd.mul (tendsto_id.mono_left nhdsWithin_le_nhds)
    simpa using this
  have t4 : Tendsto (fun x : ℝ => Real.pi^2/12 * x) (𝓝[>] (0:ℝ)) (𝓝 0) := by
    have h0 : Tendsto (fun x : ℝ => x) (𝓝[>] (0:ℝ)) (𝓝 0) :=
      tendsto_id.mono_left nhdsWithin_le_nhds
    simpa using h0.const_mul (Real.pi^2/12)
  have hcomb := ((t1.sub t2).sub t3).sub t4
  simp only [sub_zero] at hcomb
  refine hcomb.congr ?_
  intro x
  unfold quadAltV
  ring


/-- `J(-y)` is differentiable at `0` with derivative `0`, and vanishes there. -/
theorem quadAltJneg_hasDerivAt_zero :
    HasDerivAt (fun y : ℝ => quadAltJclosed (-y)) 0 0 := by
  have hF0 : quadAltFclosed 0 = 0 := by
    unfold quadAltFclosed; norm_num
  have hM : HasDerivAt quadAltMclosed
      (quadAltFclosed 0 - 2 * quadAltFclosed (-0)) 0 :=
    quadAltMclosed_hasDerivAt' (by norm_num)
  rw [neg_zero, hF0] at hM
  norm_num at hM
  have hMn : HasDerivAt (fun y : ℝ => quadAltMclosed (-y)) 0 0 := by
    have hM' : HasDerivAt quadAltMclosed 0 ((fun r : ℝ => -r) 0) := by simpa using hM
    have := HasDerivAt.comp 0 hM' (hasDerivAt_neg (0:ℝ))
    simpa using this
  have hsq : HasDerivAt (fun y : ℝ => y ^ 2) 0 0 := by
    simpa using (hasDerivAt_pow 2 (0:ℝ))
  have hdl0 : HasDerivAt dilog 1 0 := RamanujanChallenge.P26.dilog_hasDerivAt_zero26
  have hdsq : HasDerivAt (fun y : ℝ => dilog (y ^ 2)) 0 0 := by
    have hz : HasDerivAt dilog 1 ((fun y : ℝ => y ^ 2) 0) := by
      show HasDerivAt dilog 1 ((0:ℝ) ^ 2)
      rw [show ((0:ℝ)) ^ 2 = 0 from by norm_num]
      exact hdl0
    have hcomp := HasDerivAt.comp (h := fun y : ℝ => y ^ 2) (0:ℝ) hz hsq
    simpa using hcomp
  have hsum := hMn.add hdsq
  simpa [quadAltJclosed] using hsum

theorem quadAltJneg_zero : quadAltJclosed (-(0:ℝ)) = 0 := by
  rw [neg_zero]
  unfold quadAltJclosed
  rw [quadAltMclosed_zero]
  simp [dilog_zero]

/-- `F x = -2 V x · J(-x) → 0` as `x → 0⁺`. -/
theorem quadAltF_tendsto_zero_right :
    Tendsto (fun x : ℝ => -2 * quadAltV x * quadAltJclosed (-x)) (𝓝[>] (0:ℝ)) (𝓝 0) := by
  have hq : Tendsto (fun x : ℝ => quadAltJclosed (-x) / x) (𝓝[≠] (0:ℝ)) (𝓝 0) :=
    slope_tendsto_of_hasDerivAt_zero _ 0 quadAltJneg_hasDerivAt_zero quadAltJneg_zero
  have hq' : Tendsto (fun x : ℝ => quadAltJclosed (-x) / x) (𝓝[>] (0:ℝ)) (𝓝 0) :=
    hq.mono_left (nhdsWithin_mono _ (fun x hx => ne_of_gt hx))
  have hprod := (quadAltV_mul_self_tendsto.const_mul (-2)).mul hq'
  simp only [mul_zero, zero_mul] at hprod
  refine hprod.congr' ?_
  filter_upwards [self_mem_nhdsWithin] with x hx
  have hxne : x ≠ 0 := ne_of_gt hx
  field_simp


/-- Landen identity (Q6047 (4.5)): `Li2(−x) = −Li2(x/(1+x)) − log(1+x)²/2`
for `0 < x < 1`. Proved via `g' = 0` on `(0,1)` and `g(0) = 0`. -/
theorem quadAlt_dilog_landen {x : ℝ} (hx0 : 0 < x) (hx1 : x < 1) :
    dilog (-x) = -dilog (x / (1 + x)) - Real.log (1 + x) ^ 2 / 2 := by
  -- 不必自造：Problem26WeightThree.dilog_landen_half26 就是这条, 代 x/(1+x)
  have h1px : (0:ℝ) < 1 + x := by linarith
  have hu0 : 0 < x / (1 + x) := by positivity
  have huhalf : x / (1 + x) ≤ 1 / 2 := by
    rw [div_le_iff₀ h1px]; linarith
  have h1u : 1 - x / (1 + x) = 1 / (1 + x) := by
    field_simp
    ring
  have hkey := RamanujanChallenge.P26.dilog_landen_half26 hu0 huhalf
  rw [h1u] at hkey
  have harg : -(x / (1 + x)) / (1 / (1 + x)) = -x := by field_simp
  rw [harg] at hkey
  have hlog : Real.log (1 / (1 + x)) = -Real.log (1 + x) := by
    rw [one_div, Real.log_inv]
  rw [hlog] at hkey
  rw [hkey]; ring

/-- Bridge `−2V(x) = W0(2x/(1+x))` for `0 < x < 1` (Q6047 (4.6)-(4.7)+(5.2)):
the Landen form `W` becomes `W0` under `t = 2x/(1+x)` (`u = t/2`). -/
theorem quadAlt_neg2V_eq_W0 {x : ℝ} (hx0 : 0 < x) (hx1 : x < 1) :
    -2 * quadAltV x = W0 (2 * x / (1 + x)) := by
  unfold quadAltV W0
  have hd := quadAlt_dilog_landen hx0 hx1
  rw [hd]
  have hlog : Real.log (2 * x / (1 + x) / 2) = Real.log (x / (1 + x)) := by
    congr 1
    ring
  rw [hlog]
  have hxlog : Real.log (x / (1 + x)) = Real.log x - Real.log (1 + x) := by
    rw [Real.log_div (ne_of_gt hx0) (by positivity : (1 + x) ≠ 0)]
  rw [hxlog]
  ring

/-- `W0(1) = 0` (double zero at `t=1`, needed for boundary terms). -/
theorem quadAltW0_one : W0 1 = 0 := by
  unfold W0
  rw [quadAlt_dilog_half]
  have hlog : Real.log (1 / 2) ^ 2 = Real.log 2 ^ 2 := by
    rw [Real.log_div]
    · rw [Real.log_one]
      ring
    · norm_num
    · norm_num
  rw [hlog]
  ring

/-- `W0` is differentiable at `1` with derivative `0` (the hypothesis of
`quadAltW0_hasDerivAt` only ever needed `t < 2`, so `t = 1` is an interior point). -/
theorem quadAltW0_hasDerivAt_one : HasDerivAt W0 0 1 := by
  have h := quadAltW0_hasDerivAt (t := 1) (by norm_num) (by norm_num)
  norm_num at h
  exact h

/-- `W0 t / (t - 1) → 0` as `t → 1`. -/
theorem quadAltW0_slope_tendsto :
    Tendsto (fun t : ℝ => W0 t / (t - 1)) (𝓝[≠] (1:ℝ)) (𝓝 0) :=
  slope_tendsto_of_hasDerivAt_eq_zero W0 1 0 quadAltW0_hasDerivAt_one quadAltW0_one

/-! ## Right-endpoint atoms: transport through `x ↦ 1-x`, and `V`'s double zero at 1

`V 1 = 0` and `V' 1 = log 1/(1·2) = 0`, so `V` vanishes to second order at `1`.
Proving that directly would want `dilog` differentiable at the boundary point
`-1`; instead we go through `-2 V x = W0 (2x/(1+x))`, for which `1` is an
interior point (`W0` only ever needed `t < 2`). -/

theorem tendsto_one_sub_nhdsWithin :
    Tendsto (fun x : ℝ => 1 - x) (𝓝[<] (1:ℝ)) (𝓝[>] (0:ℝ)) := by
  rw [tendsto_nhdsWithin_iff]
  constructor
  · have hcont : Continuous (fun x : ℝ => 1 - x) := by fun_prop
    have hc : Tendsto (fun x : ℝ => 1 - x) (𝓝 (1:ℝ)) (𝓝 0) := by
      have h1 := hcont.tendsto (1:ℝ)
      simpa using h1
    exact hc.mono_left nhdsWithin_le_nhds
  · filter_upwards [self_mem_nhdsWithin] with x hx
    have hx1 : x < 1 := hx
    simp only [Set.mem_Ioi]
    linarith

theorem oneSub_logSq_tendsto :
    Tendsto (fun x : ℝ => Real.log (1-x) ^ 2 * (1-x)) (𝓝[<] (1:ℝ)) (𝓝 0) :=
  logSq_mul_self_tendsto.comp tendsto_one_sub_nhdsWithin

theorem oneSub_log_tendsto :
    Tendsto (fun x : ℝ => Real.log (1-x) * (1-x)) (𝓝[<] (1:ℝ)) (𝓝 0) :=
  log_mul_self_tendsto.comp tendsto_one_sub_nhdsWithin

/-- The Möbius map `x ↦ 2x/(1+x)` sends `𝓝[<] 1` into `𝓝[≠] 1`. -/
theorem tendsto_mobius_nhdsNe_one :
    Tendsto (fun x : ℝ => 2 * x / (1 + x)) (𝓝[<] (1:ℝ)) (𝓝[≠] (1:ℝ)) := by
  rw [tendsto_nhdsWithin_iff]
  constructor
  · have hc : ContinuousAt (fun x : ℝ => 2 * x / (1 + x)) 1 := by
      apply ContinuousAt.div (by fun_prop) (by fun_prop)
      norm_num
    have h1 := hc.tendsto
    norm_num at h1
    exact h1.mono_left nhdsWithin_le_nhds
  · filter_upwards [self_mem_nhdsWithin,
      (eventually_gt_nhds (show (0:ℝ) < 1 by norm_num)).filter_mono
        nhdsWithin_le_nhds] with x hx hxpos
    have hx1 : x < 1 := hx
    have h1x : (0:ℝ) < 1 + x := by linarith
    simp only [Set.mem_compl_iff, Set.mem_singleton_iff]
    intro hcon
    rw [div_eq_one_iff_eq (ne_of_gt h1x)] at hcon
    linarith

/-- `V x / (x - 1) → 0` as `x → 1⁻`: `V` has a double zero at `1`, seen through
`-2 V x = W0 (2x/(1+x))` where `1` is an interior point of `W0`'s good range. -/
theorem quadAltV_slope_tendsto_one :
    Tendsto (fun x : ℝ => quadAltV x / (x - 1)) (𝓝[<] (1:ℝ)) (𝓝 0) := by
  have hcomp := quadAltW0_slope_tendsto.comp tendsto_mobius_nhdsNe_one
  -- W0(t)/(t-1) = -2·V x·(1+x)/(x-1), so multiply by -1/(2(1+x)) → -1/4
  have hg : Tendsto (fun x : ℝ => -1 / (2 * (1 + x))) (𝓝[<] (1:ℝ)) (𝓝 (-(1/4))) := by
    have hc : ContinuousAt (fun x : ℝ => -1 / (2 * (1 + x))) 1 := by
      apply ContinuousAt.div (by fun_prop) (by fun_prop)
      norm_num
    have h1 := hc.tendsto
    norm_num at h1
    exact h1.mono_left nhdsWithin_le_nhds
  have hmul := hcomp.mul hg
  simp only [zero_mul] at hmul
  refine hmul.congr' ?_
  filter_upwards [self_mem_nhdsWithin,
    (eventually_gt_nhds (show (0:ℝ) < 1 by norm_num)).filter_mono
      nhdsWithin_le_nhds] with x hx hxpos
  have hx1 : x < 1 := hx
  have h1x : (0:ℝ) < 1 + x := by linarith
  have hxne : x - 1 ≠ 0 := by linarith
  have hkey : 2 * x / (1 + x) - 1 = (x - 1) / (1 + x) := by field_simp; ring
  rw [Function.comp_apply, hkey, ← quadAlt_neg2V_eq_W0 hxpos hx1]
  field_simp

/-- `(1-x) · J(-x) → 0` as `x → 1⁻`: `J(-x)` blows up only like `log(1-x)^2`,
which the factor `(1-x)` kills. -/
theorem oneSub_mul_quadAltJneg_tendsto :
    Tendsto (fun x : ℝ => (1 - x) * quadAltJclosed (-x)) (𝓝[<] (1:ℝ)) (𝓝 0) := by
  have hone : Tendsto (fun x : ℝ => 1 - x) (𝓝[<] (1:ℝ)) (𝓝 0) :=
    tendsto_nhdsWithin_iff.mp tendsto_one_sub_nhdsWithin |>.1
  -- bounded factors: log(1+x) → log 2
  have hlog1p : Tendsto (fun x : ℝ => Real.log (1 + x)) (𝓝[<] (1:ℝ)) (𝓝 (Real.log 2)) := by
    have hc : ContinuousAt (fun x : ℝ => Real.log (1 + x)) 1 := by
      apply ContinuousAt.log (by fun_prop); norm_num
    have h1 := hc.tendsto; norm_num at h1
    exact h1.mono_left nhdsWithin_le_nhds
  -- dilog((1-x)/2) → dilog 0 = 0
  have hdil0 : Tendsto (fun x : ℝ => dilog ((1 - x)/2)) (𝓝[<] (1:ℝ)) (𝓝 0) := by
    have hin : Tendsto (fun x : ℝ => (1 - x)/2) (𝓝[<] (1:ℝ)) (𝓝 0) := by
      simpa using hone.div_const 2
    have hcd : ContinuousWithinAt dilog (Icc (-1:ℝ) 1) 0 :=
      dilog_continuousOn_unit 0 (by norm_num)
    have hev : ∀ᶠ x in 𝓝[<] (1:ℝ), (1 - x)/2 ∈ Icc (-1:ℝ) 1 := by
      filter_upwards [self_mem_nhdsWithin,
        (eventually_gt_nhds (show (-1:ℝ) < 1 by norm_num)).filter_mono
          nhdsWithin_le_nhds] with x hx hxg
      have hx1 : x < 1 := hx
      constructor <;> [linarith; linarith]
    have := hcd.tendsto.comp (tendsto_nhdsWithin_iff.mpr ⟨hin, hev⟩)
    simpa [dilog_zero] using this
  -- dilog(x²) → dilog 1
  have hdil1 : Tendsto (fun x : ℝ => dilog (x^2)) (𝓝[<] (1:ℝ)) (𝓝 (Real.pi^2/6)) := by
    have hin : Tendsto (fun x : ℝ => x^2) (𝓝[<] (1:ℝ)) (𝓝 1) := by
      have hc : ContinuousAt (fun x : ℝ => x^2) 1 := by fun_prop
      have h1 := hc.tendsto; norm_num at h1
      exact h1.mono_left nhdsWithin_le_nhds
    have hcd : ContinuousWithinAt dilog (Icc (-1:ℝ) 1) 1 :=
      dilog_continuousOn_unit 1 (by norm_num)
    have hev : ∀ᶠ x in 𝓝[<] (1:ℝ), x^2 ∈ Icc (-1:ℝ) 1 := by
      filter_upwards [self_mem_nhdsWithin,
        (eventually_gt_nhds (show (-1:ℝ) < 1 by norm_num)).filter_mono
          nhdsWithin_le_nhds] with x hx hxg
      have hx1 : x < 1 := hx
      constructor
      · nlinarith
      · nlinarith
    have := hcd.tendsto.comp (tendsto_nhdsWithin_iff.mpr ⟨hin, hev⟩)
    simpa [dilog_one] using this
  -- assemble the seven terms
  have T1 : Tendsto (fun x : ℝ => (1-x) * (Real.log (1+x)^2/2)) (𝓝[<] (1:ℝ)) (𝓝 0) := by
    have := hone.mul ((hlog1p.pow 2).div_const 2); simpa using this
  have T2 : Tendsto (fun x : ℝ => (1-x) * Real.log (1-x)^2) (𝓝[<] (1:ℝ)) (𝓝 0) := by
    have := oneSub_logSq_tendsto; simpa [mul_comm] using this
  have T3 : Tendsto (fun x : ℝ => (1-x) * (2 * Real.log (1+x) * Real.log (1-x)))
      (𝓝[<] (1:ℝ)) (𝓝 0) := by
    have h := (hlog1p.const_mul 2).mul oneSub_log_tendsto
    simpa [mul_comm, mul_assoc, mul_left_comm] using h
  have T4 : Tendsto (fun x : ℝ => (1-x) * (2 * Real.log 2 * Real.log (1-x)))
      (𝓝[<] (1:ℝ)) (𝓝 0) := by
    have h := oneSub_log_tendsto.const_mul (2 * Real.log 2)
    simpa [mul_comm, mul_assoc, mul_left_comm] using h
  have T5 : Tendsto (fun x : ℝ => (1-x) * (Real.pi^2/6 - Real.log 2^2))
      (𝓝[<] (1:ℝ)) (𝓝 0) := by
    have := hone.mul_const (Real.pi^2/6 - Real.log 2^2); simpa using this
  have T6 : Tendsto (fun x : ℝ => (1-x) * (-2 * dilog ((1-x)/2))) (𝓝[<] (1:ℝ)) (𝓝 0) := by
    have := hone.mul (hdil0.const_mul (-2)); simpa using this
  have T7 : Tendsto (fun x : ℝ => (1-x) * dilog (x^2)) (𝓝[<] (1:ℝ)) (𝓝 0) := by
    have := hone.mul hdil1; simpa using this
  have hsum := ((((((T1.add T2).add T3).add T4).add T5).add T6).add T7)
  simp only [add_zero] at hsum
  refine hsum.congr ?_
  intro x
  unfold quadAltJclosed quadAltMclosed
  ring

/-- `F x = -2 V x · J(-x) → 0` as `x → 1⁻`. -/
theorem quadAltF_tendsto_zero_left :
    Tendsto (fun x : ℝ => -2 * quadAltV x * quadAltJclosed (-x)) (𝓝[<] (1:ℝ)) (𝓝 0) := by
  have hA := quadAltV_slope_tendsto_one
  have hB := oneSub_mul_quadAltJneg_tendsto
  have hprod := (hA.mul hB).const_mul (2:ℝ)
  simp only [mul_zero] at hprod
  refine hprod.congr' ?_
  filter_upwards [self_mem_nhdsWithin] with x hx
  have hx1 : x < 1 := hx
  have hne : x - 1 ≠ 0 := by linarith
  field_simp
  ring

/-! ## Integrability atoms for the IBP side conditions

Integrability is much weaker than continuity: `log` and `(log x)^2` are both
integrable on `[0,1]`, so the endpoint singularities of `A` and `B` need only be
DOMINATED, not removed. -/

/-- `|log x| ≤ 4 * x ^ (-(1:ℝ)/4)` for `0 < x ≤ 1`. -/
theorem abs_log_le_rpow {x : ℝ} (hx0 : 0 < x) (hx1 : x ≤ 1) :
    |Real.log x| ≤ 4 * x ^ (-(1:ℝ)/4) := by
  have hlognonpos : Real.log x ≤ 0 := Real.log_nonpos hx0.le hx1
  rw [abs_of_nonpos hlognonpos]
  have hinv : (0:ℝ) ≤ x⁻¹ := by positivity
  have h := Real.log_le_rpow_div hinv (show (0:ℝ) < 1/4 by norm_num)
  rw [Real.log_inv] at h
  have hrw : (x⁻¹) ^ ((1:ℝ)/4) = x ^ (-(1:ℝ)/4) := by
    rw [← Real.rpow_neg_one x, ← Real.rpow_mul hx0.le]
    norm_num
  rw [hrw] at h
  calc -Real.log x ≤ x ^ (-(1:ℝ)/4) / (1/4) := h
    _ = 4 * x ^ (-(1:ℝ)/4) := by ring

/-- `(log x)^2` is interval-integrable on `[0,1]`: dominated by `16 x^(-1/2)`. -/
theorem intervalIntegrable_logSq :
    IntervalIntegrable (fun x : ℝ => Real.log x ^ 2) MeasureTheory.volume 0 1 := by
  have hmaj : IntervalIntegrable (fun x : ℝ => 16 * x ^ (-(1:ℝ)/2)) MeasureTheory.volume 0 1 :=
    (intervalIntegral.intervalIntegrable_rpow' (r := -(1:ℝ)/2) (by norm_num)).const_mul 16
  rw [intervalIntegrable_iff_integrableOn_Ioc_of_le (by norm_num : (0:ℝ) ≤ 1)] at hmaj ⊢
  refine hmaj.mono' ?_ ?_
  · exact (Real.measurable_log.pow_const 2).aestronglyMeasurable
  · filter_upwards [MeasureTheory.ae_restrict_mem measurableSet_Ioc] with x hx
    have hx0 : 0 < x := hx.1
    have hx1 : x ≤ 1 := hx.2
    have h := abs_log_le_rpow hx0 hx1
    have habs : (0:ℝ) ≤ |Real.log x| := abs_nonneg _
    have hr : (0:ℝ) < x ^ (-(1:ℝ)/4) := Real.rpow_pos_of_pos hx0 _
    rw [Real.norm_eq_abs, abs_of_nonneg (sq_nonneg _)]
    calc Real.log x ^ 2 = |Real.log x| ^ 2 := by rw [sq_abs]
      _ ≤ (4 * x ^ (-(1:ℝ)/4)) ^ 2 := by nlinarith
      _ = 16 * x ^ (-(1:ℝ)/2) := by
          rw [mul_pow, ← Real.rpow_natCast (x ^ (-(1:ℝ)/4)) 2, ← Real.rpow_mul hx0.le]
          norm_num

/-- Continuity on the OPEN interval plus an integrable majorant gives
`IntervalIntegrable`: endpoint values are irrelevant, which is exactly what an
integrable endpoint singularity needs. -/
theorem intervalIntegrable_of_continuousOn_Ioo_of_le
    {f g : ℝ → ℝ} {a b : ℝ} (hab : a ≤ b)
    (hf : ContinuousOn f (Set.Ioo a b))
    (hg : IntervalIntegrable g MeasureTheory.volume a b)
    (hfg : ∀ x ∈ Set.Ioo a b, ‖f x‖ ≤ g x) :
    IntervalIntegrable f MeasureTheory.volume a b := by
  rw [intervalIntegrable_iff_integrableOn_Ioo_of_le hab] at hg ⊢
  exact hg.mono' (hf.aestronglyMeasurable measurableSet_Ioo)
    (by
      filter_upwards [MeasureTheory.ae_restrict_mem measurableSet_Ioo] with x hx
      exact hfg x hx)

/-- A function continuous on an open interval and with finite one-sided limits
at both endpoints is interval-integrable.  `extendFrom` supplies the continuous
representative on the closed interval; the endpoint changes are null. -/
theorem intervalIntegrable_of_continuousOn_Ioo_of_tendsto
    {f : ℝ → ℝ} {a b la lb : ℝ} (hab : a < b)
    (hf : ContinuousOn f (Set.Ioo a b))
    (ha : Tendsto f (𝓝[>] a) (𝓝 la))
    (hb : Tendsto f (𝓝[<] b) (𝓝 lb)) :
    IntervalIntegrable f MeasureTheory.volume a b := by
  let fext : ℝ → ℝ := extendFrom (Set.Ioo a b) f
  have hfext : ContinuousOn fext (Set.Icc a b) := by
    exact continuousOn_Icc_extendFrom_Ioo hf ha hb
  have hint : IntervalIntegrable fext MeasureTheory.volume a b :=
    hfext.intervalIntegrable_of_Icc hab.le
  apply hint.congr_ae
  filter_upwards [MeasureTheory.ae_restrict_mem measurableSet_uIoc,
    MeasureTheory.Measure.ae_ne
      (MeasureTheory.volume.restrict (Set.uIoc a b)) b] with x hx hxb
  have hx' : x ∈ Set.Ioo a b := by
    rw [Set.uIoc_of_le hab.le] at hx
    exact ⟨hx.1, lt_of_le_of_ne hx.2 hxb⟩
  exact extendFrom_extends hf x hx'

/-- Multiplication by a function with finite endpoint limits preserves interval
integrability, even when the displayed multiplier has arbitrary endpoint
values. -/
theorem IntervalIntegrable.mul_of_continuousOn_Ioo_of_tendsto
    {f g : ℝ → ℝ} {a b ga gb : ℝ} (hab : a < b)
    (hf : IntervalIntegrable f MeasureTheory.volume a b)
    (hg : ContinuousOn g (Set.Ioo a b))
    (hga : Tendsto g (𝓝[>] a) (𝓝 ga))
    (hgb : Tendsto g (𝓝[<] b) (𝓝 gb)) :
    IntervalIntegrable (fun x => f x * g x) MeasureTheory.volume a b := by
  let gext : ℝ → ℝ := extendFrom (Set.Ioo a b) g
  have hgext : ContinuousOn gext (Set.Icc a b) := by
    exact continuousOn_Icc_extendFrom_Ioo hg hga hgb
  have hint : IntervalIntegrable (fun x => f x * gext x) MeasureTheory.volume a b :=
    hf.mul_continuousOn (by simpa [Set.uIcc_of_le hab.le] using hgext)
  apply hint.congr_ae
  filter_upwards [MeasureTheory.ae_restrict_mem measurableSet_uIoc,
    MeasureTheory.Measure.ae_ne
      (MeasureTheory.volume.restrict (Set.uIoc a b)) b] with x hx hxb
  have hx' : x ∈ Set.Ioo a b := by
    rw [Set.uIoc_of_le hab.le] at hx
    exact ⟨hx.1, lt_of_le_of_ne hx.2 hxb⟩
  change f x * extendFrom (Set.Ioo a b) g x = f x * g x
  rw [extendFrom_extends hg x hx']

/-- The double zero of `W0` at one, in the quantitative form needed to
dominate the logarithmic kernels. -/
theorem quadAltW0_quadratic_tendsto :
    Tendsto (fun t : ℝ => W0 t / (t - 1) ^ 2) (𝓝[≠] (1 : ℝ)) (𝓝 (-2)) := by
  let dW : ℝ → ℝ := fun t => -2 * Real.log (t / (2 - t)) / t
  have hdW_one : HasDerivAt dW (-4) 1 := by
    have hden : HasDerivAt (fun t : ℝ => 2 - t) (-1) 1 := by
      simpa using (hasDerivAt_const (1 : ℝ) 2).sub (hasDerivAt_id (1 : ℝ))
    have hmob : HasDerivAt (fun t : ℝ => t / (2 - t)) 2 1 := by
      convert (hasDerivAt_id (1 : ℝ)).div hden (by norm_num) using 1 <;> norm_num
    have hlog : HasDerivAt (fun t : ℝ => Real.log (t / (2 - t))) 2 1 := by
      have hlogAt : HasDerivAt Real.log 1 ((fun t : ℝ => t / (2 - t)) 1) := by
        convert Real.hasDerivAt_log (by norm_num : (1 : ℝ) ≠ 0) using 1 <;>
          norm_num [div_eq_mul_inv]
      have h := HasDerivAt.comp (h := fun t : ℝ => t / (2 - t)) 1 hlogAt hmob
      convert h using 1 <;> norm_num
    unfold dW
    convert (hlog.const_mul (-2)).div (hasDerivAt_id (1 : ℝ)) (by norm_num) using 1 <;>
      norm_num
  have hdW_zero : dW 1 = 0 := by
    unfold dW
    norm_num
  have hslope : Tendsto (fun t : ℝ => dW t / (t - 1))
      (𝓝[≠] (1 : ℝ)) (𝓝 (-4)) :=
    slope_tendsto_of_hasDerivAt_eq_zero dW 1 (-4) hdW_one hdW_zero
  have hratio : Tendsto (fun t : ℝ => dW t / (2 * (t - 1)))
      (𝓝[≠] (1 : ℝ)) (𝓝 (-2)) := by
    have h := hslope.div_const 2
    norm_num at h
    refine h.congr' ?_
    filter_upwards [self_mem_nhdsWithin] with t ht
    have htne : t - 1 ≠ 0 := sub_ne_zero.mpr ht
    field_simp [htne]
  have hWderiv : ∀ᶠ t in 𝓝[≠] (1 : ℝ), HasDerivAt W0 (dW t) t := by
    filter_upwards [
      (eventually_gt_nhds (show (0 : ℝ) < 1 by norm_num)).filter_mono nhdsWithin_le_nhds,
      (eventually_lt_nhds (show (1 : ℝ) < 2 by norm_num)).filter_mono nhdsWithin_le_nhds]
      with t ht0 ht2
    exact quadAltW0_hasDerivAt ht0 ht2
  have hsqderiv : ∀ᶠ t in 𝓝[≠] (1 : ℝ),
      HasDerivAt (fun y : ℝ => (y - 1) ^ 2) (2 * (t - 1)) t := by
    filter_upwards with t
    simpa [id_eq] using ((hasDerivAt_id t).sub_const 1).pow 2
  have hsqne : ∀ᶠ t in 𝓝[≠] (1 : ℝ), 2 * (t - 1) ≠ 0 := by
    filter_upwards [self_mem_nhdsWithin] with t ht
    exact mul_ne_zero (by norm_num) (sub_ne_zero.mpr ht)
  have hWlim : Tendsto W0 (𝓝[≠] (1 : ℝ)) (𝓝 0) := by
    simpa [quadAltW0_one] using
      quadAltW0_hasDerivAt_one.continuousAt.tendsto.mono_left nhdsWithin_le_nhds
  have hsqlim : Tendsto (fun t : ℝ => (t - 1) ^ 2)
      (𝓝[≠] (1 : ℝ)) (𝓝 0) := by
    have hid : Tendsto (fun t : ℝ => t) (𝓝 (1 : ℝ)) (𝓝 1) := tendsto_id
    have hsub : Tendsto (fun t : ℝ => t - 1) (𝓝 (1 : ℝ)) (𝓝 0) := by
      simpa using (hid.sub_const (1 : ℝ))
    have hpow := hsub.pow 2
    norm_num at hpow
    exact hpow.mono_left
      (nhdsWithin_le_nhds (a := (1 : ℝ)) (s := ({1} : Set ℝ)ᶜ))
  exact HasDerivAt.lhopital_zero_nhdsNE hWderiv hsqderiv hsqne hWlim hsqlim hratio

/-- The logarithmic singularity of `W0` at zero is interval-integrable. -/
theorem quadAltW0_intervalIntegrable :
    IntervalIntegrable W0 MeasureTheory.volume 0 1 := by
  have hlog : IntervalIntegrable Real.log MeasureTheory.volume 0 1 :=
    intervalIntegral.intervalIntegrable_log'
  have hconstLog : IntervalIntegrable (fun _ : ℝ => Real.log 2 ^ 2)
      MeasureTheory.volume 0 1 := intervalIntegrable_const
  have hexpand : IntervalIntegrable
      (fun t : ℝ => Real.log t ^ 2 - 2 * Real.log 2 * Real.log t + Real.log 2 ^ 2)
      MeasureTheory.volume 0 1 :=
    (intervalIntegrable_logSq.sub (hlog.const_mul (2 * Real.log 2))).add hconstLog
  have hlogdiv : IntervalIntegrable (fun t : ℝ => Real.log (t / 2) ^ 2)
      MeasureTheory.volume 0 1 := by
    apply hexpand.congr_ae
    filter_upwards [MeasureTheory.ae_restrict_mem measurableSet_uIoc] with t ht
    have ht0 : 0 < t := by
      simpa [Set.uIoc_of_le (by norm_num : (0 : ℝ) ≤ 1)] using ht.1
    rw [Real.log_div (ne_of_gt ht0) (by norm_num : (2 : ℝ) ≠ 0)]
    ring
  have hdilog_cont : ContinuousOn (fun t : ℝ => dilog (t / 2)) (Set.Icc 0 1) := by
    apply dilog_continuousOn_unit.comp (by fun_prop)
    intro t ht
    constructor <;> linarith [ht.1, ht.2]
  have hdilog : IntervalIntegrable (fun t : ℝ => dilog (t / 2))
      MeasureTheory.volume 0 1 := hdilog_cont.intervalIntegrable_of_Icc (by norm_num)
  have hconstPi : IntervalIntegrable (fun _ : ℝ => Real.pi ^ 2 / 6)
      MeasureTheory.volume 0 1 := intervalIntegrable_const
  unfold W0
  exact (hconstPi.sub (hdilog.const_mul 2)).sub hlogdiv

theorem quadAltH1_div_self_tendsto_zero_right :
    Tendsto (fun t : ℝ => H1 t / t) (𝓝[>] (0 : ℝ)) (𝓝 1) := by
  have hinner : HasDerivAt (fun t : ℝ => 1 - t) (-1) 0 := by
    simpa using (hasDerivAt_const (0 : ℝ) 1).sub (hasDerivAt_id (0 : ℝ))
  have hH1 : HasDerivAt H1 1 0 := by
    unfold H1
    convert (hinner.log (by norm_num)).neg using 1 <;> norm_num
  simpa [H1, smul_eq_mul, div_eq_mul_inv, mul_comm] using
    hH1.tendsto_slope_zero_right

theorem quadAltH2_div_self_tendsto_zero_right :
    Tendsto (fun t : ℝ => H2 t / t) (𝓝[>] (0 : ℝ)) (𝓝 (1 / 2)) := by
  have hinner : HasDerivAt (fun t : ℝ => 1 - t / 2) (-1 / 2) 0 := by
    convert (hasDerivAt_const (0 : ℝ) 1).sub ((hasDerivAt_id (0 : ℝ)).div_const 2) using 1 <;>
      norm_num
  have hH2 : HasDerivAt H2 (1 / 2) 0 := by
    unfold H2
    convert (hinner.log (by norm_num)).neg using 1 <;> norm_num
  simpa [H2, smul_eq_mul, div_eq_mul_inv, mul_comm] using
    hH2.tendsto_slope_zero_right

theorem quadAltH1_mul_oneSub_tendsto_one :
    Tendsto (fun t : ℝ => (1 - t) * H1 t) (𝓝[<] (1 : ℝ)) (𝓝 0) := by
  simpa [H1, mul_comm] using oneSub_log_tendsto.neg

theorem quadAltH2_tendsto_one :
    Tendsto H2 (𝓝[<] (1 : ℝ)) (𝓝 (Real.log 2)) := by
  have hc : ContinuousAt H2 1 := by
    unfold H2
    apply ContinuousAt.neg
    apply ContinuousAt.log (by fun_prop)
    norm_num
  have h := hc.tendsto.mono_left
    (nhdsWithin_le_nhds (a := (1 : ℝ)) (s := Set.Iio 1))
  have hval : H2 1 = Real.log 2 := by
    unfold H2
    norm_num
    rw [Real.log_div]
    · simp
    · norm_num
    · norm_num
  rw [hval] at h
  exact h

theorem quadAltH1_continuousAt {t : ℝ} (ht : t < 1) : ContinuousAt H1 t := by
  unfold H1
  exact (ContinuousAt.log (by fun_prop) (by linarith : 1 - t ≠ 0)).neg

theorem quadAltH2_continuousAt {t : ℝ} (ht : t < 2) : ContinuousAt H2 t := by
  unfold H2
  have hinner : ContinuousAt (fun y : ℝ => 1 - y / 2) t := by fun_prop
  exact (hinner.log (by linarith : 1 - t / 2 ≠ 0)).neg

/-- A finite left endpoint multiplier and a quadratic right endpoint
vanishing factor suffice to multiply the integrable logarithmic kernel `W0`. -/
theorem quadAltW0_mul_intervalIntegrable
    {q : ℝ → ℝ} {q0 : ℝ}
    (hq : ContinuousOn q (Set.Ioo 0 1))
    (hq0 : Tendsto q (𝓝[>] (0 : ℝ)) (𝓝 q0))
    (hq1 : Tendsto (fun t : ℝ => (t - 1) ^ 2 * q t)
      (𝓝[<] (1 : ℝ)) (𝓝 0)) :
    IntervalIntegrable (fun t : ℝ => W0 t * q t) MeasureTheory.volume 0 1 := by
  have hWleft : IntervalIntegrable W0 MeasureTheory.volume 0 (1 / 2) := by
    apply quadAltW0_intervalIntegrable.mono_set
    rw [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1 / 2),
      Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)]
    intro t ht
    exact ⟨ht.1, by linarith [ht.2]⟩
  have hqmid_at : ContinuousAt q (1 / 2) :=
    (hq (1 / 2) (by norm_num)).continuousAt (Ioo_mem_nhds (by norm_num) (by norm_num))
  have hleft : IntervalIntegrable (fun t : ℝ => W0 t * q t)
      MeasureTheory.volume 0 (1 / 2) :=
    IntervalIntegrable.mul_of_continuousOn_Ioo_of_tendsto (by norm_num) hWleft
      (hq.mono (by intro t ht; exact ⟨ht.1, by linarith [ht.2]⟩)) hq0
      (hqmid_at.tendsto.mono_left nhdsWithin_le_nhds)
  have hWratio : Tendsto (fun t : ℝ => W0 t / (t - 1) ^ 2)
      (𝓝[<] (1 : ℝ)) (𝓝 (-2)) :=
    quadAltW0_quadratic_tendsto.mono_left
      (nhdsWithin_mono _ (by intro t ht; exact ne_of_lt ht))
  have hright_lim : Tendsto (fun t : ℝ => W0 t * q t)
      (𝓝[<] (1 : ℝ)) (𝓝 0) := by
    have h := hWratio.mul hq1
    norm_num at h
    refine h.congr' ?_
    filter_upwards [self_mem_nhdsWithin] with t ht
    have htne : t - 1 ≠ 0 := ne_of_lt (sub_neg.mpr ht)
    field_simp [htne]
  have hWq_cont : ContinuousOn (fun t : ℝ => W0 t * q t) (Set.Ioo (1 / 2) 1) := by
    intro t ht
    have ht0 : 0 < t := by linarith [ht.1]
    have ht2 : t < 2 := by linarith [ht.2]
    have hqt : ContinuousAt q t :=
      (hq t ⟨ht0, ht.2⟩).continuousAt (Ioo_mem_nhds ht0 ht.2)
    exact ((quadAltW0_hasDerivAt ht0 ht2).continuousAt.mul hqt).continuousWithinAt
  have hright := intervalIntegrable_of_continuousOn_Ioo_of_tendsto (by norm_num)
    hWq_cont
    (((quadAltW0_hasDerivAt (by norm_num : (0 : ℝ) < 1 / 2) (by norm_num : (1 / 2 : ℝ) < 2)).continuousAt.mul
      hqmid_at).tendsto.mono_left nhdsWithin_le_nhds)
    hright_lim
  exact hleft.trans hright

theorem quadAltI10_kernel_intervalIntegrable : IntervalIntegrable
    (fun t : ℝ => W0 t * H1 t / t) MeasureTheory.volume 0 1 := by
  have hqcont : ContinuousOn (fun t : ℝ => H1 t / t) (Set.Ioo 0 1) := by
    intro t ht
    exact ((quadAltH1_continuousAt ht.2).div continuousAt_id (ne_of_gt ht.1)).continuousWithinAt
  have hone : Tendsto (fun t : ℝ => 1 - t) (𝓝[<] (1 : ℝ)) (𝓝 0) :=
    tendsto_nhdsWithin_iff.mp tendsto_one_sub_nhdsWithin |>.1
  have hid : Tendsto (fun t : ℝ => t) (𝓝[<] (1 : ℝ)) (𝓝 1) :=
    tendsto_id.mono_left nhdsWithin_le_nhds
  have hratio : Tendsto (fun t : ℝ => (1 - t) / t)
      (𝓝[<] (1 : ℝ)) (𝓝 0) := by
    simpa using hone.div hid (by norm_num)
  have hq1 : Tendsto (fun t : ℝ => (t - 1) ^ 2 * (H1 t / t))
      (𝓝[<] (1 : ℝ)) (𝓝 0) := by
    have h := quadAltH1_mul_oneSub_tendsto_one.mul hratio
    norm_num at h
    refine h.congr' ?_
    filter_upwards with t
    ring
  convert quadAltW0_mul_intervalIntegrable hqcont
    quadAltH1_div_self_tendsto_zero_right hq1 using 1 <;> ring

theorem quadAltI11_kernel_intervalIntegrable : IntervalIntegrable
    (fun t : ℝ => W0 t * H1 t / (1 - t)) MeasureTheory.volume 0 1 := by
  have hqcont : ContinuousOn (fun t : ℝ => H1 t / (1 - t)) (Set.Ioo 0 1) := by
    intro t ht
    exact ((quadAltH1_continuousAt ht.2).div (by fun_prop)
      (sub_ne_zero.mpr (ne_of_gt ht.2))).continuousWithinAt
  have hH10 : Tendsto H1 (𝓝[>] (0 : ℝ)) (𝓝 0) := by
    convert (quadAltH1_continuousAt (by norm_num : (0 : ℝ) < 1)).tendsto.mono_left
      (show 𝓝[>] (0 : ℝ) ≤ 𝓝 0 from nhdsWithin_le_nhds) using 1 <;> simp [H1]
  have hden0 : Tendsto (fun t : ℝ => 1 - t) (𝓝[>] (0 : ℝ)) (𝓝 1) := by
    have hid : Tendsto (fun t : ℝ => t) (𝓝[>] (0 : ℝ)) (𝓝 0) :=
      tendsto_id.mono_left nhdsWithin_le_nhds
    simpa using tendsto_const_nhds.sub hid
  have hq0 : Tendsto (fun t : ℝ => H1 t / (1 - t))
      (𝓝[>] (0 : ℝ)) (𝓝 0) := by
    simpa using hH10.div hden0 (by norm_num)
  have hq1 : Tendsto (fun t : ℝ => (t - 1) ^ 2 * (H1 t / (1 - t)))
      (𝓝[<] (1 : ℝ)) (𝓝 0) := by
    refine quadAltH1_mul_oneSub_tendsto_one.congr' ?_
    filter_upwards [self_mem_nhdsWithin] with t ht
    have hne : 1 - t ≠ 0 := sub_ne_zero.mpr (ne_of_gt ht)
    field_simp [hne]
    ring
  convert quadAltW0_mul_intervalIntegrable hqcont hq0 hq1 using 1 <;> ring

theorem quadAltI12_kernel_intervalIntegrable : IntervalIntegrable
    (fun t : ℝ => W0 t * H1 t / (2 - t)) MeasureTheory.volume 0 1 := by
  have hqcont : ContinuousOn (fun t : ℝ => H1 t / (2 - t)) (Set.Ioo 0 1) := by
    intro t ht
    exact ((quadAltH1_continuousAt ht.2).div (by fun_prop)
      (sub_ne_zero.mpr (by linarith [ht.2]))).continuousWithinAt
  have hH10 : Tendsto H1 (𝓝[>] (0 : ℝ)) (𝓝 0) := by
    convert (quadAltH1_continuousAt (by norm_num : (0 : ℝ) < 1)).tendsto.mono_left
      (show 𝓝[>] (0 : ℝ) ≤ 𝓝 0 from nhdsWithin_le_nhds) using 1 <;> simp [H1]
  have hden0 : Tendsto (fun t : ℝ => 2 - t) (𝓝[>] (0 : ℝ)) (𝓝 2) := by
    have hid : Tendsto (fun t : ℝ => t) (𝓝[>] (0 : ℝ)) (𝓝 0) :=
      tendsto_id.mono_left nhdsWithin_le_nhds
    simpa using tendsto_const_nhds.sub hid
  have hq0 : Tendsto (fun t : ℝ => H1 t / (2 - t))
      (𝓝[>] (0 : ℝ)) (𝓝 0) := by
    simpa using hH10.div hden0 (by norm_num)
  have hone : Tendsto (fun t : ℝ => 1 - t) (𝓝[<] (1 : ℝ)) (𝓝 0) :=
    tendsto_nhdsWithin_iff.mp tendsto_one_sub_nhdsWithin |>.1
  have hden1 : Tendsto (fun t : ℝ => 2 - t) (𝓝[<] (1 : ℝ)) (𝓝 1) := by
    have hid : Tendsto (fun t : ℝ => t) (𝓝[<] (1 : ℝ)) (𝓝 1) :=
      tendsto_id.mono_left nhdsWithin_le_nhds
    convert tendsto_const_nhds.sub hid using 1 <;> norm_num
  have hratio : Tendsto (fun t : ℝ => (1 - t) / (2 - t))
      (𝓝[<] (1 : ℝ)) (𝓝 0) := by
    simpa using hone.div hden1 (by norm_num)
  have hq1 : Tendsto (fun t : ℝ => (t - 1) ^ 2 * (H1 t / (2 - t)))
      (𝓝[<] (1 : ℝ)) (𝓝 0) := by
    have h := quadAltH1_mul_oneSub_tendsto_one.mul hratio
    norm_num at h
    refine h.congr' ?_
    filter_upwards with t
    ring
  convert quadAltW0_mul_intervalIntegrable hqcont hq0 hq1 using 1 <;> ring

theorem quadAltI20_kernel_intervalIntegrable : IntervalIntegrable
    (fun t : ℝ => W0 t * H2 t / t) MeasureTheory.volume 0 1 := by
  have hqcont : ContinuousOn (fun t : ℝ => H2 t / t) (Set.Ioo 0 1) := by
    intro t ht
    exact ((quadAltH2_continuousAt (by linarith [ht.2])).div continuousAt_id
      (ne_of_gt ht.1)).continuousWithinAt
  have hone : Tendsto (fun t : ℝ => 1 - t) (𝓝[<] (1 : ℝ)) (𝓝 0) :=
    tendsto_nhdsWithin_iff.mp tendsto_one_sub_nhdsWithin |>.1
  have hid : Tendsto (fun t : ℝ => t) (𝓝[<] (1 : ℝ)) (𝓝 1) :=
    tendsto_id.mono_left nhdsWithin_le_nhds
  have hq1 : Tendsto (fun t : ℝ => (t - 1) ^ 2 * (H2 t / t))
      (𝓝[<] (1 : ℝ)) (𝓝 0) := by
    have h := ((hone.pow 2).mul quadAltH2_tendsto_one).div hid (by norm_num)
    norm_num at h
    refine h.congr' ?_
    filter_upwards with t
    simp only [Pi.div_apply]
    ring
  convert quadAltW0_mul_intervalIntegrable hqcont
    quadAltH2_div_self_tendsto_zero_right hq1 using 1 <;> ring

theorem quadAltI21_kernel_intervalIntegrable : IntervalIntegrable
    (fun t : ℝ => W0 t * H2 t / (1 - t)) MeasureTheory.volume 0 1 := by
  have hqcont : ContinuousOn (fun t : ℝ => H2 t / (1 - t)) (Set.Ioo 0 1) := by
    intro t ht
    exact ((quadAltH2_continuousAt (by linarith [ht.2])).div (by fun_prop)
      (sub_ne_zero.mpr (ne_of_gt ht.2))).continuousWithinAt
  have hH20 : Tendsto H2 (𝓝[>] (0 : ℝ)) (𝓝 0) := by
    convert (quadAltH2_continuousAt (by norm_num : (0 : ℝ) < 2)).tendsto.mono_left
      (show 𝓝[>] (0 : ℝ) ≤ 𝓝 0 from nhdsWithin_le_nhds) using 1 <;> simp [H2]
  have hden0 : Tendsto (fun t : ℝ => 1 - t) (𝓝[>] (0 : ℝ)) (𝓝 1) := by
    have hid : Tendsto (fun t : ℝ => t) (𝓝[>] (0 : ℝ)) (𝓝 0) :=
      tendsto_id.mono_left nhdsWithin_le_nhds
    simpa using tendsto_const_nhds.sub hid
  have hq0 : Tendsto (fun t : ℝ => H2 t / (1 - t))
      (𝓝[>] (0 : ℝ)) (𝓝 0) := by
    simpa using hH20.div hden0 (by norm_num)
  have hone : Tendsto (fun t : ℝ => 1 - t) (𝓝[<] (1 : ℝ)) (𝓝 0) :=
    tendsto_nhdsWithin_iff.mp tendsto_one_sub_nhdsWithin |>.1
  have hq1 : Tendsto (fun t : ℝ => (t - 1) ^ 2 * (H2 t / (1 - t)))
      (𝓝[<] (1 : ℝ)) (𝓝 0) := by
    have h := hone.mul quadAltH2_tendsto_one
    norm_num at h
    refine h.congr' ?_
    filter_upwards [self_mem_nhdsWithin] with t ht
    have hne : 1 - t ≠ 0 := sub_ne_zero.mpr (ne_of_gt ht)
    field_simp [hne]
    ring
  convert quadAltW0_mul_intervalIntegrable hqcont hq0 hq1 using 1 <;> ring

theorem quadAltI22_kernel_intervalIntegrable : IntervalIntegrable
    (fun t : ℝ => W0 t * H2 t / (2 - t)) MeasureTheory.volume 0 1 := by
  have hqcont : ContinuousOn (fun t : ℝ => H2 t / (2 - t)) (Set.Ioo 0 1) := by
    intro t ht
    exact ((quadAltH2_continuousAt (by linarith [ht.2])).div (by fun_prop)
      (sub_ne_zero.mpr (by linarith [ht.2]))).continuousWithinAt
  have hH20 : Tendsto H2 (𝓝[>] (0 : ℝ)) (𝓝 0) := by
    convert (quadAltH2_continuousAt (by norm_num : (0 : ℝ) < 2)).tendsto.mono_left
      (show 𝓝[>] (0 : ℝ) ≤ 𝓝 0 from nhdsWithin_le_nhds) using 1 <;> simp [H2]
  have hden0 : Tendsto (fun t : ℝ => 2 - t) (𝓝[>] (0 : ℝ)) (𝓝 2) := by
    have hid : Tendsto (fun t : ℝ => t) (𝓝[>] (0 : ℝ)) (𝓝 0) :=
      tendsto_id.mono_left nhdsWithin_le_nhds
    simpa using tendsto_const_nhds.sub hid
  have hq0 : Tendsto (fun t : ℝ => H2 t / (2 - t))
      (𝓝[>] (0 : ℝ)) (𝓝 0) := by
    simpa using hH20.div hden0 (by norm_num)
  have hone : Tendsto (fun t : ℝ => 1 - t) (𝓝[<] (1 : ℝ)) (𝓝 0) :=
    tendsto_nhdsWithin_iff.mp tendsto_one_sub_nhdsWithin |>.1
  have hden1 : Tendsto (fun t : ℝ => 2 - t) (𝓝[<] (1 : ℝ)) (𝓝 1) := by
    have hid : Tendsto (fun t : ℝ => t) (𝓝[<] (1 : ℝ)) (𝓝 1) :=
      tendsto_id.mono_left nhdsWithin_le_nhds
    convert tendsto_const_nhds.sub hid using 1 <;> norm_num
  have hq1 : Tendsto (fun t : ℝ => (t - 1) ^ 2 * (H2 t / (2 - t)))
      (𝓝[<] (1 : ℝ)) (𝓝 0) := by
    have h := ((hone.pow 2).mul quadAltH2_tendsto_one).div hden1 (by norm_num)
    norm_num at h
    refine h.congr' ?_
    filter_upwards with t
    simp only [Pi.div_apply]
    ring
  convert quadAltW0_mul_intervalIntegrable hqcont hq0 hq1 using 1 <;> ring

/-! ## Continuity of the pieces on the open interval

Each of `Mclosed`, `Jclosed`, `V`, `Dminus` is continuous away from the points
where a `log` or `dilog` argument leaves its good range. These feed the
integrability side conditions. -/

theorem quadAltMclosed_continuousAt {x : ℝ} (hx : x ∈ Set.Ioo (-1:ℝ) 1) :
    ContinuousAt quadAltMclosed x := by
  have h1 : (0:ℝ) < 1 - x := by linarith [hx.2]
  have h2 : (0:ℝ) < 1 + x := by linarith [hx.1]
  have h3 : (0:ℝ) < (1 + x)/2 := by positivity
  have hlm : ContinuousAt (fun y : ℝ => Real.log (1 - y)) x := by
    apply ContinuousAt.log (by fun_prop); linarith
  have hlp : ContinuousAt (fun y : ℝ => Real.log (1 + y)) x := by
    apply ContinuousAt.log (by fun_prop); linarith
  have hd : ContinuousAt (fun y : ℝ => dilog ((1+y)/2)) x := by
    have habs : |(1+x)/2| < 1 := by
      rw [abs_lt]
      exact ⟨by linarith [hx.1], by linarith [hx.2]⟩
    have hne : (1+x)/2 ≠ 0 := ne_of_gt h3
    have hdil : ContinuousAt dilog ((fun y : ℝ => (1+y)/2) x) :=
      (dilog_hasDerivAt_of_abs_lt_one habs hne).continuousAt
    have hin : ContinuousAt (fun y : ℝ => (1+y)/2) x := by fun_prop
    exact ContinuousAt.comp (f := fun y : ℝ => (1+y)/2) (g := dilog) hdil hin
  unfold quadAltMclosed
  exact ((((((hlm.pow 2).div_const 2).add (hlp.pow 2)).add
      ((hlm.const_mul 2).mul hlp)).add
      ((hlp.const_mul (2 * Real.log 2)))).add continuousAt_const).sub
      continuousAt_const |>.sub (hd.const_mul 2)

theorem quadAltJclosed_continuousAt {x : ℝ} (hx : x ∈ Set.Ioo (-1:ℝ) 1) :
    ContinuousAt quadAltJclosed x := by
  have hM := quadAltMclosed_continuousAt hx
  have hsq : ContinuousAt (fun y : ℝ => dilog (y^2)) x := by
    by_cases hx0 : x = 0
    · subst hx0
      have hdil : ContinuousAt dilog ((fun y : ℝ => y^2) 0) := by
        simpa using RamanujanChallenge.P26.dilog_hasDerivAt_zero26.continuousAt
      exact ContinuousAt.comp (f := fun y : ℝ => y^2) (g := dilog) hdil (by fun_prop)
    · have habs : |x^2| < 1 := by
        rw [abs_of_nonneg (sq_nonneg x)]
        nlinarith [hx.1, hx.2]
      have hne : x^2 ≠ 0 := pow_ne_zero 2 hx0
      have hdil : ContinuousAt dilog ((fun y : ℝ => y^2) x) :=
        (dilog_hasDerivAt_of_abs_lt_one habs hne).continuousAt
      exact ContinuousAt.comp (f := fun y : ℝ => y^2) (g := dilog) hdil (by fun_prop)
  unfold quadAltJclosed
  exact hM.add hsq

theorem quadAltV_continuousAt {x : ℝ} (hx0 : 0 < x) (hx1 : x < 1) :
    ContinuousAt quadAltV x := by
  have hlog : ContinuousAt Real.log x := Real.continuousAt_log (ne_of_gt hx0)
  have hlp : ContinuousAt (fun y : ℝ => Real.log (1 + y)) x := by
    apply ContinuousAt.log (by fun_prop); linarith
  have hdil : ContinuousAt (fun y : ℝ => dilog (-y)) x := by
    have habs : |(-x)| < 1 := by rw [abs_neg, abs_of_pos hx0]; exact hx1
    have hne : (-x) ≠ 0 := by simpa using ne_of_gt hx0
    have hd : ContinuousAt dilog ((fun y : ℝ => -y) x) :=
      (dilog_hasDerivAt_of_abs_lt_one habs hne).continuousAt
    exact ContinuousAt.comp (f := fun y : ℝ => -y) (g := dilog) hd (by fun_prop)
  unfold quadAltV
  exact (((hlog.pow 2).div_const 2).sub (hlog.mul hlp)).sub hdil |>.sub continuousAt_const

theorem quadAltDminus_continuousAt {x : ℝ} (hx0 : 0 < x) (hx1 : x < 1) :
    ContinuousAt quadAltDminus x := by
  have h1m : (0:ℝ) < 1 - x := by linarith
  have h1p : (0:ℝ) < 1 + x := by linarith
  have hlm : ContinuousAt (fun y : ℝ => Real.log (1 - y)) x := by
    apply ContinuousAt.log (by fun_prop); linarith
  have hlp : ContinuousAt (fun y : ℝ => Real.log (1 + y)) x := by
    apply ContinuousAt.log (by fun_prop); linarith
  unfold quadAltDminus
  have t1 : ContinuousAt (fun y : ℝ => -(Real.log (1+y) + 2 * Real.log (1-y))/(1+y)) x := by
    apply ContinuousAt.div (by exact ((hlp.add (hlm.const_mul 2)).neg)) (by fun_prop)
    linarith
  have t2 : ContinuousAt (fun y : ℝ => 2 * (Real.log (1-y) + 2 * Real.log (1+y))/(1-y)) x := by
    apply ContinuousAt.div (by exact ((hlm.add (hlp.const_mul 2)).const_mul 2)) (by fun_prop)
    linarith
  have t3 : ContinuousAt (fun y : ℝ => 2 * (Real.log (1-y) + Real.log (1+y))/y) x := by
    apply ContinuousAt.div (by exact ((hlm.add hlp).const_mul 2)) (by fun_prop)
    linarith
  exact (t1.add t2).add t3

theorem quadAltA_continuousOn : ContinuousOn
    (fun x : ℝ => -2 * (Real.log x / (x * (1 + x))) * quadAltJclosed (-x))
    (Set.Ioo 0 1) := by
  intro x hx
  have hlog : ContinuousAt Real.log x := Real.continuousAt_log (ne_of_gt hx.1)
  have hJ : ContinuousAt (fun y : ℝ => quadAltJclosed (-y)) x := by
    apply ContinuousAt.comp (f := fun y : ℝ => -y) (g := quadAltJclosed)
    · exact quadAltJclosed_continuousAt (by constructor <;> linarith [hx.1, hx.2])
    · fun_prop
  have h1x : 1 + x ≠ 0 := ne_of_gt (by linarith [hx.1] : 0 < 1 + x)
  have hden : x * (1 + x) ≠ 0 := mul_ne_zero (ne_of_gt hx.1) h1x
  exact ((hlog.div (by fun_prop) hden).const_mul (-2)).mul hJ
    |>.continuousWithinAt

theorem quadAltB_continuousOn : ContinuousOn
    (fun x : ℝ => (-2 * quadAltV x) * (-quadAltDminus x)) (Set.Ioo 0 1) := by
  intro x hx
  exact ((quadAltV_continuousAt hx.1 hx.2).const_mul (-2)).mul
    (quadAltDminus_continuousAt hx.1 hx.2).neg |>.continuousWithinAt

/-- `Dminus x → 0` at the left endpoint. -/
theorem quadAltDminus_tendsto_zero_right :
    Tendsto quadAltDminus (𝓝[>] (0 : ℝ)) (𝓝 0) := by
  have hid : Tendsto (fun x : ℝ => x) (𝓝[>] (0 : ℝ)) (𝓝 0) :=
    tendsto_id.mono_left nhdsWithin_le_nhds
  have hplus : Tendsto (fun x : ℝ => Real.log (1 + x))
      (𝓝[>] (0 : ℝ)) (𝓝 0) := by
    have hc : ContinuousAt (fun x : ℝ => Real.log (1 + x)) 0 := by
      apply ContinuousAt.log (by fun_prop)
      norm_num
    simpa using hc.tendsto.mono_left nhdsWithin_le_nhds
  have hminus : Tendsto (fun x : ℝ => Real.log (1 - x))
      (𝓝[>] (0 : ℝ)) (𝓝 0) := by
    have hc : ContinuousAt (fun x : ℝ => Real.log (1 - x)) 0 := by
      apply ContinuousAt.log (by fun_prop)
      norm_num
    simpa using hc.tendsto.mono_left nhdsWithin_le_nhds
  have hplusSlope : Tendsto (fun x : ℝ => Real.log (1 + x) / x)
      (𝓝[>] (0 : ℝ)) (𝓝 1) := by
    have hinner : HasDerivAt (fun x : ℝ => 1 + x) 1 0 := by
      simpa using (hasDerivAt_const (0 : ℝ) 1).add (hasDerivAt_id (0 : ℝ))
    have hlog : HasDerivAt (fun x : ℝ => Real.log (1 + x)) 1 0 := by
      convert hinner.log (by norm_num) using 1 <;> norm_num
    simpa [smul_eq_mul, div_eq_mul_inv, mul_comm] using hlog.tendsto_slope_zero_right
  have hminusSlope : Tendsto (fun x : ℝ => Real.log (1 - x) / x)
      (𝓝[>] (0 : ℝ)) (𝓝 (-1)) := by
    have hinner : HasDerivAt (fun x : ℝ => 1 - x) (-1) 0 := by
      simpa using (hasDerivAt_const (0 : ℝ) 1).sub (hasDerivAt_id (0 : ℝ))
    have hlog : HasDerivAt (fun x : ℝ => Real.log (1 - x)) (-1) 0 := by
      convert hinner.log (by norm_num) using 1 <;> norm_num
    simpa [smul_eq_mul, div_eq_mul_inv, mul_comm] using hlog.tendsto_slope_zero_right
  have h1p : Tendsto (fun x : ℝ => 1 + x) (𝓝[>] (0 : ℝ)) (𝓝 1) := by
    simpa using tendsto_const_nhds.add hid
  have h1m : Tendsto (fun x : ℝ => 1 - x) (𝓝[>] (0 : ℝ)) (𝓝 1) := by
    simpa using tendsto_const_nhds.sub hid
  have t1 : Tendsto
      (fun x : ℝ => -(Real.log (1 + x) + 2 * Real.log (1 - x)) / (1 + x))
      (𝓝[>] (0 : ℝ)) (𝓝 0) := by
    simpa using (hplus.add (hminus.const_mul 2)).neg.div h1p (by norm_num)
  have t2 : Tendsto
      (fun x : ℝ => 2 * (Real.log (1 - x) + 2 * Real.log (1 + x)) / (1 - x))
      (𝓝[>] (0 : ℝ)) (𝓝 0) := by
    simpa using ((hminus.add (hplus.const_mul 2)).const_mul 2).div h1m (by norm_num)
  have t3 : Tendsto
      (fun x : ℝ => 2 * (Real.log (1 - x) + Real.log (1 + x)) / x)
      (𝓝[>] (0 : ℝ)) (𝓝 0) := by
    have h := (hminusSlope.add hplusSlope).const_mul 2
    norm_num at h
    refine h.congr' ?_
    filter_upwards [self_mem_nhdsWithin] with x hx
    have hxne : x ≠ 0 := ne_of_gt hx
    field_simp [hxne]
  have hsum := (t1.add t2).add t3
  norm_num at hsum
  refine hsum.congr' ?_
  filter_upwards with x
  unfold quadAltDminus
  ring

/-- The square of the distance to one absorbs the right-endpoint singularity
of `Dminus`. -/
theorem quadAltDminus_mul_oneSub_sq_tendsto :
    Tendsto (fun x : ℝ => (1 - x) ^ 2 * quadAltDminus x)
      (𝓝[<] (1 : ℝ)) (𝓝 0) := by
  have hone : Tendsto (fun x : ℝ => 1 - x) (𝓝[<] (1 : ℝ)) (𝓝 0) :=
    tendsto_nhdsWithin_iff.mp tendsto_one_sub_nhdsWithin |>.1
  have hid : Tendsto (fun x : ℝ => x) (𝓝[<] (1 : ℝ)) (𝓝 1) :=
    tendsto_id.mono_left nhdsWithin_le_nhds
  have hplus : Tendsto (fun x : ℝ => Real.log (1 + x))
      (𝓝[<] (1 : ℝ)) (𝓝 (Real.log 2)) := by
    have hc : ContinuousAt (fun x : ℝ => Real.log (1 + x)) 1 := by
      apply ContinuousAt.log (by fun_prop)
      norm_num
    convert hc.tendsto.mono_left nhdsWithin_le_nhds using 1 <;> norm_num
  have hdenp : Tendsto (fun x : ℝ => 1 + x) (𝓝[<] (1 : ℝ)) (𝓝 2) := by
    convert tendsto_const_nhds.add hid using 1 <;> norm_num
  have hsqPlus := (hone.pow 2).mul hplus
  have hsqMinus := hone.mul oneSub_log_tendsto
  have t1num := (hsqPlus.add (hsqMinus.const_mul 2)).neg
  have t1 : Tendsto
      (fun x : ℝ => (1 - x) ^ 2 *
        (-(Real.log (1 + x) + 2 * Real.log (1 - x)) / (1 + x)))
      (𝓝[<] (1 : ℝ)) (𝓝 0) := by
    have h := t1num.div hdenp (by norm_num)
    norm_num at h
    refine h.congr' ?_
    filter_upwards [self_mem_nhdsWithin,
      (eventually_gt_nhds (show (0 : ℝ) < 1 by norm_num)).filter_mono nhdsWithin_le_nhds]
      with x hx1 hx0
    have h1x : 1 + x ≠ 0 := by linarith
    simp only [Pi.div_apply]
    field_simp [h1x]
    ring
  have t2raw := (oneSub_log_tendsto.add ((hone.mul hplus).const_mul 2)).const_mul 2
  have t2 : Tendsto
      (fun x : ℝ => (1 - x) ^ 2 *
        (2 * (Real.log (1 - x) + 2 * Real.log (1 + x)) / (1 - x)))
      (𝓝[<] (1 : ℝ)) (𝓝 0) := by
    norm_num at t2raw
    refine t2raw.congr' ?_
    filter_upwards [self_mem_nhdsWithin] with x hx
    have hne : 1 - x ≠ 0 := sub_ne_zero.mpr (ne_of_gt hx)
    field_simp [hne]
  have t3num := ((hone.mul oneSub_log_tendsto).add ((hone.pow 2).mul hplus)).const_mul 2
  have t3 : Tendsto
      (fun x : ℝ => (1 - x) ^ 2 *
        (2 * (Real.log (1 - x) + Real.log (1 + x)) / x))
      (𝓝[<] (1 : ℝ)) (𝓝 0) := by
    have h := t3num.div hid (by norm_num)
    norm_num at h
    refine h.congr' ?_
    filter_upwards [
      (eventually_gt_nhds (show (0 : ℝ) < 1 by norm_num)).filter_mono nhdsWithin_le_nhds]
      with x hx0
    have hxne : x ≠ 0 := ne_of_gt hx0
    simp only [Pi.div_apply]
    field_simp [hxne]
  have hsum := (t1.add t2).add t3
  norm_num at hsum
  simpa [quadAltDminus, mul_add] using hsum

/-- The `A` summand in the integration-by-parts derivative is integrable. -/
theorem quadAltA_intervalIntegrable : IntervalIntegrable
    (fun x : ℝ => -2 * (Real.log x / (x * (1 + x))) * quadAltJclosed (-x))
    MeasureTheory.volume 0 1 := by
  let q : ℝ → ℝ := fun x => (-2 / (1 + x)) * (quadAltJclosed (-x) / x)
  have hqcont : ContinuousOn q (Set.Ioo 0 1) := by
    intro x hx
    have hJ : ContinuousAt (fun y : ℝ => quadAltJclosed (-y)) x := by
      apply ContinuousAt.comp (f := fun y : ℝ => -y) (g := quadAltJclosed)
      · exact quadAltJclosed_continuousAt (by constructor <;> linarith [hx.1, hx.2])
      · fun_prop
    exact ((continuousAt_const.div (by fun_prop)
      (ne_of_gt (by linarith [hx.1] : 0 < 1 + x))).mul
      (hJ.div continuousAt_id (ne_of_gt hx.1))).continuousWithinAt
  have hJq : Tendsto (fun x : ℝ => quadAltJclosed (-x) / x)
      (𝓝[>] (0 : ℝ)) (𝓝 0) := by
    exact (slope_tendsto_of_hasDerivAt_zero _ 0 quadAltJneg_hasDerivAt_zero
      quadAltJneg_zero).mono_left (nhdsWithin_mono _ (fun x hx => ne_of_gt hx))
  have hfac0 : Tendsto (fun x : ℝ => -2 / (1 + x))
      (𝓝[>] (0 : ℝ)) (𝓝 (-2)) := by
    have hc : ContinuousAt (fun x : ℝ => -2 / (1 + x)) 0 := by
      apply ContinuousAt.div (by fun_prop) (by fun_prop)
      norm_num
    simpa using hc.tendsto.mono_left nhdsWithin_le_nhds
  have hq0 : Tendsto q (𝓝[>] (0 : ℝ)) (𝓝 0) := by
    simpa [q] using hfac0.mul hJq
  have hqmid_at : ContinuousAt q (1 / 2) :=
    (hqcont (1 / 2) (by norm_num)).continuousAt (Ioo_mem_nhds (by norm_num) (by norm_num))
  have hqmid : Tendsto q (𝓝[<] (1 / 2 : ℝ)) (𝓝 (q (1 / 2))) :=
    hqmid_at.tendsto.mono_left nhdsWithin_le_nhds
  have hleft0 : IntervalIntegrable (fun x : ℝ => Real.log x * q x)
      MeasureTheory.volume 0 (1 / 2) := by
    exact IntervalIntegrable.mul_of_continuousOn_Ioo_of_tendsto (by norm_num)
      (intervalIntegral.intervalIntegrable_log' :
        IntervalIntegrable Real.log MeasureTheory.volume 0 (1 / 2))
      (hqcont.mono (by intro x hx; exact ⟨hx.1, by linarith [hx.2]⟩)) hq0 hqmid
  have hleft : IntervalIntegrable
      (fun x : ℝ => -2 * (Real.log x / (x * (1 + x))) * quadAltJclosed (-x))
      MeasureTheory.volume 0 (1 / 2) := by
    apply hleft0.congr_ae
    filter_upwards [MeasureTheory.ae_restrict_mem measurableSet_uIoc] with x hx
    have hx0 : 0 < x := by
      simpa [Set.uIoc_of_le (by norm_num : (0 : ℝ) ≤ 1 / 2)] using hx.1
    unfold q
    field_simp [ne_of_gt hx0, ne_of_gt (by linarith : 0 < 1 + x)]
  have hlogSlope : Tendsto (fun x : ℝ => Real.log x / (x - 1))
      (𝓝[<] (1 : ℝ)) (𝓝 1) := by
    have h := slope_tendsto_of_hasDerivAt_eq_zero Real.log 1 1
      (by simpa using Real.hasDerivAt_log (by norm_num : (1 : ℝ) ≠ 0)) Real.log_one
    exact h.mono_left (nhdsWithin_mono _ (fun x hx => ne_of_lt hx))
  have hJone : Tendsto (fun x : ℝ => (x - 1) * quadAltJclosed (-x))
      (𝓝[<] (1 : ℝ)) (𝓝 0) := by
    have h := oneSub_mul_quadAltJneg_tendsto.neg
    norm_num at h
    refine h.congr' ?_
    filter_upwards with x
    ring
  have hfac1 : Tendsto (fun x : ℝ => -2 / (x * (1 + x)))
      (𝓝[<] (1 : ℝ)) (𝓝 (-1)) := by
    have hc : ContinuousAt (fun x : ℝ => -2 / (x * (1 + x))) 1 := by
      apply ContinuousAt.div (by fun_prop) (by fun_prop)
      norm_num
    convert hc.tendsto.mono_left nhdsWithin_le_nhds using 1 <;> norm_num
  have hAone : Tendsto
      (fun x : ℝ => -2 * (Real.log x / (x * (1 + x))) * quadAltJclosed (-x))
      (𝓝[<] (1 : ℝ)) (𝓝 0) := by
    have h := (hlogSlope.mul hJone).mul hfac1
    norm_num at h
    refine h.congr' ?_
    filter_upwards [self_mem_nhdsWithin,
      (eventually_gt_nhds (show (0 : ℝ) < 1 by norm_num)).filter_mono nhdsWithin_le_nhds]
      with x hx1 hx0
    have hxne : x ≠ 0 := ne_of_gt hx0
    have hx1ne : x - 1 ≠ 0 := ne_of_lt (sub_neg.mpr hx1)
    field_simp [hxne, hx1ne]
  have hAmid_at := (quadAltA_continuousOn (1 / 2) (by norm_num)).continuousAt
    (Ioo_mem_nhds (by norm_num) (by norm_num))
  have hright := intervalIntegrable_of_continuousOn_Ioo_of_tendsto (by norm_num)
    (quadAltA_continuousOn.mono (by intro x hx; exact ⟨by linarith [hx.1], hx.2⟩))
    (hAmid_at.tendsto.mono_left nhdsWithin_le_nhds) hAone
  exact hleft.trans hright

/-- `V` itself is integrable on the left half-interval; its only singular term
there is `(log x)^2`. -/
theorem quadAltV_intervalIntegrable_zero_half :
    IntervalIntegrable quadAltV MeasureTheory.volume 0 (1 / 2) := by
  have hsq : IntervalIntegrable (fun x : ℝ => Real.log x ^ 2)
      MeasureTheory.volume 0 (1 / 2) := by
    apply intervalIntegrable_logSq.mono_set
    rw [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1 / 2),
      Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)]
    intro x hx
    exact ⟨hx.1, by linarith [hx.2]⟩
  have hlog : IntervalIntegrable Real.log MeasureTheory.volume 0 (1 / 2) :=
    intervalIntegral.intervalIntegrable_log'
  have hlp : ContinuousOn (fun x : ℝ => Real.log (1 + x)) (Set.uIcc 0 (1 / 2)) := by
    rw [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1 / 2)]
    intro x hx
    apply ContinuousAt.continuousWithinAt
    apply ContinuousAt.log (by fun_prop)
    linarith [hx.1]
  have hprod := hlog.mul_continuousOn hlp
  have hdil_cont : ContinuousOn (fun x : ℝ => dilog (-x)) (Set.Icc 0 (1 / 2)) := by
    apply dilog_continuousOn_unit.comp (by fun_prop)
    intro x hx
    constructor <;> linarith [hx.1, hx.2]
  have hdil : IntervalIntegrable (fun x : ℝ => dilog (-x))
      MeasureTheory.volume 0 (1 / 2) := hdil_cont.intervalIntegrable_of_Icc (by norm_num)
  have hc : IntervalIntegrable (fun _ : ℝ => Real.pi ^ 2 / 12)
      MeasureTheory.volume 0 (1 / 2) := intervalIntegrable_const
  unfold quadAltV
  exact (((hsq.div_const 2).sub hprod).sub hdil).sub hc

/-- The `B` summand in the integration-by-parts derivative is integrable. -/
theorem quadAltB_intervalIntegrable : IntervalIntegrable
    (fun x : ℝ => (-2 * quadAltV x) * (-quadAltDminus x))
    MeasureTheory.volume 0 1 := by
  have hDmid_at : ContinuousAt quadAltDminus (1 / 2) :=
    quadAltDminus_continuousAt (by norm_num) (by norm_num)
  have hVD : IntervalIntegrable (fun x : ℝ => quadAltV x * quadAltDminus x)
      MeasureTheory.volume 0 (1 / 2) := by
    exact IntervalIntegrable.mul_of_continuousOn_Ioo_of_tendsto (by norm_num)
      quadAltV_intervalIntegrable_zero_half
      (fun x hx => (quadAltDminus_continuousAt hx.1 (by linarith [hx.2])).continuousWithinAt)
      quadAltDminus_tendsto_zero_right
      (hDmid_at.tendsto.mono_left nhdsWithin_le_nhds)
  have hleft : IntervalIntegrable
      (fun x : ℝ => (-2 * quadAltV x) * (-quadAltDminus x))
      MeasureTheory.volume 0 (1 / 2) := by
    convert hVD.const_mul 2 using 1
    funext x
    ring
  have hmobD : Tendsto
      (fun x : ℝ => (2 * x / (1 + x) - 1) ^ 2 * (-quadAltDminus x))
      (𝓝[<] (1 : ℝ)) (𝓝 0) := by
    have hfac : Tendsto (fun x : ℝ => 1 / (1 + x) ^ 2)
        (𝓝[<] (1 : ℝ)) (𝓝 (1 / 4)) := by
      have hc : ContinuousAt (fun x : ℝ => 1 / (1 + x) ^ 2) 1 := by
        apply ContinuousAt.div (by fun_prop) (by fun_prop)
        norm_num
      convert hc.tendsto.mono_left nhdsWithin_le_nhds using 1 <;> norm_num
    have h := quadAltDminus_mul_oneSub_sq_tendsto.neg.mul hfac
    norm_num at h
    refine h.congr' ?_
    filter_upwards [
      (eventually_gt_nhds (show (0 : ℝ) < 1 by norm_num)).filter_mono nhdsWithin_le_nhds]
      with x hx0
    have h1x : 1 + x ≠ 0 := by linarith
    field_simp [h1x]
    ring
  have hWmob := quadAltW0_quadratic_tendsto.comp tendsto_mobius_nhdsNe_one
  have hB1raw := hWmob.mul hmobD
  norm_num at hB1raw
  have hB1 : Tendsto
      (fun x : ℝ => (-2 * quadAltV x) * (-quadAltDminus x))
      (𝓝[<] (1 : ℝ)) (𝓝 0) := by
    refine hB1raw.congr' ?_
    filter_upwards [self_mem_nhdsWithin,
      (eventually_gt_nhds (show (0 : ℝ) < 1 by norm_num)).filter_mono nhdsWithin_le_nhds]
      with x hx1 hx0
    have hx1' : x < 1 := hx1
    have h1x : 1 + x ≠ 0 := by linarith
    have hmne : 2 * x / (1 + x) - 1 ≠ 0 := by
      have hdiv : 2 * x / (1 + x) < 1 :=
        (div_lt_one (by linarith : 0 < 1 + x)).mpr (by linarith [hx1'])
      exact ne_of_lt (sub_neg.mpr hdiv)
    have hnum : 2 * x - (1 + x) ≠ 0 := by linarith [hx1']
    rw [quadAlt_neg2V_eq_W0 hx0 hx1]
    field_simp [hmne, h1x, hnum]
  have hBmid_at := (quadAltB_continuousOn (1 / 2) (by norm_num)).continuousAt
    (Ioo_mem_nhds (by norm_num) (by norm_num))
  have hright := intervalIntegrable_of_continuousOn_Ioo_of_tendsto (by norm_num)
    (quadAltB_continuousOn.mono (by intro x hx; exact ⟨by linarith [hx.1], hx.2⟩))
    (hBmid_at.tendsto.mono_left nhdsWithin_le_nhds) hB1
  exact hleft.trans hright

/-- Integration by parts (Q6047 (4.9)): with `F = −2V·J(−x)`, `F(1)=F(0)=0`
(`V(1)=0`, `J(0)=0`), so `∫₀¹ (−log x)/x·Q(−x) = ∫₀¹ (−2V(x))·Dminus(x)`. -/
theorem quadAltCoeffIntegral_eq_neg2V_Dminus :
    (∫ x : ℝ in (0 : ℝ)..1, (-Real.log x) / x * quadAltQclosed (-x)) =
      ∫ x : ℝ in (0 : ℝ)..1, -2 * quadAltV x * quadAltDminus x := by
  -- 1. kernel = −2·V'·J(−x)（逐点）
  have hkernel : (∫ x : ℝ in (0 : ℝ)..1, (-Real.log x) / x * quadAltQclosed (-x))
      = ∫ x : ℝ in (0 : ℝ)..1, -2 * (Real.log x / (x * (1 + x))) * quadAltJclosed (-x) := by
    apply intervalIntegral.integral_congr_ae
    filter_upwards [MeasureTheory.Measure.ae_ne MeasureTheory.volume (1 : ℝ)] with x hxne hx
    have hx0 : 0 < x := by simpa using hx.1
    have hx1 : x < 1 := lt_of_le_of_ne (by simpa using hx.2) hxne
    exact quadAltCoeffKernel_eq_VJ hx0 hx1
  -- 2. 分部积分
  let F : ℝ → ℝ := fun x => -2 * quadAltV x * quadAltJclosed (-x)
  have hFderiv : ∀ x ∈ Set.Ioo (0 : ℝ) 1, HasDerivAt F
      (-2 * (Real.log x / (x * (1 + x))) * quadAltJclosed (-x) +
        (-2 * quadAltV x) * (-quadAltDminus x)) x := by
    intro x hx
    have hx1 : x < 1 := hx.2
    have hV' : HasDerivAt (fun y : ℝ => -2 * quadAltV y)
        (-2 * (Real.log x / (x * (1 + x)))) x := by
      convert (quadAltV_hasDerivAt hx.1 hx1).const_mul (-2) using 1
    have hJ : HasDerivAt (fun y : ℝ => quadAltJclosed (-y)) (-quadAltDminus x) x :=
      quadAltJclosed_neg_hasDerivAt hx.1 hx1
    have hprod := hV'.mul hJ
    unfold F
    convert hprod using 1
  have hFint : IntervalIntegrable
      (fun x : ℝ => -2 * (Real.log x / (x * (1 + x))) * quadAltJclosed (-x) +
        (-2 * quadAltV x) * (-quadAltDminus x)) MeasureTheory.volume 0 1 := by
    exact quadAltA_intervalIntegrable.add quadAltB_intervalIntegrable
  -- endpoint LIMITS replace ContinuousOn on the closed interval: both are 0.
  have hFlim0 : Tendsto F (𝓝[>] (0:ℝ)) (𝓝 0) := quadAltF_tendsto_zero_right
  have hFlim1 : Tendsto F (𝓝[<] (1:ℝ)) (𝓝 0) := quadAltF_tendsto_zero_left
  have hFTC := intervalIntegral.integral_eq_sub_of_hasDerivAt_of_tendsto
    (a := 0) (b := 1) (f := F)
    (f' := fun x : ℝ => -2 * (Real.log x / (x * (1 + x))) * quadAltJclosed (-x) +
      (-2 * quadAltV x) * (-quadAltDminus x))
    (by norm_num) hFderiv hFint hFlim0 hFlim1
  simp only [sub_self] at hFTC
  -- ∫₀¹ F' = 0 → ∫₀¹ (−2V')·J(−x) = −∫₀¹ (−2V)·(−Dminus) = ∫₀¹ (−2V)·Dminus
  have hFTC' : (∫ x : ℝ in (0 : ℝ)..1,
      -2 * (Real.log x / (x * (1 + x))) * quadAltJclosed (-x)) =
      -∫ x : ℝ in (0 : ℝ)..1, (-2 * quadAltV x) * (-quadAltDminus x) := by
    have hsplit : (∫ x in (0 : ℝ)..1, -2 * (Real.log x / (x * (1 + x))) * quadAltJclosed (-x) +
          (-2 * quadAltV x) * (-quadAltDminus x)) = 0 := by
      exact hFTC
    have hA_int : IntervalIntegrable
        (fun x : ℝ => -2 * (Real.log x / (x * (1 + x))) * quadAltJclosed (-x))
        MeasureTheory.volume 0 1 := by
      exact quadAltA_intervalIntegrable
    have hB_int : IntervalIntegrable
        (fun x : ℝ => (-2 * quadAltV x) * (-quadAltDminus x)) MeasureTheory.volume 0 1 := by
      exact quadAltB_intervalIntegrable
    have hsum := intervalIntegral.integral_add (μ := MeasureTheory.volume)
      (a := (0 : ℝ)) (b := 1)
      (f := fun x : ℝ => -2 * (Real.log x / (x * (1 + x))) * quadAltJclosed (-x))
      (g := fun x : ℝ => (-2 * quadAltV x) * (-quadAltDminus x))
      (hf := hA_int) (hg := hB_int)
    rw [hsum] at hsplit
    linarith
  have hval : -∫ x : ℝ in (0 : ℝ)..1, (-2 * quadAltV x) * (-quadAltDminus x)
      = ∫ x : ℝ in (0 : ℝ)..1, -2 * quadAltV x * quadAltDminus x := by
    rw [← intervalIntegral.integral_neg]
    congr 1
    funext x
    ring
  calc
    (∫ x : ℝ in (0 : ℝ)..1, (-Real.log x) / x * quadAltQclosed (-x))
        = ∫ x : ℝ in (0 : ℝ)..1, -2 * (Real.log x / (x * (1 + x))) * quadAltJclosed (-x) := hkernel
    _ = ∫ x : ℝ in (0 : ℝ)..1, -2 * quadAltV x * quadAltDminus x := by
          rw [hFTC', hval]


/-- Möbius substitution (Q6047 (5.1)+(5.4)): `x = t/(2−t)` turns
`∫₀¹ −2V(x)·Dminus(x) dx` into the `W0·(H1/H2)` combination. -/
theorem quadAltMobiusSubst :
    (∫ x : ℝ in (0 : ℝ)..1, -2 * quadAltV x * quadAltDminus x) =
      ∫ t : ℝ in (0 : ℝ)..1,
        W0 t * (H1 t * (-2 / t - 2 / (1 - t) + 2 / (2 - t)) +
                H2 t * (4 / t + 6 / (1 - t) - 5 / (2 - t))) := by
  let f : ℝ → ℝ := fun t => t / (2 - t)
  let f' : ℝ → ℝ := fun t => 2 / (2 - t) ^ 2
  let g : ℝ → ℝ := fun x => -2 * quadAltV x * quadAltDminus x
  have hf : ContinuousOn f (Set.uIcc (0 : ℝ) 1) := by
    unfold f
    apply ContinuousOn.div
    · exact continuousOn_id
    · exact (continuousOn_const.sub continuousOn_id)
    · intro t ht
      have ht1 : t ≤ 1 := by
        simpa [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)] using ht.2
      linarith
  have hff' : ∀ t ∈ Set.Ioo (0 : ℝ) 1, HasDerivAt f (f' t) t := by
    intro t ht
    have hden : HasDerivAt (fun u : ℝ => 2 - u) (-1) t := by
      simpa using (hasDerivAt_const t (2 : ℝ)).sub (hasDerivAt_id t)
    have h2mt : 2 - t ≠ 0 := by linarith [ht.2]
    have hinv : HasDerivAt (fun u : ℝ => (2 - u)⁻¹) (1 / (2 - t) ^ 2) t := by
      have hc := (hasDerivAt_inv (x := 2 - t) h2mt).comp t hden
      convert hc using 1
      field_simp [h2mt]
    have hm : HasDerivAt (fun u : ℝ => u * (2 - u)⁻¹)
        (1 * (2 - t)⁻¹ + t * (1 / (2 - t) ^ 2)) t :=
      (hasDerivAt_id t).mul hinv
    unfold f f'
    have hval : 2 / (2 - t) ^ 2 = 1 * (2 - t)⁻¹ + t * (1 / (2 - t) ^ 2) := by
      field_simp [h2mt]
      ring
    rw [hval]
    simpa [div_eq_mul_inv] using hm
  have hf' : ∀ t ∈ Set.Ioo (0 : ℝ) 1, 0 ≤ f' t := by
    intro t ht
    unfold f'
    positivity
  have hsubst := intervalIntegral.integral_comp_mul_deriv_of_deriv_nonneg
    (a := (0 : ℝ)) (b := 1) (f := f) (f' := f') (g := g) hf (by simpa using hff') (by simpa using hf')
  have hf01 : f 0 = 0 := by unfold f; ring
  have hf11 : f 1 = 1 := by unfold f; ring
  have hsubst' : (∫ x : ℝ in (0 : ℝ)..1, g x) = ∫ t : ℝ in (0 : ℝ)..1, (g ∘ f) t * f' t := by
    rw [hf01, hf11] at hsubst
    exact hsubst.symm
  have hpoint' : ∀ t ∈ Set.Ioc (0 : ℝ) 1, t ≠ 1 →
      (g ∘ f) t * f' t = W0 t * (H1 t * (-2 / t - 2 / (1 - t) + 2 / (2 - t)) +
                H2 t * (4 / t + 6 / (1 - t) - 5 / (2 - t))) := by
    intro t ht htne
    have ht0 : 0 < t := ht.1
    have ht1 : t < 1 := lt_of_le_of_ne ht.2 htne
    have h2mt : 0 < 2 - t := by linarith
    have hft0 : 0 < t / (2 - t) := by positivity
    have hft1 : t / (2 - t) < 1 := by
      rw [div_lt_one h2mt]
      linarith
    have hW := quadAlt_neg2V_eq_W0 hft0 hft1
    have hmob := quadAltMobius_identity ht0 ht1
    unfold g f f'
    calc
      -2 * quadAltV (t / (2 - t)) * quadAltDminus (t / (2 - t)) * (2 / (2 - t) ^ 2)
          = W0 t * (quadAltDminus (t / (2 - t)) * (2 / (2 - t) ^ 2)) := by
              have hnorm : 2 * (t / (2 - t)) / (1 + t / (2 - t)) = t := by
                field_simp [h2mt.ne']
                ring
              rw [quadAlt_neg2V_eq_W0 hft0 hft1, hnorm]
              ring
      _ = W0 t * (H1 t * (-2 / t - 2 / (1 - t) + 2 / (2 - t)) +
                  H2 t * (4 / t + 6 / (1 - t) - 5 / (2 - t))) := by
            simpa [mul_div_assoc] using congrArg (fun z => W0 t * z) hmob
  have hpoint : (∀ᵐ x ∂MeasureTheory.volume, x ∈ Set.uIoc (0 : ℝ) 1 →
      (g ∘ f) x * f' x = W0 x * (H1 x * (-2 / x - 2 / (1 - x) + 2 / (2 - x)) +
                H2 x * (4 / x + 6 / (1 - x) - 5 / (2 - x)))) := by
    filter_upwards [MeasureTheory.Measure.ae_ne MeasureTheory.volume (1 : ℝ)] with t htne ht
    have htI : t ∈ Set.Ioc (0 : ℝ) 1 := by
      simpa [Set.uIoc_of_le (by norm_num : (0 : ℝ) ≤ 1)] using ht
    exact hpoint' t htI htne
  rw [hsubst']
  apply intervalIntegral.integral_congr_ae
  exact hpoint


/-- Six-integral reduction (Q6047 (5.6)): the `W0·(H1/H2)` integrand splits
linearly into the six `I` integrals. -/
theorem quadAltSixIntegralLinear :
    (∫ t : ℝ in (0 : ℝ)..1, W0 t * (H1 t * (-2 / t - 2 / (1 - t) + 2 / (2 - t)) +
        H2 t * (4 / t + 6 / (1 - t) - 5 / (2 - t)))) =
      -2 * I10 - 2 * I11 + 2 * I12 + 4 * I20 + 6 * I21 - 5 * I22 := by
  -- integrand 逐点展开为线性组合（ring）
  have hpt : (∫ t : ℝ in (0 : ℝ)..1, W0 t * (H1 t * (-2 / t - 2 / (1 - t) + 2 / (2 - t)) +
        H2 t * (4 / t + 6 / (1 - t) - 5 / (2 - t)))) =
      ∫ t : ℝ in (0 : ℝ)..1,
        (-2 * (W0 t * H1 t / t) - 2 * (W0 t * H1 t / (1 - t)) +
          2 * (W0 t * H1 t / (2 - t)) + 4 * (W0 t * H2 t / t) +
          6 * (W0 t * H2 t / (1 - t)) - 5 * (W0 t * H2 t / (2 - t))) := by
    apply intervalIntegral.integral_congr
    intro t ht
    unfold W0 H1 H2
    ring
  -- 积分线性拆分（每项可积性 stub）
  have hint : ∀ k : Fin 6, IntervalIntegrable
      (fun t : ℝ => match k with
        | 0 => -2 * (W0 t * H1 t / t)
        | 1 => -2 * (W0 t * H1 t / (1 - t))
        | 2 => 2 * (W0 t * H1 t / (2 - t))
        | 3 => 4 * (W0 t * H2 t / t)
        | 4 => 6 * (W0 t * H2 t / (1 - t))
        | _ => -5 * (W0 t * H2 t / (2 - t)))
      MeasureTheory.volume 0 1 := by
    intro k
    fin_cases k
    · simpa using quadAltI10_kernel_intervalIntegrable.const_mul (-2)
    · simpa using quadAltI11_kernel_intervalIntegrable.const_mul (-2)
    · simpa using quadAltI12_kernel_intervalIntegrable.const_mul 2
    · simpa using quadAltI20_kernel_intervalIntegrable.const_mul 4
    · simpa using quadAltI21_kernel_intervalIntegrable.const_mul 6
    · simpa using quadAltI22_kernel_intervalIntegrable.const_mul (-5)
  rw [hpt]
  -- 线性拆分：∫(Σ) = Σ∫（逐步 integral_add/sub/const_mul）
  have hlin : (∫ t : ℝ in (0 : ℝ)..1,
        (-2 * (W0 t * H1 t / t) - 2 * (W0 t * H1 t / (1 - t)) +
          2 * (W0 t * H1 t / (2 - t)) + 4 * (W0 t * H2 t / t) +
          6 * (W0 t * H2 t / (1 - t)) - 5 * (W0 t * H2 t / (2 - t)))) =
      (∫ t : ℝ in (0 : ℝ)..1, -2 * (W0 t * H1 t / t)) +
        (∫ t : ℝ in (0 : ℝ)..1, -2 * (W0 t * H1 t / (1 - t))) +
        (∫ t : ℝ in (0 : ℝ)..1, 2 * (W0 t * H1 t / (2 - t))) +
        (∫ t : ℝ in (0 : ℝ)..1, 4 * (W0 t * H2 t / t)) +
        (∫ t : ℝ in (0 : ℝ)..1, 6 * (W0 t * H2 t / (1 - t))) -
        (∫ t : ℝ in (0 : ℝ)..1, 5 * (W0 t * H2 t / (2 - t))) := by
    have h0 : IntervalIntegrable (fun t : ℝ => -2 * (W0 t * H1 t / t))
        MeasureTheory.volume 0 1 := by simpa using hint (0 : Fin 6)
    have h1 : IntervalIntegrable (fun t : ℝ => 2 * (W0 t * H1 t / (1 - t)))
        MeasureTheory.volume 0 1 := by
      simpa using quadAltI11_kernel_intervalIntegrable.const_mul 2
    have h2 : IntervalIntegrable (fun t : ℝ => 2 * (W0 t * H1 t / (2 - t)))
        MeasureTheory.volume 0 1 := by simpa using hint (2 : Fin 6)
    have h3 : IntervalIntegrable (fun t : ℝ => 4 * (W0 t * H2 t / t))
        MeasureTheory.volume 0 1 := by simpa using hint (3 : Fin 6)
    have h4 : IntervalIntegrable (fun t : ℝ => 6 * (W0 t * H2 t / (1 - t)))
        MeasureTheory.volume 0 1 := by simpa using hint (4 : Fin 6)
    have h5 : IntervalIntegrable (fun t : ℝ => 5 * (W0 t * H2 t / (2 - t)))
        MeasureTheory.volume 0 1 := by
      simpa using quadAltI22_kernel_intervalIntegrable.const_mul 5
    have h01 := h0.sub h1
    have h012 := h01.add h2
    have h0123 := h012.add h3
    have h01234 := h0123.add h4
    have hneg1 : -(∫ t : ℝ in (0 : ℝ)..1, 2 * (W0 t * H1 t / (1 - t))) =
        ∫ t : ℝ in (0 : ℝ)..1, -2 * (W0 t * H1 t / (1 - t)) := by
      rw [← intervalIntegral.integral_neg]
      congr 1
      funext t
      ring
    rw [intervalIntegral.integral_sub h01234 h5,
      intervalIntegral.integral_add h0123 h4,
      intervalIntegral.integral_add h012 h3,
      intervalIntegral.integral_add h01 h2,
      intervalIntegral.integral_sub h0 h1]
    linear_combination hneg1
  rw [hlin]
  -- 每项 → I（常数提取 + 定义）
  simp [I10, I11, I12, I20, I21, I22,
    intervalIntegral.integral_const_mul, intervalIntegral.integral_mul_const]
  ring


/-- Layer D assembled (Q6047 (4.9)+(5.1)+(5.4)+(5.6)):
`∫₀¹ (−log x)/x·Q(−x) dx = −2I10 − 2I11 + 2I12 + 4I20 + 6I21 − 5I22`. -/
theorem quadAltCoeffIntegral_eq_six :
    (∫ x : ℝ in (0 : ℝ)..1, (-Real.log x) / x * quadAltQclosed (-x)) =
      -2 * I10 - 2 * I11 + 2 * I12 + 4 * I20 + 6 * I21 - 5 * I22 := by
  calc
    (∫ x : ℝ in (0 : ℝ)..1, (-Real.log x) / x * quadAltQclosed (-x))
        = ∫ x : ℝ in (0 : ℝ)..1, -2 * quadAltV x * quadAltDminus x := quadAltCoeffIntegral_eq_neg2V_Dminus
    _ = ∫ t : ℝ in (0 : ℝ)..1,
          W0 t * (H1 t * (-2 / t - 2 / (1 - t) + 2 / (2 - t)) +
                  H2 t * (4 / t + 6 / (1 - t) - 5 / (2 - t))) := quadAltMobiusSubst
    _ = -2 * I10 - 2 * I11 + 2 * I12 + 4 * I20 + 6 * I21 - 5 * I22 := quadAltSixIntegralLinear


/-- `r(t) = log(t/(2−t))` (Q6047 §5). -/
def quadAltR (t : ℝ) : ℝ := Real.log (t / (2 - t))

/-- `W0 t · (H1 t)²/2 → 0` as `t → 0⁺`: `H1 t ~ t` beats `W0 t ~ -(log(t/2))²`. -/
theorem quadAltW0_mul_g11_tendsto_right :
    Tendsto (fun t : ℝ => W0 t * (H1 t ^ 2 / 2)) (𝓝[>] (0:ℝ)) (𝓝 0) := by
  -- H1 t / t → 1, so H1 t ^2/2 = t * (t * (H1 t/t)^2/2)
  have hH1slope : Tendsto (fun t : ℝ => H1 t / t) (𝓝[≠] (0:ℝ)) (𝓝 1) := by
    have hd : HasDerivAt H1 1 0 := by
      have hc : HasDerivAt (fun y : ℝ => 1 - y) (-1) 0 := by
        simpa using (hasDerivAt_const (0:ℝ) (1:ℝ)).sub (hasDerivAt_id (0:ℝ))
      have hlog : HasDerivAt Real.log (1 / (1 - (0:ℝ))) ((fun y : ℝ => 1 - y) 0) := by
        norm_num
        simpa using Real.hasDerivAt_log (by norm_num : (1:ℝ) ≠ 0)
      have hcomp := HasDerivAt.comp (h := fun y : ℝ => 1 - y) (0:ℝ) hlog hc
      unfold H1
      simpa using hcomp.neg
    have h0 : H1 0 = 0 := by unfold H1; simp
    exact slope_tendsto_of_hasDerivAt_zero H1 1 hd h0
  have hH1slope' : Tendsto (fun t : ℝ => H1 t / t) (𝓝[>] (0:ℝ)) (𝓝 1) :=
    hH1slope.mono_left (nhdsWithin_mono _ (fun x hx => ne_of_gt hx))
  -- W0 t · t → 0
  have hW0t : Tendsto (fun t : ℝ => W0 t * t) (𝓝[>] (0:ℝ)) (𝓝 0) := by
    have hid : Tendsto (fun t : ℝ => t) (𝓝[>] (0:ℝ)) (𝓝 0) :=
      tendsto_id.mono_left nhdsWithin_le_nhds
    have hdil : Tendsto (fun t : ℝ => dilog (t/2) * t) (𝓝[>] (0:ℝ)) (𝓝 0) := by
      have hc : ContinuousAt (fun t : ℝ => dilog (t/2)) 0 := by
        have hd0 : ContinuousAt dilog ((fun t : ℝ => t/2) 0) := by
          simpa using RamanujanChallenge.P26.dilog_hasDerivAt_zero26.continuousAt
        exact ContinuousAt.comp (f := fun t : ℝ => t/2) (g := dilog) hd0 (by fun_prop)
      have h1 := hc.tendsto.mono_left (nhdsWithin_le_nhds (a := (0:ℝ)) (s := Set.Ioi 0))
      simpa [dilog_zero] using h1.mul hid
    have hlogsq : Tendsto (fun t : ℝ => Real.log (t/2) ^ 2 * t) (𝓝[>] (0:ℝ)) (𝓝 0) := by
      have hA := logSq_mul_self_tendsto
      have hB := log_mul_self_tendsto
      have := ((hA.sub (hB.const_mul (2 * Real.log 2))).add (hid.const_mul (Real.log 2 ^ 2)))
      simp only [sub_zero, add_zero, mul_zero] at this
      refine this.congr' ?_
      filter_upwards [self_mem_nhdsWithin] with t ht
      have ht0 : (0:ℝ) < t := ht
      rw [Real.log_div (ne_of_gt ht0) (by norm_num)]
      ring
    have := ((hid.const_mul (Real.pi ^ 2 / 6)).sub (hdil.const_mul 2)).sub hlogsq
    simp only [mul_zero, sub_zero] at this
    refine this.congr ?_
    intro t
    unfold W0
    ring
  -- t · (H1 t/t)^2/2 → 0
  have hrest : Tendsto (fun t : ℝ => t * (H1 t / t) ^ 2 / 2) (𝓝[>] (0:ℝ)) (𝓝 0) := by
    have hid : Tendsto (fun t : ℝ => t) (𝓝[>] (0:ℝ)) (𝓝 0) :=
      tendsto_id.mono_left nhdsWithin_le_nhds
    have := ((hid.mul (hH1slope'.pow 2)).div_const 2)
    simpa using this
  have hmul := hW0t.mul hrest
  simp only [mul_zero] at hmul
  refine hmul.congr' ?_
  filter_upwards [self_mem_nhdsWithin] with t ht
  have ht0 : (0:ℝ) < t := ht
  field_simp

/-- `W0 t · (H1 t)²/2 → 0` as `t → 1⁻`: `W0`'s double zero at `1` beats the
`log(1-t)²` growth of `H1²`. -/
theorem quadAltW0_mul_g11_tendsto_left :
    Tendsto (fun t : ℝ => W0 t * (H1 t ^ 2 / 2)) (𝓝[<] (1:ℝ)) (𝓝 0) := by
  -- W0 t = (W0 t/(t-1)) * (t-1), and (t-1)*log(1-t)^2 → 0
  have hslope : Tendsto (fun t : ℝ => W0 t / (t - 1)) (𝓝[<] (1:ℝ)) (𝓝 0) :=
    quadAltW0_slope_tendsto.mono_left
      (nhdsWithin_mono _ (fun x hx => by
        simp only [Set.mem_compl_iff, Set.mem_singleton_iff]
        exact ne_of_lt hx))
  have hlog : Tendsto (fun t : ℝ => (t - 1) * (Real.log (1 - t) ^ 2 / 2))
      (𝓝[<] (1:ℝ)) (𝓝 0) := by
    have h := oneSub_logSq_tendsto
    have := (h.const_mul (-(1:ℝ)/2))
    simp only [mul_zero] at this
    refine this.congr ?_
    intro t
    ring
  have hmul := hslope.mul hlog
  simp only [mul_zero] at hmul
  refine hmul.congr' ?_
  filter_upwards [self_mem_nhdsWithin] with t ht
  have ht1 : t < 1 := ht
  have hne : t - 1 ≠ 0 := by linarith
  unfold H1
  field_simp

/-- The first integration-by-parts summand tends to zero at the left endpoint. -/
theorem quadAltA11_tendsto_right :
    Tendsto (fun t : ℝ => W0 t * (H1 t / (1 - t))) (𝓝[>] (0:ℝ)) (𝓝 0) := by
  -- W0·t → 0 已有(藏在 quadAltW0_mul_g11_tendsto_right 的证明里), 这里独立重建
  have hW0t : Tendsto (fun t : ℝ => W0 t * t) (𝓝[>] (0:ℝ)) (𝓝 0) := by
    have hid : Tendsto (fun t : ℝ => t) (𝓝[>] (0:ℝ)) (𝓝 0) :=
      tendsto_id.mono_left nhdsWithin_le_nhds
    have hdil : Tendsto (fun t : ℝ => dilog (t/2) * t) (𝓝[>] (0:ℝ)) (𝓝 0) := by
      have hd0 : ContinuousAt dilog ((fun t : ℝ => t/2) 0) := by
        simpa using RamanujanChallenge.P26.dilog_hasDerivAt_zero26.continuousAt
      have hc : ContinuousAt (fun t : ℝ => dilog (t/2)) 0 :=
        ContinuousAt.comp (f := fun t : ℝ => t/2) (g := dilog) hd0 (by fun_prop)
      have h1 := hc.tendsto.mono_left (nhdsWithin_le_nhds (a := (0:ℝ)) (s := Set.Ioi 0))
      simpa [dilog_zero] using h1.mul hid
    have hlogsq : Tendsto (fun t : ℝ => Real.log (t/2) ^ 2 * t) (𝓝[>] (0:ℝ)) (𝓝 0) := by
      have hh := ((logSq_mul_self_tendsto.sub
        (log_mul_self_tendsto.const_mul (2 * Real.log 2))).add
        (hid.const_mul (Real.log 2 ^ 2)))
      simp only [sub_zero, add_zero, mul_zero] at hh
      refine hh.congr' ?_
      filter_upwards [self_mem_nhdsWithin] with t ht
      rw [Real.log_div (ne_of_gt ht) (by norm_num)]
      ring
    have hh := ((hid.const_mul (Real.pi ^ 2 / 6)).sub (hdil.const_mul 2)).sub hlogsq
    simp only [mul_zero, sub_zero] at hh
    refine hh.congr ?_
    intro t; unfold W0; ring
  -- H1 t/(t(1-t)) → 1
  have hq : Tendsto (fun t : ℝ => H1 t / (t * (1 - t))) (𝓝[>] (0:ℝ)) (𝓝 1) := by
    have hslope : Tendsto (fun t : ℝ => H1 t / t) (𝓝[≠] (0:ℝ)) (𝓝 1) := by
      have hd : HasDerivAt H1 1 0 := by
        have hc : HasDerivAt (fun y : ℝ => 1 - y) (-1) 0 := by
          simpa using (hasDerivAt_const (0:ℝ) (1:ℝ)).sub (hasDerivAt_id (0:ℝ))
        have hlog : HasDerivAt Real.log (1 / (1 - (0:ℝ))) ((fun y : ℝ => 1 - y) 0) := by
          norm_num
          simpa using Real.hasDerivAt_log (by norm_num : (1:ℝ) ≠ 0)
        have hcomp := HasDerivAt.comp (h := fun y : ℝ => 1 - y) (0:ℝ) hlog hc
        unfold H1; simpa using hcomp.neg
      exact slope_tendsto_of_hasDerivAt_zero H1 1 hd (by unfold H1; simp)
    have hs' : Tendsto (fun t : ℝ => H1 t / t) (𝓝[>] (0:ℝ)) (𝓝 1) :=
      hslope.mono_left (nhdsWithin_mono _ (fun x hx => ne_of_gt hx))
    have hinv : Tendsto (fun t : ℝ => (1 - t)⁻¹) (𝓝[>] (0:ℝ)) (𝓝 1) := by
      have hc : ContinuousAt (fun t : ℝ => (1 - t)⁻¹) 0 := by
        exact ContinuousAt.inv₀ (by fun_prop) (by norm_num)
      have h1 := hc.tendsto; norm_num at h1
      exact h1.mono_left nhdsWithin_le_nhds
    have := hs'.mul hinv
    simp only [mul_one] at this
    refine this.congr' ?_
    filter_upwards [self_mem_nhdsWithin] with t ht
    have h0 : t ≠ 0 := ne_of_gt ht
    field_simp
  have hmul := hW0t.mul hq
  simp only [zero_mul] at hmul
  refine hmul.congr' ?_
  filter_upwards [self_mem_nhdsWithin] with t ht
  have h0 : t ≠ 0 := ne_of_gt ht
  field_simp

/-- The first summand in the `I11` integration-by-parts derivative is
integrable.  Near one it is dominated using only the first-order quotient
`W0 t / (t - 1)`. -/
theorem quadAltA11_intervalIntegrable : IntervalIntegrable
    (fun y : ℝ => W0 y * (H1 y / (1 - y))) MeasureTheory.volume 0 1 := by
  let A : ℝ → ℝ := fun y => W0 y * (H1 y / (1 - y))
  have hAcont : ContinuousOn A (Set.Ioo 0 1) := by
    intro y hy
    exact ((quadAltW0_hasDerivAt hy.1 (by linarith [hy.2])).continuousAt.mul
      ((quadAltH1_continuousAt hy.2).div (by fun_prop)
        (sub_ne_zero.mpr (ne_of_gt hy.2)))).continuousWithinAt
  have hAleft : ContinuousOn A (Set.Ioo 0 (3 / 4 : ℝ)) :=
    hAcont.mono (by intro y hy; exact ⟨hy.1, by linarith [hy.2]⟩)
  have hAmid : Tendsto A (𝓝[<] (3 / 4 : ℝ)) (𝓝 (A (3 / 4))) := by
    have hAt : ContinuousAt A (3 / 4 : ℝ) :=
      (hAcont (3 / 4) (by norm_num)).continuousAt
        (Ioo_mem_nhds (by norm_num) (by norm_num))
    exact hAt.tendsto.mono_left nhdsWithin_le_nhds
  let Aleft : ℝ → ℝ := extendFrom (Set.Ioo 0 (3 / 4 : ℝ)) A
  have hAleft_cont : ContinuousOn Aleft (Set.Icc 0 (3 / 4 : ℝ)) := by
    exact continuousOn_Icc_extendFrom_Ioo hAleft quadAltA11_tendsto_right hAmid
  obtain ⟨C0, hC0⟩ := isCompact_Icc.exists_bound_of_continuousOn hAleft_cont
  let q : ℝ → ℝ := fun y => W0 y / (y - 1)
  have hqcont : ContinuousOn q (Set.Ioo (1 / 2 : ℝ) 1) := by
    intro y hy
    have hW : ContinuousAt W0 y :=
      (quadAltW0_hasDerivAt (by linarith [hy.1]) (by linarith [hy.2])).continuousAt
    have hden : ContinuousAt (fun z : ℝ => z - 1) y :=
      continuousAt_id.sub continuousAt_const
    exact (hW.div hden (sub_ne_zero.mpr (ne_of_lt hy.2))).continuousWithinAt
  have hqmid : Tendsto q (𝓝[>] (1 / 2 : ℝ)) (𝓝 (q (1 / 2))) := by
    have hqt : ContinuousAt q (1 / 2 : ℝ) := by
      dsimp [q]
      have hW : ContinuousAt W0 (1 / 2 : ℝ) :=
        (quadAltW0_hasDerivAt (by norm_num) (by norm_num)).continuousAt
      exact hW.div (by fun_prop) (by norm_num)
    exact hqt.tendsto.mono_left nhdsWithin_le_nhds
  have hqone : Tendsto q (𝓝[<] (1 : ℝ)) (𝓝 0) := by
    exact quadAltW0_slope_tendsto.mono_left
      (nhdsWithin_mono _ (fun y hy => ne_of_lt hy))
  let qright : ℝ → ℝ := extendFrom (Set.Ioo (1 / 2 : ℝ) 1) q
  have hqright_cont : ContinuousOn qright (Set.Icc (1 / 2 : ℝ) 1) := by
    exact continuousOn_Icc_extendFrom_Ioo hqcont hqmid hqone
  obtain ⟨C1, hC1⟩ := isCompact_Icc.exists_bound_of_continuousOn hqright_cont
  have hH1 : IntervalIntegrable H1 MeasureTheory.volume 0 1 := by
    have hlog : IntervalIntegrable (fun y : ℝ => Real.log (1 - y))
        MeasureTheory.volume 0 1 := by
      simpa using ((intervalIntegral.intervalIntegrable_log' :
        IntervalIntegrable Real.log MeasureTheory.volume 0 1).comp_sub_left 1).symm
    simpa [H1] using hlog.neg
  let C : ℝ := max |C0| |C1|
  have hC0C : |C0| ≤ C := le_max_left _ _
  have hC1C : |C1| ≤ C := le_max_right _ _
  have hCnonneg : 0 ≤ C := (abs_nonneg C0).trans hC0C
  let g : ℝ → ℝ := fun y => C * (1 + ‖H1 y‖)
  have hg : IntervalIntegrable g MeasureTheory.volume 0 1 := by
    have hOne : IntervalIntegrable (fun _ : ℝ => (1 : ℝ))
        MeasureTheory.volume 0 1 := intervalIntegrable_const
    exact (hOne.add hH1.norm).const_mul C
  apply intervalIntegrable_of_continuousOn_Ioo_of_le (by norm_num) hAcont hg
  intro y hy
  by_cases hyl : y < 3 / 4
  · have hyleft : y ∈ Set.Ioo (0 : ℝ) (3 / 4) := ⟨hy.1, hyl⟩
    have heq : Aleft y = A y := extendFrom_extends hAleft y hyleft
    calc
      ‖A y‖ = ‖Aleft y‖ := by rw [heq]
      _ ≤ C0 := hC0 y ⟨hy.1.le, hyl.le⟩
      _ ≤ |C0| := le_abs_self C0
      _ ≤ C := hC0C
      _ ≤ C * (1 + ‖H1 y‖) := by
        have hfac : 1 ≤ 1 + ‖H1 y‖ := le_add_of_nonneg_right (norm_nonneg _)
        simpa using (mul_le_mul_of_nonneg_left hfac hCnonneg)
  · have hyr : (1 / 2 : ℝ) < y := by linarith
    have hyright : y ∈ Set.Ioo (1 / 2 : ℝ) 1 := ⟨hyr, hy.2⟩
    have heq : qright y = q y := extendFrom_extends hqcont y hyright
    have hqbound : ‖q y‖ ≤ |C1| := by
      calc
        ‖q y‖ = ‖qright y‖ := by rw [heq]
        _ ≤ C1 := hC1 y ⟨hyr.le, hy.2.le⟩
        _ ≤ |C1| := le_abs_self C1
    have hnorm : ‖A y‖ = ‖q y‖ * ‖H1 y‖ := by
      have hy1 : y - 1 ≠ 0 := sub_ne_zero.mpr (ne_of_lt hy.2)
      have h1y : 1 - y ≠ 0 := sub_ne_zero.mpr (ne_of_gt hy.2)
      have hAq : A y = -(q y * H1 y) := by
        dsimp [A, q]
        field_simp [hy1, h1y]
        ring
      rw [hAq, norm_neg, norm_mul]
    calc
      ‖A y‖ = ‖q y‖ * ‖H1 y‖ := hnorm
      _ ≤ |C1| * ‖H1 y‖ := mul_le_mul_of_nonneg_right hqbound (norm_nonneg _)
      _ ≤ C * ‖H1 y‖ := mul_le_mul_of_nonneg_right hC1C (norm_nonneg _)
      _ ≤ C * (1 + ‖H1 y‖) := by
        exact mul_le_mul_of_nonneg_left (le_add_of_nonneg_left (by norm_num)) hCnonneg

/-- The second summand in the `I11` integration-by-parts derivative is
integrable.  Its endpoint limits use `t * log t → 0` at zero and the simple
zero of `quadAltR` at one; a single log-square majorant then works globally. -/
theorem quadAltB11_intervalIntegrable : IntervalIntegrable
    (fun y : ℝ => -2 * (quadAltR y / y) * (H1 y ^ 2 / 2))
    MeasureTheory.volume 0 1 := by
  let g11 : ℝ → ℝ := fun t => H1 t ^ 2 / 2
  let B : ℝ → ℝ := fun y => -2 * (quadAltR y / y) * g11 y
  have hBcont : ContinuousOn B (Set.Ioo 0 1) := by
    intro y hy
    have hratio : ContinuousAt (fun z : ℝ => z / (2 - z)) y :=
      continuousAt_id.div (by fun_prop) (sub_ne_zero.mpr (by linarith [hy.2]))
    have hR : ContinuousAt quadAltR y := by
      unfold quadAltR
      exact hratio.log (div_ne_zero (ne_of_gt hy.1)
        (sub_ne_zero.mpr (by linarith [hy.2])))
    have hg : ContinuousAt g11 y := by
      dsimp [g11]
      exact ((quadAltH1_continuousAt hy.2).pow 2).div_const 2
    exact ((continuousAt_const.mul
      (hR.div continuousAt_id (ne_of_gt hy.1))).mul hg).continuousWithinAt
  have hB0 : Tendsto B (𝓝[>] (0 : ℝ)) (𝓝 0) := by
    have hid : Tendsto (fun y : ℝ => y) (𝓝[>] (0 : ℝ)) (𝓝 0) :=
      tendsto_id.mono_left nhdsWithin_le_nhds
    have hlogden : Tendsto (fun y : ℝ => Real.log (2 - y))
        (𝓝[>] (0 : ℝ)) (𝓝 (Real.log 2)) := by
      have hc : ContinuousAt (fun y : ℝ => Real.log (2 - y)) 0 := by
        exact (show ContinuousAt (fun y : ℝ => 2 - y) 0 by fun_prop).log (by norm_num)
      simpa using hc.tendsto.mono_left nhdsWithin_le_nhds
    have hRt : Tendsto (fun y : ℝ => quadAltR y * y)
        (𝓝[>] (0 : ℝ)) (𝓝 0) := by
      have h := log_mul_self_tendsto.sub (hlogden.mul hid)
      simp only [mul_zero, sub_zero] at h
      refine h.congr' ?_
      filter_upwards [self_mem_nhdsWithin,
        (eventually_lt_nhds (show (0 : ℝ) < 2 by norm_num)).filter_mono
          nhdsWithin_le_nhds] with y hy hy2
      unfold quadAltR
      rw [Real.log_div (ne_of_gt hy) (by linarith : 2 - y ≠ 0)]
      ring
    have h := (hRt.mul (quadAltH1_div_self_tendsto_zero_right.pow 2)).neg
    simp only [zero_mul, neg_zero] at h
    refine h.congr' ?_
    filter_upwards [self_mem_nhdsWithin] with y hy
    have hy0 : y ≠ 0 := ne_of_gt hy
    dsimp [B, g11]
    field_simp [hy0]
  have hRslope : Tendsto (fun y : ℝ => quadAltR y / (y - 1))
      (𝓝[<] (1 : ℝ)) (𝓝 2) := by
    have hden : HasDerivAt (fun y : ℝ => 2 - y) (-1) 1 := by
      simpa using (hasDerivAt_const (1 : ℝ) 2).sub (hasDerivAt_id (1 : ℝ))
    have hmob : HasDerivAt (fun y : ℝ => y / (2 - y)) 2 1 := by
      convert (hasDerivAt_id (1 : ℝ)).div hden (by norm_num) using 1 <;> norm_num
    have hlog : HasDerivAt quadAltR 2 1 := by
      unfold quadAltR
      have hlogAt : HasDerivAt Real.log 1 ((fun y : ℝ => y / (2 - y)) 1) := by
        convert Real.hasDerivAt_log (by norm_num : (1 : ℝ) ≠ 0) using 1 <;>
          norm_num [div_eq_mul_inv]
      have h := HasDerivAt.comp (h := fun y : ℝ => y / (2 - y)) 1 hlogAt hmob
      convert h using 1 <;> norm_num
    have hR1 : quadAltR 1 = 0 := by unfold quadAltR; norm_num
    exact (slope_tendsto_of_hasDerivAt_eq_zero quadAltR 1 2 hlog hR1).mono_left
      (nhdsWithin_mono _ (fun y hy => ne_of_lt hy))
  have hHsq : Tendsto (fun y : ℝ => (1 - y) * H1 y ^ 2)
      (𝓝[<] (1 : ℝ)) (𝓝 0) := by
    simpa [H1, mul_comm] using oneSub_logSq_tendsto
  have hinv : Tendsto (fun y : ℝ => y⁻¹) (𝓝[<] (1 : ℝ)) (𝓝 1) := by
    have hc : ContinuousAt (fun y : ℝ => y⁻¹) 1 := continuousAt_id.inv₀ (by norm_num)
    simpa using hc.tendsto.mono_left nhdsWithin_le_nhds
  have hB1 : Tendsto B (𝓝[<] (1 : ℝ)) (𝓝 0) := by
    have h := (hRslope.mul hHsq).mul hinv
    simp only [mul_zero, zero_mul] at h
    refine h.congr' ?_
    filter_upwards [self_mem_nhdsWithin,
      (eventually_gt_nhds (show (0 : ℝ) < 1 by norm_num)).filter_mono
        nhdsWithin_le_nhds] with y hy1 hy0
    have hyne : y ≠ 0 := ne_of_gt hy0
    have hy1ne : y - 1 ≠ 0 := ne_of_lt (sub_neg.mpr hy1)
    dsimp [B, g11]
    field_simp [hyne, hy1ne]
    ring
  let Bext : ℝ → ℝ := extendFrom (Set.Ioo 0 1) B
  have hBext_cont : ContinuousOn Bext (Set.Icc 0 1) := by
    exact continuousOn_Icc_extendFrom_Ioo hBcont hB0 hB1
  obtain ⟨C, hC⟩ := isCompact_Icc.exists_bound_of_continuousOn hBext_cont
  have hOneSubLogSq : IntervalIntegrable (fun y : ℝ => Real.log (1 - y) ^ 2)
      MeasureTheory.volume 0 1 := by
    simpa using (intervalIntegrable_logSq.comp_sub_left 1).symm
  let g : ℝ → ℝ := fun y => |C| *
    (1 + Real.log y ^ 2 + Real.log (1 - y) ^ 2)
  have hg : IntervalIntegrable g MeasureTheory.volume 0 1 := by
    have hOne : IntervalIntegrable (fun _ : ℝ => (1 : ℝ))
        MeasureTheory.volume 0 1 := intervalIntegrable_const
    exact ((hOne.add intervalIntegrable_logSq).add
      hOneSubLogSq).const_mul |C|
  have hBint : IntervalIntegrable B MeasureTheory.volume 0 1 := by
    apply intervalIntegrable_of_continuousOn_Ioo_of_le (by norm_num) hBcont hg
    intro y hy
    have heq : Bext y = B y := extendFrom_extends hBcont y hy
    calc
      ‖B y‖ = ‖Bext y‖ := by rw [heq]
      _ ≤ C := hC y ⟨hy.1.le, hy.2.le⟩
      _ ≤ |C| := le_abs_self C
      _ ≤ |C| * (1 + Real.log y ^ 2 + Real.log (1 - y) ^ 2) := by
        have hfac : 1 ≤ 1 + Real.log y ^ 2 + Real.log (1 - y) ^ 2 := by
          nlinarith [sq_nonneg (Real.log y), sq_nonneg (Real.log (1 - y))]
        simpa using mul_le_mul_of_nonneg_left hfac (abs_nonneg C)
  simpa [B, g11] using hBint

/-- `I11` via the derivative certificate (Q6047 (6.6) with `g11 = H1²/2`):
`I11 = ∫₀¹ r(t)·H1(t)²/t dt`. -/
theorem quadAltI11_eq_integral :
    I11 = ∫ t : ℝ in (0 : ℝ)..1, quadAltR t * H1 t ^ 2 / t := by
  unfold I11
  let g11 : ℝ → ℝ := fun t => H1 t ^ 2 / 2
  have hg11' : ∀ t ∈ Set.Ioo (0 : ℝ) 1, HasDerivAt g11 (H1 t / (1 - t)) t := by
    intro t ht
    have h1mt : 1 - t ≠ 0 := by linarith [ht.2]
    have hd : HasDerivAt (fun y : ℝ => Real.log (1 - y)) (-1 / (1 - t)) t := by
      have hc : HasDerivAt (fun y : ℝ => 1 - y) (-1) t := by
        simpa using (hasDerivAt_const t (1 : ℝ)).sub (hasDerivAt_id t)
      have hlog := Real.hasDerivAt_log h1mt
      have hcomp := HasDerivAt.comp (h := fun y : ℝ => 1 - y) t hlog hc
      convert hcomp using 1
      field_simp [h1mt]
    have hH1' : HasDerivAt H1 (1 / (1 - t)) t := by
      unfold H1
      convert hd.neg using 1
      field_simp [h1mt]
    have hpow := hH1'.pow 2
    unfold g11
    convert hpow.div_const 2 using 1
    field_simp [h1mt]
    ring
  have hg11_0 : g11 0 = 0 := by
    unfold g11 H1
    simp
  have hg11_1 : g11 1 = 0 := by
    unfold g11 H1
    simp
  have hprod' : ∀ t ∈ Set.Ioo (0 : ℝ) 1,
      HasDerivAt (fun y : ℝ => W0 y * g11 y)
        ((W0 t * (H1 t / (1 - t))) + (-2 * (quadAltR t / t) * g11 t)) t := by
    intro t ht
    have hW0' := quadAltW0_hasDerivAt ht.1 (by linarith [ht.2])
    -- hW0' : HasDerivAt W0 (-2 * log(t/(2-t)) / t) t = -2·r/t
    have hg := hg11' t ht
    have hprod := hW0'.mul hg
    unfold g11
    convert hprod using 1
    unfold quadAltR
    ring
  have hlim0 : Tendsto (fun y : ℝ => W0 y * g11 y) (𝓝[>] (0:ℝ)) (𝓝 0) := by
    unfold g11; exact quadAltW0_mul_g11_tendsto_right
  have hlim1 : Tendsto (fun y : ℝ => W0 y * g11 y) (𝓝[<] (1:ℝ)) (𝓝 0) := by
    unfold g11; exact quadAltW0_mul_g11_tendsto_left
  have hA : IntervalIntegrable (fun y : ℝ => W0 y * (H1 y / (1 - y)))
      MeasureTheory.volume 0 1 := quadAltA11_intervalIntegrable
  have hB : IntervalIntegrable (fun y : ℝ => -2 * (quadAltR y / y) * g11 y)
      MeasureTheory.volume 0 1 := by
    simpa [g11] using quadAltB11_intervalIntegrable
  have hint : IntervalIntegrable
      (fun y : ℝ => (W0 y * (H1 y / (1 - y))) + (-2 * (quadAltR y / y) * g11 y))
      MeasureTheory.volume 0 1 := by
    simpa using hA.add hB
  have hFTC := intervalIntegral.integral_eq_sub_of_hasDerivAt_of_tendsto
    (a := 0) (b := 1)
    (f := fun y : ℝ => W0 y * g11 y)
    (f' := fun y : ℝ => (W0 y * (H1 y / (1 - y))) + (-2 * (quadAltR y / y) * g11 y))
    (by norm_num) hprod' hint hlim0 hlim1
  simp only [sub_self] at hFTC
  -- ∫(W0·H1/(1−t)) + ∫(−2r/t·g11) = 0 → ∫W0·H1/(1−t) = ∫ 2r/t·g11 = ∫ r·H1²/t
  have hsplit : (∫ y : ℝ in (0 : ℝ)..1, (W0 y * (H1 y / (1 - y))) + (-2 * (quadAltR y / y) * g11 y)) = 0 := by
    exact hFTC
  have hsum := intervalIntegral.integral_add (μ := MeasureTheory.volume)
    (a := (0 : ℝ)) (b := 1)
    (f := fun y : ℝ => W0 y * (H1 y / (1 - y)))
    (g := fun y : ℝ => -2 * (quadAltR y / y) * g11 y)
    (hf := hA) (hg := hB)
  rw [hsum] at hsplit
  have hmul : ∫ y : ℝ in (0 : ℝ)..1, -2 * (quadAltR y / y) * g11 y
      = -(∫ y : ℝ in (0 : ℝ)..1, 2 * (quadAltR y / y) * g11 y) := by
    rw [← intervalIntegral.integral_neg]
    congr 1
    funext y
    ring
  rw [hmul] at hsplit
  have hmain : (∫ y : ℝ in (0 : ℝ)..1, W0 y * (H1 y / (1 - y)))
      = ∫ y : ℝ in (0 : ℝ)..1, 2 * (quadAltR y / y) * g11 y := by
    linarith
  have hcongr : (∫ y : ℝ in (0 : ℝ)..1, 2 * (quadAltR y / y) * g11 y)
      = ∫ y : ℝ in (0 : ℝ)..1, quadAltR y * H1 y ^ 2 / y := by
    apply intervalIntegral.integral_congr
    intro y hy
    unfold g11
    ring
  simpa [mul_div_assoc] using (hmain.trans hcongr)


/-! ## Layer E, step 1: the IBP that removes `K`'s termwise integration

`d[log²x · log²(1+x)/2] = log x·log²(1+x)/x + log²x·log(1+x)/(1+x)` and BOTH
boundary terms vanish — at `1` because `log 1 = 0`, at `0` because
`log(1+x) = O(x)` beats `log²x`.  So `K` equals the negative of the other piece,
and the expensive route through
`MeasureTheory.hasSum_integral_of_summable_integral_norm` is not needed here. -/

/-- `log(1+x) ≤ x` for `x ≥ 0`. -/
theorem log_one_add_le_self {x : ℝ} (hx : 0 ≤ x) : Real.log (1 + x) ≤ x := by
  have h := Real.log_le_sub_one_of_pos (show (0:ℝ) < 1 + x by linarith)
  simpa using h

/-- `0 ≤ log(1+x)` for `x ≥ 0`. -/
theorem log_one_add_nonneg {x : ℝ} (hx : 0 ≤ x) : 0 ≤ Real.log (1 + x) :=
  Real.log_nonneg (by linarith)

/-- The IBP antiderivative `F = log²x · log²(1+x) / 2`. -/
noncomputable def quadAltKF (x : ℝ) : ℝ :=
  Real.log x ^ 2 * Real.log (1 + x) ^ 2 / 2

theorem quadAltKF_hasDerivAt {x : ℝ} (hx0 : 0 < x) :
    HasDerivAt quadAltKF
      (Real.log x * Real.log (1 + x) ^ 2 / x
        + Real.log x ^ 2 * Real.log (1 + x) / (1 + x)) x := by
  have hx1 : (0:ℝ) < 1 + x := by linarith
  have hlx : HasDerivAt (fun y : ℝ => Real.log y) (1 / x) x := by
    simpa [one_div] using Real.hasDerivAt_log (ne_of_gt hx0)
  have hlp : HasDerivAt (fun y : ℝ => Real.log (1 + y)) (1 / (1 + x)) x := by
    have hinner : HasDerivAt (fun y : ℝ => 1 + y) 1 x := by
      simpa using (hasDerivAt_id x).const_add (1 : ℝ)
    have := (Real.hasDerivAt_log (ne_of_gt hx1)).comp x hinner
    simpa [one_div] using this
  have hmain := (((hlx.pow 2).mul (hlp.pow 2)).div_const 2)
  convert hmain using 1
  simp only [Pi.pow_apply]
  have hxne : x ≠ 0 := ne_of_gt hx0
  have hx1ne : (1:ℝ) + x ≠ 0 := ne_of_gt hx1
  field_simp
  ring

/-- Left endpoint: `F → 0` as `x → 0⁺`, because `log(1+x) = O(x)` beats `log²x`. -/
theorem quadAltKF_tendsto_zero_right :
    Tendsto quadAltKF (𝓝[>] (0:ℝ)) (𝓝 0) := by
  have hsq : Tendsto (fun x : ℝ => (Real.log x ^ 2 * x) * x / 2) (𝓝[>] (0:ℝ)) (𝓝 0) := by
    have h1 := logSq_mul_self_tendsto
    have h2 : Tendsto (fun x : ℝ => x) (𝓝[>] (0:ℝ)) (𝓝 0) :=
      tendsto_id.mono_left nhdsWithin_le_nhds
    simpa using ((h1.mul h2).div_const 2)
  refine squeeze_zero_norm' ?_ hsq
  filter_upwards [self_mem_nhdsWithin] with x hx
  have hx0 : (0:ℝ) < x := hx
  have hle : Real.log (1 + x) ≤ x := log_one_add_le_self (le_of_lt hx0)
  have hnn : 0 ≤ Real.log (1 + x) := log_one_add_nonneg (le_of_lt hx0)
  have hsqle : Real.log (1 + x) ^ 2 ≤ x ^ 2 := by nlinarith
  have : ‖quadAltKF x‖ = Real.log x ^ 2 * Real.log (1 + x) ^ 2 / 2 := by
    unfold quadAltKF
    rw [Real.norm_eq_abs, abs_of_nonneg (by positivity)]
  rw [this]
  have hpos : (0:ℝ) ≤ Real.log x ^ 2 := sq_nonneg _
  have : Real.log x ^ 2 * Real.log (1 + x) ^ 2 ≤ Real.log x ^ 2 * x ^ 2 :=
    mul_le_mul_of_nonneg_left hsqle hpos
  nlinarith [this]

/-- Right endpoint: `F → 0` as `x → 1⁻`, because `log 1 = 0`. -/
theorem quadAltKF_tendsto_zero_left :
    Tendsto quadAltKF (𝓝[<] (1:ℝ)) (𝓝 0) := by
  have hc : ContinuousAt quadAltKF 1 := by
    unfold quadAltKF
    have h1 : ContinuousAt Real.log 1 := Real.continuousAt_log (by norm_num)
    have h2 : ContinuousAt (fun y : ℝ => Real.log (1 + y)) 1 := by
      apply ContinuousAt.log (by fun_prop); norm_num
    exact ((h1.pow 2).mul (h2.pow 2)).div_const 2
  have hval : quadAltKF 1 = 0 := by unfold quadAltKF; simp
  have := hc.tendsto
  rw [hval] at this
  exact this.mono_left nhdsWithin_le_nhds

/-- `1 + log²x` is interval-integrable on `[0,1]`; it dominates both IBP pieces. -/
theorem intervalIntegrable_one_add_logSq :
    IntervalIntegrable (fun x : ℝ => 1 + Real.log x ^ 2) MeasureTheory.volume 0 1 :=
  (intervalIntegrable_const (c := (1:ℝ))).add intervalIntegrable_logSq

/-- The `A` piece `log x · log²(1+x) / x`, dominated by `|log x| ≤ 1 + log²x`. -/
theorem quadAltKA_intervalIntegrable :
    IntervalIntegrable (fun x : ℝ => Real.log x * Real.log (1 + x) ^ 2 / x)
      MeasureTheory.volume 0 1 := by
  refine intervalIntegrable_of_continuousOn_Ioo_of_le (by norm_num) ?_
    intervalIntegrable_one_add_logSq ?_
  · intro x hx
    have hx0 : (0:ℝ) < x := hx.1
    refine ContinuousAt.continuousWithinAt ?_
    have h1 : ContinuousAt (fun y : ℝ => Real.log y) x := Real.continuousAt_log (ne_of_gt hx0)
    have h2 : ContinuousAt (fun y : ℝ => Real.log (1 + y)) x := by
      apply ContinuousAt.log (by fun_prop); linarith
    exact (h1.mul (h2.pow 2)).div continuousAt_id (ne_of_gt hx0)
  · intro x hx
    have hx0 : (0:ℝ) < x := hx.1
    have hx1 : x < 1 := hx.2
    have hle : Real.log (1 + x) ≤ x := log_one_add_le_self (le_of_lt hx0)
    have hnn : 0 ≤ Real.log (1 + x) := log_one_add_nonneg (le_of_lt hx0)
    have hsq : Real.log (1 + x) ^ 2 ≤ x ^ 2 := by nlinarith
    have habs : |Real.log x * Real.log (1 + x) ^ 2 / x|
        = |Real.log x| * (Real.log (1 + x) ^ 2 / x) := by
      rw [abs_div, abs_mul, abs_of_pos hx0,
        abs_of_nonneg (sq_nonneg (Real.log (1 + x)))]
      ring
    rw [Real.norm_eq_abs, habs]
    have hfrac : Real.log (1 + x) ^ 2 / x ≤ 1 := by
      rw [div_le_one hx0]; nlinarith
    have hfnn : 0 ≤ Real.log (1 + x) ^ 2 / x := by positivity
    have hlog : |Real.log x| ≤ 1 + Real.log x ^ 2 := by
      nlinarith [sq_abs (Real.log x), sq_nonneg (|Real.log x| - 1)]
    nlinarith [abs_nonneg (Real.log x)]

/-- The `B` piece is the `K` integrand itself, dominated by `log 2 · log²x`. -/
theorem quadAltKB_intervalIntegrable :
    IntervalIntegrable (fun x : ℝ => Real.log x ^ 2 * Real.log (1 + x) / (1 + x))
      MeasureTheory.volume 0 1 := by
  refine intervalIntegrable_of_continuousOn_Ioo_of_le (by norm_num) ?_
    intervalIntegrable_one_add_logSq ?_
  · intro x hx
    have hx0 : (0:ℝ) < x := hx.1
    refine ContinuousAt.continuousWithinAt ?_
    have h1 : ContinuousAt (fun y : ℝ => Real.log y) x := Real.continuousAt_log (ne_of_gt hx0)
    have h2 : ContinuousAt (fun y : ℝ => Real.log (1 + y)) x := by
      apply ContinuousAt.log (by fun_prop); linarith
    exact ((h1.pow 2).mul h2).div (by fun_prop) (by linarith)
  · intro x hx
    have hx0 : (0:ℝ) < x := hx.1
    have hx1 : x < 1 := hx.2
    have hnn : 0 ≤ Real.log (1 + x) := log_one_add_nonneg (le_of_lt hx0)
    have hle : Real.log (1 + x) ≤ x := log_one_add_le_self (le_of_lt hx0)
    have h1x : (0:ℝ) < 1 + x := by linarith
    have habs : |Real.log x ^ 2 * Real.log (1 + x) / (1 + x)|
        = Real.log x ^ 2 * (Real.log (1 + x) / (1 + x)) := by
      rw [abs_div, abs_mul, abs_of_nonneg (sq_nonneg (Real.log x)),
        abs_of_nonneg hnn, abs_of_pos h1x]
      ring
    rw [Real.norm_eq_abs, habs]
    have hfrac : Real.log (1 + x) / (1 + x) ≤ 1 := by
      rw [div_le_one h1x]; linarith
    have hfnn : 0 ≤ Real.log (1 + x) / (1 + x) := by positivity
    nlinarith [sq_nonneg (Real.log x)]

/-- **The IBP step for `K`.**  Both boundary terms of
`d[log²x · log²(1+x)/2] = log x log²(1+x)/x + log²x log(1+x)/(1+x)`
vanish, so the two pieces are negatives of each other. -/
theorem quadAltK_eq_neg_integral :
    quadAltK = -∫ x in (0:ℝ)..1, Real.log x * Real.log (1 + x) ^ 2 / x := by
  have hsum : IntervalIntegrable
      (fun x : ℝ => Real.log x * Real.log (1 + x) ^ 2 / x
        + Real.log x ^ 2 * Real.log (1 + x) / (1 + x)) MeasureTheory.volume 0 1 :=
    quadAltKA_intervalIntegrable.add quadAltKB_intervalIntegrable
  have hzero := intervalIntegral.integral_eq_sub_of_hasDerivAt_of_tendsto
    (by norm_num : (0:ℝ) < 1)
    (fun x hx => quadAltKF_hasDerivAt hx.1) hsum
    quadAltKF_tendsto_zero_right quadAltKF_tendsto_zero_left
  rw [intervalIntegral.integral_add quadAltKA_intervalIntegrable
    quadAltKB_intervalIntegrable] at hzero
  simp only [sub_zero] at hzero
  unfold quadAltK
  linarith [hzero]

/-- The second logarithmic moment `∫₀¹ xⁿ log²x dx = 2/(n+1)³`.

Companion to `P26.integral_pow_mul_log26` (`∫₀¹ xⁿ log x = -1/(n+1)²`), proved
the same way: an explicit antiderivative plus the tendsto form of FTC. -/
theorem integral_pow_mul_logSq (n : ℕ) :
    (∫ x : ℝ in (0:ℝ)..1, x ^ n * Real.log x ^ 2) = 2 / ((n : ℝ) + 1) ^ 3 := by
  have hn1 : (0:ℝ) < (n : ℝ) + 1 := by positivity
  have hn1ne : ((n : ℝ) + 1) ≠ 0 := ne_of_gt hn1
  set c : ℝ := (n : ℝ) + 1 with hc
  let F : ℝ → ℝ := fun x =>
    x ^ (n + 1) * (Real.log x ^ 2 / c - 2 * Real.log x / c ^ 2 + 2 / c ^ 3)
  have hInt : IntervalIntegrable (fun x : ℝ => x ^ n * Real.log x ^ 2)
      MeasureTheory.volume 0 1 :=
    intervalIntegrable_logSq.continuousOn_mul (continuousOn_pow n)
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt_of_tendsto
    (f := F) (fa := 0) (fb := 2 / c ^ 3) (hint := hInt)]
  · ring
  · norm_num
  · intro x hx
    have hx0 : x ≠ 0 := ne_of_gt hx.1
    have hpow : HasDerivAt (fun y : ℝ => y ^ (n + 1)) (c * x ^ n) x := by
      convert hasDerivAt_pow (n + 1) x using 1
      simp [hc]
    have hlog : HasDerivAt Real.log x⁻¹ x := Real.hasDerivAt_log hx0
    have hrest : HasDerivAt
        (fun y : ℝ => Real.log y ^ 2 / c - 2 * Real.log y / c ^ 2 + 2 / c ^ 3)
        ((2 * Real.log x * x⁻¹) / c - (2 * x⁻¹) / c ^ 2) x := by
      have h1 : HasDerivAt (fun y : ℝ => Real.log y ^ 2 / c)
          ((2 * Real.log x ^ (2 - 1) * x⁻¹) / c) x := (hlog.pow 2).div_const c
      have h2 : HasDerivAt (fun y : ℝ => 2 * Real.log y / c ^ 2)
          ((2 * x⁻¹) / c ^ 2) x := (hlog.const_mul 2).div_const (c ^ 2)
      have := (h1.sub h2).add_const (2 / c ^ 3)
      simpa using this
    dsimp only [F]
    convert hpow.mul hrest using 1
    field_simp
    ring
  · -- left endpoint: xⁿ⁺¹·(log²x/c − 2log x/c² + 2/c³) → 0
    have hx0 : Tendsto (fun x : ℝ => x ^ n) (𝓝[>] (0:ℝ)) (𝓝 ((0:ℝ) ^ n)) :=
      ((continuous_pow n).continuousAt.tendsto).mono_left nhdsWithin_le_nhds
    have hA : Tendsto (fun x : ℝ => (Real.log x ^ 2 * x) * x ^ n / c)
        (𝓝[>] (0:ℝ)) (𝓝 0) := by
      simpa using ((logSq_mul_self_tendsto.mul hx0).div_const c)
    have hB : Tendsto (fun x : ℝ => 2 * ((Real.log x * x) * x ^ n) / c ^ 2)
        (𝓝[>] (0:ℝ)) (𝓝 0) := by
      simpa using (((log_mul_self_tendsto.mul hx0).const_mul 2).div_const (c ^ 2))
    have hC : Tendsto (fun x : ℝ => 2 * (x * x ^ n) / c ^ 3)
        (𝓝[>] (0:ℝ)) (𝓝 0) := by
      have h : Tendsto (fun x : ℝ => x) (𝓝[>] (0:ℝ)) (𝓝 0) :=
        tendsto_id.mono_left nhdsWithin_le_nhds
      simpa using (((h.mul hx0).const_mul 2).div_const (c ^ 3))
    have := (hA.sub hB).add hC
    simp only [sub_zero, add_zero] at this
    refine this.congr ?_
    intro x
    dsimp only [F]
    ring
  · have hcont : ContinuousAt F 1 := by
      dsimp only [F]
      have h1 : ContinuousAt Real.log 1 := Real.continuousAt_log (by norm_num)
      exact (continuousAt_pow _ _).mul
        ((((h1.pow 2).div_const c).sub ((h1.const_mul 2).div_const (c ^ 2))).add
          continuousAt_const)
    have hval : F 1 = 2 / c ^ 3 := by dsimp only [F]; simp
    have := hcont.tendsto
    rw [hval] at this
    exact this.mono_left nhdsWithin_le_nhds

/-! ## Layer E, step 2: `K` as an alternating linear Euler series

The repo already proves `harmonicNumber_generating_hasSum`, and at `x = -t` that
IS the generating function of the `K` integrand:
`log(1+t)/(1+t) = ∑_{n≥0} (-1)ⁿ H_{n+1} tⁿ⁺¹`.  Integrating termwise against
`∫₀¹ x^{n+1} log²x = 2/(n+2)³` therefore gives `K` with no new machinery — the
same `hasSum_integral_of_summable_integral_norm` template this development
already uses three times. -/

/-- The moment family whose termwise integral produces `K`. -/
noncomputable def quadAltKMoment (n : ℕ) (x : ℝ) : ℝ :=
  (-1 : ℝ) ^ n * harmonicNumber (n + 1) * (x ^ (n + 1) * Real.log x ^ 2)

theorem quadAltKMoment_intervalIntegrable (n : ℕ) :
    IntervalIntegrable (quadAltKMoment n) MeasureTheory.volume 0 1 := by
  have h : IntervalIntegrable (fun x : ℝ => x ^ (n + 1) * Real.log x ^ 2)
      MeasureTheory.volume 0 1 :=
    intervalIntegrable_logSq.continuousOn_mul (continuousOn_pow (n + 1))
  unfold quadAltKMoment
  simpa [mul_assoc] using h.const_mul ((-1 : ℝ) ^ n * harmonicNumber (n + 1))

theorem quadAltKMoment_integral (n : ℕ) :
    (∫ x : ℝ in (0:ℝ)..1, quadAltKMoment n x)
      = 2 * (-1 : ℝ) ^ n * harmonicNumber (n + 1) / ((n : ℝ) + 2) ^ 3 := by
  have hm := integral_pow_mul_logSq (n + 1)
  have hcast : ((n + 1 : ℕ) : ℝ) + 1 = (n : ℝ) + 2 := by push_cast; ring
  rw [hcast] at hm
  unfold quadAltKMoment
  rw [intervalIntegral.integral_const_mul, hm]
  ring

theorem quadAltKMoment_hasSum_pointwise {x : ℝ} (hx0 : 0 < x) (hx1 : x < 1) :
    HasSum (fun n : ℕ => quadAltKMoment n x)
      (Real.log x ^ 2 * Real.log (1 + x) / (1 + x)) := by
  have hx1p : (0:ℝ) < 1 + x := by linarith
  have habs : |(-x)| < 1 := by rw [abs_neg, abs_of_pos hx0]; exact hx1
  have hne : (-x) ≠ 0 := by simpa using ne_of_gt hx0
  have hgen := harmonicNumber_generating_hasSum habs hne
  have hmul := hgen.mul_left (x * Real.log x ^ 2)
  have hval : (x * Real.log x ^ 2) *
      (-Real.log (1 - -x) / ((-x) * (1 - -x)))
      = Real.log x ^ 2 * Real.log (1 + x) / (1 + x) := by
    rw [show (1 : ℝ) - -x = 1 + x by ring]
    field_simp
  rw [hval] at hmul
  refine hmul.congr_fun ?_
  intro n
  unfold quadAltKMoment
  rw [neg_pow]
  ring

theorem quadAltKMoment_integral_norm_summable :
    Summable (fun n : ℕ => ∫ x : ℝ in (0:ℝ)..1, ‖quadAltKMoment n x‖) := by
  have hval : ∀ n : ℕ, (∫ x : ℝ in (0:ℝ)..1, ‖quadAltKMoment n x‖)
      = 2 * harmonicNumber (n + 1) / ((n : ℝ) + 2) ^ 3 := by
    intro n
    have hcongr : ∀ x ∈ Set.uIcc (0:ℝ) 1,
        ‖quadAltKMoment n x‖ = harmonicNumber (n + 1) * (x ^ (n + 1) * Real.log x ^ 2) := by
      intro x hx
      rw [Set.uIcc_of_le (by norm_num : (0:ℝ) ≤ 1)] at hx
      unfold quadAltKMoment
      rw [Real.norm_eq_abs, abs_mul, abs_mul, abs_pow, abs_neg, abs_one, one_pow,
        one_mul, abs_of_nonneg (harmonicNumber_nonneg (n + 1)),
        abs_of_nonneg (mul_nonneg (pow_nonneg hx.1 (n + 1)) (sq_nonneg _))]
    rw [intervalIntegral.integral_congr hcongr, intervalIntegral.integral_const_mul]
    have hm := integral_pow_mul_logSq (n + 1)
    have hcast : ((n + 1 : ℕ) : ℝ) + 1 = (n : ℝ) + 2 := by push_cast; ring
    rw [hcast] at hm
    rw [hm]; ring
  refine Summable.congr ?_ (fun n => (hval n).symm)
  have hbd : ∀ n : ℕ, ‖2 * harmonicNumber (n + 1) / ((n : ℝ) + 2) ^ 3‖
      ≤ 2 * (harmonicNumber (n + 1) ^ 2 / ((n : ℝ) + 1) ^ 2) := by
    intro n
    have hH : 0 ≤ harmonicNumber (n + 1) := harmonicNumber_nonneg (n + 1)
    have hHsq : harmonicNumber (n + 1) ≤ harmonicNumber (n + 1) ^ 2 :=
      harmonicNumber_succ_le_sq n
    rw [Real.norm_eq_abs, abs_of_nonneg (by positivity)]
    rw [div_le_iff₀ (by positivity)]
    have h1 : ((n : ℝ) + 1) ^ 2 ≤ ((n : ℝ) + 2) ^ 3 := by
      nlinarith [Nat.cast_nonneg (α := ℝ) n]
    have hq : (0:ℝ) ≤ 2 * (harmonicNumber (n + 1) ^ 2 / ((n : ℝ) + 1) ^ 2) := by
      positivity
    have step1 : 2 * (harmonicNumber (n + 1) ^ 2 / ((n : ℝ) + 1) ^ 2) * ((n : ℝ) + 1) ^ 2
        ≤ 2 * (harmonicNumber (n + 1) ^ 2 / ((n : ℝ) + 1) ^ 2) * ((n : ℝ) + 2) ^ 3 :=
      mul_le_mul_of_nonneg_left h1 hq
    have step2 : 2 * (harmonicNumber (n + 1) ^ 2 / ((n : ℝ) + 1) ^ 2) * ((n : ℝ) + 1) ^ 2
        = 2 * harmonicNumber (n + 1) ^ 2 := by
      field_simp
    linarith
  exact (summable_harmonicNumber_succ_sq_div.mul_left 2).of_norm_bounded hbd

/-- **The `K` series.**  `K = 2 ∑ (-1)ⁿ H_{n+1}/(n+2)³`. -/
theorem quadAltK_hasSum :
    HasSum (fun n : ℕ => 2 * (-1 : ℝ) ^ n * harmonicNumber (n + 1) / ((n : ℝ) + 2) ^ 3)
      quadAltK := by
  have hInt : ∀ n : ℕ, MeasureTheory.Integrable (quadAltKMoment n)
      (MeasureTheory.volume.restrict (Set.Ioc (0:ℝ) 1)) :=
    fun n => (quadAltKMoment_intervalIntegrable n).1
  have hNorm : Summable (fun n : ℕ =>
      ∫ x : ℝ in Set.Ioc (0:ℝ) 1, ‖quadAltKMoment n x‖) := by
    simpa only [← intervalIntegral.integral_of_le (by norm_num : (0:ℝ) ≤ 1)] using
      quadAltKMoment_integral_norm_summable
  have h := MeasureTheory.hasSum_integral_of_summable_integral_norm
    (μ := MeasureTheory.volume.restrict (Set.Ioc (0:ℝ) 1)) hInt hNorm
  have h' : HasSum (fun n : ℕ =>
      2 * (-1 : ℝ) ^ n * harmonicNumber (n + 1) / ((n : ℝ) + 2) ^ 3)
      (∫ x : ℝ in Set.Ioc (0:ℝ) 1, ∑' n : ℕ, quadAltKMoment n x) := by
    convert h using 1
    funext n
    rw [← intervalIntegral.integral_of_le (by norm_num : (0:ℝ) ≤ 1)]
    exact (quadAltKMoment_integral n).symm
  convert h' using 1
  unfold quadAltK
  rw [intervalIntegral.integral_of_le (by norm_num : (0:ℝ) ≤ 1)]
  apply MeasureTheory.setIntegral_congr_ae measurableSet_Ioc
  filter_upwards [MeasureTheory.Measure.ae_ne MeasureTheory.volume (1:ℝ)] with x hxne hx
  have hxlt : x < 1 := lt_of_le_of_ne hx.2 hxne
  exact (quadAltKMoment_hasSum_pointwise hx.1 hxlt).tsum_eq.symm

/-- **Layer E, step 3: the value of `K`.**

Purely a series rearrangement: with
`a n = (-1)^{n+1} H_{n+1}/(n+1)³` and `b n = (-1)^{n+1}/(n+1)⁴`, both of which
this development already evaluates in closed form, the `K` summand is exactly
`2(a(n+1) - b(n+1))` — the `1/(n+2)⁴` produced by `H_{n+2} = H_{n+1} + 1/(n+2)`
is cancelled by `b`.  No new analysis. -/
theorem quadAltK_eq :
    quadAltK = (1 / 5) * (Real.pi ^ 2 / 6) ^ 2
      - 2 * alternatingCubicLinearEulerValue24 - 2 * cubicLinearEulerValue24 := by
  have hHsucc : ∀ m : ℕ, harmonicNumber (m + 1) = harmonicNumber m + 1 / ((m : ℝ) + 1) := by
    intro m; simp [harmonicNumber, Finset.sum_range_succ]
  set Va : ℝ := 2 * polylog4 (1 / 2) + (1 / 12 : ℝ) * Real.log 2 ^ 4 -
      (1 / 2 : ℝ) * Real.log 2 ^ 2 * (Real.pi ^ 2 / 6) +
      (7 / 4 : ℝ) * Real.log 2 * zeta3_24 -
      (11 / 4 : ℝ) * (Real.pi ^ 4 / 90) with hVa
  set Vb : ℝ := -(7 / 8 : ℝ) * (Real.pi ^ 4 / 90) with hVb
  have ha : HasSum
      (fun n : ℕ => (-1 : ℝ) ^ (n + 1) * harmonicNumber (n + 1) / ((n : ℝ) + 1) ^ 3) Va :=
    alternatingHarmonicCubic_hasSum24
  have hb : HasSum (fun n : ℕ => (-1 : ℝ) ^ (n + 1) / ((n : ℝ) + 1) ^ 4) Vb :=
    alternatingZetaFour_hasSum24
  have hc : HasSum (fun n : ℕ =>
      (-1 : ℝ) ^ (n + 1) * harmonicNumber (n + 1) / ((n : ℝ) + 1) ^ 3
        - (-1 : ℝ) ^ (n + 1) / ((n : ℝ) + 1) ^ 4) (Va - Vb) := ha.sub hb
  have hc1 := (hasSum_nat_add_iff' (f := fun n : ℕ =>
      (-1 : ℝ) ^ (n + 1) * harmonicNumber (n + 1) / ((n : ℝ) + 1) ^ 3
        - (-1 : ℝ) ^ (n + 1) / ((n : ℝ) + 1) ^ 4) 1).mpr hc
  have hzero : (∑ i ∈ Finset.range 1,
      ((-1 : ℝ) ^ (i + 1) * harmonicNumber (i + 1) / ((i : ℝ) + 1) ^ 3
        - (-1 : ℝ) ^ (i + 1) / ((i : ℝ) + 1) ^ 4)) = 0 := by
    norm_num [harmonicNumber]
  rw [hzero, sub_zero] at hc1
  have hcomb := hc1.mul_left 2
  have hkey : HasSum
      (fun n : ℕ => 2 * (-1 : ℝ) ^ n * harmonicNumber (n + 1) / ((n : ℝ) + 2) ^ 3)
      (2 * (Va - Vb)) := by
    refine hcomb.congr_fun ?_
    intro n
    have hH : harmonicNumber (n + 1 + 1) = harmonicNumber (n + 1) + 1 / ((n : ℝ) + 2) := by
      rw [hHsucc (n + 1)]; push_cast; ring_nf
    have hsign : (-1 : ℝ) ^ (n + 1 + 1) = (-1 : ℝ) ^ n := by
      rw [pow_succ, pow_succ]; ring
    have hne : ((n : ℝ) + 2) ≠ 0 := by positivity
    rw [show ((n + 1 : ℕ) : ℝ) + 1 = (n : ℝ) + 2 by push_cast; ring, hH, hsign]
    field_simp
    ring
  have := quadAltK_hasSum.unique hkey
  rw [this, hVa, hVb]
  unfold alternatingCubicLinearEulerValue24 cubicLinearEulerValue24
  ring

/-! ## Layer E, row `I11`

Reflection turns the post-IBP kernel into the difference between the ordinary
quartic core and its alternating complement.  Both values are supplied by the
series library in `Problem24Euler`. -/

theorem quadAltI11_reflected :
    I11 = ∫ x : ℝ in (0 : ℝ)..1,
      (Real.log x ^ 2 * Real.log (1 - x) / (1 - x) -
        Real.log x ^ 2 * Real.log (1 + x) / (1 - x)) := by
  rw [quadAltI11_eq_integral]
  let f : ℝ → ℝ := fun t => quadAltR t * H1 t ^ 2 / t
  have hreflect := intervalIntegral.integral_comp_sub_left
    (a := (0 : ℝ)) (b := 1) f 1
  calc
    (∫ t : ℝ in (0 : ℝ)..1, quadAltR t * H1 t ^ 2 / t) =
        ∫ x : ℝ in (0 : ℝ)..1, f (1 - x) := by
      simpa [f] using hreflect.symm
    _ = ∫ x : ℝ in (0 : ℝ)..1,
        (Real.log x ^ 2 * Real.log (1 - x) / (1 - x) -
          Real.log x ^ 2 * Real.log (1 + x) / (1 - x)) := by
      apply intervalIntegral.integral_congr
      intro x hx
      rw [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)] at hx
      by_cases hxzero : x = 0
      · subst x
        simp [f, quadAltR, H1]
      by_cases hxone : x = 1
      · subst x
        simp [f, quadAltR, H1]
      have hx0 : 0 < x := lt_of_le_of_ne hx.1 (Ne.symm hxzero)
      have hx1 : x < 1 := lt_of_le_of_ne hx.2 hxone
      have h1xne : 1 - x ≠ 0 := ne_of_gt (sub_pos.mpr hx1)
      have h1pxne : 1 + x ≠ 0 := by linarith
      unfold f quadAltR H1
      dsimp only
      rw [show 2 - (1 - x) = 1 + x by ring,
        show 1 - (1 - x) = x by ring,
        Real.log_div h1xne h1pxne]
      field_simp [h1xne]

theorem quadAltI11_eq :
    I11 = -(7 / 2) * Real.log 2 * zeta3_24
      + (3 / 4) * (Real.pi ^ 2 / 6) ^ 2 := by
  let core : ℝ → ℝ := fun x =>
    Real.log x * Real.log (1 - x) ^ 2 / x
  let first : ℝ → ℝ := fun x =>
    Real.log x ^ 2 * Real.log (1 - x) / (1 - x)
  let alt : ℝ → ℝ := fun x =>
    Real.log x ^ 2 * Real.log (1 + x) / (1 - x)
  have hFirstInt : IntervalIntegrable first MeasureTheory.volume 0 1 := by
    have h :=
      (quarticCoreIntervalIntegrable24_export.comp_sub_left 1).symm
    convert h using 1
    funext x
    unfold first
    rw [show 1 - (1 - x) = x by ring]
    ring
    all_goals norm_num
  have hAltInt : IntervalIntegrable alt MeasureTheory.volume 0 1 := by
    exact quarticAlternatingComplementIntervalIntegrable24
  have hFirst : (∫ x : ℝ in (0 : ℝ)..1, first x) =
      -(1 / 2 : ℝ) * (Real.pi ^ 4 / 90) := by
    have hreflect := intervalIntegral.integral_comp_sub_left
      (a := (0 : ℝ)) (b := 1) core 1
    calc
      (∫ x : ℝ in (0 : ℝ)..1, first x) =
          ∫ x : ℝ in (0 : ℝ)..1, core (1 - x) := by
        apply intervalIntegral.integral_congr
        intro x _
        unfold first core
        dsimp only
        rw [show 1 - (1 - x) = x by ring]
        ring
      _ = ∫ x : ℝ in (0 : ℝ)..1, core x := by
        simpa using hreflect
      _ = -(1 / 2 : ℝ) * (Real.pi ^ 4 / 90) := by
        exact quarticCoreIntegral24_export
  rw [quadAltI11_reflected]
  change (∫ x : ℝ in (0 : ℝ)..1, first x - alt x) = _
  rw [intervalIntegral.integral_sub hFirstInt hAltInt,
    hFirst, quarticAlternatingComplementIntegral24]
  ring

/-! ## Layer E, row `I20` -/

theorem quadAltI20_eq :
    I20 = -(1 / 2) * quadAltK -
      (1 / 2) * (Real.pi ^ 2 / 6) ^ 2 := by
  let value : ℝ :=
    -2 * polylog4 (1 / 2) -
      (1 / 12 : ℝ) * Real.log 2 ^ 4 +
      (1 / 2 : ℝ) * Real.log 2 ^ 2 * (Real.pi ^ 2 / 6) -
      (7 / 4 : ℝ) * Real.log 2 * zeta3_24 +
      (1 / 4 : ℝ) * (Real.pi ^ 2 / 6) ^ 2
  have hint : IntervalIntegrable
      (fun t : ℝ =>
        (Real.pi ^ 2 / 6 - 2 * dilog (t / 2) -
          Real.log (t / 2) ^ 2) *
        (-Real.log (1 - t / 2)) / t)
      MeasureTheory.volume 0 1 := by
    simpa [W0, H2] using quadAltI20_kernel_intervalIntegrable
  have hvalue :
      (∫ t : ℝ in (0 : ℝ)..1,
        (Real.pi ^ 2 / 6 - 2 * dilog (t / 2) -
          Real.log (t / 2) ^ 2) *
        (-Real.log (1 - t / 2)) / t) = value := by
    rw [intervalIntegral.integral_eq_sub_of_hasDerivAt_of_tendsto
      (f := quadAltI20Primitive24)
      (fa := (0 : ℝ)) (fb := value)
      (by norm_num)
      (fun t ht => quadAltI20Primitive24_hasDerivAt ht.1 ht.2)
      hint quadAltI20Primitive24_tendsto_zero
      (by simpa [value] using quadAltI20Primitive24_tendsto_one)]
    ring
  calc
    I20 = value := by
      unfold I20 W0 H2
      exact hvalue
    _ = -(1 / 2) * quadAltK -
        (1 / 2) * (Real.pi ^ 2 / 6) ^ 2 := by
      rw [quadAltK_eq]
      unfold value alternatingCubicLinearEulerValue24
        cubicLinearEulerValue24
      ring

/-! ## Layer E, row `I22` -/

private theorem quadAltI22_radial_intervalIntegrable :
    IntervalIntegrable
      (fun t : ℝ => quadAltR t * H2 t ^ 2 / t)
      MeasureTheory.volume 0 1 := by
  let half : ℝ → ℝ := fun u =>
    (Real.log u - Real.log (1 - u)) *
      Real.log (1 - u) ^ 2 / u
  have hhalf : IntervalIntegrable half MeasureTheory.volume 0 (1 / 2) := by
    exact quadAltI22HalfRadialIntervalIntegrable24
  have hcomp := hhalf.comp_mul_left (c := (1 / 2 : ℝ))
  have hcomp' : IntervalIntegrable (fun t : ℝ => half (t / 2))
      MeasureTheory.volume 0 1 := by
    convert hcomp using 1 <;> norm_num
    funext t
    ring
  apply IntervalIntegrable.congr
    (f := fun t : ℝ => (1 / 2 : ℝ) * half (t / 2)) ?_
    (hcomp'.const_mul (1 / 2 : ℝ))
  intro t ht
  rw [Set.uIoc_of_le (by norm_num : (0 : ℝ) ≤ 1)] at ht
  have htzero : t ≠ 0 := ne_of_gt ht.1
  have hthalf : t / 2 ≠ 0 := div_ne_zero htzero (by norm_num)
  have h2t : 2 - t ≠ 0 := by linarith [ht.2]
  have hhalfden : 1 - t / 2 ≠ 0 := by linarith [ht.2]
  unfold half quadAltR H2
  dsimp only
  rw [show t / (2 - t) = (t / 2) / (1 - t / 2) by
    field_simp [h2t, hhalfden],
    Real.log_div hthalf hhalfden]
  field_simp [htzero]

private theorem quadAltI22_radial_eq_half :
    (∫ t : ℝ in (0 : ℝ)..1,
      quadAltR t * H2 t ^ 2 / t) =
      ∫ u : ℝ in (0 : ℝ)..(1 / 2),
        (Real.log u - Real.log (1 - u)) *
          Real.log (1 - u) ^ 2 / u := by
  let radial : ℝ → ℝ := fun t => quadAltR t * H2 t ^ 2 / t
  have hscale := intervalIntegral.smul_integral_comp_mul_left
    (a := (0 : ℝ)) (b := (1 / 2 : ℝ)) radial (2 : ℝ)
  calc
    (∫ t : ℝ in (0 : ℝ)..1, quadAltR t * H2 t ^ 2 / t) =
        (2 : ℝ) * ∫ u : ℝ in (0 : ℝ)..(1 / 2), radial (2 * u) := by
      simpa [radial, smul_eq_mul] using hscale.symm
    _ = ∫ u : ℝ in (0 : ℝ)..(1 / 2), 2 * radial (2 * u) := by
      rw [← intervalIntegral.integral_const_mul]
    _ = ∫ u : ℝ in (0 : ℝ)..(1 / 2),
        (Real.log u - Real.log (1 - u)) *
          Real.log (1 - u) ^ 2 / u := by
      apply intervalIntegral.integral_congr
      intro u hu
      rw [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1 / 2)] at hu
      by_cases huzero : u = 0
      · subst u
        simp [radial, quadAltR, H2]
      · have hu0 : 0 < u := lt_of_le_of_ne hu.1 (Ne.symm huzero)
        have h1u : 1 - u ≠ 0 := by linarith [hu.2]
        unfold radial quadAltR H2
        dsimp only
        rw [show (2 * u) / (2 - 2 * u) = u / (1 - u) by
          field_simp [h1u],
          Real.log_div (ne_of_gt hu0) h1u]
        field_simp [huzero]

private theorem quadAltW0_mul_g22_tendsto_right :
    Tendsto (fun t : ℝ => W0 t * (H2 t ^ 2 / 2))
      (𝓝[>] (0 : ℝ)) (𝓝 0) := by
  have hratio : Tendsto
      (fun t : ℝ => (H2 t / t) / (H1 t / t))
      (𝓝[>] (0 : ℝ)) (𝓝 (1 / 2)) := by
    simpa using quadAltH2_div_self_tendsto_zero_right.div
      quadAltH1_div_self_tendsto_zero_right (by norm_num)
  have h := quadAltW0_mul_g11_tendsto_right.mul (hratio.pow 2)
  norm_num at h
  refine h.congr' ?_
  filter_upwards [self_mem_nhdsWithin,
    (eventually_lt_nhds (show (0 : ℝ) < 1 by norm_num)).filter_mono
      nhdsWithin_le_nhds] with t ht0 ht1
  have htne : t ≠ 0 := ne_of_gt ht0
  have hH1pos : 0 < H1 t := by
    unfold H1
    have htpos : 0 < t := ht0
    have hsubpos : 0 < 1 - t := sub_pos.mpr ht1
    have hsublt : 1 - t < 1 := sub_lt_self 1 htpos
    exact neg_pos.mpr (Real.log_neg hsubpos hsublt)
  have hH1ne : H1 t ≠ 0 := ne_of_gt hH1pos
  field_simp [htne, hH1ne]

private theorem quadAltW0_mul_g22_tendsto_left :
    Tendsto (fun t : ℝ => W0 t * (H2 t ^ 2 / 2))
      (𝓝[<] (1 : ℝ)) (𝓝 0) := by
  have hW : Tendsto W0 (𝓝[<] (1 : ℝ)) (𝓝 0) := by
    simpa [quadAltW0_one] using
      quadAltW0_hasDerivAt_one.continuousAt.tendsto.mono_left
        nhdsWithin_le_nhds
  have h := hW.mul ((quadAltH2_tendsto_one.pow 2).div_const 2)
  simpa using h

private theorem quadAltI22_eq_radial :
    I22 = ∫ t : ℝ in (0 : ℝ)..1,
      quadAltR t * H2 t ^ 2 / t := by
  let g22 : ℝ → ℝ := fun t => H2 t ^ 2 / 2
  have hg22' : ∀ t ∈ Set.Ioo (0 : ℝ) 1,
      HasDerivAt g22 (H2 t / (2 - t)) t := by
    intro t ht
    have hden : 1 - t / 2 ≠ 0 := by linarith [ht.2]
    have hinner : HasDerivAt (fun y : ℝ => 1 - y / 2)
        (-1 / 2) t := by
      convert (hasDerivAt_const t 1).sub
        ((hasDerivAt_id t).div_const 2) using 1 <;> norm_num
    have hH2 : HasDerivAt H2 (1 / (2 - t)) t := by
      unfold H2
      convert (hinner.log hden).neg using 1
      field_simp [hden]
    unfold g22
    convert (hH2.pow 2).div_const 2 using 1
    field_simp
    ring
  have hprod' : ∀ t ∈ Set.Ioo (0 : ℝ) 1,
      HasDerivAt (fun y : ℝ => W0 y * g22 y)
        (W0 t * (H2 t / (2 - t)) +
          (-2 * (quadAltR t / t) * g22 t)) t := by
    intro t ht
    have hprod :=
      (quadAltW0_hasDerivAt ht.1 (by linarith [ht.2])).mul
        (hg22' t ht)
    unfold g22
    convert hprod using 1
    unfold quadAltR
    ring
  have hA : IntervalIntegrable
      (fun t : ℝ => W0 t * (H2 t / (2 - t)))
      MeasureTheory.volume 0 1 := by
    apply IntervalIntegrable.congr
      (f := fun t : ℝ => W0 t * H2 t / (2 - t)) ?_
      quadAltI22_kernel_intervalIntegrable
    intro t _
    ring
  have hB : IntervalIntegrable
      (fun t : ℝ => -2 * (quadAltR t / t) * g22 t)
      MeasureTheory.volume 0 1 := by
    have hrad := quadAltI22_radial_intervalIntegrable.neg
    apply IntervalIntegrable.congr
      (f := fun t : ℝ => -(quadAltR t * H2 t ^ 2 / t)) ?_ hrad
    intro t _
    unfold g22
    ring
  have hFTC := intervalIntegral.integral_eq_sub_of_hasDerivAt_of_tendsto
    (a := (0 : ℝ)) (b := 1)
    (f := fun t : ℝ => W0 t * g22 t)
    (f' := fun t : ℝ => W0 t * (H2 t / (2 - t)) +
      (-2 * (quadAltR t / t) * g22 t))
    (by norm_num) hprod' (hA.add hB)
    (by simpa [g22] using quadAltW0_mul_g22_tendsto_right)
    (by simpa [g22] using quadAltW0_mul_g22_tendsto_left)
  simp only [sub_self] at hFTC
  rw [intervalIntegral.integral_add hA hB] at hFTC
  have hBint :
      (∫ t : ℝ in (0 : ℝ)..1,
        -2 * (quadAltR t / t) * g22 t) =
        -(∫ t : ℝ in (0 : ℝ)..1,
          quadAltR t * H2 t ^ 2 / t) := by
    rw [← intervalIntegral.integral_neg]
    apply intervalIntegral.integral_congr
    intro t _
    unfold g22
    ring
  rw [hBint] at hFTC
  have hmain :
      (∫ t : ℝ in (0 : ℝ)..1, W0 t * (H2 t / (2 - t))) =
        ∫ t : ℝ in (0 : ℝ)..1, quadAltR t * H2 t ^ 2 / t := by
    linarith
  unfold I22
  calc
    (∫ t : ℝ in (0 : ℝ)..1, W0 t * H2 t / (2 - t)) =
        ∫ t : ℝ in (0 : ℝ)..1, W0 t * (H2 t / (2 - t)) := by
      apply intervalIntegral.integral_congr
      intro t _
      ring
    _ = ∫ t : ℝ in (0 : ℝ)..1,
        quadAltR t * H2 t ^ 2 / t := hmain

theorem quadAltI22_eq :
    I22 = -(3 / 2) * quadAltK +
      (1 / 20) * (Real.pi ^ 2 / 6) ^ 2 := by
  rw [quadAltI22_eq_radial, quadAltI22_radial_eq_half,
    quadAltI22HalfRadialIntegral24, quadAltK_eq]
  unfold alternatingCubicLinearEulerValue24 cubicLinearEulerValue24
  ring

private noncomputable def i10PowerH1Primitive (n : ℕ) (x : ℝ) : ℝ :=
  H1 x * ((x ^ (n + 1) - 1) / ((n : ℝ) + 1)) +
    (1 / ((n : ℝ) + 1)) *
      ∑ k ∈ Finset.range (n + 1), x ^ (k + 1) / ((k : ℝ) + 1)

private theorem i10PowerH1Primitive_hasDerivAt (n : ℕ)
    {x : ℝ} (hx0 : 0 < x) (hx1 : x < 1) :
    HasDerivAt (i10PowerH1Primitive n) (x ^ n * H1 x) x := by
  have hxne : x ≠ 0 := ne_of_gt hx0
  have h1xne : 1 - x ≠ 0 := ne_of_gt (sub_pos.mpr hx1)
  have hH1 : HasDerivAt H1 (1 / (1 - x)) x := by
    have hinner : HasDerivAt (fun y : ℝ => 1 - y) (-1) x := by
      convert (hasDerivAt_const x 1).sub (hasDerivAt_id x) using 1 <;> norm_num
    unfold H1
    convert (hinner.log h1xne).neg using 1
    field_simp [h1xne]
  have hpow : HasDerivAt
      (fun y : ℝ => (y ^ (n + 1) - 1) / ((n : ℝ) + 1))
      (x ^ n) x := by
    convert ((hasDerivAt_pow (n + 1) x).sub_const 1).div_const
      ((n : ℝ) + 1) using 1
    simp only [Nat.cast_add, Nat.cast_one]
    field_simp
    simp
  have hsum : HasDerivAt
      (fun y : ℝ => (1 / ((n : ℝ) + 1)) *
        ∑ k ∈ Finset.range (n + 1), y ^ (k + 1) / ((k : ℝ) + 1))
      ((1 / ((n : ℝ) + 1)) *
        ∑ k ∈ Finset.range (n + 1), x ^ k) x := by
    have heach : ∀ k ∈ Finset.range (n + 1), HasDerivAt
        (fun y : ℝ => y ^ (k + 1) / ((k : ℝ) + 1)) (x ^ k) x := by
      intro k hk
      convert (hasDerivAt_pow (k + 1) x).div_const ((k : ℝ) + 1) using 1
      simp only [Nat.cast_add, Nat.cast_one]
      field_simp
      simp
    have hs := HasDerivAt.sum heach
    convert hs.const_mul (1 / ((n : ℝ) + 1)) using 1 <;>
      simp only [Finset.sum_apply]
  unfold i10PowerH1Primitive
  convert (hH1.mul hpow).add hsum using 1
  have hxone : x ≠ 1 := ne_of_lt hx1
  have hxsub : x - 1 ≠ 0 := sub_ne_zero.mpr hxone
  rw [geom_sum_eq hxone]
  simp only [pow_succ]
  field_simp [h1xne, hxone, hxsub]
  ring_nf

private theorem i10PowerH1Primitive_tendsto_zero (n : ℕ) :
    Tendsto (i10PowerH1Primitive n) (𝓝[>] (0 : ℝ)) (𝓝 0) := by
  have hH : Tendsto H1 (𝓝[>] (0 : ℝ)) (𝓝 0) := by
    have hc : ContinuousAt H1 0 := by
      unfold H1
      exact ((continuousAt_const.sub continuousAt_id).log (by norm_num)).neg
    simpa [H1] using hc.tendsto.mono_left nhdsWithin_le_nhds
  have hpow : Tendsto
      (fun x : ℝ => (x ^ (n + 1) - 1) / ((n : ℝ) + 1))
      (𝓝[>] (0 : ℝ)) (𝓝 (-(1 / ((n : ℝ) + 1)))) := by
    have hc : ContinuousAt
        (fun x : ℝ => (x ^ (n + 1) - 1) / ((n : ℝ) + 1)) 0 := by
      fun_prop
    convert hc.tendsto.mono_left nhdsWithin_le_nhds using 1
    simp [zero_pow (Nat.succ_ne_zero n)]
    ring
  have hsum : Tendsto
      (fun x : ℝ => (1 / ((n : ℝ) + 1)) *
        ∑ k ∈ Finset.range (n + 1), x ^ (k + 1) / ((k : ℝ) + 1))
      (𝓝[>] (0 : ℝ)) (𝓝 0) := by
    have hc : ContinuousAt
        (fun x : ℝ => (1 / ((n : ℝ) + 1)) *
          ∑ k ∈ Finset.range (n + 1), x ^ (k + 1) / ((k : ℝ) + 1)) 0 := by
      fun_prop
    simpa [zero_pow (Nat.succ_ne_zero _)] using
      hc.tendsto.mono_left nhdsWithin_le_nhds
  unfold i10PowerH1Primitive
  simpa using (hH.mul hpow).add hsum

private theorem i10PowerH1Primitive_tendsto_one (n : ℕ) :
    Tendsto (i10PowerH1Primitive n) (𝓝[<] (1 : ℝ))
      (𝓝 (harmonicNumber (n + 1) / ((n : ℝ) + 1))) := by
  have hgeom : Tendsto
      (fun x : ℝ =>
        (∑ k ∈ Finset.range (n + 1), x ^ k) / ((n : ℝ) + 1))
      (𝓝[<] (1 : ℝ)) (𝓝 1) := by
    have hc : ContinuousAt
        (fun x : ℝ =>
          (∑ k ∈ Finset.range (n + 1), x ^ k) / ((n : ℝ) + 1)) 1 := by
      fun_prop
    convert hc.tendsto.mono_left nhdsWithin_le_nhds using 1
    simp
    field_simp
  have hfirst := oneSub_log_tendsto.mul hgeom
  simp only [zero_mul] at hfirst
  have hfirst' : Tendsto
      (fun x : ℝ => H1 x * ((x ^ (n + 1) - 1) / ((n : ℝ) + 1)))
      (𝓝[<] (1 : ℝ)) (𝓝 0) := by
    apply hfirst.congr'
    filter_upwards [self_mem_nhdsWithin] with x hx
    have hxone : x ≠ 1 := ne_of_lt hx
    have hxsub : x - 1 ≠ 0 := sub_ne_zero.mpr hxone
    unfold H1
    rw [geom_sum_eq hxone]
    field_simp [hxone, hxsub]
    ring
  have hsecond : Tendsto
      (fun x : ℝ => (1 / ((n : ℝ) + 1)) *
        ∑ k ∈ Finset.range (n + 1), x ^ (k + 1) / ((k : ℝ) + 1))
      (𝓝[<] (1 : ℝ))
      (𝓝 (harmonicNumber (n + 1) / ((n : ℝ) + 1))) := by
    have hc : ContinuousAt
        (fun x : ℝ => (1 / ((n : ℝ) + 1)) *
          ∑ k ∈ Finset.range (n + 1), x ^ (k + 1) / ((k : ℝ) + 1)) 1 := by
      fun_prop
    convert hc.tendsto.mono_left nhdsWithin_le_nhds using 1
    simp [harmonicNumber]
    ring
  unfold i10PowerH1Primitive
  simpa using hfirst'.add hsecond

private theorem i10PowerH1_intervalIntegrable (n : ℕ) :
    IntervalIntegrable (fun x : ℝ => x ^ n * H1 x)
      MeasureTheory.volume 0 1 := by
  have hlog := (intervalIntegral.intervalIntegrable_log'
    (a := (0 : ℝ)) (b := 1)).comp_sub_left 1
  have hH1 : IntervalIntegrable H1 MeasureTheory.volume 0 1 := by
    simpa [H1] using hlog.symm.neg
  exact hH1.continuousOn_mul (continuousOn_pow n)

theorem i10PowerH1_integral (n : ℕ) :
    (∫ x : ℝ in 0..1, x ^ n * H1 x) =
      harmonicNumber (n + 1) / ((n : ℝ) + 1) := by
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt_of_tendsto
    (f := i10PowerH1Primitive n)
    (fa := (0 : ℝ))
    (fb := harmonicNumber (n + 1) / ((n : ℝ) + 1))
    (by norm_num)
    (fun x hx => i10PowerH1Primitive_hasDerivAt n hx.1 hx.2)
    (i10PowerH1_intervalIntegrable n)
    (i10PowerH1Primitive_tendsto_zero n)
    (i10PowerH1Primitive_tendsto_one n)]
  ring

private noncomputable def i10HalfHarmonicTerm (n : ℕ) : ℝ :=
  harmonicNumber (n + 1) /
    ((2 : ℝ) ^ (n + 1) * ((n : ℝ) + 1) ^ 3)

private theorem i10HalfHarmonicTerm_summable :
    Summable i10HalfHarmonicTerm := by
  apply summable_harmonicNumber_succ_sq_div.of_norm_bounded
  intro n
  have hH : 0 ≤ harmonicNumber (n + 1) := harmonicNumber_nonneg (n + 1)
  have hHsq : harmonicNumber (n + 1) ≤ harmonicNumber (n + 1) ^ 2 :=
    harmonicNumber_succ_le_sq n
  have hpow : (1 : ℝ) ≤ 2 ^ (n + 1) := one_le_pow₀ (by norm_num)
  rw [Real.norm_eq_abs, abs_of_nonneg (by
    unfold i10HalfHarmonicTerm
    positivity)]
  unfold i10HalfHarmonicTerm
  calc
    harmonicNumber (n + 1) /
        ((2 : ℝ) ^ (n + 1) * ((n : ℝ) + 1) ^ 3) ≤
        harmonicNumber (n + 1) / ((n : ℝ) + 1) ^ 3 := by
      apply div_le_div_of_nonneg_left hH (by positivity)
      exact le_mul_of_one_le_left (by positivity) hpow
    _ ≤ harmonicNumber (n + 1) ^ 2 / ((n : ℝ) + 1) ^ 3 := by
      gcongr
    _ ≤ harmonicNumber (n + 1) ^ 2 / ((n : ℝ) + 1) ^ 2 := by
      apply div_le_div_of_nonneg_left (sq_nonneg _) (by positivity)
      have hn0 : (0 : ℝ) ≤ n := Nat.cast_nonneg n
      have hn : (1 : ℝ) ≤ (n : ℝ) + 1 := by linarith
      nlinarith [sq_nonneg ((n : ℝ) + 1)]

private noncomputable def i10CrossMoment (n : ℕ) (x : ℝ) : ℝ :=
  ((x / 2) ^ (n + 1) / ((n : ℝ) + 1) ^ 2) * (H1 x / x)

private theorem i10CrossMoment_intervalIntegrable (n : ℕ) :
    IntervalIntegrable (i10CrossMoment n) MeasureTheory.volume 0 1 := by
  apply IntervalIntegrable.congr
    (f := fun x : ℝ =>
      (1 / ((2 : ℝ) ^ (n + 1) * ((n : ℝ) + 1) ^ 2)) *
        (x ^ n * H1 x)) ?_
    ((i10PowerH1_intervalIntegrable n).const_mul
      (1 / ((2 : ℝ) ^ (n + 1) * ((n : ℝ) + 1) ^ 2)))
  intro x hx
  rw [Set.uIoc_of_le (by norm_num : (0 : ℝ) ≤ 1)] at hx
  have hxne : x ≠ 0 := ne_of_gt hx.1
  unfold i10CrossMoment
  rw [div_pow]
  field_simp [hxne]
  ring

private theorem i10CrossMoment_integral (n : ℕ) :
    (∫ x : ℝ in 0..1, i10CrossMoment n x) =
      i10HalfHarmonicTerm n := by
  calc
    (∫ x : ℝ in 0..1, i10CrossMoment n x) =
        ∫ x : ℝ in 0..1,
          (1 / ((2 : ℝ) ^ (n + 1) * ((n : ℝ) + 1) ^ 2)) *
            (x ^ n * H1 x) := by
      apply intervalIntegral.integral_congr
      intro x hx
      by_cases hxzero : x = 0
      · subst x
        simp [i10CrossMoment, H1]
      · unfold i10CrossMoment
        rw [div_pow]
        field_simp [hxzero]
        ring
    _ = i10HalfHarmonicTerm n := by
      rw [intervalIntegral.integral_const_mul, i10PowerH1_integral]
      unfold i10HalfHarmonicTerm
      rw [pow_succ]
      field_simp

private theorem i10CrossMoment_hasSum_pointwise
    {x : ℝ} (hx0 : 0 < x) (hx1 : x < 1) :
    HasSum (fun n : ℕ => i10CrossMoment n x)
      (dilog (x / 2) * H1 x / x) := by
  have hs : HasSum
      (fun n : ℕ => (x / 2) ^ (n + 1) / ((n : ℝ) + 1) ^ 2)
      (dilog (x / 2)) := by
    have habs : |x / 2| ≤ 1 := by
      rw [abs_of_pos (by positivity : 0 < x / 2)]
      linarith
    simpa only [dilog, Nat.cast_add, Nat.cast_one] using
      (dilog_summable habs).hasSum
  have hmul := hs.mul_right (H1 x / x)
  convert hmul using 1
  · ring

private theorem i10CrossMoment_integral_norm (n : ℕ) :
    (∫ x : ℝ in 0..1, ‖i10CrossMoment n x‖) =
      i10HalfHarmonicTerm n := by
  calc
    (∫ x : ℝ in 0..1, ‖i10CrossMoment n x‖) =
        ∫ x : ℝ in 0..1, i10CrossMoment n x := by
      apply intervalIntegral.integral_congr
      intro x hx
      rw [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)] at hx
      change |i10CrossMoment n x| = i10CrossMoment n x
      rw [abs_of_nonneg]
      unfold i10CrossMoment
      have hx0 : 0 ≤ x := hx.1
      have hx1 : x ≤ 1 := hx.2
      have hH1 : 0 ≤ H1 x := by
        unfold H1
        exact neg_nonneg.mpr (Real.log_nonpos
          (by linarith : 0 ≤ 1 - x) (by linarith : 1 - x ≤ 1))
      exact mul_nonneg (div_nonneg (pow_nonneg (by positivity) _)
        (sq_nonneg _)) (div_nonneg hH1 hx0)
    _ = i10HalfHarmonicTerm n := i10CrossMoment_integral n

private theorem i10CrossIntegral_hasSum :
    HasSum i10HalfHarmonicTerm
      (∫ x : ℝ in 0..1, dilog (x / 2) * H1 x / x) := by
  have hInt : ∀ n : ℕ, MeasureTheory.Integrable (i10CrossMoment n)
      (MeasureTheory.volume.restrict (Set.Ioc (0 : ℝ) 1)) :=
    fun n => (i10CrossMoment_intervalIntegrable n).1
  have hNorm : Summable (fun n : ℕ =>
      ∫ x : ℝ in Set.Ioc (0 : ℝ) 1, ‖i10CrossMoment n x‖) := by
    simpa only [← intervalIntegral.integral_of_le
      (by norm_num : (0 : ℝ) ≤ 1)] using
      (i10HalfHarmonicTerm_summable.congr fun n =>
        (i10CrossMoment_integral_norm n).symm)
  have h := MeasureTheory.hasSum_integral_of_summable_integral_norm
    (μ := MeasureTheory.volume.restrict (Set.Ioc (0 : ℝ) 1)) hInt hNorm
  have h' : HasSum i10HalfHarmonicTerm
      (∫ x : ℝ in Set.Ioc (0 : ℝ) 1,
        ∑' n : ℕ, i10CrossMoment n x) := by
    convert h using 1
    funext n
    rw [← intervalIntegral.integral_of_le (by norm_num : (0 : ℝ) ≤ 1)]
    exact (i10CrossMoment_integral n).symm
  convert h' using 1
  rw [intervalIntegral.integral_of_le (by norm_num : (0 : ℝ) ≤ 1)]
  apply MeasureTheory.setIntegral_congr_ae measurableSet_Ioc
  filter_upwards [MeasureTheory.Measure.ae_ne MeasureTheory.volume (1 : ℝ)]
    with x hxne hx
  have hxlt : x < 1 := lt_of_le_of_ne hx.2 hxne
  exact (i10CrossMoment_hasSum_pointwise hx.1 hxlt).tsum_eq.symm

private noncomputable def i10HalfMoment (n : ℕ) (x : ℝ) : ℝ :=
  (1 / 2 : ℝ) * harmonicNumber (n + 1) *
    x ^ n * Real.log (2 * x) ^ 2

private theorem i10HalfMoment_intervalIntegrable (n : ℕ) :
    IntervalIntegrable (i10HalfMoment n)
      MeasureTheory.volume 0 (1 / 2) := by
  let base : ℝ → ℝ := fun y => y ^ n * Real.log y ^ 2
  have hbase : IntervalIntegrable base MeasureTheory.volume 0 1 := by
    exact intervalIntegrable_logSq.continuousOn_mul (continuousOn_pow n)
  have hcomp := hbase.comp_mul_left (c := (2 : ℝ))
  have hcomp' : IntervalIntegrable (fun x : ℝ => base (2 * x))
      MeasureTheory.volume 0 (1 / 2) := by
    convert hcomp using 1 <;> norm_num
  apply IntervalIntegrable.congr
    (f := fun x : ℝ =>
      (harmonicNumber (n + 1) / (2 : ℝ) ^ (n + 1)) * base (2 * x)) ?_
    (hcomp'.const_mul
      (harmonicNumber (n + 1) / (2 : ℝ) ^ (n + 1)))
  intro x _
  unfold i10HalfMoment
  dsimp only [base]
  rw [mul_pow, pow_succ]
  field_simp

private theorem i10HalfMoment_integral (n : ℕ) :
    (∫ x : ℝ in 0..(1 / 2), i10HalfMoment n x) =
      i10HalfHarmonicTerm n := by
  let base : ℝ → ℝ := fun y => y ^ n * Real.log y ^ 2
  have hscale := intervalIntegral.smul_integral_comp_mul_left
    (a := (0 : ℝ)) (b := (1 / 2 : ℝ)) base (2 : ℝ)
  have hbase := integral_pow_mul_logSq n
  have hcomp : (∫ x : ℝ in 0..(1 / 2), base (2 * x)) =
      1 / ((n : ℝ) + 1) ^ 3 := by
    have hs : (2 : ℝ) * (∫ x : ℝ in 0..(1 / 2), base (2 * x)) =
        2 / ((n : ℝ) + 1) ^ 3 := by
      have hscale' : (2 : ℝ) * (∫ x : ℝ in 0..(1 / 2), base (2 * x)) =
          ∫ x : ℝ in 0..1, base x := by
        simpa [smul_eq_mul] using hscale
      have hbase' : (∫ x : ℝ in 0..1, base x) =
          2 / ((n : ℝ) + 1) ^ 3 := by
        simpa [base] using hbase
      exact hscale'.trans hbase'
    calc
      (∫ x : ℝ in 0..(1 / 2), base (2 * x)) =
          (1 / 2 : ℝ) *
            (2 * ∫ x : ℝ in 0..(1 / 2), base (2 * x)) := by ring
      _ = (1 / 2 : ℝ) * (2 / ((n : ℝ) + 1) ^ 3) := by rw [hs]
      _ = 1 / ((n : ℝ) + 1) ^ 3 := by ring
  calc
    (∫ x : ℝ in 0..(1 / 2), i10HalfMoment n x) =
        ∫ x : ℝ in 0..(1 / 2),
          (harmonicNumber (n + 1) / (2 : ℝ) ^ (n + 1)) *
            base (2 * x) := by
      apply intervalIntegral.integral_congr
      intro x _
      unfold i10HalfMoment
      dsimp only [base]
      rw [mul_pow]
      simp only [pow_succ]
      field_simp
    _ = i10HalfHarmonicTerm n := by
      rw [intervalIntegral.integral_const_mul, hcomp]
      unfold i10HalfHarmonicTerm
      rw [pow_succ]
      field_simp

private theorem i10HalfMoment_hasSum_pointwise
    {x : ℝ} (hx0 : 0 < x) (hxhalf : x < 1 / 2) :
    HasSum (fun n : ℕ => i10HalfMoment n x)
      (-(1 / 2 : ℝ) * Real.log (1 - x) * Real.log (2 * x) ^ 2 /
        (x * (1 - x))) := by
  have hs := harmonicNumber_generating_hasSum
    (x := x) (by rw [abs_of_pos hx0]; linarith) (ne_of_gt hx0)
  have hmul := hs.mul_left ((1 / 2 : ℝ) * Real.log (2 * x) ^ 2)
  convert hmul using 1
  · funext n
    unfold i10HalfMoment
    ring
  · field_simp [ne_of_gt hx0, ne_of_gt (by linarith : 0 < 1 - x)]

private theorem i10HalfMoment_integral_norm (n : ℕ) :
    (∫ x : ℝ in 0..(1 / 2), ‖i10HalfMoment n x‖) =
      i10HalfHarmonicTerm n := by
  calc
    (∫ x : ℝ in 0..(1 / 2), ‖i10HalfMoment n x‖) =
        ∫ x : ℝ in 0..(1 / 2), i10HalfMoment n x := by
      apply intervalIntegral.integral_congr
      intro x hx
      rw [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1 / 2)] at hx
      change |i10HalfMoment n x| = i10HalfMoment n x
      rw [abs_of_nonneg]
      unfold i10HalfMoment
      have hH : 0 ≤ harmonicNumber (n + 1) := harmonicNumber_nonneg (n + 1)
      have hxpow : 0 ≤ x ^ n := pow_nonneg hx.1 n
      have hlog : 0 ≤ Real.log (2 * x) ^ 2 := sq_nonneg _
      positivity
    _ = i10HalfHarmonicTerm n := i10HalfMoment_integral n

private theorem i10HalfIntegral_hasSum :
    HasSum i10HalfHarmonicTerm
      (∫ x : ℝ in 0..(1 / 2),
        -(1 / 2 : ℝ) * Real.log (1 - x) * Real.log (2 * x) ^ 2 /
          (x * (1 - x))) := by
  have hInt : ∀ n : ℕ, MeasureTheory.Integrable (i10HalfMoment n)
      (MeasureTheory.volume.restrict (Set.Ioc (0 : ℝ) (1 / 2))) :=
    fun n => (i10HalfMoment_intervalIntegrable n).1
  have hNorm : Summable (fun n : ℕ =>
      ∫ x : ℝ in Set.Ioc (0 : ℝ) (1 / 2), ‖i10HalfMoment n x‖) := by
    simpa only [← intervalIntegral.integral_of_le
      (by norm_num : (0 : ℝ) ≤ 1 / 2)] using
      (i10HalfHarmonicTerm_summable.congr fun n =>
        (i10HalfMoment_integral_norm n).symm)
  have h := MeasureTheory.hasSum_integral_of_summable_integral_norm
    (μ := MeasureTheory.volume.restrict (Set.Ioc (0 : ℝ) (1 / 2))) hInt hNorm
  have h' : HasSum i10HalfHarmonicTerm
      (∫ x : ℝ in Set.Ioc (0 : ℝ) (1 / 2),
        ∑' n : ℕ, i10HalfMoment n x) := by
    convert h using 1
    funext n
    rw [← intervalIntegral.integral_of_le
      (by norm_num : (0 : ℝ) ≤ 1 / 2)]
    exact (i10HalfMoment_integral n).symm
  convert h' using 1
  rw [intervalIntegral.integral_of_le
    (by norm_num : (0 : ℝ) ≤ 1 / 2)]
  apply MeasureTheory.setIntegral_congr_ae measurableSet_Ioc
  filter_upwards [MeasureTheory.Measure.ae_ne MeasureTheory.volume (1 / 2 : ℝ)]
    with x hxne hx
  have hxlt : x < 1 / 2 := lt_of_le_of_ne hx.2 hxne
  exact (i10HalfMoment_hasSum_pointwise hx.1 hxlt).tsum_eq.symm

private theorem i10Cross_eq_halfIntegral :
    (∫ x : ℝ in 0..1, dilog (x / 2) * H1 x / x) =
      ∫ x : ℝ in 0..(1 / 2),
        -(1 / 2 : ℝ) * Real.log (1 - x) * Real.log (2 * x) ^ 2 /
          (x * (1 - x)) := by
  exact i10CrossIntegral_hasSum.unique i10HalfIntegral_hasSum

private theorem i10LogPlusSquare_intervalIntegrable :
    IntervalIntegrable
      (fun x : ℝ => (Real.log x + Real.log 2) ^ 2)
      MeasureTheory.volume 0 (1 / 2) := by
  have hsq : IntervalIntegrable (fun x : ℝ => Real.log x ^ 2)
      MeasureTheory.volume 0 (1 / 2) := by
    apply intervalIntegrable_logSq.mono_set
    rw [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1 / 2),
      Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)]
    intro x hx
    exact ⟨hx.1, by linarith [hx.2]⟩
  have hlog : IntervalIntegrable Real.log MeasureTheory.volume 0 (1 / 2) :=
    intervalIntegral.intervalIntegrable_log'
  have hraw := (hsq.add (hlog.const_mul (2 * Real.log 2))).add
    (intervalIntegrable_const (c := Real.log 2 ^ 2))
  apply IntervalIntegrable.congr
    (f := fun x : ℝ =>
      Real.log x ^ 2 + 2 * Real.log 2 * Real.log x + Real.log 2 ^ 2) ?_ hraw
  intro x _
  ring

private theorem i10Ax_intervalIntegrable :
    IntervalIntegrable
      (fun x : ℝ =>
        Real.log (1 - x) * (Real.log x + Real.log 2) ^ 2 / x)
      MeasureTheory.volume 0 (1 / 2) := by
  let q : ℝ → ℝ := fun x => Real.log (1 - x) / x
  have hq : ContinuousOn q (Set.Ioo (0 : ℝ) (1 / 2)) := by
    intro x hx
    have hxlt : x < 1 / 2 := hx.2
    unfold q
    exact ((continuousAt_const.sub continuousAt_id).log
      (ne_of_gt (by linarith : 0 < 1 - x))).div
      continuousAt_id (ne_of_gt hx.1) |>.continuousWithinAt
  have hq0 : Tendsto q (𝓝[>] (0 : ℝ)) (𝓝 (-1)) := by
    have h := quadAltH1_div_self_tendsto_zero_right.neg
    convert h using 1 <;> simp only [q, H1] <;> ring
  have hqhalf : Tendsto q (𝓝[<] (1 / 2 : ℝ))
      (𝓝 (Real.log (1 - (1 / 2 : ℝ)) / (1 / 2))) := by
    have hc : ContinuousAt q (1 / 2) := by
      unfold q
      exact ((continuousAt_const.sub continuousAt_id).log
        (by norm_num : 1 - (1 / 2 : ℝ) ≠ 0)).div
          continuousAt_id (by norm_num)
    exact hc.tendsto.mono_left nhdsWithin_le_nhds
  have hmul := IntervalIntegrable.mul_of_continuousOn_Ioo_of_tendsto
    (by norm_num : (0 : ℝ) < 1 / 2)
    i10LogPlusSquare_intervalIntegrable hq hq0 hqhalf
  apply IntervalIntegrable.congr
    (f := fun x : ℝ => (Real.log x + Real.log 2) ^ 2 * q x) ?_ hmul
  intro x _
  unfold q
  ring

private noncomputable def i10AxPrimitive (x : ℝ) : ℝ :=
  -(Real.log x + Real.log 2) ^ 2 * dilog x +
    2 * (Real.log x + Real.log 2) * RamanujanChallenge.P26.trilog26 x -
    2 * polylog4 x

private theorem i10AxPrimitive_hasDerivAt
    {x : ℝ} (hx0 : 0 < x) (hxhalf : x < 1 / 2) :
    HasDerivAt i10AxPrimitive
      (Real.log (1 - x) * (Real.log x + Real.log 2) ^ 2 / x) x := by
  have hxne : x ≠ 0 := ne_of_gt hx0
  have hx1 : x < 1 := hxhalf.trans (by norm_num)
  have hlog := Real.hasDerivAt_log hxne
  have hd := dilog_hasDerivAt hx0 hx1
  have ht := RamanujanChallenge.P26.trilog26_hasDerivAt_of_abs_lt_one
    (by rw [abs_of_pos hx0]; exact hx1) hxne
  have hq := polylog4_hasDerivAt24_export
    (by rw [abs_of_pos hx0]; exact hx1) hxne
  unfold i10AxPrimitive
  have htotal :=
    (((hlog.add_const (Real.log 2)).pow 2).mul hd).neg.add
      (((hlog.add_const (Real.log 2)).mul ht).const_mul 2) |>.sub
      (hq.const_mul 2)
  convert htotal using 1
  · funext y
    simp only [Pi.add_apply, Pi.neg_apply, Pi.mul_apply, Pi.pow_apply,
      Pi.sub_apply]
    ring
  · field_simp [hxne]
    simp only [Pi.add_apply, Pi.pow_apply]
    ring

private theorem i10AxPrimitive_tendsto_zero :
    Tendsto i10AxPrimitive (𝓝[>] (0 : ℝ)) (𝓝 0) := by
  have hdSlope :
      Tendsto (fun x : ℝ => x⁻¹ * dilog x)
        (𝓝[>] (0 : ℝ)) (𝓝 1) := by
    simpa [dilog_zero] using
      RamanujanChallenge.P26.dilog_hasDerivAt_zero26.tendsto_slope_zero_right
  have htSlope :
      Tendsto (fun x : ℝ => x⁻¹ * RamanujanChallenge.P26.trilog26 x)
        (𝓝[>] (0 : ℝ)) (𝓝 1) := by
    simpa [RamanujanChallenge.P26.trilog26_zero] using
      trilog26_hasDerivAt_zero24_export.tendsto_slope_zero_right
  have hd : Tendsto dilog (𝓝[>] (0 : ℝ)) (𝓝 0) := by
    simpa [dilog_zero] using
      RamanujanChallenge.P26.dilog_hasDerivAt_zero26.continuousAt.tendsto.mono_left
        nhdsWithin_le_nhds
  have hq : Tendsto polylog4 (𝓝[>] (0 : ℝ)) (𝓝 0) := by
    simpa using
      polylog4_continuousAt_zero24_export.tendsto.mono_left nhdsWithin_le_nhds
  have hlogSqDilog :
      Tendsto (fun x : ℝ => Real.log x ^ 2 * dilog x)
        (𝓝[>] (0 : ℝ)) (𝓝 0) := by
    have hraw := logSq_mul_self_tendsto.mul hdSlope
    have hraw' : Tendsto (fun x : ℝ => Real.log x ^ 2 * dilog x)
        (𝓝[>] (0 : ℝ)) (𝓝 (0 * 1)) := by
      apply hraw.congr'
      filter_upwards [self_mem_nhdsWithin] with x hx
      have hx0 : 0 < x := hx
      field_simp [ne_of_gt hx0]
    simpa using hraw'
  have hlogDilog :
      Tendsto (fun x : ℝ => Real.log x * dilog x)
        (𝓝[>] (0 : ℝ)) (𝓝 0) := by
    have hraw := log_mul_self_tendsto.mul hdSlope
    have hraw' : Tendsto (fun x : ℝ => Real.log x * dilog x)
        (𝓝[>] (0 : ℝ)) (𝓝 (0 * 1)) := by
      apply hraw.congr'
      filter_upwards [self_mem_nhdsWithin] with x hx
      have hx0 : 0 < x := hx
      field_simp [ne_of_gt hx0]
    simpa using hraw'
  have hlogTrilog :
      Tendsto (fun x : ℝ => Real.log x * RamanujanChallenge.P26.trilog26 x)
        (𝓝[>] (0 : ℝ)) (𝓝 0) := by
    have hraw := log_mul_self_tendsto.mul htSlope
    have hraw' : Tendsto
        (fun x : ℝ => Real.log x * RamanujanChallenge.P26.trilog26 x)
        (𝓝[>] (0 : ℝ)) (𝓝 (0 * 1)) := by
      apply hraw.congr'
      filter_upwards [self_mem_nhdsWithin] with x hx
      have hx0 : 0 < x := hx
      field_simp [ne_of_gt hx0]
    simpa using hraw'
  have ht : Tendsto RamanujanChallenge.P26.trilog26
      (𝓝[>] (0 : ℝ)) (𝓝 0) := by
    simpa [RamanujanChallenge.P26.trilog26_zero] using
      trilog26_hasDerivAt_zero24_export.continuousAt.tendsto.mono_left
        nhdsWithin_le_nhds
  have hfirst :=
    ((hlogSqDilog.add (hlogDilog.const_mul (2 * Real.log 2))).add
      (hd.const_mul (Real.log 2 ^ 2))).neg
  have hsecond :=
    ((hlogTrilog.add (ht.const_mul (Real.log 2))).const_mul 2)
  have htotal := (hfirst.add hsecond).sub (hq.const_mul 2)
  simp only [neg_zero, add_zero, mul_zero, sub_zero] at htotal
  refine htotal.congr' ?_
  filter_upwards with x
  unfold i10AxPrimitive
  ring

private theorem i10AxPrimitive_tendsto_half :
    Tendsto i10AxPrimitive (𝓝[<] (1 / 2 : ℝ))
      (𝓝 (-2 * polylog4 (1 / 2))) := by
  have hd : ContinuousAt dilog (1 / 2 : ℝ) :=
    dilog_continuousOn_unit.continuousAt
      (Icc_mem_nhds (by norm_num : (-1 : ℝ) < 1 / 2)
        (by norm_num : (1 / 2 : ℝ) < 1))
  have ht : ContinuousAt RamanujanChallenge.P26.trilog26 (1 / 2 : ℝ) :=
    RamanujanChallenge.P26.trilog26_continuousOn_unit.continuousAt
      (Icc_mem_nhds (by norm_num : (-1 : ℝ) < 1 / 2)
        (by norm_num : (1 / 2 : ℝ) < 1))
  have hq : ContinuousAt polylog4 (1 / 2 : ℝ) :=
    polylog4_continuousOn_unit24_export.continuousAt
      (Icc_mem_nhds (by norm_num : (-1 : ℝ) < 1 / 2) (by norm_num))
  have hc : ContinuousAt i10AxPrimitive (1 / 2 : ℝ) := by
    unfold i10AxPrimitive
    have hlog : ContinuousAt Real.log (1 / 2 : ℝ) :=
      Real.continuousAt_log (by norm_num)
    have htotal := (((hlog.add_const (Real.log 2)).pow 2).mul hd).neg.add
      (((hlog.add_const (Real.log 2)).mul ht).const_mul 2) |>.sub
      (hq.const_mul 2)
    convert htotal using 1
    funext x
    simp only [Pi.add_apply, Pi.neg_apply, Pi.mul_apply, Pi.pow_apply,
      Pi.sub_apply]
    ring
  have hloghalf : Real.log (1 / 2 : ℝ) = -Real.log 2 := by
    rw [one_div, Real.log_inv]
  simpa [i10AxPrimitive, hloghalf] using
    hc.tendsto.mono_left nhdsWithin_le_nhds

private theorem i10Ax_integral :
    (∫ x : ℝ in 0..(1 / 2),
      Real.log (1 - x) * (Real.log x + Real.log 2) ^ 2 / x) =
      -2 * polylog4 (1 / 2) := by
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt_of_tendsto
    (f := i10AxPrimitive) (fa := (0 : ℝ))
    (fb := -2 * polylog4 (1 / 2)) (by norm_num)
    (fun x hx => i10AxPrimitive_hasDerivAt hx.1 hx.2)
    i10Ax_intervalIntegrable i10AxPrimitive_tendsto_zero
    i10AxPrimitive_tendsto_half]
  ring

private theorem i10A1_intervalIntegrable :
    IntervalIntegrable
      (fun x : ℝ =>
        Real.log (1 - x) * (Real.log x + Real.log 2) ^ 2 / (1 - x))
      MeasureTheory.volume 0 (1 / 2) := by
  let g : ℝ → ℝ := fun x => x / (1 - x)
  have hg : ContinuousOn g (Set.uIcc (0 : ℝ) (1 / 2)) := by
    rw [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1 / 2)]
    intro x hx
    have hxle : x ≤ 1 / 2 := hx.2
    unfold g
    exact continuousAt_id.div (continuousAt_const.sub continuousAt_id)
      (ne_of_gt (by linarith : 0 < 1 - x)) |>.continuousWithinAt
  have hprod := i10Ax_intervalIntegrable.mul_continuousOn hg
  apply IntervalIntegrable.congr
    (f := fun x : ℝ =>
      (Real.log (1 - x) * (Real.log x + Real.log 2) ^ 2 / x) * g x) ?_
    hprod
  intro x hx
  rw [Set.uIoc_of_le (by norm_num : (0 : ℝ) ≤ 1 / 2)] at hx
  have hxle : x ≤ 1 / 2 := hx.2
  by_cases hxzero : x = 0
  · subst x
    simp [g]
  · have hden : 1 - x ≠ 0 := ne_of_gt (by linarith : 0 < 1 - x)
    unfold g
    field_simp [hxzero, hden]

private noncomputable def i10A1Primitive (x : ℝ) : ℝ :=
  -(1 / 2 : ℝ) * Real.log (1 - x) ^ 2 *
    (Real.log x + Real.log 2) ^ 2

private theorem i10A1Primitive_hasDerivAt
    {x : ℝ} (hx0 : 0 < x) (hxhalf : x < 1 / 2) :
    HasDerivAt i10A1Primitive
      (Real.log (1 - x) * (Real.log x + Real.log 2) ^ 2 / (1 - x) -
        Real.log x * Real.log (1 - x) ^ 2 / x -
        Real.log 2 * (Real.log (1 - x) ^ 2 / x)) x := by
  have hxne : x ≠ 0 := ne_of_gt hx0
  have hden : 1 - x ≠ 0 := ne_of_gt (by linarith : 0 < 1 - x)
  have hinner : HasDerivAt (fun y : ℝ => 1 - y) (-1) x := by
    convert (hasDerivAt_const x 1).sub (hasDerivAt_id x) using 1 <;> norm_num
  have ha : HasDerivAt (fun y : ℝ => Real.log (1 - y))
      (-1 / (1 - x)) x := by
    convert hinner.log hden using 1
  have hl := Real.hasDerivAt_log hxne
  unfold i10A1Primitive
  have htotal := ((ha.pow 2).mul ((hl.add_const (Real.log 2)).pow 2)).const_mul
    (-(1 / 2 : ℝ))
  convert htotal using 1
  · funext y
    simp only [Pi.mul_apply, Pi.pow_apply, Pi.neg_apply]
    ring
  · simp only [Pi.add_apply, Pi.pow_apply]
    field_simp [hxne, hden]
    ring

private theorem i10A1Primitive_tendsto_zero :
    Tendsto i10A1Primitive (𝓝[>] (0 : ℝ)) (𝓝 0) := by
  have hH := quadAltH1_div_self_tendsto_zero_right
  have hid : Tendsto (fun x : ℝ => x) (𝓝[>] (0 : ℝ)) (𝓝 0) :=
    tendsto_id.mono_left nhdsWithin_le_nhds
  have hv := log_mul_self_tendsto.add (hid.const_mul (Real.log 2))
  have hraw := ((hH.pow 2).mul (hv.pow 2)).const_mul (-(1 / 2 : ℝ))
  have hraw' : Tendsto
      (fun x : ℝ => -((1 / 2 : ℝ) * ((H1 x / x) ^ 2 *
        (Real.log x * x + Real.log 2 * x) ^ 2)))
      (𝓝[>] (0 : ℝ)) (𝓝 0) := by
    simpa using hraw
  refine hraw'.congr' ?_
  filter_upwards [self_mem_nhdsWithin] with x hx
  have hx0 : 0 < x := hx
  unfold i10A1Primitive H1
  field_simp [ne_of_gt hx0]

private theorem i10A1Primitive_tendsto_half :
    Tendsto i10A1Primitive (𝓝[<] (1 / 2 : ℝ)) (𝓝 0) := by
  have hlog : ContinuousAt Real.log (1 / 2 : ℝ) :=
    Real.continuousAt_log (by norm_num)
  have hlogSub : ContinuousAt (fun x : ℝ => Real.log (1 - x)) (1 / 2) :=
    (continuousAt_const.sub continuousAt_id).log (by norm_num)
  have hc : ContinuousAt i10A1Primitive (1 / 2 : ℝ) := by
    unfold i10A1Primitive
    have htotal := (((hlogSub.pow 2).mul
      ((hlog.add_const (Real.log 2)).pow 2)).const_mul (-(1 / 2 : ℝ)))
    convert htotal using 1
    funext x
    simp only [Pi.mul_apply, Pi.pow_apply]
    ring
  have hloghalf : Real.log (1 / 2 : ℝ) = -Real.log 2 := by
    rw [one_div, Real.log_inv]
  simpa [i10A1Primitive, hloghalf] using
    hc.tendsto.mono_left nhdsWithin_le_nhds

private theorem i10A1_integral :
    (∫ x : ℝ in 0..(1 / 2),
      Real.log (1 - x) * (Real.log x + Real.log 2) ^ 2 / (1 - x)) =
      -(1 / 12 : ℝ) * Real.log 2 ^ 4 +
        (1 / 4 : ℝ) * Real.log 2 * zeta3_24 -
        (1 / 4 : ℝ) * (Real.pi ^ 4 / 90) := by
  let f : ℝ → ℝ := fun x =>
    Real.log (1 - x) * (Real.log x + Real.log 2) ^ 2 / (1 - x)
  let j : ℝ → ℝ := fun x => Real.log x * Real.log (1 - x) ^ 2 / x
  let b : ℝ → ℝ := fun x => Real.log (1 - x) ^ 2 / x
  have hf : IntervalIntegrable f MeasureTheory.volume 0 (1 / 2) := by
    simpa [f] using i10A1_intervalIntegrable
  have hj : IntervalIntegrable j MeasureTheory.volume 0 (1 / 2) := by
    simpa [j] using quadAltHalfQuarticCoreIntervalIntegrable24
  have hb : IntervalIntegrable b MeasureTheory.volume 0 (1 / 2) := by
    simpa [b] using quadAltHalfLogSquareIntervalIntegrable24
  have hder : IntervalIntegrable
      (fun x : ℝ => f x - j x - Real.log 2 * b x)
      MeasureTheory.volume 0 (1 / 2) :=
    (hf.sub hj).sub (hb.const_mul (Real.log 2))
  have hzero : (∫ x : ℝ in 0..(1 / 2),
      (f x - j x - Real.log 2 * b x)) = 0 := by
    rw [intervalIntegral.integral_eq_sub_of_hasDerivAt_of_tendsto
      (f := i10A1Primitive) (fa := (0 : ℝ)) (fb := (0 : ℝ))
      (by norm_num)
      (fun x hx => by
        simpa [f, j, b] using i10A1Primitive_hasDerivAt hx.1 hx.2)
      hder i10A1Primitive_tendsto_zero i10A1Primitive_tendsto_half]
    ring
  rw [intervalIntegral.integral_sub (hf.sub hj) (hb.const_mul (Real.log 2)),
    intervalIntegral.integral_sub hf hj,
    intervalIntegral.integral_const_mul] at hzero
  have hjval : (∫ x : ℝ in 0..(1 / 2), j x) =
      (1 / 4 : ℝ) * Real.log 2 ^ 4 -
        (1 / 4 : ℝ) * (Real.pi ^ 4 / 90) := by
    simpa [j] using quadAltHalfQuarticCoreIntegral24
  have hbval : (∫ x : ℝ in 0..(1 / 2), b x) =
      (1 / 4 : ℝ) * zeta3_24 -
        (1 / 3 : ℝ) * Real.log 2 ^ 3 := by
    simpa [b] using quadAltHalfLogSquareIntegral24
  rw [hjval, hbval] at hzero
  dsimp only [f] at hzero
  linarith

private theorem i10HalfIntegral_eq :
    (∫ x : ℝ in 0..(1 / 2),
      -(1 / 2 : ℝ) * Real.log (1 - x) * Real.log (2 * x) ^ 2 /
        (x * (1 - x))) =
      polylog4 (1 / 2) + (1 / 24 : ℝ) * Real.log 2 ^ 4 -
        (1 / 8 : ℝ) * Real.log 2 * zeta3_24 +
        (1 / 8 : ℝ) * (Real.pi ^ 4 / 90) := by
  calc
    (∫ x : ℝ in 0..(1 / 2),
      -(1 / 2 : ℝ) * Real.log (1 - x) * Real.log (2 * x) ^ 2 /
        (x * (1 - x))) =
        ∫ x : ℝ in 0..(1 / 2), -(1 / 2 : ℝ) *
          (Real.log (1 - x) * (Real.log x + Real.log 2) ^ 2 / x +
           Real.log (1 - x) * (Real.log x + Real.log 2) ^ 2 / (1 - x)) := by
      apply intervalIntegral.integral_congr
      intro x hx
      rw [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1 / 2)] at hx
      have hxle : x ≤ 1 / 2 := hx.2
      by_cases hxzero : x = 0
      · subst x
        simp
      · have hden : 1 - x ≠ 0 := ne_of_gt (by linarith : 0 < 1 - x)
        dsimp only
        rw [Real.log_mul (by norm_num : (2 : ℝ) ≠ 0) hxzero]
        field_simp [hxzero, hden]
        ring
    _ = polylog4 (1 / 2) + (1 / 24 : ℝ) * Real.log 2 ^ 4 -
        (1 / 8 : ℝ) * Real.log 2 * zeta3_24 +
        (1 / 8 : ℝ) * (Real.pi ^ 4 / 90) := by
      rw [intervalIntegral.integral_const_mul,
        intervalIntegral.integral_add i10Ax_intervalIntegrable i10A1_intervalIntegrable,
        i10Ax_integral, i10A1_integral]
      ring

private noncomputable def i10Q (x : ℝ) : ℝ := H1 x / x

private theorem i10Q_continuousOn_half :
    ContinuousOn i10Q (Set.Ioo (0 : ℝ) (1 / 2)) := by
  intro x hx
  have hxlt : x < 1 / 2 := hx.2
  unfold i10Q H1
  exact ((continuousAt_const.sub continuousAt_id).log
    (ne_of_gt (by linarith : 0 < 1 - x))).neg.div
      continuousAt_id (ne_of_gt hx.1) |>.continuousWithinAt

private theorem i10Q_tendsto_zero :
    Tendsto i10Q (𝓝[>] (0 : ℝ)) (𝓝 1) := by
  exact quadAltH1_div_self_tendsto_zero_right

private theorem i10Q_tendsto_half :
    Tendsto i10Q (𝓝[<] (1 / 2 : ℝ))
      (𝓝 (H1 (1 / 2) / (1 / 2))) := by
  have hc : ContinuousAt i10Q (1 / 2 : ℝ) := by
    unfold i10Q H1
    exact ((continuousAt_const.sub continuousAt_id).log (by norm_num)).neg.div
      continuousAt_id (by norm_num)
  exact hc.tendsto.mono_left nhdsWithin_le_nhds

private theorem i10H1_intervalIntegrable :
    IntervalIntegrable H1 MeasureTheory.volume (1 / 2) 1 := by
  have hlog := (intervalIntegral.intervalIntegrable_log'
    (a := (0 : ℝ)) (b := (1 / 2))).comp_sub_left 1
  convert hlog.symm.neg using 1 <;> norm_num [H1]

private theorem i10H1Div_intervalIntegrable :
    IntervalIntegrable (fun x : ℝ => H1 x / x)
      MeasureTheory.volume 0 1 := by
  have hleft : IntervalIntegrable i10Q MeasureTheory.volume 0 (1 / 2) :=
    intervalIntegrable_of_continuousOn_Ioo_of_tendsto (by norm_num)
      i10Q_continuousOn_half i10Q_tendsto_zero i10Q_tendsto_half
  have hg : ContinuousOn (fun x : ℝ => 1 / x)
      (Set.uIcc (1 / 2 : ℝ) 1) := by
    rw [Set.uIcc_of_le (by norm_num : (1 / 2 : ℝ) ≤ 1)]
    intro x hx
    have hxlo : 1 / 2 ≤ x := hx.1
    exact continuousAt_const.div continuousAt_id
      (ne_of_gt (by linarith : 0 < x)) |>.continuousWithinAt
  have hright := i10H1_intervalIntegrable.continuousOn_mul hg
  apply IntervalIntegrable.trans (b := (1 / 2 : ℝ))
  · simpa [i10Q] using hleft
  · apply IntervalIntegrable.congr
      (f := fun x : ℝ => (1 / x) * H1 x) ?_ hright
    intro x _
    ring

private theorem i10LogH1Div_intervalIntegrable :
    IntervalIntegrable (fun x : ℝ => Real.log x * H1 x / x)
      MeasureTheory.volume 0 1 := by
  have hlog : IntervalIntegrable Real.log MeasureTheory.volume 0 (1 / 2) :=
    intervalIntegral.intervalIntegrable_log'
  have hleft := IntervalIntegrable.mul_of_continuousOn_Ioo_of_tendsto
    (by norm_num : (0 : ℝ) < 1 / 2) hlog i10Q_continuousOn_half
      i10Q_tendsto_zero i10Q_tendsto_half
  have hg : ContinuousOn (fun x : ℝ => Real.log x / x)
      (Set.uIcc (1 / 2 : ℝ) 1) := by
    rw [Set.uIcc_of_le (by norm_num : (1 / 2 : ℝ) ≤ 1)]
    intro x hx
    have hxlo : 1 / 2 ≤ x := hx.1
    exact (Real.continuousAt_log (ne_of_gt (by linarith : 0 < x))).div
      continuousAt_id (ne_of_gt (by linarith : 0 < x)) |>.continuousWithinAt
  have hright := i10H1_intervalIntegrable.continuousOn_mul hg
  apply IntervalIntegrable.trans (b := (1 / 2 : ℝ))
  · apply IntervalIntegrable.congr
      (f := fun x : ℝ => Real.log x * i10Q x) ?_ hleft
    intro x _
    unfold i10Q
    ring
  · apply IntervalIntegrable.congr
      (f := fun x : ℝ => (Real.log x / x) * H1 x) ?_ hright
    intro x _
    ring

private theorem i10LogSqH1Div_intervalIntegrable :
    IntervalIntegrable (fun x : ℝ => Real.log x ^ 2 * H1 x / x)
      MeasureTheory.volume 0 1 := by
  have hsq : IntervalIntegrable (fun x : ℝ => Real.log x ^ 2)
      MeasureTheory.volume 0 (1 / 2) := by
    apply intervalIntegrable_logSq.mono_set
    rw [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1 / 2),
      Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)]
    intro x hx
    exact ⟨hx.1, by linarith [hx.2]⟩
  have hleft := IntervalIntegrable.mul_of_continuousOn_Ioo_of_tendsto
    (by norm_num : (0 : ℝ) < 1 / 2) hsq i10Q_continuousOn_half
      i10Q_tendsto_zero i10Q_tendsto_half
  have hg : ContinuousOn (fun x : ℝ => Real.log x ^ 2 / x)
      (Set.uIcc (1 / 2 : ℝ) 1) := by
    rw [Set.uIcc_of_le (by norm_num : (1 / 2 : ℝ) ≤ 1)]
    intro x hx
    have hxlo : 1 / 2 ≤ x := hx.1
    exact ((Real.continuousAt_log (ne_of_gt (by linarith : 0 < x))).pow 2).div
      continuousAt_id (ne_of_gt (by linarith : 0 < x)) |>.continuousWithinAt
  have hright := i10H1_intervalIntegrable.continuousOn_mul hg
  apply IntervalIntegrable.trans (b := (1 / 2 : ℝ))
  · apply IntervalIntegrable.congr
      (f := fun x : ℝ => Real.log x ^ 2 * i10Q x) ?_ hleft
    intro x _
    unfold i10Q
    ring
  · apply IntervalIntegrable.congr
      (f := fun x : ℝ => (Real.log x ^ 2 / x) * H1 x) ?_ hright
    intro x _
    ring

private theorem i10_dilog_tendsto_one :
    Tendsto dilog (𝓝[<] (1 : ℝ)) (𝓝 (Real.pi ^ 2 / 6)) := by
  have hWithin : ContinuousWithinAt dilog (Iio (1 : ℝ)) 1 :=
    (dilog_continuousOn_unit 1 (by norm_num)).mono_of_mem_nhdsWithin
      (Icc_mem_nhdsLT (show (-1 : ℝ) < 1 by norm_num))
  simpa [dilog_one] using hWithin.tendsto

private theorem i10_trilog_tendsto_one :
    Tendsto RamanujanChallenge.P26.trilog26
      (𝓝[<] (1 : ℝ)) (𝓝 zeta3_24) := by
  have hWithin : ContinuousWithinAt RamanujanChallenge.P26.trilog26
      (Iio (1 : ℝ)) 1 :=
    (RamanujanChallenge.P26.trilog26_continuousOn_unit
      1 (by norm_num)).mono_of_mem_nhdsWithin
      (Icc_mem_nhdsLT (show (-1 : ℝ) < 1 by norm_num))
  simpa [RamanujanChallenge.P26.trilog26_one,
    RamanujanChallenge.P26.zeta3, zeta3_24] using hWithin.tendsto

private theorem i10_polylog4_tendsto_one :
    Tendsto polylog4 (𝓝[<] (1 : ℝ)) (𝓝 (Real.pi ^ 4 / 90)) := by
  have hWithin : ContinuousWithinAt polylog4 (Iio (1 : ℝ)) 1 :=
    (polylog4_continuousOn_unit24_export
      1 (by norm_num)).mono_of_mem_nhdsWithin
      (Icc_mem_nhdsLT (show (-1 : ℝ) < 1 by norm_num))
  simpa [polylog4_one24_export] using hWithin.tendsto

private theorem i10_log_dilog_tendsto_zero :
    Tendsto (fun x : ℝ => Real.log x * dilog x)
      (𝓝[>] (0 : ℝ)) (𝓝 0) := by
  have hslope : Tendsto (fun x : ℝ => x⁻¹ * dilog x)
      (𝓝[>] (0 : ℝ)) (𝓝 1) := by
    simpa [dilog_zero] using
      RamanujanChallenge.P26.dilog_hasDerivAt_zero26.tendsto_slope_zero_right
  have hraw := log_mul_self_tendsto.mul hslope
  have hraw' : Tendsto (fun x : ℝ => Real.log x * dilog x)
      (𝓝[>] (0 : ℝ)) (𝓝 (0 * 1)) := by
    apply hraw.congr'
    filter_upwards [self_mem_nhdsWithin] with x hx
    have hx0 : 0 < x := hx
    field_simp [ne_of_gt hx0]
  simpa using hraw'

private theorem i10_logSq_dilog_tendsto_zero :
    Tendsto (fun x : ℝ => Real.log x ^ 2 * dilog x)
      (𝓝[>] (0 : ℝ)) (𝓝 0) := by
  have hslope : Tendsto (fun x : ℝ => x⁻¹ * dilog x)
      (𝓝[>] (0 : ℝ)) (𝓝 1) := by
    simpa [dilog_zero] using
      RamanujanChallenge.P26.dilog_hasDerivAt_zero26.tendsto_slope_zero_right
  have hraw := logSq_mul_self_tendsto.mul hslope
  have hraw' : Tendsto (fun x : ℝ => Real.log x ^ 2 * dilog x)
      (𝓝[>] (0 : ℝ)) (𝓝 (0 * 1)) := by
    apply hraw.congr'
    filter_upwards [self_mem_nhdsWithin] with x hx
    have hx0 : 0 < x := hx
    field_simp [ne_of_gt hx0]
  simpa using hraw'

private theorem i10_log_trilog_tendsto_zero :
    Tendsto (fun x : ℝ =>
      Real.log x * RamanujanChallenge.P26.trilog26 x)
      (𝓝[>] (0 : ℝ)) (𝓝 0) := by
  have hslope : Tendsto
      (fun x : ℝ => x⁻¹ * RamanujanChallenge.P26.trilog26 x)
      (𝓝[>] (0 : ℝ)) (𝓝 1) := by
    simpa [RamanujanChallenge.P26.trilog26_zero] using
      trilog26_hasDerivAt_zero24_export.tendsto_slope_zero_right
  have hraw := log_mul_self_tendsto.mul hslope
  have hraw' : Tendsto
      (fun x : ℝ => Real.log x * RamanujanChallenge.P26.trilog26 x)
      (𝓝[>] (0 : ℝ)) (𝓝 (0 * 1)) := by
    apply hraw.congr'
    filter_upwards [self_mem_nhdsWithin] with x hx
    have hx0 : 0 < x := hx
    field_simp [ne_of_gt hx0]
  simpa using hraw'

private theorem i10H1Div_integral :
    (∫ x : ℝ in 0..1, H1 x / x) = Real.pi ^ 2 / 6 := by
  have hzero : Tendsto dilog (𝓝[>] (0 : ℝ)) (𝓝 0) := by
    simpa [dilog_zero] using
      RamanujanChallenge.P26.dilog_hasDerivAt_zero26.continuousAt.tendsto.mono_left
        nhdsWithin_le_nhds
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt_of_tendsto
    (f := dilog) (fa := (0 : ℝ)) (fb := Real.pi ^ 2 / 6)
    (by norm_num) (fun x hx => by
      simpa [H1] using dilog_hasDerivAt hx.1 hx.2)
    i10H1Div_intervalIntegrable hzero i10_dilog_tendsto_one]
  ring

private noncomputable def i10LogDilogPrimitive (x : ℝ) : ℝ :=
  Real.log x * dilog x - RamanujanChallenge.P26.trilog26 x

private theorem i10LogDilogPrimitive_hasDerivAt
    {x : ℝ} (hx0 : 0 < x) (hx1 : x < 1) :
    HasDerivAt i10LogDilogPrimitive (Real.log x * H1 x / x) x := by
  have hxne : x ≠ 0 := ne_of_gt hx0
  have hlog := Real.hasDerivAt_log hxne
  have hd := dilog_hasDerivAt hx0 hx1
  have ht := RamanujanChallenge.P26.trilog26_hasDerivAt_of_abs_lt_one
    (by rw [abs_of_pos hx0]; exact hx1) hxne
  unfold i10LogDilogPrimitive
  convert (hlog.mul hd).sub ht using 1
  unfold H1
  field_simp [hxne]
  ring

private theorem i10LogDilogPrimitive_tendsto_zero :
    Tendsto i10LogDilogPrimitive (𝓝[>] (0 : ℝ)) (𝓝 0) := by
  have ht : Tendsto RamanujanChallenge.P26.trilog26
      (𝓝[>] (0 : ℝ)) (𝓝 0) := by
    simpa [RamanujanChallenge.P26.trilog26_zero] using
      trilog26_hasDerivAt_zero24_export.continuousAt.tendsto.mono_left
        nhdsWithin_le_nhds
  simpa [i10LogDilogPrimitive] using i10_log_dilog_tendsto_zero.sub ht

private theorem i10LogDilogPrimitive_tendsto_one :
    Tendsto i10LogDilogPrimitive (𝓝[<] (1 : ℝ)) (𝓝 (-zeta3_24)) := by
  have hlog : Tendsto Real.log (𝓝[<] (1 : ℝ)) (𝓝 0) := by
    simpa using (Real.continuousAt_log (by norm_num : (1 : ℝ) ≠ 0)).tendsto.mono_left
      nhdsWithin_le_nhds
  have htotal := (hlog.mul i10_dilog_tendsto_one).sub i10_trilog_tendsto_one
  simpa [i10LogDilogPrimitive] using htotal

private theorem i10LogH1Div_integral :
    (∫ x : ℝ in 0..1, Real.log x * H1 x / x) = -zeta3_24 := by
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt_of_tendsto
    (f := i10LogDilogPrimitive) (fa := (0 : ℝ)) (fb := -zeta3_24)
    (by norm_num) (fun x hx => i10LogDilogPrimitive_hasDerivAt hx.1 hx.2)
    i10LogH1Div_intervalIntegrable i10LogDilogPrimitive_tendsto_zero
    i10LogDilogPrimitive_tendsto_one]
  ring

private noncomputable def i10LogSqDilogPrimitive (x : ℝ) : ℝ :=
  Real.log x ^ 2 * dilog x -
    2 * Real.log x * RamanujanChallenge.P26.trilog26 x +
    2 * polylog4 x

private theorem i10LogSqDilogPrimitive_hasDerivAt
    {x : ℝ} (hx0 : 0 < x) (hx1 : x < 1) :
    HasDerivAt i10LogSqDilogPrimitive (Real.log x ^ 2 * H1 x / x) x := by
  have hxne : x ≠ 0 := ne_of_gt hx0
  have hlog := Real.hasDerivAt_log hxne
  have hd := dilog_hasDerivAt hx0 hx1
  have ht := RamanujanChallenge.P26.trilog26_hasDerivAt_of_abs_lt_one
    (by rw [abs_of_pos hx0]; exact hx1) hxne
  have hq := polylog4_hasDerivAt24_export
    (by rw [abs_of_pos hx0]; exact hx1) hxne
  unfold i10LogSqDilogPrimitive
  have htotal := (((hlog.pow 2).mul hd).sub
    ((hlog.mul ht).const_mul 2)).add (hq.const_mul 2)
  convert htotal using 1
  · funext y
    simp only [Pi.add_apply, Pi.sub_apply, Pi.mul_apply, Pi.pow_apply]
    ring
  · field_simp [hxne]
    unfold H1
    simp only [Pi.pow_apply]
    ring

private theorem i10LogSqDilogPrimitive_tendsto_zero :
    Tendsto i10LogSqDilogPrimitive (𝓝[>] (0 : ℝ)) (𝓝 0) := by
  have hq : Tendsto polylog4 (𝓝[>] (0 : ℝ)) (𝓝 0) := by
    simpa using polylog4_continuousAt_zero24_export.tendsto.mono_left
      nhdsWithin_le_nhds
  have htotal := (i10_logSq_dilog_tendsto_zero.sub
    (i10_log_trilog_tendsto_zero.const_mul 2)).add (hq.const_mul 2)
  have htotal' : Tendsto
      (fun x : ℝ => Real.log x ^ 2 * dilog x -
        2 * (Real.log x * RamanujanChallenge.P26.trilog26 x) +
        2 * polylog4 x) (𝓝[>] (0 : ℝ)) (𝓝 0) := by
    simpa using htotal
  refine htotal'.congr' ?_
  filter_upwards with x
  unfold i10LogSqDilogPrimitive
  ring

private theorem i10LogSqDilogPrimitive_tendsto_one :
    Tendsto i10LogSqDilogPrimitive (𝓝[<] (1 : ℝ))
      (𝓝 (2 * (Real.pi ^ 4 / 90))) := by
  have hlog : Tendsto Real.log (𝓝[<] (1 : ℝ)) (𝓝 0) := by
    simpa using (Real.continuousAt_log (by norm_num : (1 : ℝ) ≠ 0)).tendsto.mono_left
      nhdsWithin_le_nhds
  have htotal := (((hlog.pow 2).mul i10_dilog_tendsto_one).sub
    ((hlog.mul i10_trilog_tendsto_one).const_mul 2)).add
      (i10_polylog4_tendsto_one.const_mul 2)
  have htotal' : Tendsto
      (fun x : ℝ => Real.log x ^ 2 * dilog x -
        2 * (Real.log x * RamanujanChallenge.P26.trilog26 x) +
        2 * polylog4 x) (𝓝[<] (1 : ℝ))
        (𝓝 (2 * (Real.pi ^ 4 / 90))) := by
    simpa using htotal
  refine htotal'.congr' ?_
  filter_upwards with x
  unfold i10LogSqDilogPrimitive
  ring

private theorem i10LogSqH1Div_integral :
    (∫ x : ℝ in 0..1, Real.log x ^ 2 * H1 x / x) =
      2 * (Real.pi ^ 4 / 90) := by
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt_of_tendsto
    (f := i10LogSqDilogPrimitive) (fa := (0 : ℝ))
    (fb := 2 * (Real.pi ^ 4 / 90)) (by norm_num)
    (fun x hx => i10LogSqDilogPrimitive_hasDerivAt hx.1 hx.2)
    i10LogSqH1Div_intervalIntegrable i10LogSqDilogPrimitive_tendsto_zero
    i10LogSqDilogPrimitive_tendsto_one]
  ring

private theorem i10Cross_intervalIntegrable :
    IntervalIntegrable (fun x : ℝ => dilog (x / 2) * H1 x / x)
      MeasureTheory.volume 0 1 := by
  have hg : ContinuousOn (fun x : ℝ => dilog (x / 2)) (Set.uIcc (0 : ℝ) 1) := by
    rw [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)]
    apply dilog_continuousOn_unit.comp (by fun_prop)
    intro x hx
    constructor <;> linarith [hx.1, hx.2]
  have hprod := i10H1Div_intervalIntegrable.continuousOn_mul hg
  apply IntervalIntegrable.congr
    (f := fun x : ℝ => dilog (x / 2) * (H1 x / x)) ?_ hprod
  intro x _
  ring

private theorem i10Cross_integral :
    (∫ x : ℝ in 0..1, dilog (x / 2) * H1 x / x) =
      polylog4 (1 / 2) + (1 / 24 : ℝ) * Real.log 2 ^ 4 -
        (1 / 8 : ℝ) * Real.log 2 * zeta3_24 +
        (1 / 8 : ℝ) * (Real.pi ^ 4 / 90) := by
  exact i10Cross_eq_halfIntegral.trans i10HalfIntegral_eq

private theorem quadAltI10_raw :
    I10 = (Real.pi ^ 2 / 6) ^ 2 -
      2 * polylog4 (1 / 2) - (1 / 12 : ℝ) * Real.log 2 ^ 4 -
      (7 / 4 : ℝ) * Real.log 2 * zeta3_24 -
      Real.log 2 ^ 2 * (Real.pi ^ 2 / 6) -
      (9 / 4 : ℝ) * (Real.pi ^ 4 / 90) := by
  let q : ℝ → ℝ := fun x => H1 x / x
  let c : ℝ → ℝ := fun x => dilog (x / 2) * H1 x / x
  let a1 : ℝ → ℝ := fun x => Real.log x * H1 x / x
  let a2 : ℝ → ℝ := fun x => Real.log x ^ 2 * H1 x / x
  let e : ℝ → ℝ := fun x =>
    ((Real.pi ^ 2 / 6) * q x - 2 * c x) -
      (a2 x - 2 * Real.log 2 * a1 x + Real.log 2 ^ 2 * q x)
  have hq : IntervalIntegrable q MeasureTheory.volume 0 1 := by
    simpa [q] using i10H1Div_intervalIntegrable
  have hc : IntervalIntegrable c MeasureTheory.volume 0 1 := by
    simpa [c] using i10Cross_intervalIntegrable
  have ha1 : IntervalIntegrable a1 MeasureTheory.volume 0 1 := by
    simpa [a1] using i10LogH1Div_intervalIntegrable
  have ha2 : IntervalIntegrable a2 MeasureTheory.volume 0 1 := by
    simpa [a2] using i10LogSqH1Div_intervalIntegrable
  have hfirst := (hq.const_mul (Real.pi ^ 2 / 6)).sub (hc.const_mul 2)
  have hsecond := (ha2.sub (ha1.const_mul (2 * Real.log 2))).add
    (hq.const_mul (Real.log 2 ^ 2))
  have he : IntervalIntegrable e MeasureTheory.volume 0 1 := by
    exact hfirst.sub hsecond
  calc
    I10 = ∫ x : ℝ in 0..1, e x := by
      unfold I10
      apply intervalIntegral.integral_congr
      intro x hx
      rw [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)] at hx
      by_cases hxzero : x = 0
      · subst x
        simp [e, q, c, a1, a2, W0, H1]
      · dsimp only [e, q, c, a1, a2]
        unfold W0
        rw [Real.log_div hxzero (by norm_num : (2 : ℝ) ≠ 0)]
        ring
    _ = (Real.pi ^ 2 / 6) ^ 2 -
      2 * polylog4 (1 / 2) - (1 / 12 : ℝ) * Real.log 2 ^ 4 -
      (7 / 4 : ℝ) * Real.log 2 * zeta3_24 -
      Real.log 2 ^ 2 * (Real.pi ^ 2 / 6) -
      (9 / 4 : ℝ) * (Real.pi ^ 4 / 90) := by
      dsimp only [e]
      rw [intervalIntegral.integral_sub hfirst hsecond,
        intervalIntegral.integral_sub (hq.const_mul (Real.pi ^ 2 / 6))
          (hc.const_mul 2),
        intervalIntegral.integral_add
          (ha2.sub (ha1.const_mul (2 * Real.log 2)))
          (hq.const_mul (Real.log 2 ^ 2)),
        intervalIntegral.integral_sub ha2 (ha1.const_mul (2 * Real.log 2))]
      repeat rw [intervalIntegral.integral_const_mul]
      dsimp only [q, c, a1, a2]
      rw [i10H1Div_integral, i10Cross_integral, i10LogH1Div_integral,
        i10LogSqH1Div_integral]
      ring

theorem quadAltI10_eq :
    I10 = -(1 / 2) * quadAltK
      - (3 / 2) * Real.log 2 ^ 2 * (Real.pi ^ 2 / 6)
      - (13 / 20) * (Real.pi ^ 2 / 6) ^ 2 := by
  rw [quadAltI10_raw, quadAltK_eq]
  unfold alternatingCubicLinearEulerValue24 cubicLinearEulerValue24
  ring


private def testAltZetaTwoTerm (n : ℕ) : ℝ :=
  (-1 : ℝ) ^ (n + 1) / ((n : ℝ) + 1) ^ 2

private theorem testAltZetaTwoTerm_hasSum :
    HasSum testAltZetaTwoTerm
      (-(1 / 2 : ℝ) * (Real.pi ^ 2 / 6)) := by
  convert RamanujanChallenge.P26.alternatingZeta2Term26_hasSum using 1
  · funext n
    unfold testAltZetaTwoTerm
      RamanujanChallenge.P26.alternatingZeta2Term26
      RamanujanChallenge.P26.zeta2Term26
    push_cast
    ring
  · ring

private def testAltZetaTwoTriangleTerm (n : ℕ) : ℝ :=
  testAltZetaTwoTerm n *
    ∑ k ∈ Finset.range n, testAltZetaTwoTerm k

private theorem testAltZetaTwoTerm_abs (n : ℕ) :
    |testAltZetaTwoTerm n| = 1 / ((n : ℝ) + 1) ^ 2 := by
  unfold testAltZetaTwoTerm
  rw [abs_div, abs_pow]
  norm_num

private theorem testAltZetaTwoTriangleTerm_summable :
    Summable testAltZetaTwoTriangleTerm := by
  apply (shifted_zeta_two_hasSum.summable.mul_left
    (Real.pi ^ 2 / 6)).of_norm_bounded
  intro n
  unfold testAltZetaTwoTriangleTerm
  rw [Real.norm_eq_abs, abs_mul, testAltZetaTwoTerm_abs]
  have hinner :
      |∑ k ∈ Finset.range n, testAltZetaTwoTerm k| ≤
        Real.pi ^ 2 / 6 := by
    calc
      |∑ k ∈ Finset.range n, testAltZetaTwoTerm k| ≤
          ∑ k ∈ Finset.range n, |testAltZetaTwoTerm k| :=
        Finset.abs_sum_le_sum_abs _ _
      _ = harmonicSquare24 n := by
        unfold harmonicSquare24
        apply Finset.sum_congr rfl
        intro k _
        rw [testAltZetaTwoTerm_abs]
      _ ≤ Real.pi ^ 2 / 6 := harmonicSquare24_le_zeta_two n
  change
    1 / ((n : ℝ) + 1) ^ 2 *
        |∑ k ∈ Finset.range n, testAltZetaTwoTerm k| ≤
      Real.pi ^ 2 / 6 * (1 / ((n : ℝ) + 1) ^ 2)
  nlinarith [show 0 ≤ 1 / ((n : ℝ) + 1) ^ 2 by positivity]

private theorem testAltZetaTwoTerm_sq_hasSum :
    HasSum (fun n : ℕ => testAltZetaTwoTerm n ^ 2)
      (Real.pi ^ 4 / 90) := by
  convert shifted_zeta_four_hasSum24 using 1
  funext n
  unfold testAltZetaTwoTerm
  rw [div_pow]
  have hsign : ((-1 : ℝ) ^ (n + 1)) ^ 2 = 1 :=
    quadAlt_neg_one_pow_sq (n + 1)
  rw [hsign]
  push_cast
  ring

private theorem testAltZetaTwoTriangle_partial (N : ℕ) :
    2 * ∑ n ∈ Finset.range N, testAltZetaTwoTriangleTerm n =
      (∑ n ∈ Finset.range N, testAltZetaTwoTerm n) ^ 2 -
        ∑ n ∈ Finset.range N, testAltZetaTwoTerm n ^ 2 := by
  induction N with
  | zero => simp
  | succ N ih =>
      rw [Finset.sum_range_succ, Finset.sum_range_succ,
        Finset.sum_range_succ]
      unfold testAltZetaTwoTriangleTerm at *
      rw [mul_add, ih]
      ring_nf

private theorem testAltZetaTwoTriangleTerm_hasSum :
    HasSum testAltZetaTwoTriangleTerm
      (-(3 / 16 : ℝ) * (Real.pi ^ 4 / 90)) := by
  have hleft := testAltZetaTwoTerm_hasSum.tendsto_sum_nat
  have hdiag := testAltZetaTwoTerm_sq_hasSum.tendsto_sum_nat
  have hraw := (hleft.pow 2).sub hdiag
  have htwice :
      Tendsto
        (fun N : ℕ =>
          2 * ∑ n ∈ Finset.range N, testAltZetaTwoTriangleTerm n)
        atTop
        (𝓝 ((-(1 / 2 : ℝ) * (Real.pi ^ 2 / 6)) ^ 2 -
          Real.pi ^ 4 / 90)) := by
    apply hraw.congr'
    filter_upwards with N
    exact (testAltZetaTwoTriangle_partial N).symm
  have hhalf := htwice.const_mul (1 / 2 : ℝ)
  have hlim :
      Tendsto
        (fun N : ℕ =>
          ∑ n ∈ Finset.range N, testAltZetaTwoTriangleTerm n)
        atTop
        (𝓝 (-(3 / 16 : ℝ) * (Real.pi ^ 4 / 90))) := by
    convert hhalf using 1 <;> ring
  exact
    (testAltZetaTwoTriangleTerm_summable.hasSum_iff_tendsto_nat).2 hlim

private theorem test_neg_one_pow_sub {k n : ℕ} (hkn : k ≤ n) :
    (-1 : ℝ) ^ (n - k) = (-1 : ℝ) ^ n * (-1 : ℝ) ^ k := by
  have hsplit :
      (-1 : ℝ) ^ n = (-1 : ℝ) ^ (n - k) * (-1 : ℝ) ^ k := by
    rw [← pow_add, Nat.sub_add_cancel hkn]
  calc
    (-1 : ℝ) ^ (n - k) =
        (-1 : ℝ) ^ (n - k) * (((-1 : ℝ) ^ k) ^ 2) := by
          rw [quadAlt_neg_one_pow_sq]
          ring
    _ = ((-1 : ℝ) ^ (n - k) * (-1 : ℝ) ^ k) *
        (-1 : ℝ) ^ k := by ring
    _ = (-1 : ℝ) ^ n * (-1 : ℝ) ^ k := by rw [hsplit]

private def testDilogPlusCoeff (n : ℕ) : ℝ :=
  ∑ k ∈ Finset.range (n + 1),
    (1 / ((k : ℝ) + 1) ^ 2) * (-1 : ℝ) ^ (n - k)

private theorem testDilogPlusCoeff_abs_le (n : ℕ) :
    |testDilogPlusCoeff n| ≤ Real.pi ^ 2 / 6 := by
  calc
    |testDilogPlusCoeff n| ≤
        ∑ k ∈ Finset.range (n + 1),
          |(1 / ((k : ℝ) + 1) ^ 2) * (-1 : ℝ) ^ (n - k)| := by
      unfold testDilogPlusCoeff
      exact Finset.abs_sum_le_sum_abs _ _
    _ = harmonicSquare24 (n + 1) := by
      unfold harmonicSquare24
      apply Finset.sum_congr rfl
      intro k _
      rw [abs_mul, abs_pow]
      norm_num
    _ ≤ Real.pi ^ 2 / 6 := harmonicSquare24_le_zeta_two (n + 1)

private theorem testDilogPlusCoeff_integral_eq_triangle (n : ℕ) :
    -(testDilogPlusCoeff n) / ((n : ℝ) + 2) ^ 2 =
      testAltZetaTwoTriangleTerm (n + 1) := by
  unfold testDilogPlusCoeff testAltZetaTwoTriangleTerm
  rw [Finset.mul_sum, ← Finset.sum_neg_distrib,
    Finset.sum_div]
  apply Finset.sum_congr rfl
  intro k hk
  have hkn : k ≤ n := Nat.le_of_lt_succ (Finset.mem_range.mp hk)
  rw [test_neg_one_pow_sub hkn]
  unfold testAltZetaTwoTerm
  push_cast
  have hn : (n : ℝ) + 2 ≠ 0 := by positivity
  have hk' : (k : ℝ) + 1 ≠ 0 := by positivity
  field_simp [hn, hk']
  have hs1 := quadAlt_neg_one_pow_sq (n + 1)
  have hs2 := quadAlt_neg_one_pow_sq (k + 1)
  ring_nf at hs1 hs2 ⊢

private def testDilogPlusMoment (n : ℕ) (x : ℝ) : ℝ :=
  testDilogPlusCoeff n * (x ^ (n + 1) * Real.log x)

private theorem testDilogPlusMoment_intervalIntegrable (n : ℕ) :
    IntervalIntegrable (testDilogPlusMoment n)
      MeasureTheory.volume 0 1 := by
  unfold testDilogPlusMoment
  exact
    ((intervalIntegral.intervalIntegrable_log').continuousOn_mul
      (continuousOn_pow (n + 1))).const_mul (testDilogPlusCoeff n)

private theorem testDilogPlusMoment_integral (n : ℕ) :
    (∫ x : ℝ in 0..1, testDilogPlusMoment n x) =
      testAltZetaTwoTriangleTerm (n + 1) := by
  unfold testDilogPlusMoment
  rw [intervalIntegral.integral_const_mul,
    RamanujanChallenge.P26.integral_pow_mul_log26 (n + 1)]
  convert testDilogPlusCoeff_integral_eq_triangle n using 1 <;>
    push_cast <;> ring

private theorem testDilogPlusMoment_integral_norm (n : ℕ) :
    (∫ x : ℝ in 0..1, ‖testDilogPlusMoment n x‖) =
      |testDilogPlusCoeff n| / ((n : ℝ) + 2) ^ 2 := by
  calc
    (∫ x : ℝ in 0..1, ‖testDilogPlusMoment n x‖) =
        ∫ x : ℝ in 0..1,
          |testDilogPlusCoeff n| *
            (-(x ^ (n + 1) * Real.log x)) := by
      apply intervalIntegral.integral_congr
      intro x hx
      rw [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)] at hx
      unfold testDilogPlusMoment
      change
        |testDilogPlusCoeff n * (x ^ (n + 1) * Real.log x)| = _
      have hnonpos : x ^ (n + 1) * Real.log x ≤ 0 :=
        mul_nonpos_of_nonneg_of_nonpos
        (pow_nonneg hx.1 (n + 1))
        (Real.log_nonpos hx.1 hx.2)
      rw [abs_mul, abs_of_nonpos hnonpos]
    _ = |testDilogPlusCoeff n| / ((n : ℝ) + 2) ^ 2 := by
      rw [intervalIntegral.integral_const_mul,
        intervalIntegral.integral_neg,
        RamanujanChallenge.P26.integral_pow_mul_log26 (n + 1)]
      push_cast
      ring

private theorem testDilogPlusMoment_integral_norm_summable :
    Summable (fun n : ℕ =>
      ∫ x : ℝ in 0..1, ‖testDilogPlusMoment n x‖) := by
  have hmajor : Summable
      (fun n : ℕ => (Real.pi ^ 2 / 6) /
        ((n : ℝ) + 2) ^ 2) := by
    have htail : Summable
        (fun n : ℕ => 1 / (((n + 1 : ℕ) : ℝ) + 1) ^ 2) := by
      have hraw :=
        shifted_zeta_two_hasSum.summable.comp_injective Nat.succ_injective
      exact hraw.congr fun n => by
        dsimp only [Function.comp_apply]
    exact (htail.mul_left (Real.pi ^ 2 / 6)).congr fun n => by
      push_cast
      ring
  apply hmajor.of_nonneg_of_le
  · intro n
    exact intervalIntegral.integral_nonneg (by norm_num) fun _ _ => norm_nonneg _
  · intro n
    rw [testDilogPlusMoment_integral_norm]
    calc
      |testDilogPlusCoeff n| / ((n : ℝ) + 2) ^ 2 ≤
          (Real.pi ^ 2 / 6) / ((n : ℝ) + 2) ^ 2 := by
        exact mul_le_mul_of_nonneg_right
          (testDilogPlusCoeff_abs_le n) (by positivity)
      _ = (Real.pi ^ 2 / 6) / ((n : ℝ) + 2) ^ 2 := rfl

private theorem testDilogPlusMoment_hasSum_pointwise
    {x : ℝ} (hx0 : 0 < x) (hx1 : x < 1) :
    HasSum (fun n : ℕ => testDilogPlusMoment n x)
      (dilog x * Real.log x / (1 + x)) := by
  let f : ℕ → ℝ := fun n => x ^ (n + 1) / ((n : ℝ) + 1) ^ 2
  let g : ℕ → ℝ := fun n => (-x) ^ n
  have hf : HasSum f (dilog x) := by
    simpa only [f, dilog, Nat.cast_add, Nat.cast_one] using
      (dilog_summable (by rw [abs_of_pos hx0]; exact hx1.le)).hasSum
  have hg : HasSum g (1 + x)⁻¹ := by
    simpa [g] using hasSum_geometric_of_norm_lt_one
      (show ‖-x‖ < 1 by simpa [Real.norm_eq_abs, abs_of_pos hx0] using hx1)
  have hproduct :
      HasSum
        (fun n : ℕ =>
          ∑ k ∈ Finset.range (n + 1), f k * g (n - k))
        ((∑' n, f n) * ∑' n, g n) := by
    apply hasSum_sum_range_mul_of_summable_norm
    · exact hf.summable.norm
    · exact hg.summable.norm
  rw [hf.tsum_eq, hg.tsum_eq] at hproduct
  have hcoeff (n : ℕ) :
      (∑ k ∈ Finset.range (n + 1), f k * g (n - k)) =
        testDilogPlusCoeff n * x ^ (n + 1) := by
    unfold testDilogPlusCoeff f g
    rw [Finset.sum_mul]
    apply Finset.sum_congr rfl
    intro k hk
    have hkn : k ≤ n := Nat.le_of_lt_succ (Finset.mem_range.mp hk)
    rw [neg_pow]
    have hpow : x ^ (k + 1) * x ^ (n - k) = x ^ (n + 1) := by
      rw [← pow_add]
      congr 1
      omega
    calc
      x ^ (k + 1) / ((k : ℝ) + 1) ^ 2 *
          ((-1 : ℝ) ^ (n - k) * x ^ (n - k)) =
          (1 / ((k : ℝ) + 1) ^ 2) * (-1 : ℝ) ^ (n - k) *
            (x ^ (k + 1) * x ^ (n - k)) := by ring
      _ = _ := by rw [hpow]
  have hseries :
      HasSum (fun n : ℕ => testDilogPlusCoeff n * x ^ (n + 1))
        (dilog x * (1 + x)⁻¹) := by
    convert hproduct using 1
    funext n
    exact (hcoeff n).symm
  have hlog := hseries.mul_right (Real.log x)
  convert hlog using 1
  · funext n
    unfold testDilogPlusMoment
    ring
  · field_simp [show 1 + x ≠ 0 by linarith]

private theorem testDilogPlusIntegral_hasSum :
    HasSum (fun n : ℕ => testAltZetaTwoTriangleTerm (n + 1))
      (∫ x : ℝ in 0..1, dilog x * Real.log x / (1 + x)) := by
  have hInt : ∀ n : ℕ, MeasureTheory.Integrable (testDilogPlusMoment n)
      (MeasureTheory.volume.restrict (Set.Ioc (0 : ℝ) 1)) :=
    fun n => (testDilogPlusMoment_intervalIntegrable n).1
  have hNorm : Summable (fun n : ℕ =>
      ∫ x : ℝ in Set.Ioc (0 : ℝ) 1, ‖testDilogPlusMoment n x‖) := by
    simpa only [← intervalIntegral.integral_of_le
      (by norm_num : (0 : ℝ) ≤ 1)] using
      testDilogPlusMoment_integral_norm_summable
  have h := MeasureTheory.hasSum_integral_of_summable_integral_norm
    (μ := MeasureTheory.volume.restrict (Set.Ioc (0 : ℝ) 1)) hInt hNorm
  have h' : HasSum (fun n : ℕ => testAltZetaTwoTriangleTerm (n + 1))
      (∫ x : ℝ in Set.Ioc (0 : ℝ) 1,
        ∑' n : ℕ, testDilogPlusMoment n x) := by
    convert h using 1
    funext n
    rw [← intervalIntegral.integral_of_le (by norm_num : (0 : ℝ) ≤ 1)]
    exact (testDilogPlusMoment_integral n).symm
  convert h' using 1
  rw [intervalIntegral.integral_of_le (by norm_num : (0 : ℝ) ≤ 1)]
  apply MeasureTheory.setIntegral_congr_ae measurableSet_Ioc
  filter_upwards [MeasureTheory.Measure.ae_ne MeasureTheory.volume (1 : ℝ)]
    with x hxne hx
  have hxlt : x < 1 := lt_of_le_of_ne hx.2 hxne
  exact (testDilogPlusMoment_hasSum_pointwise hx.1 hxlt).tsum_eq.symm

private theorem testDilogPlusIntegral :
    (∫ x : ℝ in 0..1, dilog x * Real.log x / (1 + x)) =
      -(3 / 16 : ℝ) * (Real.pi ^ 4 / 90) := by
  have htail :
      HasSum (fun n : ℕ => testAltZetaTwoTriangleTerm (n + 1))
        (-(3 / 16 : ℝ) * (Real.pi ^ 4 / 90)) := by
    apply (hasSum_nat_add_iff (f := testAltZetaTwoTriangleTerm) 1).mpr
    convert testAltZetaTwoTriangleTerm_hasSum using 1
    simp [testAltZetaTwoTriangleTerm]
  exact testDilogPlusIntegral_hasSum.unique htail

private def testDilogRadialMoment (n : ℕ) (x : ℝ) : ℝ :=
  (x ^ (n + 1) / ((n : ℝ) + 1) ^ 2) * (Real.log x / x)

private theorem testDilogRadialMoment_intervalIntegrable (n : ℕ) :
    IntervalIntegrable (testDilogRadialMoment n)
      MeasureTheory.volume 0 1 := by
  apply IntervalIntegrable.congr
    (f := fun x : ℝ =>
      (1 / ((n : ℝ) + 1) ^ 2) * (x ^ n * Real.log x)) ?_
    (((intervalIntegral.intervalIntegrable_log').continuousOn_mul
      (continuousOn_pow n)).const_mul (1 / ((n : ℝ) + 1) ^ 2))
  intro x hx
  have hx' : x ∈ Ioc (0 : ℝ) 1 := by
    simpa [Set.uIoc_of_le (by norm_num : (0 : ℝ) ≤ 1)] using hx
  unfold testDilogRadialMoment
  field_simp [ne_of_gt hx'.1]
  ring

private theorem testDilogRadialMoment_integral (n : ℕ) :
    (∫ x : ℝ in 0..1, testDilogRadialMoment n x) =
      -(1 / ((n : ℝ) + 1) ^ 4) := by
  calc
    (∫ x : ℝ in 0..1, testDilogRadialMoment n x) =
        ∫ x : ℝ in 0..1,
          (1 / ((n : ℝ) + 1) ^ 2) *
            (x ^ n * Real.log x) := by
      apply intervalIntegral.integral_congr
      intro x hx
      by_cases hxzero : x = 0
      · subst x
        simp [testDilogRadialMoment]
      · unfold testDilogRadialMoment
        field_simp [hxzero]
        ring
    _ = -(1 / ((n : ℝ) + 1) ^ 4) := by
      rw [intervalIntegral.integral_const_mul,
        RamanujanChallenge.P26.integral_pow_mul_log26 n]
      have hn : (n : ℝ) + 1 ≠ 0 := by positivity
      field_simp [hn]

private theorem testDilogRadialMoment_nonpos
    (n : ℕ) {x : ℝ} (hx0 : 0 ≤ x) (hx1 : x ≤ 1) :
    testDilogRadialMoment n x ≤ 0 := by
  by_cases hxzero : x = 0
  · subst x
    simp [testDilogRadialMoment]
  · unfold testDilogRadialMoment
    exact mul_nonpos_of_nonneg_of_nonpos
      (div_nonneg (pow_nonneg hx0 _) (sq_nonneg _))
      (div_nonpos_of_nonpos_of_nonneg
        (Real.log_nonpos hx0 hx1) hx0)

private theorem testDilogRadialMoment_integral_norm (n : ℕ) :
    (∫ x : ℝ in 0..1, ‖testDilogRadialMoment n x‖) =
      1 / ((n : ℝ) + 1) ^ 4 := by
  calc
    (∫ x : ℝ in 0..1, ‖testDilogRadialMoment n x‖) =
        ∫ x : ℝ in 0..1, -(testDilogRadialMoment n x) := by
      apply intervalIntegral.integral_congr
      intro x hx
      rw [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)] at hx
      change |testDilogRadialMoment n x| = -testDilogRadialMoment n x
      rw [abs_of_nonpos
        (testDilogRadialMoment_nonpos n hx.1 hx.2)]
    _ = 1 / ((n : ℝ) + 1) ^ 4 := by
      rw [intervalIntegral.integral_neg,
        testDilogRadialMoment_integral]
      ring

private theorem testDilogRadialMoment_hasSum_pointwise
    {x : ℝ} (hx0 : 0 < x) (hx1 : x < 1) :
    HasSum (fun n : ℕ => testDilogRadialMoment n x)
      (dilog x * Real.log x / x) := by
  have hd : HasSum
      (fun n : ℕ => x ^ (n + 1) / ((n : ℝ) + 1) ^ 2)
      (dilog x) := by
    simpa only [dilog, Nat.cast_add, Nat.cast_one] using
      (dilog_summable (by rw [abs_of_pos hx0]; exact hx1.le)).hasSum
  have hmul := hd.mul_right (Real.log x / x)
  convert hmul using 1
  ring

private theorem testDilogRadialIntegral_hasSum :
    HasSum (fun n : ℕ => -(1 / ((n : ℝ) + 1) ^ 4))
      (∫ x : ℝ in 0..1, dilog x * Real.log x / x) := by
  have hInt : ∀ n : ℕ, MeasureTheory.Integrable (testDilogRadialMoment n)
      (MeasureTheory.volume.restrict (Set.Ioc (0 : ℝ) 1)) :=
    fun n => (testDilogRadialMoment_intervalIntegrable n).1
  have hNorm : Summable (fun n : ℕ =>
      ∫ x : ℝ in Set.Ioc (0 : ℝ) 1, ‖testDilogRadialMoment n x‖) := by
    simpa only [← intervalIntegral.integral_of_le
      (by norm_num : (0 : ℝ) ≤ 1)] using
      (shifted_zeta_four_hasSum24.summable.congr fun n =>
        (testDilogRadialMoment_integral_norm n).symm)
  have h := MeasureTheory.hasSum_integral_of_summable_integral_norm
    (μ := MeasureTheory.volume.restrict (Set.Ioc (0 : ℝ) 1)) hInt hNorm
  have h' : HasSum (fun n : ℕ => -(1 / ((n : ℝ) + 1) ^ 4))
      (∫ x : ℝ in Set.Ioc (0 : ℝ) 1,
        ∑' n : ℕ, testDilogRadialMoment n x) := by
    convert h using 1
    funext n
    rw [← intervalIntegral.integral_of_le (by norm_num : (0 : ℝ) ≤ 1)]
    exact (testDilogRadialMoment_integral n).symm
  convert h' using 1
  rw [intervalIntegral.integral_of_le (by norm_num : (0 : ℝ) ≤ 1)]
  apply MeasureTheory.setIntegral_congr_ae measurableSet_Ioc
  filter_upwards [MeasureTheory.Measure.ae_ne MeasureTheory.volume (1 : ℝ)]
    with x hxne hx
  have hxlt : x < 1 := lt_of_le_of_ne hx.2 hxne
  exact (testDilogRadialMoment_hasSum_pointwise hx.1 hxlt).tsum_eq.symm

private theorem testDilogRadialIntegral :
    (∫ x : ℝ in 0..1, dilog x * Real.log x / x) =
      -(Real.pi ^ 4 / 90) := by
  have hzeta :
      HasSum (fun n : ℕ => -(1 / ((n : ℝ) + 1) ^ 4))
        (-(Real.pi ^ 4 / 90)) := by
    convert shifted_zeta_four_hasSum24.neg using 1
  exact testDilogRadialIntegral_hasSum.unique hzeta

private theorem testMobiusIntegral (g : ℝ → ℝ) :
    (∫ x : ℝ in 0..1, g x) =
      ∫ t : ℝ in 0..1,
        g (t / (2 - t)) * (2 / (2 - t) ^ 2) := by
  let f : ℝ → ℝ := fun t => t / (2 - t)
  let f' : ℝ → ℝ := fun t => 2 / (2 - t) ^ 2
  have hf : ContinuousOn f (Set.uIcc (0 : ℝ) 1) := by
    unfold f
    apply ContinuousOn.div
    · exact continuousOn_id
    · exact continuousOn_const.sub continuousOn_id
    · intro t ht
      have ht1 : t ≤ 1 := by
        simpa [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)] using ht.2
      linarith
  have hff' : ∀ t ∈ Set.Ioo (0 : ℝ) 1,
      HasDerivAt f (f' t) t := by
    intro t ht
    have hden : HasDerivAt (fun u : ℝ => 2 - u) (-1) t := by
      simpa using (hasDerivAt_const t (2 : ℝ)).sub (hasDerivAt_id t)
    have h2mt : 2 - t ≠ 0 := by linarith [ht.2]
    have hinv : HasDerivAt (fun u : ℝ => (2 - u)⁻¹)
        (1 / (2 - t) ^ 2) t := by
      have hc := (hasDerivAt_inv (x := 2 - t) h2mt).comp t hden
      convert hc using 1
      field_simp [h2mt]
    have hm : HasDerivAt (fun u : ℝ => u * (2 - u)⁻¹)
        (1 * (2 - t)⁻¹ + t * (1 / (2 - t) ^ 2)) t :=
      (hasDerivAt_id t).mul hinv
    unfold f f'
    have hval : 2 / (2 - t) ^ 2 =
        1 * (2 - t)⁻¹ + t * (1 / (2 - t) ^ 2) := by
      field_simp [h2mt]
      ring
    rw [hval]
    simpa [div_eq_mul_inv] using hm
  have hf' : ∀ t ∈ Set.Ioo (0 : ℝ) 1, 0 ≤ f' t := by
    intro t _
    unfold f'
    positivity
  have hsubst := intervalIntegral.integral_comp_mul_deriv_of_deriv_nonneg
    (a := (0 : ℝ)) (b := 1) (f := f) (f' := f') (g := g)
    hf (by simpa using hff') (by simpa using hf')
  have hf01 : f 0 = 0 := by unfold f; ring
  have hf11 : f 1 = 1 := by unfold f; ring
  rw [hf01, hf11] at hsubst
  simpa [f, f', Function.comp_apply] using hsubst.symm

private theorem testI12_substitution :
    I12 - I22 =
      2 * ∫ x : ℝ in 0..1,
        quadAltV x * Real.log (1 - x) / (1 + x) := by
  let g : ℝ → ℝ := fun x =>
    2 * quadAltV x * Real.log (1 - x) / (1 + x)
  have hsub := testMobiusIntegral g
  have hrows : I12 - I22 =
      ∫ t : ℝ in 0..1,
        W0 t * H1 t / (2 - t) -
          W0 t * H2 t / (2 - t) := by
    unfold I12 I22
    rw [← intervalIntegral.integral_sub
      quadAltI12_kernel_intervalIntegrable
      quadAltI22_kernel_intervalIntegrable]
  have hright :
      2 * ∫ x : ℝ in 0..1,
          quadAltV x * Real.log (1 - x) / (1 + x) =
        ∫ x : ℝ in 0..1, g x := by
    rw [← intervalIntegral.integral_const_mul]
    apply intervalIntegral.integral_congr
    intro x _
    dsimp only [g]
    ring
  rw [hrows, hright, hsub]
  apply intervalIntegral.integral_congr_ae
  filter_upwards [MeasureTheory.Measure.ae_ne MeasureTheory.volume (1 : ℝ)]
    with t htne ht
  have ht' : t ∈ Ioc (0 : ℝ) 1 := by
    simpa [Set.uIoc_of_le (by norm_num : (0 : ℝ) ≤ 1)] using ht
  have ht0 : 0 < t := ht'.1
  have ht1 : t < 1 := lt_of_le_of_ne ht'.2 htne
  have h2mt : 0 < 2 - t := by linarith
  have hx0 : 0 < t / (2 - t) := by positivity
  have hx1 : t / (2 - t) < 1 := by
    rw [div_lt_one h2mt]
    linarith
  have hnorm :
      2 * (t / (2 - t)) / (1 + t / (2 - t)) = t := by
    field_simp [ne_of_gt h2mt]
    ring
  have hW := quadAlt_neg2V_eq_W0 hx0 hx1
  rw [hnorm] at hW
  have hlog :
      Real.log (1 - t / (2 - t)) = H2 t - H1 t := by
    have harg :
        1 - t / (2 - t) = (1 - t) / (1 - t / 2) := by
      field_simp [ne_of_gt h2mt]
      ring
    rw [harg, Real.log_div (by linarith : 1 - t ≠ 0)
      (by linarith : 1 - t / 2 ≠ 0)]
    unfold H1 H2
    ring
  dsimp only [g]
  rw [← hW, hlog]
  have hplus : 1 + t / (2 - t) ≠ 0 := by positivity
  field_simp [ne_of_gt h2mt, hplus]
  ring

private theorem testI10I20_substitution :
    I10 - I20 =
      2 * ∫ x : ℝ in 0..1,
        quadAltV x * Real.log (1 - x) / (x * (1 + x)) := by
  let g : ℝ → ℝ := fun x =>
    2 * quadAltV x * Real.log (1 - x) / (x * (1 + x))
  have hsub := testMobiusIntegral g
  have hrows : I10 - I20 =
      ∫ t : ℝ in 0..1,
        W0 t * H1 t / t - W0 t * H2 t / t := by
    unfold I10 I20
    rw [← intervalIntegral.integral_sub
      quadAltI10_kernel_intervalIntegrable
      quadAltI20_kernel_intervalIntegrable]
  have hright :
      2 * ∫ x : ℝ in 0..1,
          quadAltV x * Real.log (1 - x) / (x * (1 + x)) =
        ∫ x : ℝ in 0..1, g x := by
    rw [← intervalIntegral.integral_const_mul]
    apply intervalIntegral.integral_congr
    intro x _
    dsimp only [g]
    ring
  rw [hrows, hright, hsub]
  apply intervalIntegral.integral_congr_ae
  filter_upwards [MeasureTheory.Measure.ae_ne MeasureTheory.volume (1 : ℝ)]
    with t htne ht
  have ht' : t ∈ Ioc (0 : ℝ) 1 := by
    simpa [Set.uIoc_of_le (by norm_num : (0 : ℝ) ≤ 1)] using ht
  have ht0 : 0 < t := ht'.1
  have ht1 : t < 1 := lt_of_le_of_ne ht'.2 htne
  have h2mt : 0 < 2 - t := by linarith
  have hx0 : 0 < t / (2 - t) := by positivity
  have hx1 : t / (2 - t) < 1 := by
    rw [div_lt_one h2mt]
    linarith
  have hnorm :
      2 * (t / (2 - t)) / (1 + t / (2 - t)) = t := by
    field_simp [ne_of_gt h2mt]
    ring
  have hW := quadAlt_neg2V_eq_W0 hx0 hx1
  rw [hnorm] at hW
  have hlog :
      Real.log (1 - t / (2 - t)) = H2 t - H1 t := by
    have harg :
        1 - t / (2 - t) = (1 - t) / (1 - t / 2) := by
      field_simp [ne_of_gt h2mt]
      ring
    rw [harg, Real.log_div (by linarith : 1 - t ≠ 0)
      (by linarith : 1 - t / 2 ≠ 0)]
    unfold H1 H2
    ring
  dsimp only [g]
  rw [← hW, hlog]
  have hplus : 1 + t / (2 - t) ≠ 0 := by positivity
  field_simp [ne_of_gt ht0, ne_of_gt h2mt, hplus]
  ring

private def testDilogSlope (x : ℝ) : ℝ :=
  Function.update (fun y : ℝ => dilog y / y) 0 1 x

private theorem testDilogSlope_continuousOn :
    ContinuousOn testDilogSlope (Icc (0 : ℝ) 1) := by
  intro x hx
  by_cases hxzero : x = 0
  · subst x
    have hc := RamanujanChallenge.P26.dilog_hasDerivAt_zero26.continuousAt_div
    have hc' : ContinuousAt testDilogSlope 0 := by
      convert hc using 1
      funext y
      by_cases hy : y = 0
      · subst y
        simp [testDilogSlope]
      · simp [testDilogSlope, hy, dilog_zero]
    exact hc'.continuousWithinAt
  · have hd : ContinuousWithinAt dilog (Icc (0 : ℝ) 1) x :=
      (dilog_continuousOn_unit.mono (by
        intro y hy
        constructor <;> linarith [hy.1, hy.2])) x hx
    have hbase : ContinuousWithinAt (fun y : ℝ => dilog y / y)
        (Icc (0 : ℝ) 1) x :=
      hd.div continuousWithinAt_id hxzero
    have heq : testDilogSlope =ᶠ[𝓝[Icc (0 : ℝ) 1] x]
        (fun y : ℝ => dilog y / y) := by
      filter_upwards [
        (eventually_ne_nhds hxzero).filter_mono nhdsWithin_le_nhds] with y hy
      simp [testDilogSlope, hy]
    exact hbase.congr_of_eventuallyEq heq (by
      simp [testDilogSlope, hxzero])

private theorem testDilogRadialKernel_intervalIntegrable :
    IntervalIntegrable
      (fun x : ℝ => dilog x * Real.log x / x)
      MeasureTheory.volume 0 1 := by
  have hslope : ContinuousOn testDilogSlope (Set.uIcc (0 : ℝ) 1) := by
    simpa [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)] using
      testDilogSlope_continuousOn
  have hprod := (intervalIntegral.intervalIntegrable_log').continuousOn_mul
    hslope
  apply IntervalIntegrable.congr
    (f := fun x : ℝ => testDilogSlope x * Real.log x) ?_ hprod
  intro x hx
  have hx' : x ∈ Ioc (0 : ℝ) 1 := by
    simpa [Set.uIoc_of_le (by norm_num : (0 : ℝ) ≤ 1)] using hx
  simp [testDilogSlope, ne_of_gt hx'.1]
  ring

private theorem testDilogDerivativeKernel_intervalIntegrable :
    IntervalIntegrable
      (fun x : ℝ => dilog x * Real.log x / (x * (1 + x)))
      MeasureTheory.volume 0 1 := by
  have hc : ContinuousOn (fun x : ℝ => 1 / (1 + x))
      (Set.uIcc (0 : ℝ) 1) := by
    rw [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)]
    intro x hx
    have hden : (1 : ℝ) + x ≠ 0 := by linarith [hx.1]
    exact continuousAt_const.div (continuousAt_const.add continuousAt_id)
      hden |>.continuousWithinAt
  have hprod := testDilogRadialKernel_intervalIntegrable.continuousOn_mul hc
  apply IntervalIntegrable.congr
    (f := fun x : ℝ => (1 / (1 + x)) *
      (dilog x * Real.log x / x)) ?_ hprod
  intro x hx
  have hx' : x ∈ Ioc (0 : ℝ) 1 := by
    simpa [Set.uIoc_of_le (by norm_num : (0 : ℝ) ≤ 1)] using hx
  have hxne : x ≠ 0 := ne_of_gt hx'.1
  have hplus : 1 + x ≠ 0 := by linarith [hx'.1]
  field_simp [hxne, hplus]

private def testVDilogPrimitive (x : ℝ) : ℝ :=
  quadAltV x * dilog x

private theorem testVDilogPrimitive_hasDerivAt
    {x : ℝ} (hx0 : 0 < x) (hx1 : x < 1) :
    HasDerivAt testVDilogPrimitive
      (dilog x * Real.log x / (x * (1 + x)) +
        quadAltV x * H1 x / x) x := by
  have hprod := (quadAltV_hasDerivAt hx0 hx1).mul
    (dilog_hasDerivAt hx0 hx1)
  unfold testVDilogPrimitive H1
  convert hprod using 1
  field_simp [ne_of_gt hx0, show 1 + x ≠ 0 by linarith]

private theorem testVDilogPrimitive_tendsto_zero :
    Tendsto testVDilogPrimitive (𝓝[>] (0 : ℝ)) (𝓝 0) := by
  have hslope : Tendsto (fun x : ℝ => dilog x / x)
      (𝓝[>] (0 : ℝ)) (𝓝 1) := by
    simpa [dilog_zero, div_eq_mul_inv, mul_comm] using
      RamanujanChallenge.P26.dilog_hasDerivAt_zero26.tendsto_slope_zero_right
  have hraw := quadAltV_mul_self_tendsto.mul hslope
  have hraw' : Tendsto testVDilogPrimitive
      (𝓝[>] (0 : ℝ)) (𝓝 (0 * 1)) :=
    hraw.congr' (by
      filter_upwards [self_mem_nhdsWithin] with x hx
      have hx0 : 0 < x := hx
      unfold testVDilogPrimitive
      field_simp [ne_of_gt hx0])
  simpa using hraw'

private theorem testVDilogPrimitive_tendsto_one :
    Tendsto testVDilogPrimitive (𝓝[<] (1 : ℝ)) (𝓝 0) := by
  have hsub : Tendsto (fun x : ℝ => x - 1)
      (𝓝[<] (1 : ℝ)) (𝓝 0) := by
    have hc : ContinuousAt (fun x : ℝ => x - 1) 1 := by fun_prop
    simpa using hc.tendsto.mono_left nhdsWithin_le_nhds
  have hVraw := quadAltV_slope_tendsto_one.mul hsub
  have hV : Tendsto quadAltV (𝓝[<] (1 : ℝ)) (𝓝 0) := by
    simpa only [zero_mul] using hVraw.congr' (by
      filter_upwards [self_mem_nhdsWithin] with x hx
      have hxlt : x < 1 := hx
      have hxne : x - 1 ≠ 0 := by linarith
      field_simp [hxne])
  have hd : Tendsto dilog (𝓝[<] (1 : ℝ))
      (𝓝 (Real.pi ^ 2 / 6)) := by
    have hc : ContinuousWithinAt dilog (Iio (1 : ℝ)) 1 :=
      (dilog_continuousOn_unit 1 (by norm_num)).mono_of_mem_nhdsWithin
        (Icc_mem_nhdsLT (show (-1 : ℝ) < 1 by norm_num))
    simpa [dilog_one] using hc.tendsto
  unfold testVDilogPrimitive
  simpa using hV.mul hd

private theorem testVH1Div_intervalIntegrable :
    IntervalIntegrable
      (fun x : ℝ => quadAltV x * H1 x / x)
      MeasureTheory.volume 0 1 := by
  have hlogPlus : ContinuousOn (fun x : ℝ => Real.log (1 + x))
      (Set.uIcc (0 : ℝ) 1) := by
    rw [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)]
    intro x hx
    exact (ContinuousAt.log (by fun_prop)
      (by linarith [hx.1] : 1 + x ≠ 0)).continuousWithinAt
  have hdilogNeg : ContinuousOn (fun x : ℝ => dilog (-x))
      (Set.uIcc (0 : ℝ) 1) := by
    rw [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)]
    apply dilog_continuousOn_unit.comp (by fun_prop)
    intro x hx
    constructor <;> linarith [hx.1, hx.2]
  let q : ℝ → ℝ := fun x =>
    (Real.log x ^ 2 * H1 x / x) / 2 -
      Real.log (1 + x) * (Real.log x * H1 x / x) -
      dilog (-x) * (H1 x / x) -
      (Real.pi ^ 2 / 12) * (H1 x / x)
  have hq : IntervalIntegrable q MeasureTheory.volume 0 1 := by
    exact (((i10LogSqH1Div_intervalIntegrable.div_const 2).sub
      (i10LogH1Div_intervalIntegrable.continuousOn_mul hlogPlus)).sub
      (i10H1Div_intervalIntegrable.continuousOn_mul hdilogNeg)).sub
      (i10H1Div_intervalIntegrable.const_mul (Real.pi ^ 2 / 12))
  apply IntervalIntegrable.congr (f := q) ?_ hq
  intro x hx
  have hx' : x ∈ Ioc (0 : ℝ) 1 := by
    simpa [Set.uIoc_of_le (by norm_num : (0 : ℝ) ≤ 1)] using hx
  have hxne : x ≠ 0 := ne_of_gt hx'.1
  dsimp only [q]
  unfold quadAltV
  field_simp [hxne]

private theorem testDilogPlusKernel_intervalIntegrable :
    IntervalIntegrable
      (fun x : ℝ => dilog x * Real.log x / (1 + x))
      MeasureTheory.volume 0 1 := by
  have hc : ContinuousOn (fun x : ℝ => x / (1 + x))
      (Set.uIcc (0 : ℝ) 1) := by
    rw [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)]
    intro x hx
    exact continuousAt_id.div (continuousAt_const.add continuousAt_id)
      (by linarith [hx.1] : 1 + x ≠ 0) |>.continuousWithinAt
  have hprod := testDilogRadialKernel_intervalIntegrable.continuousOn_mul hc
  apply IntervalIntegrable.congr
    (f := fun x : ℝ => (x / (1 + x)) *
      (dilog x * Real.log x / x)) ?_ hprod
  intro x hx
  have hx' : x ∈ Ioc (0 : ℝ) 1 := by
    simpa [Set.uIoc_of_le (by norm_num : (0 : ℝ) ≤ 1)] using hx
  have hxne : x ≠ 0 := ne_of_gt hx'.1
  have hplus : 1 + x ≠ 0 := by linarith [hx'.1]
  field_simp [hxne, hplus]

private theorem testVLogIntegral_eq_dilogDerivative :
    (∫ x : ℝ in 0..1,
      quadAltV x * Real.log (1 - x) / x) =
      ∫ x : ℝ in 0..1,
        dilog x * Real.log x / (x * (1 + x)) := by
  have hderivInt :=
    testDilogDerivativeKernel_intervalIntegrable.add
      testVH1Div_intervalIntegrable
  have hFTC := intervalIntegral.integral_eq_sub_of_hasDerivAt_of_tendsto
    (a := (0 : ℝ)) (b := 1) (f := testVDilogPrimitive)
    (f' := fun x : ℝ =>
      dilog x * Real.log x / (x * (1 + x)) +
        quadAltV x * H1 x / x)
    (by norm_num)
    (fun x hx => testVDilogPrimitive_hasDerivAt hx.1 hx.2)
    hderivInt testVDilogPrimitive_tendsto_zero
    testVDilogPrimitive_tendsto_one
  have hsplit := intervalIntegral.integral_add
    testDilogDerivativeKernel_intervalIntegrable
    testVH1Div_intervalIntegrable
  rw [hsplit] at hFTC
  norm_num at hFTC
  calc
    (∫ x : ℝ in 0..1,
      quadAltV x * Real.log (1 - x) / x) =
        -(∫ x : ℝ in 0..1, quadAltV x * H1 x / x) := by
      rw [← intervalIntegral.integral_neg]
      apply intervalIntegral.integral_congr
      intro x _
      unfold H1
      ring
    _ = ∫ x : ℝ in 0..1,
        dilog x * Real.log x / (x * (1 + x)) := by linarith

private theorem testDilogDerivativeIntegral :
    (∫ x : ℝ in 0..1,
      dilog x * Real.log x / (x * (1 + x))) =
      -(13 / 16 : ℝ) * (Real.pi ^ 4 / 90) := by
  calc
    (∫ x : ℝ in 0..1,
      dilog x * Real.log x / (x * (1 + x))) =
        (∫ x : ℝ in 0..1, dilog x * Real.log x / x) -
          ∫ x : ℝ in 0..1, dilog x * Real.log x / (1 + x) := by
      rw [← intervalIntegral.integral_sub
        testDilogRadialKernel_intervalIntegrable
        testDilogPlusKernel_intervalIntegrable]
      apply intervalIntegral.integral_congr_ae
      filter_upwards [MeasureTheory.Measure.ae_ne MeasureTheory.volume (0 : ℝ)]
        with x hxzero hx
      have hx' : x ∈ Ioc (0 : ℝ) 1 := by
        simpa [Set.uIoc_of_le (by norm_num : (0 : ℝ) ≤ 1)] using hx
      have hplus : 1 + x ≠ 0 := by linarith [hx'.1]
      field_simp [hxzero, hplus]
      ring
    _ = -(13 / 16 : ℝ) * (Real.pi ^ 4 / 90) := by
      rw [testDilogRadialIntegral, testDilogPlusIntegral]
      ring

private theorem testEndpointDifferenceSum :
    (I12 - I22) + (I10 - I20) =
      2 * ∫ x : ℝ in 0..1,
        quadAltV x * Real.log (1 - x) / x := by
  let c : ℝ → ℝ := fun x =>
    quadAltV x * Real.log (1 - x) / x
  have hc : IntervalIntegrable c MeasureTheory.volume 0 1 := by
    apply IntervalIntegrable.congr
      (f := fun x : ℝ => -(quadAltV x * H1 x / x)) ?_
      testVH1Div_intervalIntegrable.neg
    intro x _
    dsimp only [c]
    unfold H1
    ring
  have hfacA : ContinuousOn (fun x : ℝ => x / (1 + x))
      (Set.uIcc (0 : ℝ) 1) := by
    rw [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)]
    intro x hx
    exact continuousAt_id.div (continuousAt_const.add continuousAt_id)
      (by linarith [hx.1] : 1 + x ≠ 0) |>.continuousWithinAt
  have hfacB : ContinuousOn (fun x : ℝ => 1 / (1 + x))
      (Set.uIcc (0 : ℝ) 1) := by
    rw [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)]
    intro x hx
    exact continuousAt_const.div (continuousAt_const.add continuousAt_id)
      (by linarith [hx.1] : 1 + x ≠ 0) |>.continuousWithinAt
  let a : ℝ → ℝ := fun x =>
    quadAltV x * Real.log (1 - x) / (1 + x)
  let b : ℝ → ℝ := fun x =>
    quadAltV x * Real.log (1 - x) / (x * (1 + x))
  have ha : IntervalIntegrable a MeasureTheory.volume 0 1 := by
    have hprod := hc.continuousOn_mul hfacA
    apply IntervalIntegrable.congr
      (f := fun x : ℝ => (x / (1 + x)) * c x) ?_ hprod
    intro x hx
    have hx' : x ∈ Ioc (0 : ℝ) 1 := by
      simpa [Set.uIoc_of_le (by norm_num : (0 : ℝ) ≤ 1)] using hx
    have hxne : x ≠ 0 := ne_of_gt hx'.1
    have hplus : 1 + x ≠ 0 := by linarith [hx'.1]
    dsimp only [a, c]
    field_simp [hxne, hplus]
  have hb : IntervalIntegrable b MeasureTheory.volume 0 1 := by
    have hprod := hc.continuousOn_mul hfacB
    apply IntervalIntegrable.congr
      (f := fun x : ℝ => (1 / (1 + x)) * c x) ?_ hprod
    intro x hx
    have hx' : x ∈ Ioc (0 : ℝ) 1 := by
      simpa [Set.uIoc_of_le (by norm_num : (0 : ℝ) ≤ 1)] using hx
    have hxne : x ≠ 0 := ne_of_gt hx'.1
    have hplus : 1 + x ≠ 0 := by linarith [hx'.1]
    dsimp only [b, c]
    field_simp [hxne, hplus]
  have hab :
      (∫ x : ℝ in 0..1, a x) + (∫ x : ℝ in 0..1, b x) =
        ∫ x : ℝ in 0..1, c x := by
    rw [← intervalIntegral.integral_add ha hb]
    apply intervalIntegral.integral_congr_ae
    filter_upwards [MeasureTheory.Measure.ae_ne MeasureTheory.volume (0 : ℝ)]
      with x hxzero hx
    have hx' : x ∈ Ioc (0 : ℝ) 1 := by
      simpa [Set.uIoc_of_le (by norm_num : (0 : ℝ) ≤ 1)] using hx
    have hplus : 1 + x ≠ 0 := by linarith [hx'.1]
    dsimp only [a, b, c]
    field_simp [hxzero, hplus]
    ring
  rw [testI12_substitution, testI10I20_substitution]
  dsimp only [a, b, c] at hab
  linarith

theorem quadAltI12_eq :
    I12 = -(3 / 2) * quadAltK
      + (3 / 2) * Real.log 2 ^ 2 * (Real.pi ^ 2 / 6)
      - (9 / 20) * (Real.pi ^ 2 / 6) ^ 2 := by
  have hsum := testEndpointDifferenceSum
  rw [testVLogIntegral_eq_dilogDerivative,
    testDilogDerivativeIntegral, quadAltI10_eq, quadAltI20_eq,
    quadAltI22_eq] at hsum
  have hzeta : Real.pi ^ 4 / 90 =
      (2 / 5 : ℝ) * (Real.pi ^ 2 / 6) ^ 2 := by ring
  rw [hzeta] at hsum
  linarith

private theorem testI21_substitution :
    I21 - I22 =
      -2 * ∫ x : ℝ in 0..1,
        quadAltV x * Real.log (1 + x) / (1 - x) := by
  let g : ℝ → ℝ := fun x =>
    -2 * quadAltV x * Real.log (1 + x) / (1 - x)
  have hsub := testMobiusIntegral g
  have hrows : I21 - I22 =
      ∫ t : ℝ in 0..1,
        W0 t * H2 t / (1 - t) -
          W0 t * H2 t / (2 - t) := by
    unfold I21 I22
    rw [← intervalIntegral.integral_sub
      quadAltI21_kernel_intervalIntegrable
      quadAltI22_kernel_intervalIntegrable]
  have hright :
      -2 * ∫ x : ℝ in 0..1,
          quadAltV x * Real.log (1 + x) / (1 - x) =
        ∫ x : ℝ in 0..1, g x := by
    rw [← intervalIntegral.integral_const_mul]
    apply intervalIntegral.integral_congr
    intro x _
    dsimp only [g]
    ring
  rw [hrows, hright, hsub]
  apply intervalIntegral.integral_congr_ae
  filter_upwards [MeasureTheory.Measure.ae_ne MeasureTheory.volume (1 : ℝ)]
    with t htne ht
  have ht' : t ∈ Ioc (0 : ℝ) 1 := by
    simpa [Set.uIoc_of_le (by norm_num : (0 : ℝ) ≤ 1)] using ht
  have ht0 : 0 < t := ht'.1
  have ht1 : t < 1 := lt_of_le_of_ne ht'.2 htne
  have h2mt : 0 < 2 - t := by linarith
  have hx0 : 0 < t / (2 - t) := by positivity
  have hx1 : t / (2 - t) < 1 := by
    rw [div_lt_one h2mt]
    linarith
  have hnorm :
      2 * (t / (2 - t)) / (1 + t / (2 - t)) = t := by
    field_simp [ne_of_gt h2mt]
    ring
  have hW := quadAlt_neg2V_eq_W0 hx0 hx1
  rw [hnorm] at hW
  have hlog : Real.log (1 + t / (2 - t)) = H2 t := by
    have harg :
        1 + t / (2 - t) = (1 - t / 2)⁻¹ := by
      field_simp [ne_of_gt h2mt, (by linarith : 1 - t / 2 ≠ 0)]
      ring
    rw [harg, Real.log_inv]
    unfold H2
    rfl
  have honeSub :
      1 - t / (2 - t) = 2 * (1 - t) / (2 - t) := by
    field_simp [ne_of_gt h2mt]
    ring
  dsimp only [g]
  rw [← hW, hlog, honeSub]
  field_simp [ne_of_gt h2mt, (by linarith : 1 - t ≠ 0)]
  ring

private theorem testVLogMinusPlus_intervalIntegrable :
    IntervalIntegrable
      (fun x : ℝ => quadAltV x * Real.log (1 - x) / (1 + x))
      MeasureTheory.volume 0 1 := by
  let c : ℝ → ℝ := fun x =>
    quadAltV x * Real.log (1 - x) / x
  have hc : IntervalIntegrable c MeasureTheory.volume 0 1 := by
    apply IntervalIntegrable.congr
      (f := fun x : ℝ => -(quadAltV x * H1 x / x)) ?_
      testVH1Div_intervalIntegrable.neg
    intro x _
    dsimp only [c]
    unfold H1
    ring
  have hfac : ContinuousOn (fun x : ℝ => x / (1 + x))
      (Set.uIcc (0 : ℝ) 1) := by
    rw [Set.uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)]
    intro x hx
    exact continuousAt_id.div (continuousAt_const.add continuousAt_id)
      (by linarith [hx.1] : 1 + x ≠ 0) |>.continuousWithinAt
  have hprod := hc.continuousOn_mul hfac
  apply IntervalIntegrable.congr
    (f := fun x : ℝ => (x / (1 + x)) * c x) ?_ hprod
  intro x hx
  have hx' : x ∈ Ioc (0 : ℝ) 1 := by
    simpa [Set.uIoc_of_le (by norm_num : (0 : ℝ) ≤ 1)] using hx
  have hxne : x ≠ 0 := ne_of_gt hx'.1
  have hplus : 1 + x ≠ 0 := by linarith [hx'.1]
  dsimp only [c]
  field_simp [hxne, hplus]

private theorem testVLogPlusMinus_intervalIntegrable :
    IntervalIntegrable
      (fun x : ℝ => quadAltV x * Real.log (1 + x) / (1 - x))
      MeasureTheory.volume 0 1 := by
  have hcont : ContinuousOn
      (fun x : ℝ => quadAltV x * Real.log (1 + x) / (1 - x))
      (Set.Ioo (0 : ℝ) 1) := by
    intro x hx
    exact ((quadAltV_continuousAt hx.1 hx.2).mul
      (ContinuousAt.log (by fun_prop)
        (by linarith [hx.1] : 1 + x ≠ 0))).div
      (by fun_prop) (by linarith [hx.2] : 1 - x ≠ 0) |>.continuousWithinAt
  have hlogDeriv : HasDerivAt (fun x : ℝ => Real.log (1 + x)) 1 0 := by
    have hinner : HasDerivAt (fun x : ℝ => 1 + x) 1 0 := by
      convert (hasDerivAt_const (0 : ℝ) 1).add (hasDerivAt_id 0) using 1 <;>
        norm_num
    have hlog := Real.hasDerivAt_log
      (show 1 + (0 : ℝ) ≠ 0 by norm_num)
    have hcomp := HasDerivAt.comp (h := fun x : ℝ => 1 + x)
      (0 : ℝ) hlog hinner
    simpa using hcomp
  have hlogSlope : Tendsto (fun x : ℝ => Real.log (1 + x) / x)
      (𝓝[>] (0 : ℝ)) (𝓝 1) :=
    (slope_tendsto_of_hasDerivAt_zero _ 1 hlogDeriv (by norm_num)).mono_left
      (nhdsWithin_mono _ (fun x hx => ne_of_gt hx))
  have hden0 : Tendsto (fun x : ℝ => 1 - x)
      (𝓝[>] (0 : ℝ)) (𝓝 1) := by
    have hid : Tendsto (fun x : ℝ => x) (𝓝[>] (0 : ℝ)) (𝓝 0) :=
      tendsto_id.mono_left nhdsWithin_le_nhds
    simpa using tendsto_const_nhds.sub hid
  have hzero : Tendsto
      (fun x : ℝ => quadAltV x * Real.log (1 + x) / (1 - x))
      (𝓝[>] (0 : ℝ)) (𝓝 0) := by
    have h := (quadAltV_mul_self_tendsto.mul hlogSlope).div hden0
      (by norm_num)
    norm_num at h
    refine h.congr' ?_
    filter_upwards [self_mem_nhdsWithin] with x hx
    have hxne : x ≠ 0 := ne_of_gt hx
    simp only [Pi.div_apply]
    field_simp [hxne]
  have hlogOne : Tendsto (fun x : ℝ => Real.log (1 + x))
      (𝓝[<] (1 : ℝ)) (𝓝 (Real.log 2)) := by
    have hc : ContinuousAt (fun x : ℝ => Real.log (1 + x)) 1 := by
      exact ContinuousAt.log (by fun_prop) (by norm_num)
    convert hc.tendsto.mono_left nhdsWithin_le_nhds using 1 <;> norm_num
  have hone : Tendsto
      (fun x : ℝ => quadAltV x * Real.log (1 + x) / (1 - x))
      (𝓝[<] (1 : ℝ)) (𝓝 0) := by
    have h := (quadAltV_slope_tendsto_one.mul hlogOne).neg
    simp only [zero_mul, neg_zero] at h
    refine h.congr' ?_
    filter_upwards [self_mem_nhdsWithin] with x hx
    have hne : x - 1 ≠ 0 := sub_ne_zero.mpr (ne_of_lt hx)
    have hden : 1 - x ≠ 0 := sub_ne_zero.mpr (ne_of_gt hx)
    field_simp [hne, hden]
    ring
  exact intervalIntegrable_of_continuousOn_Ioo_of_tendsto
    (by norm_num) hcont hzero hone

private noncomputable def testMixedLogPrimitive (x : ℝ) : ℝ :=
  quadAltV x * Real.log (1 - x) * Real.log (1 + x)

private theorem testMixedLogPrimitive_hasDerivAt
    {x : ℝ} (hx0 : 0 < x) (hx1 : x < 1) :
    HasDerivAt testMixedLogPrimitive
      (Real.log x * Real.log (1 - x) * Real.log (1 + x) /
          (x * (1 + x)) -
        quadAltV x * Real.log (1 + x) / (1 - x) +
        quadAltV x * Real.log (1 - x) / (1 + x)) x := by
  have hminus : 0 < 1 - x := by linarith
  have hplus : 0 < 1 + x := by linarith
  have hm : HasDerivAt (fun y : ℝ => Real.log (1 - y))
      (-1 / (1 - x)) x := by
    have hinner : HasDerivAt (fun y : ℝ => 1 - y) (-1) x := by
      convert (hasDerivAt_const x (1 : ℝ)).sub (hasDerivAt_id x) using 1 <;>
        norm_num
    have hlog := Real.hasDerivAt_log (ne_of_gt hminus)
    convert HasDerivAt.comp (h := fun y : ℝ => 1 - y) x hlog hinner using 1
    field_simp [ne_of_gt hminus]
  have hp : HasDerivAt (fun y : ℝ => Real.log (1 + y))
      (1 / (1 + x)) x := by
    have hinner : HasDerivAt (fun y : ℝ => 1 + y) 1 x := by
      convert (hasDerivAt_const x (1 : ℝ)).add (hasDerivAt_id x) using 1 <;>
        norm_num
    have hlog := Real.hasDerivAt_log (ne_of_gt hplus)
    convert HasDerivAt.comp (h := fun y : ℝ => 1 + y) x hlog hinner using 1
    field_simp [ne_of_gt hplus]
  have hprod := ((quadAltV_hasDerivAt hx0 hx1).mul hm).mul hp
  unfold testMixedLogPrimitive
  convert hprod using 1
  simp only [Pi.mul_apply]
  field_simp [ne_of_gt hx0, ne_of_gt hminus, ne_of_gt hplus]
  ring

private theorem testMixedLogPrimitive_tendsto_zero :
    Tendsto testMixedLogPrimitive (𝓝[>] (0 : ℝ)) (𝓝 0) := by
  have hminus : Tendsto (fun x : ℝ => Real.log (1 - x))
      (𝓝[>] (0 : ℝ)) (𝓝 0) := by
    have hc : ContinuousAt (fun x : ℝ => Real.log (1 - x)) 0 := by
      exact ContinuousAt.log (by fun_prop) (by norm_num)
    simpa using hc.tendsto.mono_left nhdsWithin_le_nhds
  have hlogDeriv : HasDerivAt (fun x : ℝ => Real.log (1 + x)) 1 0 := by
    have hinner : HasDerivAt (fun x : ℝ => 1 + x) 1 0 := by
      convert (hasDerivAt_const (0 : ℝ) 1).add (hasDerivAt_id 0) using 1 <;>
        norm_num
    have hlog := Real.hasDerivAt_log
      (show 1 + (0 : ℝ) ≠ 0 by norm_num)
    have hcomp := HasDerivAt.comp (h := fun x : ℝ => 1 + x)
      (0 : ℝ) hlog hinner
    simpa using hcomp
  have hplusSlope : Tendsto (fun x : ℝ => Real.log (1 + x) / x)
      (𝓝[>] (0 : ℝ)) (𝓝 1) :=
    (slope_tendsto_of_hasDerivAt_zero _ 1 hlogDeriv (by norm_num)).mono_left
      (nhdsWithin_mono _ (fun x hx => ne_of_gt hx))
  have h := (quadAltV_mul_self_tendsto.mul hminus).mul hplusSlope
  simp only [zero_mul] at h
  refine h.congr' ?_
  filter_upwards [self_mem_nhdsWithin] with x hx
  have hxne : x ≠ 0 := ne_of_gt hx
  unfold testMixedLogPrimitive
  field_simp [hxne]

private theorem testMixedLogPrimitive_tendsto_one :
    Tendsto testMixedLogPrimitive (𝓝[<] (1 : ℝ)) (𝓝 0) := by
  have hminus : Tendsto (fun x : ℝ => (x - 1) * Real.log (1 - x))
      (𝓝[<] (1 : ℝ)) (𝓝 0) := by
    have h := oneSub_log_tendsto.neg
    simp only [neg_zero] at h
    refine h.congr' ?_
    filter_upwards with x
    ring
  have hplus : Tendsto (fun x : ℝ => Real.log (1 + x))
      (𝓝[<] (1 : ℝ)) (𝓝 (Real.log 2)) := by
    have hc : ContinuousAt (fun x : ℝ => Real.log (1 + x)) 1 := by
      exact ContinuousAt.log (by fun_prop) (by norm_num)
    convert hc.tendsto.mono_left nhdsWithin_le_nhds using 1 <;> norm_num
  have h := (quadAltV_slope_tendsto_one.mul hminus).mul hplus
  simp only [zero_mul] at h
  refine h.congr' ?_
  filter_upwards [self_mem_nhdsWithin] with x hx
  have hne : x - 1 ≠ 0 := sub_ne_zero.mpr (ne_of_lt hx)
  unfold testMixedLogPrimitive
  field_simp [hne]

private theorem testMixedLogIntegralRelation :
    (∫ x : ℝ in 0..1,
      quadAltV x * Real.log (1 - x) / (1 + x)) -
        ∫ x : ℝ in 0..1,
          quadAltV x * Real.log (1 + x) / (1 - x) =
      -(∫ x : ℝ in 0..1,
        Real.log x * Real.log (1 - x) * Real.log (1 + x) /
          (x * (1 + x))) := by
  have hderivInt :=
    (quadAltMixedLogIntervalIntegrable24.sub
      testVLogPlusMinus_intervalIntegrable).add
      testVLogMinusPlus_intervalIntegrable
  have hFTC := intervalIntegral.integral_eq_sub_of_hasDerivAt_of_tendsto
    (a := (0 : ℝ)) (b := 1) (f := testMixedLogPrimitive)
    (f' := fun x : ℝ =>
      Real.log x * Real.log (1 - x) * Real.log (1 + x) /
          (x * (1 + x)) -
        quadAltV x * Real.log (1 + x) / (1 - x) +
        quadAltV x * Real.log (1 - x) / (1 + x))
    (by norm_num)
    (fun x hx => testMixedLogPrimitive_hasDerivAt hx.1 hx.2)
    hderivInt testMixedLogPrimitive_tendsto_zero
    testMixedLogPrimitive_tendsto_one
  rw [intervalIntegral.integral_add
      (quadAltMixedLogIntervalIntegrable24.sub
        testVLogPlusMinus_intervalIntegrable)
      testVLogMinusPlus_intervalIntegrable,
    intervalIntegral.integral_sub quadAltMixedLogIntervalIntegrable24
      testVLogPlusMinus_intervalIntegrable] at hFTC
  norm_num at hFTC
  linarith

theorem quadAltI21_eq :
    I21 = -(3 / 2) * quadAltK
      - 3 * Real.log 2 ^ 2 * (Real.pi ^ 2 / 6)
      + (7 / 4) * Real.log 2 * zeta3_24
      + (3 / 10) * (Real.pi ^ 2 / 6) ^ 2 := by
  have h21 := testI21_substitution
  have h12 := testI12_substitution
  have hmix := testMixedLogIntegralRelation
  have hrow : I21 = 2 * I22 - I12 -
      2 * ∫ x : ℝ in 0..1,
        Real.log x * Real.log (1 - x) * Real.log (1 + x) /
          (x * (1 + x)) := by
    linarith
  rw [quadAltI22_eq, quadAltI12_eq, quadAltMixedLogIntegral24] at hrow
  linarith

/-- The missing P2.4 hypothesis: the outer-alternating quadratic Euler sum
evaluates to `alternatingQuadraticEulerValue24`. -/
theorem alternatingQuadraticEulerTerm24_hasSum :
    HasSum alternatingQuadraticEulerTerm24 alternatingQuadraticEulerValue24 := by
  have htsum : (∑' n : ℕ, alternatingQuadraticEulerTerm24 n) =
      alternatingQuadraticEulerValue24 := by
    calc
      (∑' n : ℕ, alternatingQuadraticEulerTerm24 n) =
          ∫ x : ℝ in (0 : ℝ)..1,
            (-Real.log x) / x * quadAltQclosed (-x) :=
        quadAlt_tsum_eq_coeff_integral
      _ = -2 * I10 - 2 * I11 + 2 * I12 + 4 * I20 + 6 * I21 - 5 * I22 :=
        quadAltCoeffIntegral_eq_six
      _ = sixIntegralCombination := rfl
      _ = bridgeValue := quadAltSixIntegral_eq_bridgeValue
        quadAltK_eq quadAltI10_eq quadAltI11_eq quadAltI12_eq
        quadAltI20_eq quadAltI21_eq quadAltI22_eq
      _ = alternatingQuadraticEulerValue24 :=
        bridgeValue_eq_alternatingQuadraticEulerValue24
  rw [← htsum]
  exact summable_quadAlt.hasSum

end RamanujanChallenge.P24QuadAlt
