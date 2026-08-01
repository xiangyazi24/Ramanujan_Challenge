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

/-! ## Layer D: integration-by-parts objects (Q6047 §4-5)

`W(x) = Z2 − 2 Li2(u) − log(u)²` with `u = x/(1+x)` (Landen form, (4.7));
after the Möbius map `t = 2x/(1+x)` this becomes `W0(t)`. -/

/-- `W0(t) = Z2 − 2 Li2(t/2) − log(t/2)²` (Q6047 (5.2)). -/
def W0 (t : ℝ) : ℝ :=
  Real.pi ^ 2 / 6 - 2 * dilog (t / 2) - Real.log (t / 2) ^ 2

/-- `H1(t) = -log(1-t)`. -/
def H1 (t : ℝ) : ℝ := -Real.log (1 - t)

/-- `H2(t) = -log(1-t/2)`. -/
def H2 (t : ℝ) : ℝ := -Real.log (1 - t / 2)

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

/-! ## The target theorem (open) -/

/-- Summability of the outer-alternating quadratic terms (already proved in
`Problem24`). -/
theorem summable_quadAlt : Summable alternatingQuadraticEulerTerm24 :=
  summable_alternatingQuadraticEulerTerm24

/-- The missing P2.4 hypothesis: the outer-alternating quadratic Euler sum
evaluates to `alternatingQuadraticEulerValue24`.  This is the target of the
Q6047 six-integral certificate; the endpoint lemmas I10..I22 remain to be
filled in. -/
theorem alternatingQuadraticEulerTerm24_hasSum :
    HasSum alternatingQuadraticEulerTerm24 alternatingQuadraticEulerValue24 := by
  -- route (Q6047 Layer F):
  --   rw [← summable_quadAlt.hasSum.tsum_eq]
  --   tsum = ∫₀¹ (-log x)/x * Q(-x) dx   (coefficient integration)
  --        = sixIntegralCombination        (IBP + Mobius)
  --        = bridgeValue                   (endpoint evaluations)
  --        = alternatingQuadraticEulerValue24
  sorry

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

end RamanujanChallenge.P24QuadAlt

