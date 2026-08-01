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

end RamanujanChallenge.P24QuadAlt

