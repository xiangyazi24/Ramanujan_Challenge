import Ramanujan31.Dilog.FiveTerm

/-!
# Finite endpoint transformation certificates

This file replays finite linear combinations of Abel five-term and Euler
complement relations.  The algebraic hypotheses are deliberately explicit:
later endpoint modules discharge them in the appropriate number field.
-/

namespace Real

private theorem rogers_complement_rewrite
    {x y : ℝ} (hx0 : 0 < x) (hx1 : x < 1) (hy : y = 1 - x) :
    rogers x + rogers y = Real.pi ^ 2 / 6 := by
  rw [hy]
  exact rogers_add_rogers_one_sub hx0 hx1

private theorem rogers_five_term_rewrite
    {x y p q r : ℝ}
    (hx0 : 0 < x) (hx1 : x < 1) (hy0 : 0 < y) (hy1 : y < 1)
    (hp : x * y = p)
    (hq : x * (1 - y) / (1 - x * y) = q)
    (hr : y * (1 - x) / (1 - x * y) = r) :
    rogers x + rogers y = rogers p + rogers q + rogers r := by
  have h := rogers_five_term hx0 hx1 hy0 hy1
  rw [hq, hr, hp] at h
  exact h

/--
The 22-relation certificate reducing the three non-common beta chart
arguments to the first four standard `π/17` arguments.

The certificate is a rational linear combination of twelve Euler complement
relations and ten Abel five-term relations.  All cross-identifications are
listed as hypotheses so that they can be checked by polynomial normalization
at the algebraic endpoint.
-/
theorem beta_chart_to_standard_certificate
    {b c d s1 s2 s3 s4 : ℝ}
    {u10 u11 u12 u13 u14 u15 u17 u26 u27 u30 u32 u33 u34 u35 u36 u37
      u38 u39 u51 u54 u56 u69 u70 u71 u73 u74 : ℝ}
    (hunit :
      ∀ x ∈
        [b, c, d, s1, s2, s3, s4, u10, u11, u12, u13, u14, u15,
          u17, u27, u30, u32, u33, u36, u51, u39, u35, u54],
        0 < x ∧ x < 1)
    (hc10 : u10 = 1 - c)
    (hc11 : u11 = 1 - d)
    (hc12 : u12 = 1 - s1)
    (hc13 : u13 = 1 - s2)
    (hc14 : u14 = 1 - s3)
    (hc15 : u15 = 1 - s4)
    (hc17 : u37 = 1 - u17)
    (hc27 : u51 = 1 - u27)
    (hc30 : u54 = 1 - u30)
    (hc32 : u39 = 1 - u32)
    (hc33 : u56 = 1 - u33)
    (hc36 : u38 = 1 - u36)
    (hf1p : b * b = s4)
    (hf1q : b * (1 - b) / (1 - b * b) = u17)
    (hf1r : b * (1 - b) / (1 - b * b) = u17)
    (hf2p : s4 * u11 = c)
    (hf2q : s4 * (1 - u11) / (1 - s4 * u11) = u26)
    (hf2r : u11 * (1 - s4) / (1 - s4 * u11) = u27)
    (hf3p : u10 * u11 = u33)
    (hf3q : u10 * (1 - u11) / (1 - u10 * u11) = u34)
    (hf3r : u11 * (1 - u10) / (1 - u10 * u11) = u35)
    (hf4p : u12 * u13 = u36)
    (hf4q : u12 * (1 - u13) / (1 - u12 * u13) = s1)
    (hf4r : u13 * (1 - u12) / (1 - u12 * u13) = u36)
    (hf5p : u13 * u14 = u37)
    (hf5q : u13 * (1 - u14) / (1 - u13 * u14) = u38)
    (hf5r : u14 * (1 - u13) / (1 - u13 * u14) = u39)
    (hf6p : u14 * u15 = u34)
    (hf6q : u14 * (1 - u15) / (1 - u14 * u15) = u32)
    (hf6r : u15 * (1 - u14) / (1 - u14 * u15) = u30)
    (hf7p : b * u51 = u69)
    (hf7q : b * (1 - u51) / (1 - b * u51) = u70)
    (hf7r : u51 * (1 - b) / (1 - b * u51) = u71)
    (hf8p : c * u39 = u73)
    (hf8q : c * (1 - u39) / (1 - c * u39) = u26)
    (hf8r : u39 * (1 - c) / (1 - c * u39) = u74)
    (hf9p : s1 * u35 = u73)
    (hf9q : s1 * (1 - u35) / (1 - s1 * u35) = u71)
    (hf9r : u35 * (1 - s1) / (1 - s1 * u35) = u69)
    (hf10p : u11 * u54 = u56)
    (hf10q : u11 * (1 - u54) / (1 - u11 * u54) = u74)
    (hf10r : u54 * (1 - u11) / (1 - u11 * u54) = u70) :
    rogers b + rogers c + rogers d =
      rogers s1 + rogers s2 + rogers s3 + rogers s4 := by
  have hb := hunit b (by simp)
  have hc := hunit c (by simp)
  have hd := hunit d (by simp)
  have hs1 := hunit s1 (by simp)
  have hs2 := hunit s2 (by simp)
  have hs3 := hunit s3 (by simp)
  have hs4 := hunit s4 (by simp)
  have hu10 := hunit u10 (by simp)
  have hu11 := hunit u11 (by simp)
  have hu12 := hunit u12 (by simp)
  have hu13 := hunit u13 (by simp)
  have hu14 := hunit u14 (by simp)
  have hu15 := hunit u15 (by simp)
  have hu17 := hunit u17 (by simp)
  have hu27 := hunit u27 (by simp)
  have hu30 := hunit u30 (by simp)
  have hu32 := hunit u32 (by simp)
  have hu33 := hunit u33 (by simp)
  have hu36 := hunit u36 (by simp)
  have hu51 := hunit u51 (by simp)
  have hu39 := hunit u39 (by simp)
  have hu35 := hunit u35 (by simp)
  have hu54 := hunit u54 (by simp)

  have c1 := rogers_complement_rewrite hc.1 hc.2 hc10
  have c2 := rogers_complement_rewrite hd.1 hd.2 hc11
  have c3 := rogers_complement_rewrite hs1.1 hs1.2 hc12
  have c4 := rogers_complement_rewrite hs2.1 hs2.2 hc13
  have c5 := rogers_complement_rewrite hs3.1 hs3.2 hc14
  have c6 := rogers_complement_rewrite hs4.1 hs4.2 hc15
  have c7 := rogers_complement_rewrite hu17.1 hu17.2 hc17
  have c8 := rogers_complement_rewrite hu27.1 hu27.2 hc27
  have c9 := rogers_complement_rewrite hu30.1 hu30.2 hc30
  have c10 := rogers_complement_rewrite hu32.1 hu32.2 hc32
  have c11 := rogers_complement_rewrite hu33.1 hu33.2 hc33
  have c12 := rogers_complement_rewrite hu36.1 hu36.2 hc36

  have f1 := rogers_five_term_rewrite
    hb.1 hb.2 hb.1 hb.2 hf1p hf1q hf1r
  have f2 := rogers_five_term_rewrite
    hs4.1 hs4.2 hu11.1 hu11.2 hf2p hf2q hf2r
  have f3 := rogers_five_term_rewrite
    hu10.1 hu10.2 hu11.1 hu11.2 hf3p hf3q hf3r
  have f4 := rogers_five_term_rewrite
    hu12.1 hu12.2 hu13.1 hu13.2 hf4p hf4q hf4r
  have f5 := rogers_five_term_rewrite
    hu13.1 hu13.2 hu14.1 hu14.2 hf5p hf5q hf5r
  have f6 := rogers_five_term_rewrite
    hu14.1 hu14.2 hu15.1 hu15.2 hf6p hf6q hf6r
  have f7 := rogers_five_term_rewrite
    hb.1 hb.2 hu51.1 hu51.2 hf7p hf7q hf7r
  have f8 := rogers_five_term_rewrite
    hc.1 hc.2 hu39.1 hu39.2 hf8p hf8q hf8r
  have f9 := rogers_five_term_rewrite
    hs1.1 hs1.2 hu35.1 hu35.2 hf9p hf9q hf9r
  have f10 := rogers_five_term_rewrite
    hu11.1 hu11.2 hu54.1 hu54.2 hf10p hf10q hf10r

  linarith [c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12,
    f1, f2, f3, f4, f5, f6, f7, f8, f9, f10]

end Real
