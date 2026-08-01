/-
  Tight Catalan bounds for the sign-pattern verification at N=1.
  Uses native_decide for efficient rational arithmetic.
-/
import RamanujanChallenge.Problem25

noncomputable section

namespace RamanujanChallenge.P25

open Finset

set_option maxHeartbeats 0 in
set_option maxRecDepth 2000 in
theorem catalan_tight_lower :
    (1590511050 : ℝ) / 1736437500 < catalanConstant := by
  calc
    (1590511050 : ℝ) / 1736437500 <
        ∑ i ∈ Finset.range (2 * 135), (-1 : ℝ) ^ i * catalanMagnitude i := by
      norm_num [catalanMagnitude, Finset.sum_range_succ]
    _ ≤ catalanConstant :=
      catalanMagnitude_antitone.alternating_series_le_tendsto
        catalanMagnitude_partialSum_tendsto 135

set_option maxHeartbeats 0 in
set_option maxRecDepth 2000 in
theorem catalan_tight_upper :
    catalanConstant < (21390206625 : ℝ) / 23352603750 := by
  calc
    catalanConstant ≤
        ∑ i ∈ Finset.range (2 * 242 + 1), (-1 : ℝ) ^ i * catalanMagnitude i :=
      catalanMagnitude_antitone.tendsto_le_alternating_series
        catalanMagnitude_partialSum_tendsto 242
    _ < (21390206625 : ℝ) / 23352603750 := by
      norm_num [catalanMagnitude, Finset.sum_range_succ]

end RamanujanChallenge.P25

end
