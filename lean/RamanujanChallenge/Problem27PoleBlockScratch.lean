import RamanujanChallenge.Problem27BarnesNormalization
import Mathlib.Analysis.Calculus.ParametricIntegral
import Mathlib.Analysis.PSeries
import Mathlib.Analysis.SpecialFunctions.Gamma.Basic
import Mathlib.MeasureTheory.Integral.DominatedConvergence
import Mathlib.MeasureTheory.Integral.ExpDecay
import Mathlib.MeasureTheory.Integral.Prod
import Mathlib.NumberTheory.ZetaValues
import Mathlib.Topology.Algebra.InfiniteSum.NatInt

open Filter Set MeasureTheory

noncomputable section

namespace RamanujanChallenge.P27

#check MeasureTheory.integral_integral_swap
#check MeasureTheory.integral_tsum_of_summable_integral_norm
#check MeasureTheory.hasSum_integral_of_summable_integral_norm
#check MeasureTheory.Integrable.mul_prod
#check MeasureTheory.Integrable.smul_prod
#check MeasureTheory.Integrable.mono'
#check integrableOn_exp_mul_complex_Ioi
#check integral_exp_mul_complex_Ioi
#check MeasureTheory.integral_Ioi_of_hasDerivAt_of_tendsto
#check MeasureTheory.integral_Ioi_mul_deriv_eq_deriv_mul
#check Real.GammaIntegral_convergent
#check Real.integral_rpow_mul_exp_neg_mul_Ioi
#check Real.Gamma_nat_eq_factorial
#check hasSum_geometric_of_norm_lt_one
#check summable_geometric_of_norm_lt_one
#check Summable.sum_add_tsum_nat_add
#check sum_add_tsum_nat_add
#check hasSum_nat_add_iff'
#check summable_nat_add_iff
#check hasSum_zeta_two
#check MeasureTheory.integral_const_mul
#check MeasureTheory.integral_mul_const
#check MeasureTheory.integral_smul
#check MeasureTheory.integral_indicator
#check MeasureTheory.IntegrableOn
#check MeasureTheory.integrableOn_Ioi_iff_integrableOn_Ici
#check integrableOn_Ici_iff_integrableOn_Ioi
#check Real.norm_exp
#check Complex.norm_exp
#check Complex.norm_real
#check Complex.exp_add
#check Complex.exp_mul_I
#check Complex.exp_neg_mul_I

end RamanujanChallenge.P27
