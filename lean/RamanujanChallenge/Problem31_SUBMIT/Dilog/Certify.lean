import Ramanujan31.Dilog.RealBounds

/-!
# A certified evaluation of `Li₂`, end to end

This file demonstrates that the machinery of `Dilog/RealBounds.lean` actually
closes a numerical goal inside Lean, with no `sorry`, no `native_decide`, and no
floating point anywhere: the enclosure is a rational inequality and `norm_num`
discharges it.

## What this establishes, and what it does not

It establishes that "certified evaluation of `Li₂` in Lean" is not a research
problem — the analytic content is the tail bound, which is proved, and the rest
is rational arithmetic.

It does **not** yet establish that the *particular* eight evaluations Problem 3.1
needs are tractable this way.  Here is the obstruction, stated plainly, because
it determines what has to be built next.

## The tractability finding

The four tetrahedron shapes have `3`-cycle representatives in `(0,1)`; at the
first endpoint they are

    0.86365088,  0.85328205,  0.98210210,  0.83144898.

For the first, second and fourth, a few hundred terms give `10⁻¹²` and exact
rational arithmetic is fine.  For `0.9821021` it is not: the tail bound
`x^{N+1}/((N+1)²(1-x))` needs `N ≈ 1000`, and the exact rational partial sum then
carries a denominator divisible by `lcm(1,…,1000)²`, a number with roughly `870`
digits before the powers of the argument are even accounted for.  `norm_num` on
that is not a good plan.

The fix is not more arithmetic, it is **Euler's reflection formula**

    Li₂(x) + Li₂(1-x) = π²/6 - log x · log(1-x),

which sends `0.9821021` to `0.0178979`, where **eight terms** suffice.

So the reflection formula is the keystone twice over: it is what makes
`D(1-z) = -D(z)` — hence `D(1/(1-z)) = D(z)`, the functional equation the
Bloch–Wigner argument runs on — and it is also what makes the numerical
certificate arithmetically feasible.  One lemma, both gaps.

`Li₂(1) = π²/6` (`Complex.dilog_one`, from Mathlib's `hasSum_zeta_two`) is the
anchor that reflection is proved against, and it is already done.
-/

namespace Real

/-- **A certified evaluation.**  `Li₂(1/10)` agrees with its `12`-term partial sum
to better than `10⁻¹²`.

The proof is: apply the tail bound, then `norm_num`.  Nothing else. -/
example : |Real.dilog (1/10) - Real.dilogPartial (1/10) 12| ≤ 1/10 ^ 12 := by
  refine Real.dilog_enclosure_width (by norm_num) (by norm_num) 12 ?_
  norm_num [Real.dilogTerm]

/-- The partial sum is an explicit rational — nothing is hidden inside a
`tsum`. -/
example : Real.dilogPartial (1/10) 12
    = 1/10 + (1/10) ^ 2 / 4 + (1/10) ^ 3 / 9 + (1/10) ^ 4 / 16 + (1/10) ^ 5 / 25
      + (1/10) ^ 6 / 36 + (1/10) ^ 7 / 49 + (1/10) ^ 8 / 64 + (1/10) ^ 9 / 81
      + (1/10) ^ 10 / 100 + (1/10) ^ 11 / 121 + (1/10) ^ 12 / 144 := by
  simp [Real.dilogPartial, Finset.sum_range_succ, Real.dilogTerm]
  norm_num

/-- The same at an argument close to the ones the problem actually produces
(`0.86` is the largest of the three benign representatives).  Here `120` terms
give better than `10⁻⁸`, which is already past the `3·10⁻⁸` the reconstruction
needs. -/
example : |Real.dilog (86/100) - Real.dilogPartial (86/100) 120| ≤ 1/10 ^ 8 := by
  refine Real.dilog_enclosure_width (by norm_num) (by norm_num) 120 ?_
  norm_num

end Real
