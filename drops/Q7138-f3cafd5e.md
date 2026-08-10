ANSWER Q7138 f3cafd5e

# Valuation-level primitive/nonprimitive decomposition audit

## Verdict

The proposed exact Mobius inversion on return levels does **not** hold in the required form for valuations of the adjacent fully-deflated content. The obstruction is that the content is an intersection-multiplicity invariant of a zero-dimensional algebra, not a sum over geometric roots. A root which lies on several return loci contributes with local Artin length, and inclusion-exclusion on the underlying sets does not preserve that length.

A useful replacement is a one-sided valuation inequality using saturated ideals and Fitting ideals.

## 1. Saturated algebra formulation

For a gap triple `(a,b,c)` write

```
F=N_a^o(x), G=N_b^o(x+a), J=N_c^o(x+a+b).
```

Let `I=(F,G,J)` in `Z[x]`. Let `E` be the product of all forbidden factors:

```
E = (center factors) * product_{d in D(a,b,c)} N_d(x)
```

where `D(a,b,c)` is the set of intermediate return gaps. The primitive algebra is the zero-dimensional algebra

```
A_prim = (Q[x]/I) : E^infinity.
```

The fully deflated content controls the Fitting ideal of this saturated algebra only. It does not canonically split into Fitting ideals of all smaller return strata.

## 2. Why exact Mobius inversion fails

A concrete local counterexample:

Take the Artin algebra

```
A = k[t]/(t^4).
```

Let two return divisors correspond locally to `(t)` and `(t^2)`. The lengths are

```
length A/(t)=1,
length A/(t^2)=2,
length A/(t^3)=3.
```

A set-theoretic Mobius subtraction sees only the same closed point. It cannot recover the multiplicity contribution 4 from the local algebra. Therefore a formula of the form

```
v_p(C_fd)=primitive + sum(return-chain charges)
```
with exact coefficients determined only by the return poset is false in general.

The missing datum is the nilpotent filtration of each local component.

## 3. Surviving inequality

The correct statement is:

```
v_p(C_fd(a,b,c))
<= v_p(C_prim(a,b,c))
 + sum_{D} m_D(a,b,c) v_p(C_D)
 + v_p(E_content),
```

where:

* `C_prim` is the Fitting generator of the saturated algebra;
* `C_D` is the adjacent content attached to a smaller return chain;
* `m_D` is the maximum local intersection multiplicity with that return stratum.

This follows from the exact sequence induced by saturation:

```
0 -> I^sat/I -> Z[x]/I -> Z[x]/I^sat -> 0
```

and multiplicativity of zeroth Fitting ideals up to containment:

```
Fitt_0(M_extension) subset Fitt_0(M_1)Fitt_0(M_2).
```

Taking `p`-adic valuations gives the inequality above.

## 4. Consequence for H-mass

The desired `H^{3+o(1)}` compression is therefore conditional on proving two new inputs:

1. a uniform bound

```
sum_{a+b+c<=H} log C_prim(a,b,c) << H^{3+o(1)};
```

2. a return-poset multiplicity bound

```
sum_D m_D log C_D << H^{3+o(1)}.
```

The existing adjacent resultant carrier handles the first type of content, but it does not automatically bound the second because repeated returns can increase Artin length without increasing the number of geometric roots.

## 5. Practical recommendation

Do not replace the current AC^tr hypothesis by an exact primitive/nonprimitive decomposition. The safe publication statement is:

* fully-deflated adjacent content removes the center degeneracy;
* residual skipped chains are candidates for a secondary charging argument;
* an exact valuation Mobius inversion is unavailable without extra local-length control;
* the strongest unconditional replacement is the Fitting-ideal inequality above.

Any future compression theorem should explicitly introduce local intersection multiplicity weights, not only the return-time poset.
