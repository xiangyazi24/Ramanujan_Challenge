# Codex Spec: Prove aperyB_recurrence_int (WZ Identity)

## Target

Prove in Lean 4 (with Mathlib) the following theorem from
`lean/RamanujanChallenge/Problem32/Main.lean`:

```lean
theorem aperyB_recurrence_int (n : ℕ) (hn : n ≥ 1) :
    ((n : ℤ) + 1) ^ 3 * aperyB (n + 1) =
      (34 * (n : ℤ) ^ 3 + 51 * (n : ℤ) ^ 2 + 27 * (n : ℤ) + 5) * aperyB n -
      (n : ℤ) ^ 3 * aperyB (n - 1)
```

where `aperyB` is defined in `lean/RamanujanChallenge/Problem32/AperyDef.lean`:

```lean
def aperyB (n : ℕ) : ℤ :=
  ∑ k ∈ Finset.range (n + 1),
    (Nat.choose n k : ℤ) ^ 2 * (Nat.choose (n + k) k : ℤ) ^ 2
```

## Mathematical Background

This is the Apéry recurrence for the numbers b_n = Σ_{k=0}^n C(n,k)²C(n+k,k)².
It is a standard hypergeometric identity proved by Zeilberger's algorithm.

## WZ Certificate (computed by Sage/Maxima)

The Zeilberger algorithm gives the certificate:

```
R(n,k) = 4k⁴(2n+3)(4n²+12n-2k²+3k+8) / ((n-k+1)²(n-k+2)²)
```

with recurrence coefficients (in the shifted convention n, n+1, n+2):

```
a₀(n) = -(n+1)³
a₁(n) = (2n+3)(17n²+51n+39)
a₂(n) = -(n+2)³
```

The recurrence a₀ b_n + a₁ b_{n+1} + a₂ b_{n+2} = 0 is equivalent to
our target: (m+1)³ b_{m+1} = P(m) b_m - m³ b_{m-1} with m = n.

## Key Polynomial Identity (VERIFIED)

Define f(n,k) = C(n,k)² C(n+k,k)². For 0 ≤ k < n (where f(n,k) > 0):

The WZ equation divided by f(n,k) and multiplied by (n-k+1)²(n-k+2)² gives:

```
LHS(n,k) = -(n+1)³(n+1-k)²(n+2-k)²
          + (2n+3)(17n²+51n+39)(n+1+k)²(n+2-k)²
          - (n+2)³(n+2+k)²(n+1+k)²

RHS(n,k) = 4(2n+3)(4n²+12n-2(k+1)²+3(k+1)+8)(n+k+1)²(n-k+2)²
          - 4k⁴(2n+3)(4n²+12n-2k²+3k+8)
```

**LHS(n,k) = RHS(n,k) as polynomials in ℤ[n,k].** (Verified by sympy: expand(LHS-RHS) = 0.)

This identity can be verified by `ring` in Lean 4.

## Proof Strategy

### Option A: Direct sum manipulation (RECOMMENDED)

Prove the recurrence directly using the Finset.sum API:

1. Express LHS - RHS of the recurrence as a single sum:
   (n+1)³ Σ_{k≤n+1} f(n+1,k) - P(n) Σ_{k≤n} f(n,k) + n³ Σ_{k≤n-1} f(n-1,k)

2. Reindex and combine into a single sum over k = 0 to n+1.

3. For each summand, show the combined expression equals a telescoping difference.

4. The telescope collapses to boundary terms which are 0.

### Option B: WZ telescoping (MORE PRINCIPLED)

1. **Define the anti-difference** G(n,k) for 0 ≤ k ≤ n:
   G(n,k) = R(n,k) · f(n,k)
   = 4k⁴(2n+3)(4n²+12n-2k²+3k+8) · C(n,k)² · C(n+k,k)² / ((n-k+1)²(n-k+2)²)

   Note: G(n,0) = 0 since k⁴|_{k=0} = 0.

2. **Prove the WZ equation for 0 ≤ k ≤ n-1:**
   a₀ f(n,k) + a₁ f(n+1,k) + a₂ f(n+2,k) = G(n,k+1) - G(n,k)

   This requires: after dividing by f(n,k) > 0, the equation becomes the
   polynomial identity LHS(n,k) = RHS(n,k) verified above.

   In Lean: express the identity in multiplicative form (no divisions), as:
   [a₀ + a₁ · f(n+1,k)/f(n,k) + a₂ · f(n+2,k)/f(n,k)] · (n-k+1)²(n-k+2)²
   = [R(n,k+1) · f(n,k+1)/f(n,k) - R(n,k)] · (n-k+1)²(n-k+2)²

   Both sides are polynomials in n,k multiplied by f(n,k). Use ring after
   clearing the binomial coefficients via ratio identities:
   - f(n+1,k)/f(n,k) = (n+1+k)²/(n+1-k)²
   - f(n,k+1)/f(n,k) = (n-k)²(n+k+1)²/(k+1)⁴

3. **Telescope the sum from k=0 to n-1:**
   Σ_{k=0}^{n-1} [a₀ f(n,k) + a₁ f(n+1,k) + a₂ f(n+2,k)] = G(n,n) - G(n,0) = G(n,n)

4. **Handle boundary terms (k = n, n+1, n+2):**
   Complete the partial sums to the full sums b_n, b_{n+1}, b_{n+2}.
   The remaining terms are:
   - a₀ · f(n,n)
   - a₁ · [f(n+1,n) + f(n+1,n+1)]
   - a₂ · [f(n+2,n) + f(n+2,n+1) + f(n+2,n+2)]

5. **Show G(n,n) + boundary_terms = 0.**
   This is a concrete computation involving C(n,n)=1, C(2n,n), C(2n+1,n),
   C(2n+2,n), C(2n+2,n+1), C(2n+4,n+2), etc.
   This final identity can be verified by ring after expressing everything
   in terms of C(2n,n) and using binomial coefficient identities.

### Option C: Equality with recurrence-defined sequence (CLEANEST for Lean)

The file already defines `aperyBQ : ℕ → ℚ` by the recurrence.
If we can prove `(aperyB n : ℚ) = aperyBQ n` for all n, the recurrence follows.

Proof: by strong induction. Base: n=0,1 (already proved). Step requires the
WZ identity to show that the closed form at n+2 matches the recurrence value.
THIS IS CIRCULAR — the WZ identity IS what we're proving.

## Hard Rules

- No `sorry`, no custom `axiom`.
- No `native_decide` on unbounded terms (OK for specific small values).
- All new declarations must compile in `lake build`.
- `#print axioms` must show only `{propext, Classical.choice, Quot.sound}`.

## Verification

`lake build` must succeed. The `sorry` at line 48 of Main.lean must be replaced.

## Stall Protocol

If stuck, deliver what compiles + exact goal state at the stuck point +
which lemma/API is missing. Do NOT fake with `sorry`.

## Key Mathlib APIs likely needed

- `Finset.sum_range_succ`
- `Finset.sum_congr`
- `Nat.choose_succ_succ`
- `Nat.choose_self`
- `Nat.choose_zero_right`
- `Int.mul_comm`, `Int.mul_assoc`
- `ring`, `omega`, `norm_num`
- `Finset.sum_sub_distrib` or similar for combining sums

## Notes

The WZ approach is the most principled, but the boundary term handling is the
hardest part. If the boundary terms prove too difficult, an alternative is to
verify the identity for n=1 through n=50 by `norm_num` (covering the base cases)
and then use the WZ equation only for the inductive step on the inner sum.
However, the inductive step still requires the WZ polynomial identity.

The polynomial identity LHS(n,k) = RHS(n,k) has total degree ~13 in (n,k).
Lean's `ring` tactic should handle this.
