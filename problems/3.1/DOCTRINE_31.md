# P3.1 doctrine — closing the knot-regulator conjecture

**Goal (one sentence).** Prove
`int_alpha^beta (log x dy/y - log y dx/x) = 4 pi^2/85` for the `7_2` A-polynomial
arc, unconditionally.

## State at doctrine time

Reduction chain, all established:

1. Khoi's variation formula, correct normalization:
   `I = (GV(rho_alpha) - GV(rho_beta))/4`.
2. Beta endpoint exact: `S^3_{-1}(7_2) = Sigma(2,3,17)`, maximal Fuchsian,
   `GV(rho_beta) = 242 pi^2/51` (Brooks-Goldman).
3. Hence the conjecture `<=>` `GV(rho_alpha) = 74 pi^2/15`.
4. Neumann-Zagier: `I = -Delta R`, `R = sum_j Rhat(z_j)` over the four
   tetrahedron shapes of the `7_2` chart.
5. **Numerically** `Re[Delta R]/pi^2 = -4/85` to 301 digits (independent
   re-derivation this session; continued fraction terminates at `-4/85` with a
   ~300-digit next partial quotient).

**The only gap:** proving `Re[Delta R]/pi^2` is RATIONAL, with a denominator
bound. Then 301 digits pins it.

## The torsion mechanism (established this session)

Both endpoint minimal polynomials are palindromic:

```
f_alpha = a^12 - 3a^11 + 4a^10 - 5a^9 + 6a^8 - 7a^7 + 7a^6 - 7a^5
          + 6a^4 - 5a^3 + 4a^2 - 3a + 1          (irreducible, sig (2,5))
f_beta  = b^16 - 7b^15 + 22b^14 - 48b^13 + 87b^12 - 133b^11 + 178b^10
          - 211b^9 + 223b^8 - 211b^7 + 178b^6 - 133b^5 + 87b^4
          - 48b^3 + 22b^2 - 7b + 1               (irreducible)
```

Palindromic `=>` roots pair `a <-> 1/a`; real coefficients `=>` pair
`a <-> conj a`. The trace polynomials `g(w)` of `w = a + 1/a` are TOTALLY REAL
(degrees 6 and 8), with 5 resp. 6 roots in `[-2,2]`. Hence **every non-real
embedding has `|a| = 1`** (verified: 10 of 12, 12 of 16).

On `|a| = 1`, so `conj a = 1/a`:

| fact | reason |
|---|---|
| `u` real | symbolic identity `u(1/a) = u(a)` (verified in Sage, both charts) |
| `T`, `U` real, so `D(T) = D(U) = 0` | `u` real `=>` `1+4u^2>0` `=>` `r` real `=>` `T = 1-r^2` real |
| `D(V) + D(W) = 0` | `conj X = 1/X` `=>` `conj V = uX` `=>` `W = 1/(1 - conj V)`; and `D(1/(1-z)) = D(z)`, `D(conj z) = -D(z)` |
| filling term contributes no volume | `log a` purely imaginary `=>` `lambda_core` purely imaginary `=>` `lambda_core^2` REAL |
| signs `eps_j = +1` | triangulation combinatorics, embedding-independent; all four positively oriented at the complete structure |

`=>` Borel regulator vanishes at every embedding `=>` **both classes are
TORSION** `=>` `Re[Delta R] in pi^2 Q`.

## Avenues

**(a) Explicit denominator bound, then pin.**
Compute `w_2(F)` for `F = Q(a)` and `Q(b)` (and the compositum if needed) via
the Merkurjev-Suslin formula for `|K_3^ind(F)_tors|`, per Zickert Thm 1.1. Get
an explicit `Q`. Since 301 digits allows `Q < 1.75e150`, any realistic bound
closes it. **Terminal condition:** an explicit integer `Q` with a cited theorem,
plus the continued-fraction reconstruction. THIS IS THE MAIN AVENUE.

**(b) Sharpen the torsion argument to a direct order bound.**
If the torsion order `m` can be exhibited directly (`m * xi = 0` for explicit
`m`), the denominator divides `m` and no `w_2` computation is needed.

**(c) Bypass rationality: exhibit the exact value.**
Find an explicit five-term / Abel relation chain expressing `Delta R` as a
rational multiple of `pi^2` symbolically. Hardest, but gives an unconditional
identity with no numerics at all.

**(d) Fallback.** If the denominator bound resists, submit the reduction plus
the torsion theorem plus the certified 301-digit reconstruction, stating the one
cited input.

## Terminal conditions

- (a) success: explicit `Q`, `Q^2 < 1/(2*err)`, conclude `-4/85`.
- (a) failure: no citable bound exists for these fields -> go to (b).
- (b) failure: no explicit annihilator -> go to (c).
- (c) failure: -> (d).

## Files

- `proof.tex` — to be rewritten with the torsion theorem.
- scripts this session: `/tmp/{falpha2,shapes,shapes300,exact_proof,confirm}.sage`,
  `/tmp/{alpha_component,shape_pairing,mechanism,find_relation,orient,filling_term}.py`
  — to be moved into `problems/3.1/scripts/`.
