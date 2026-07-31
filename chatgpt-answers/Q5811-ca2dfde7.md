ANSWER Q5811 ca2dfde7

# Q5811 audit response: correction of Q5804 and current status

## 1. Withdrawal of the unsupported fixed-R obstruction claim

The Q5804 statement that a fixed-R exterior-lattice ansatz already gave a rigorous \Omega(n) primitive-height obstruction was not justified. In particular:

- no explicit coefficient matrix was provided;

- no Smith normal form or invariant factors were computed;

- no specialization argument converted a rank statement into a lower bound on the distinguished Apéry integer minors;

- the Vandermonde/jet argument only controls generic lattice rank and does not exclude exceptional low-height integer specializations.

Therefore the claimed theorem is withdrawn. The correct status is: no fixed-R impossibility theorem has been established.

The new exact gcd data you supplied are precisely the kind of exceptional specialization phenomenon that invalidates that argument. Large initial heights can collapse after the actual shell normalization.

## 2. What can be concluded from the hostile data

The observations

- n=321: running gcd bit lengths

- n=200: hostile core

show that the endpoint-state cancellation is much stronger than a naive dimension count predicts.

They do not prove a polynomial carrier exists, but they rule out the previous argument.

## 3. Correct reduced problem

For lengths

[

L,L+1,ldots,L+R,

]

let

[

D_L=det(E_{d,L},E_{d+1,L})

]

(or the corresponding primitive exterior minors after the exact universal node/Smith content division).

The only admissible question is now whether there is a fixed R combination

[

C=sum_{i=0}^{R} c_i(n)D_{L+i}

]

with:

1. all common target primes retained with multiplicity at least one;

1. all universal content removed exactly;

1. primitive height o(n).

The previous rank argument does not decide this.

## 4. Required computation for R=2,3,4,...

The decisive audit should compute, for each fixed R:

[

G_R(n)=gcd(D_L,D_{L+1},ldots,D_{L+R})

]

after dividing only by certified universal factors.

Then factor

[

G_R(n)=U_R(n)V_R(n),

]

where U_R is the predicted universal polynomial/content contribution and V_R is the residual factor.

A positive result would require proving that V_R(n) divides a fixed polynomial in the shell parameters (n,D,N), or another explicitly bounded-height family.

A negative result would require an actual Smith obstruction:

- explicit matrix M_R(n);

- invariant factors;

- specialization proof;

- proof that every admissible coefficient vector leaves a surviving factor of linear height.

The previous answer did not provide these items and therefore cannot stand.

## 5. Current smallest open ansatz

After removing the invalid fixed-R obstruction, the smallest genuinely open ansatz is:

[

R=2 quad (three adjacent lengths).

]

If R=2 fails, test R=3,4,\ldots until either:

- an explicit bounded-degree residual carrier appears, or

- a real finite-rank obstruction with matrices and Smith data is obtained.

The next deliverable should therefore be a computational certificate, not a rank heuristic.

## 6. Checker skeleton

```python
for n in range(1,1001):
    shell = apery_first_cell_shell(n)
    minors = []
    for L in range(L0, L0+R+1):
        E = primitive_cross_minor(shell, L)
        E //= certified_universal_content(E, n, L)
        minors.append(E)
    g = gcd_list(minors)
    record(n, R, factor(g), bit_length(g))
```

The checker must additionally record:

- exact removed content;

- target-prime support;

- residual quotient after dividing candidate polynomial factors;

- Smith invariants of the coefficient lattice if a Bezout relation is attempted.

## Final verdict

The Q5804 fixed-R \Omega(n) obstruction is withdrawn. The actual Apéry shell data demonstrate that short translated length windows can have massive exceptional cancellation, so a proof must work with the normalized shell arithmetic itself. The correct next step is the exact R=2,3,4 gcd/residual-factor audit through n\ge1000 and only then a search for either a polynomial carrier or a genuine Smith obstruction.