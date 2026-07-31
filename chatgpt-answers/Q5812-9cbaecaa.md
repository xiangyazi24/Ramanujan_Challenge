ANSWER Q5812 9cbaecaa

# Status of the one-step Desnanot–Jacobi ansatz

## 1. Exact Pascal identity

For

[

G_{d,L}(F)=sum_{i=0}^{L}(-1)^ibinom{d-1+i}{i}binom{d+L}{L-i}F_{d-1+i},

]

the adjacent-parameter Pascal relations give a two-term boundary relation. Writing (E=d-1+i), one obtains

[

G_{d,L}(F)-G_{d+1,L-1}(F)=alpha_{d,L}F_{d-1}-beta_{d,L}F_{d+L},

]

with explicit rational binomial factors (alpha_{d,L},beta_{d,L}). The exact factors depend on the normalization convention for the Newton projector; after clearing denominators they are integers.

Therefore the determinant of two rows of such projectors satisfies the Desnanot–Jacobi reduction

[

det(U_m,V_m)-det(U_{m-1},V_{m-1})

]

The right hand side contains only endpoint shell entries.

## 2. What follows rigorously

For the actual vector shell (F=(Y,W)), the bilinearity of the determinant gives

[

H_m=frac{det(G_{d_m,L_m}(F),G_{D,L_m}(F))}{g_m}.

]

The Pascal relation implies that a one-step relation can only have the form

[

A_m H_m+B_m H_{m-1}=Q_m,

]

where (A_m,B_m) are the cleared Pascal multipliers and (Q_m) is a determinant involving the newly exposed boundary packets.

The maximality condition on (M) gives two possible boundary collapses:

- lower boundary: (d_M-1=lfloor(n-1)/2rfloor), so Cartier reduction applies to the first cell;

- upper boundary: (D+L_M=n-1), so the last shell index is the Apéry terminal index.

Thus the remaining packets are not arbitrary: they reduce to finitely many shell coefficients or, in the doubled-period gauge, to Apéry Casoratian terms.

## 3. Missing step

The requested conclusion would require proving that the entire boundary determinant equals a fixed scalar of logarithmic height O(log n):

[

Q_m=Q(n),qquad log |Q(n)|=O(log n).

]

The current identities only show that (Q_m) is supported on boundary packets. They do not show cancellation of all interior shell dependence.

In particular, the Apéry Casoratian controls

[

a_{r+1}b_r-a_rb_{r+1}=6/(r+1)^3,

]

but it controls only the two-dimensional Apéry module coordinates. The boundary packet determinant has extra shell coordinates unless an additional gauge identity identifies them with the Apéry basis.

## 4. Finite obstruction to the exact ansatz

A countermodel to the pure one-step claim is obtained by replacing the actual shell vector by

[

F=(F_1,F_2)

]

where both components satisfy the same certified shell recurrences but the endpoint jets are independent. The Pascal determinant identity still holds, but the boundary term contains an independent endpoint Wronskian. Hence no universal scalar (Q(n)) follows from Pascal algebra alone.

For the actual Apéry shell, this obstruction is not yet a disproof: it reduces the problem to a finite symbolic test:

1. Compute the endpoint packet basis produced by the two Pascal boundary operators.

1. Express every packet through the doubled-period gauge basis ((b_r,b_{r+1})) plus residual shell directions.

1. Check whether the residual coefficient matrix vanishes identically.

The decisive CAS target is therefore not a recurrence for (H_m), but the exact residual matrix

[

R_m=operatorname{Proj}_{text{shell}/text{Apéry module}}(Q_m).

]

If (R_m=0), then the desired identity follows and (Q) is a rational multiple of the Apéry Casoratian. If (R_mneq0), the one-step Desnanot–Jacobi route cannot yield an O(log n) scalar without adding additional carriers.

## Conclusion

The Pascal/Christoffel–Darboux mechanism genuinely reduces (H_M,H_{M-1}) to boundary data. The remaining theorem needed is a finite-dimensional boundary-gauge collapse, not a new gcd argument. The next exact symbolic computation should be the residual endpoint matrix above; it decides the ansatz completely.