ANSWER Q7259 3ce2efee

# Sym² square-root audit of the Apéry Hasse polynomial

## Verdict

The strongest proved statement is a square-root factorization in the Sym²/elliptic description, but the actual square root does not give a new zero restriction beyond the existing Apéry recurrence data.

The factorization

\[
H_p(t)=\sum_{j=0}^{p-1}b_jt^j=\Delta(t)^\epsilon B_p(t)^2
\]

is structural: the Hasse polynomial is a Frobenius trace object coming from a symmetric square, so the odd discriminant factor accounts for the determinant character and the remaining factor is a square. However, the square root polynomial \(B_p\) is not an arbitrary square root; it is the Hasse polynomial of the underlying rank-two Franel/elliptic family after the required pullback and normalization.

The available theorem gives a geometric square root, not an independent recurrence constraint on the coefficients \(b_j\).

# 1. Exact square-root normalization

Let the elliptic/Franel family be written in its rank-two form with period vector

\[
\Pi(t)=\binom{\omega_1(t)}{\omega_2(t)}.
\]

The Apéry Picard--Fuchs object is the symmetric square:

\[
\operatorname{Sym}^2(\Pi).
\]

After the pullback used in `toric_mellin_square.tex`, the Hasse polynomial has the form

\[
H_p(t)=\chi_p(t)\,B_p(t)^2,
\]

where

\[
\chi_p(t)=\Delta(t)^{\epsilon_p}
\]

is the determinant/discriminant character contribution.

Thus the actual square root is the rank-two Frobenius polynomial

\[
\boxed{B_p(t)=\operatorname{Hasse}_{E_p}(\phi(t))\cdot u_p(t)}
\]

with:

- \(\phi\): the explicit pullback from the elliptic parameter to the Apéry parameter;
- \(u_p(t)\): the normalization factor removing the determinant character;
- the equality understood in \(\mathbf F_p[t]\).

The important point is that \(B_p\) is determined only up to multiplication by a scalar square unit unless the normalization is fixed.

# 2. Mellin coefficients

Write

\[
B_p(t)=\sum_{r=0}^{d_p}\beta_r t^r.
\]

Then

\[
H_p(t)=\Delta(t)^\epsilon\sum_{r,s}\beta_r\beta_s t^{r+s}.
\]

Therefore the Apéry coefficients satisfy

\[
\boxed{
 b_j=\sum_{r+s+e=j}d_e\beta_r\beta_s
}
\]

where

\[
\Delta(t)^\epsilon=\sum_e d_et^e.
\]

If the discriminant factor is removed first, the reduced coefficients are simply the quadratic convolution

\[
\boxed{
 h_j=\sum_{r+s=j}\beta_r\beta_s.
}
\]

In Mellin language, if

\[
\widehat B(\chi)=\sum_t B_p(t)\chi(t),
\]

then multiplication of coefficients corresponds to additive convolution of Mellin modes:

\[
\widehat{B^2}(\chi)=\sum_{\chi_1\chi_2=\chi}
\widehat B(\chi_1)\widehat B(\chi_2).
\]

This is exact, but it is only a change of coordinates.

# 3. Does the square root force extra vanishing restrictions?

No new restriction is proved.

A reciprocal square polynomial can have arbitrary coefficient cancellations. The square-root structure only says that the coefficient vector lies in the image of the quadratic Veronese map

\[
(\beta_r)\mapsto(\sum_{r+s=j}\beta_r\beta_s).
\]

It does not imply:

- fixed support restrictions;
- rank deficiency of the Hankel matrix;
- Jacobi-factor restrictions;
- p-adic valuation constraints on individual coefficients.

The reason is that quadratic convolution is not injective in a way compatible with coefficient vanishing. The equation

\[
\sum_{r+s=j}\beta_r\beta_s=0
\]

is a single quadratic relation among many Mellin coefficients.

# 4. Small-prime audit

Direct computation for small primes is consistent with the square factorization:

- the residual factor after removing \(\Delta^\epsilon\) is a square;
- the square root degree matches the rank-two family degree;
- zeros of \(H_p\) occur with the expected parity except at discriminant zeros.

However the observed zero patterns do not imply a theorem: they are compatible with generic quadratic convolution cancellation.

# 5. Hankel/Gram determinant attempt

Let

\[
M_j=\sum_r\beta_r\beta_{r+j}.
\]

Then Hankel matrices built from \(M_j\) have rank constraints depending on the finite degree of \(B_p\).

But these are constraints on the square-root coefficients \(\beta_r\), not directly on the Apéry coefficient zeros.

The map

\[
\beta\mapsto\beta*\beta
\]

forgets the linear structure. A zero determinant of a Hankel matrix does not imply a zero/nonzero pattern for the sequence \(b_j\).

# 6. Scoped no-go for the Mellin route

The Mellin transform is an invertible Fourier transform on \(\mathbf F_p^\times\). Therefore replacing \(B_p\) by its Mellin coefficients does not reduce the arithmetic complexity.

The only new statement is the exact convolution identity

\[
\widehat H(\chi)=\widehat\Delta^\epsilon*\widehat B*\widehat B.
\]

To obtain a stronger zero theorem one would need an additional input such as:

1. a support theorem for \(\widehat B\);
2. a rank drop in the Mellin convolution map;
3. a p-adic valuation theorem for the convolution terms;
4. an independent nonconcentration theorem for the Frobenius parameters.

The square-root factorization alone provides none of these.

Therefore the exact conclusion is:

\[
\boxed{
\text{Sym}^2\text{ squareness explains the factorization, but the actual }B_p
\text{ gives no presently proved improvement on Apéry zero-gap bounds.}
}
\]

D-finiteness, symmetric-square closure, or Mellin invertibility alone do not imply bounded conductor or finite-field zero nonconcentration.