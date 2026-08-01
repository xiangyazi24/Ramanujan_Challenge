ANSWER Q6660 f9aa6bef

# Hostile audit verdict

Both exponents survive. Neither theorem is broken. One statement-level repair is needed for the small-window range, and the absolute-index shift must be made explicit in the proof. The corrected leading constant is

```text
3^(4/3)/2,
```

not `3/2`.

The strongest bankable forms are the exact min-bounds

```text
d_D(r)
 <= min_{1<=Y<=D}
      { (3/2)Y(Y-1) + floor(D/(Y+1)) },
```

and, for an interval J of H consecutive regular indices,

```text
#{n in J: pi(n)=v}
 <= 1 + min_{1<=Y<=H-1}
          { (3/2)Y(Y-1) + floor((H-1)/(Y+1)) }
```

when `H>=2`, with the trivial bound `<=1` when `H=1`.

# Attack 1: the shifted restart chain

Verdict: PROVED.

The correct object is the absolute-index homogeneous solution

```text
y_n^(r) = b_r c_n - c_r b_n.
```

Its zeros in the interval `[r,r+D]` are exactly the returns to `pi(r)`. The recurrence state is always normalized with the absolute index:

```text
X_n(y) = [n^3 y_{n-1}:y_n].
```

It satisfies

```text
X_{n+1}(y) = M_n X_n(y)
```

with the global matrix `M_n`, not a reindexed matrix. If two consecutive zeros of `y^(r)` in the return row occur at

```text
n_i = r+e_i,
n_{i+1} = r+e_{i+1},
h = e_{i+1}-e_i,
```

then

```text
X_{n_i+1}(y^(r)) = [0:1]
```

and

```text
X_{n_{i+1}}(y^(r))
 = G_{h-1}(n_i+1)[0:1].
```

The bottom coordinate is therefore

```text
U_{h-1}(n_i+1)
 = U_{h-1}(r+e_i+1)
 = N_h(r+e_i),
```

so

```text
N_h(r+e_i)=0.
```

Thus the relative start `e_i` is a root of the translated polynomial

```text
E -> N_h(r+E),
```

not literally of `N_h(E)`. Translation is an automorphism of `F_p[E]`, so this polynomial is nonzero whenever `N_h` is nonzero, and it has the same degree `3(h-1)`. This is exactly what the degree count needs. The global indexing causes no loss.

The companion-at-r theorem does not require the rational identity (D1). It follows directly from the homogeneous solution, the Casoratian initialization, and the transfer formula.

# Attack 2: boundaries, first zero, and small D

Verdict: FIXABLE.

There is no substantive boundary failure, but the statement should be repaired as follows.

Use only intervals contained in

```text
I_p = {1,...,p-2}.
```

For the base theorem require

```text
1 <= r,
r+D <= p-2.
```

A terminal return at `p-2` is allowed. To detect it, the last transfer matrix used is `M_{p-3}`. No value at `p-1`, no matrix `M_{p-2}`, and no quantity at `z+h+1` is needed. The reset after the terminal zero is also not needed.

For the base row, `d=0` is a known zero of `y^(r)` but is not counted in `d_D(r)`. If the zeros are

```text
0=e_0<e_1<...<e_q<=D,
```

then `q=d_D(r)` is exactly the number of inter-zero gaps. Hence there is no additive `+1` in the base-return theorem.

For a general projective fiber in an arbitrary interval J, the interval need not begin at a hit. If the fiber contains `k>=1` points, then `k` is one plus the number of consecutive gaps. This is the unique source of the additive `+1`.

The earlier restriction `2<=Y` unnecessarily excluded `D=1`, `H=1`, and `H=2`. Allowing

```text
Y=1
```

fixes all small cases. Since consecutive zeros are impossible, there are no short gaps when `Y=1`, and every gap has length at least `2=Y+1`.

# Attack 3: all projective fibers and coordinate order

Verdict: PROVED.

Write the projective point in the same coordinate order as the orbit:

```text
v=[alpha:beta],
pi(n)=[b_n:c_n].
```

Then define

```text
y_n(v)=alpha c_n-beta b_n.
```

This gives

```text
y_n(v)=0
 <=> alpha c_n-beta b_n=0
 <=> [b_n:c_n]=[alpha:beta].
```

There is no `[B:A]` reversal once the coordinate convention is stated this way.

The extreme fibers are included:

```text
v=[1:0] gives y=c,
v=[0:1] gives y=-b.
```

The solution is never identically zero. Indeed,

```text
y_0=-beta,
y_1=alpha-5beta,
```

so `y_0=y_1=0` would force `alpha=beta=0`, impossible for a projective point. Also `pi(n)` is defined throughout `I_p`, since the nonzero Casoratian forbids `b_n=c_n=0`.

The target is fixed during each fiber count. In the base-return theorem it may depend on the row through `v=pi(r)`, but the estimate is uniform in `v`, so this row dependence is harmless.

# Attack 4: the min-formula and the constant

Verdict: PROVED for `3^(4/3)/2`; the advertised coefficient `3/2` is BROKEN as an output of this optimization.

For the base theorem the exact bound is

```text
B(D)=min_{1<=Y<=D}
       { (3/2)Y(Y-1)+floor(D/(Y+1)) }.
```

Put

```text
x=(D/3)^(1/3),
Y=ceil(x).
```

Then `1<=Y<=D`, `Y(Y-1)<=x(x+1)`, and `Y+1>=x+1>x`. Hence

```text
B(D)
 <= (3/2)x(x+1)+D/x
 = (9/2)x^2+(3/2)x
 = (3^(4/3)/2)D^(2/3)
   +(3^(2/3)/2)D^(1/3).
```

This is valid for every integer `D>=1`. For a general interval replace `D` by `H-1` and add `1`.

The leading constant from the restart-degree balance is therefore

```text
3^(4/3)/2 = 2.163374... .
```

The smaller coefficient `3/2` is not produced by the displayed short-gap and long-gap terms. It could only come from additional arithmetic input beyond this restart proof.

# Attack 5: bad primes and leading-coefficient drops

Verdict: PROVED for every prime `p>=7`, with no extra exclusions.

Every gap used in these theorems has

```text
2<=h<=p-3,
```

so it lies inside the proved range in which `N_h` is not the zero polynomial modulo `p`. A leading coefficient may drop modulo `p`, but that only decreases the degree and strengthens the root bound. The proof needs neither the exact leading coefficient, irreducibility, nor squarefreeness.

All transfer determinants used are products of nonzero factors `(n+1)^6` with absolute indices inside the regular nonwrapping range. Thus there are no singular-transfer primes beyond the already imposed `p>=7` and interval restrictions.

# Merge-ready ASCII LaTeX

```latex
\begin{theorem}[Uniform return multiplicity from a fixed base]
\label{thm:br-base-return}
Let $p\geq 7$ be prime, let
\[
 I_p=\{1,\ldots,p-2\},
 \qquad
 \pi(n)=[b_n:c_n]\in\mathbf P^1(\mathbf F_p),
\]
and suppose that
\[
 1\leq r,
 \qquad
 1\leq D,
 \qquad
 r+D\leq p-2.
\]
Put
\[
 d_D(r)=\#\{1\leq d\leq D:\pi(r+d)=\pi(r)\}.
\]
Then
\begin{equation}
 d_D(r)
 \leq
 \min_{1\leq Y\leq D}
 \left\{
   \frac32Y(Y-1)
   +\left\lfloor\frac{D}{Y+1}\right\rfloor
 \right\}.
 \label{eq:br-base-min}
\end{equation}
In particular,
\begin{equation}
 d_D(r)
 \leq
 \frac{3^{4/3}}2D^{2/3}
 +\frac{3^{2/3}}2D^{1/3},
 \label{eq:br-base-explicit}
\end{equation}
and hence
\[
 \max_{1\leq r\leq p-2-D}d_D(r)\ll D^{2/3}
\]
with an absolute implied constant.
\end{theorem}

\begin{proof}
For the fixed base $r$, define an absolute-index solution
\[
 y^{(r)}_n=b_r c_n-c_r b_n.
\]
It satisfies the Apery recurrence because it is a fixed linear
combination of the two solutions $b$ and $c$.  Moreover,
\[
 y^{(r)}_r=0,
 \qquad
 y^{(r)}_{r+1}
 =b_r c_{r+1}-b_{r+1}c_r
 =\frac1{(r+1)^3}\neq0,
\]
where the last equality is the Casoratian identity.  Thus
$y^{(r)}$ is not the zero solution.  Since $(b_n,c_n)$ is nonzero on
$I_p$, one also has
\begin{equation}
 y^{(r)}_n=0
 \quad\Longleftrightarrow\quad
 \pi(n)=\pi(r).
 \label{eq:br-return-zero}
\end{equation}

For any nonzero homogeneous solution $y$, put
\[
 X_n(y)=[n^3y_{n-1}:y_n].
\]
A direct use of the recurrence gives
\begin{equation}
 X_{n+1}(y)=M_nX_n(y),
 \qquad
 M_n=
 \begin{pmatrix}
 0&(n+1)^6\\
 -1&P(n)
 \end{pmatrix}.
 \label{eq:br-state-transfer}
\end{equation}
Indeed,
\[
 M_n
 \binom{n^3y_{n-1}}{y_n}
 =
 \binom{(n+1)^6y_n}{-n^3y_{n-1}+P(n)y_n}
 =
 \binom{(n+1)^6y_n}{(n+1)^3y_{n+1}},
\]
which represents $[(n+1)^3y_n:y_{n+1}]$.
All matrices used below are invertible because their determinants are
$(n+1)^6\neq0$ in the regular nonwrapping range.  Consequently a
nonzero solution cannot have two consecutive zeros.  If $y_z=0$ and
there is a later zero in the regular range, then
\[
 X_z(y)=[1:0],
 \qquad
 X_{z+1}(y)=M_z[1:0]=[0:1].
\]
Thus every zero resets the state to the fixed point $[0:1]$ one step
later.

Let two consecutive zeros occur at $z$ and $z+h$.  Then $h\geq2$, and
using the proved transfer formula with the absolute start $z+1$ gives
\begin{align*}
 X_{z+h}(y)
 &=G_{h-1}(z+1)[0:1]\\
 &=\bigl[(z+h)^6U_{h-2}(z+1):U_{h-1}(z+1)\bigr].
\end{align*}
The left side is $[1:0]$.  Since the transfer product is invertible,
the displayed top coordinate is nonzero, and therefore
\begin{equation}
 U_{h-1}(z+1)=0.
 \label{eq:br-U-gap}
\end{equation}
By the identity $N_h(X)=U_{h-1}(X+1)$, this is
\begin{equation}
 N_h(z)=0.
 \label{eq:br-N-gap}
\end{equation}
Notice that all indices here are absolute indices.

List the zeros of $y^{(r)}$ in $[r,r+D]$ as
\[
 r=n_0<n_1<\cdots<n_q,
\]
so that $q=d_D(r)$ by \eqref{eq:br-return-zero}.  Write
\[
 e_i=n_i-r,
 \qquad
 h_i=n_{i+1}-n_i=e_{i+1}-e_i.
\]
For a fixed gap length $h$, equation \eqref{eq:br-N-gap} says that each
corresponding relative start $e_i$ is a root of
\[
 F_{r,h}(E)=N_h(r+E).
\]
Translation is an automorphism of $\mathbf F_p[E]$.  Hence the proved
nonvanishing of $N_h$ in the range $h<p$ implies that $F_{r,h}$ is not
the zero polynomial, and
\[
 \deg F_{r,h}\leq3(h-1).
\]
Thus, for each fixed $h$, at most $3(h-1)$ of the gaps can have length
$h$.

Fix $Y$ with $1\leq Y\leq D$.  The number of gaps with $h_i\leq Y$ is
at most
\[
 \sum_{h=2}^Y3(h-1)=\frac32Y(Y-1).
\]
Every remaining gap has length at least $Y+1$, while
\[
 \sum_{i=0}^{q-1}h_i=n_q-r\leq D.
\]
Hence the number of remaining gaps is at most
\[
 \left\lfloor\frac{D}{Y+1}\right\rfloor.
\]
Since $q$ is exactly the total number of gaps, this proves
\eqref{eq:br-base-min}.  This argument includes $D=1$: take $Y=1$,
which gives $d_1(r)=0$.

For the explicit estimate, put
\[
 x=(D/3)^{1/3},
 \qquad
 Y=\lceil x\rceil.
\]
Then $1\leq Y\leq D$, $Y(Y-1)\leq x(x+1)$, and $Y+1>x$.  Therefore
\begin{align*}
 d_D(r)
 &\leq\frac32x(x+1)+\frac{D}{x}\\
 &=\frac92x^2+\frac32x\\
 &=\frac{3^{4/3}}2D^{2/3}
   +\frac{3^{2/3}}2D^{1/3}.
\end{align*}
This proves \eqref{eq:br-base-explicit}.
\end{proof}

\begin{theorem}[Uniform projective-fiber bound in a window]
\label{thm:br-projective-fiber}
Let $p\geq7$ be prime, let $J$ be an interval of $H$ consecutive
integers contained in $I_p$, and let $v\in\mathbf P^1(\mathbf F_p)$.
If $H=1$, then
\[
 \#\{n\in J:\pi(n)=v\}\leq1.
\]
If $H\geq2$, then
\begin{equation}
 \#\{n\in J:\pi(n)=v\}
 \leq
 1+
 \min_{1\leq Y\leq H-1}
 \left\{
   \frac32Y(Y-1)
   +\left\lfloor\frac{H-1}{Y+1}\right\rfloor
 \right\}.
 \label{eq:br-fiber-min}
\end{equation}
Consequently,
\begin{equation}
 \#\{n\in J:\pi(n)=v\}
 \leq
 1+rac{3^{4/3}}2(H-1)^{2/3}
  +\frac{3^{2/3}}2(H-1)^{1/3},
 \label{eq:br-fiber-explicit}
\end{equation}
and in particular the fiber size is $O(H^{2/3})$, uniformly in
$p$, $J$, and $v$.
\end{theorem}

\begin{proof}
Write the projective point in the same coordinate order as the orbit:
\[
 v=[\alpha:\beta],
 \qquad
 \pi(n)=[b_n:c_n].
\]
Define
\[
 y_n=\alpha c_n-\beta b_n.
\]
This is a homogeneous solution of the Apery recurrence.  It is not the
zero solution, since
\[
 y_0=-\beta,
 \qquad
 y_1=\alpha-5\beta,
\]
and simultaneous vanishing would imply $\alpha=\beta=0$.  Also the
Casoratian identity implies $(b_n,c_n)\neq(0,0)$ on $I_p$.  Therefore
\begin{equation}
 y_n=0
 \quad\Longleftrightarrow\quad
 \alpha c_n-\beta b_n=0
 \quad\Longleftrightarrow\quad
 \pi(n)=v.
 \label{eq:br-fiber-zero}
\end{equation}
This includes the two extreme fibers: $v=[1:0]$ gives $y=c$, while
$v=[0:1]$ gives $y=-b$.

List the zeros of $y$ in $J$ as
\[
 z_1<z_2<\cdots<z_k.
\]
If $k\leq1$, the result is immediate.  For each $1\leq i<k$, put
\[
 h_i=z_{i+1}-z_i.
\]
The absolute-index reset computation in the proof of
Theorem~\ref{thm:br-base-return} applies to this arbitrary nonzero
homogeneous solution and gives
\[
 h_i\geq2,
 \qquad
 N_{h_i}(z_i)=0.
\]
For a fixed $h$, the starts $z_i$ are distinct roots of the nonzero
polynomial $N_h$, whose degree is at most $3(h-1)$.  Thus at most
$3(h-1)$ gaps have length $h$.

Fix $Y$ with $1\leq Y\leq H-1$.  The number of gaps of length at most
$Y$ is at most
\[
 \sum_{h=2}^Y3(h-1)=\frac32Y(Y-1).
\]
Every other gap has length at least $Y+1$, and
\[
 \sum_{i=1}^{k-1}h_i=z_k-z_1\leq H-1.
\]
Hence the number of long gaps is at most
\[
 \left\lfloor\frac{H-1}{Y+1}\right\rfloor.
\]
Since $k$ is one plus the number of gaps, this proves
\eqref{eq:br-fiber-min}.  The argument does not require the left endpoint
of $J$ to be a zero; the additive $1$ accounts for the first zero of the
fiber wherever it occurs in $J$.  A zero at the right endpoint of $J$
or at $p-2$ is harmless, since no reset after the terminal zero is used.

Finally apply the explicit estimate from
Theorem~\ref{thm:br-base-return} with $D=H-1$ and add $1$.  This gives
\eqref{eq:br-fiber-explicit}.
\end{proof}
```

# Final referee decision

[THM-BASE-RETURN-2/3]: ACCEPT AFTER MINOR REPAIR.

[THM-PROJECTIVE-FIBER-2/3]: ACCEPT AFTER MINOR REPAIR.

The repairs are statement hygiene, not a change of exponent or scope:

```text
- expose the absolute shift N_h(r+E);
- permit Y=1;
- distinguish the base row's no-+1 count from the general fiber's +1;
- use 3^(4/3)/2 as the leading constant.
```

No exceptional primes beyond `p>=7` are needed.