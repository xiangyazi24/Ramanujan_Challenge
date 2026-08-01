# Quarter-point value search for the Apéry square-root branch

## Verdict

No closed form was found in the families specified in
`CODEX_SPEC_CRON_quarterpoint_value.md`.

The computation covers every prime

\[
\mathcal P=\{p<20000:p\text{ prime},\ p\equiv1\pmod {24}\},
\qquad |\mathcal P|=267.
\]

The first and last primes are 73 and 19993.  Every claimed fit below was
checked on all 267 primes.  The only candidate reaching ten hits is a dense
LLL interpolation that fits the ten probe primes and fails on all other 257
primes.  None of the a-priori monomial or quartic-twist candidates reaches ten
hits.

The computation is reproducible with

```bash
python3 problems/3.2/CRON_quarterpoint_value.py
python3 problems/3.2/CRON_quarterpoint_value.py --dump-csv
```

The first command runs every fit and prints expanded failure-prime lists.  The
second emits all 267 data rows.  The canonical data digest printed by the
script is

```text
SHA-256 311912d0c0eb6058ad3e8854b02339e029b8a0a558f671c55c551ece5fb4d729
```

## Data and normalization

For every \(p\in\mathcal P\), the script computes the recurrence through
\((p-1)/2\) and records

\[
v_p=\tau_{(p-1)/4}\pmod p.
\]

It also records all requested auxiliary data:

- \(p=x^2+6y^2\), with the unique \(x>0,y>0\);
- \(p=a^2+b^2\), with \(b>0\) and the signed odd coordinate
  \(a\equiv1\pmod4\);
- \(p=c^2+2d^2\), with \(c,d>0\) (this exists for every prime in the scan);
- \(G_1=\binom{(p-1)/2}{(p-1)/4}\) and
  \(G_2=\binom{(p-1)/4}{(p-1)/8}\) modulo \(p\);
- the eight integer Jacobsthal sums
  \[
  \phi(k)=\sum_{t\bmod p}\left(\frac{t(t^2+k)}p\right),
  \quad k\in\{1,2,3,6,-1,-2,-3,-6\};
  \]
- \(q_r=r^{(p-1)/4}\pmod p\) for \(r=2,3,6\), represented by
  \(+1\) or \(-1\);
- the midpoint \(\tau_{(p-1)/2}\).

Sanity checks are exact assertions in the program:

| Check | Result |
|---|---:|
| recurrence against direct \(\sqrt{\sum b_nz^n}\), all class-1 primes below 250 | 4/4 |
| \(v_p\ne0\) | 267/267 |
| \(\tau_{(p-1)/2}=(-2\mid p)=1\) | 267/267 |
| Gauss normalization \(G_1=2a\pmod p\) | 267/267 |
| \(q_6=q_2q_3\) | 267/267 |

The quartic-symbol distributions are

| symbol | \(-1\) | \(+1\) |
|---|---:|---:|
| \(q_2\) | 139 | 128 |
| \(q_3\) | 135 | 132 |
| \(q_6\) | 138 | 129 |

## 1. Bounded-rational monomial search

The twelve basic variables are

\[
x,y,a,b,\phi(1),\phi(2),\phi(3),\phi(6),
\phi(-1),\phi(-2),\phi(-3),\phi(-6).
\]

The search uses all 91 monomials of total degree at most two: the constant,
12 degree-one monomials, and 78 quadratic monomials.  For each monomial \(M\),
it tests

\[
v_p=uG_1^eM\pmod p,
\qquad e\in\{-1,0,1\},
\]

for every reduced rational \(u=m/n\) with

\[
0<|m|\le64,\qquad 1\le n\le64.
\]

Thus 273 scaled-monomial families are checked.  There are **zero** candidates
with at least ten hits.  The best support is only \(8/267\), attained by

\[
v_p=6G_1a.
\]

## 2. Small-coefficient polynomial search for \(v^2\) and \(v^4\)

Terms containing \(p\) vanish modulo \(p\).  The representation identities

\[
x^2+6y^2=p,\qquad a^2+b^2=p
\]

also give \(x^2\equiv-6y^2\) and \(a^2\equiv-b^2\).  The program therefore
uses the following nonredundant quotient basis for degree-at-most-two
polynomials:

\[
1,x,y,a,b,xy,xa,xb,y^2,ya,yb,ab,b^2.
\]

It performs three searches for each of \(v^2\) and \(v^4\):

1. exhaustive one- and two-term polynomials with nonzero coefficients in
   \([-32,32]\);
2. exhaustive three- and four-term polynomials with nonzero coefficients in
   \([-8,8]\), using meet-in-the-middle signatures on ten probes;
3. dense CRT/LLL probes on the first 10, 12, and 16 primes, retaining only
   coefficient vectors bounded by 32.

For each power, 340288 sparse coefficient/pair signatures are screened.

- For \(v^2\), no candidate matches the ten probe primes.
- For \(v^4\), LLL produces exactly one bounded candidate:

\[
\begin{aligned}
v_p^4={}&19+7x+15y+13a-14b+9xy+2xa+15xb+10y^2\\
        &{}-9ya+13yb+24ab+24b^2\pmod p.
\end{aligned}
\]

It has hit/miss count **10/257**.  Its exact hit set is

\[
H_{\rm LLL}=\{73,97,193,241,313,337,409,433,457,577\}.
\]

Its exact failure set is therefore

\[
\mathcal P\setminus H_{\rm LLL}.
\]

This is precisely the first ten-prime training set, so it is an interpolation
artifact rather than a near-identity.  The program prints the 257 failure
primes explicitly.

## 3. Quartic-twisted linear search

The targets are both \(v/G_1\) and \(vG_1\).  The algebraic factors are all
624 nonzero forms

\[
\alpha x+\beta y+\gamma a+\delta b,
\qquad \alpha,\beta,\gamma,\delta\in\{-2,-1,0,1,2\}.
\]

The character factors include:

- all eight Walsh products of
  \[
  \epsilon_{16}(p),q_2(p),q_3(p),
  \quad
  \epsilon_{16}(p)=
  \begin{cases}1,&p\equiv1\pmod {16},\\-1,&p\equiv9\pmod {16};\end{cases}
  \]
- \(\chi_4(r)=r^{(p-1)/4}\pmod p\) for every
  \(-12\le r\le12\), \(r\notin\{-1,0,1\}\);
- products of each such \(\chi_4(r)\) with every Walsh sign factor.

After identical factor vectors are merged, 32 distinct quartic factors remain.
There are **zero** candidates with at least ten hits.  The best support is
only \(4/267\):

\[
v/G_1=a,
\]

which holds only for

\[
\{6481,9601,9697,19441\}.
\]

## Three best a-priori near-misses

The LLL interpolation above is excluded from this ranking because it was
constructed from its ten hits.  Among the fixed monomial and quartic-twist
families, and requiring distinct support sets, the three best are:

1. \(v=6G_1a\), hit/miss count **8/259**, with exact hit set
   \[
   H_1=\{937,2833,4801,6073,6793,7321,11329,14281\};
   \]
2. \(v=-4G_1b\), hit/miss count **7/260**, with exact hit set
   \[
   H_2=\{409,2377,3217,7177,10369,14449,19081\};
   \]
3. \(v=-2G_1\phi(2)\), hit/miss count **6/261**, with exact hit set
   \[
   H_3=\{193,241,1153,10753,13633,17569\}.
   \]

Their exact failure sets are, respectively,

\[
\mathcal P\setminus H_1,\qquad
\mathcal P\setminus H_2,\qquad
\mathcal P\setminus H_3.
\]

This complement notation is exact because \(\mathcal P\) is the fully specified
267-prime set above; the script also emits each expanded list.

## Distribution after a quartic sign guess

To make the requested empirical distribution unambiguous, the script chooses
among the eight Walsh sign factors by maximal absolute empirical mean of the
centered normalized residues.  It then orients the overall sign to make that
mean nonnegative.  This descriptive rule selects

\[
s_p=-\epsilon_{16}(p)q_2(p).
\]

For

\[
z_p=\frac{1}{p}\operatorname{center}_p
\left(\frac{v_p}{2a s_p}\right)\in(-1/2,1/2],
\]

the empirical summary is

| statistic | value |
|---|---:|
| minimum | -0.499721 |
| first quartile | -0.150912 |
| median | 0.038315 |
| third quartile | 0.240010 |
| maximum | 0.499246 |
| mean | 0.033003 |
| population standard deviation | 0.261132 |

The counts in ten equal bins on \([-1/2,1/2)\), from left to right, are

```text
[17, 18, 19, 34, 28, 40, 32, 27, 28, 24]
```

There is no visible concentration at a fixed residue after this sign
normalization.  Together with the very low fit counts, the data argue against
a degree-two expression in the elementary quadratic-form/Jacobsthal data.
If a closed form exists, it likely needs a genuinely non-monomial finite-field
hypergeometric or Jacobi-sum quantity rather than only its low-degree CM trace
parameters.
