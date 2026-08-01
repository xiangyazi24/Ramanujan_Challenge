# Jacobsthal/apparition recon for the Apéry square-root branches

## Scope and conventions

This is an empirical scan of every prime \(5\le p<3000\) (428 primes).  I use
the natural finite branch and its full truncation window

\[
D_p=\frac{p-1}{2},\quad a_j=\tau_j\qquad
(p\bmod24\in\{1,5,7,11\}),
\]

and

\[
D_p=\frac{p-3}{2},\quad a_j=\sigma_j\qquad
(p\bmod24\in\{13,17,19,23\}).
\]

Thus \(Y_p=\{0\le j\le D_p:a_j=0\pmod p\}\).  The implementation is
`research/scripts/codex_high_jacobsthal_recon.py`.  Its recurrence output was
checked coefficientwise against the defining series for both branches for
every prime below 200.

The derived sigma recurrence is

\[
4(n+2)^2\sigma_{n+2}
=2(68n^2+238n+209)\sigma_{n+1}-(2n+3)^2\sigma_n,
\]

with \(\sigma_0=1,\ \sigma_1=39/2\).  This complements the given tau
recurrence.

## Main findings

### 1. The quarter law is a reflection-sign law

For all 428 primes and every \(0\le j\le D_p\), the complete coefficient
arrays obey

\[
\boxed{\ a_{D_p-j}=\left(\frac{-2}{p}\right)a_j\pmod p\ }.
\]

The mod-24 signs are

| \(p\bmod24\) | 1 | 5 | 7 | 11 | 13 | 17 | 19 | 23 |
|---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| \((-2\mid p)\) | + | - | - | + | - | + | + | - |

Consequently, a forced center zero occurs exactly when the sign is negative
and \(D_p\) is even.  This gives exactly the updated table:

* \(p\equiv5\pmod{24}\):
  \(\tau_{(p-1)/4}=0\), verified \(57/57\).
* \(p\equiv23\pmod{24}\):
  \(\sigma_{(p-3)/4}=0\), verified \(55/55\).
* The other six natural-branch center/floor-center coefficients are nonzero,
  verified \(316/316\).

This also resolves the requested \(11\) and \(17\) classes.  Their truncation
degree is odd and their reflection sign is \(+1\), so the structured object is
an equal, nonzero central pair rather than a zero:

\[
\begin{aligned}
p\equiv11 &: \quad
\tau_{(p-3)/4}=\tau_{(p+1)/4}\ne0 &&(54/54),\\
p\equiv17 &: \quad
\sigma_{(p-5)/4}=\sigma_{(p-1)/4}\ne0 &&(55/55).
\end{aligned}
\]

No other universal fixed rational point was found in either class.  Among
all positions obtained by floor, ceiling, or nearest-integer rounding of
\(kD_p/q\), \(2\le q\le24\), every candidate hit at most one prime: \(1/54\)
in class 11 and \(1/55\) in class 17.

### 2. A new sigma octant skeleton

There is a second universal zero law, one scale below the center:

\[
\begin{array}{ll}
p\equiv13\pmod{24}: &
\displaystyle
\sigma_{(p-5)/8}=\sigma_{(3p-7)/8}=0 \qquad(53/53),\\[2mm]
p\equiv23\pmod{24}: &
\displaystyle
\sigma_{(p-7)/8}=\sigma_{(p-3)/4}
=\sigma_{(3p-5)/8}=0 \qquad(55/55).
\end{array}
\]

These are respectively the reflected quarter-window pair and the
quarter-window/center/three-quarter-window triple.  No analogous universal
octant zero appeared in the tau classes.

The sharp finite-range description is therefore

\[
\begin{aligned}
p\equiv5:\quad
Y_p&=\{D_p/2\}\ \sqcup\ \text{zero or one sporadic reflected pair},\\
p\equiv23:\quad
Y_p&=\{(p-7)/8,D_p/2,(3p-5)/8\}\\
&\quad\sqcup\ \text{zero, one, or two sporadic reflected pairs}.
\end{aligned}
\]

### 3. Full-window zero statistics

| \(p\bmod24\) | branch | primes | quarter zeros | primes with any zero | histogram of \(|Y_p|\) |
|---:|:---:|---:|---:|---:|:---|
| 1 | tau | 46 | 0 | 16 | \(0:30,\ 2:12,\ 4:4\) |
| 5 | tau | 57 | 57 | 57 | \(1:47,\ 3:10\) |
| 7 | tau | 54 | 0 | 14 | \(0:40,\ 2:12,\ 4:2\) |
| 11 | tau | 54 | 0 | 14 | \(0:40,\ 2:12,\ 4:2\) |
| 13 | sigma | 53 | 0 | 53 | \(2:42,\ 4:9,\ 6:2\) |
| 17 | sigma | 55 | 0 | 13 | \(0:42,\ 2:8,\ 4:5\) |
| 19 | sigma | 54 | 0 | 11 | \(0:43,\ 2:10,\ 4:1\) |
| 23 | sigma | 55 | 55 | 55 | \(3:40,\ 5:12,\ 7:3\) |

Every observed zero set is reflection-stable.  The parity of its size is
completely explained by the reflection sign and whether the reflection has an
integer fixed point.

### 4. The quadratic-form premise needs correction

All 57 class-5 primes have a unique enumerated positive representation
\(p=2x^2+3y^2\).  But no class-23 prime can have this representation: for odd
\(p=2x^2+3y^2\), reduction modulo 8 gives \(p\equiv3\) or \(5\pmod8\), whereas
\(p\equiv23\pmod{24}\) gives \(p\equiv7\pmod8\).  The scan found \(0/55\), as
it must.

Likewise, \(x^2+6y^2\) represents prime classes 1 and 7 modulo 24, not class
19.  The scan found a positive representation for all \(46/46\) class-1
primes and none for the class-19 primes.

Therefore the statement that both quarter-zero classes are exactly the
non-principal discriminant-\(-24\) genus is false as written.  The form
\(2x^2+3y^2\) covers classes 5 and 11, while the quarter-zero classes are 5
and 23.  The observed quarter law instead matches the reflection character
\((-2\mid p)\) together with the parity of the truncation degree.

### 5. Jacobsthal position and value hunts were negative beyond the skeleton

For the ten class-5 primes with an extra reflected zero pair, write the lower
zero as \(J-\delta\), \(J=(p-1)/4\).  Exhaustive fitting of
\(\delta=|Ax+By|\) for \(|A|,|B|\le64\) hit at most \(1/10\).  The rounded
position candidates

\[
p\frac{x}{x+y},\quad p\frac{y}{x+y},\quad
p\frac{x^2}{x^2+y^2},\quad p\frac{y^2}{x^2+y^2},\quad
2x^2,\quad3y^2
\]

hit no extra zero in any of the ten cases.  In class 11, the larger search
\(J\pm(Ax+By)\), \(|A|,|B|\le64\), hit at most \(2/54\), and 40 of the 54
primes have no zero at all.

For the class-5 near-quarter values at offsets \(-2,-1,+1,+2\), I tested

\[
c x,\quad c y,\quad cxy,\quad c x/y,\quad c(4x^2),\quad c(3y^2)\pmod p
\]

for every fixed rational \(c=a/b\) with \(|a|,b\le48\).  Excluding \(p=5\),
where the \(-2\) index is outside the window, the best support of any candidate
at each offset was only \(3/56\).  There is no fitted constant law.  Reflection
does give the exact universal relation

\[
\tau_{J+d}=-\tau_{J-d}\qquad(p\equiv5\pmod{24}).
\]

For class 1, where \(p=x^2+6y^2\), the same candidate grid for
\(\tau_{(p-1)/4}\) had best support \(2/46\).  The proposed class-19
quadratic-form test is not arithmetically defined, and its natural sigma
center value is nonzero for all \(54/54\) primes.

## Sharpest conjecture supported by this scan

For every odd prime away from the singular small characteristics, the natural
tau/sigma truncation satisfies the exact reflection law

\[
a_{D_p-j}=(-2\mid p)a_j.
\]

The deterministic zero skeleton consists of:

1. the reflection-fixed center precisely in classes 5 (tau) and 23 (sigma);
2. a sigma octant pair in classes 13 and 23;
3. only reflection-paired sporadic zeros beyond that skeleton.

Up to \(p<3000\), the sporadic pairs show no Jacobsthal parametrization by the
positive representation \(p=2x^2+3y^2\), no small linear position law in
\((x,y)\), and no classical fixed-constant law for the adjacent values.  The
data support a reflection-character apparition law, not a discriminant-\(-24\)
class-group parametrization of the residual zero sets.

## Complete zero sets in the two quarter-zero classes

The following are all primes below 3000, not selected examples.

### \(p\equiv5\pmod{24}\): tau, with \(p=2x^2+3y^2\)

```text
p=5    (x,y)=(1,1)   Y=[1]
p=29   (x,y)=(1,3)   Y=[7]
p=53   (x,y)=(5,1)   Y=[13]
p=101  (x,y)=(7,1)   Y=[6,25,44]
p=149  (x,y)=(1,7)   Y=[37]
p=173  (x,y)=(7,5)   Y=[43]
p=197  (x,y)=(5,7)   Y=[49]
p=269  (x,y)=(11,3)  Y=[67]
p=293  (x,y)=(5,9)   Y=[73]
p=317  (x,y)=(11,5)  Y=[79]
p=389  (x,y)=(11,7)  Y=[97]
p=461  (x,y)=(7,11)  Y=[115]
p=509  (x,y)=(1,13)  Y=[127]
p=557  (x,y)=(5,13)  Y=[139]
p=653  (x,y)=(17,5)  Y=[163]
p=677  (x,y)=(1,15)  Y=[169]
p=701  (x,y)=(13,11) Y=[41,175,309]
p=773  (x,y)=(7,15)  Y=[193]
p=797  (x,y)=(19,5)  Y=[199]
p=821  (x,y)=(17,9)  Y=[205]
p=941  (x,y)=(17,11) Y=[235]
p=1013 (x,y)=(13,15) Y=[253]
p=1061 (x,y)=(23,1)  Y=[83,265,447]
p=1109 (x,y)=(11,17) Y=[277]
p=1181 (x,y)=(7,19)  Y=[161,295,429]
p=1229 (x,y)=(19,13) Y=[307]
p=1277 (x,y)=(25,3)  Y=[122,319,516]
p=1301 (x,y)=(23,9)  Y=[325]
p=1373 (x,y)=(5,21)  Y=[321,343,365]
p=1493 (x,y)=(25,9)  Y=[19,373,727]
p=1613 (x,y)=(25,11) Y=[403]
p=1637 (x,y)=(5,23)  Y=[409]
p=1709 (x,y)=(29,3)  Y=[427]
p=1733 (x,y)=(23,15) Y=[433]
p=1877 (x,y)=(1,25)  Y=[469]
p=1901 (x,y)=(17,21) Y=[475]
p=1949 (x,y)=(31,3)  Y=[21,487,953]
p=1973 (x,y)=(7,25)  Y=[493]
p=1997 (x,y)=(31,5)  Y=[499]
p=2069 (x,y)=(31,7)  Y=[517]
p=2141 (x,y)=(23,19) Y=[535]
p=2213 (x,y)=(13,25) Y=[553]
p=2237 (x,y)=(5,27)  Y=[559]
p=2309 (x,y)=(19,23) Y=[577]
p=2333 (x,y)=(25,19) Y=[583]
p=2357 (x,y)=(29,15) Y=[589]
p=2381 (x,y)=(23,21) Y=[595]
p=2477 (x,y)=(35,3)  Y=[619]
p=2549 (x,y)=(29,17) Y=[637]
p=2621 (x,y)=(7,29)  Y=[334,655,976]
p=2693 (x,y)=(35,9)  Y=[673]
p=2741 (x,y)=(37,1)  Y=[685]
p=2789 (x,y)=(31,17) Y=[575,697,819]
p=2837 (x,y)=(25,23) Y=[709]
p=2861 (x,y)=(13,29) Y=[715]
p=2909 (x,y)=(19,27) Y=[727]
p=2957 (x,y)=(35,13) Y=[739]
```

### \(p\equiv23\pmod{24}\): sigma

There is no positive \(2x^2+3y^2\) representation for any row.

```text
p=23   Y=[2,5,8]
p=47   Y=[5,11,17]
p=71   Y=[4,8,17,26,30]
p=167  Y=[20,41,62]
p=191  Y=[23,37,47,57,71]
p=239  Y=[29,59,89]
p=263  Y=[32,65,98]
p=311  Y=[38,77,116]
p=359  Y=[44,89,134]
p=383  Y=[47,82,95,108,143]
p=431  Y=[49,53,88,107,126,161,165]
p=479  Y=[59,119,179]
p=503  Y=[62,125,188]
p=599  Y=[74,149,224]
p=647  Y=[80,106,161,216,242]
p=719  Y=[59,89,119,179,239,269,299]
p=743  Y=[74,92,185,278,296]
p=839  Y=[104,209,314]
p=863  Y=[107,215,323]
p=887  Y=[110,221,332]
p=911  Y=[46,113,214,227,240,341,408]
p=983  Y=[122,245,368]
p=1031 Y=[118,128,257,386,396]
p=1103 Y=[137,275,413]
p=1151 Y=[143,287,431]
p=1223 Y=[152,305,458]
p=1319 Y=[164,329,494]
p=1367 Y=[170,341,512]
p=1439 Y=[179,359,539]
p=1487 Y=[185,371,557]
p=1511 Y=[188,377,566]
p=1559 Y=[194,364,389,414,584]
p=1583 Y=[197,395,593]
p=1607 Y=[200,248,401,554,602]
p=1823 Y=[227,455,683]
p=1847 Y=[230,310,461,612,692]
p=1871 Y=[233,438,467,496,701]
p=2039 Y=[254,509,764]
p=2063 Y=[257,515,773]
p=2087 Y=[61,260,521,782,981]
p=2111 Y=[263,527,791]
p=2207 Y=[275,551,827]
p=2351 Y=[293,587,881]
p=2399 Y=[299,428,599,770,899]
p=2423 Y=[302,605,908]
p=2447 Y=[305,611,917]
p=2543 Y=[317,635,953]
p=2591 Y=[323,647,971]
p=2663 Y=[332,665,998]
p=2687 Y=[335,671,1007]
p=2711 Y=[338,677,1016]
p=2879 Y=[359,719,1079]
p=2903 Y=[362,725,1088]
p=2927 Y=[365,731,1097]
p=2999 Y=[374,749,1124]
```
