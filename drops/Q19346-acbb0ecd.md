ANSWER Q19346 acbb0ecd

# R1 verdict

I do **not** see an existing theorem in the literature that gives the needed
worst-case-in-\(n\) upper bound for an Apéry/P-recursive term. The closest
results split into three classes: (i) constant-coefficient recurrence results,
which control existence/size of prime divisors or averages over indices;
(ii) factorial/binomial results, whose pointwise strength comes from
Kummer/Legendre carry structure unavailable for a sum with cancellation; and
(iii) diagonal/Lucas/Dwork/automatic results, which are local in one
characteristic and do not correlate the defining primes \(p\asymp n\).

For the deep attack, the route I would pursue is a **characteristic-zero
anti-diagonal carrier followed by a short-cluster, two-row,
three-characteristic resultant**. I can prove the carrier and the cluster
reduction exactly. The first tempting “extra equation” after reflection is
false in the finite data: after removing the structural central factor of the
reflected gap polynomial, the quotient does *not* vanish at actual Apéry
zeros. The smallest surviving premise is therefore genuinely a three-column
identity/resultant before scalarizing the rows.

This is option **(iii)**: an exact new reduction, the first failed implication,
and a narrowly scoped replacement lemma.

---

# Part A. Literature / technique audit

I searched specifically for a bound of the shape

\[
\#\{p\in(cn,n]:p\mid u_n\}
\quad\hbox{or}\quad
\sum_{\substack{p\in(cn,n]\\p\mid u_n}}\log p
\]

holding for **every** \(n\), with \(u_n\) holonomic/P-recursive,
binomial-sum, or Apéry-like. I did not find a theorem that transfers.

## A1. Erdős--Graham--Ruzsa--Straus and recent binomial work

P. Erdős, R. L. Graham, I. Z. Ruzsa, E. G. Straus,
*On the Prime Factors of \(\binom{2n}{n}\)*,
Math. Comp. **29** (1975), 83--92, obtains quantitative information on
small-prime factors of the central binomial coefficient, including
worst-case and average statements. Its pointwise strength is based on the
Kummer/Legendre carry criterion.

That criterion has no analogue for
\[
b_r=\sum_k\binom rk^2\binom{r+k}{k}^2,
\]
because \(p\mid b_r\) is cancellation among \(p\)-adic units. The
Gross--Koblitz audit already in this repository makes this exact: in the
top-half regime there is a whole minimum-slope packet, and its unit sum is
again the reflected Apéry coefficient. Thus the EGRS mechanism does not
transfer.

A recent related paper is Hung M. Bui--Kyle Pratt--Alexandru Zaharescu,
*Binomial coefficients with divisors avoiding an interval*,
arXiv:2605.21221 (20 May 2026). They prove strong “divisor near \(n\)”
results for sufficiently large \(k\), and under GRH construct small-\(k\)
counterexamples to the old Erdős--Graham expectation. Again the mechanism
is factorial/carry/sieve structure, not cancellation in a holonomic sum.

**Verdict:** relevant analogy and warning, no transferable theorem.

## A2. Luca on binary holonomic sequences

Florian Luca, *Prime divisors of binary holonomic sequences*,
Adv. Appl. Math. **40** (2008), 168--179,
DOI 10.1016/j.aam.2006.12.001.

For a rational sequence satisfying a second-order recurrence with polynomial
coefficients, subject to a non-eventual-binary-recurrence condition, Luca
proves that the product of numerators and denominators of the nonzero terms
through \(N\) has at least \(c\log N\) distinct prime factors.

This is the closest theorem I found by recurrence class, but its conclusion
has the wrong quantifiers and direction: it is a **lower bound** on distinct
prime factors of a product of many terms, not an upper bound for
proportional-window primes dividing one prescribed term.

**Verdict:** closest class match, wrong observable.

## A3. Schinzel--Stewart / primitive-divisor theory

Schinzel--Stewart and later Stewart results for nondegenerate binary linear
recurrences give primitive divisors and lower bounds for the greatest prime
factor of a term. These are constant-coefficient recurrences, whereas
Apéry's coefficients depend cubically on \(n\). They also prove existence or
largeness of prime factors, while we need an upper bound on how many primes
in a moving interval divide one term.

**Verdict:** wrong recurrence class and wrong direction.

## A4. Results on \(n\mid u_n\)

Juan José Alba González--Florian Luca--Carl Pomerance--Igor Shparlinski,
*On numbers \(n\) dividing the \(n\)-th term of a linear recurrence*,
Proc. Edinburgh Math. Soc. **55** (2012), 271--289,
DOI 10.1017/S0013091510001355, gives upper and lower bounds for
\[
\#\{m\le x:m\mid u_m\}
\]
for nondegenerate constant-coefficient linear recurrences with simple roots.
Related work of Sanna and others studies densities involving
\(\gcd(m,u_m)\) and ranks of apparition.

Our event
\[
p\mid b_{n-p},\qquad p\in(n/2,n]
\]
has both modulus and coefficient index moving together. There is no known
divisibility-sequence law or rank-of-apparition parametrization converting
these primes to divisors of a fixed recurrence index.

**Verdict:** no transfer.

## A5. Delaygue / Dwork congruences

Éric Delaygue, *Arithmetic properties of Apéry-like numbers*,
Compos. Math. **154** (2018), 249--274,
DOI 10.1112/S0010437X17007552, proves \(p\)-adic valuation results for
Apéry-like factorial-ratio multisums and an effective criterion for the
\(p\)-Lucas property for almost all primes. Delaygue--Rivoal--Roques'
Dwork formal congruence machinery gives stronger local congruence and
integrality packages for broad hypergeometric families.

This is directly relevant to why
\[
b_{qp+r}\equiv b_qb_r\pmod p
\]
exists. It does **not** correlate the zero events for different
characteristics. Our problem begins after Lucas has done its job.

**Verdict:** supplies the local reduction, not the missing horizontal theorem.

## A6. Automaticity / Christol / rational diagonals

Eric Rowland--Reem Yassawi,
*Automatic congruences for diagonals of rational functions*,
J. Théorie des Nombres de Bordeaux **27** (2015), 245--288,
DOI 10.5802/jtnb.901, constructs automata modulo \(p^\alpha\) for diagonals
of rational functions and treats Apéry-type congruences.

For each **fixed** \(p\), this can describe
\(\{r:p\mid b_r\}\) automatically. But the automaton depends on \(p\), and
for \(0\le r<p\) the target is exactly the seed row \(Z_p\). Automaticity
therefore does not couple \(Z_p,Z_q,Z_\ell\) for three characteristics.
Derksen-type characteristic-\(p\) zero-set structure has the same limitation.

**Verdict:** exact finite-\(p\) structure, no worst-case cross-\(p\) control.

### Part A conclusion

As of 22 September 2026, I found no black-box theorem proving
\[
N_q(n)=o(n/\log n)
\]
for every \(n\), even for \(q=1\). The literature audit supports the
repository diagnosis: the missing input must be horizontal/cross-
characteristic, not another per-prime marginal estimate.

---

# Part B. One coherent attack

## B1. A single characteristic-zero array carries every top-half event

Use the proved diagonal representation
\[
b_r\equiv[x^ry^rz^rw^r]F^{p-1}\pmod p,
\qquad
F=(1-x-y)(1-z-w)-xyzw.
\]

Set
\[
G=xyzw\,F
\]
and define the integral two-parameter array
\[
C_{s,m}:=[x^sy^sz^sw^s]G^m\in\mathbb Z.
\tag{B.1}
\]

Since \(G^m=(xyzw)^mF^m\),
\[
C_{m+r,m}=[x^ry^rz^rw^r]F^m.
\tag{B.2}
\]

Now fix **\(n\)**, not \(p\), and form
\[
H_n(T)
=\sum_{m\ge0}C_{n-1,m}T^m
=[x^{n-1}y^{n-1}z^{n-1}w^{n-1}]
  \frac1{1-TG}.
\tag{B.3}
\]

Every monomial of \(G\) has each coordinate exponent equal to \(1\) or \(2\).
Hence
\[
C_{s,m}=0\quad\text{unless}\quad m\le s\le2m,
\tag{B.4}
\]
so
\[
\left\lceil\frac{n-1}{2}\right\rceil\le m\le n-1.
\tag{B.5}
\]
After \(p=m+1\), this is exactly the q=1 top window.

### Carrier lemma

For prime \(p\in(n/2,n]\), put
\[
m=p-1,\qquad r=n-p=n-1-m.
\]
Then
\[
[T^{p-1}]H_n
=C_{n-1,p-1}
=[x^ry^rz^rw^r]F^{p-1},
\]
and therefore
\[
\boxed{
p\mid b_{n-p}
\iff
p\mid[T^{p-1}]H_n.
}
\tag{B.6}
\]

Thus
\[
\boxed{
N_1(n)=
\#\{p\in(n/2,n]\text{ prime}:p\mid[T^{p-1}]H_n\}.
}
\tag{B.7}
\]

This is a genuine quantifier change: for fixed \(n\), all the changing
characteristics now test coefficients of **one integral polynomial
\(H_n(T)\)**, itself a fixed rational diagonal.

No equality over \(\mathbb Z\) between
\(C_{n-1,p-1}\) and \(b_{n-p}\) is used; only the congruence modulo the
defining prime is asserted.

For symbolic elimination there is also the explicit one-sum formula. For
\(0\le r\le m\),
\[
\boxed{
[x^ry^rz^rw^r]F^m
=
\sum_{j=\max(0,\,2r-m)}^r
(-1)^j\binom mj
\left(
\frac{(m-j)!}
{(r-j)!^2(m-2r+j)!}
\right)^2.
}
\tag{B.8}
\]
So the array \(C_{s,m}\), with \(r=s-m\), is a proper-hypergeometric
two-parameter object suitable for Ore/creative-telescoping elimination.

I searched the repository and did not find this exact anti-diagonal carrier
already recorded.

## B2. Reflection gives a second row, but not a second independent congruence

For \(p=m+1\) prime, the proved reflection
\[
b_r\equiv b_{p-1-r}\pmod p
\]
becomes
\[
\boxed{
C_{s,m}\equiv C_{3m-s,m}\pmod{m+1}.
}
\tag{B.9}
\]

At \(s=n-1\), a bad column therefore satisfies
\[
m+1\mid C_{n-1,m},
\qquad
m+1\mid C_{3m-n+1,m}.
\tag{B.10}
\]

This is **not** codimension two: the second congruence is reflection of the
first. Its value is that the two integers lie in the same
characteristic-zero array, so we can form cross-column determinants before
reducing.

Define
\[
V_n(m)=
\binom{C_{n-1,m}}{C_{3m-n+1,m}}.
\tag{B.11}
\]

The new object to seek is a determinant/resultant involving
\(V_n(m),V_n(m-h),V_n(m-h-k)\), not a stronger statistic of one \(Z_p\).

## B3. A positive-proportion bad set forces a short three-prime cluster

Let
\[
R_n=\{r<n/2:n-r\text{ prime},\ n-r\mid b_r\},
\qquad M=|R_n|=N_1(n),
\]
and write
\[
0\le r_1<\cdots<r_M<n/2.
\]

Assume, for some fixed \(\varepsilon>0\),
\[
M\ge\varepsilon\frac n{\log n}.
\tag{B.12}
\]
Put
\[
H=\left\lceil\frac{4\log n}{\varepsilon}\right\rceil.
\tag{B.13}
\]
Since the total adjacent-gap length is \(<n/2\), fewer than
\[
\frac{n}{2H}\le\frac M8
\]
adjacent gaps exceed \(H\). A long gap destroys at most two consecutive
triples, so at least \(M/2\) consecutive triples have both gaps at most
\(H\), for large \(n\).

There are at most \(H^2\) ordered gap pairs. Hence some
\(1\le h,k\le H\) occurs at least
\[
\frac{M}{2H^2}
\gg_\varepsilon
\frac n{\log^3n}
\tag{B.14}
\]
times.

Equivalently, for \(\gg_\varepsilon n/\log^3 n\) starting primes,
\[
p,\quad p-h,\quad p-h-k
\tag{B.15}
\]
are all prime and all bad for the **same \(n\)**.

In terms of \(m=p-1\), each occurrence gives the exact two-row,
three-characteristic packet
\[
\begin{array}{c|ccc}
&m&m-h&m-h-k\\ \hline
\text{top}
&C_{n-1,m}
&C_{n-1,m-h}
&C_{n-1,m-h-k}\\
\text{reflected}
&C_{3m-n+1,m}
&C_{3m-n+1-3h,m-h}
&C_{3m-n+1-3h-3k,m-h-k}.
\end{array}
\tag{B.16}
\]
The three columns are divisible by the three distinct characteristics
\[
m+1,\quad m-h+1,\quad m-h-k+1,
\]
respectively.

This is exactly the shape demanded by the current interface no-go:
two rows, three characteristics, now with only logarithmic inter-column
separation.

## B4. The narrowly scoped missing lemma

Let \(\mathcal C_n(h,k)\) count primes \(p\in(n/2,n]\) such that
\[
p,\ p-h,\ p-h-k
\]
are prime and
\[
p\mid b_{n-p},\quad
p-h\mid b_{n-p+h},\quad
p-h-k\mid b_{n-p+h+k}.
\tag{B.17}
\]

The following suffices for q=1.

> **(L3C) Logarithmic three-characteristic codegree.**
> For every fixed \(A>0\),
> \[
> \max_{1\le h,k\le A\log n}
> \mathcal C_n(h,k)
> =
> o_A\left(\frac n{\log^3 n}\right).
> \tag{B.18}
> \]

Indeed, (B.12) produces \(h,k\le(4/\varepsilon+o(1))\log n\) violating
(B.18). Since \(\varepsilon\) is arbitrary, L3C gives
\[
N_1(n)=o(n/\log n).
\]

A fixed power saving
\[
\mathcal C_n(h,k)\ll_A n^{1-\delta_A}
\tag{B.19}
\]
would be much stronger than necessary.

This is a strict shrinkage of the target: no global \(Z_p\) theorem, no
high factorial moment, no averaging over \(n\), and no long gaps. The whole
load-bearing problem is reduced to a fixed logarithmic gap pair.

## B5. The first obvious “second equation” fails exactly

One might hope that reflection gives an extra vanishing from the existing
continuant gap polynomial \(N_d(x)\). It does not.

For a reflection pair
\[
x,\quad p-1-x,
\]
choose the smaller \(x\), and let
\[
d=p-1-2x.
\]
The proved symmetry
\[
N_d(-X-d-1)=(-1)^{d-1}N_d(X)
\tag{B.20}
\]
has \(d\) even. Therefore \(N_d\) contains the structural central factor
\[
2X+d+1.
\tag{B.21}
\]
At \(X=x\), that factor equals \(p\). Thus
\[
p\mid N_d(x)
\]
is automatic for **every reflected projective pair**; it does not know that
the distinguished coordinate \(b_x\) vanishes.

The smallest tempting repair is
\[
p\mid b_x
\Longrightarrow
p\mid\frac{N_d(x)}{2x+d+1},
\tag{B.22}
\]
equivalently
\[
p^2\mid N_d(x).
\tag{B.23}
\]

I tested (B.23) with exact integer modular arithmetic for every prime
\(7\le p\le10000\):

\[
\begin{array}{c|r}
\text{active primes }(Z_p\ne\varnothing)&470\\
\text{noncentral reflection zero-pairs}&585\\
\text{central zeros}&2\\
\text{noncentral pairs with }p^2\mid N_d(x)&\mathbf 0.
\end{array}
\]

The first example is
\[
p=17,\quad x=3,\quad d=10,
\qquad
N_{10}(3)\equiv34\pmod{17^2}.
\]
After removing the forced factor \(17\), the quotient is \(2\pmod{17}\).

So the deflated-reflection/first-jet implication is not just unproved; it
is empirically false in every tested noncentral target. This agrees with the
repository's previous “no automatic second vanishing congruence” audits.

## B6. What would actually repair the route

The repair must couple **different columns before reduction**.

The first object I would compute is
\[
\Delta_h(n,m)
=
\det\bigl(V_n(m),V_n(m-h)\bigr),
\tag{B.24}
\]
followed by a three-column elimination
\[
\mathcal R_{h,k}(n,m)
\quad\text{built from}\quad
V_n(m),V_n(m-h),V_n(m-h-k).
\tag{B.25}
\]

Badness of two columns gives divisibility by two distinct defining primes of
a common integral determinant; the third column adds a second adjacent minor.
That is genuinely cross-characteristic.

The decisive question is not whether such determinants exist; they do.
It is whether the explicit hypergeometric/Ore structure (B.8) forces a
**nonzero primitive** minor/resultant whose algebraic complexity is only
\(\operatorname{poly}(h,k)\), or at worst \(n^{o(1)}\) when
\(h,k=O(\log n)\), rather than paying \(\Theta(n)\) logarithmic height per
row.

This matches the strongest surviving algebraic class already isolated in
the repository's scale-sensitive cross-row audit: use adjacent/banded
identities **before** multiplying one exponential integer carrier per row.

My concrete next symbolic program would be:

1. derive a bivariate Ore system for \(C_{s,m}\) from (B.8);
2. restrict it to \(s=n-1\) and \(s=3m-n+1\);
3. eliminate the hypergeometric certificate across
   \(m,m-h,m-h-k\);
4. primitive-saturate the resulting \(2\times2\) and \(2\times3\) minors;
5. measure their degree and height as functions of \(n,h,k\).

If, even for fixed \(h,k\), every primitive minor has logarithmic height
\(\Theta(n)\), this route dies cleanly by the existing height-capacity
barrier. I have **not** proved the desired low-height minor; that is the
exact missing premise.

---

# Numerical test of the narrowed target

I computed the exact q=1 incidence table from the recurrence for every
\(n\le20000\). The distribution is

\[
\begin{array}{c|rrrr}
N_1(n)&0&1&2&3\\ \hline
\#\{n\le20000\}&18486&1455&56&3.
\end{array}
\]

The maximum \(3\) occurs at
\[
n=321,\quad11576,\quad18444.
\]

For \(n=321\),
\[
r=(110,128,142),\qquad
p=(211,193,179),
\]
with short gaps
\[
(h,k)=(18,14).
\]
This is the same three-prime packet visible in the repository's
\(n=321\) horizontal Mellin audit.

Taking
\[
H(n)=\lceil8\log n\rceil,
\]
the other two triple-hit indices have no triple with both adjacent gaps
\(\le H(n)\). Across the full scan, the maximum multiplicity of one fixed
short pair \((h,k)\) at one fixed \(n\) is **1**.

This is useful falsification evidence for L3C, not asymptotic evidence.

Here is a standalone Python script reproducing both the q=1 scan and the
deflated-reflection \(p^2\) test.

~~~python
from math import log

def primes_upto(N):
    isprime = bytearray(b"\x01") * (N + 1)
    if N >= 0:
        isprime[0] = 0
    if N >= 1:
        isprime[1] = 0
    for p in range(2, int(N**0.5) + 1):
        if isprime[p]:
            isprime[p*p:N+1:p] = b"\x00" * (((N - p*p)//p) + 1)
    return [p for p in range(2, N + 1) if isprime[p]]

def P(n, mod):
    return ((2*n + 1) * (17*n*n + 17*n + 5)) % mod

def zero_row(p):
    # r in [0,p-1] with b_r == 0 mod p
    if p < 7:
        return []
    z = []
    bm1, b = 1 % p, 5 % p
    if bm1 == 0:
        z.append(0)
    if b == 0:
        z.append(1)
    for m in range(1, p - 1):
        rhs = (P(m, p) * b - pow(m, 3, p) * bm1) % p
        bn = rhs * pow(pow(m + 1, 3, p), -1, p) % p
        if bn == 0:
            z.append(m + 1)
        bm1, b = b, bn
    return z

def q1_scan(N=20000):
    hits = [[] for _ in range(N + 1)]
    for p in primes_upto(N):
        if p < 7:
            continue
        for r in zero_row(p):
            n = p + r
            if n <= N and p > n/2:
                hits[n].append(r)

    hist = {}
    maxM = 0
    argmax = []
    max_short_codegree = 0
    max_record = None

    for n in range(1, N + 1):
        rs = sorted(hits[n])
        M = len(rs)
        hist[M] = hist.get(M, 0) + 1
        if M > maxM:
            maxM, argmax = M, [n]
        elif M == maxM and M:
            argmax.append(n)

        H = int(8 * log(max(n, 3)) + 0.999999)
        codeg = {}
        for i in range(M):
            for j in range(i + 1, M):
                h = rs[j] - rs[i]
                if h > H:
                    break
                for k0 in range(j + 1, M):
                    k = rs[k0] - rs[j]
                    if k > H:
                        break
                    key = (h, k)
                    codeg[key] = codeg.get(key, 0) + 1

        for key, v in codeg.items():
            if v > max_short_codegree:
                max_short_codegree = v
                max_record = (n, key, rs, H)

    print("hist =", hist)
    print("max N1 =", maxM, "at", argmax)
    print("max short (h,k)-codegree =", max_short_codegree)
    print("record =", max_record)

# N_1(x)=1, N_2(x)=P(x+1),
# N_{h+1}(x)=P(x+h)N_h(x)-(x+h)^6 N_{h-1}(x).
def N_gap_mod(x, h, mod):
    if h == 1:
        return 1 % mod
    a = 1 % mod
    b = P(x + 1, mod)
    if h == 2:
        return b
    for d in range(2, h):
        c = (P(x + d, mod) * b - pow(x + d, 6, mod) * a) % mod
        a, b = b, c
    return b

def reflected_p2_scan(PMAX=10000):
    active = pairs = central = p2 = 0
    first_nonzero = []
    for p in primes_upto(PMAX):
        if p < 7:
            continue
        Z = set(zero_row(p))
        if Z:
            active += 1
        for r in sorted(Z):
            rr = p - 1 - r
            if r == rr:
                central += 1
                continue
            if r > rr:
                continue
            assert rr in Z
            pairs += 1
            d = rr - r
            v = N_gap_mod(r, d, p*p)
            if v == 0:
                p2 += 1
            elif len(first_nonzero) < 8:
                first_nonzero.append((p, r, d, v))
    print("active primes =", active)
    print("noncentral reflection pairs =", pairs)
    print("central zeros =", central)
    print("p^2 | N_d(r) =", p2)
    print("first nonzero residues =", first_nonzero)

if __name__ == "__main__":
    q1_scan(20000)
    reflected_p2_scan(10000)
~~~

Expected key output:

~~~text
hist = {0: 18486, 1: 1455, 2: 56, 3: 3}
max N1 = 3 at [321, 11576, 18444]
max short (h,k)-codegree = 1
active primes = 470
noncentral reflection pairs = 585
central zeros = 2
p^2 | N_d(r) = 0
~~~

---

# Part C. Adversarial audit

## C1. Reflection is not independent

This is the biggest possible false breakthrough. Equation (B.9) is the same
reflection congruence already banked. Treating top and bottom as two
independent equations would violate the interface no-go. The second row is
useful only if an **integral cross-column determinant/syzygy** uses both rows
before reducing modulo three characteristics.

## C2. The common carrier \(H_n\) is not itself a sieve

Writing all events as
\[
p\mid[T^{p-1}]H_n
\]
does not make self-divisibility rare. D-finiteness/rational-diagonal
structure alone does not control coefficients at their own varying prime
indices. The carrier matters only because it supplies one characteristic-zero
object on which a three-column identity could act.

## C3. No integer equality was smuggled in

Only
\[
C_{n-1,p-1}\equiv b_{n-p}\pmod p
\]
is used. Replacing this congruence by an equality would be circular and false.

## C4. The clustering argument is only at contradiction scale

The pigeonhole lemma assumes
\[
N_1(n)\ge\varepsilon n/\log n.
\]
It does not say typical bad primes cluster. To get little-\(o\), L3C must hold
for every fixed \(A\), with constants allowed to depend on \(A\). There is no
hidden uniformity in \(\varepsilon\).

## C5. No multiplicity loss

\(N_1(n)\) counts distinct primes, so \(p\)-adic valuation multiplicities are
irrelevant here. Since \(p\asymp n\), prime count and logarithmic mass are
equivalent at the scale of the target. Central reflection zeros must be
treated separately; I excluded them from the \(p^2\) test.

## C6. Where this escapes -- and could fall back into -- the interface no-go

The proved carrier uses the **same integer array** \(C_{s,m}\) for all
characteristics. An arbitrary adversarial family of independent \(Z_p\)'s
need not lift to such an array.

But if the next step uses only row cardinalities, reflection, per-column gap
bounds, or separate resultants, it immediately collapses back to the old
interface. The genuinely new datum must be a syzygy/minor involving at least
the three nearby \(m\)-columns simultaneously.

## C7. The load-bearing algebra is still unproved

I have not derived the minimal Ore ideal of \(C_{s,m}\), factored
\(\Delta_h(n,m)\), or proved that a three-column primitive minor has
\(\operatorname{poly}(h,k)\) height. This is the exact open step. It is
possible that every such minor still has \(\Theta(n)\) logarithmic height
for fixed \(h,k\); that would kill the route.

## C8. Numerical and literature limits

The q=1 scan stops at \(n=20000\); the reflected first-jet audit stops at
\(p=10000\). They are falsification tests, not asymptotic evidence. I did
not test q>1 here.

Likewise, a literature search cannot prove nonexistence of a theorem. My
claim is narrower: among the named families and recent relevant work through
September 2026, I found no result with both the correct recurrence class and
the correct worst-case proportional-window quantifiers.

---

# Recommendation for R2

Do **not** spend R2 on another marginal model for \(Z_p\). Derive the
Ore/creative-telescoping relations for (B.8) and ask whether the two-row
vectors \(V_n(m)\) admit a low-complexity Plücker/Desnanot identity across
\(m,m-h,m-h-k\).

The first binary test should be:

\[
\boxed{
\text{For fixed }h,k,\ 
\text{does the primitive three-column elimination have }
\log\operatorname{height}=O_{h,k}(\log n)
\text{ (or }n^{o(1)}\text{), rather than }\Theta(n)?
}
\]

A positive answer would be the first recurrence-specific
cross-characteristic input I can see that escapes the current interface
no-go. A negative answer would kill this attack cleanly and early.
