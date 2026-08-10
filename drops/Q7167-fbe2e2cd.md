ANSWER Q7167 fbe2e2cd

# Audit of P32 periodic extremizer R5

## Result

PASS. The displayed branches are numerically consistent with the actual cap \(A_Z(h)\le 3(h-1)\). Periodicity produces repeated differences, but the repetitions are below the linear cap at the relevant scales.

For \(N=mR\), \(R=m/100\), and
\[
S=\sum_{i=0}^{m-1}(3m+i)=m(7m-1)/2,
\]
we have
\[
L=RS=\frac{m^2(7m-1)}{200}.
\]
Hence \(4L+1\) is cubic in \(m\), and any chosen prime \(p>4L+1\) satisfies \(p>H^2\) for \(H=12m\) and sufficiently large \(m\).

## Geometry

The reflection
\[
x\mapsto p-1-x
\]
fixes \(Z\). All internal gaps are
\[
3m+(i\bmod m)\ge3m>1,
\]
so there is no adjacency. Three consecutive gaps have total length at most
\[
3(4m-1)=12m-3<H.
\]
After removing the two endpoint windows in each half, the number of short consecutive four-point windows is
\[
2(N-2)=\Omega(H^2).
\]
They are off-center because their span is \(<H\) and they lie entirely in one reflected half.

## Pair-cap cases

### Non-period differences

For a within-half difference with block length \(\ell=qm+r\), \(0<r<m\), the value is
\[
qS+F_r(a).
\]
The partial-period sums have at most two starts for a fixed target. The number of possible remainders is bounded by the growth of \(F_r\), giving the stated estimate
\[
W(h)\le 2(R+1)(h/(12m)+2).
\]
With \(R=m/100\), this is much smaller than \(3(h-1)\) throughout the admissible range.

### Full-period differences

If \(h=qS\), the only extra repetition comes from period shifts. There are at most
\[
N+1=mR+1
\]
starts. Since
\[
N+1=O(m^2),\qquad S=\Theta(m^2),
\]
and \(h\ge S\), we get
\[
N+1<3(h-1).
\]

### Cross-half differences

Cross-half pairs satisfy an equation of the form
\[
P_a+P_b=p-1-h.
\]
The stated bound
\[
W_{cross}(h)\le N+1
\]
is sufficient. These differences lie in the separate cross-half range \(>2L\), so they cannot mix with within-half pairs. Again \(N+1=O(m^2)\) is below the cap because these \(h\)'s are of order \(L\).

## Final audit

No displayed inequality has a counterexample. The earlier objection from raw periodic multiplicity compared \(O(mR)\) to a constant instead of to \(3(h-1)\). At the relevant scales the cap is much larger. The construction passes the requested checks: pair cap, palindrome, no adjacency, \(p>H^2\), short-window count, maximum span, and off-center condition.