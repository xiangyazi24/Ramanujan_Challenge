ANSWER Q7901 eed6cfc4

# Audit of the paired real-kernel formula and ACF3

## 1. Exact reflection factorization

For odd p, reflection invariance gives

$$r\in Z_p\iff p-1-r\in Z_p.$$

Pair the terms r and p-1-r. With

$$F_p(a)=\sum_{r\in Z_p}e_p(ar),$$

we have

$$e_p(ar)+e_p(a(p-1-r))=e_p(ar)+e_p(-a-a r)=e_p(-a/2)(e_p(a(r+1/2))+e_p(-a(r+1/2))).$$

Therefore

$$F_p(a)=e_p(-a/2)R_p(a),$$

where

$$R_p(a)=\sum_{r\in Z_p}e_p(a(r+1/2))$$

is real. Reflection also gives R_p(-a)=R_p(a).

The product phase in the cubic term is

$$e_p(-h\overline{q\ell}/2)e_q(-h\overline{p\ell}/2)e_\ell(-h\overline{pq}/2).$$

Because inverses are taken modulo the respective primes, this combines to the CRT half-shift phase

$$e_Q(-h/2),$$

not a cancellation. It remains coupled to the kernel.

## 2. Pairing h and -h

For

$$W_M(x)=\sum_{m\le M}e(mx)=e((M+1)x/2)\frac{\sin(\pi Mx)}{\sin(\pi x)},$$

one must keep the phase convention. Pairing h and Q-h gives the real contribution

$$Q^{-1}\sum_{1\le h\le(Q-1)/2,(h,Q)=1} \frac{\sin(2\pi M h/Q)}{\sin(\pi h/Q)}\times \text{(phase sign)}\times \prod R_p.$$

The claimed factor $(-1)^h$ is not automatic from reflection alone. It depends on replacing the original Fourier kernel by a centered interval with an integer half-shift. For general M=X^2 the exact sign must be recomputed from the chosen W_M normalization. Reflection removes the irrational-looking half-shift but does not by itself create $(-1)^h$.

## 3. Can oscillation plus Parseval save ACF3?

The available information is

$$|Z_p|\ll p^{2/3},$$

and

$$\sum_a |R_p(a)|^2=p|Z_p|.$$

Applying Parseval independently gives only

$$\sum_{a}|R_p(a)|^2\ll X^{5/3}.$$

For three primes the diagonal bound gives a size comparable to

$$\prod X^{5/6}=X^{5/2}$$

for square-root-normalized transforms, which is not enough to beat the required tuple scale after the number of prime triples is included.

The h-kernel has size

$$\left|\frac{\sin(2\pi Mh/Q)}{\sin(\pi h/Q)}\right|\ll \min(M,Q/h).$$

Since Q~X^3 and M=X^2, the range h\lesssim X contributes a kernel of size X^2. Summing this range with only Cauchy-Schwarz loses the full cancellation factor needed: the h-sum is too large by roughly a factor X.

## 4. Reflection-symmetric obstruction model

Take each R_p to be concentrated on an arbitrary symmetric set of frequencies of size p^{1/3} with amplitudes chosen to saturate Parseval. Then

$$\|R_p\|_2^2\asymp p|Z_p|$$

but the supports can align with the large-kernel h range. The product correlation then has size comparable to the Cauchy bound. Reflection and no-adjacent-zero conditions do not prevent this frequency alignment.

Hence Parseval plus reflection cannot prove a power saving.

## 5. Minimal missing estimate

A sufficient replacement is a genuinely mixed inequality for centered transforms:

$$\sum_{p,q,\ell\sim X}\left|\sum_{h\le Q/2}K_M(h,Q)R_p(h\bar{q\ell})R_q(h\bar{p\ell})R_\ell(h\bar{pq})\right|
\ll X^{2+o(1)}\lambda_X^3.$$

The required input is therefore not a one-prime energy estimate but a trilinear dispersion estimate controlling frequency alignment across different primes.

## Conclusion

The paired real-kernel reduction is valid only after fixing the exact Fourier kernel convention; the $(-1)^h$ factor requires an additional normalization check. More importantly, the oscillatory kernel together with reflection, spacing, cardinality, and Parseval does not force the ACF3 power saving. A new mixed-prime trilinear estimate is genuinely required.