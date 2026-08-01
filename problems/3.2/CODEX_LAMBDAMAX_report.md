# Direct power-iteration measurement of the centered Gram matrix

## 1. Histograms and finite-field arithmetic

For each prime, `M=p-2`, `I_h={1,...,M-h}`, and `H=ceil(p^(2/3))` was evaluated as the least integer with `H^3>=p^2`. The Apéry recurrence produced `b_n,c_n`; every table value used `Delta_(r,h)=b_r c_(r+h)-b_(r+h)c_r (mod p)`. Simultaneously, the numerator recurrence checked `Delta_(r,h) prod_(j=1)^h(r+j)^3=N_h(r)` at every histogram entry. Counts and all displayed integer invariants are exact; only eigendata and displayed ratios use binary64 arithmetic.

| p | H | sum_h card(I_h) | D1 checks | max n_h(a) | max n_1(a) |
|---:|---:|---:|---:|---:|---:|
| 997 | 100 | 94450 | 94450 | 7 | 3 |
| 1999 | 159 | 304803 | 304803 | 9 | 3 |
| 4001 | 253 | 979616 | 979616 | 8 | 1 |

All D1 checks passed. The independent identity `Delta_(r,1)=(r+1)^(-3)` also passed at every admissible `r` for all three primes.

## 2. Numerical Parseval check at p=199

Using `H=35`, all 1225 ordered row pairs were checked against FFT evaluations of `S_h(t)=sum_a n_h(a) exp(-2 pi i t a/p)`. The comparison was

`p <q_h,q_k> = p sum_a n_h(a)n_k(a)-|I_h||I_k| = sum_(t!=0) S_h(t) conjugate(S_k(t))`.

The maximum real-part error was `2.183e-11` and the maximum spurious imaginary part was `4.396e-12`; the largest compared exact entry had magnitude `77800`. All 6265 underlying D1 checks passed.

## 3. Matrix-free iteration and exact checks

No campaign-size Gram matrix was formed. With `Q_(h,a)=q_h(a)`, every application of `Gamma/p^2` used exactly the two table passes `v -> Q^T v -> Q(Q^T v)/p`. For each block, one seeded random start was run for exactly 40 steps; its Rayleigh quotient gave `lambda_1`. The same fixed procedure on the orthogonal complement of its vector gave `lambda_2`. There was no convergence-based early stop.

The exact checks use `Gamma_hh=p sum_a n_h(a)^2-|I_h|^2`, `Tr(Gamma)=sum_h Gamma_hh`, and

`1^T Gamma 1 = p ||sum_h q_h||^2 = p N_coinc-S^2`,

where `S=sum_h |I_h|` and `N_coinc=sum_a(sum_h n_h(a))^2`. The next table gives the exact integer values; `(h=1 removed)` recomputes every quantity on the principal block indexed by `2,...,H`. The norm identity was checked a second exact way as `sum_a(p sum_h n_h(a)-S)^2/p`; it agreed in every block.

| p | block | S | N_coinc | D_max (argmax h) | 1^T Gamma 1 | Tr(Gamma) |
|---:|:---|---:|---:|---:|---:|---:|
| 997 | all h | 94450 | 9051178 | 1981030 (h=1) | 103221966 | 94383966 |
| 997 | h=1 removed | 93456 | 8863272 | 1066466 (h=16) | 102658248 | 92402936 |
| 1999 | all h | 304803 | 46811323 | 7978000 (h=1) | 670965868 | 611529304 |
| 1999 | h=1 removed | 302807 | 46198471 | 4313806 (h=4) | 658664280 | 603551304 |
| 4001 | all h | 979616 | 240920516 | 16824084 (h=9) | 4275477060 | 3895811424 |
| 4001 | h=1 removed | 975618 | 238967006 | 16824084 (h=9) | 4276509082 | 3895799430 |

## 4. Top two eigenvalues and normalized diagnostics

Here `rows` is `H` for the full block and `H-1` after peeling. Thus the all-ones column is normalized by the actual number of retained rows. Residuals are for the normalized matrix `Gamma/p^2`.

| p | block | rows | 40-step power lambda_1/p^2 | 40-step deflated lambda_2/p^2 | D_max/p^2 | (1^T Gamma 1)/(p^2 rows) | Tr/p^2 | residuals (1,2) |
|---:|:---|---:|---:|---:|---:|---:|---:|:---|
| 997 | all h | 100 | 2.41528009 | 1.85490387 | 1.99296988 | 1.03844096 | 94.95282839 | 2.20e-05, 1.49e-02 |
| 997 | h=1 removed | 99 | 1.90078950 | 1.85503677 | 1.07289371 | 1.04320182 | 92.95985851 | 1.29e-02, 2.76e-02 |
| 1999 | all h | 159 | 2.30635333 | 1.80103581 | 1.99649600 | 1.05603355 | 153.03532306 | 5.66e-05, 2.65e-02 |
| 1999 | h=1 removed | 158 | 1.79773356 | 1.77476494 | 1.07953076 | 1.04323329 | 151.03882707 | 1.63e-02, 3.91e-02 |
| 4001 | all h | 253 | 1.74444143 | 1.71068284 | 1.05097969 | 1.05566703 | 243.36651553 | 3.02e-02, 2.84e-02 |
| 4001 | h=1 removed | 252 | 1.76265094 | 1.72575512 | 1.05097969 | 1.06011201 | 243.36576628 | 8.16e-03, 8.88e-03 |

The largest measured deflation overlap `|<u_1,u_2>|` was `3.963e-17`.

Because 40 scalar power steps need not resolve a clustered spectral edge, the same matrix-free operator was independently audited with 40-step fully reorthogonalized symmetric Lanczos. Only the resulting 40 by 40 tridiagonal matrix was diagonalized. These are the preferred spectral estimates; in particular, they restore the required principal-submatrix interlacing in the `p=4001` cell.

| p | block | audited lambda_1/p^2 | audited lambda_2/p^2 | power shortfall (1,2) |
|---:|:---|---:|---:|:---|
| 997 | all h | 2.41528009 | 1.89064138 | 9.067e-10, 3.574e-02 |
| 997 | h=1 removed | 1.90409636 | 1.88462855 | 3.307e-03, 2.959e-02 |
| 1999 | all h | 2.30635334 | 1.83140730 | 6.203e-09, 3.037e-02 |
| 1999 | h=1 removed | 1.83149552 | 1.80520680 | 3.376e-02, 3.044e-02 |
| 4001 | all h | 1.76339557 | 1.72687912 | 1.895e-02, 1.620e-02 |
| 4001 | h=1 removed | 1.76339483 | 1.72687820 | 7.439e-04, 1.123e-03 |

## 5. Off-diagonal operator

For `O=Gamma-diag(Gamma)`, one seeded random start was run for exactly 40 symmetric-matrix power steps. The Rayleigh-quotient sign was retained; this distinguishes a dominant negative eigenvalue from a positive one. Residuals refer to `O/p^2`.

| p | block | signed power Rayleigh/p^2 | 40-step abs(Rayleigh)/p^2 | 40-step random-scale ratio | residual |
|---:|:---|---:|---:|---:|---:|
| 997 | all h | +1.04303156 | 1.04303156 | 3.29340415 | 1.61e-03 |
| 997 | h=1 removed | +0.94240408 | 0.94240408 | 2.99066069 | 1.27e-02 |
| 1999 | all h | +0.87232307 | 0.87232307 | 3.09303757 | 9.05e-03 |
| 1999 | h=1 removed | +0.84434117 | 0.84434117 | 3.00327998 | 1.64e-02 |
| 4001 | all h | +0.77353686 | 0.77353686 | 3.07613244 | 1.94e-02 |
| 4001 | h=1 removed | +0.78249122 | 0.78249122 | 3.11790932 | 1.73e-02 |

The symmetric Lanczos audit also computed both algebraic edges of `O/p^2`; the positive edge dominates in every cell.

| p | block | lambda_max(O)/p^2 | lambda_min(O)/p^2 | audited norm(O)/p^2 | audited random-scale ratio |
|---:|:---|---:|---:|---:|---:|
| 997 | all h | +1.04305563 | -0.67034029 | 1.04305563 | 3.29348015 |
| 997 | h=1 removed | +0.95751923 | -0.65152687 | 0.95751923 | 3.03862764 |
| 1999 | all h | +0.87529359 | -0.62905872 | 0.87529359 | 3.10357027 |
| 1999 | h=1 removed | +0.86453408 | -0.60482083 | 0.86453408 | 3.07510515 |
| 4001 | all h | +0.78878873 | -0.54692936 | 0.78878873 | 3.13678471 |
| 4001 | h=1 removed | +0.78878699 | -0.55661167 | 0.78878699 | 3.14299540 |

## 6. Verdict

**[OP-OFF-0] finite test: supported at these three sizes.** The full `lambda_1/p^2` values lie in `[1.7634, 2.4153]`; after removing the exceptional cube row they lie in `[1.7634, 1.9041]`. There is no observed growth with `p` in this range. This is direct finite evidence for the bounded scenario, not an asymptotic proof.

**Random-matrix off-diagonal scale: supported.** The ratios `||O||/(sqrt(rows) p^(3/2))` lie in `[3.1036, 3.2935]` for the full blocks and `[3.0386, 3.1430]` after peeling. They remain constant-sized, so the data fit the proposed `sqrt(H) p^(3/2)` scale. Peeling changes the cube-row outlier in the PSD spectrum for primes `p=1 mod 3`, but does not expose a growing off-diagonal norm.

The scalar-power rows are finite-iteration lower estimates in absolute Rayleigh quotient. The displayed residuals and shortfalls quantify their unresolved spectral error; the numerical ranges in this verdict use the fixed-step symmetric Lanczos audit.
