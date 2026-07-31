# Q3.2: Zero Fiber Theory and the Horizontal Decorrelation Problem

**Date**: 2026-07-29  
**Status**: |Z_p| boundedness is INSUFFICIENT and WRONG; need horizontal decorrelation  
**Key revision**: Q5523 (ChatGPT, 07-29) demolished the |Z_p|≤D route

---

## 1. Setup and Definitions

The Apéry numbers for ζ(3):
$$b_n = \sum_{k=0}^{n} \binom{n}{k}^2 \binom{n+k}{k}^2$$

For each prime p, define the **zero fiber**:
$$Z_p = \{ r \in \{0, 1, \ldots, p-1\} : b_r \equiv 0 \pmod{p} \}$$

The **target primes** for a given n:
$$T_n = \{ p \in (n/2, n] : p \text{ prime}, \; p \mid b_{n-p} \}$$

The **radical**: $R_n = \prod_{p \in T_n} p$.

**P3.2 asks**: prove $\log G_n = o(n)$ for all $n$, where $G_n = \gcd(d_n a_n, d_n b_n)$.

Since $\log G_n \leq \log R_n + O(\log n)$, it suffices to prove $\log R_n = o(n)$.

---

## 2. Reflection Symmetry (PROVED)

**Theorem**. For every prime $p \geq 3$ and every $r \in \{0, \ldots, p-1\}$:
$$b_{p-1-r} \equiv b_r \pmod{p}$$

**Proof**. We use three binomial coefficient identities mod p.

**Step 1**: For $0 \leq k \leq p-1-r$:
$$\binom{p-1-r}{k} \equiv (-1)^k \binom{r+k}{k} \pmod{p}$$

*Proof of Step 1*: 
$$\binom{p-1-r}{k} = \frac{(p-1-r)(p-2-r)\cdots(p-r-k)}{k!} \equiv \frac{(-1)^k(r+1)(r+2)\cdots(r+k)}{k!} = (-1)^k \binom{r+k}{k} \pmod{p}$$

**Step 2**: For $0 \leq k \leq r$:
$$\binom{p-1-r+k}{k} \equiv (-1)^k \binom{r}{k} \pmod{p}$$

*Proof of Step 2*:
$$\binom{p-1-r+k}{k} = \frac{(p-1-r+k)(p-2-r+k)\cdots(p-r)}{k!} \equiv \frac{(-1)^k(r-k+1)(r-k+2)\cdots r}{k!} = (-1)^k \binom{r}{k} \pmod{p}$$

**Step 3**: For $k > r$ (with $k \leq p-1-r$ so that $p-1-r+k < 2p$):
$$\binom{p-1-r+k}{k} \equiv 0 \pmod{p}$$

*Proof of Step 3*: Write $p-1-r+k = p + (k-r-1)$ in base $p$: digits $(k-r-1, 1)$. And $k$ in base $p$: digits $(k, 0)$. By Lucas's theorem, $\binom{p-1-r+k}{k} \equiv \binom{1}{0} \cdot \binom{k-r-1}{k} \pmod{p}$. Since $k > r \geq 0$, we have $k > k-r-1 \geq 0$, so $\binom{k-r-1}{k} = 0$. ∎

**Combining** (for $r \leq (p-1)/2$, so $p-1-r \geq r$):

$$b_{p-1-r} = \sum_{k=0}^{p-1-r} \binom{p-1-r}{k}^2 \binom{p-1-r+k}{k}^2$$

Split at $k = r$: terms with $k > r$ vanish by Step 3. For $k \leq r$, apply Steps 1 and 2:

$$b_{p-1-r} \equiv \sum_{k=0}^{r} \left[(-1)^k \binom{r+k}{k}\right]^2 \left[(-1)^k \binom{r}{k}\right]^2 = \sum_{k=0}^{r} \binom{r+k}{k}^2 \binom{r}{k}^2 = b_r \pmod{p}$$

(The $(-1)^{2k} = 1$ factors cancel.)

For $r > (p-1)/2$: let $s = p-1-r < (p-1)/2$. By the case above, $b_{p-1-s} \equiv b_s \pmod{p}$, i.e., $b_r \equiv b_{p-1-r} \pmod{p}$. ∎

**Verified computationally**: all 94 primes $p \leq 500$, all 21,534 pairs $(p, r)$.

---

## 3. Consequences of Reflection

**Corollary 1**: $Z_p$ is invariant under $r \mapsto p-1-r$.

**Corollary 2**: $|Z_p|$ is even unless $(p-1)/2 \in Z_p$ (i.e., unless $p \mid b_{(p-1)/2}$).

**Corollary 3**: Define $Z_p^- = \{ s \in [0, (p-3)/2] : p \mid b_s \}$. Then
$$|Z_p| = 2|Z_p^-| + \varepsilon, \quad \varepsilon = \mathbf{1}[p \mid b_{(p-1)/2}]$$

---

## 4. Empirical Statistics (p ≤ 50,000; 5,132 primes)

| |Z_p| | Count | Frequency | Poisson(1/2) prediction |
|-------|-------|-----------|------------------------|
| 0     | 3109  | 0.6060    | 0.6065                 |
| 1     | 2     | 0.0004    | (midpoint correction)  |
| 2     | 1536  | 0.2994    | 0.3033                 |
| 4     | 416   | 0.0811    | 0.0758                 |
| 6     | 60    | 0.0117    | 0.0126                 |
| 8     | 9     | 0.0018    | 0.0016                 |

**Max |Z_p| = 8** (at 9 primes: 3727, 6841, 13463, 16451, 28933, 32063, 38197, 43151, 48109).

**|Z_p| = 1 primes**: Only p = 11 and p = 3137 (both at the midpoint $r = (p-1)/2$).

---

## 5. The Poisson(1/2) Model

The distribution of $|Z_p^-|$ is strikingly close to Poisson(1/2). Explanation:

For each $s \in [0, (p-3)/2]$, the event "$p \mid b_s$" requires $p$ to be a prime factor of the integer $b_s \approx (1+\sqrt{2})^{4s}$. Heuristically, this occurs with probability $\approx 1/p$.

The total:
$$E[|Z_p^-|] \approx \sum_{s=0}^{(p-3)/2} \frac{1}{p} \approx \frac{1}{2}$$

By the Poisson limit theorem (sum of many rare independent events), $|Z_p^-| \sim \text{Poisson}(1/2)$.

**Prediction**: $|Z_p| = 10$ should first appear near $p \sim 5 \times 10^5$ (since $P(|Z_p^-| \geq 5) \approx 1.6 \times 10^{-4}$ and $\pi(5 \times 10^5) \approx 41,538$, giving expected count $\approx 6.6$).

---

## 6. The Reduction: |Z_p| Bounded ⟹ P3.2

**Assume** $|Z_p| \leq D$ for all primes $p$ (absolute constant).

Then for each $n$:
$$|T_n| = \#\{ p \in (n/2, n] : p \mid b_{n-p} \} = \sum_{p \in (n/2,n]} \mathbf{1}[(n-p) \in Z_p]$$

Each prime $p$ contributes at most 1. The "expected" value:
$$E[|T_n|] \approx \sum_{p \in (n/2,n]} \frac{|Z_p|}{p} \leq D \sum_{p \in (n/2,n]} \frac{1}{p} \approx \frac{D \ln 2}{\ln n} \to 0$$

**Key insight**: Even without a hard bound on $|Z_p|$, the Poisson model gives $|T_n| \sim \text{Poisson}(\ln 2 / \ln n)$, which is almost surely 0 or 1 for large $n$.

---

## 7. P3.2 Under the Poisson Model (No Hard Bound Needed)

Under the Poisson model for $|Z_p^-|$:

1. $E[|Z_p|] \approx 1$ for all primes $p$
2. $E[|T_n|] \approx \ln 2 / \ln n \to 0$
3. $|T_n| \sim \text{Poisson}(\ln 2 / \ln n)$
4. $\max_{n \leq N} |T_n| = O(\ln N / \ln \ln N)$ (by Chernoff + union bound)
5. $\log R_n \leq |T_n| \cdot \log n = O((\log n)^2 / \log \log n) = o(n)$

**This proves P3.2** — assuming the Poisson model (or any model where $E[|Z_p|] = O(1)$).

---

## 8. What Remains for a Rigorous Proof

The Poisson model is a heuristic. **UPDATE (2026-07-29)**: Routes (A)-(D) below are ALL insufficient per se (see §12). The planted diagonal counterexample shows no per-prime estimate proves P3.2. The actual target is the pair energy $H_2$ or a cross-prime theorem.

**(A) Hard bound**: ~~Prove $|Z_p| \leq D$.~~ Even $|Z_p| = O(1)$ doesn't prove P3.2 (planted diagonal).

**(B) Average bound**: $\sum_{p \leq x} |Z_p| = O(x/\ln x)$ is true (from $|Z_p| \ll p^{2/3}$) but gives only $E_1 \ll N^{8/3}/\log^2 N$, NOT the needed $o(N^2/\log^2 N)$.

**(C) Large sieve**: Blocked by modulus-square barrier. Spacing $N^{-2}$ forces operator norm $\gg N^2$.

**(D) Direct counting**: This IS the correct target, but needs CROSS-PRIME information.

---

## 9. Connection to Known Results

- **Codegree amplification** (known): $\#\{n \leq N : \log G_n > \varepsilon n\} = O_\varepsilon((\log N)^2)$. The polylog exceptional set comes from a simple counting argument; our goal is to eliminate it entirely.

- **Beukers (1987)**: $b_{p-1} \equiv a_p \pmod{p}$ where $a_p$ is the $p$-th coefficient of a weight-4 modular form. This connects $|Z_p|$ to the arithmetic of the modular form.

- **Stienstra-Beukers (1985)**: The Apéry family is a family of K3 surfaces; the Picard-Fuchs equation has monodromy group $\subseteq \text{Sp}(4)$.

---

## 10. ⚡ LUCAS COLLAPSE (Q5546, 2026-07-29 major revision)

**Theorem.** For every prime $p \geq 7$ and $n = p + r$ with $0 \leq r < p$:
$$p \mid b_{n-p} \iff p \mid b_n$$

Proof: Gessel's Lucas congruence gives $b_{p+r} \equiv 5 b_r \pmod{p}$. Since $5$ is a unit mod $p \geq 7$, $b_r \equiv 0 \iff b_{p+r} \equiv 0$.

**Consequence**: ALL target primes divide ONE integer:
$$T_n = \{p \in (n/2, n] : p \mid b_n\}, \quad K_n = \omega_{(n/2,n]}(b_n)$$

This completely changes the problem's nature. It is NOT about cross-prime "alignment" of independent zero fibers. It is about the prime factorization of a single P-recursive integer $b_n$.

---

## 11. H₂ SUFFICES for P3.2 (Q5550)

**Theorem.** $H_2(N) = \sum_{N < n \leq 2N} (K_n)_2 = o(N^2/\log^2 N) \implies M_N = o(N/\log N)$.

Proof: $M_N \leq 1 + \sqrt{H_2(N)}$, so $H_2 = o(N^2/\log^2 N) \implies M_N = o(N/\log N)$.

Even $H_2 \ll N^{2-\delta}$ for any $\delta > 0$ suffices. No growing-$k$ moments needed.

**Current unconditional**: $E_1(N) = H_2(N) \ll N^{8/3}/\log^2 N$ (from $|Z_p| \ll p^{2/3}$). Gap: $N^{2/3}$.

---

## 12. ALL STANDARD METHODS FAIL (Q5546, Q5549)

### Planted diagonal counterexample
For any $n_0$, define $Z_p^* = \{n_0 - p, p-1-(n_0-p)\}$ (reflection pair). Then:
- $|Z_p^*| \leq 2$, exact reflection, bounded additive energy, good Fourier
- Yet $K^*(n_0) = \pi(n_0) - \pi(n_0/2) \asymp n_0/\log n_0$

**No per-prime estimate can prove P3.2.** This includes:
1. BV theorem: range mismatch (moduli $\sim x$ but BV only to $x^{1/2}$)
2. Large sieve: spacing $N^{-2}$, operator norm $\gg N^2/\log N$
3. Turán-Kubilius: $f(r) = 1_{Z_p}(r)$ not additive in $r$
4. CRT: interval $N$ < one period $pq \sim N^2$
5. Even $|Z_p| = O(1)$ still misses by $\log$ factor

### Spike-to-close-pairs reduction (Q5549 §6-7)
If $K_n \geq \varepsilon n/\log n$, Kummer pruning + pigeonhole forces $\Omega_\varepsilon(n)$ same-kernel close pairs with gap $\leq \Delta = n^{o(1)}$.

**Sufficient condition (CP/SK.5)**: $\mathcal{C}_\varepsilon(N) = o_\varepsilon(N) \implies$ P3.2.

---

## 13. ⛔ HEIGHT CONTRADICTION BLOCKED (Q5554)

The global relation $A = b_n/R$ with $A \equiv b_n/p_i \cdot (P'(p_i))^{-1} \cdot (-1)^{K-1} \pmod{p_i}$ is just CRT reconstruction.

1. **Pair moduli collapse**: $\text{lcm}(\text{pair moduli}) \leq R$, NOT $R^2$
2. **$R^2 > A$ is consistent**: CRT determines $A$ uniquely, no contradiction
3. **mod-$p^2$ jet parametrizes**: $K$ equations on $K$ local CRT components, no overdetermination
4. **Close-pair gap is a unit**: small $h = q - p$ doesn't reduce local entropy

**Critical threshold** (if small-lift representative existed at order $m$):
$$K_n \leq \frac{cn}{m \log(n/2)} + o(n/\log n)$$
At $m = 7$: $K \leq 0.50n/\log n$. But NO small representative exists yet.

---

## 14. CFVZ SQUARE FACTORIZATION (Q5555)

$A_p(t) = B_p(t)^2$ (pure classes) or $(t^2 - 34t + 1) B_p(t)^2$.

**Does NOT bound $|Z_p|$**: a reciprocal square can have $p - O(1)$ zero coefficients. Squareness constrains polynomial roots, not coefficient zeros. Cauchy-Davenport gives support bounds but not cancellation control.

**Genuine gain**: $B_p$ comes from a FIXED second-order recurrence independent of $p$:
$$(2n-1)^2 s_n = \tilde{P}(n) s_{n-1} - n^2 s_{n-2}$$
This is the symmetric square root of the Apéry operator.

---

## 15. VIABLE PATHS (as of 2026-07-29 18:00 CT)

### Path A: Cross-prime divided-Frobenius container
Seek invariant $D$ and integer $M_D(n)$ with: compression ($n^{o(1)}$ values), container ($p \mid M_D$), height ($o(n)$), nondegeneracy.
**Status**: No container found. Free quotient digits block Bezout.

### Path D: Holonomic GCD bypass ⭐
$\log \gcd(b_n, \binom{n}{\lfloor n/2 \rfloor}) = o(n)$ would prove top-half directly.

**Empirical data (n ≤ 2000)**:
- Max $\log(\gcd)/n$ by dyadic block: 0.25 → 0.025 (monotonically decreasing)
- $K_n$ distribution: 89.9% = 0, 9.6% = 1, 0.5% = 2, max = 3
- Poisson fit: $E[K^2]/(E[K] + E[K]^2) = 1.005$ (perfect)
- $H_2$ per block: negligible compared to $N^2/\log^2 N$

**Theory needed**: BCZ/Corvaja-Zannier subspace theorem for P-recursive sequences.
No existing theorem covers this case (only constant-coefficient linear recurrences / S-units).

### Unbounded-order tower (Life window)
CRT tower: $A \bmod R^m$ determined at $m = 7$ (claimed, unverified).
For P3.2: need $m \to \infty$ with small-lift/non-stabilization theorem.
Threshold: $\rho \leq dc/(mn) + H_n/(mn)$; need $m/d \to \infty$.

---

## 16. EMPIRICAL DATA SUMMARY

### Close pairs for $K_n \geq 2$ (n ≤ 2000)
| n | K | Primes | Gaps | Same-kernel |
|---|---|--------|------|-------------|
| 200 | 2 | 139, 181 | 42 | Yes (g=1) |
| 272 | 2 | 191, 233 | 42 | Yes (g=1) |
| 300 | 2 | 191, 227 | 36 | Yes (g=1) |
| 321 | 3 | 179, 193, 211 | 14,32,18 | No |
| 755 | 2 | 593, 733 | 140 | Yes (g=2) |
| 1041 | 2 | 577, 769 | 192 | Yes (g=16) |
| 1103 | 2 | 761, 821 | 60 | No |
| 1265 | 2 | 769, 1069 | 300 | No |
| 1300 | 2 | 701, 937 | 236 | No |
| 1773 | 2 | 1229, 1427 | 198 | No |

5/12 pairs are same-kernel (42%). Gaps range 14-300 (mean 109).

### Answer files (2026-07-29)
| Q# | Source | Topic | Key finding |
|----|--------|-------|-------------|
| Q5546 | dm4 | Lucas collapse + methods fail | T_n = {p \| b_n} |
| Q5549 | dm1 | Collision energy | E₁ ≪ N^{8/3}, Path D |
| Q5550 | dm3 | Higher independence | H₂ suffices |
| Q5551 | dm5 | Continuant spacing | Aggregate only |
| Q5554 | dm2 | Global relation | Height approach dead |
| Q5555 | dm4 | Hasse square | No |Z_p| improvement |
| Q5538 | dm2 | Level 6 vs 8 forms | Midpoint only |
| Q5540 | dm2 | Modular form details | η(2τ)⁴η(4τ)⁴ |
