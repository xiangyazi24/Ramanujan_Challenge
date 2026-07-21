# Fable oracle report (2026-07-20) — horizontal dispersion route for (mh2)

## Lemma 0 (Lucas collapse) — VERIFIED, restructures everything
For m∈(N,2N], N<p≤m: base-p digits of m are (1,m−p), Gessel-Lucas gives
b_m ≡ 5·b_{m−p} mod p. So X_p(m)=1_{p|b_{m−p}} = 1_{p|b_m}, and
L_N(m) = #{p∈(N,m]: p|b_m}. Hence (mh2) is NOT cross-prime — it follows from
pointwise anti-concentration (♦): max_m #{p∈(N,m]: p|b_m} ≪ N^{o(1)}.
Size calibration: log b_2N≈7N vs Σlog p≈N ⟹ adversary size-permitted ⟹ needs
non-archimedean input on prime factorization of individual b_m.

## Q1 — character→geometry: ℓ-adic NO, crystalline YES
j∈Z_p is 𝔭|S_j, INVISIBLE to ℓ-adic topology. Gross-Koblitz/Stickelberger:
v_𝔭(g(ω^{−j}))=s_p(j)/(p−1) — mod-𝔭 vanishing governed by DIGIT combinatorics.
Character aspect mod p ≡ Stickelberger aspect, NOT monodromy. FKM/AFT in j = DEAD.
Positive: S_j = alternating sum of ≤B Frobenius eigenvalues; mod 𝔭 only slope-0
survive; exactly-one-unit ⟹ j∉Z_p. So Z_p ⊆ {j: unit-count ≠1} = E_p(digit-defined,
Hodge non-ordinary) ∪ A_p(≥2 units collide mod 𝔭, rate 1/p). Matches Poisson(1/2);
explains j=(p−1)/2 landmark (weight-4 level-8 form ordinarity).

## Q2 — two-characteristic coupling: DEAD + dissolved by Lemma 0. Zero effort.

## Q3 — bypass sheaf: only through Lemma 0; classical bypasses all die
Fourier detection Σe_p(ar+s·b_r) = Heilbronn-type, outside all machinery. DEAD.
Dispersion/vdC killed by size-permitted adversary. DEAD. Hooley-Δ DEAD. Lucas-
amplification propagates consistently, no contradiction. DEAD. Survivor: after
crystalline dichotomy, moving-target half = "digit-defined E_p uncorrelated with
affine j=m−p" (Mauduit-Rivat) + Poisson collision count on A_p.

## Q4 — LEMMA A (spec for Codex): Jacobi-sum skeleton + Stickelberger valuations
(A1) unfold c_{p,j}=−Σω^j(Λ(u)) to b_j ≡ −Σ_k c(j,k)Π J_i(j,k) mod 𝔭 (≤5 Jacobi/
Gauss sums, args affine in (j,k)); verify exact p≤200 all j; anchor j=(p−1)/2 ↔ γ_p.
(A2) Gross-Koblitz ⟹ v(k;j)=Σ⟨(αj+βk)/(p−1)⟩, U(j)=#{k:v=0}; U(j)=1⟹p∤b_j, so
Z_p⊆{j:U(j)≠1} — first structural container beyond size.
(A3) decision gate p≤2000-5000: measure |{U≠1}|. S1 polylog(alive)/S3 p^θ digit-
explicit(alive, Mauduit-Rivat)/S2 unstructured collisions(crystalline alone
insufficient — report, itself publishable). Secondary: any |Z_p|≪p^θ,θ<1 = new
Apéry theorem + down payment on (♦).

Dispatched to Codex gpt-5.6-sol xhigh as CODEX_SPEC_lemmaA.md.
