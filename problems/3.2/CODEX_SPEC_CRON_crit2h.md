# CODEX SPEC: [CRIT-2H] critical-value audit for f_h = N_h/A_h (implements Q6563 §12)

## Mission

Execute the exact computation plan of `chatgpt-answers/Q6563.md` §12 (READ THAT FILE FIRST, especially §1, §5, §12 and the gates). This is the machine side of the new [CRIT-2H] route to absolute irreducibility of the collision curves G_h — it replaces bivariate factorization with a one-variable critical-value audit.

## Objects

```
P(X)=34X^3+51X^2+27X+5
N_1(X)=1, N_2(X)=P(X+1), N_{h+1}(X)=P(X+h)N_h(X)-(X+h)^6 N_{h-1}(X)
A_h(X)=prod_{j=1..h}(X+j)^3
f_h = N_h/A_h
B_h = A_h N_h' - A_h' N_h
C_h = B_h / prod_{j=1..h}(X+j)^2     (exact division — this is Gate 2)
V_h(T) = primitive_part( Res_X( C_h(X), N_h(X) - T*A_h(X) ) )
```

## Gates (run first; report PASS/FAIL each; a FAIL is a finding, not a reason to stop the rest)

- Gate 0: verify the h=1 absolute factorization: G_1(X,Y) = -[(X+1)^2+(X+1)(Y+1)+(Y+1)^2] and that it splits over Q(omega).
- Gate 1: verify the pole-value identity N_h(-j) = (-1)^(j-1) * ((j-1)!)^3 * ((h-j)!)^3 * b[j-1] * b[h-j] for all 1<=j<=h, for h<=12 (b = Apery numbers: (m+1)^3 b_{m+1} = P(m) b_m - m^3 b_{m-1}, b_0=1, b_1=5).
- Gate 2: C_h is an exact polynomial quotient with deg(C_h)=4h-4 and lc(C_h)=-3*lc(N_h).
- Gate 3: mirror laws C_h(-h-1-X)=C_h(X) and V_h(-T)=V_h(T) up to scalar; V_h = const*W_h(T^2) with deg W_h = 2h-2.
- Gate 4: for h=2..6, cross-check [CRIT-2H] verdict against direct absolute factorization of G_h(X,Y) (over Q, plus a modular absolute-irreducibility check as in Q6563 §9.4).

## Main audit, for h = 2..H (push H as far as runtime allows; H>=30 target; switch to modular arithmetic per Q6563 §12 for large h)

```
g      = gcd(V_h, dV_h/dT)
rad    = V_h/g
repeat = gcd(rad, g)
simple = rad/repeat
if simple(0)==0: simple = simple/T
s_h    = deg(simple)
record: deg(V_h), deg(g), v_T(V_h at 0), s_h,
        full_Morse := (V_h(0)!=0 and deg(g)==0),
        crit_2h    := (s_h >= 2h-1)
```

For h beyond exact-arithmetic reach: compute everything mod one or several auxiliary primes ell > 3h (distinct ell as independent certificates); a modular s_{h,ell} >= 2h-1 with degree/content preservation is a characteristic-zero certificate for that h (Q6563 §9.4/§12) — CHECK and RECORD the degree-preservation conditions explicitly per (h, ell).

## Analysis outputs

1. Per-h table: s_h vs the threshold 2h-1 vs the Morse maximum 4h-4. Is full Morse holding? If s_h ever dips: report the colliding critical values exactly (which T-roots repeat, their multiplicity pattern) — a hidden symmetry there would be a major structural finding.
2. The sequence s_h: constant offset from 4h-4? Parity structure (mirror pairing a<->-a is expected — quantify beyond it)?
3. Store V_h coefficient data for h<=12 (JSON) for later p-adic/Newton-polygon analysis.

## Deliverables

- `CRON_crit2h.py` (progress prints per h, every <=10s during long h).
- `CODEX_CRIT2H_report.md`: gates, per-h table, analysis, verdict on [CRIT-2H] empirical status.
- `crit2h_results.json`.

## Rules

- Exact integer/rational arithmetic (sympy or flint) for h<=12; modular for larger h with ell recorded.
- No effort cap. If a gate fails, triple-check the conventions against Q6563 §1 (N_1=1 vs N_h(r)=K_{h-1}(r+1) indexing) before reporting — an indexing off-by-one is the most likely false alarm.
