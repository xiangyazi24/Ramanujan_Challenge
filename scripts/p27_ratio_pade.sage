#!/usr/bin/env sage
"""P2.7: Search for rational function r_{n+1}/r_n where r_n = q_n/W_n."""

# Problem 2.7 recurrence coefficients
def A27(n):
    return 1024*(2*n+5)^4*(2*n+7)^3*(2*n+9)^3*(946*n^2+6407*n+10860)

def B27(n):
    P6 = 104060*n^6+1745370*n^5+12145238*n^4+44886481*n^3+92943995*n^2+102256019*n+46709052
    return 128*(2*n+7)^3*(2*n+9)^3*P6

def C27(n):
    P5 = 3784*n^5+57792*n^4+351019*n^3+1059230*n^2+1587211*n+944620
    return 16*(n+3)^4*(2*n+9)^3*P5

def D27(n):
    return (n+3)^4*(n+4)^6*(946*n^2+4515*n+5399)

# Cooper T_k
T = [QQ(1), QQ(4), QQ(28)]
for k in range(2, 250):
    t = (2*(2*k+1)*(5*k^2+5*k+2)*T[k] - 8*k*(7*k^2+1)*T[k-1] + 22*k*(2*k-1)*(k-1)*T[k-2]) / (k+1)^3
    T.append(t)

# W_n
NMAX = 60
W = []
for n in range(NMAX):
    val = QQ(0)
    for j in range(min(2*n+1, len(T))):
        val += binomial(2*n, j) * QQ(-2)^(2*n-j) * T[j]
    W.append(val / QQ(256)^n)

# q_n
q = [QQ(-215040420000), QQ(-167282265043404)/QQ(905), QQ(-964185327658080)/QQ(6071)]
for n in range(2, NMAX-1):
    q_next = (B27(n)*q[n] - C27(n)*q[n-1] + D27(n)*q[n-2]) / A27(n)
    q.append(q_next)

# r_n = q_n / W_n
r = [q[n]/W[n] for n in range(min(len(q), len(W)))]

# Consecutive ratios rr_n = r_{n+1}/r_n
rr = [r[n+1]/r[n] for n in range(len(r)-1)]

print("=== r_{n+1}/r_n ===", flush=True)
for n in range(min(20, len(rr))):
    print(f"  n={n}: {float(rr[n]):.15f}")

# Padé search: rr_n = P(n)/Q(n)
R_poly.<x> = PolynomialRing(QQ)

for d in range(1, 12):
    n_unk = 2*d + 1
    n_pts = min(n_unk + 5, len(rr) - 2)
    if n_pts < n_unk:
        print(f"  d={d}: not enough data"); break

    M = matrix(QQ, n_pts, n_unk)
    b = vector(QQ, n_pts)
    for i in range(n_pts):
        nn = QQ(i + 1)
        for k in range(d+1):
            M[i, k] = -nn^k
        for k in range(d):
            M[i, d+1+k] = rr[i+1] * nn^k
        b[i] = -rr[i+1] * nn^d

    try:
        sol = M.solve_right(b)
    except ValueError:
        print(f"  d={d}: no solution"); continue

    # Verify holdout
    ok = True
    for i in range(n_pts, min(n_pts+10, len(rr)-1)):
        nn = QQ(i+1)
        P_val = sum(sol[k]*nn^k for k in range(d+1))
        Q_val = nn^d + sum(sol[d+1+k]*nn^k for k in range(d))
        if Q_val == 0 or P_val/Q_val != rr[i+1]:
            ok = False; break

    if ok:
        P_poly = sum(sol[k]*x^k for k in range(d+1))
        Q_poly = x^d + sum(sol[d+1+k]*x^k for k in range(d))
        print(f"\n*** MATCH at degree {d}: r_{{n+1}}/r_n = P(n)/Q(n) ***")
        print(f"  P(n) = {P_poly}")
        print(f"  Q(n) = {Q_poly}")
        print(f"  P factors: {P_poly.factor()}")
        print(f"  Q factors: {Q_poly.factor()}")

        # Full verification
        all_ok = True
        for i in range(1, len(rr)):
            nn = QQ(i)
            P_val = P_poly(nn)
            Q_val = Q_poly(nn)
            if Q_val == 0 or P_val/Q_val != rr[i]:
                all_ok = False
                print(f"  FAIL at n={i}")
                break
        if all_ok:
            print(f"  VERIFIED for all n=1..{len(rr)-1}")
        break
    else:
        if d <= 3:
            print(f"  d={d}: holdout fails")

# Also try: is r_n itself a polynomial/rational function of n?
print("\n=== Is r_n = P(n)/Q(n)? ===")
for d in range(1, 8):
    n_unk = 2*d + 1
    n_pts = min(n_unk + 5, len(r) - 2)
    if n_pts < n_unk: break

    M = matrix(QQ, n_pts, n_unk)
    b = vector(QQ, n_pts)
    for i in range(n_pts):
        nn = QQ(i)
        for k in range(d+1):
            M[i, k] = -nn^k
        for k in range(d):
            M[i, d+1+k] = r[i] * nn^k
        b[i] = -r[i] * nn^d

    try:
        sol = M.solve_right(b)
    except ValueError:
        continue

    ok = True
    for i in range(n_pts, min(n_pts+5, len(r))):
        nn = QQ(i)
        P_val = sum(sol[k]*nn^k for k in range(d+1))
        Q_val = nn^d + sum(sol[d+1+k]*nn^k for k in range(d))
        if Q_val == 0 or P_val/Q_val != r[i]:
            ok = False; break
    if ok:
        P_poly = sum(sol[k]*x^k for k in range(d+1))
        Q_poly = x^d + sum(sol[d+1+k]*x^k for k in range(d))
        print(f"\n*** r_n = P(n)/Q(n) at degree {d} ***")
        print(f"  P(n) = {P_poly}")
        print(f"  Q(n) = {Q_poly}")
        break

# Check: is q_n / (h_n * W_n) a polynomial for specific h_n?
# Try h_n = product of (a_i + j) for j=0..n-1
print("\n=== Checking q_n / (W_n * Gamma-ratio) ===")
# Try h_n such that h_{n+1}/h_n matches certain patterns
# From the exponent gap analysis: need h_{n+1}/h_n ~ 1 + 3/(2n) + O(1/n^2)
# Central binomial: C(2n,n)/4^n has ratio (2n+1)/(2n+2) ~ 1 - 1/(2n)
# (3/2)_n/n! has ratio (n+1/2)/(n+1) ~ 1 - 1/(2n)
# Need opposite sign: (n+a)/(n+b) ~ 1 + (a-b)/n

# For 3/(2n): need (a-b) = 3/2 in ratio (n+a)/(n+b)
# E.g., h_n = (n+5/2)!/(n+1)! = (5/2)_n * ...
# Let's just compute h_n = prod_{j=0}^{n-1} (2j+5)/(2j+2) = (5/2)_n * 2^n / (2)_n
print("  h_n = prod (2j+5)/(2j+2) = (5/2)_n * 2^n / (2)_n")
h_vals = [QQ(1)]
for n in range(1, NMAX):
    h_vals.append(h_vals[-1] * QQ(2*(n-1)+5) / QQ(2*(n-1)+2))

rh = [q[n]/(W[n]*h_vals[n]) for n in range(min(len(q),len(W),len(h_vals)))]
print("  q/(W*h) ratios:")
for n in range(min(10, len(rh))):
    print(f"    n={n}: {float(rh[n]):.10f}")

# Check finite diffs
print("  Diffs:")
for n in range(1, min(8, len(rh))):
    print(f"    Δ(n={n}) = {float(rh[n]-rh[n-1]):.10f}")

print("\nDone.")
