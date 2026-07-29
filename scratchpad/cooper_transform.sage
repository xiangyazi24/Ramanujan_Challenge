#!/usr/bin/env sage
"""
Test: does the centered binomial transform of Cooper level-11 solutions
produce P2.7 solutions?

Cooper: (n+1)^3 T_{n+1} = 2(2n+1)(5n^2+5n+2) T_n - 8n(7n^2+1) T_{n-1}
        + 22n(2n-1)(n-1) T_{n-2}
P2.7:   A(n) q_{n+1} = B(n) q_n - C(n) q_{n-1} + D(n) q_{n-2}

Centered binomial transform: (KT)_n = (1/256^n) Σ_{k=0}^{2n} C(2n,k)(-2)^{2n-k} T_k

Spectral identity: 16 P27((t-2)^2/4) = -C11(t) C11(4-t)
"""
import mpmath
mpmath.mp.dps = 100

# ============================================================
# Cooper level-11 sequence
# ============================================================
N_COOP = 100
T = [QQ(0)] * N_COOP
T[0] = QQ(1)
T[1] = QQ(4)
T[2] = QQ(28)

for n in range(2, N_COOP - 1):
    T[n+1] = (2*(2*n+1)*(5*n^2+5*n+2)*T[n] - 8*n*(7*n^2+1)*T[n-1]
              + 22*n*(2*n-1)*(n-1)*T[n-2]) / (n+1)^3

print("Cooper T_k first values:")
for k in range(10):
    print("  T_%d = %s" % (k, T[k]))

# Verify Cooper recurrence
print("\nVerification of Cooper recurrence:")
for n in range(2, 8):
    lhs = (n+1)^3 * T[n+1]
    rhs = 2*(2*n+1)*(5*n^2+5*n+2)*T[n] - 8*n*(7*n^2+1)*T[n-1] + 22*n*(2*n-1)*(n-1)*T[n-2]
    print("  n=%d: match=%s" % (n, lhs == rhs))

# ============================================================
# Centered even binomial transform
# ============================================================
def centered_transform(T_seq, n, normalize=True):
    """(KT)_n = (1/256^n) Σ_{k=0}^{2n} C(2n,k) (-2)^{2n-k} T_k"""
    val = QQ(0)
    for k in range(2*n + 1):
        val += binomial(2*n, k) * (-2)^(2*n - k) * T_seq[k]
    if normalize:
        val /= QQ(256)^n
    return val

print("\n=== Centered binomial transform (KT)_n ===")
KT = []
for n in range(20):
    val = centered_transform(T, n)
    KT.append(val)
    print("  (KT)_%d = %s" % (n, float(val)))

# ============================================================
# Check if (KT)_n satisfies P2.7
# ============================================================
def A_p27(n):
    return QQ(1024) * (2*n+5)^4 * (2*n+7)^3 * (2*n+9)^3 * (946*n^2+6407*n+10860)
def B_p27(n):
    return QQ(128) * (2*n+7)^3 * (2*n+9)^3 * (104060*n^6 + 1745370*n^5 +
        12145238*n^4 + 44886481*n^3 + 92943995*n^2 + 102256019*n + 46709052)
def C_p27(n):
    return QQ(16) * (n+3)^4 * (2*n+9)^3 * (3784*n^5 + 57792*n^4 +
        351019*n^3 + 1059230*n^2 + 1587211*n + 944620)
def D_p27(n):
    return QQ(1) * (n+3)^4 * (n+4)^6 * (946*n^2 + 4515*n + 5399)

print("\n=== P2.7 recurrence residual for (KT) ===")
for n in range(2, min(15, len(KT)-3)):
    res = A_p27(n) * KT[n+1] - B_p27(n) * KT[n] + C_p27(n-1)/A_p27(n-1) * KT[n-1]
    # Wait, the recurrence is: u_{n+1} = B/A u_n - C_{n-1}/A_{n-1} u_{n-1} + D_{n-2}/A_{n-2} u_{n-2}
    # Written as: A_n u_{n+1} - B_n u_n + C_{n-1} u_{n-1}/... no.
    #
    # The MONIC form: u_{n+1} = (B/A) u_n - (C_{n-1}/A_{n-1}) u_{n-1} + (D_{n-2}/A_{n-2}) u_{n-2}
    # The STANDARD form: P3(n) u_{n+3} + P2(n) u_{n+2} + P1(n) u_{n+1} + P0(n) u_n = 0
    # From the summary: A_n u_{n+1} = B_n u_n - C_{n-1}/A_{n-1} u_{n-1}...
    # Actually the recurrence is:
    # u_{n+1} = (B(n)/A(n)) u_n - (C(n-1)/A(n-1)) u_{n-1} + (D(n-2)/A(n-2)) u_{n-2}
    #
    # So the residual is:
    # u_{n+1} - (B(n)/A(n)) u_n + (C(n-1)/A(n-1)) u_{n-1} - (D(n-2)/A(n-2)) u_{n-2}
    pass

# Let me redo with the correct form
print("=== P2.7 residual: u_{n+1} - (B/A)u_n + (C_{n-1}/A_{n-1})u_{n-1} - (D_{n-2}/A_{n-2})u_{n-2} ===")
for n in range(2, min(15, len(KT)-1)):
    pred = B_p27(n)/A_p27(n) * KT[n] - C_p27(n-1)/A_p27(n-1) * KT[n-1] + D_p27(n-2)/A_p27(n-2) * KT[n-2]
    actual = KT[n+1]
    res = actual - pred
    if actual != 0:
        rel = float(res / actual)
    else:
        rel = float(res)
    print("  n=%d: residual = %e (relative: %e)" % (n, float(res), rel))

# ============================================================
# Also try the UN-normalized transform: ΣC(2n,k)(-2)^{2n-k}T_k (no 256^n division)
# ============================================================
print("\n=== Un-normalized transform ===")
KT_unnorm = []
for n in range(20):
    val = centered_transform(T, n, normalize=False)
    KT_unnorm.append(val)

print("=== P2.7 residual for un-normalized ===")
for n in range(2, min(15, len(KT_unnorm)-1)):
    pred = B_p27(n)/A_p27(n) * KT_unnorm[n] - C_p27(n-1)/A_p27(n-1) * KT_unnorm[n-1] + D_p27(n-2)/A_p27(n-2) * KT_unnorm[n-2]
    actual = KT_unnorm[n+1]
    res = actual - pred
    if actual != 0:
        rel = float(res / actual)
    else:
        rel = float(res)
    print("  n=%d: rel residual = %e" % (n, rel))

# ============================================================
# Try with Pochhammer gauge: v_n = 256^n * h_n * (KT)_n
# ============================================================
def gamma_gauge(n):
    if n == 0:
        return QQ(1)
    num = QQ(2)^(-20*n)
    for j in range(n):
        num *= (3+j)^4 * (4+j)^6
    den = QQ(1)
    for j in range(n):
        den *= (QQ(5)/2+j)^4 * (QQ(7)/2+j)^3 * (QQ(9)/2+j)^3
    return num / den

print("\n=== Gauged transform: g_n * (KT)_n ===")
# Try different gauge factors
for gauge_name, gauge_fn in [
    ("γ_n", gamma_gauge),
    ("256^n·γ_n", lambda n: QQ(256)^n * gamma_gauge(n)),
    ("64^n·γ_n", lambda n: QQ(64)^n * gamma_gauge(n)),
]:
    print("\nGauge: %s" % gauge_name)
    gauged = [gauge_fn(n) * KT[n] for n in range(min(15, len(KT)))]
    for n in range(2, min(12, len(gauged)-1)):
        pred = B_p27(n)/A_p27(n) * gauged[n] - C_p27(n-1)/A_p27(n-1) * gauged[n-1] + D_p27(n-2)/A_p27(n-2) * gauged[n-2]
        actual = gauged[n+1]
        res = actual - pred
        if actual != 0:
            rel = float(res / actual)
        else:
            rel = float(res)
        if abs(rel) < 1e-10:
            print("  n=%d: *** MATCH! *** rel = %e" % (n, rel))
        else:
            print("  n=%d: rel = %e" % (n, rel))

# ============================================================
# Check what recurrence (KT)_n DOES satisfy
# ============================================================
print("\n=== Guess recurrence for (KT)_n ===")
from ore_algebra import OreAlgebra, guess
Rn = PolynomialRing(QQ, 'nn')
OS = OreAlgebra(Rn, 'Snn')

try:
    rec = guess(KT[:18], OS)
    print("Found recurrence! Order %d" % rec.order())
    for j in range(rec.order()+1):
        print("  P_%d: degree %d" % (j, rec[j].degree()))
    # Check Poincaré polynomial
    P_leading = [rec[j].leading_coefficient() for j in range(rec.order()+1)]
    print("Poincaré coefficients: %s" % P_leading)
    poin = sum(P_leading[j] * var('mu')^j for j in range(len(P_leading)))
    print("Poincaré polynomial: %s" % poin)
    print("P27: 4μ³ - 220μ² + 8μ - 1")
except Exception as e:
    print("No relation found: %s" % e)

# Also try with explicit order
for order in [3, 4, 5, 6]:
    try:
        rec = guess(KT[:18], OS, order=order)
        print("Order %d: found!" % order)
        for j in range(rec.order()+1):
            print("  P_%d: degree %d" % (j, rec[j].degree()))
        break
    except:
        print("Order %d: no relation" % order)

# ============================================================
# Compare (KT) ratios with P2.7 ratios
# ============================================================
print("\n=== Growth comparison ===")
print("(KT)_{n+1}/(KT)_n vs q_{n+1}/q_n:")

q_p27 = [QQ(0)] * 20
q_p27[0] = QQ(-215040420000)
q_p27[1] = QQ(-167282265043404) / QQ(905)
q_p27[2] = QQ(-964185327658080) / QQ(6071)
for n in range(2, 19):
    q_p27[n+1] = B_p27(n)/A_p27(n)*q_p27[n] - C_p27(n-1)/A_p27(n-1)*q_p27[n-1] + D_p27(n-2)/A_p27(n-2)*q_p27[n-2]

for n in range(min(12, len(KT)-1)):
    if KT[n] != 0:
        r_kt = float(KT[n+1]/KT[n])
        r_p27 = float(q_p27[n+1]/q_p27[n])
        print("  n=%d: KT ratio = %.10f, P2.7 ratio = %.10f" % (n, r_kt, r_p27))
