# P2.7 Cooper level-11 trace basis decomposition
# 256^n * q_n = sum_{k=0}^n f_k * Phi11(n,k)
# Phi11(n,k) = 4^{n-k} C(2n,2k) T_{2k} - 2*4^{n-k-1} C(2n,2k+1) T_{2k+1}
# T_j = Cooper level-11 sequence.

import sys

NMAX = 400          # number of q_n terms
TMAX = 2*NMAX + 2   # need T up to 2n+1

# ---------- 1. Cooper T_j ----------
print("=== Step 1: Cooper level-11 T_j ===")
T = [QQ(0)]*(TMAX+1)
T[0] = QQ(1); T[1] = QQ(4); T[2] = QQ(28)
for j in range(2, TMAX):
    # (j+1)^3 T_{j+1} = 2(2j+1)(5j^2+5j+2) T_j - 8j(7j^2+1) T_{j-1} + 22j(2j-1)(j-1) T_{j-2}
    rhs = 2*(2*j+1)*(5*j^2+5*j+2)*T[j] - 8*j*(7*j^2+1)*T[j-1] + 22*j*(2*j-1)*(j-1)*T[j-2]
    T[j+1] = rhs / (j+1)^3

print("T_0..T_9:", [T[j] for j in range(10)])
assert T[3] == 268, "T_3 check failed: got %s" % T[3]
print("T_3 = 268 verified.")
allint = all(t in ZZ for t in T)
print("All T_j integers up to j=%d: %s" % (TMAX, allint))
sys.stdout.flush()

# ---------- 2. q_n via P2.7 recurrence ----------
print("\n=== Step 2: q_n (exact QQ) ===")
def Apoly(n): return 1024*(2*n+5)^4*(2*n+7)^3*(2*n+9)^3*(946*n^2+6407*n+10860)
def Bpoly(n): return 128*(2*n+7)^3*(2*n+9)^3*(104060*n^6+1745370*n^5+12145238*n^4+44886481*n^3+92943995*n^2+102256019*n+46709052)
def Cpoly(n): return 16*(n+3)^4*(2*n+9)^3*(3784*n^5+57792*n^4+351019*n^3+1059230*n^2+1587211*n+944620)
def Dpoly(n): return (n+3)^4*(n+4)^6*(946*n^2+4515*n+5399)

def run_recurrence(u0, u1, u2, N):
    u = [QQ(u0), QQ(u1), QQ(u2)]
    for n in range(2, N):
        # u_{n+1} = (B_n/A_n) u_n - (C_{n-1}/A_{n-1}) u_{n-1} + (D_{n-2}/A_{n-2}) u_{n-2}
        val = QQ(Bpoly(n))/QQ(Apoly(n)) * u[n] \
            - QQ(Cpoly(n-1))/QQ(Apoly(n-1)) * u[n-1] \
            + QQ(Dpoly(n-2))/QQ(Apoly(n-2)) * u[n-2]
        u.append(val)
        if (n+1) % 50 == 0:
            print("  computed u_%d" % (n+1)); sys.stdout.flush()
    return u

q = run_recurrence(QQ(-215040420000),
                   QQ(-167282265043404)/905,
                   QQ(-964185327658080)/6071, NMAX)
print("q_0..q_2:", q[0], q[1], q[2])
sys.stdout.flush()

# ---------- 3. U_n = 256^n q_n ----------
U = [QQ(256)^n * q[n] for n in range(NMAX+1)]

# ---------- 4. basis inversion ----------
print("\n=== Step 4: invert triangular system for f_k ===")
def Phi11(n, k):
    val = QQ(4)^(n-k) * binomial(2*n, 2*k) * T[2*k]
    if k < n:
        val -= 2 * QQ(4)^(n-k-1) * binomial(2*n, 2*k+1) * T[2*k+1]
    return val

def decompose(Uvals, N, label):
    f = []
    for k in range(N+1):
        s = sum(f[j]*Phi11(k, j) for j in range(k))
        fk = (Uvals[k] - s) / T[2*k]     # Phi11(k,k) = T_{2k}
        f.append(fk)
        if (k+1) % 50 == 0:
            print("  [%s] computed f_%d" % (label, k)); sys.stdout.flush()
    return f

f = decompose(U, NMAX, "f")

# ---------- 5. inspect f_k ----------
print("\n=== Step 5: first 20 f_k ===")
for k in range(20):
    print("f_%d = %s" % (k, f[k]))
ints = [k for k in range(NMAX+1) if f[k] in ZZ]
nonints = [k for k in range(NMAX+1) if f[k] not in ZZ]
print("Integer f_k count: %d / %d" % (len(ints), NMAX+1))
if nonints:
    print("First non-integer indices:", nonints[:10])
    dens = [f[k].denominator() for k in nonints[:20]]
    print("Their denominators:", dens)
    maxden = max(f[k].denominator() for k in range(NMAX+1))
    print("Max denominator over all k:", maxden if maxden < 10^40 else "%.3e" % float(maxden))
# growth rate
RRF = RealField(200)
print("\nGrowth: log|f_k|/k for sample k:")
for k in [10, 50, 100, 200, 300, 400]:
    if k <= NMAX and f[k] != 0:
        lg = RRF(abs(f[k])).log(10)
        root = (RRF(abs(f[k])).log()/k).exp()
        print("  k=%d: |f_k| ~ 10^%.2f, (|f_k|)^(1/k) ~ %.4f" % (k, float(lg), float(root)))
sys.stdout.flush()

# ---------- 6. guess recurrence for f_k ----------
messy = nonints and max(f[k].denominator() for k in range(min(50,NMAX+1))) > 10^6
if not messy:
    print("\n=== Step 6: ore_algebra.guess on f_k ===")
    from ore_algebra import OreAlgebra, guess
    R.<kk> = QQ['k']
    A = OreAlgebra(R, 'Sk')
    # clear denominators if small, else pass rationals directly
    if nonints:
        L = lcm([f[k].denominator() for k in range(NMAX+1)])
        print("Common denominator lcm bits: %d" % L.nbits())
        if L.nbits() < 64:
            flist = [ZZ(f[k]*L) for k in range(NMAX+1)]
        else:
            flist = [f[k] for k in range(NMAX+1)]
    else:
        flist = [ZZ(f[k]) for k in range(NMAX+1)]
    rec = None
    for cut in [120, 250, len(flist)]:
        try:
            rec = guess(flist[:cut], A)
            print("GUESS SUCCESS with %d terms" % cut)
            break
        except Exception as e:
            print("guess with %d terms failed: %s" % (cut, e))
    if rec is not None:
        print("order  =", rec.order())
        d = max(c.degree() for c in rec.coefficients())
        print("degree =", d)
        print("Recurrence operator:")
        print(rec)
        # Poincare polynomial: leading behavior — substitute Sk -> x, take top-degree coeff in k
        Px = PolynomialRing(QQ, 'x'); x = Px.gen()
        coeffs = rec.coefficients(sparse=False)
        poincare = sum(QQ(coeffs[i][d]) * x^i for i in range(len(coeffs)))
        print("Poincare polynomial (coeff of k^%d):" % d, poincare)
        print("Poincare roots:", poincare.roots(CC))
    sys.stdout.flush()
else:
    print("\n=== Step 6 skipped: f_k denominators explode ===")
    print("Trying ALTERNATIVE normalization 64^n * q_n")
    U64 = [QQ(64)^n * q[n] for n in range(NMAX+1)]
    f64 = decompose(U64, min(60, NMAX), "f64")
    print("first 20 f64_k:")
    for k in range(20):
        print("f64_%d = %s" % (k, f64[k]))

# ---------- 7. same for p_n ----------
print("\n=== Step 7: p_n decomposition ===")
# Canonical P2.7 p-sequence initial values (from PDF, see comprehensive_2_7.py):
p = run_recurrence(QQ(-612218384750),
                   QQ(-9525021973931919)/QQ(18100),
                   QQ(-29561828382772029)/QQ(65380), NMAX)
target = zeta(2) + zeta(3)
rat = RRF(p[NMAX])/RRF(q[NMAX])
print("p_400/q_400 = %.15f, zeta(2)+zeta(3) = %.15f" % (float(rat), float(target)))
Up = [QQ(256)^n * p[n] for n in range(NMAX+1)]
g = decompose(Up, NMAX, "g")
print("first 10 g_k:")
for k in range(10):
    print("g_%d = %s" % (k, g[k]))
print("\ng_k/f_k for k=0..20:")
for k in range(21):
    if f[k] != 0:
        r = g[k]/f[k]
        rstr = str(r) if max(r.numerator().nbits(), r.denominator().nbits()) < 200 else "(big rational)"
        print("k=%d: g/f = %s ~ %.10f" % (k, rstr, float(RRF(r))))
    else:
        print("k=%d: f_k = 0, g_k = %s" % (k, g[k]))
print("\nDONE")
