# p27_basis_decomp.sage
# Problem 2.7 (Ramanujan Challenge): q_n -> zeta(2)+zeta(3) via 4-term recurrence.
# Replicate the P2.5 strategy: decompose q_n = sum_k f(k) * Phi(n,k) over several
# triangular bases, invert triangularly for f(k), and try ore_algebra.guess on f.
#
# Run (from the repository root): sage scripts/p27_basis_decomp.sage

import time
from ore_algebra import OreAlgebra, guess

# ----------------------------------------------------------------------------
# Recurrence coefficients (degree 12 in n)
# ----------------------------------------------------------------------------
def Acoef(n): return 1024*(2*n+5)**4*(2*n+7)**3*(2*n+9)**3*(946*n**2+6407*n+10860)
def Bcoef(n): return 128*(2*n+7)**3*(2*n+9)**3*(104060*n**6+1745370*n**5+12145238*n**4+44886481*n**3+92943995*n**2+102256019*n+46709052)
def Ccoef(n): return 16*(n+3)**4*(2*n+9)**3*(3784*n**5+57792*n**4+351019*n**3+1059230*n**2+1587211*n+944620)
def Dcoef(n): return (n+3)**4*(n+4)**6*(946*n**2+4515*n+5399)

NMAX = 400   # number of terms (indices 0..NMAX)

# ----------------------------------------------------------------------------
# Step 1: compute q_n exactly over QQ
# ----------------------------------------------------------------------------
print("="*78)
print("Step 1: computing q_n exactly (QQ) for n = 0..%d" % NMAX)
print("="*78)
t0 = time.time()

q = [QQ(0)]*(NMAX+1)
q[0] = QQ(-215040420000)
q[1] = QQ(-167282265043404)/905
q[2] = QQ(-964185327658080)/6071

for n in range(2, NMAX):
    # u_{n+1} = (B_n/A_n) u_n - (C_{n-1}/A_{n-1}) u_{n-1} + (D_{n-2}/A_{n-2}) u_{n-2}
    q[n+1] = (QQ(Bcoef(n))/QQ(Acoef(n)))*q[n] \
           - (QQ(Ccoef(n-1))/QQ(Acoef(n-1)))*q[n-1] \
           + (QQ(Dcoef(n-2))/QQ(Acoef(n-2)))*q[n-2]
    if (n+1) % 50 == 0:
        print("  computed q_%d  (%.1fs; numerator digits ~%d, denominator digits ~%d)"
              % (n+1, time.time()-t0,
                 len(str(q[n+1].numerator())), len(str(q[n+1].denominator()))))

print("q computation done in %.1fs" % (time.time()-t0))

# sanity: numeric limit check against zeta(2)+zeta(3)
target = zeta(2) + zeta(3)
print("\nSanity check of limit:")
print("  q_%d (numeric)      = %s" % (NMAX, numerical_approx(q[NMAX], digits=30)))
print("  zeta(2)+zeta(3)     = %s" % numerical_approx(target, digits=30))
print("  |q_N - target|      = %s" % numerical_approx(abs(q[NMAX] - target), digits=10))

# ----------------------------------------------------------------------------
# Step 4 (early): look at the arithmetic of q_n / natural normalization
# ----------------------------------------------------------------------------
print("\n" + "="*78)
print("Step 4: arithmetic of q_n (normalization probe)")
print("="*78)
print("q_0 = %s = %s" % (q[0], factor(q[0])))
print("\nDenominators of q_n, n = 0..12 (factored):")
for i in range(13):
    d = q[i].denominator()
    print("  den(q_%-2d) = %-20s = %s" % (i, d, factor(d)))
L20 = lcm([q[i].denominator() for i in range(21)])
print("\nlcm(den(q_0..q_20)) = %s" % factor(L20))
print("\ngcd of numerators of q_0..q_12: %s" % factor(gcd([q[i].numerator() for i in range(13)])))

# ----------------------------------------------------------------------------
# Step 2/3: basis decompositions + guessing
# ----------------------------------------------------------------------------
BASES = [
    ("(a) Binomial            C(n,k)",
        lambda nn, kk: binomial(nn, kk)),
    ("(b) Apery-2   C(n,k)^2 C(n+k,k)",
        lambda nn, kk: binomial(nn, kk)**2 * binomial(nn+kk, kk)),
    ("(c) Apery-3   C(n,k)^2 C(n+k,k)^2",
        lambda nn, kk: binomial(nn, kk)**2 * binomial(nn+kk, kk)**2),
    ("(d) Delannoy  2^k C(2k,k) C(n,k) C(n+k,k)",
        lambda nn, kk: 2**kk * binomial(2*kk, kk) * binomial(nn, kk) * binomial(nn+kk, kk)),
    ("(e) CDelannoy^2  4^k C(2k,k)^2 C(n,k) C(n+k,k)",
        lambda nn, kk: 4**kk * binomial(2*kk, kk)**2 * binomial(nn, kk) * binomial(nn+kk, kk)),
]

Rk = PolynomialRing(QQ, 'k')
Aop = OreAlgebra(Rk, 'Sk')
lam = polygen(QQ, 'lam')

def poincare_poly(rec):
    """Characteristic polynomial at infinity: sum of top-degree coefficients."""
    coeffs = rec.coefficients(sparse=False)   # p_0(k), ..., p_r(k)
    polys = [Rk(c) for c in coeffs]
    d = max(p.degree() for p in polys if p != 0)
    return sum(QQ(p[d]) * lam**i for i, p in enumerate(polys))

def describe_niceness(f, label):
    head = f[:10]
    print("  first 10 f(k): ")
    for i, v in enumerate(head):
        print("    f(%d) = %s" % (i, v))
    n_int = sum(1 for v in f if v.denominator() == 1)
    print("  integrality: %d / %d values are integers" % (n_int, len(f)))
    if n_int < len(f):
        dens = [f[i].denominator() for i in range(min(30, len(f)))]
        print("  denominators of f(0..29): %s" % dens)
        big = f[min(30, len(f)-1)].denominator()
        print("  den(f(30)) factored: %s" % factor(big))
    sizes = [len(str(f[i].numerator())) for i in (10, 50, 100, min(200, len(f)-1))]
    print("  numerator digit sizes at k=10,50,100,~200: %s" % sizes)

def try_guess(f, label):
    """Try ore_algebra.guess with 400 terms, then 200, then 120."""
    for cut in (len(f), 200, 120):
        data = [QQ(v) for v in f[:cut]]
        try:
            t = time.time()
            rec = guess(data, Aop)
            dt = time.time()-t
            coeffs = rec.coefficients(sparse=False)
            degs = [Rk(c).degree() for c in coeffs]
            print("  GUESS SUCCESS with %d terms (%.1fs): order %d, max coeff degree %d"
                  % (cut, dt, rec.order(), max(degs)))
            print("  coefficient degrees: %s" % degs)
            pp = poincare_poly(rec)
            print("  Poincare polynomial (char poly at infinity):")
            print("    %s" % pp)
            print("    factored: %s" % factor(pp))
            try:
                rts = pp.roots(QQbar, multiplicities=False)
                print("    roots (QQbar->CC): %s" % [numerical_approx(CC(r), digits=8) for r in rts])
            except Exception as e:
                print("    (root computation failed: %s)" % e)
            if rec.order() <= 6 and max(degs) <= 20:
                print("  operator: %s" % rec)
            return rec, cut
        except Exception as e:
            print("  guess with %d terms failed: %s" % (cut, str(e)[:120]))
    return None, None

summary = []
for label, Phi in BASES:
    print("\n" + "="*78)
    print("Basis %s" % label)
    print("="*78)
    t0 = time.time()
    # triangular inversion: f(k) = (q_k - sum_{j<k} f(j) Phi(k,j)) / Phi(k,k)
    kmax = NMAX
    f = []
    ok = True
    for k in range(kmax+1):
        s = q[k]
        for j in range(k):
            s -= f[j] * Phi(k, j)
        pkk = Phi(k, k)
        if pkk == 0:
            print("  Phi(%d,%d) = 0 -- basis not triangular-invertible, skipping" % (k, k))
            ok = False
            break
        f.append(s / pkk)
        if k > 0 and k % 100 == 0:
            print("  inverted through k=%d (%.1fs)" % (k, time.time()-t0))
        # fallback: if inversion is dragging, stop at 200
        if k == 200 and time.time()-t0 > 240:
            print("  inversion slow; falling back to 200 terms")
            kmax = 200
            break
    if not ok:
        summary.append((label, "not invertible", None, None, None))
        continue
    print("  inversion done: %d terms in %.1fs" % (len(f), time.time()-t0))

    describe_niceness(f, label)
    rec, used = try_guess(f, label)
    if rec is not None:
        coeffs = rec.coefficients(sparse=False)
        degs = [Rk(c).degree() for c in coeffs]
        summary.append((label, "recurrence found (%d terms)" % used,
                        rec.order(), max(degs), poincare_poly(rec)))
    else:
        summary.append((label, "no recurrence found", None, None, None))

# ----------------------------------------------------------------------------
# Step 5: summary
# ----------------------------------------------------------------------------
print("\n" + "="*78)
print("SUMMARY")
print("="*78)
for label, status, order, deg, pp in summary:
    print("%s" % label)
    print("    status: %s" % status)
    if order is not None:
        print("    order = %d, max coeff degree = %d" % (order, deg))
        print("    Poincare poly: %s  = %s" % (pp, factor(pp)))
print("\nDone.")
