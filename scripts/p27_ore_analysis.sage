#!/usr/bin/env sage
"""P2.7: Find the polynomial-coefficient recurrence operator via ore_algebra,
then analyze its Ore structure: hypergeometric solutions, right factors, LCLM/GCRD."""
import sys, time
from ore_algebra import OreAlgebra, guess

t0 = time.time()

# Recurrence coefficients
def A(n): return 1024*(2*n+5)**4*(2*n+7)**3*(2*n+9)**3*(946*n**2+6407*n+10860)
def B(n): return 128*(2*n+7)**3*(2*n+9)**3*(104060*n**6+1745370*n**5+12145238*n**4+44886481*n**3+92943995*n**2+102256019*n+46709052)
def C(n): return 16*(n+3)**4*(2*n+9)**3*(3784*n**5+57792*n**4+351019*n**3+1059230*n**2+1587211*n+944620)
def D(n): return (n+3)**4*(n+4)**6*(946*n**2+4515*n+5399)

# Compute q_n exactly using QQ arithmetic
NMAX = 300
q = [QQ(0)] * (NMAX+1)
q[0] = QQ(-215040420000)
q[1] = QQ(-167282265043404) / QQ(905)
q[2] = QQ(-964185327658080) / QQ(6071)

print("Computing q_n for n=0..%d..." % NMAX, flush=True)
for n in range(2, NMAX):
    An = QQ(A(n)); Bn = QQ(B(n))
    Cn1 = QQ(C(n-1)); An1 = QQ(A(n-1))
    Dn2 = QQ(D(n-2)); An2 = QQ(A(n-2))
    q[n+1] = (Bn/An)*q[n] - (Cn1/An1)*q[n-1] + (Dn2/An2)*q[n-2]
    if n % 50 == 0:
        print("  n=%d done (%.1fs)" % (n, time.time()-t0), flush=True)

print("q_n computed (%.1fs)" % (time.time()-t0), flush=True)

# Step 1: Use ore_algebra.guess to find minimal polynomial recurrence
R = QQ['n']; n = R.gen()
Ore = OreAlgebra(R, 'Sn')
print("\nGuessing minimal polynomial-coefficient recurrence...", flush=True)
try:
    L = guess(q[:NMAX+1], Ore)
    print("Found operator L:")
    print("  order = %d" % L.order())
    print("  degree = %d" % L.degree())

    # Print the operator coefficients
    for i in range(L.order()+1):
        c = L[i]
        print("  L[%d] has degree %d" % (i, c.degree()))

    # Verify
    test = sum(L[i](n=k)*q[k+i] for i in range(L.order()+1) for k in [5])[0] if False else None
    # manual verify
    ok = True
    for k in range(NMAX - L.order()):
        val = sum(L[i](n=k)*q[k+i] for i in range(L.order()+1))
        if val != 0:
            print("  VERIFICATION FAILED at n=%d: val=%s" % (k, val))
            ok = False
            break
    if ok:
        print("  Verified for n=0..%d ✓" % (NMAX - L.order()))

except Exception as e:
    print("  guess failed: %s" % e)
    L = None

if L is None:
    print("\nFallback: trying with fewer terms...", flush=True)
    for nterms in [200, 150, 100]:
        try:
            L = guess(q[:nterms], Ore)
            print("  Found with %d terms: order=%d, degree=%d" % (nterms, L.order(), L.degree()))
            break
        except:
            continue

if L is None:
    print("FAILED to find recurrence operator. Exiting.")
    sys.exit(1)

print("\n=== Operator analysis (%.1fs) ===" % (time.time()-t0), flush=True)

# Step 2: Poincaré polynomial
print("\nPoincaré polynomial (leading coefficients of L[i]):")
poincare_coeffs = []
max_deg = max(L[i].degree() for i in range(L.order()+1))
for i in range(L.order()+1):
    c = L[i]
    lc = c.leading_coefficient() if c.degree() == max_deg else 0
    poincare_coeffs.append(lc)
    print("  L[%d]: degree %d, leading coeff = %s" % (i, c.degree(), lc))

Rx = QQ['x']
ppoly = Rx(poincare_coeffs)
print("Poincaré polynomial: %s" % ppoly)
print("Poincaré roots: %s" % [r[0] for r in ppoly.roots(QQbar)])

# Step 3: Search for hypergeometric (Petkovšek) solutions
print("\n=== Hypergeometric solution search (%.1fs) ===" % (time.time()-t0), flush=True)

# Method: try to find right factors of order 1
# A right factor of order 1 means L = Q * L1 where L1 = Sn - r(n) for some r(n) in QQ(n)
# This is equivalent to a hypergeometric solution v_n with v_{n+1}/v_n = r(n)
try:
    factors = L.right_factors(1)
    if factors:
        print("Found %d right factor(s) of order 1:" % len(factors))
        for f in factors:
            print("  %s" % f)
            # Extract the rational function r(n) = coefficient of Sn^0, negated
            r_n = -f[0] / f[1]
            print("  => v_{n+1}/v_n = %s" % r_n)
    else:
        print("No right factors of order 1 found (no hypergeometric solutions).")
except Exception as e:
    print("right_factors failed: %s" % e)

# Step 4: Try to factor the operator
print("\n=== Operator factorization (%.1fs) ===" % (time.time()-t0), flush=True)
try:
    factored = L.factor()
    print("Factorization: %s" % factored)
    print("Number of factors: %d" % len(factored))
    for i, (fac, mult) in enumerate(factored):
        print("  Factor %d (mult %d): order %d, degree %d" % (i, mult, fac.order(), fac.degree()))
except Exception as e:
    print("Factorization failed: %s" % e)

# Step 5: Check if the operator is self-adjoint or has special structure
print("\n=== Adjoint analysis ===" , flush=True)
try:
    Ladj = L.adjoint()
    print("L adjoint: order %d, degree %d" % (Ladj.order(), Ladj.degree()))
    # Check if L is self-adjoint (up to scalar)
    if L.order() == Ladj.order():
        # Compare leading coefficients
        ratio = Ladj[L.order()] / L[L.order()]
        is_const = ratio.is_constant()
        print("L_adj / L leading ratio is constant: %s" % is_const)
except Exception as e:
    print("Adjoint analysis failed: %s" % e)

# Step 6: Compute p_n and check its relation
print("\n=== Computing p_n (%.1fs) ===" % (time.time()-t0), flush=True)
p = [QQ(0)] * (NMAX+1)
p[0] = QQ(-612218384750)
p[1] = QQ(-9525021973931919) / QQ(18100)
p[2] = QQ(-29561828382772029) / QQ(65380)

for nn in range(2, NMAX):
    An = QQ(A(nn)); Bn = QQ(B(nn))
    Cn1 = QQ(C(nn-1)); An1 = QQ(A(nn-1))
    Dn2 = QQ(D(nn-2)); An2 = QQ(A(nn-2))
    p[nn+1] = (Bn/An)*p[nn] - (Cn1/An1)*p[nn-1] + (Dn2/An2)*p[nn-2]

# Check that p_n also satisfies the same operator L
print("Checking p_n under operator L...", flush=True)
ok = True
for k in range(min(50, NMAX - L.order())):
    val = sum(L[i](n=k)*p[k+i] for i in range(L.order()+1))
    if val != 0:
        print("  p_n FAILS at n=%d" % k)
        ok = False
        break
if ok:
    print("  p_n satisfies L ✓ (checked n=0..49)")

# Step 7: Look for an operator annihilating the "error" e_n = p_n - L*q_n
# Actually, e_n is also annihilated by L since both p_n and q_n are.
# But maybe e_n satisfies a LOWER order operator?
# e_n = p_n - (zeta(2)+zeta(3)) * q_n — this requires a specific value of L.
# Instead, check: does the operator have a right GCRD with a shifted version?

print("\n=== Summary ===" , flush=True)
print("Operator L: order %d, degree %d" % (L.order(), L.degree()))
print("Total time: %.1fs" % (time.time()-t0))
