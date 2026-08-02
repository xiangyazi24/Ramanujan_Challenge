#!/usr/bin/env sage
"""
P2.5: full k-recurrence analysis for the Delannoy-basis coefficients.

f_e1(k): coefficients of Qhat^{e1}_N = (prod M_H)_{11}  (f(0)=1, f(1)=5749/3136)
f_q(k):  coefficients of the combined q-row Qhat_N, q = (33750,-36000,9000)
g_p(k):  coefficients of the combined p-row Phat_N, p = (30921,-32972,8240)
         (g_p/f_q -> Catalan G geometrically)

For each: guess minimal recurrence with ore_algebra, verify EXACTLY on all
computed values, do Poincare analysis, and compare the three operators.  The
script also computes LCLM(L_f,L_g) and writes all primitive integral
coefficients to lean/specs/p25_k_recurrence_coefficients.txt.
"""
import sys, time
from pathlib import Path
from ore_algebra import OreAlgebra, guess

KMAX = 400

def M_entries(n):
    m11 = (-2*n-5)*(n+3)**2*(136*n**4+1424*n**3+5548*n**2+9551*n+6141)
    m12 = 384*n**6+6384*n**5+44168*n**4+162698*n**3+336377*n**2+369933*n+169011
    m13 = -(480*n**4+4980*n**3+19210*n**2+32690*n+20730)
    m21 = (n+2)**2*(n+3)**2*(4*n+10)*(48*n**3+386*n**2+1017*n+879)
    m22 = (n+2)**2*(-272*n**5-3848*n**4-21732*n**3-61184*n**2-85761*n-47808)
    m23 = (n+2)**2*(320*n**3+2540*n**2+6610*n+5640)
    m31 = (-4*n-10)*(n+2)**2*(n+3)**2*(32*n**4+302*n**3+1037*n**2+1530*n+813)
    m32 = (n+2)**2*(192*n**6+2984*n**5+19116*n**4+64452*n**3+120256*n**2+117279*n+46476)
    m33 = (n+2)**2*(-16*n**5-408*n**4-2912*n**3-8884*n**2-12254*n-6240)
    return [[m11,m12,m13],[m21,m22,m23],[m31,m32,m33]]

def delta_H(n):
    return -2*(n+2)**2*(n+3)**2*(2*n+5)*(2*n+7)**2

t0 = time.time()
E = [vector(QQ, [1,0,0]), vector(QQ, [0,1,0]), vector(QQ, [0,0,1])]
first_col = [[QQ(E[j][0]) for j in range(3)]]  # (prod M)_{j,1} history
cur = [vector(v) for v in E]
for N in range(KMAX):
    M = matrix(QQ, M_entries(N)) / delta_H(N)
    cur = [v * M for v in cur]
    first_col.append([cur[j][0] for j in range(3)])
q = vector(QQ, [33750, -36000, 9000])
p = vector(QQ, [30921, -32972, 8240])
Qe1_vals = [first_col[N][0] for N in range(KMAX+1)]
Q_vals = [sum(q[j]*first_col[N][j] for j in range(3)) for N in range(KMAX+1)]
P_vals = [sum(p[j]*first_col[N][j] for j in range(3)) for N in range(KMAX+1)]
print("trajectories done (%.1fs)" % (time.time()-t0)); sys.stdout.flush()

def Bsummand(N, k):
    return 2**k * binomial(2*k, k) * binomial(N, k) * binomial(N+k, k)

def decompose(vals):
    c = []
    for K in range(len(vals)):
        rhs = vals[K]
        for k in range(K):
            rhs -= c[k] * Bsummand(K, k)
        c.append(rhs / Bsummand(K, K))
    return c

f_e1 = decompose(Qe1_vals)
f_q  = decompose(Q_vals)
g_p  = decompose(P_vals)
assert f_e1[1] == 5749/3136 and f_e1[2] == 16811771/4572288
print("decompositions done (%.1fs)" % (time.time()-t0)); sys.stdout.flush()

R = QQ['k']; kk = R.gen()
A = OreAlgebra(R, 'Sk')

def apply_op(L, seq):
    """Return first bad k, or None if L annihilates seq everywhere."""
    r = L.order()
    cs = [L[i] for i in range(r+1)]
    for k in range(len(seq)-r):
        if sum(cs[i](k)*seq[k+i] for i in range(r+1)) != 0:
            return k
    return None

def analyze(seq, name):
    print("="*70)
    print("[%s]" % name)
    t1 = time.time()
    try:
        L = guess(seq, A)
    except ValueError:
        print("  no recurrence found with %d terms" % len(seq))
        return None
    r, d = L.order(), L.degree()
    print("  order = %d, degree = %d  (guessed in %.1fs)" % (r, d, time.time()-t1))
    bad = apply_op(L, seq)
    print("  EXACT verification on all k=0..%d: %s"
          % (len(seq)-1-r, "PASS" if bad is None else "FAIL at k=%d" % bad))
    lc = L.leading_coefficient()
    print("  leading coefficient c_%d(k) factored:" % r)
    print("   ", lc.factor())
    print("  trailing coefficient c_0(k) factored:")
    print("   ", L[0].factor())
    cs = [L[i] for i in range(r+1)]
    D = max(c.degree() for c in cs)
    xi = polygen(QQ, 'xi')
    poin = sum(c[D]*xi**i for i, c in enumerate(cs))
    poin /= poin.content() if hasattr(poin, 'content') else 1
    print("  Poincare polynomial:", poin.factor())
    print("  Poincare roots:", sorted(poin.roots(CC, multiplicities=True),
                                      key=lambda t: abs(t[0])))
    sys.stdout.flush()
    return L

Le1 = analyze(f_e1, "f_e1  (e1-trajectory Qhat)")
Lg  = analyze(g_p,  "g_p   (combined P-row)")
Lf  = analyze(f_q,  "f_q   (combined Q-row)")

print("="*70)
print("CROSS CHECKS")
if Le1 is not None:
    for nm, seq in (("g_p", g_p), ("f_q", f_q)):
        bad = apply_op(Le1, seq)
        print("  does %s satisfy f_e1's operator?  %s"
              % (nm, "YES (all k)" if bad is None else "NO (fails at k=%d)" % bad))
if Lf is not None and Lg is not None:
    # compare normalized operators
    same = (Lf.order() == Lg.order() and
            all(Lf[i]*Lg[Lg.order()] == Lg[i]*Lf[Lf.order()]
                for i in range(Lf.order()+1)))
    print("  Lf == Lg (up to normalization)?  %s" % ("YES" if same else "NO"))


def primitive_integral_operator(L):
    """Clear rational denominators/content and fix the leading sign.

    The normalization is global across all shift coefficients, so it changes
    an operator only by a nonzero rational scalar and keeps the recurrence
    exactly the same.
    """
    denominator = ZZ.one()
    for i in range(L.order()+1):
        for coefficient in R(L[i]).coefficients():
            denominator = lcm(denominator, coefficient.denominator())
    polys = [R(denominator * L[i]) for i in range(L.order()+1)]
    content = ZZ.zero()
    for poly in polys:
        for coefficient in poly.coefficients():
            content = gcd(content, abs(ZZ(coefficient)))
    polys = [R(poly / content) for poly in polys]
    if polys[-1].leading_coefficient() < 0:
        polys = [-poly for poly in polys]
    return A(polys)


def poincare_polynomial(L):
    cs = [R(L[i]) for i in range(L.order()+1)]
    degree = max(c.degree() for c in cs)
    xi = polygen(QQ, 'xi')
    return sum(c[degree]*xi**i for i, c in enumerate(cs))


def append_operator(lines, label, description, L):
    L = primitive_integral_operator(L)
    lines.extend([
        "",
        "== %s ==" % label,
        "# %s" % description,
        "order = %d" % L.order(),
        "degree = %d" % L.degree(),
        "normalization = primitive integral coefficients; leading coefficient of c_%d is positive"
        % L.order(),
        "Poincare polynomial = %s" % poincare_polynomial(L).factor(),
    ])
    for i in range(L.order()+1):
        lines.append("c_%d(k) = %s" % (i, R(L[i])))
    return L


if Lf is not None and Lg is not None:
    print("="*70)
    t1 = time.time()
    Lcommon = Lf.lclm(Lg)
    print("LCLM(L_f,L_g): order %d, degree %d  (%.1fs)"
          % (Lcommon.order(), Lcommon.degree(), time.time()-t1))
    badf = apply_op(Lcommon, f_q)
    badg = apply_op(Lcommon, g_p)
    print("  exact verification: f %s; g %s"
          % ("PASS" if badf is None else "FAIL at k=%d" % badf,
             "PASS" if badg is None else "FAIL at k=%d" % badg))
    assert Lcommon.order() == 9
    assert badf is None and badg is None

    lines = [
        "# Problem 2.5: exact k-recurrences for the Delannoy coefficients",
        "# Generated by scripts/p25_k_recurrence_guess2.sage with SageMath 10.9",
        "# and ore_algebra from 401 exact CMF terms (k = 0,...,400).",
        "#",
        "# B(N,k) = 2^k binomial(2k,k) binomial(N,k) binomial(N+k,k).",
        "# Qhat_N = sum_{k=0}^N f(k) B(N,k), with q = (33750,-36000,9000).",
        "# Phat_N = sum_{k=0}^N g(k) B(N,k), with p = (30921,-32972,8240).",
        "#",
        "# Convention: Sk*s(k) = s(k+1).  For L = sum_i c_i(k) Sk^i,",
        "# the recurrence is sum_i c_i(k) s(k+i) = 0.",
        "# Every displayed operator is normalized to primitive integral",
        "# coefficients with positive leading coefficient in its highest shift.",
        "# L_f and L_g were returned by ore_algebra.guess; each was checked",
        "# exactly on every available recurrence row k = 0,...,392.",
        "# The LCLM was checked exactly on f and g for k = 0,...,391.",
        "",
        "f(0..2) = %s" % [f_q[i] for i in range(3)],
        "g(0..2) = %s" % [g_p[i] for i in range(3)],
    ]
    Lf_out = append_operator(
        lines, "L_f", "minimal guessed annihilator of f(k)", Lf)
    Lg_out = append_operator(
        lines, "L_g", "minimal guessed annihilator of g(k)", Lg)
    Lcommon_out = append_operator(
        lines, "LCLM(L_f,L_g)", "common order-9 annihilator of f(k) and g(k)",
        Lcommon)

    assert apply_op(Lf_out, f_q) is None
    assert apply_op(Lg_out, g_p) is None
    assert apply_op(Lcommon_out, f_q) is None
    assert apply_op(Lcommon_out, g_p) is None

    output = (Path(__file__).resolve().parent.parent /
              "lean/specs/p25_k_recurrence_coefficients.txt")
    output.write_text("\n".join(lines) + "\n")
    print("wrote exact coefficients to", output)
print("total %.1fs" % (time.time()-t0))
