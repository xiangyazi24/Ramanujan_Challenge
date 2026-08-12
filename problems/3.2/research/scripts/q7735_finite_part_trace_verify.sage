#!/usr/bin/env sage -python
"""Exact symbolic verifier for Q7735 finite-part transport.

Run with Sage's Python mode:

    sage -python problems/3.2/research/scripts/q7735_finite_part_trace_verify.sage

This does NOT fit a recurrence.  It constructs R_n(X) in QQ(X), extracts every
regular finite part by differentiation, independently reconstructs the T/q
principal parts, checks the universal multiplier transport at every old pole
and the new-pole boundary formula, checks a quotient-ring trace representative,
and verifies the stated recurrence as exact rational arithmetic on a requested
initial range.

Important normalization guard: the harmonic summand displayed in Q7735 is not
minus the finite part at the same pole.  Only the SUM of those harmonic
summands equals minus the SUM of the regular finite parts.  The verifier checks
that trace identity separately and never uses the false termwise statement.

It also constructs the unique degree < 2 deg(D_n) Hermite polynomial P_n with
  P_n(a)  = D_n''(a)/(2 D_n'(a)),
  P_n'(a) = D_n'''(a)/(3 D_n'(a)) - D_n''(a)^2/(4 D_n'(a)^2),
so that W_n=D_n'/D_n-P_n has Laurent expansion 1/(X-a)+O((X-a)^2)
at every root a.  This verifies the finite-part extraction kernel and its
Riccati congruence modulo D_n^2 without a harmonic WZ certificate.
"""

from sage.all import QQ, PolynomialRing, binomial
import argparse


def P_apery(n):
    return 34*n**3 + 51*n**2 + 27*n + 5


def harmonic(n, power=1):
    return sum(QQ(1) / QQ(k**power) for k in range(1, n+1))


def data(n, PR, x):
    D = PR.one()
    for a in range(n+1):
        D *= x-a
    U = PR.one()
    for u in range(1, n+1):
        U *= x+u
    return U, D


def Rfun(n, K, x):
    U, D = data(n, x.parent(), x)
    return K(U) / K(D)


def finite_part_from_rational(F, a, K, x):
    """Regular finite part at a known double pole a."""
    z2F = (x-a)**2 * F
    return z2F.derivative(x, 2)(a) / 2


def local_T_q(n, a):
    T = QQ(binomial(n,a)**2 * binomial(n+a,a)**2)
    q = harmonic(n+a) + harmonic(n-a) - 2*harmonic(a)
    return T, q


def harmonic_trace_summand(n, a):
    """The displayed Q7735 summand; only its total is a finite-part trace."""
    T, q = local_T_q(n, a)
    bracket = (-harmonic(a,2) - harmonic(n-a,2)
               + 2*q*(harmonic(a)-harmonic(n-a)))
    return T*bracket


def hermite_P(D, PR, x):
    """Unique P of degree < 2 deg D with the required value/derivative jets."""
    mod = D**2
    residues = []
    moduli = []
    n = D.degree()-1
    D1, D2, D3 = D.derivative(), D.derivative(2), D.derivative(3)
    for a in range(n+1):
        da = D1(a)
        h = D2(a)/(2*da)
        k = D3(a)/(3*da) - D2(a)**2/(4*da**2)
        residues.append(PR(h + k*(x-a)))
        moduli.append((x-a)**2)
    Pcur = PR.zero()
    Mcur = PR.one()
    for res, m in zip(residues, moduli):
        rhs = (res-Pcur) % m
        inv = (Mcur % m).inverse_mod(m)
        t = (rhs*inv) % m
        Pcur += Mcur*t
        Mcur *= m
        Pcur %= Mcur
    assert Mcur == mod
    return Pcur % mod


def check_one(n, verbose=False):
    PR = PolynomialRing(QQ, 'x')
    x = PR.gen()
    K = PR.fraction_field()
    Rn = Rfun(n, K, x)
    Fn = Rn**2
    U, D = data(n, PR, x)

    # 1. The displayed harmonic formula equals the NEGATIVE TOTAL finite-part
    # trace.  It is deliberately not asserted pole-by-pole.
    harmonic_e = sum(harmonic_trace_summand(n,a) for a in range(n+1))
    fps = [finite_part_from_rational(Fn,a,K,x) for a in range(n+1)]
    e = -sum(fps)
    assert harmonic_e == e, ("global finite-part equivalence",n,harmonic_e,e)
    if n == 1:
        # Explicit regression against the tempting false termwise statement.
        assert [-v for v in fps] == [QQ(-8),QQ(-5)]
        assert [harmonic_trace_summand(1,a) for a in range(2)] == [QQ(-5),QQ(-8)]

    # 2. Quotient-ring trace representative for the finite parts themselves.
    D1, D2, D3 = D.derivative(), D.derivative(2), D.derivative(3)
    U1, U2 = U.derivative(), U.derivative(2)
    invU = U.inverse_mod(D)
    invD1 = D1.inverse_mod(D)
    Au = (U1*invU) % D
    Bd = (D2*invD1) % D
    fpClass = ((U*invD1)**2 *
        ((U2*invU) + Au**2 - 2*Au*Bd + QQ(3)/4*Bd**2
         - QQ(1)/3*(D3*invD1))) % D
    trace_fp = sum(fpClass(a) for a in range(n+1))
    assert trace_fp == sum(fps) == -e, ("quotient trace",n,trace_fp,e)
    for a in range(n+1):
        assert fpClass(a) == fps[a], ("quotient local finite part",n,a)

    # 3. Hermite finite-part extraction kernel.
    HP = hermite_P(D, PR, x)
    W = K(D1)/K(D) - K(HP)
    for a in range(n+1):
        reg = W - K.one()/(x-a)
        assert reg(a) == 0, ("Hermite constant jet",n,a)
        assert reg.derivative(x)(a) == 0, ("Hermite linear jet",n,a)
        G = Fn*W
        residue = (((x-a)**3 * G).derivative(x,2)(a))/2
        assert residue == fps[a], ("Hermite extraction",n,a)
    ric = D2 - 2*HP*D1 + D*(HP**2-HP.derivative())
    assert ric % (D**2) == 0, ("Riccati congruence",n,ric % (D**2))

    # 4. Exact n -> n+1 and n -> n+2 multiplier transport at every OLD pole.
    # eps_{n,a} := -FP_a(F_n).  These epsilons sum to e_n, but are NOT the
    # displayed harmonic summands term-by-term.
    c = QQ(n+1); d = QQ(n+2)
    Ac = K((x+c)**2) / K((x-c)**2)
    Ad = K((x+d)**2) / K((x-d)**2)
    Bmul = Ac*Ad
    F1 = Ac*Fn
    F2 = Bmul*Fn
    eps0 = [-v for v in fps]
    for a in range(n+1):
        T,q = local_T_q(n,a)
        E1transport = (Ac(a)*eps0[a] - 2*Ac.derivative(x)(a)*q*T
                       - Ac.derivative(x,2)(a)*T/2)
        E2transport = (Bmul(a)*eps0[a] - 2*Bmul.derivative(x)(a)*q*T
                       - Bmul.derivative(x,2)(a)*T/2)
        eps1 = -finite_part_from_rational(F1,a,K,x)
        eps2 = -finite_part_from_rational(F2,a,K,x)
        assert E1transport == eps1, ("n+1 finite-part transport",n,a)
        assert E2transport == eps2, ("n+2 finite-part transport",n,a)

    # 5. New-pole boundary formula: if F is analytic at c then
    # FP_c[((x+c)/(x-c))^2 F] = F(c)+4c F'(c)+2c^2 F''(c).
    fp1c = finite_part_from_rational(F1, n+1, K, x)
    boundary1 = Fn(c)+4*c*Fn.derivative(x)(c)+2*c**2*Fn.derivative(x,2)(c)
    assert fp1c == boundary1
    fp2d = finite_part_from_rational(F2, n+2, K, x)
    boundary2 = F1(d)+4*d*F1.derivative(x)(d)+2*d**2*F1.derivative(x,2)(d)
    assert fp2d == boundary2

    # The NEW double-pole coefficient is four central-binomial squares.
    T1c,_ = local_T_q(n+1,n+1)
    central = QQ(binomial(2*n+1,n)**2)
    assert T1c == 4*central
    assert Fn(c) == central/c**2

    # 6. Exact old-pole trace normal form for the requested recurrence.
    Kmul = c**3 - QQ(P_apery(n+1))*Ac + d**3*Bmul
    old_trace = QQ.zero()
    direct_old = QQ.zero()
    for a in range(n+1):
        T,q = local_T_q(n,a)
        old_trace += (Kmul(a)*eps0[a]
            - 2*Kmul.derivative(x)(a)*q*T
            - Kmul.derivative(x,2)(a)*T/2)
        direct_old += (c**3*eps0[a]
            - QQ(P_apery(n+1))*(-finite_part_from_rational(F1,a,K,x))
            + d**3*(-finite_part_from_rational(F2,a,K,x)))
    assert old_trace == direct_old

    eps1c = -fp1c
    eps2c = -finite_part_from_rational(F2,n+1,K,x)
    eps2d = -fp2d
    rhs = -(QQ(287*n**2+813*n+578)/QQ(n+2))*central
    assert old_trace - QQ(P_apery(n+1))*eps1c + d**3*(eps2c+eps2d) == rhs, (
        "full recurrence",n)

    # Independent total check from the harmonic definition, not from the old
    # trace decomposition.
    e1 = sum(harmonic_trace_summand(n+1,a) for a in range(n+2))
    e2 = sum(harmonic_trace_summand(n+2,a) for a in range(n+3))
    assert c**3*e - QQ(P_apery(n+1))*e1 + d**3*e2 == rhs

    if verbose:
        print("n",n,"e",e,"old_trace",old_trace,"rhs",rhs)
    return e


def symbolic_multiplier_checks():
    PR = PolynomialRing(QQ, names=('c','x'))
    c,x = PR.gens(); K=PR.fraction_field(); d=c+1
    A = K((x+c)**2)/K((x-c)**2)
    logder = A.derivative(x)/A
    assert logder == K(4*c)/(c**2-x**2)
    second = A.derivative(x,2)/A
    assert second == K(8*c*(2*c+x))/(c**2-x**2)**2

    B = A * K((x+d)**2)/K((x-d)**2)
    Kap = c**3 - (34*c**3+51*c**2+27*c+5)*A + d**3*B
    num,den = Kap.numerator(),Kap.denominator()
    assert num.degree(x) == den.degree(x)
    lc_ratio = num.coefficient({x:num.degree(x)}) / den.coefficient({x:den.degree(x)})
    assert lc_ratio == -4*(2*c+1)**3
    print("SYMBOLIC_MULTIPLIER PASS")


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--N',type=int,default=10)
    ap.add_argument('--verbose',action='store_true')
    args=ap.parse_args()
    symbolic_multiplier_checks()
    vals=[]
    for n in range(args.N+1):
        vals.append(check_one(n,args.verbose))
    print("EXACT_N_RANGE",0,args.N)
    print("E_PREFIX",vals[:min(len(vals),6)])
    print("Q7735_SYMBOLIC_VERIFY PASS")


if __name__=='__main__':
    main()
