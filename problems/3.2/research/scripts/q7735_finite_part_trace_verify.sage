#!/usr/bin/env sage -python
"""Exact symbolic verifier for Q7735 finite-part transport.

Run with Sage's Python mode:

    sage -python problems/3.2/research/scripts/q7735_finite_part_trace_verify.sage

The verifier does not fit a recurrence.  It constructs R_n(X) in QQ(X),
extracts every finite part by differentiation, independently reconstructs the
T/q partial-fraction data, checks the universal multiplier transport at old
poles and the new-pole boundary formula, checks the quotient-ring/Hermite
trace representatives, and finally checks the stated inhomogeneous recurrence
as an exact rational identity for a configurable initial range of n.

It also constructs the unique degree < 2 deg(D_n) Hermite polynomial P_n with
  P_n(a)  = D_n''(a)/(2 D_n'(a)),
  P_n'(a) = D_n'''(a)/(3 D_n'(a)) - D_n''(a)^2/(4 D_n'(a)^2),
so that W_n=D_n'/D_n-P_n has Laurent expansion 1/(X-a)+O((X-a)^2)
at every root a.  This verifies the finite-part extraction kernel and its
Riccati congruence modulo D_n^2 without any harmonic WZ certificate.
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
    """Finite part at a known double pole a."""
    z2F = (x-a)**2 * F
    # If z2F is analytic, FP is half its second derivative at a.
    return z2F.derivative(x, 2)(a) / 2


def local_T_q(n, a):
    T = QQ(binomial(n,a)**2 * binomial(n+a,a)**2)
    q = harmonic(n+a) + harmonic(n-a) - 2*harmonic(a)
    return T, q


def E_formula(n, a):
    T, q = local_T_q(n, a)
    bracket = (-harmonic(a,2) - harmonic(n-a,2)
               + 2*q*(harmonic(a)-harmonic(n-a)))
    return T*bracket


def hermite_P(D, PR, x):
    """Unique P of degree < 2 deg D with the required value/derivative jets."""
    mod = D**2
    # Chinese remainder construction modulo (x-a)^2 for roots 0,...,n.
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
    # Pairwise-coprime CRT by accumulation.
    Pcur = PR.zero()
    Mcur = PR.one()
    for res, m in zip(residues, moduli):
        # Pcur + Mcur*t == res mod m.
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

    # 1. Finite part = negative of the supplied E summand.
    e = QQ.zero()
    for a in range(n+1):
        fp = finite_part_from_rational(Fn, a, K, x)
        Ef = E_formula(n,a)
        assert fp == -Ef, ("finite part formula", n, a, fp, Ef)
        e += Ef

    # 2. Quotient-ring trace representative, with all inverses taken mod D.
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
    assert trace_fp == -e, ("quotient trace", n, trace_fp, e)

    # 3. Hermite finite-part extraction kernel.
    HP = hermite_P(D, PR, x)
    W = K(D1)/K(D) - K(HP)
    for a in range(n+1):
        # W - 1/(x-a) has zero value and derivative at a.
        reg = W - K.one()/(x-a)
        assert reg(a) == 0, ("Hermite constant jet",n,a)
        assert reg.derivative(x)(a) == 0, ("Hermite linear jet",n,a)
        # Therefore Res(F*W,a)=FP(F,a).  Extract residue by one derivative
        # after multiplying by (x-a)^3, since F*W has order <=3.
        G = Fn*W
        residue = (((x-a)**3 * G).derivative(x,2)(a))/2
        assert residue == finite_part_from_rational(Fn,a,K,x), (
            "Hermite extraction",n,a)
    ric = D2 - 2*HP*D1 + D*(HP**2-HP.derivative())
    assert ric % (D**2) == 0, ("Riccati congruence",n,ric % (D**2))

    # 4. Exact n -> n+1 and n -> n+2 multiplier transport at every old pole.
    c = QQ(n+1); d = QQ(n+2)
    Ac = K((x+c)**2) / K((x-c)**2)
    Ad = K((x+d)**2) / K((x-d)**2)
    Bmul = Ac*Ad
    F1 = Ac*Fn
    F2 = Bmul*Fn
    for a in range(n+1):
        T,q = local_T_q(n,a)
        E0 = E_formula(n,a)
        E1transport = (Ac(a)*E0 - 2*Ac.derivative(x)(a)*q*T
                       - Ac.derivative(x,2)(a)*T/2)
        E2transport = (Bmul(a)*E0 - 2*Bmul.derivative(x)(a)*q*T
                       - Bmul.derivative(x,2)(a)*T/2)
        assert E1transport == E_formula(n+1,a), ("n+1 transport",n,a)
        assert E2transport == E_formula(n+2,a), ("n+2 transport",n,a)

    # 5. New-pole boundary formula: if F is analytic at c then
    # FP_c[((x+c)/(x-c))^2 F] = F(c)+4c F'(c)+2c^2 F''(c).
    fp1c = finite_part_from_rational(F1, n+1, K, x)
    boundary1 = Fn(c)+4*c*Fn.derivative(x)(c)+2*c**2*Fn.derivative(x,2)(c)
    assert fp1c == boundary1
    # and analogously at d for F_{n+1}.
    fp2d = finite_part_from_rational(F2, n+2, K, x)
    boundary2 = F1(d)+4*d*F1.derivative(x)(d)+2*d**2*F1.derivative(x,2)(d)
    assert fp2d == boundary2

    # The new double-pole coefficient is 4 times the central-binomial square.
    T1c,_ = local_T_q(n+1,n+1)
    central = QQ(binomial(2*n+1,n)**2)
    assert T1c == 4*central
    assert Fn(c) == central/c**2

    # 6. The exact old-pole trace normal form for the requested recurrence.
    Kmul = c**3 - QQ(P_apery(n+1))*Ac + d**3*Bmul
    old_trace = QQ.zero()
    for a in range(n+1):
        T,q = local_T_q(n,a)
        old_trace += (Kmul(a)*E_formula(n,a)
            - 2*Kmul.derivative(x)(a)*q*T
            - Kmul.derivative(x,2)(a)*T/2)
    direct_old = sum(c**3*E_formula(n,a)
        - QQ(P_apery(n+1))*E_formula(n+1,a)
        + d**3*E_formula(n+2,a) for a in range(n+1))
    assert old_trace == direct_old

    E1c = E_formula(n+1,n+1)
    E2c = E_formula(n+2,n+1)
    E2d = E_formula(n+2,n+2)
    rhs_boundary = -(QQ(287*n**2+813*n+578)/QQ(n+2))*central
    assert old_trace - QQ(P_apery(n+1))*E1c + d**3*(E2c+E2d) == rhs_boundary, (
        "full recurrence",n)

    if verbose:
        print("n",n,"e",e,"old_trace",old_trace,"boundary_rhs",rhs_boundary)
    return e


def symbolic_multiplier_checks():
    # Pure rational-function identities with symbolic c,x, independent of any
    # finite-n recurrence data.
    PR = PolynomialRing(QQ, names=('c','x'))
    c,x = PR.gens(); K=PR.fraction_field(); d=c+1
    A = K((x+c)**2)/K((x-c)**2)
    logder = A.derivative(x)/A
    assert logder == K(4*c)/(c**2-x**2)
    second = A.derivative(x,2)/A
    assert second == K(8*c*(2*c+x))/(c**2-x**2)**2

    B = A * K((x+d)**2)/K((x-d)**2)
    Kap = c**3 - (34*c**3+51*c**2+27*c+5)*A + d**3*B
    # The value at infinity is the compact cube -4(2c+1)^3.
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
