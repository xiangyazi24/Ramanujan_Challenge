from sage.all import *
from ore_algebra import OreAlgebra

# ============================================================
# 0. Rings and helpers
# ============================================================
Rn.<n> = PolynomialRing(QQ)
Kn = Rn.fraction_field()
Rt.<t> = PolynomialRing(QQ)
Kt = Rt.fraction_field()
OD.<Dt> = OreAlgebra(Kt, 'Dt')
theta = Kt(t)*Dt


def eval_poly_at_op(f, T):
    out = OD.zero()
    for a in reversed(Rn(f).list()):
        out = out*T + Kt(a)
    return out


def primitive_coefficients(c):
    c = [Kn(a) for a in c]
    den = Rn.one()
    for a in c:
        den = lcm(den, Rn(a.denominator()))
    nums = [Rn(den*a) for a in c]
    g = nums[0]
    for a in nums[1:]:
        g = gcd(g, a)
    nums = [a.quo_rem(g)[0] for a in nums]
    scl = QQ(1)/QQ(nums[-1].leading_coefficient())
    return [Rn(scl*a) for a in nums]


def right_quo_rem(P, Q):
    P = OD(P)
    Q = OD(Q)
    quo = OD.zero()
    rem = P
    while rem != 0 and rem.order() >= Q.order():
        k = rem.order()-Q.order()
        a = rem[rem.order()]/Q[Q.order()]
        term = Kt(a)*Dt^k
        quo += term
        rem -= term*Q
    return quo, rem


def apply_op(L, f):
    f = Kt(f)
    ans = Kt(0)
    der = f
    for k in range(L.order()+1):
        if k == 0:
            der = f
        elif k == 1:
            der = f.derivative()
        else:
            der = der.derivative()
        ans += Kt(L[k])*der
    return ans


def conjugate_by_log_derivative(L, ell):
    """Return g^(-1) L g when g'/g = ell."""
    P = Dt + Kt(ell)
    out = OD.zero()
    power = OD.one()
    for k in range(L.order()+1):
        if k > 0:
            power *= P
        out += Kt(L[k])*power
    return out


# ============================================================
# 1. Original P2.7 rational recurrence
# ============================================================
def A(z):
    return 1024*(2*z+5)^4*(2*z+7)^3*(2*z+9)^3 \
           *(946*z^2+6407*z+10860)

def B(z):
    return 128*(2*z+7)^3*(2*z+9)^3 * (
        104060*z^6+1745370*z^5+12145238*z^4
        +44886481*z^3+92943995*z^2+102256019*z
        +46709052
    )

def C(z):
    return 16*(z+3)^4*(2*z+9)^3 * (
        3784*z^5+57792*z^4+351019*z^3
        +1059230*z^2+1587211*z+944620
    )

def D(z):
    return (z+3)^4*(z+4)^6*(946*z^2+4515*z+5399)


crat = [
    Kn(-D(n)/A(n)),
    Kn(C(n+1)/A(n+1)),
    Kn(-B(n+2)/A(n+2)),
    Kn(1),
]
p = primitive_coefficients(crat)
assert [f.degree() for f in p] == [18]*4

rel18 = [p[j][17]/p[j][18] for j in range(4)]
assert rel18 == [QQ(1503)/22]*4
print('relative n^17 coefficients of true recurrence =', rel18)

Ldisplay_coeffs = [Rn(-D(n)), Rn(C(n+1)), Rn(-B(n+2)), Rn(A(n+2))]
assert any(
    Ldisplay_coeffs[j]*p[3] != Ldisplay_coeffs[3]*p[j]
    for j in range(3)
)
print('non-associate display form is DIFFERENT from primitive form: CONFIRMED')

# ============================================================
# 2. Mellin adjoint and exact Euler quotient
# ============================================================
Mdag = OD.zero()
for j in range(4):
    Mdag += Kt(t^j)*eval_poly_at_op(p[j], -theta-j-1)
assert Mdag.order() == 18

Rx.<x> = PolynomialRing(QQ)
shifted = [Rx(p[j](x-j)) for j in range(4)]
gshift = shifted[0]
for f in shifted[1:]:
    gshift = gcd(gshift, f)
gshift = gshift.monic()
assert gshift == x^2 + QQ(105)/22*x + QQ(5399)/946
print('Euler factor g(n) =', gshift)

Rpoly = []
for f in shifted:
    q, rem = f.quo_rem(gshift)
    assert rem == 0
    Rpoly.append(q)

Euler2 = eval_poly_at_op(gshift, -theta-1)
assert Euler2 == theta^2-QQ(61)/22*theta+QQ(915)/473

M16 = OD.zero()
for j in range(4):
    M16 += Kt(t^j)*eval_poly_at_op(Rpoly[j], -theta-1)
assert M16.order() == 16
assert Mdag == M16*Euler2
print('M_dag = M16 * E: VERIFIED')

gn = n^2+QQ(105)/22*n+QQ(5399)/946
rrec = []
for j in range(4):
    q, rem = p[j].quo_rem(Rn(gn(n=n+j)))
    assert rem == 0
    rrec.append(q)
assert [f.degree() for f in rrec] == [16]*4
rel16 = [rrec[j][15]/rrec[j][16] for j in range(4)]
assert rel16 == [QQ(699)/11-2*j for j in range(4)]
print('relative n^15 coefficients of quotient recurrence =', rel16)

# ============================================================
# 3. Singularities and exact indicial data
# ============================================================
chi = t^3-QQ(55)/64*t^2+QQ(1)/2048*t-QQ(1)/2^20
assert M16[16] == Kt(t^16*chi)
print('\nleading coefficient =', factor(M16[16]))

Rrho.<rho> = PolynomialRing(QQ)
I0 = Rrho(Rpoly[0](-rho-1))
Iinf = Rrho(Rpoly[3](rho-1))
print('I_0 =', factor(I0))
print('I_infinity =', factor(Iinf))

# Finite cubic singularity over K=Q(mu).
RX.<X> = PolynomialRing(QQ)
Pmu = 4*X^3-220*X^2+8*X-1
Kmu.<mu> = NumberField(Pmu, embedding=RR(54.9))
lam = mu/64

def eval_Kt_at(f, value):
    f = Kt(f)
    num = Rt(f.numerator())
    den = Rt(f.denominator())
    return Kmu(num(value))/Kmu(den(value))

special = Kmu(15) - eval_Kt_at(M16[15], lam) / eval_Kt_at(Kt(M16[16]).derivative(), lam)
assert special == -3
Ilam = prod(rho-k for k in range(15))*(rho+3)
print('finite special exponent =', special)
print('I_lambda =', factor(Ilam))

# Fuchs sum
assert -Iinf[15]/Iinf[16] == QQ(71)/11
assert QQ(523)/11 + QQ(71)/11 + 3*QQ(102) == 360
print('Fuchs sum = 360: VERIFIED')

# ============================================================
# 4. No additional pure Euler factor, right or left
# ============================================================
gright = Rpoly[0]
for f in Rpoly[1:]:
    gright = gcd(gright, f)
assert gright.degree() == 0

left_shifted = [Rx(Rpoly[j](x+j)) for j in range(4)]
gleft = left_shifted[0]
for f in left_shifted[1:]:
    gleft = gcd(gleft, f)
assert gleft.degree() == 0
print('no further right or left Euler factor')

# ============================================================
# 5. Full exact differential factorization
# ============================================================
print('\n=== FULL FACTORIZATION OF M16 ===')
try:
    fac = M16.factor()
    print('factorization =', fac)
except Exception as err:
    fac = None
    print('M16.factor() failed:', repr(err))

def try_right_factors(L, order):
    for name in ('right_factors', 'right_factors_of_order'):
        if hasattr(L, name):
            meth = getattr(L, name)
            try:
                ans = meth(order)
                print('%s(order=%s) ->' % (name, order), ans)
                return ans
            except TypeError:
                try:
                    ans = meth(order=order)
                    print('%s(order=%s) ->' % (name, order), ans)
                    return ans
                except Exception as err:
                    print(name, 'failed:', repr(err))
            except Exception as err:
                print(name, 'failed:', repr(err))
    return None

right_factor_data = {}
for ord0 in [1,2,3,4]:
    right_factor_data[ord0] = try_right_factors(M16, ord0)

def verify_right_factor(L, F):
    Q, R = right_quo_rem(L, F)
    assert R == 0
    assert L == Q*F
    return Q

reported_right_factors = []
for F in reported_right_factors:
    Q = verify_right_factor(M16, F)
    print('verified right factor order', F.order(), 'left quotient order', Q.order())

# ============================================================
# 6. Independent hyperexponential search
# ============================================================
def polynomial_solutions_of_gauged_operator(Lg, maxdeg):
    """Find P(t), deg P<=maxdeg, with Lg(P)=0 by exact linear algebra."""
    vals = [apply_op(Lg, Kt(t^i)) for i in range(maxdeg+1)]
    den = Rt.one()
    for v in vals:
        den = lcm(den, Rt(Kt(v).denominator()))
    polys = [Rt(den*v) for v in vals]
    md = max((f.degree() if f else -1) for f in polys)
    M = matrix(QQ, md+1, maxdeg+1,
               lambda row, col: polys[col][row] if row <= polys[col].degree() else 0)
    ans = []
    for vec in M.right_kernel().basis():
        P = Rt(sum(vec[i]*t^i for i in range(maxdeg+1)))
        if P and apply_op(Lg, Kt(P)) == 0:
            ans.append(P)
    return ans

print('\n=== HYPEREXPONENTIAL ANSATZ SEARCH ===')
origin_rational_exponents = [QQ(2), QQ(3), QQ(5)/2, QQ(7)/2]
finite_exponents = [-3] + list(range(0,15))
MAX_POLY_DEG = 12

for rho0 in origin_rational_exponents:
    for kappa in finite_exponents:
        ell = Kt(rho0/t + kappa*chi.derivative()/chi)
        Lg = conjugate_by_log_derivative(M16, ell)
        sols = polynomial_solutions_of_gauged_operator(Lg, MAX_POLY_DEG)
        for Psol in sols:
            logder = ell + Kt(Psol.derivative()/Psol)
            F1 = Dt-Kt(logder)
            Qleft, rem = right_quo_rem(M16, F1)
            assert rem == 0
            print('FOUND order-1 factor: rho0=', rho0,
                  'finite exponent=', kappa, 'P=', factor(Psol))
            print('factor =', F1)

# ============================================================
# 7. Inverted AESZ #209 comparison
# ============================================================
def U1(T):
    return 1902*T^4+3708*T^3+2789*T^2+935*T+119

def U2(T):
    return 62408*T^4+68576*T^3-10029*T^2-24106*T-5661

def U3(T):
    return 66180*T^4+33048*T^3+20785*T^2+17799*T+4794

def U4(T):
    return 196*T^3+498*T^2+487*T+169

def P209_terms(T):
    return [
        17^2*T^4,
        -34*U1(T),
        4*U2(T),
        -4*U3(T),
        128*(2*T+1)*U4(T),
        -4096*(T+1)^2*(2*T+1)*(2*T+3),
    ]

AESZ_inv = OD.zero()
for j, term in enumerate(P209_terms(-theta)):
    AESZ_inv += Kt((256*t)^(5-j))*term
print('\ninverted AESZ order =', AESZ_inv.order())
print('inverted AESZ leading coefficient =', factor(AESZ_inv[AESZ_inv.order()]))

def try_gcrd(Aop, Bop):
    for name in ('gcrd', 'greatest_common_right_divisor'):
        if hasattr(Aop, name):
            try:
                G = getattr(Aop, name)(Bop)
                print(name, 'order =', G.order(), 'operator =', G)
                return G
            except Exception as err:
                print(name, 'failed:', repr(err))
    return None

G209 = try_gcrd(M16, AESZ_inv)
if G209 is not None and G209.order() > 0:
    verify_right_factor(M16, G209)
    verify_right_factor(AESZ_inv, G209)

# ============================================================
# 8. Recurrence-level Pochhammer gauges
# ============================================================
def rh(z):
    return Kn(2)^(-20)*(z+3)^4*(z+4)^6 / (
        (z+QQ(5)/2)^4*(z+QQ(7)/2)^3*(z+QQ(9)/2)^3
    )

def hquot(j):
    out = Kn(1)
    for k in range(j):
        out *= rh(n+k)
    return out

def mellin_from_rational_coeffs(c):
    pp = primitive_coefficients(c)
    MM = OD.zero()
    for j in range(4):
        MM += Kt(t^j)*eval_poly_at_op(pp[j], -theta-j-1)
    return pp, MM

for label, coeffs in [
    ('u_n = h_n v_n', [crat[j]*hquot(j) for j in range(4)]),
    ('u_n = v_n/h_n', [crat[j]/hquot(j) for j in range(4)]),
]:
    pp, MM = mellin_from_rational_coeffs(coeffs)
    print('\n===', label, '===')
    print('primitive degrees =', [f.degree() for f in pp])
    print('Mellin order =', MM.order())
    print('leading coefficient =', factor(MM[MM.order()]))
    try:
        print('factorization =', MM.factor())
    except Exception as err:
        print('factorization failed:', repr(err))

print('\nDone.')
