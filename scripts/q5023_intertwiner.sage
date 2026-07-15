# Q5023: exact Ore-intertwiner search for Problem 2.7 versus Zudilin Thm. 2
# Run with: sage scripts/q5023_intertwiner.sage

R.<n> = PolynomialRing(QQ)
K = R.fraction_field()
RT.<T> = PolynomialRing(QQ)
SNAME = 'Sn'


def shift(f, k):
    f = K(f)
    if k == 0:
        return f
    return K(f.numerator()(n + k)) / K(f.denominator()(n + k))


def trim(op):
    op = [K(c) for c in op]
    while len(op) > 1 and op[-1] == 0:
        op.pop()
    return op


def ore_mul(P, Q):
    P = trim(P); Q = trim(Q)
    out = [K(0)] * (len(P) + len(Q) - 1)
    for i, a in enumerate(P):
        if a == 0:
            continue
        for j, b in enumerate(Q):
            if b != 0:
                out[i+j] += a * shift(b, i)
    return trim(out)


def right_divmod(F, G):
    F = trim(F[:]); G = trim(G)
    m = len(G) - 1
    if G[-1] == 0:
        raise ValueError('zero leading coefficient')
    Q = [K(0)] * max(1, len(F) - m)
    while len(F) - 1 >= m and not (len(F) == 1 and F[0] == 0):
        ell = len(F) - 1
        k = ell - m
        q = F[-1] / shift(G[-1], k)
        Q[k] += q
        for j in range(m + 1):
            F[k+j] -= q * shift(G[j], k)
        F = trim(F)
    return trim(Q), trim(F)


def ore_sub(P, Q):
    L = max(len(P), len(Q))
    out = [K(0)] * L
    for i in range(L):
        out[i] = (P[i] if i < len(P) else 0) - (Q[i] if i < len(Q) else 0)
    return trim(out)


def op_string(op):
    terms = []
    for i, c in enumerate(op):
        if c == 0:
            continue
        num = factor(c.numerator())
        den = factor(c.denominator())
        cs = str(num) if c.denominator() == 1 else '(' + str(num) + ')/(' + str(den) + ')'
        if i == 0:
            terms.append('(' + cs + ')')
        elif i == 1:
            terms.append('(' + cs + ')*' + SNAME)
        else:
            terms.append('(' + cs + ')*' + SNAME + '^' + str(i))
    return ' + '.join(terms) if terms else '0'


# Target coefficients.
def A(x):
    return 1024*(2*x+5)^4*(2*x+7)^3*(2*x+9)^3*(946*x^2+6407*x+10860)

def B(x):
    return 128*(2*x+7)^3*(2*x+9)^3*(104060*x^6+1745370*x^5+12145238*x^4+44886481*x^3+92943995*x^2+102256019*x+46709052)

def C(x):
    return 16*(x+3)^4*(2*x+9)^3*(3784*x^5+57792*x^4+351019*x^3+1059230*x^2+1587211*x+944620)

def D(x):
    return (x+3)^4*(x+4)^6*(946*x^2+4515*x+5399)

L27 = [K(-D(n)), K(C(n+1)), K(-B(n+2)), K(A(n+2))]

# Zudilin Theorem 2, shifted as in the question.
def zp3(x):
    return 2*(946*x^2-731*x+153)*(2*x+1)*(x+1)^3

def zp2(x):
    return -2*(104060*x^6+127710*x^5+12788*x^4-34525*x^3-8482*x^2+3298*x+1071)

def zp1(x):
    return 2*(3784*x^5-1032*x^4-1925*x^3+853*x^2+328*x-184)*x

def zp0(x):
    return -(946*x^2+1161*x+368)*x*(x-1)^3

LZ = [K(zp0(n+2)), K(zp1(n+2)), K(zp2(n+2)), K(zp3(n+2))]
# If z_n satisfies LZ, then y_n=64^{-n}z_n satisfies sum 64^k a_k(n)y_{n+k}=0.
LZ64 = [K((64^i) * LZ[i]) for i in range(4)]


def characteristic_at_infinity(L):
    deg = max(c.numerator().degree() - c.denominator().degree() for c in L)
    coeffs = []
    for c in L:
        num, den = c.numerator(), c.denominator()
        if den.degree() != 0:
            raise ValueError('unexpected rational coefficient in input operator')
        coeffs.append(num[deg] / den[0] if num.degree() == deg else QQ(0))
    return sum(coeffs[i]*T^i for i in range(len(coeffs)))


def limit_at_infinity(f):
    f = K(f)
    dn = f.numerator().degree(); dd = f.denominator().degree()
    if dn < dd:
        return QQ(0)
    if dn > dd:
        return Infinity if f.numerator().leading_coefficient()/f.denominator().leading_coefficient() > 0 else -Infinity
    return f.numerator().leading_coefficient()/f.denominator().leading_coefficient()


print('=== Exact input diagnostics ===')
chi27 = characteristic_at_infinity(L27)
chiZ = characteristic_at_infinity(LZ)
chiZ64 = characteristic_at_infinity(LZ64)
print('chi27(T)  =', factor(chi27))
print('chiZ(T)   =', factor(chiZ))
print('chiZ64(T) =', factor(chiZ64))
print('gcd(chi27,chiZ)   =', gcd(chi27, chiZ))
print('gcd(chi27,chiZ64) =', gcd(chi27, chiZ64))
print('chi27 == chiZ(64*T):', bool(chi27 == chiZ(64*T)))

# Determinant ratio for companion systems. A rational matrix gauge U would
# require ratio = det(U(n+1))/det(U(n)), whose limit at infinity is 1.
det27 = -L27[0]/L27[3]
detZ = -LZ[0]/LZ[3]
detZ64 = -LZ64[0]/LZ64[3]
ratio_raw = K(det27/detZ)
ratio_scaled = K(det27/detZ64)
print('det27/detZ =', factor(ratio_raw))
print('limit det27/detZ at infinity =', limit_at_infinity(ratio_raw))
print('det27/detZ64 =', factor(ratio_scaled))
print('limit det27/detZ64 at infinity =', limit_at_infinity(ratio_scaled))
print()


def basis_remainders(source, d, Dmax):
    rems = []
    labels = []
    for i in range(d+1):
        for j in range(Dmax+1):
            E = [K(0)]*(i+1)
            E[i] = K(n^j)
            F = ore_mul(L27, E)
            q, r = right_divmod(F, source)
            r += [K(0)]*(3-len(r))
            rems.append(r[:3])
            labels.append((i,j))
    return labels, rems


def exact_remainder_matrix(rems):
    cols = len(rems)
    rows = []
    denominator_degrees = []
    for k in range(3):
        den = R(1)
        for b in range(cols):
            den = lcm(den, rems[b][k].denominator())
        denominator_degrees.append(den.degree())
        polys = [R(rems[b][k] * den) for b in range(cols)]
        nonzero = [p for p in polys if p != 0]
        maxdeg = max(p.degree() for p in nonzero) if nonzero else -1
        for e in range(maxdeg+1):
            row = [p[e] for p in polys]
            if any(c != 0 for c in row):
                rows.append(row)
    return matrix(QQ, rows), denominator_degrees


def vector_to_operator(v, labels, d):
    op = [K(0)]*(d+1)
    for c, (i,j) in zip(v, labels):
        if c:
            op[i] += K(c*n^j)
    return trim(op)


def verify_intertwiner(source, Rop):
    product = ore_mul(L27, Rop)
    Qop, rem = right_divmod(product, source)
    rem += [K(0)]*(3-len(rem))
    okrem = all(x == 0 for x in rem[:3])
    diff = ore_sub(product, ore_mul(Qop, source))
    okid = all(x == 0 for x in diff)
    return okrem and okid, Qop, rem, diff


def search(source, label, dmax=3, Dmax=20):
    print('=== Polynomial intertwiner search:', label, '===')
    for d in range(dmax+1):
        labels, rems = basis_remainders(source, d, Dmax)
        M, dd = exact_remainder_matrix(rems)
        rank = M.rank()
        nullity = M.ncols() - rank
        print('d=%d, D<=%d: matrix %d x %d, rank=%d, nullity=%d, common-den degrees=%s' %
              (d,Dmax,M.nrows(),M.ncols(),rank,nullity,dd))
        if nullity == 0:
            continue
        ker = M.right_kernel().basis()
        for idx, v in enumerate(ker):
            Rop = vector_to_operator(v, labels, d)
            ok, Qop, rem, diff = verify_intertwiner(source, Rop)
            print('  kernel vector', idx, 'symbolic verification:', ok)
            if ok and Rop != [0]:
                print('FOUND R of order', len(Rop)-1)
                print('R =', op_string(Rop))
                print('Q =', op_string(Qop))
                print('Exact identity verified:', all(x == 0 for x in diff))
                return Rop, Qop
    print('NO polynomial-coefficient intertwiner in the requested box.')
    return None

raw = search(LZ, 'raw LZ', dmax=3, Dmax=20)
print()
scaled = search(LZ64, 'exponentially normalized LZ64 (z_n -> 64^{-n} z_n)', dmax=3, Dmax=20)
print()

if raw is None:
    print('RAW CONCLUSION: no R with polynomial coefficients, order <= 3 and coefficient degree <= 20.')
if scaled is None:
    print('SCALED CONCLUSION: no polynomial R in the same requested box after the necessary 64^{-n} normalization.')
else:
    print('SCALED CONCLUSION: exact rational Ore intertwiner found after the necessary exponential normalization.')
