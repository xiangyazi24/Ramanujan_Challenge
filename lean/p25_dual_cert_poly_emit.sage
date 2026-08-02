from sage.all import *

Kn.<n> = PolynomialRing(QQ)
K = Kn.fraction_field()
R.<p,q,v> = PolynomialRing(K)
lines = open('/tmp/p25_dual_cert.out').read().splitlines()
certs = []
for i,line in enumerate(lines):
    if line.startswith(('Pp terms', 'Pq terms', 'Pv terms')):
        certs.append(R(lines[i+1]))
assert len(certs) == 9

S.<nn,pp,qq,vv> = PolynomialRing(QQ)
delta = 4 * (2*n+3) * (n+2)

def n_horner(poly):
    poly = Kn(poly)
    if poly == 0:
        return '0'
    coeffs = poly.list()
    ans = str(coeffs[-1])
    for c in reversed(coeffs[:-1]):
        ans = '(' + ans + ' * dualCertN + (' + str(c) + '))'
    return ans

def lean_terms(poly):
    poly = S(poly)
    terms = []
    # Group the common polynomial in n for every spatial monomial.
    spatial = {}
    for (an,ap,aq,av), coeff in poly.dict().items():
        spatial.setdefault((ap,aq,av), Kn.zero())
        spatial[(ap,aq,av)] += coeff*n**an
    for (ap,aq,av), coeff in sorted(spatial.items(), reverse=True):
        cs = ', '.join(str(c) for c in coeff.list())
        terms.append('{ nCoeffs := [%s], pExp := %d, qExp := %d, vExp := %d }' %
                     (cs, ap, aq, av))
    return ',\n    '.join(terms)

names = ['Pp0','Pq0','Pv0','Pp1','Pq1','Pv1','Pp2','Pq2','Pv2']
for name,poly in zip(names, certs):
    scaled = S.zero()
    for (ap,aq,av), coeff in (delta*poly).dict().items():
        coeff = Kn(coeff)
        coeffS = sum(S(c) * nn**i for i,c in enumerate(coeff.list()))
        scaled += coeffS * pp**ap * qq**aq * vv**av
    print('def dualCert%sTerms : List DualCertTerm := [' % name)
    print('    ' + lean_terms(scaled))
    print('  ]')
    print()
    print('def dualCert%sPoly : DualCertPoly :=' % name)
    print('  dualCertOfTerms dualCert%sTerms' % name)
    print()
