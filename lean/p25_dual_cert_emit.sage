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
delta = 1

def lean_poly(P):
    P = R(delta*P)
    terms = []
    for (a,b,c), coeff in sorted(P.dict().items(), reverse=True):
        cs = str(coeff).replace('^', ' ^ ')
        factors = ['(' + cs + ')']
        if a: factors.append('p' if a == 1 else 'p ^ %d' % a)
        if b: factors.append('q' if b == 1 else 'q ^ %d' % b)
        if c: factors.append('v' if c == 1 else 'v ^ %d' % c)
        terms.append(' * '.join(factors))
    return ' +\n    '.join(terms) if terms else '0'

names = ['Pp0','Pq0','Pv0','Pp1','Pq1','Pv1','Pp2','Pq2','Pv2']
for name,P in zip(names,certs):
    print('private def dualCert%s (n p q v : ℝ) : ℝ :=' % name)
    print('  ' + lean_poly(P))
    print()

derivatives = []
for i in range(3):
    derivatives.extend([certs[3*i].derivative(p), certs[3*i+1].derivative(q),
                        certs[3*i+2].derivative(v)])
dn = ['DPp0','DPq0','DPv0','DPp1','DPq1','DPv1','DPp2','DPq2','DPv2']
for name,P in zip(dn,derivatives):
    print('private def dualCert%s (n p q v : ℝ) : ℝ :=' % name)
    print('  ' + lean_poly(P))
    print()
