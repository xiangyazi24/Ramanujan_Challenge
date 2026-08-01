#!/usr/bin/env sage-python
"""Guess scalar recurrences for the two exact Lima bracket margins modulo p."""

from sage.all import GF, matrix


prime = 2305843009213693951
field = GF(prime)
terms = 900


def inv(value):
    return field(value) ** -1


def challenge(n):
    raw = [
        (2*n+5)*(n+3)**2*(136*n**4+1424*n**3+5548*n**2+9551*n+6141),
        384*n**6+6384*n**5+44168*n**4+162698*n**3+336377*n**2+369933*n+169011,
        480*n**4+4980*n**3+19210*n**2+32690*n+20730,
        (n+2)**2*(n+3)**2*(4*n+10)*(48*n**3+386*n**2+1017*n+879),
        (n+2)**2*(272*n**5+3848*n**4+21732*n**3+61184*n**2+85761*n+47808),
        (n+2)**2*(320*n**3+2540*n**2+6610*n+5640),
        (4*n+10)*(n+2)**2*(n+3)**2*(32*n**4+302*n**3+1037*n**2+1530*n+813),
        (n+2)**2*(192*n**6+2984*n**5+19116*n**4+64452*n**3+120256*n**2+117279*n+46476),
        (n+2)**2*(16*n**5+408*n**4+2912*n**3+8884*n**2+12254*n+6240),
    ]
    gauge = 2*(n+2)**2*(n+3)**2*(2*n+5)*(2*n+7)**2
    return matrix(field, 3, 3, [field(value)/field(gauge) for value in raw])


def ratio(m):
    return -field((m+1)**3*(3*m+5))/field((2*m+3)**3*(3*m+2))


def lima_term(m):
    # Recurrence is used below; only a few seed terms are needed.
    answer = field(1)
    for index in range(m):
        answer *= ratio(index)
    return answer


partial = field(0)
for m in range(4):
    partial += lima_term(m)
c4 = lima_term(4)

p = [field(30921), field(32972), field(8240)]
q = [field(33750), field(36000), field(9000)]
lower = [q[j]*partial-p[j] for j in range(3)]
tail = [q[j]*c4 for j in range(3)]
upper = [-lower[j]-tail[j] for j in range(3)]

lower_values = []
upper_values = []
for n in range(terms):
    lower_values.append(lower[2])
    upper_values.append(upper[0])
    C = challenge(n)
    m = 2*n+4
    a = ratio(m)
    s = a*ratio(m+1)
    r = 1+a
    t = -(a+s)
    lower = [sum((lower[i]+r*tail[i])*C[i, j] for i in range(3)) for j in range(3)]
    upper = [sum((upper[i]+t*tail[i])*C[i, j] for i in range(3)) for j in range(3)]
    tail = [s*sum(tail[i]*C[i, j] for i in range(3)) for j in range(3)]


def search(values, name):
    print("search", name, flush=True)
    order = 6
    candidates = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    for degree in candidates:
        unknowns = (order+1)*(degree+1)
        rows = []
        for n in range(unknowns+8):
            powers = [field(1)]
            for _ in range(degree):
                powers.append(powers[-1]*field(n))
            rows.append([values[n+shift]*powers[d]
                         for shift in range(order+1) for d in range(degree+1)])
        kernel = matrix(field, rows).right_kernel()
        print("degree", degree, "nullity", kernel.dimension(), flush=True)
        if kernel.dimension():
            vector = kernel.basis()[0]
            # Verify against every unused term.
            good = True
            for n in range(unknowns+8, len(values)-order):
                total = field(0)
                for shift in range(order+1):
                    polynomial = field(0)
                    power = field(1)
                    for d in range(degree+1):
                        polynomial += vector[shift*(degree+1)+d]*power
                        power *= field(n)
                    total += polynomial*values[n+shift]
                if total:
                    good = False
                    print("failed", n, flush=True)
                    break
            print("verified", good, flush=True)
            if good:
                return degree, vector
    return None


for values, name in [(lower_values, "lower"), (upper_values, "upper")]:
    answer = search(values, name)
    if answer is not None:
        degree, vector = answer
        print("FOUND", name, degree, flush=True)
        for shift in range(7):
            coefficients = list(vector[shift*(degree+1):(shift+1)*(degree+1)])
            while coefficients and not coefficients[-1]:
                coefficients.pop()
            print(shift, len(coefficients)-1, [int(value) for value in coefficients], flush=True)
