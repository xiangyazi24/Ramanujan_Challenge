#!/usr/bin/env python3
"""Temporary exact pointwise recurrence coefficients for a cross-product."""

from fractions import Fraction as F


def wilson_coefficients(m):
    a4 = 12265 + 29296*m + 26176*m**2 + 10368*m**3 + 1536*m**4
    c4 = 313 + 1904*m + 4288*m**2 + 4224*m**3 + 1536*m**4
    b = (111992515 + 1144683736*m + 5147619352*m**2
         + 13412393984*m**3 + 22433518592*m**4
         + 25185342464*m**5 + 19235018752*m**6
         + 9876373504*m**7 + 3265527808*m**8
         + 628359168*m**9 + 53477376*m**10)
    a = 4*(m+1)**2*(4*m+1)**2*(4*m+3)**2*a4
    c = 4*(m+2)**2*(4*m+5)**2*(4*m+7)**2*c4
    return F(a), F(b), F(c)


def challenge(n):
    e = [
        (-2*n-5)*(n+3)**2*(136*n**4+1424*n**3+5548*n**2+9551*n+6141),
        384*n**6+6384*n**5+44168*n**4+162698*n**3+336377*n**2+369933*n+169011,
        -(480*n**4+4980*n**3+19210*n**2+32690*n+20730),
        (n+2)**2*(n+3)**2*(4*n+10)*(48*n**3+386*n**2+1017*n+879),
        (n+2)**2*(-272*n**5-3848*n**4-21732*n**3-61184*n**2-85761*n-47808),
        (n+2)**2*(320*n**3+2540*n**2+6610*n+5640),
        (-4*n-10)*(n+2)**2*(n+3)**2*(32*n**4+302*n**3+1037*n**2+1530*n+813),
        (n+2)**2*(192*n**6+2984*n**5+19116*n**4+64452*n**3+120256*n**2+117279*n+46476),
        (n+2)**2*(-16*n**5-408*n**4-2912*n**3-8884*n**2-12254*n-6240),
    ]
    d = F(-2*(n+2)**2*(n+3)**2*(2*n+5)*(2*n+7)**2)
    return [[F(e[3*i+j], d) for j in range(3)] for i in range(3)]


def transition(n):
    a, b, c = wilson_coefficients(n + 2)
    u = [[F(0), F(1)], [-a/c, b/c]]
    m = challenge(n)
    # row-major vectorization: (U C M)_(r,j)
    return [[sum(u[r][s] * m[i][j] if old == 3*s+i else 0
                 for s in range(2) for i in range(3))
             for old in range(6)]
            for r in range(2) for j in range(3)]


def matvec(matrix, vector):
    return [sum(matrix[i][j]*vector[j] for j in range(len(vector)))
            for i in range(len(matrix))]


def matmul(left, right):
    return [[sum(left[i][k]*right[k][j] for k in range(len(right)))
             for j in range(len(right[0]))] for i in range(len(left))]


def solve(matrix, rhs):
    a = [row[:] + [rhs[i]] for i, row in enumerate(matrix)]
    for col in range(len(a)):
        pivot = next(i for i in range(col, len(a)) if a[i][col])
        a[col], a[pivot] = a[pivot], a[col]
        scale = a[col][col]
        a[col] = [x/scale for x in a[col]]
        for i in range(len(a)):
            if i != col and a[i][col]:
                scale = a[i][col]
                a[i] = [a[i][j]-scale*a[col][j] for j in range(len(a)+1)]
    return [a[i][-1] for i in range(len(a))]


def coefficients(n):
    # Column state update is z(n+1)=L(n)z(n); collect e_0^T products.
    product = [[F(int(i == j)) for j in range(6)] for i in range(6)]
    vectors = [product[0][:]]
    for shift in range(6):
        product = matmul(transition(n + shift), product)
        vectors.append(product[0][:])
    basis = [[vectors[column][row] for column in range(6)] for row in range(6)]
    return solve(basis, [-value for value in vectors[6]]) + [F(1)]


for n in list(range(15)) + [20, 30, 50, 100]:
    c = coefficients(n)
    print(n, ['+' if x > 0 else '-' if x < 0 else '0' for x in c],
          [float(x) for x in c])
