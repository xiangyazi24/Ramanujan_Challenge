#!/usr/bin/env sage-python
"""Fast coefficient-linear search for a Wilson quotient gauge."""

from sage.all import QQ, PolynomialRing, matrix

R = PolynomialRing(QQ, "n")
n = R.gen()

M = matrix(R, [
    [(-2*n-5)*(n+3)**2*(136*n**4+1424*n**3+5548*n**2+9551*n+6141),
     384*n**6+6384*n**5+44168*n**4+162698*n**3+336377*n**2+369933*n+169011,
     -(480*n**4+4980*n**3+19210*n**2+32690*n+20730)],
    [(n+2)**2*(n+3)**2*(4*n+10)*(48*n**3+386*n**2+1017*n+879),
     (n+2)**2*(-272*n**5-3848*n**4-21732*n**3-61184*n**2-85761*n-47808),
     (n+2)**2*(320*n**3+2540*n**2+6610*n+5640)],
    [(-4*n-10)*(n+2)**2*(n+3)**2*(32*n**4+302*n**3+1037*n**2+1530*n+813),
     (n+2)**2*(192*n**6+2984*n**5+19116*n**4+64452*n**3+120256*n**2+117279*n+46476),
     (n+2)**2*(-16*n**5-408*n**4-2912*n**3-8884*n**2-12254*n-6240)],
])
dA = -2*(n+2)**2*(n+3)**2*(2*n+5)*(2*n+7)**2
d0 = 8*(n+3)**2*(4*n+11)**2
dT = d0*(4*n+13)**2
Tnum = matrix(R, [
    [(4*n+13)**2*(4*n+9)**2*(40*n**2+228*n+325),
     (4*n+13)**2*(1536*n**4+16512*n**3+66496*n**2+118896*n+79641)],
    [(4*n+9)**2*(1536*n**4+18432*n**3+82880*n**2+165504*n+123841),
     59392*n**6+1012736*n**5+7184384*n**4+27140352*n**3+
     57583336*n**2+65059404*n+30580677],
])


def shifted(poly):
    return poly(n=n+1)


def residual_basis(row, column, power, denominator):
    P = matrix(R, 3, 2, 0)
    P[row, column] = n**power
    Pnext = P.apply_map(shifted)
    return dT*M*Pnext*denominator-dA*P*Tnum.transpose()*shifted(denominator)


def search(denominator, maximum_degree):
    denominator = R(denominator)
    print("denominator", denominator.factor(), flush=True)
    all_columns = []
    maximum_residual_degree = -1
    for degree in range(maximum_degree+1):
        for row in range(3):
            for column in range(2):
                residual = residual_basis(row, column, degree, denominator)
                all_columns.append(residual.list())
                maximum_residual_degree = max(
                    maximum_residual_degree,
                    max(value.degree() for value in residual.list()),
                )
        coefficient_matrix = matrix(
            QQ, 6*(maximum_residual_degree+1), len(all_columns),
            lambda equation, unknown: all_columns[unknown][equation //
                (maximum_residual_degree+1)][equation %
                (maximum_residual_degree+1)],
        )
        kernel = coefficient_matrix.right_kernel_matrix()
        print(" degree", degree, "shape", coefficient_matrix.dimensions(),
              "nullity", kernel.nrows(), flush=True)
        if kernel.nrows():
            for solution in kernel.rows():
                gauge = matrix(R, 3, 2, lambda row, column: sum(
                    solution[6*power+2*row+column]*n**power
                    for power in range(degree+1)
                ))
                print(" rank", gauge.rank(), flush=True)
                for gauge_row in gauge.rows():
                    print([value.factor()/denominator.factor()
                           for value in gauge_row], flush=True)
                if gauge.rank() == 2:
                    return gauge/denominator
    return None


candidates = [
    1,
    (n+1)*(n+2)*(2*n+3)**2*(2*n+5)**2,
    (n+1)*(n+2)*(2*n+3)**2*(2*n+5)**2*(4*n+9)**2,
    (n+1)**2*(n+2)**2*(2*n+3)**2*(2*n+5)**2*(4*n+9)**2,
    (n+1)**2*(n+2)**2*(2*n+3)**4*(2*n+5)**4*(4*n+9)**2,
]
for candidate in candidates:
    answer = search(candidate, 16)
    if answer is not None:
        print("FOUND", flush=True)
        raise SystemExit(0)
print("not found", flush=True)
