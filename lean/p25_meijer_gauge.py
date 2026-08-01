#!/usr/bin/env python3
"""Temporary modular search for a Meijer-G coboundary to P2.5."""

from fractions import Fraction


PRIME = 1_000_000_007


def inv(value):
    return pow(value % PRIME, PRIME - 2, PRIME)


def rat(value):
    value = Fraction(value)
    return value.numerator % PRIME * inv(value.denominator) % PRIME


def challenge(n):
    n = Fraction(n)
    matrix = [
        [
            (-2*n-5)*(n+3)**2*(136*n**4+1424*n**3+5548*n**2+9551*n+6141),
            384*n**6+6384*n**5+44168*n**4+162698*n**3+336377*n**2+369933*n+169011,
            -(480*n**4+4980*n**3+19210*n**2+32690*n+20730),
        ],
        [
            (n+2)**2*(n+3)**2*(4*n+10)*(48*n**3+386*n**2+1017*n+879),
            (n+2)**2*(-272*n**5-3848*n**4-21732*n**3-61184*n**2-85761*n-47808),
            (n+2)**2*(320*n**3+2540*n**2+6610*n+5640),
        ],
        [
            (-4*n-10)*(n+2)**2*(n+3)**2*(32*n**4+302*n**3+1037*n**2+1530*n+813),
            (n+2)**2*(192*n**6+2984*n**5+19116*n**4+64452*n**3+120256*n**2+117279*n+46476),
            (n+2)**2*(-16*n**5-408*n**4-2912*n**3-8884*n**2-12254*n-6240),
        ],
    ]
    delta = -2*(n+2)**2*(n+3)**2*(2*n+5)*(2*n+7)**2
    return [[rat(value / delta) for value in row] for row in matrix]


def meijer(n):
    n = Fraction(n)
    matrix = [
        [
            (1024*n**6-4096*n**5-5440*n**4-1024*n**3+688*n**2+260*n+21) / (2*n*(2*n+1)*(4*n+1)**2*(4*n+3)**2),
            (8*n+3)*(2048*n**5+3520*n**4+1216*n**3-412*n**2-216*n-19) / (4*n*(2*n+1)*(4*n+1)**2*(4*n+3)**2),
            -(8*n+3)*(3072*n**6+8704*n**5+7584*n**4+1296*n**3-1006*n**2-365*n-29) / (4*n*(2*n+1)*(4*n+1)**2*(4*n+3)**2),
        ],
        [
            -2*(8*n+3)*(384*n**4+160*n**3-156*n**2-77*n-7) / (n*(2*n+1)*(4*n+1)**2*(4*n+3)**2),
            (8704*n**6+15104*n**5+3328*n**4-6256*n**3-4086*n**2-863*n-57) / (n*(2*n+1)*(4*n+1)**2*(4*n+3)**2),
            -(8*n+3)*(1536*n**6+3840*n**5+2080*n**4-1368*n**3-1522*n**2-395*n-29) / (n*(2*n+1)*(4*n+1)**2*(4*n+3)**2),
        ],
        [
            2*(8*n+3)*(256*n**3+280*n**2+86*n+7) / (n*(2*n+1)*(4*n+1)**2*(4*n+3)**2),
            -(8*n+3)*(768*n**4+1600*n**3+1064*n**2+258*n+19) / (n*(2*n+1)*(4*n+1)**2*(4*n+3)**2),
            (8704*n**6+31488*n**5+43264*n**4+28624*n**3+9538*n**2+1507*n+87) / (n*(2*n+1)*(4*n+1)**2*(4*n+3)**2),
        ],
    ]
    return [[rat(value) for value in row] for row in matrix]


def matrix_rank(rows, columns):
    rows = [row[:] for row in rows]
    rank = 0
    for column in range(columns):
        pivot = next((r for r in range(rank, len(rows)) if rows[r][column]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        scale = inv(rows[rank][column])
        rows[rank] = [value * scale % PRIME for value in rows[rank]]
        for r in range(rank + 1, len(rows)):
            if rows[r][column]:
                scale = rows[r][column]
                rows[r] = [
                    (rows[r][c] - scale * rows[rank][c]) % PRIME
                    for c in range(columns)
                ]
        rank += 1
        if rank == len(rows):
            break
    return rank


def inverse_transpose(matrix):
    augmented = [
        matrix[row][:] + [1 if row == column else 0 for column in range(3)]
        for row in range(3)
    ]
    for column in range(3):
        pivot = next(row for row in range(column, 3) if augmented[row][column])
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        scale = inv(augmented[column][column])
        augmented[column] = [value * scale % PRIME for value in augmented[column]]
        for row in range(3):
            if row != column and augmented[row][column]:
                scale = augmented[row][column]
                augmented[row] = [
                    (augmented[row][c] - scale * augmented[column][c]) % PRIME
                    for c in range(6)
                ]
    inverse = [row[3:] for row in augmented]
    return [[inverse[j][i] for j in range(3)] for i in range(3)]


def system(degree, shift, samples, dual):
    columns = 9 * (degree + 1)
    rows = []
    for n in samples:
        left = challenge(n)
        right = meijer(Fraction(n) + shift)
        if dual:
            right = inverse_transpose(right)
        powers = [pow(n % PRIME, r, PRIME) for r in range(degree + 1)]
        next_powers = [pow((n + 1) % PRIME, r, PRIME) for r in range(degree + 1)]
        for i in range(3):
            for j in range(3):
                row = [0] * columns
                for k in range(3):
                    for r in range(degree + 1):
                        # A(n) U(n+1)
                        index = (k * 3 + j) * (degree + 1) + r
                        row[index] = (row[index] + left[i][k] * next_powers[r]) % PRIME
                        # -U(n) T(n+s)
                        index = (i * 3 + k) * (degree + 1) + r
                        row[index] = (row[index] - powers[r] * right[k][j]) % PRIME
                rows.append(row)
    return rows, columns


def main():
    for dual in (False, True):
        print("DUAL", dual)
        for shift in [Fraction(value, 2) for value in range(-4, 13)]:
            for degree in range(0, 13):
                rows, columns = system(degree, shift, range(10, 10 + degree + 8), dual)
                nullity = columns - matrix_rank(rows, columns)
                if nullity:
                    print("HIT", "shift", shift, "degree", degree, "nullity", nullity)
                    break
            else:
                print("none shift", shift)


if __name__ == "__main__":
    main()
