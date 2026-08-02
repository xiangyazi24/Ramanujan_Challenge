#!/usr/bin/env python3
"""Modular polynomial-intertwiner scan for the derived 3F2 P2.5 orbit."""

from fractions import Fraction as F

from p25_meijer_gauge import PRIME, inv, rat, challenge, matrix_rank, inverse_transpose


def derived(value):
    n = F(value)
    d0 = n**2 * (2*n-1)**2 * (4*n+3) * (4*n+5) * (16*n-12)
    d1 = n * (2*n-1)**2 * (4*n+3) * (4*n+5) * (16*n-12)
    matrix = [
        [
            (4*n+1)*(256*n**6-400*n**4+304*n**3+252*n**2-64*n+24)/d0,
            (4*n+1)*(8*n+4)*(32*n**3+24*n**2-2*n+3)/d0,
            (4*n+1)*(12*n+6)*(32*n**4+64*n**3+22*n**2+2*n+3)/d0,
        ],
        [
            96*(n+1)*(2*n+1)*(4*n-3)*(4*n-1)*(4*n+1)/d1,
            4*(4*n-3)*(4*n-1)*(4*n+1)*(68*n**3+140*n**2+81*n+12)/d1,
            24*(n+1)*(2*n+1)*(4*n-3)*(4*n-1)*(4*n+1)*(8*n**2+16*n+3)/d1,
        ],
        [
            (4*n+1)*(32*n+16)*(64*n**3+32*n**2-20*n+3)/d0,
            (4*n+1)**2*(16*n+8)*(48*n**3+36*n**2-20*n+3)/d0,
            (4*n+1)*(4352*n**6+12800*n**5+9328*n**4+208*n**3-756*n**2+24*n+36)/d0,
        ],
    ]
    return [[rat(entry) for entry in row] for row in matrix]


def transpose(matrix):
    return [[matrix[j][i] for j in range(3)] for i in range(3)]


def system(degree, shift, samples, variant):
    columns = 9 * (degree + 1)
    rows = []
    for value in samples:
        left = challenge(value)
        right = derived(F(value) + shift)
        if variant & 1:
            left = transpose(left)
            right = transpose(right)
        if variant & 2:
            right = inverse_transpose(right)
        if variant & 4:
            left, right = right, left
        powers = [pow(value % PRIME, r, PRIME) for r in range(degree + 1)]
        next_powers = [pow((value + 1) % PRIME, r, PRIME) for r in range(degree + 1)]
        for i in range(3):
            for j in range(3):
                row = [0] * columns
                for k in range(3):
                    for r in range(degree + 1):
                        # left(n) U(n+1) = U(n) right(n)
                        index = (k*3+j)*(degree+1)+r
                        row[index] = (row[index] + left[i][k]*next_powers[r]) % PRIME
                        index = (i*3+k)*(degree+1)+r
                        row[index] = (row[index] - powers[r]*right[k][j]) % PRIME
                rows.append(row)
    return rows, columns


def main():
    shifts = [F(value, 2) for value in range(-8, 17)]
    for variant in range(8):
        print("VARIANT", variant, flush=True)
        for shift in shifts:
            hit = None
            for degree in range(13):
                start = 20
                rows, columns = system(
                    degree, shift, range(start, start + degree + 8), variant
                )
                nullity = columns - matrix_rank(rows, columns)
                if nullity:
                    hit = (degree, nullity)
                    print("HIT", shift, degree, nullity, flush=True)
                    break
            if hit is None:
                print("none", shift, flush=True)


if __name__ == "__main__":
    main()
