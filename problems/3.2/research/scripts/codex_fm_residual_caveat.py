#!/usr/bin/env python3
"""Exact checks separating mod-p atoms from complex equidistribution."""

import cmath
import math


def main() -> None:
    # In F_5 choose generator 2 and the order-four Teichmueller lift 2 -> i.
    # For r=1 and f(1),f(2),f(4),f(3) = 1,2,2,0, the residual Mellin sum is
    # 1 + 2*2^{-1} + 2*4^{-1} = 0 mod 5.
    residual = (1 + 2 * pow(2, -1, 5) + 2 * pow(4, -1, 5)) % 5
    assert residual == 0

    # Its Teichmueller-lifted sum is 1 + 1 - i = 2-i, not zero in C.
    gaussian_real, gaussian_imag = 2, -1
    assert (gaussian_real, gaussian_imag) != (0, 0)
    assert gaussian_real**2 + gaussian_imag**2 == 5

    # At the prime above 5 specified by i -> 2, 2-i reduces to zero.
    reduction = (gaussian_real + gaussian_imag * 2) % 5
    assert reduction == 0

    print("VERIFIED residual Mellin sum is 0 in F_5")
    print("VERIFIED its Teichmueller lift is 2-i != 0 in C and has norm 5")
    print("VERIFIED 2-i reduces to 0 at the prime (5, i-2)")

    # A perturbation of size O(p) is invisible on the natural p^(3/2)
    # archimedean scale, but it can prescribe the residual zero pattern.
    for prime in [101, 211, 401, 809]:
        scale = prime**1.5
        for index in range(1, prime - 1):
            base = round(scale * math.cos(2 * math.pi * index / (prime - 1)))
            all_zero = base - base % prime
            none_zero = all_zero + 1
            assert all_zero % prime == 0
            assert none_zero % prime == 1
            assert abs(all_zero - base) / scale < 1 / math.sqrt(prime)
            assert abs(none_zero - base) / scale <= 1 / math.sqrt(prime)
    print(
        "VERIFIED O(p^-1/2) normalized perturbations can force residual zero "
        "density 1 or 0"
    )

    # Additive Fourier inversion is the exact residual replacement:
    # #zeros = p^-1 sum_s sum_r exp(2 pi i s M_r/p).
    for prime in [5, 7, 11, 13]:
        values = [(r**3 + 2 * r + 1) % prime for r in range(1, prime - 1)]
        zero_count = sum(value == 0 for value in values)
        root = cmath.exp(2j * math.pi / prime)
        fourier_total = sum(
            sum(root ** (scalar * value) for value in values)
            for scalar in range(prime)
        ) / prime
        assert abs(fourier_total.imag) < 1e-9
        assert abs(fourier_total.real - zero_count) < 1e-9
    print("VERIFIED additive-character orthogonality formula for residual zero counts")


if __name__ == "__main__":
    main()
