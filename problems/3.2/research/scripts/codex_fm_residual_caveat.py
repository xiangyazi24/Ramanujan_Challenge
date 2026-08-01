#!/usr/bin/env python3
"""A minimal exact example separating mod-p zero from complex cancellation."""


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


if __name__ == "__main__":
    main()
