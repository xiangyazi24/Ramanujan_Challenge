#!/usr/bin/env python3
"""Check the two CFVZ square-root branches and their rank-two recurrences."""

from math import comb


def primes_up_to(bound: int) -> list[int]:
    return [
        candidate
        for candidate in range(2, bound + 1)
        if all(candidate % divisor for divisor in range(2, int(candidate**0.5) + 1))
    ]


def apery(index: int) -> int:
    return sum(
        comb(index, k) ** 2 * comb(index + k, k) ** 2
        for k in range(index + 1)
    )


def multiply(left: list[int], right: list[int], prime: int) -> list[int]:
    product = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            product[i + j] = (product[i + j] + a * b) % prime
    while len(product) > 1 and product[-1] == 0:
        product.pop()
    return product


def divide_by_q(polynomial: list[int], prime: int) -> tuple[list[int], list[int]]:
    # q(t)=1-34t+t^2 is monic, so ordinary high-to-low division is exact.
    remainder = polynomial[:]
    quotient = [0] * (len(polynomial) - 2)
    for degree in range(len(polynomial) - 1, 1, -1):
        coefficient = remainder[degree] % prime
        quotient[degree - 2] = coefficient
        remainder[degree] = 0
        remainder[degree - 1] = (remainder[degree - 1] + 34 * coefficient) % prime
        remainder[degree - 2] = (remainder[degree - 2] - coefficient) % prime
    return quotient, remainder[:2]


def square_root(polynomial: list[int], degree: int, prime: int) -> list[int]:
    assert polynomial[0] % prime == 1
    root = [1]
    inverse_two = pow(2, -1, prime)
    for index in range(1, degree + 1):
        interior = sum(root[k] * root[index - k] for k in range(1, index))
        root.append((polynomial[index] - interior) * inverse_two % prime)
    return root


def check_rank_two_recurrence(root: list[int], epsilon: int, prime: int) -> None:
    for m in range(1, len(root) - 1):
        if epsilon == 0:
            residual = (
                4 * (m + 1) ** 2 * root[m + 1]
                - (136 * m * m + 68 * m + 10) * root[m]
                + (2 * m - 1) ** 2 * root[m - 1]
            )
        else:
            residual = (
                4 * (m + 1) ** 2 * root[m + 1]
                - (136 * m * m + 204 * m + 78) * root[m]
                + (2 * m + 1) ** 2 * root[m - 1]
            )
        assert residual % prime == 0


def main() -> None:
    checked = 0
    class_checks = set()
    for prime in primes_up_to(199):
        if prime < 5:
            continue
        ap = [apery(n) % prime for n in range(prime)]
        chi_minus_six = 1 if pow(-6 % prime, (prime - 1) // 2, prime) == 1 else -1
        epsilon = (1 - chi_minus_six) // 2
        class_checks.add((prime % 24, epsilon))
        if epsilon == 0:
            radicand = ap
            degree = (prime - 1) // 2
        else:
            radicand, remainder = divide_by_q(ap, prime)
            assert remainder == [0, 0]
            degree = (prime - 3) // 2
        root = square_root(radicand, degree, prime)
        reconstructed = multiply(root, root, prime)
        if epsilon:
            reconstructed = multiply([1, -34 % prime, 1], reconstructed, prime)
        assert reconstructed == ap
        check_rank_two_recurrence(root, epsilon, prime)
        checked += 1

    expected = {
        (1, 0), (5, 0), (7, 0), (11, 0),
        (13, 1), (17, 1), (19, 1), (23, 1),
    }
    assert class_checks == expected
    print(
        "VERIFIED CFVZ A_p=q^epsilon B_p^2 and the corresponding S_+/S_- "
        f"rank-two recurrence for {checked} primes 5<=p<=199"
    )
    print("VERIFIED epsilon=0 on p mod 24 in {1,5,7,11} and epsilon=1 on {13,17,19,23}")


if __name__ == "__main__":
    main()
