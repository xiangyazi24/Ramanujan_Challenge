import numpy as np
from sympy import primitive_root


def gauss_table(p: int):
    """Return primitive root g and G[m] = sum_x omega^m(x) exp(2*pi*i*x/p)."""
    if p < 3:
        raise ValueError("p must be an odd prime")

    N = p - 1
    g = int(primitive_root(p))

    phase = np.empty(N, dtype=np.complex128)
    x = 1
    for r in range(N):
        phase[r] = np.exp(2j * np.pi * x / p)
        x = (x * g) % p

    # numpy.ifft has the positive Fourier sign and includes 1/N.
    G = N * np.fft.ifft(phase)
    return g, G


def apery_trace_one(p: int, z: int, G=None):
    """Compute the complex hypergeometric lift T_{p,z}, 1 <= z <= p-2."""
    N = p - 1
    if not (1 <= z <= p - 2):
        raise ValueError("The nondegenerate formula requires 1 <= z <= p-2")

    if G is None:
        _, G = gauss_table(p)

    j = np.arange(N)
    numerator = (
        G[(j + z) % N] ** 2
        * G[(j - z) % N] ** 2
        * G[(-j) % N] ** 4
    )
    denominator = G[z] ** 2 * G[(-z) % N] ** 2
    return -np.mean(numerator / denominator)


def apery_trace_all(p: int):
    """Compute T_{p,z} for all nondegenerate z. Cost after FFT: O(p^2)."""
    _, G = gauss_table(p)
    N = p - 1
    j = np.arange(N)
    Gminusj4 = G[(-j) % N] ** 4

    T = np.full(p, np.nan + 0j, dtype=np.complex128)
    for z in range(1, N):
        numerator = (
            G[(j + z) % N] ** 2
            * G[(j - z) % N] ** 2
            * Gminusj4
        )
        denominator = G[z] ** 2 * G[(-z) % N] ** 2
        T[z] = -np.mean(numerator / denominator)

    return T

p = 1009
T = apery_trace_all(p)
z = np.arange(1, p - 1)
normalized = T[z].real / (p ** 1.5)

print("max imaginary residual:", np.max(np.abs(T[z].imag)))
print("range:", normalized.min(), normalized.max())
print("central trace:", T[(p - 1) // 2])
