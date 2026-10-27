"""Quartic and quadratic Gauss sums over Z[i] at primary Gaussian primes.

Definitions follow the paper (s2_notation.tex, s1_intro.tex):
  e(z) = exp(2 pi i Tr z), Tr z = z + conj(z);  c is primary if c = 1 mod lambda^3, lambda = 1 + i;
  g_4(r, c) = sum_{x mod c} (x/c)_4 e(r x / c)  for r in lambda^{-2} O and odd c.
An index r is passed as rnum with r = rnum / 2 (rnum a Gaussian integer), which covers lambda^{-2} O.
Gaussian integers are (real, imag) pairs of Python ints; residue arrays are int64 numpy arrays.
"""
import numpy as np

UNITS = [(1, 0), (0, 1), (-1, 0), (0, -1)]  # 1, i, -1, -i  ->  1j**k for k = 0..3


def gmul(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def norm(a):
    return a[0] * a[0] + a[1] * a[1]


def is_primary(c):
    # (c - 1) / lambda^3 = (c - 1)(-2 - 2i) / 8 must be a Gaussian integer
    x, y = c[0] - 1, c[1]
    return (-2 * x + 2 * y) % 8 == 0 and (-2 * x - 2 * y) % 8 == 0


def primary_associate(c):
    hits = [gmul(u, c) for u in UNITS if is_primary(gmul(u, c))]
    assert len(hits) == 1, c
    return hits[0]


def rational_primes(X):
    s = np.ones(X + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(X ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return [int(p) for p in np.nonzero(s)[0]]


def two_squares(p):
    for a in range(1, int(p ** 0.5) + 1):
        b2 = p - a * a
        b = int(round(b2 ** 0.5))
        if b * b == b2:
            return a, b
    raise ValueError(p)


def gaussian_primes(X):
    """Primary Gaussian primes c with N(c) <= X, lambda excluded. Split primes give two (conjugate) entries."""
    out = []
    for p in rational_primes(X):
        if p == 2:
            continue
        if p % 4 == 1:
            a, b = two_squares(p)
            out.append(primary_associate((a, b)))
            out.append(primary_associate((a, -b)))
        elif p * p <= X:
            out.append(primary_associate((p, 0)))
    return out


def residues(c):
    """A complete residue system mod the prime c: 0..p-1 if N(c) = p is prime (split), else s + t i, 0 <= s, t < q."""
    n = norm(c)
    if c[0] != 0 and c[1] != 0:  # split prime: Z/p -> O/c is a bijection
        xr = np.arange(n, dtype=np.int64)
        return xr, np.zeros_like(xr)
    q = abs(c[0]) + abs(c[1])  # inert prime: c = unit * q, N(c) = q^2
    assert q * q == n
    s, t = np.meshgrid(np.arange(q, dtype=np.int64), np.arange(q, dtype=np.int64), indexing="ij")
    return s.ravel(), t.ravel()


def reduce(zr, zi, c, n):
    """z mod c: subtract c times the Gaussian integer nearest to z / c = z conj(c) / n."""
    a, b = c
    nr = zr * a + zi * b
    ni = zi * a - zr * b
    qr = np.floor_divide(2 * nr + n, 2 * n)
    qi = np.floor_divide(2 * ni + n, 2 * n)
    return zr - (qr * a - qi * b), zi - (qr * b + qi * a)


def powmod(zr, zi, e, c, n):
    rr, ri = np.ones_like(zr), np.zeros_like(zr)
    br, bi = reduce(zr, zi, c, n)
    while e:
        if e & 1:
            rr, ri = reduce(rr * br - ri * bi, rr * bi + ri * br, c, n)
        e >>= 1
        if e:
            br, bi = reduce(br * br - bi * bi, 2 * br * bi, c, n)
    return rr, ri


def quartic_symbol(zr, zi, c):
    """(z/c)_4 for a primary prime c, by Euler's criterion z^((N c - 1)/4) = (z/c)_4 mod c."""
    zr = np.asarray(zr, dtype=np.int64)
    zi = np.asarray(zi, dtype=np.int64)
    n = norm(c)
    pr, pi_ = powmod(zr, zi, (n - 1) // 4, c, n)
    out = np.zeros(zr.shape, dtype=complex)
    matched = np.zeros(zr.shape, dtype=int)
    for k, (ur, ui) in enumerate(UNITS):
        dr, di = reduce(pr - ur, pi_ - ui, c, n)
        hit = (dr == 0) & (di == 0)
        out[hit] = (1, 1j, -1, -1j)[k]
        matched += hit
    r0, i0 = reduce(zr, zi, c, n)
    zero = (r0 == 0) & (i0 == 0)
    assert np.all(matched[~zero] == 1) and np.all(matched[zero] == 0), c
    return out


def symbol_composite(zr, zi, prime_factors):
    """(z/k)_4 for squarefree k = product of the given primary primes (multiplicative in the denominator)."""
    out = np.ones(np.shape(zr), dtype=complex)
    for p in prime_factors:
        out = out * quartic_symbol(zr, zi, p)
    return out


def gauss_sum(rnum, c, quadratic=False):
    """g(r, c) with r = rnum / 2; the quartic symbol, or its square (the quadratic symbol) if quadratic."""
    n = norm(c)
    xr, xi = residues(c)
    chi = quartic_symbol(xr, xi, c)
    if quadratic:
        chi = chi * chi
    Rr, Ri = rnum
    a, b = c
    # Tr(r x / c) = 2 Re(rnum x conj(c) / (2 n)) = Re(rnum x conj(c)) / n
    pr = Rr * xr - Ri * xi
    pim = Rr * xi + Ri * xr
    t = pr * a + pim * b
    phase = np.exp(2j * np.pi * (np.mod(t, n) / n))
    return complex(np.sum(chi * phase))
