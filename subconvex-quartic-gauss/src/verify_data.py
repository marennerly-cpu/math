"""Independent checks of figdata.npz and of the symbol code; prints every figure that a caption will quote."""
import os
import random
import sys
import numpy as np
from gauss4 import gaussian_primes, norm, quartic_symbol, is_primary, gmul, rational_primes
from compute_data import K_FACTORS

HERE = os.path.dirname(os.path.abspath(__file__))
d = np.load(os.path.join(HERE, "figdata.npz"))


def ref_symbol(z, c):
    """(z/c)_4 for a Gaussian integer z, by plain Python arithmetic independent of gauss4's reduce/powmod.
    Split c = a+bi: map z = s+ti to s + t*u in F_p with u = -a/b = i mod c.  Inert c = unit*q: O/q has
    componentwise reduction mod q, and Euler's criterion is evaluated with tuple multiplication mod q."""
    a, b = int(c[0]), int(c[1])
    n = a * a + b * b
    s, t = int(z[0]), int(z[1])
    if a != 0 and b != 0:
        p = n
        u = (-a * pow(b, -1, p)) % p
        val = pow((s + t * u) % p, (p - 1) // 4, p)
        return {0: 0, 1: 1, u: 1j, p - 1: -1, p - u: -1j}[val]
    q = abs(a) + abs(b)
    base, acc, e = (s % q, t % q), (1, 0), (n - 1) // 4
    while e:
        if e & 1:
            acc = ((acc[0] * base[0] - acc[1] * base[1]) % q, (acc[0] * base[1] + acc[1] * base[0]) % q)
        base = ((base[0] * base[0] - base[1] * base[1]) % q, (2 * base[0] * base[1]) % q)
        e >>= 1
    return {(0, 0): 0, (1, 0): 1, (0, 1): 1j, (q - 1, 0): -1, (0, q - 1): -1j}[acc]


def ref_gauss(r, c):
    """(quartic g(r,c), quadratic g(r,c)) by direct summation in complex floating point: Tr(w) = 2 Re(w), w = r x / c,
    reduced mod 1 before exponentiating."""
    n = norm(c)
    if c[0] != 0 and c[1] != 0:
        res = [(x, 0) for x in range(n)]
    else:
        q = abs(c[0]) + abs(c[1])
        res = [(s, t) for s in range(q) for t in range(q)]
    cc = complex(*c)
    t4 = t2 = 0j
    for x in res:
        chi = ref_symbol(x, c)
        ph = np.exp(2j * np.pi * np.mod(2 * (r * complex(*x) / cc).real, 1.0))
        t4 += chi * ph
        t2 += chi * chi * ph
    return t4, t2


# 1. quartic symbol against the independent reference, on random complex residues of 40 split and all inert primes
random.seed(1)
allp = gaussian_primes(16000)
split = [c for c in allp if c[0] != 0 and c[1] != 0]
inert = [c for c in allp if c[0] == 0 or c[1] == 0]
bad = tested = 0
for c in random.sample(split, 40) + inert:
    zs = [(random.randint(-300, 300), random.randint(-300, 300)) for _ in range(300)]
    ours = quartic_symbol(np.array([z[0] for z in zs]), np.array([z[1] for z in zs]), c)
    ref = np.array([ref_symbol(z, c) for z in zs])
    bad += int(np.sum(np.abs(ours - ref) > 1e-12))
    tested += len(zs)
print(f"[1] symbol vs independent reference, 40 split + {len(inert)} inert primes, {tested} complex arguments: mismatches = {bad}")

# 2. quartic reciprocity for distinct primary primes: (p/q)_4 = (q/p)_4 (-1)^(((Np-1)/4)((Nq-1)/4))
P = gaussian_primes(3000)
viol = 0
for _ in range(2000):
    p, q = random.sample(P, 2)
    pq = quartic_symbol(np.array([p[0]]), np.array([p[1]]), q)[0]
    qp = quartic_symbol(np.array([q[0]]), np.array([q[1]]), p)[0]
    sign = (-1) ** (((norm(p) - 1) // 4) * ((norm(q) - 1) // 4))
    viol += int(abs(pq - sign * qp) > 1e-12)
print(f"[2] quartic reciprocity on 2000 random pairs: violations = {viol}")

# 3. all listed primes primary, distinct, and prime count reconciles with rational primes
P1 = [tuple(int(v) for v in c) for c in d["P1"]]
X1 = int(d["X1"])
rp = rational_primes(X1)
expected = 2 * sum(1 for p in rp if p % 4 == 1) + sum(1 for p in rp if p % 4 == 3 and p * p <= X1)
print(f"[3] fig1 primes: {len(P1)} (expected from rational primes: {expected}), all primary: {all(is_primary(c) for c in P1)}, "
      f"distinct: {len(set(P1)) == len(P1)}, max N = {max(norm(c) for c in P1)}")

# 4. |g| = sqrt(N) (normalised modulus 1), and the quadratic sums are +-1
g2, g4 = d["g2"], d["g4"]
print(f"[4] max | |g4/sqrtN| - 1 | = {np.max(np.abs(np.abs(g4) - 1)):.2e};  max | |g2/sqrtN| - 1 | = {np.max(np.abs(np.abs(g2) - 1)):.2e}")
dist_pm1 = np.minimum(np.abs(g2 - 1), np.abs(g2 + 1))
print(f"    max distance of g2/sqrtN from {{+1,-1}} = {dist_pm1.max():.2e};  count +1: {int(np.sum(np.abs(g2 - 1) < 1e-6))}, count -1: {int(np.sum(np.abs(g2 + 1) < 1e-6))}")
ang = np.angle(g4)
quad = np.histogram(np.mod(ang, 2 * np.pi), bins=4, range=(0, 2 * np.pi))[0]
print(f"    quartic phases by quarter-turn: {quad.tolist()};  |mean of g4/sqrtN| = {abs(g4.mean()):.4f}")
dist4 = np.min(np.abs(g4[:, None] - np.array([1, 1j, -1, -1j])[None, :]), axis=1)
print(f"    quartic phases within 0.01 of a unit: {int(np.sum(dist4 < 0.01))} of {len(g4)}")

# 5. figure 2 identity: g4(-k^3/2, c) = (k/c)_4 g4(-1/2, c), and (k/c)_4 = (c/k)_4 since N(k) = 1 mod 8
a0, ak, kc, ck = d["a0"], d["ak"], d["k_over_c"], d["c_over_k"]
print(f"[5] fig2 primes: {len(a0)};  max |ak - (k/c)a0| = {np.max(np.abs(ak - kc * a0)):.2e};  (k/c)_4 == (c/k)_4 for all: {bool(np.all(np.abs(kc - ck) < 1e-12))}")
k = tuple(int(v) for v in d["k"])
print(f"    k = {k}, N(k) = {norm(k)}, N(k) mod 8 = {norm(k) % 8}, k^3 = {tuple(int(v) for v in d['k3'])}, gmul check = {gmul(gmul(k, k), k)}")

# 6. quilt: period k (and ik), zeros exactly at multiples of a prime factor of k, value counts
q = d["quilt"]
W = int(d["W"])
per_bad = 0
for s in range(-W, W + 1):
    for t in range(-W, W + 1):
        for (dr, di) in (k, (-k[1], k[0])):
            s2, t2 = s + dr, t + di
            if -W <= s2 <= W and -W <= t2 <= W:
                per_bad += int(abs(q[t + W, s + W] - q[t2 + W, s2 + W]) > 1e-12)
vals, counts = np.unique(np.round(q, 6), return_counts=True)
zero_bad = 0
for t in range(-W, W + 1):
    for s in range(-W, W + 1):
        divisible = any(((s * f[0] + t * f[1]) % norm(f) == 0 and (t * f[0] - s * f[1]) % norm(f) == 0) for f in K_FACTORS)
        zero_bad += int((abs(q[t + W, s + W]) < 1e-12) != divisible)
print(f"[6] quilt period violations (shifts by k, ik): {per_bad};  zero exactly where a prime of k divides c: violations = {zero_bad};  "
      f"values: {dict(zip([str(v) for v in vals], counts.tolist()))}")

# 7. Gauss sums recomputed by the independent reference (floating-point trace, plain-Python symbol)
P2 = [tuple(int(v) for v in c) for c in d["P2"]]
k3 = complex(*(int(v) for v in d["k3"]))
e2 = e4 = 0.0
plus1 = 0
for j, c in enumerate(P1):
    r4, r2 = ref_gauss(-0.5, c)
    s = np.sqrt(norm(c))
    e4 = max(e4, abs(r4 / s - g4[j]))
    e2 = max(e2, abs(r2 / s - g2[j]))
    plus1 += int(abs(r2 / s - 1) < 1e-9)
ea0 = eak = 0.0
for j, c in enumerate(P2):
    s = np.sqrt(norm(c))
    ea0 = max(ea0, abs(ref_gauss(-0.5, c)[0] / s - a0[j]))
    eak = max(eak, abs(ref_gauss(-k3 / 2, c)[0] / s - ak[j]))
print(f"[7] reference recomputation, all {len(P1)} fig-1 primes: max err g2 {e2:.1e}, g4 {e4:.1e}, reference g2 = +1 at {plus1};  "
      f"all {len(P2)} fig-2 primes: max err a0 {ea0:.1e}, ak {eak:.1e}")

failures = [name for name, ok in (
    ("symbol", bad == 0), ("reciprocity", viol == 0), ("prime count", len(P1) == expected),
    ("primary/distinct", all(is_primary(c) for c in P1) and len(set(P1)) == len(P1)),
    ("modulus", np.max(np.abs(np.abs(g4) - 1)) < 1e-9 and np.max(np.abs(np.abs(g2) - 1)) < 1e-9),
    ("identity", np.max(np.abs(ak - kc * a0)) < 1e-9 and bool(np.all(np.abs(kc - ck) < 1e-12))),
    ("quilt", per_bad == 0 and zero_bad == 0),
    ("reference", max(e2, e4, ea0, eak) < 1e-8 and plus1 == len(P1)))
    if not ok]
print("ALL CHECKS PASS" if not failures else f"FAILED: {failures}")
sys.exit(1 if failures else 0)
