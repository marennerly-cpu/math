"""Computes the data for figures 1 and 2 and writes figdata.npz next to this script."""
import os
import numpy as np
from gauss4 import gaussian_primes, gauss_sum, gmul, norm, quartic_symbol, symbol_composite

HERE = os.path.dirname(os.path.abspath(__file__))
X1 = 16000            # figure 1: all primary primes with N(c) <= X1
W = 30                # figure 2: window |Re c|, |Im c| <= W
K_FACTORS = [(-1, 2), (3, 2)]   # primary primes of norms 5 and 13
k = gmul(*K_FACTORS)            # k = -7 + 4i, N(k) = 65 = 1 mod 8
R0 = (-1, 0)                    # r0 = -1/2
k3 = gmul(gmul(k, k), k)
RK = (-k3[0], -k3[1])           # -k^3 / 2


def main():
    # figure 1
    P1 = gaussian_primes(X1)
    n1 = np.array([norm(c) for c in P1], dtype=float)
    g2 = np.array([gauss_sum(R0, c, quadratic=True) for c in P1]) / np.sqrt(n1)
    g4 = np.array([gauss_sum(R0, c) for c in P1]) / np.sqrt(n1)

    # figure 2: primes in the window, k's own prime factors excluded
    P2 = [c for c in gaussian_primes(2 * W * W)
          if abs(c[0]) <= W and abs(c[1]) <= W and c not in K_FACTORS]
    n2 = np.array([norm(c) for c in P2], dtype=float)
    a0 = np.array([gauss_sum(R0, c) for c in P2]) / np.sqrt(n2)
    ak = np.array([gauss_sum(RK, c) for c in P2]) / np.sqrt(n2)
    k_over_c = np.array([quartic_symbol(np.array([k[0]]), np.array([k[1]]), c)[0] for c in P2])
    c_over_k = np.array([symbol_composite(np.array([c[0]]), np.array([c[1]]), K_FACTORS)[0] for c in P2])

    s, t = np.meshgrid(np.arange(-W, W + 1), np.arange(-W, W + 1), indexing="xy")
    quilt = symbol_composite(s.ravel(), t.ravel(), K_FACTORS).reshape(s.shape)  # quilt[row=t+W, col=s+W]

    np.savez(os.path.join(HERE, "figdata.npz"),
             P1=np.array(P1), g2=g2, g4=g4, P2=np.array(P2), a0=a0, ak=ak,
             k_over_c=k_over_c, c_over_k=c_over_k, quilt=quilt, k=np.array(k), k3=np.array(k3),
             X1=X1, W=W)
    print(f"fig1 primes: {len(P1)}; fig2 primes: {len(P2)}; k = {k}, N(k) = {norm(k)}, k^3 = {k3}")


if __name__ == "__main__":
    main()
