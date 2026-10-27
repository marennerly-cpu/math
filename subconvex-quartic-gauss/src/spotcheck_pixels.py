"""Spot-checks rendered figure 2: at several primes, the pixel colour in panels (a), quilt and (c) against the colour
expected from figdata.npz, and checks that the (c) phase equals the (a) phase turned by the quilt value."""
import os
import numpy as np
import matplotlib.image as mpimg
from palette import phase_rgb, UNIT_HEX
import matplotlib.colors as mc

HERE = os.path.dirname(os.path.abspath(__file__))
d = np.load(os.path.join(HERE, "figdata.npz"))
img = mpimg.imread(os.path.join(os.path.dirname(HERE), "figures", "new_fig2_cubic_index_trick.png"))[:, :, :3]
H, Wpx = img.shape[:2]
W = int(d["W"])
side, y0 = 0.272, 0.205
hgt = side * 16 / 8.4
xs = (0.03, 0.364, 0.698)


def pix(panel, x, y):
    fx = xs[panel] + (x + W + 0.5) / (2 * W + 1) * side
    fy = y0 + (y + W + 0.5) / (2 * W + 1) * hgt
    return int(round(fy * H)), int(round(fx * Wpx))


P2, a0, ak, quilt = d["P2"], d["a0"], d["ak"], d["quilt"]
lut = {1: UNIT_HEX[0], 1j: UNIT_HEX[1], -1: UNIT_HEX[2], -1j: UNIT_HEX[3]}
rng = np.random.default_rng(3)
worst = 0.0
for j in rng.choice(len(P2), 6, replace=False):
    x, y = int(P2[j, 0]), int(P2[j, 1])
    q = quilt[y + W, x + W]
    turn = np.angle(ak[j]) - np.angle(a0[j]) - np.angle(q)
    turn = (turn + np.pi) % (2 * np.pi) - np.pi
    r, c = pix(0, x, y)
    gotA = img[H - 1 - r, c]
    r, c = pix(2, x, y)
    gotC = img[H - 1 - r, c]
    r, c = pix(1, x, y)
    gotQ = img[H - 1 - r, c]
    expA, expC = phase_rgb(np.angle(a0[j])), phase_rgb(np.angle(ak[j]))
    expQ = np.array(mc.to_rgb(lut[complex(np.round(q.real), np.round(q.imag))]))
    errs = [np.abs(gotA - expA).max(), np.abs(gotC - expC).max(), np.abs(gotQ - expQ).max()]
    worst = max(worst, max(errs))
    print(f"c = {x:+d}{y:+d}i  quilt {q:.0f}  residual turn {turn:+.1e}  pixel error a/c/quilt:" + ", ".join(f"{e:.3f}" for e in errs))
print("worst pixel channel error:", round(float(worst), 3))
