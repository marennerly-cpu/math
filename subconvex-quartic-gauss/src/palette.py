"""Cyclic phase colours on a dark surface.

Hue runs once round the OKLCH wheel (phase 0 = the value +1 = gold, i = teal, -1 = indigo, -i = crimson); lightness is
interpolated between four anchors so that the four units stay distinct under colour-vision deficiency
(validated with the dataviz skill's validate_palette.py, dark mode, all pairs: worst CVD dE 14.1, normal 16.6).
"""
import numpy as np

C = 0.12
H0 = 100.0
L_ANCHORS = np.array([0.66, 0.66, 0.54, 0.50])     # at phases 0, pi/2, pi, 3pi/2
BG = "#0b1020"
PANEL = "#111830"
INK = "#e9edf5"
MUTED = "#9aa4b8"
FAINT = "#2c3652"


def oklch_to_srgb(Lv, Cv, hdeg):
    h = np.deg2rad(hdeg)
    a, b = Cv * np.cos(h), Cv * np.sin(h)
    l_ = Lv + 0.3963377774 * a + 0.2158037573 * b
    m_ = Lv - 0.1055613458 * a - 0.0638541728 * b
    s_ = Lv - 0.0894841775 * a - 1.2914855480 * b
    l, m, s = l_ ** 3, m_ ** 3, s_ ** 3
    r = 4.0767416621 * l - 3.3077115913 * m + 0.2309699292 * s
    g = -1.2684380046 * l + 2.6097574011 * m - 0.3413193965 * s
    bb = -0.0041960863 * l - 0.7034186147 * m + 1.7076147010 * s
    lin = np.stack([r, g, bb], axis=-1)
    clipped = float(np.max(np.abs(lin - np.clip(lin, 0, 1))))
    lin = np.clip(lin, 0, 1)
    srgb = np.where(lin <= 0.0031308, 12.92 * lin, 1.055 * lin ** (1 / 2.4) - 0.055)
    return srgb, clipped


def lightness(theta):
    u = np.mod(np.asarray(theta, dtype=float), 2 * np.pi) / (np.pi / 2)   # 0..4
    j = np.floor(u).astype(int) % 4
    f = u - np.floor(u)
    return (1 - f) * L_ANCHORS[j] + f * L_ANCHORS[(j + 1) % 4]


def phase_rgb(theta):
    """theta in radians (array) -> sRGB rows."""
    theta = np.asarray(theta, dtype=float)
    rgb, _ = oklch_to_srgb(lightness(theta), C, H0 + np.rad2deg(theta))
    return rgb


def hexes(thetas):
    rgb = phase_rgb(np.asarray(thetas))
    return ["#%02x%02x%02x" % tuple(int(round(255 * v)) for v in row) for row in np.atleast_2d(rgb)]


UNIT_HEX = hexes(np.array([0, np.pi / 2, np.pi, 3 * np.pi / 2]))

if __name__ == "__main__":
    th = np.deg2rad(np.arange(0, 360, 1.0))
    _, clip = oklch_to_srgb(lightness(th), C, H0 + np.rad2deg(th))
    print("max gamut clip over the wheel:", round(clip, 4))
    print("units 1, i, -1, -i:", ",".join(UNIT_HEX))
