"""Renders the three article figures from figdata.npz (figures 1, 2) and from text taken from the paper (figure 3)."""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import FancyBboxPatch, Polygon, Circle
from palette import BG, PANEL, INK, MUTED, FAINT, UNIT_HEX, phase_rgb, hexes

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "figures")
DPI = 200

for f in ("cambria.ttc", "cambriab.ttf", "calibri.ttf", "calibrib.ttf", "calibrii.ttf"):
    p = os.path.join(r"C:\Windows\Fonts", f)
    if os.path.exists(p):
        font_manager.fontManager.addfont(p)
plt.rcParams.update({
    "font.family": "Calibri", "mathtext.fontset": "stix", "text.color": INK,
    "axes.edgecolor": FAINT, "savefig.facecolor": BG, "figure.facecolor": BG,
})
SERIF = "Cambria"
GOLD, TEAL, INDIGO, CRIMSON = UNIT_HEX

d = np.load(os.path.join(HERE, "figdata.npz"))


def title(fig, main, sub, y=0.955):
    fig.text(0.03, y, main, fontsize=26, family=SERIF, weight="bold", color=INK, va="top")
    fig.text(0.03, y - 0.062, sub, fontsize=14.5, color=MUTED, va="top")


def plane_axes(ax, lim):
    ax.set_facecolor(BG)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.axhline(0, color=FAINT, lw=0.8, zorder=0)
    ax.axvline(0, color=FAINT, lw=0.8, zorder=0)
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


def glow_scatter(ax, x, y, rgb, core=7, halo=46):
    ax.scatter(x, y, s=halo, c=rgb, alpha=0.16, linewidths=0, zorder=2)
    ax.scatter(x, y, s=core, c=rgb, alpha=1.0, linewidths=0, zorder=3)


def wheel_key(ax, label_size=13):
    """A ring coloured by phase, labelled at the four units."""
    th = np.linspace(0, 2 * np.pi, 361)
    r = np.array([0.72, 1.0])
    T, R = np.meshgrid(th, r)
    cols = phase_rgb(0.5 * (th[:-1] + th[1:]))
    ax.pcolormesh(T, R, np.arange(360)[None, :], cmap=matplotlib.colors.ListedColormap(cols), shading="flat")
    ax.set_ylim(0, 1.0)
    ax.set_axis_off()
    for ang, lab in ((0, "1"), (90, "$i$"), (180, "$-1$"), (270, "$-i$")):
        a = np.deg2rad(ang)
        ax.text(a, 1.16, lab, ha="center", va="center", fontsize=label_size, color=INK, clip_on=False)
    ax.text(0, 0, "direction", ha="center", va="center", fontsize=11, color=MUTED)


# ---------------------------------------------------------------------------------------------------------------
def figure1():
    P1, g2, g4 = d["P1"], d["g2"], d["g4"]
    X1 = int(d["X1"])
    n = len(P1)
    lim = np.sqrt(X1) * 1.03
    fig = plt.figure(figsize=(15, 8.6), dpi=DPI)
    title(fig, "One direction, or every direction",
          "Each dot is a Gaussian prime $c$, placed at its own position in the complex plane and coloured by the direction\n"
          "of its normalized Gauss sum $g(-\\frac{1}{2},c)/\\sqrt{N(c)}$, which always has length exactly 1.")
    W1 = 0.27
    axL = fig.add_axes([0.03, 0.13, W1, W1 * 15 / 8.6])
    W2 = 0.42
    axR = fig.add_axes([0.335, 0.035, W2, W2 * 15 / 8.6 * 0.98])
    for ax in (axL, axR):
        plane_axes(ax, lim)
        ax.add_patch(Circle((0, 0), np.sqrt(X1), fill=False, ec=FAINT, lw=0.8, ls=(0, (2, 3))))
    glow_scatter(axL, P1[:, 0], P1[:, 1], phase_rgb(np.angle(g2)), core=3.2, halo=20)
    glow_scatter(axR, P1[:, 0], P1[:, 1], phase_rgb(np.angle(g4)))

    yL = 0.13 + W1 * 15 / 8.6
    fig.text(0.03, yL + 0.035, "Quadratic (Gauss's squares)", fontsize=17, family=SERIF, weight="bold", color=INK)
    fig.text(0.03, yL + 0.008, "$(x/c)_4^{\\,2}$ in place of $(x/c)_4$", fontsize=12.5, color=MUTED)
    fig.text(0.03, 0.105, f"All {n:,} sums equal exactly $+1$.\nThe arrows never turn.", fontsize=14, color=GOLD, va="top",
             weight="bold")
    yR = 0.035 + W2 * 15 / 8.6 * 0.98
    fig.text(0.355, yR - 0.005, "Quartic (fourth powers)", fontsize=17, family=SERIF, weight="bold", color=INK)
    fig.text(0.355, yR - 0.037, "$g_4(-\\frac{1}{2},c)/\\sqrt{N(c)}$", fontsize=12.5, color=MUTED)

    axW = fig.add_axes([0.805, 0.53, 0.15, 0.15 * 15 / 8.6], projection="polar")
    wheel_key(axW)
    fig.text(0.88, 0.475, "colour = direction of the arrow", ha="center", fontsize=12, color=MUTED)

    axH = fig.add_axes([0.805, 0.13, 0.15, 0.15 * 15 / 8.6], projection="polar")
    bins = 72
    h, e = np.histogram(np.mod(np.angle(g4), 2 * np.pi), bins=bins, range=(0, 2 * np.pi))
    ctr = 0.5 * (e[:-1] + e[1:])
    axH.bar(ctr, h, width=2 * np.pi / bins * 0.86, color=phase_rgb(ctr), bottom=0)
    axH.set_facecolor(BG)
    axH.set_xticks([])
    axH.set_yticks([])
    axH.spines["polar"].set_color(FAINT)
    axH.set_ylim(0, h.max() * 1.08)
    fig.text(0.88, 0.105, f"How the {n:,} quartic arrows point:\nevery direction, roughly equally",
             ha="center", va="top", fontsize=12, color=MUTED)

    fig.text(0.03, 0.018,
             f"All primary Gaussian primes $c$ with $N(c)\\leq${X1:,} ($c\\neq 1+i$); index $-\\frac{{1}}{{2}}$ is the base index of the paper. "
             "$g(r,c)=\\sum_{x\\,\\mathrm{mod}\\,c}\\,\\chi(x)\\,e^{2\\pi i\\,\\mathrm{Tr}(rx/c)}$ with $\\chi(x)=(x/c)_4$ (quartic) or $(x/c)_4^{\\,2}$ (quadratic), "
             "computed directly from the definition.",
             fontsize=10.5, color=MUTED)
    fig.savefig(os.path.join(OUT, "new_fig1_one_direction.png"), dpi=DPI)
    plt.close(fig)


# ---------------------------------------------------------------------------------------------------------------
def figure2():
    P2, a0, ak, quilt = d["P2"], d["a0"], d["ak"], d["quilt"]
    W = int(d["W"])
    k = tuple(int(v) for v in d["k"])
    k3 = tuple(int(v) for v in d["k3"])
    lim = W + 0.5
    fig = plt.figure(figsize=(16, 8.4), dpi=DPI)
    title(fig, "The cubic-index trick",
          "Moving the index from $-\\frac{1}{2}$ to $-k^3/2$ does not create new chaos: it recolours the old chaos by a perfectly periodic quilt.")
    side = 0.272
    hgt = side * 16 / 8.4
    y0 = 0.205
    xs = (0.03, 0.364, 0.698)
    axA = fig.add_axes([xs[0], y0, side, hgt])
    axQ = fig.add_axes([xs[1], y0, side, hgt])
    axC = fig.add_axes([xs[2], y0, side, hgt])
    for ax in (axA, axQ, axC):
        plane_axes(ax, lim)
    glow_scatter(axA, P2[:, 0], P2[:, 1], phase_rgb(np.angle(a0)), core=44, halo=190)
    glow_scatter(axC, P2[:, 0], P2[:, 1], phase_rgb(np.angle(ak)), core=44, halo=190)

    img = np.zeros(quilt.shape + (3,))
    lut = {1: GOLD, 1j: TEAL, -1: INDIGO, -1j: CRIMSON}
    for val, hx in lut.items():
        mask = np.abs(quilt - val) < 1e-9
        img[mask] = matplotlib.colors.to_rgb(hx)
    img[np.abs(quilt) < 1e-9] = matplotlib.colors.to_rgb("#151c33")
    axQ.imshow(img, origin="lower", extent=(-W - 0.5, W + 0.5, -W - 0.5, W + 0.5), interpolation="nearest", zorder=1)
    ik = (-k[1], k[0])
    # the period lattice k Z[i] (the multiples of k) as white points, one cell outlined
    jj, ll = np.meshgrid(np.arange(-10, 11), np.arange(-10, 11))
    lx = jj * k[0] + ll * ik[0]
    ly = jj * k[1] + ll * ik[1]
    inside = (np.abs(lx) <= W) & (np.abs(ly) <= W)
    axQ.plot(lx[inside], ly[inside], "o", ms=3.2, color=INK, zorder=3)
    cell = [(0, 0), k, (k[0] + ik[0], k[1] + ik[1]), ik]
    axQ.add_patch(Polygon(cell, closed=True, fill=False, ec=INK, lw=2.2, zorder=4))
    axQ.plot([0], [0], "o", ms=5, color=INK, zorder=5)
    axQ.set_xlim(-lim, lim)
    axQ.set_ylim(-lim, lim)

    top = y0 + hgt
    heads = (("Gauss-sum phases at index $-\\frac{1}{2}$", "$g_4(-\\frac{1}{2},c)/\\sqrt{N(c)}$ at each prime $c$"),
             ("The quilt $(c/k)_4$,  $k=%d+%di$" % k, "every Gaussian integer $c$; it repeats with period $k$ (white points)"),
             ("Phases at index $-k^3/2$", "$g_4(-k^3/2,c)/\\sqrt{N(c)}$, with $k^3=%d+%di$" % k3))
    for x, (h1, h2) in zip(xs, heads):
        fig.text(x, top + 0.052, h1, fontsize=16.5, family=SERIF, weight="bold", color=INK)
        fig.text(x, top + 0.018, h2, fontsize=12, color=MUTED)
    ymid = y0 + hgt / 2
    fig.text(0.3485, ymid, "×", fontsize=46, color=INK, ha="center", va="center")
    fig.text(0.6825, ymid, "=", fontsize=46, color=INK, ha="center", va="center")

    lx = 0.364
    for j, (hx, lab) in enumerate(zip(UNIT_HEX + ["#151c33"], ["$1$", "$i$", "$-1$", "$-i$", "0 (shares a prime with $k$)"])):
        xx = lx + 0.052 * j
        fig.patches.append(FancyBboxPatch((xx, y0 - 0.045), 0.014, 0.024, boxstyle="round,pad=0.002",
                                          transform=fig.transFigure, fc=hx, ec=FAINT, lw=0.6))
        fig.text(xx + 0.018, y0 - 0.033, lab, fontsize=12, color=INK, va="center")

    fig.text(0.03, 0.118,
             "$g_4(-k^3/2,\\,c)\\;=\\;(c/k)_4\\cdot g_4(-\\frac{1}{2},\\,c)$",
             fontsize=19, color=INK)
    fig.text(0.03, 0.072,
             f"Substituting $x\\to xk^{{-3}}$ in the Gauss sum pulls out $(k/c)_4$, and quartic reciprocity turns it into $(c/k)_4$ "
             f"(no sign here, since $N(k)={k[0]**2 + k[1]**2}$ leaves remainder 1 on division by 8).\n"
             f"Verified at all {len(P2)} primes shown by summing both sides directly. Bounding the series at $-k^3/2$ thus asks how far "
             "the irregular phases can fall into step with the quilt.",
             fontsize=12.5, color=MUTED, va="top", linespacing=1.5)
    fig.savefig(os.path.join(OUT, "new_fig2_cubic_index_trick.png"), dpi=DPI)
    plt.close(fig)


# ---------------------------------------------------------------------------------------------------------------
def figure3():
    fig = plt.figure(figsize=(16, 9.4), dpi=DPI)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 16)
    ax.set_ylim(0, 9.4)
    ax.axis("off")
    title(fig, "Same blueprint, new engine",
          "Holowinsky and Nelson's 2018 argument for $\\mathrm{GL}_3$, and the same argument run for quartic Gauss sums. "
          "Teal marks the parts that needed new work.")
    xs = np.array([1.45, 4.2, 6.95, 9.7, 12.45])
    yH, yP = 6.7, 4.05
    hn_col, new_col = "#7b869c", TEAL

    stage = ["The target", "Zero frequency\n+ amplifier", "Reciprocity", "Voronoi summation", "Cauchy–Schwarz\nand the endgame"]
    hn = ["$L(\\pi\\otimes\\chi,\\frac{1}{2})$: a fixed cusp\nform $\\pi$ on $\\mathrm{GL}_3$, twisted by\n$\\chi$ of prime conductor $q$",
          "write $\\chi$ as the zero-frequency\nterm of a Poisson-type formula;\namplify, using that $\\chi$ is\nmultiplicative",
          "reciprocity",
          "the $\\mathrm{GL}_3$ Voronoi formula",
          "Cauchy–Schwarz\nand the Weil bound"]
    paper = ["Kubota's series at $-k^3/2$,\nturned into sums twisted by\n$\\chi=(\\cdot/k)_4$, $k$ squarefree",
             "the same over $\\mathbb{Z}[i]$: amplifier\nprimes of norms near $K^{1/9}$ and\n$K^{5/18}$; Poisson over auxiliary\nmoduli $r$",
             "additive twists modulo $c=rs$:\nodd, but not necessarily\nsquarefree",
             "Theorem B: a Voronoi formula\nat every odd modulus, derived\nfrom Eckhardt–Patterson's\nfunctional equation",
             "Kloosterman correlations at\nprimes dividing $c$ once; prime-\npower losses absorbed by an\naverage over squarefull parts"]
    is_new = [True, False, False, True, True]

    xe = 14.75
    for y, col in ((yH, hn_col), (yP, new_col)):
        ax.plot([xs[0], xe], [y, y], color=col, lw=5, solid_capstyle="round", alpha=0.55, zorder=1)
        ax.add_patch(Circle((xe, y), 0.24, fc=col, ec="none", zorder=3))
    ax.text(xe, 7.62, "Result", ha="center", va="bottom", fontsize=15, family=SERIF, weight="bold", color=INK)
    for j, x in enumerate(xs):
        ax.text(x, 7.62, stage[j], ha="center", va="bottom", fontsize=15, family=SERIF, weight="bold", color=INK, linespacing=1.1)
        ax.add_patch(Circle((x, yH), 0.17, fc=BG, ec=hn_col, lw=3, zorder=3))
        ax.text(x, yH - 0.33, hn[j], ha="center", va="top", fontsize=12, color=MUTED, linespacing=1.35)
        ec = new_col if is_new[j] else "#c9d1df"
        if is_new[j]:
            ax.add_patch(Circle((x, yP), 0.36, fc=new_col, ec="none", alpha=0.18, zorder=2))
        ax.add_patch(Circle((x, yP), 0.17, fc=BG, ec=ec, lw=3, zorder=3))
        ax.text(x, yP - 0.33, paper[j], ha="center", va="top", fontsize=12, color=INK, linespacing=1.35)
        ax.plot([x, x], [yP + 0.25, yH - 0.25], color=FAINT, lw=1, ls=(0, (2, 3)), zorder=0)

    ax.text(0.22, yH + 0.3, "Holowinsky–Nelson (2018)", fontsize=13.5, color=hn_col, weight="bold")
    ax.text(0.22, yP + 0.3, "This paper", fontsize=13.5, color=new_col, weight="bold")
    ax.text(xe, yH - 0.36, "$q^{3/4-\\delta}$\nfor every\n$\\delta<1/36$", ha="center", va="top", fontsize=14, color=MUTED,
            linespacing=1.3)
    ax.text(xe, yP - 0.36, "$K^{3/4-1/36+\\varepsilon}$", ha="center", va="top", fontsize=15.5, color=INK)
    ax.text(xe, yP - 0.82, "Theorem A", ha="center", va="top", fontsize=12.5, color=new_col, weight="bold")
    ax.text(0.22, 2.32,
            "Also new: the approximate functional equation on the critical line (Appendix A), proved from the analytic continuation "
            "and functional equation due to Diaconu,\nin the form stated by David, Dunn, Hamieh and Lin.",
            fontsize=11.5, color=MUTED, va="top", linespacing=1.4)

    # the ruler: exponent of K, in 36ths
    y0, xa, xb = 0.42, 1.0, 15.0
    lo, hi = 24, 28
    X = lambda n36: xa + (n36 - lo) / (hi - lo) * (xb - xa)
    ax.plot([xa, xb], [y0, y0], color=MUTED, lw=1.2)
    for n36 in range(lo, hi + 1):
        ax.plot([X(n36), X(n36)], [y0, y0 + 0.16], color=MUTED, lw=1.2)
        ax.text(X(n36), y0 - 0.12, f"{n36}/36", ha="center", va="top", fontsize=11, color=MUTED)
    ax.add_patch(FancyBboxPatch((X(26), y0 + 0.02), X(27) - X(26), 0.42, boxstyle="square,pad=0", fc=new_col, alpha=0.22, ec="none"))
    for n36, lab, col in ((27, "$3/4$: convexity, the bound symmetry gives", "#c9d1df"),
                          (26, "$13/18$ $(+\\varepsilon)$: Theorem A", new_col)):
        ax.plot([X(n36), X(n36)], [y0, y0 + 0.62], color=col, lw=2.6)
        ax.text(X(n36) + (0.1 if n36 == 27 else -0.1), y0 + 0.45, lab, ha="left" if n36 == 27 else "right",
                va="center", fontsize=13, color=col, weight="bold")
    ax.text(0.5 * (X(26) + X(27)), y0 + 0.23, "1/36", ha="center", va="center", fontsize=13, color=INK, weight="bold")
    ax.text(0.22, y0 + 1.0, "Exponent of $K=N(k)$ in the bound for the series at $-k^3/2$, height $t$ fixed  (smaller is stronger):",
            fontsize=12.5, color=MUTED, va="bottom")
    fig.savefig(os.path.join(OUT, "new_fig3_blueprint.png"), dpi=DPI)
    plt.close(fig)


if __name__ == "__main__":
    figure1()
    figure2()
    figure3()
    print("written:", [f for f in sorted(os.listdir(OUT)) if f.startswith("new_fig")])
