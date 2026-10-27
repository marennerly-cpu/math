# Figures for the lay-reader article

Three figures (`../figures/new_fig1_one_direction.png`, `../figures/new_fig2_cubic_index_trick.png`, `../figures/new_fig3_blueprint.png`), 200 dpi.

Requires Python 3 with numpy and matplotlib. The figures use the Cambria and Calibri fonts (installed with Windows/Office); without
them matplotlib substitutes another font and the text may not fit the layout exactly. `figdata.npz` is the output of
`compute_data.py`, included so the figures can be re-rendered without recomputing.

Rebuild, in order (each step must exit 0 before the next):

    python compute_data.py     # Gauss sums from the paper's definitions -> figdata.npz (~20 s)
    python verify_data.py      # independent checks; prints ALL CHECKS PASS, exits 1 on any failure (~65 s)
    python make_figures.py     # renders the PNGs
    python spotcheck_pixels.py # rendered pixel colours in figure 2 vs the data

Facts the figures state, and the check behind each (verify_data.py, run 2026-10-26):

- All 1,855 primary Gaussian primes with N(c) <= 16,000 (count reconciled with rational primes); the normalized quadratic
  sums g(-1/2, c)/sqrt(N c) all equal +1, the quartic ones all have modulus 1: checks [3], [4], and [7] (an independent
  plain-Python recomputation of every plotted sum, max error 8e-14).
- Quartic phases by quarter-turn: 506 / 464 / 427 / 458.
- g4(-k^3/2, c) = (c/k)_4 g4(-1/2, c) at all 187 primes in the window, k = -7+4i, N(k) = 65 = 1 mod 8: checks [5], [7].
- The quilt (c/k)_4 has period k and vanishes exactly where a prime of k divides c: check [6] (exhaustive).
- Figure 3 text: the abstract of Holowinsky-Nelson (arXiv 1801.08593) and the introduction of `../paper.pdf` (Context and
  Outline). It is a diagram of text, not computed, and has not been checked by an outside referee.

Colours: `palette.py`; the four units pass the dataviz validator on the dark surface (worst colour-blind dE 14.1;
crimson contrast 2.95:1, relieved by labels).
