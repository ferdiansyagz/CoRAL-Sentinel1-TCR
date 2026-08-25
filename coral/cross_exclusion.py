"""
Extension module for CoRAL implementing the cross-exclusion clutter
window definition from Garthwaite et al. (2015), Geoscience Australia
Record 2015/03, Figure 4.11 / Table 4.5.

This module does NOT modify coral/corner_reflector.py or coral/dataio.py.
It re-uses calc_integrated_energy, calc_total_energy, calc_scr, calc_rcs
and get_win_bounds from corner_reflector.py unchanged, and only replaces
the clutter-region definition (originally calc_clutter_intensity, which
computes clutter as a simple square annulus: clutter window minus
target window, both centred on the CR).

Cross-exclusion instead removes a cross/plus-shaped region of width
`cross_width` pixels, extending the full length of the clutter window
in both the range and azimuth directions, before summing the remaining
four quadrants as clutter. This matches Figure 4.11 of Garthwaite et
al. (2015).

Matches the signature of the modified corner_reflector.py / dataio.py
in this project, which support a `band=` argument for SNAP BEAM-DIMAP
(.dim) inputs (e.g. Sentinel-1 IW products calibrated in SNAP).
"""
import re
import numpy as np
from datetime import datetime

from coral.dataio import readfile
from coral.corner_reflector import (
    get_win_bounds,
    calc_integrated_energy,
    calc_total_energy,
    calc_scr,
    calc_rcs,
)


def cross_exclusion_mask(n_rows, n_cols, cross_width):
    """Boolean mask (True = clutter pixel) for a window of shape
    (n_rows, n_cols) centred on the target, after excluding a
    cross/plus-shaped region of width `cross_width` pixels that spans
    the full window in both the horizontal and vertical directions.

    n_rows/n_cols should be the ACTUAL cropped window size returned by
    get_win_bounds() -- note that for an even winsz, get_win_bounds()
    (as implemented in corner_reflector.py) returns a window one pixel
    larger than winsz (e.g. winsz=24 -> a 25x25 crop), so this must not
    assume n_rows == n_cols == clt_win_sz.

    Reproduces the cross exclusion shown in Figure 4.11 of Garthwaite
    et al. (2015) -- the remaining True pixels are the four corner
    quadrants of the clutter window.
    """
    center_y = (n_rows - 1) / 2.0
    center_x = (n_cols - 1) / 2.0
    half_cross = cross_width / 2.0

    yy, xx = np.mgrid[0:n_rows, 0:n_cols]

    in_cross = (np.abs(yy - center_y) <= half_cross) | (np.abs(xx - center_x) <= half_cross)

    return ~in_cross


def calc_clutter_intensity_cross(d, cr_pos, clt_win_sz, cross_width):
    """Calculate clutter energy/count using cross-exclusion instead of
    the simple square-annulus approach in
    corner_reflector.calc_clutter_intensity.

    Parameters
    ----------
    d : ndarray, shape (n_scenes, H, W)
        Stack of image windows (same as passed into corner_reflector.loop).
    cr_pos : array-like, length 2
        Pixel position [x, y] of the CR within each window (same
        convention as corner_reflector.loop, typically [sub_im, sub_im]).
    clt_win_sz : int
        Size (pixels) of the square clutter window, centred on cr_pos.
    cross_width : int
        Width (pixels) of the cross/plus-shaped exclusion region.

    Returns
    -------
    Avg_clt : list of float
        Average clutter intensity in dB, per scene.
    Eclt : list of float
        Total integrated clutter energy (sum over the 4 quadrants), per scene.
    Nclt : list of int
        Number of clutter pixels (4 quadrants), per scene.
    """
    xmin, xmax, ymin, ymax = get_win_bounds(cr_pos, clt_win_sz)
    n_rows = ymax - ymin
    n_cols = xmax - xmin
    mask = cross_exclusion_mask(n_rows, n_cols, cross_width)

    Avg_clt = []
    Eclt = []
    Nclt = []

    for i in range(d.shape[0]):
        subd = d[i, ymin:ymax, xmin:xmax]

        clutter_pixels = subd[mask]
        E = float(clutter_pixels.sum())
        N = int(clutter_pixels.size)

        Avg_clt.append(10 * np.log10(E / N))
        Eclt.append(E)
        Nclt.append(N)

    return Avg_clt, Eclt, Nclt


def loop_cross_exclusion(files, sub_im, cr, targ_win_sz, clt_win_sz, cross_width, band=None):
    """Same as corner_reflector.loop(), except the clutter region is
    computed via cross-exclusion (calc_clutter_intensity_cross) instead
    of the simple square-annulus method (calc_clutter_intensity).

    All other steps -- target energy, total energy correction (Garthwaite
    2017 Eq. 7), SCR, and RCS (Eq. 2 & 8) -- are unchanged and reused
    directly from corner_reflector.py.

    band : str, optional
        Same as corner_reflector.loop() -- passed through to readfile()
        for SNAP BEAM-DIMAP (.dim) inputs, e.g. band='Sigma0_VV'.
    """
    d = np.empty((len(files), sub_im * 2, sub_im * 2))
    t = []

    for i, g in enumerate(files):
        m = re.search(r'\d{8}', g)
        if m:
            t.append(datetime.strptime(m.group(0), "%Y%m%d"))

        d[i], rho_r, rho_a, theta = readfile(g, sub_im, cr, band=band)

    avgI = 10 * np.log10(np.mean(d, axis=0))

    cr_pos = np.array([sub_im, sub_im])

    # Target energy: unchanged, reused from corner_reflector.py
    En, Ncr = calc_integrated_energy(d, cr_pos, targ_win_sz)

    # Clutter energy: cross-exclusion instead of square annulus
    Avg_clt, Eclt, Nclt = calc_clutter_intensity_cross(d, cr_pos, clt_win_sz, cross_width)

    # Total target energy correction, SCR, RCS: unchanged, reused from corner_reflector.py
    Ecr = calc_total_energy(Ncr, Nclt, Eclt, En)
    scr = calc_scr(Ecr, Eclt, Nclt)
    rcs = calc_rcs(Ecr, rho_r, rho_a, theta)

    return avgI, rcs, scr, Avg_clt, t
