"""
This python module contains functions for plotting CoRAL output
"""
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.patches import RegularPolygon, Rectangle
from matplotlib.colors import LogNorm
import numpy as np


def _win_extent(pos, winsz):
    """Get (xmin, xmax, ymin, ymax) extent of a window for plotting.

    Mirrors the logic in corner_reflector_adp.get_win_bounds(), so the
    box drawn on the plot always matches the window actually used in
    the energy/clutter calculation.

    Parameters
    ----------
    pos : (x, y)
        Centre position of the window, in pixel coordinates.
    winsz : float/int, or 4-tuple/list (left, right, top, bottom)
        Window size. A single number gives a window symmetric around
        `pos`; a 4-tuple/list allows each side to have a different
        (adaptive) size.
    """
    if np.isscalar(winsz):
        left = right = top = bottom = winsz
    else:
        if len(winsz) != 4:
            raise ValueError(
                "winsz must be a single number (symmetric window) or a "
                "4-tuple/list of (left, right, top, bottom) sizes; got "
                f"{winsz!r}"
            )
        left, right, top, bottom = winsz

    xmin = pos[0] - left / 2
    xmax = pos[0] + right / 2
    ymin = pos[1] - top / 2
    ymax = pos[1] + bottom / 2

    return xmin, xmax, ymin, ymax


def plot_mean_intensity(avgI, cr_pos, targ_win_sz, clt_win_sz, name):
    '''Plot image of mean SAR intensity

    Parameters
    ----------
    avgI : ndarray
        Mean intensity image (dB).
    cr_pos : (x, y)
        Corner reflector centre position, in pixel coordinates.
    targ_win_sz, clt_win_sz : float/int, or 4-tuple/list (left, right, top, bottom)
        Target/clutter window size, matching whatever was passed to
        corner_reflector_adp.loop(). A single number draws a
        symmetric box around cr_pos (the original behaviour); a
        4-tuple/list draws an adaptive box whose sides can differ
        (e.g. shrunk on the side facing a nearby target/scene edge).
    name : str
        Site name, used in the plot title and output filename.
    '''
    # set black/white colormap for plots
    cmap = plt.set_cmap('gist_gray')

    # draw new plot
    #fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
    fig, ax1 = plt.subplots(1, 1, sharey=True)
    #ax = fig.add_subplot(1,1,1)
    cax = ax1.matshow(avgI, vmin=-20, vmax=10, cmap=cmap)
    #cax = ax2.matshow(avgI_d, vmin=-20, vmax=10, cmap=cmap)

    # define target window (small padding added purely so the red box
    # remains visible outside the target pixels themselves)
    tx_min, tx_max, ty_min, ty_max = _win_extent(cr_pos, targ_win_sz)
    p1 = Rectangle((tx_min - 1, ty_min - 1),
                   (tx_max - tx_min) + 2, (ty_max - ty_min) + 2,
                   linewidth=1, edgecolor='r', facecolor='none')

    # define clutter window (padding so the yellow box sits just
    # outside the target box, matching the old drawing convention)
    cx_min, cx_max, cy_min, cy_max = _win_extent(cr_pos, clt_win_sz)
    p2 = Rectangle((cx_min - 2, cy_min - 2),
                   (cx_max - cx_min) + 4, (cy_max - cy_min) + 4,
                   linewidth=1, edgecolor='y', facecolor='none')

    # add windows to plot
    ax1.add_patch(p1)
    ax1.add_patch(p2)

    # add text labels
    #ax1.text(45, 42, name, color='w', fontsize=10)

    # plot labels
    ax1.set_xlabel('Range')
    ax1.set_ylabel('Azimuth')

    # add colorbar
    cbar = fig.colorbar(cax)
    cbar.set_label('dB')

    # add title
    ax1.set_title('Mean intensity at %s' % name)

    # x-axis labels at bottom
    ax1.xaxis.set_tick_params(labeltop='False', labelbottom='True')

    # fit subplots and save fig
    fig.tight_layout()
    #fig.set_size_inches(w=6,h=4)

    # save PNG file
    fig.savefig('adp_mean_intensity_%s.png' % name, dpi=300, bbox_inches='tight')

    return


def plot_rcs_scr(t, rcs, scr, start, end, name):
    '''Plot RCS and SCR time series'''
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    plt.plot(t, rcs, 'ro-', label='RCS')
    plt.plot(t, scr, 'bo-', label='SCR')
    plt.xlim(start, end)
    plt.ylim(0, 40)
    plt.xlabel('Date')
    plt.ylabel('RCS / SCR (dB)')
    plt.legend(loc=4)
    plt.grid(True)
    plt.title('Corner Reflector response at %s' % name)
    for label in ax.get_xticklabels():
        label.set_rotation(90)
    fig.savefig('rcs_scr_%s.png' % name, dpi=300, bbox_inches='tight')

    return

def plot_clutter(t, clt, start, end, name):
    '''Plot average clutter time series'''
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    plt.plot(t, clt, 'bo-', label='Clutter')
    plt.xlim(start, end)
    plt.ylim(-16, -2)
    plt.xlabel('Date')
    plt.ylabel('Average Clutter (dB)')
    plt.legend(loc=1)
    plt.grid(True)
    plt.title('Average Clutter at %s' % name)
    for label in ax.get_xticklabels():
        label.set_rotation(90)
    fig.savefig('clutter_%s.png' % name, dpi=300, bbox_inches='tight')

    return
