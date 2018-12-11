import numpy as np
import matplotlib.pyplot as plt

from scipy import stats
from scipy.optimize import curve_fit
import uncertainties as un
from uncertainties.unumpy import nominal_values as nom
from uncertainties.unumpy import std_devs as err

from cbsyst import Csys

from .helpers import isolate_constant_conditions

def spreadm(x, y, x_tol, y_tol, x_offset=None, y_offset=None, offset_mult=0.2):
    """
    spreadm reditributes a set of overlapping x/y points around their mean for display.

    Parameters
    ----------
    x, y : array-like
        The x and y arrays containing overlapping points
    x_tol, y_tol : float
        The overlap tolerance. Points are redistributed if any
        of them are closer to the mean than either tolerance.
    x_offset, y_offset : float
        Absolute x/y offsets to redistribute the points by.
        If None, offsets are calculated based on the number
        of points, and the sizes of the tolerances, following:
        offset = tolerance * offset_mult * n
        Defaults to None.
    offset_mult : float
        Used to automatically calculate displacement offsets. See above.

    Returns
    -------
    x_new, y_new, x_mean, y_mean

    If no points are outside the tolerances, or there is only one
    point, x_mean and y_mean are None.
    """
    x_mean = np.nanmean(x)
    y_mean = np.nanmean(y)

    x_diff = abs(x - x_mean)
    y_diff = abs(y - y_mean)

    if any(x_diff > x_tol) or any(y_diff > y_tol) or (len(x) >= 2):
        n = len(x)
        rad = 2 * np.pi / n

        rads = np.arange(n) * rad

        if x_offset is None:
            x_offset = x_tol * offset_mult * n
        if y_offset is None:
            y_offset = y_tol * offset_mult * n

        x_out = []
        y_out = []
        for r in rads:
            y_off = y_offset * np.cos(r)
            x_off = x_offset * np.sin(r)

            x_out.append(x_off)
            y_out.append(y_off)
        return x_mean + x_out, y_mean + y_out, x_mean, y_mean
    else:
        return x, y, None, None

def compare(obs, pred, errs=None, c=None, clab=None, cmap=None, figax=None, show_rmsd=False, **kwargs):
    
    if figax is None:
        fig, ax = plt.subplots(1,1)
    else:
        fig, ax = figax

    m = ax.scatter(obs, pred, c=c, cmap=cmap, **kwargs)
    if errs is not None:
        ax.errorbar(obs, pred,
                    xerr=errs,
                    lw=0, elinewidth=2, capsize=0, color=(.7, .7, .7), zorder=-1)

    ax.set_xlabel('Measured')
    ax.set_ylabel('Predicted')
    if c is not None:
        if figax is None:
            fig.colorbar(m, label=clab)
        else:
            fig.colorbar(m, label=clab, ax=ax)
    
    lims = np.concatenate([ax.get_xlim(), ax.get_ylim()])
    lims = np.array([lims.min(), lims.max()])
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.plot(lims, lims, color='k', ls='dashed', alpha=0.7, zorder=-2)

    # if errs is not None:
    #     err = np.percentile(errs, 50)
    #     ax.fill_between(lims, lims+err, lims-err, color=(0,0,0,0.2), lw=0, zorder=-5)

    if show_rmsd:
        rmsd = np.sqrt(np.sum((obs - pred)**2) / len(pred))
        ax.text(.05,.95,'RMSD: {:.2f}'.format(rmsd),
                ha='left', va='top', transform=ax.transAxes)

    if figax is None:
        return fig, ax

def lbl_position_checker(lbls, ax, pad=0.05):
    """
    Function for shifting contour labels inside minimum axis limits.
    """
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    xr = xlim[1] - xlim[0]
    yr = ylim[1] - ylim[0]
    x_min = xlim[0] + pad * xr
    y_min = ylim[0] + pad * yr

    for lb in lbls:
        lx, ly = lb.get_position()
        if lx < x_min:
            lx = x_min
        if ly < y_min:
            ly = y_min
        lb.set_position((lx, ly))


def panel_coord(n, npanels=3, limits=(.1, 0, .9, 1), pad=(.01,.05)):
    """
    Function for making sub-panels for gradient plos

    """
    
    bwidth = (limits[2] - limits[0]) / npanels
    bleft = limits[0] + (n - 1) * bwidth
    box = (bleft + pad[0], pad[1], bwidth - 2 * pad[0], 1 - 2 * pad[1])
    
    top_axis = (box[0], .56, box[2], .4)
    bottom_axis = (box[0], .14, box[2], .4)
    
    # cbwidth = 0.45
    # colorbar_left = (box[0], .09, bwidth * cbwidth, .03)
    # colorbar_right = (box[0] + box[2] - bwidth * cbwidth, .09, bwidth * cbwidth, .03)
    
    return top_axis, bottom_axis #, colorbar_left, colorbar_right

def contour_labels(cs, ax, hshift=0, vshift=0, **kwargs):
    """
    From here: https://goo.gl/XKetmJ
    """

    mid = (np.mean(ax.get_xlim()), np.mean(ax.get_ylim()))
    xrng = abs(np.diff(ax.get_xlim()))[0] * 0.5
    yrng = abs(np.diff(ax.get_ylim()))[0] * 0.5

    target = np.array((mid[0] + xrng * hshift, mid[1] + yrng * vshift))

    label_pos = []
    for line in cs.collections:
        for path in line.get_paths():
            vert = path.vertices

            # find closest point to center
            dist = np.linalg.norm(vert - target, ord=2, axis=1)
            min_ind = np.argmin(dist)
            label_pos.append(vert[min_ind, :])
    
    ax.clabel(cs, inline=True, inline_spacing=3, rightside_up=True, manual=label_pos, **kwargs)

def gradpanel(x, y, Temp, mTemp, Tsens, fig, cm, cmg, Tlim, Tgradlim, n, npanes=3, contours=[15,20,25,30]):
    # mgax, gax = (fig.add_axes(p) for p in panel_coord(n, npanes, (.07,0,.9,1), pad=(.015,.01)))
    
    mgax = fig.add_axes(panel_coord(n, npanes, (.07,0,.9,1), pad=(.015,.01))[0])

    # mgax.set_xticklabels([])    

    # temperature contours
    

    # mg/ca plot
    vmin, vmax = Tlim
    cma = mgax.pcolormesh(x, y, mTemp, vmin=vmin, vmax=vmax, cmap=cm, zorder=-1)
    csl = mgax.contour(x, y, Temp, contours, colors=[(0,0,0,0.3)], linewidths=1)
    contour_labels(csl, mgax, colors=[(0,0,0,0.8)], fmt='%.0f $^{\circ}C$')
    mgax.set_ylabel('$Mg/Ca_{cc}$ (mmol/mol)')

    # # gradient plot
    # vmin, vmax = Tgradlim
    # cmx = gax.pcolormesh(x, y, Tsens, vmin=vmin, vmax=vmax, cmap=cmg, zorder=-1)
    # # temperature contours
    # csl = gax.contour(x, y, Temp, contours, colors=[(0,0,0,0.3)], linewidths=1)
    # contour_labels(csl, gax, colors=[(0,0,0,0.8)], fmt='%.0f $^{\circ}C$')
    # gax.set_ylabel('$Mg/Ca_{cc}$ (mmol/mol)')
    
    return mgax # mgax, gax #, cbx, cbxg

def angle_of_line(x, y, ax):
    """
    Calculate angle of line on plot, in degrees.
    """
    dx = np.diff(x)[0]
    dy = np.diff(y)[0]
    
    tdx = (dx) / (ax.get_xlim()[1] - ax.get_xlim()[0])
    tdy = (dy) / (ax.get_ylim()[1] - ax.get_ylim()[0])
    
    fx, fy = ax.figure.get_size_inches()
    pos = ax.get_position()
    x_size = fx * (pos.x1 - pos.x0)
    y_size = fy * (pos.y1 - pos.y0)
    
    return np.rad2deg(np.arctan((tdy * y_size) / (tdx * x_size)))

def polypred(x, y, xnew, order):
    p = np.polyfit(x, y, order)
    return np.polyval(p, xnew)

def subset_trendline(ax, sub, xvar, color, line_order, yvar=('Measured', 'Mg/Caf'), xpad=0.2, label=None, label_x=None, label_v_offset_frac=0.05, lsargs={}, **kwargs):
    """
    Function to plot a trendline through a subset of some data
    """
    x = sub.loc[:, xvar]
    y = sub.loc[:, yvar]
    
    p = np.polyfit(x, y, line_order)
    xnew = np.linspace(x.min() - x.ptp() * xpad, x.max() + x.ptp() * xpad, 20)
    ax.plot(xnew, np.polyval(p, xnew), color=color, alpha=0.3, zorder=-2, **lsargs)
    
    if label is not None:
        dx = np.ptp(xnew) * 0.1

        if label_x > xnew.max():
            label_x_draw = xnew.max() + xnew.ptp() * 0.05
            label_x = xnew.max() - dx
            yoffset = 0
        else:
            label_x_draw = label_x
            yoffset = np.diff(ax.get_ylim()) * label_v_offset_frac

        
        xlo = label_x
        xhi = label_x + dx
        ylo = np.polyval(p, xlo)
        yhi = np.polyval(p, xhi)
        
        label_y = np.polyval(p, label_x_draw) 
        
        # tdx = (xhi - xlo) / (ax.get_xlim()[1] - ax.get_xlim()[0])
        # tdy = (yhi - ylo) / (ax.get_ylim()[1] - ax.get_ylim()[0])
        
        # # get axes absolute size
        # fx, fy = ax.figure.get_size_inches()
        # pos = ax.get_position()
        # x_size = fx * (pos.x1 - pos.x0)
        # y_size = fy * (pos.y1 - pos.y0)
        
        # label_angle = np.rad2deg(np.arctan((tdy * y_size) / (tdx * x_size)))

        label_angle = angle_of_line([xlo, xhi], [ylo, yhi], ax)
        
        ax.text(label_x_draw, label_y + yoffset, label,        
                rotation=label_angle, va='center', rotation_mode='anchor', color=color, alpha=0.4, **kwargs)

def emphasise_subset(ax, sub, xvar, color, m, s=None, yvar=('Measured', 'Mg/Caf')):
    x = sub.loc[:, xvar]
    y = sub.loc[:, yvar]
    if s is None:
        s = sub.loc[:, ('Measured', 'numberforams')]
    
    ax.scatter(x, y, marker=m, s=s, color=color, lw=0.5, edgecolor='k')
    