import numpy as np

def limits(value, tolerance):
    return value - value * tolerance, value + value * tolerance

def isolate_constant_conditions(dat, Temp=None, Mg=None, Ca=None, MgCa=None, DIC=None, pH=None,
                                tolerance=0.1, pH_tolerance=0.2):
    ind = np.ones(dat.shape[0], dtype=bool)
    
    if Temp is not None:
        lo, hi = limits(Temp, tolerance)
        ind = ind & (dat.loc[:, ('Measured', 'Temp')] >= lo) & (dat.loc[:, ('Measured', 'Temp')] <= hi)
    if Mg is not None:
        lo, hi = limits(Mg, tolerance)
        ind = ind & (dat.loc[:, ('Measured', '[Mg]sw')] >= lo) & (dat.loc[:, ('Measured', '[Mg]sw')] <= hi)
    if Ca is not None:
        lo, hi = limits(Ca, tolerance)
        ind = ind & (dat.loc[:, ('Measured', '[Ca]sw')] >= lo) & (dat.loc[:, ('Measured', '[Ca]sw')] <= hi)
    if MgCa is not None:
        lo, hi = limits(MgCa, tolerance)
        ind = ind & (dat.loc[:, ('Measured', 'Mg/Casw')] >= lo) & (dat.loc[:, ('Measured', 'Mg/Casw')] <= hi)
    if DIC is not None:
        lo, hi = limits(DIC, tolerance)
        ind = ind & (dat.loc[:, ('csys_mid', 'DIC')] >= lo) & (dat.loc[:, ('csys_mid', 'DIC')] <= hi)
    if pH is not None:
#         hi, lo = -np.log10(limits(10**-pH, tolerance * 2))
        hi, lo = pH + pH_tolerance, pH - pH_tolerance
        ind = ind & (dat.loc[:, ('csys_mid', 'pHtot')] >= lo) & (dat.loc[:, ('csys_mid', 'pHtot')] <= hi)
    
    return dat.loc[ind, :]

def gauss(x, *p):
    """ Gaussian function.

    Parameters
    ----------
    x : array_like
        Independent variable.
    *p : parameters unpacked to A, mu, sigma
        A: area
        mu: centre
        sigma: width

    Return
    ------
    array_like
        gaussian descriped by *p.
    """
    A, mu, sigma = p
    return A * np.exp(-0.5 * (-mu + x)**2 / sigma**2)

def weighted_moving_average(x, y, x_new, fwhm=300):
    """
    Calculate gaussian weigted moving mean, SD and SE.

    Parameters
    ----------
    x, y : array - like
        The x and y data to smooth
    x_new : array - like
        The new x - scale to interpolate the data

    """
    bin_avg = np.zeros(len(x_new))
    bin_std = np.zeros(len(x_new))
    bin_se = np.zeros(len(x_new))

    # Gaussian function as weights
    sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))

    for index, xn in enumerate(x_new):
        weights = gauss(x, 1, xn, sigma)
        weights /= sum(weights)
        # weighted mean
        bin_avg[index] = np.average(y, weights=weights)
        # weighted standard deviation
        bin_std[index] = np.sqrt(np.average((y - bin_avg[index])**2, weights=weights))
        # weighted standard error (mean / sqrt(n_points_in_gaussian))
        bin_se[index] = np.sqrt(np.average((y - bin_avg[index])**2, weights=weights)) / \
            np.sqrt(sum((x > xn - 2 * sigma) & (x < xn + 2 * sigma)))

    return {'mean': bin_avg,
            'std': bin_std,
            'stderr': bin_se}