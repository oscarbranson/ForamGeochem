"""
Transfer functions used for converting between Mg/Ca and environmental parameters.
"""

from uncertainties.unumpy import log, exp

# exponential, as in traditional palaeothermometer
def exp_mgca_2_temp(mgca_f, p=None):
    """
    Calculate foram Mg/Ca from temperature using the 'classic' exponential function.

    Parameters
    ----------
    mgca_f : array_like
        Foraminiferal Mg/Ca in mmol/mol
    p : array_like
        Parameters for the model, in the order [A, B].
    
    Returns
    -------
    Temperature in Celcius. : array_like
    """
    A, B = p
    return log(mgca_f / A) / B 

def exp_temp_2_mgca(temp, p=None):
    """
    Calculate temperature from foram Mg/Ca using the 'classic' exponential function.

    Parameters
    ----------
    temp : array_like
        Temperature in Celcius.
    p : array_like
        Parameters for the model, in the order [A, B].
    
    Returns
    -------
    Foraminiferal Mg/Ca in mmol/mol : array_like
    """
    A, B = p
    return A * exp(temp * B)


# Evans and Muller (2012)
# https://doi.org/10.1029/2012PA002315
def evans2012_temp_2_mgca(temp, mgca_sw_modern=5.17, mgca_sw_ancient=5.17, p=None):
    """
    Calculate foram Mg/Ca from temperature, corrected for past change in seawater Mg/Ca.

    Following Evans & Muller (2012) https://doi.org/10.1029/2012PA002315

    Parameters
    ----------
    temp : array_like
        Temperature in Celcius.
    mgca_sw_modern : array_like
        Mg/Ca of modern seawater in mol/mol.
    mgca_sw_ancient : array_like
        Mg/Ca of ancient seawater in mol/mol.

    Returns
    -------
    Foraminiferal Mg/Ca in mmol/mol : array_like
    """
    A, B, H = p
    return B * (mgca_sw_ancient**H / mgca_sw_modern**H) * exp(A * temp)

def evans2012_mgca_2_temp(mgca_f, mgca_sw_modern=5.17, mgca_sw_ancient=5.17, p=None):
    """
    Calculate temperature from foram Mg/Ca, corrected for past change in seawater Mg/Ca.

    Following Evans & Muller (2012) https://doi.org/10.1029/2012PA002315

    Parameters
    ----------
    mgca_f : array_like
        Foraminiferal Mg/Ca in mmol/mol
    mgca_sw_modern : array_like
        Mg/Ca of modern seawater in mol/mol.
    mgca_sw_ancient : array_like
        Mg/Ca of ancient seawater in mol/mol.

    Returns
    -------
    Temperature in Celcius. : array_like
    """
    A, B, H = p
    return log((mgca_f * mgca_sw_modern**H) / (mgca_sw_ancient**H * B)) / A

def evans2012_mgca_2_mgca_sw(mgca_f, temp, mgca_sw_modern=5.17, p=None):
    """
    Calculate past seawater Mg/Ca from temperature and foram Mg/Ca.

    Following Evans & Muller (2012) https://doi.org/10.1029/2012PA002315

    Parameters
    ----------
    mgca_f : array_like
        Foraminiferal Mg/Ca in mmol/mol
    temp : array_like
        Temperature in Celcius.
    mgca_sw_modern : array_like
        Mg/Ca of modern seawater in mol/mol.
    
    Returns
    -------
    Mg/Ca of ancient seawater in mol/mol : array_like
    """
    A, B, H = p
    return ((mgca_f * mgca_sw_modern**H) / (B * exp(A * temp)))**(1/H)


# Evans et al (2015)
# https://doi.org/10.1016/j.gca.2014.09.039
def evans2015_temp_2_mgca(temp, mgca_sw_modern=5.17, mgca_sw_ancient=5.17, p=(0.0168, 94.8, 2.0, 1.98, 37.0)):
    """
    Calculate foram Mg/Ca from temperature, corrected for past change in seawater Mg/Ca.

    Following Evans et al (2015) https://doi.org/10.1016/j.gca.2014.09.039

    Parameters
    ----------
    temp : array_like
        Temperature in Celcius.
    mgca_sw_modern : array_like
        Mg/Ca of modern seawater in mol/mol.
    mgca_sw_ancient : array_like
        Mg/Ca of ancient seawater in mol/mol.

    Returns
    -------
    Foraminiferal Mg/Ca in mmol/mol : array_like
    """
    A, B, H, C, D = p
    return B * exp(A * temp) * (C * mgca_sw_ancient**H + D * mgca_sw_ancient) / (C * mgca_sw_modern**H + D * mgca_sw_modern)

def evans2015_mgca_2_temp(mgca_f, mgca_sw_modern=5.17, mgca_sw_ancient=5.17, p=(0.0168, 94.8, 2.0, 1.98, 37.0)):
    """
    Calculate temperature from foram Mg/Ca, corrected for past change in seawater Mg/Ca.

    Following Evans et al (2015) https://doi.org/10.1016/j.gca.2014.09.039

    Parameters
    ----------
    mgca_f : array_like
        Foraminiferal Mg/Ca in mmol/mol
    mgca_sw_modern : array_like
        Mg/Ca of modern seawater in mol/mol.
    mgca_sw_ancient : array_like
        Mg/Ca of ancient seawater in mol/mol.

    Returns
    -------
    Temperature in Celcius. : array_like
    """
    A, B, H, C, D = p
    return log(mgca_f * (C * mgca_sw_modern**H + D * mgca_sw_modern) / (B * (C * mgca_sw_ancient**H + D * mgca_sw_ancient))) / A


# Holland et al, 2020
def holland2020_calc_mgca(temp, mgca_sw=5.17, ca_sw=10.2e-3, carb_sw=2000e-6, p=None):
    """
    Calculate foram Mg/Ca from Temperature, and seawater Mg/Ca, [Ca] and carbon chemistry.

    Parameters
    ----------
    temp : array_like
        Temperature in Celcius.
    mgca_sw : array_like
        Seawater Mg/Ca, in mol/mol.
    ca_sw : array_like
        Seawater calcium concentration, in mol kg-1.
    carb_sw : array_like
        Seawater carbon parameter - either DIC, CO3 or pH, depending on the species.
    p : array_like
        Parameters for the model, in the order [A, B, C, D, E].

    Returns
    -------
    Foraminiferal Mg/Ca in mmol/mol : array_like
    """
    A, B, C, D, E = p
    return mgca_sw**A * carb_sw**B * exp(C * ca_sw + D * temp + E)

def holland2020_calc_temp(mgca, mgca_sw=5.17, ca_sw=10.2e-3, carb_sw=2000e-6, p=None):
    """
    Calculate Temperature from foram Mg/Ca, and seawater Mg/Ca, [Ca] and carbon chemistry.

    Parameters
    ----------
    mgca : array_like
        Foraminiferal Mg/Ca in mmol/mol
    mgca_sw : array_like
        Seawater Mg/Ca, in mol/mol.
    ca_sw : array_like
        Seawater calcium concentration, in mol kg-1.
    carb_sw : array_like
        Seawater carbon parameter - either DIC, CO3 or pH, depending on the species.
    p : array_like
        Parameters for the model, in the order [A, B, C, D, E].

    Returns
    -------
    Temperature in Celcius. : array_like
    """
    A, B, C, D, E = p
    return (log(mgca / (mgca_sw**A * carb_sw**B)) - ca_sw * C - E) / D

def holland2020_calc_mgca_sw(temp, mgca, ca_sw=10.2e-3, carb_sw=2000e-6, p=None):
    """
    Calculate seawater Mg/Ca from Temperature, foram Mg/Ca, and seawater [Ca] and carbon chemistry.

    Parameters
    ----------
    temp : array_like
        Temperature in Celcius.
    mgca : array_like
        Foraminiferal Mg/Ca in mmol/mol
    ca_sw : array_like
        Seawater calcium concentration, in mol kg-1.
    carb_sw : array_like
        Seawater carbon parameter - either DIC, CO3 or pH, depending on the species.
    p : array_like
        Parameters for the model, in the order [A, B, C, D, E].

    Returns
    -------
    Seawater Mg/Ca, in mol/mol. : array_like
    """
    A, B, C, D, E = p
    return (mgca / (carb_sw**B * exp(C * ca_sw + D * temp + E)))**(1/A)

def holland2020_calc_carb_sw(temp, mgca, mgca_sw=5.17, ca_sw=10.2e-3, p=None):
    """
    Calculate seawater carbon chemistry from Temperature, foram Mg/Ca, and seawater Mg/Ca and [Ca].

    Parameters
    ----------
    temp : array_like
        Temperature in Celcius.
    mgca : array_like
        Foraminiferal Mg/Ca in mmol/mol
    mgca_sw : array_like
        Seawater Mg/Ca, in mol/mol.
    ca_sw : array_like
        Seawater calcium concentration, in mol kg-1.
    p : array_like
        Parameters for the model, in the order [A, B, C, D, E].

    Returns
    -------
    Seawater carbon parameter - either DIC, CO3 or pH, depending on the species. : array_like
    """
    A, B, C, D, E = p
    return (mgca / (mgca_sw**A * exp(C * ca_sw + D * temp + E)))**(1/B)

def holland2020_calc_Ca_sw(temp, mgca, mgca_sw=5.17, carb_sw=2000e-6, p=None):
    """
    Calculate seawater calcium concentration from Temperature, foram Mg/Ca, and seawater Mg/Ca and carbon chemistry.

    Parameters
    ----------
    temp : array_like
        Temperature in Celcius.
    mgca : array_like
        Foraminiferal Mg/Ca in mmol/mol
    mgca_sw : array_like
        Seawater Mg/Ca, in mol/mol.
    carb_sw : array_like
        Seawater carbon parameter - either DIC, CO3 or pH, depending on the species.
    p : array_like
        Parameters for the model, in the order [A, B, C, D, E].

    Returns
    -------
    Seawater calcium concentration, in mol kg-1. : array_like
    """
    A, B, C, D, E = p
    return (log(mgca / (mgca_sw**A * carb_sw**B)) - temp * D - E) / C