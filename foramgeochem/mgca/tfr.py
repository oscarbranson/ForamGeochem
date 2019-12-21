"""
Transfer functions used for converting between Mg/Ca and environmental parameters.
"""

from foramgeochem.mgca.cv_holland2020 import holland2020_calc_Ca_sw, holland2020_calc_carb_sw, holland2020_calc_mgca, holland2020_calc_mgca_sw, holland2020_calc_temp
from foramgeochem.mgca.cv_exponential import exp_mgca_2_temp, exp_temp_2_mgca

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


