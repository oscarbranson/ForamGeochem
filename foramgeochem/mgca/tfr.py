"""
Transfer functions used for converting between Mg/Ca and environmental parameters.
"""

from uncertainties.unumpy import log, exp

# exponential, as in traditional palaeothermometer
def exp_mgca_2_temp(mgca, p=None):
    A, B = p
    return log(mgca / A) / B 

def exp_temp_2_mgca(temp, p=None):
    A, B = p
    return A * exp(temp * B)


# linear, for very simple cases
def lin_mgca_2_temp(mgca, p=None):
    A, B = p
    return B + mgca * A

def lin_temp_2_mgca(temp, p=None):
    A, B = p
    return (temp - B) / A


# Evans and Muller (2012)
# https://doi.org/10.1029/2012PA002315
def evans2012_temp_2_mgca(temp, A, B, H, mgca_sw_modern=5.17, mgca_sw_ancient=5.17):
    return B * (mgca_sw_ancient**H / mgca_sw_modern**H) * exp(A * temp)

def evans2012_mgca_2_temp(mgca, A, B, H, mgca_sw_modern=5.17, mgca_sw_ancient=5.17):
    return log((mgca * mgca_sw_modern**H) / (mgca_sw_ancient**H * B)) / A

def evans2012_mgca_2_mgca_sw(mgca, temp, A, B, H, mgca_sw_modern=5.17):
    return ((mgca * mgca_sw_modern**H) / (B * exp(A * temp)))**(1/H)


# Evans et al (2015)
# https://doi.org/10.1016/j.gca.2014.09.039
def evans2015_temp_2_mgca(temp, A=0.0168, B=94.8, H=2.0, C=-1.98, D=37.0, mgca_sw_modern=5.17, mgca_sw_ancient=5.17):
    return B * exp(A * temp) * (C * mgca_sw_ancient**H + D * mgca_sw_ancient) / (C * mgca_sw_modern**H + D * mgca_sw_modern)

def evans2015_mgca_2_temp(mgca, A=0.0168, B=94.8, H=2.0, C=-1.98, D=37.0, mgca_sw_modern=5.17, mgca_sw_ancient=5.17):
    return log(mgca * (C * mgca_sw_modern**H + D * mgca_sw_modern) / (B * (C * mgca_sw_ancient**H + D * mgca_sw_ancient))) / A

# def evans2015_mgca_2_mgca_sw(mgca, temp, A=0.0168, B=94.8, H=2.0, C=-1.98, D=37.0, mgca_sw_modern=5.17)
#     return


# Holland et al, 2018
def holland2018_calc_mgca(temp, mgca_sw=5.17, ca_sw=10.2e-3, dic_sw=2000e-6, p=None):
    A, B, C1, C2, D = p
    return mgca_sw**A * B * exp((C1 * ca_sw + C2 * dic_sw + D) * temp)

def holland2018_calc_temp(mgca, mgca_sw=5.17, ca_sw=10.2e-3, dic_sw=2000e-6, p=None):
    A, B, C1, C2, D = p
    return log(mgca / (mgca_sw**A * B)) / (C1 * ca_sw + C2 * dic_sw + D)

def holland2018_calc_mgca_sw(temp, mgca, ca_sw=10.2e-3, dic_sw=2000e-6, p=None):
    A, B, C1, C2, D = p
    return (mgca / (B * exp((C1 * ca_sw + C2 * dic_sw + D) * temp)))**(1/A)

def holland2018_calc_dic_sw(temp, mgca, mgca_sw=5.17, ca_sw=10.2e-3, p=None):
    A, B, C1, C2, D = p
    return (log(mgca / (mgca_sw**A * B)) - temp * (C1 * ca_sw + D)) / (C2 * temp)

def holland2018_calc_Ca_sw(temp, mgca, mgca_sw=5.17, dic_sw=2000e-6, p=None):
    A, B, C1, C2, D = p
    return (log(mgca / (mgca_sw**A * B)) - temp * (C2 * dic_sw + D)) / (C1 * temp)
