"""
Transfer functions used for converting Mg/Ca to temperature.
"""

import numpy as np
from foramgeochem.transfer import generic

# exponential, as in traditional palaeothermometer
def exp_mgca_2_temp(mgca, A, B):
    return generic.exponential(mgca, A, B)

def exp_temp_2_mgca(temp, A, B):
    return generic.exponential_inverse(temp, A, B)


# linear, for very simple caes
def lin_mgca_2_temp(mgca, A, B):
    return generic.linear(mgca, A, B)

def lin_temp_2_mgca(temp, A, B):
    return generic.linear_inverse(temp, A, B)


# Evans and Muller (2012)
# https://doi.org/10.1029/2012PA002315
def evans2012_temp_2_mgca(temp, A, B, H, mgca_sw_modern=5.17, mgca_sw_ancient=5.17):
    return B * (mgca_sw_ancient**H / mgca_sw_modern**H) * np.exp(A * temp)

def evans2012_mgca_2_temp(mgca, A, B, H, mgca_sw_modern=5.17, mgca_sw_ancient=5.17):
    return np.log((mgca * mgca_sw_modern**H) / (mgca_sw_ancient**H * B)) / A

def evans2012_mgca_2_mgca_sw(mgca, temp, A, B, H, mgca_sw_modern=5.17):
    return ((mgca * mgca_sw_modern**H) / (B * np.exp(A * temp)))**(1/H)


# Evans et al (2015)
# https://doi.org/10.1016/j.gca.2014.09.039
def evans2015_temp_2_mgca(temp, A=0.0168, B=94.8, H=2.0, C=-1.98, D=37.0, mgca_sw_modern=5.17, mgca_sw_ancient=5.17):
    return B * np.exp(A * temp) * (C * mgca_sw_ancient**H + D * mgca_sw_ancient) / (C * mgca_sw_modern**H + D * mgca_sw_modern)

def evans2015_mgca_2_temp(mgca, A=0.0168, B=94.8, H=2.0, C=-1.98, D=37.0, mgca_sw_modern=5.17, mgca_sw_ancient=5.17):
    return np.log(mgca * (C * mgca_sw_modern**H + D * mgca_sw_modern) / (B * (C * mgca_sw_ancient**H + D * mgca_sw_ancient))) / A

# def evans2015_mgca_2_mgca_sw(mgca, temp, A=0.0168, B=94.8, H=2.0, C=-1.98, D=37.0, mgca_sw_modern=5.17)
#     return


# Holland et al, 2018
def holland2018_calc_mgca(temp, mgca_sw=5.17, ca_sw=10.2, dic_sw=2000.0,
                          A=6.40680626e-01, B=5.28277336e-01,
                          C1=2.73373145e-04, C2=4.28061422e-06,
                          D=5.38367374e-02):
    return mgca_sw**A * B * np.exp((C1 * ca_sw + C2 * dic_sw + D) * temp)

def holland2018_calc_temp(mgca, mgca_sw=5.17, ca_sw=10.2, dic_sw=2000,
                          A=6.40680626e-01, B=5.28277336e-01,
                          C1=2.73373145e-04, C2=4.28061422e-06,
                          D=5.38367374e-02):
    return np.log(mgca / (mgca_sw**A * B)) / (C1 * ca_sw + C2 * dic_sw + D)

def holland2018_calc_mgca_sw(temp, mgca, ca_sw=10.2, dic_sw=2000,
                             A=6.40680626e-01, B=5.28277336e-01,
                             C1=2.73373145e-04, C2=4.28061422e-06,
                             D=5.38367374e-02):
    return (mgca / (B * np.exp((C1 * ca_sw + C2 * dic_sw + D) * temp)))**(1/A)

def holland2018_calc_dic_sw(temp, mgca, mgca_sw=5.17, ca_sw=10.2,
                            A=6.40680626e-01, B=5.28277336e-01,
                            C1=2.73373145e-04, C2=4.28061422e-06,
                            D=5.38367374e-02):
    return (np.log(mgca / (mgca_sw**A * B)) - temp * (C1 * ca_sw + D)) / (C2 * temp)

def holland2018_calc_Ca_sw(temp, mgca, mgca_sw=5.17, dic_sw=2000,
                           A=6.40680626e-01, B=5.28277336e-01,
                           C1=2.73373145e-04, C2=4.28061422e-06,
                           D=5.38367374e-02):
    return (np.log(mgca / (mgca_sw**A * B)) - temp * (C2 * dic_sw + D)) / (C1 * temp)
