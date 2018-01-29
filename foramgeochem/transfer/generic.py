"""
Generic transfer functions used in palaeoproxies.
"""

import numpy as np

def exponential(x, A, B):
    return np.log(x / A) / B 

def exponential_inverse(y, A, B):
    return A * np.exp(y * B)

def linear(x, A, B):
    return A + x * B

def linear_inverse(y, A, B):
    return (y - A) / B