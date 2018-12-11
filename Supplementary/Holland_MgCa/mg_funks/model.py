import numpy as np
import uncertainties as un
from uncertainties import unumpy as unp

def T_fn(x, A, B, C1, C2, D):
    """
    Calculate temperature from Mg/Casw, [Ca]sw, [C] and Mg/Caforam.
    """
    if any([isinstance(p, un.core.AffineScalarFunc) for p in [A, B, C1, C2, D]]):
        log = unp.log
    else:
        log = np.log
    MgCasw, Casw, Carb, MgCacc = x
    return log(MgCacc / (MgCasw**A * B)) / (Casw * C1 + C2 * Carb + D)