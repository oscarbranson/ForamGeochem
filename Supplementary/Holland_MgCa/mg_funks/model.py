import numpy as np
import uncertainties as un
from uncertainties import unumpy as unp

def T_fn(x, A, B, C, D, E):
    """
    Calculate temperature from Mg/Casw, [Ca]sw, [C] and Mg/Caforam.
    """
    if any([isinstance(p, un.core.AffineScalarFunc) for p in [A, B, C, D, E]]):
        log = unp.log
    else:
        log = np.log
    MgCasw, Ca, DIC, MgCaf = x
    return (log(MgCaf / (MgCasw**A * DIC**B)) - Ca * C - E) / D
