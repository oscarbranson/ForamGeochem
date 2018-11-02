import uncertainties as un
import uncertainties.unumpy as unp

def ucheck(v):
    """
    Checks if v (or any element of v) is an `uncertainties` object.
    """
    if hasattr(v, '__iter__'):
        return any([isinstance(vi, (un.core.AffineScalarFunc, un.core.Variable)) for vi in v])
    else:
        return isinstance(v, (un.core.AffineScalarFunc, un.core.Variable))
    