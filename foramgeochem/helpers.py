import pkg_resources as pkgrs
import json

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
    
def load_params(json_path=None):
    """
    Loads parameters .json file.
    """
    if json_path is None:
        json_path = pkgrs.ResourceManager().resource_filename("foramgeochem", "resources/params.json")
    with open(json_path) as f:
        ps = json.load(f)
    return ps