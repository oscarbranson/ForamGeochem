"""
Functions and classes used for all proxy systems.
"""

import numpy as np
import uncertainties as un
import uncertainties.unumpy as unp
from .helpers import ucheck

import pkg_resources as pkgrs
import json

class params(object):
    """
    Container for parameters and their associated uncertainties.
    """
    def __init__(self, values, stds=None, cov=None, info=None):
        """
        Store parameters and associated uncertainties.
        """
        if isinstance(values, dict):
            vd = values.copy()
            values = vd['values']
            if 'stds' in vd:
                stds = vd['stds']
            if 'cov' in vd:
                cov = vd['cov']
            if 'info' in vd:
                info = vd['info']
        
        if info is None:
            info = 'No information given.'
        self.info = info

        if ucheck(values):
            self.values = values
        else:
            if cov is not None:
                self.values = un.correlated_values(values, cov)
            elif stds is not None:
                self.values = unp.uarray(values, stds)
            else:
                self.values = values

        if ucheck(self.values):
            self.nom_values = unp.nominal_values(self.values)
            self.stds = unp.std_devs(self.values)
            self.cov = un.covariance_matrix(self.values)
        else:
            self.nom_values = np.asanyarray(self.values)
            self.stds = np.full(len(self.values), np.nan)
            self.stds = np.full((len(self.values),len(self.values)), np.nan)
    
    def __repr__(self):
        outstr = self.info + '\n\n'
        outstr += 'Parameter Values:\n  '
        outstr += '\n  '.join('{}'.format(v) for v in self.values)
        return outstr
    
    @staticmethod
    def load(json_path=None, proxy=None, mode=None, parameters=None):
        if json_path is None:
            json_path = pkgrs.ResourceManager().resource_filename("foramgeochem", "resources/params.json")
        with open(json_path) as f:
            ps = json.load(f)
            
        try:
            ps = ps[proxy]
        except:
            raise KeyError("`proxy` incorrect. Try one of: ['{}']".format("', '".join(ps.keys())))

        try:
            ps = ps[mode]
        except:
            raise KeyError("`mode` incorrect. Try one of: ['{}']".format("', '".join(ps.keys())))
        
        try:
            ps = ps[parameters]
        except:
            raise KeyError("`parameters` incorrect. Try one of: ['{}']".format("', '".join(ps.keys())))
        
        return params(values=ps)

class proxy(object):
    """
    General framework and helpers for proxy calculation.
    """
    def __init__(self):
        self.variables = set()
        self.missing = set()
        self.last_calc = None
        
    def _var_check(self):
        """
        Check which of the required variables is missing.
        """
        missing = set()
        for v in self.variables:
            if getattr(self, v) is None:
                missing.add(v)
        self.missing = missing
                
    def _var_update(self, **kwargs):
        """
        Update the values of all variables specified as varname=value.
        """
        for k, v in kwargs.items():
            if not hasattr(self, k):
                setattr(self, k, v)
            elif v is not None:
                setattr(self, k, v)
        
        self._var_check()
    
    
