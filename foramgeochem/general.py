"""
Functions and classes used for all proxy systems.
"""

import numpy as np
import uncertainties as un
import uncertainties.unumpy as unp
from textwrap import wrap, TextWrapper

from .helpers import ucheck, load_params

linewidth = 70

class params(object):
    """
    Container for parameters and their associated uncertainties.

    Attributes
    ----------
    values : uncertainties.AffineScalarFunc or uncertainties.ufloat
        The parameter values and their packaged uncertainties.
    nom_values : array_like
        The nominal values of the parameters (without uncertainties)
    cov : array_like
        The covariance matrix of the parameter uncertainties.
    stds : array_like
        Standard deviations of the parameters (**warning**: ignores covariance)
    """
    def __init__(self, values, cov=None, stds=None, info=None, labels=None):
        """
        Store parameters and associated uncertainties.

        Parameters
        ----------
        values : array_like, dict or uncertainties object
            Either the nominal values of the parameters, a dictionary
            containing `values` and `stds` or `cov` used to create
            a parameter object, or values with associated uncertainties
            created by the `uncertainties` package.
        cov : array_like
            The covariance matrix for the specified parameters used for
            propagating errors. Only used if `values` is array_like, in
            which case i should be an array of shape (len(values), len(values)).
        stds : array_like
            Standard deviations of the values. This does not account for
            correllated uncertainties - use `cov` for that. Only used if
            `values` is array_like, in which case `stds` should be the
            same length as `values`.
        info : str
            A description of the parameteters, to be stored alongside the
            values.
        labels : array_like
            Parameter names.
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
            if 'labels' in vd:
                labels = vd['labels']
        
        if info is None:
            info = 'No information given.'
        self.info = info        
        self.labels = labels

        self.param_name = ''

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

        self.wrap = TextWrapper(linewidth).wrap
    
    def __repr__(self):
        outstr = 'Parameter Set: {}\n  '.format(self.param_name)
        outstr += '\n  '.join(self.wrap(self.info)) + '\n\n'
        outstr += 'Parameter Values:\n  '
        if self.labels is not None:
            outstr += '\n  '.join('{}: {}'.format(lab, v) for lab, v in zip(self.labels, self.values))
        else:
            outstr += '\n  '.join('{}'.format(v) for v in self.values)
        return outstr
    
    @staticmethod
    def load(proxy=None, mode=None, parameters=None, json_path=None):
        """
        Load one of the pre-defined parameter sets.

        To see a list of available options and their information, use
        `available_parameters()`.

        Parameters
        ----------
        proxy : str
            The name of the proxy system you're working in.
        mode : str
            The 'mode' of the proxy system - e.g. 'exponential' for Mg/Ca
        parameters : str
            The name of the parameter set you want to use.
        json_path : str
            The location of the json file containing defined parameters.
            If not specified, uses the built-in file included in `foramgeochem`.

        Returns
        -------
        Parameters and their associated uncertainties. : `foramgeochem.general.params` object
        """
        ps = load_params(json_path)
            
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
        
        p = params(values=ps)
        p.param_name = parameters
        
        return p

    @staticmethod
    def available_parameters(proxy=None, mode=None, parameters=None, json_path=None):
        """
        View all available parameter sets, and associated info.

        Arguments can optionally be provided if you want to see a subset of available
        parameters. For example, specify `proxy='mgca'` to just see modes and parameters
        for the Mg/Ca proxy.

        Parameters
        ----------
        proxy : str
            The name of the proxy system you're working in.
        mode : str
            The 'mode' of the proxy system - e.g. 'exponential' for Mg/Ca
        parameters : str
            The name of the parameter set you want to use.
        json_path : str
            The location of the json file containing defined parameters.
            If not specified, uses the built-in file included in `foramgeochem`.

        Returns
        -------
        None
        """
        bullet = chr(135)

        if proxy is not None:
            if isinstance(proxy, str):
                proxy = [proxy]
        if mode is not None:
            if isinstance(mode, str):
                mode = [mode]
        if parameters is not None:
            if isinstance(parameters, str):
                parameters = [parameters]

        ps = load_params(json_path)

        for k0, v0 in ps.items():
            if proxy is not None:
                if k0 not in proxy:
                    continue
            pname = '{} {} [proxy]'.format(bullet, k0)
            print('-' * (linewidth + 5))
            print(pname)
            print('-' * len(pname))

            for k1, v1 in v0.items():
                if mode is not None:
                    if k1 not in mode:
                        continue
                mname = ' {} {} [mode]'.format(bullet, k1)
                print(mname)
                print(' ' + '-' * (len(mname) - 1))

                for k2, v2 in v1.items():
                    if parameters is not None:
                        if k2 not in parameters:
                            continue
                    print('  {} {} [parameters]'.format(bullet, k2))
                    print('     ' + '\n     '.join(wrap(v2['info'], linewidth)))

        print('-' * (linewidth + 5))

class proxy(object):
    """
    General framework and helpers for proxy calculation.
    """
    def __init__(self):
        self.variables = set()
        self.missing = set()
        self.last_calc = None
        
        self.wrap = TextWrapper(linewidth).wrap

        self.fn_name = ''
        self.fn_text = ''

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
            if v is not None:
                v = np.asanyarray(v)

            if not hasattr(self, k):
                setattr(self, k, v)
            elif v is not None:
                setattr(self, k, v)
        
        self._var_check()
    
    
