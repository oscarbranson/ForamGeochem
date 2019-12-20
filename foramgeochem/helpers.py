import pkg_resources as pkgrs
import json
import os
import warnings
import numpy as np

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

def save_params(ps, json_path=None, overwrite=False):
    """
    Saves parameters .json file.
    """
    if json_path is None:
        json_path = pkgrs.ResourceManager().resource_filename("foramgeochem", "resources/params.json")
    if os.path.exists(json_path) and not overwrite:
        raise ValueError(f'{json_path} exists. Set overwrite=True if you really want to lose the old one.')
    else:
        with open(json_path, 'w') as f:
            json.dump(ps, f, indent=2)
        print(f'Parameters saved to {json_path}')

def update_params(proxy, calib_ref, calib_type, p_values=None, p_cov=None, p_labels=None, p_info=None, create_missing=False, remove_empty=False, param_file=None, save_file=False, overwrite=False):
    """
    Update or add proxy parameters to the database.

    Parameters
    ==========
    proxy : str
        The name of the proxy that the parameters belong to, for example 'mgca'.
    calib_ref : str
        The reference that the parameters come from, for example 'holland_2020'.
    calib_type : str
        The type of calibration that the parameters are derived from, for example 'O_universa'.
        This should be a short string that readily identifies the parameter set.
    p_values : array-like
        The parameter values [required]. 
    p_cov : array-like
        The covariance matrix for the parameter values. This is not required, but is
        *strongly* recommended. Without this, parameter uncertainties will not be
        propagated.
    p_labels : array-like
        The labels of the parameter values. Ideally, these should be the same as in
        the original paper, so values may be readily referenced to the literature.
    p_info : str
        A short description of the parameter set, which tells the user a bit about
        the parameters and where they come from. You might include important notes
        here such as "not suitable for application to benthics"
    create_missing : bool
        If True, new entries will be created in the database if they are not already present.
    remove_empty : bool
        If True, p_cov, p_labels and p_info are removed from an existing entry if not specified
        in the function.
    param_file : str or None
        The parameter file to modify. If None, the default parameter file is used.
    save_file : bool or str
        If True, parameters are written to the same file as the input. If a string,
        parameters are written to the new file described by the string. If False,
        the function returns the modified parameter dict.
    overwrite : bool
        If True, the parameter file will be overwritten.

    Returns
    =======
    dict : Parameter dictionary, if not saved.
    """
    params = load_params(param_file)

    if proxy not in params:
        valid_list = "\n   - ".join(params.keys())
        if create_missing:
            params[proxy] = {}
        else:
            raise ValueError(f'{proxy} is not currently in the parameter database.\nChoose one of: \n   - {valid_list}\nOr create a new entry by setting create_missing=True.')
    if calib_ref not in params[proxy]:
        valid_list = "\n   - ".join(params[proxy].keys())
        if create_missing:
            params[proxy][calib_ref] = {}
        else:
            raise ValueError(f'{calib_ref} is not currently in the parameter database for {proxy}.\nChoose one of: \n   - {valid_list}\nOr create a new entry by setting create_missing=True.')
    if calib_type not in params[proxy][calib_ref]:
        valid_list = "\n   - ".join(params[proxy][calib_ref].keys())
        if create_missing:
            params[proxy][calib_ref][calib_type] = {}
        else:
            raise ValueError(f'{calib_type} is not currently in the parameter database for {proxy}:{calib_ref}.\nChoose one of: \n   - {valid_list}\nOr create a new entry by setting create_missing=True.')
    
    
    if p_values is None:
        raise ValueError('At minimum, p_values must be given')
    else:
        params[proxy][calib_ref][calib_type]['values'] = np.asanyarray(p_values).tolist()
    
    if p_cov is None:
        warnings.warn('p_cov is not specified. Uncertainties will not be propagated for this proxy!!')
        if not remove_empty:
            warnings.warn('Old p_cov retained. Set remove_empty=True to remove them.')
        else:
            del params[proxy][calib_ref][calib_type]['cov']
    else:
        params[proxy][calib_ref][calib_type]['cov'] = np.asanyarray(p_cov).tolist()
    
    if p_info is None:
        warnings.warn('info is not specified. How will you know where these parameters are from in future?')
        if not remove_empty:
            warnings.warn('Old info retained. Set remove_empty=True to remove them.')
        else:
            del params[proxy][calib_ref][calib_type]['info']
    else:
        params[proxy][calib_ref][calib_type]['info'] = p_info
    
    if p_labels is None:
        warnings.warn('p_labels is not specified. How will you know what these parameters are in future?')
        if not remove_empty:
            warnings.warn('Old p_labels retained. Set remove_empty=True to remove them.')
        else:
            del params[proxy][calib_ref][calib_type]['labels']
    else:
        params[proxy][calib_ref][calib_type]['labels'] = np.asanyarray(p_labels).tolist()
    
    if isinstance(save_file, str):
        save_params(params, json_file=save_file, overwrite=overwrite)
    elif save_file:
        save_params(params, overwrite=overwrite)
    else:
        return params