"""
Multi-parameter relationship between Mg/Ca, [Ca]sw, Mg/Casw, Carbon Chemistry and Temperature of Holland et al. (2020).
"""
import foramgeochem
from foramgeochem.general import proxy, params
from foramgeochem.helpers import load_params
from uncertainties.unumpy import log, exp

# Holland et al, 2020
def holland2020_calc_mgca(temp, mgca_sw=5.17, ca_sw=10.2e-3, carb_sw=2000e-6, p=None):
    """
    Calculate foram Mg/Ca from Temperature, and seawater Mg/Ca, [Ca] and carbon chemistry.

    Parameters
    ----------
    temp : array_like
        Temperature in Celcius.
    mgca_sw : array_like
        Seawater Mg/Ca, in mol/mol.
    ca_sw : array_like
        Seawater calcium concentration, in mol kg-1.
    carb_sw : array_like
        Seawater carbon parameter - either DIC, CO3 or pH, depending on the species.
    p : array_like
        Parameters for the model, in the order [A, B, C, D, E].

    Returns
    -------
    Foraminiferal Mg/Ca in mmol/mol : array_like
    """
    A, B, C, D, E = p
    return mgca_sw**A * carb_sw**B * exp(C * ca_sw + D * temp + E)

def holland2020_calc_temp(mgca, mgca_sw=5.17, ca_sw=10.2e-3, carb_sw=2000e-6, p=None):
    """
    Calculate Temperature from foram Mg/Ca, and seawater Mg/Ca, [Ca] and carbon chemistry.

    Parameters
    ----------
    mgca : array_like
        Foraminiferal Mg/Ca in mmol/mol
    mgca_sw : array_like
        Seawater Mg/Ca, in mol/mol.
    ca_sw : array_like
        Seawater calcium concentration, in mol kg-1.
    carb_sw : array_like
        Seawater carbon parameter - either DIC, CO3 or pH, depending on the species.
    p : array_like
        Parameters for the model, in the order [A, B, C, D, E].

    Returns
    -------
    Temperature in Celcius. : array_like
    """
    A, B, C, D, E = p
    return (log(mgca / (mgca_sw**A * carb_sw**B)) - ca_sw * C - E) / D

def holland2020_calc_mgca_sw(temp, mgca, ca_sw=10.2e-3, carb_sw=2000e-6, p=None):
    """
    Calculate seawater Mg/Ca from Temperature, foram Mg/Ca, and seawater [Ca] and carbon chemistry.

    Parameters
    ----------
    temp : array_like
        Temperature in Celcius.
    mgca : array_like
        Foraminiferal Mg/Ca in mmol/mol
    ca_sw : array_like
        Seawater calcium concentration, in mol kg-1.
    carb_sw : array_like
        Seawater carbon parameter - either DIC, CO3 or pH, depending on the species.
    p : array_like
        Parameters for the model, in the order [A, B, C, D, E].

    Returns
    -------
    Seawater Mg/Ca, in mol/mol. : array_like
    """
    A, B, C, D, E = p
    return (mgca / (carb_sw**B * exp(C * ca_sw + D * temp + E)))**(1/A)

def holland2020_calc_carb_sw(temp, mgca, mgca_sw=5.17, ca_sw=10.2e-3, p=None):
    """
    Calculate seawater carbon chemistry from Temperature, foram Mg/Ca, and seawater Mg/Ca and [Ca].

    Parameters
    ----------
    temp : array_like
        Temperature in Celcius.
    mgca : array_like
        Foraminiferal Mg/Ca in mmol/mol
    mgca_sw : array_like
        Seawater Mg/Ca, in mol/mol.
    ca_sw : array_like
        Seawater calcium concentration, in mol kg-1.
    p : array_like
        Parameters for the model, in the order [A, B, C, D, E].

    Returns
    -------
    Seawater carbon parameter - either DIC, CO3 or pH, depending on the species. : array_like
    """
    A, B, C, D, E = p
    return (mgca / (mgca_sw**A * exp(C * ca_sw + D * temp + E)))**(1/B)

def holland2020_calc_Ca_sw(temp, mgca, mgca_sw=5.17, carb_sw=2000e-6, p=None):
    """
    Calculate seawater calcium concentration from Temperature, foram Mg/Ca, and seawater Mg/Ca and carbon chemistry.

    Parameters
    ----------
    temp : array_like
        Temperature in Celcius.
    mgca : array_like
        Foraminiferal Mg/Ca in mmol/mol
    mgca_sw : array_like
        Seawater Mg/Ca, in mol/mol.
    carb_sw : array_like
        Seawater carbon parameter - either DIC, CO3 or pH, depending on the species.
    p : array_like
        Parameters for the model, in the order [A, B, C, D, E].

    Returns
    -------
    Seawater calcium concentration, in mol kg-1. : array_like
    """
    A, B, C, D, E = p
    return (log(mgca / (mgca_sw**A * carb_sw**B)) - temp * D - E) / C

class Holland2020(proxy):
    """
    The multi-factor modified exponential relationship between Mg/Ca, temperature, carbon, [Ca] and Mg/Casw of Holland et al (2020)
    """
    def __init__(self, mgca=None, temperature=None, carb_sw=2100e-6, ca_sw=10.2e-3, mgca_sw=5.0, parameters=None):
        """
        Create an object for converting between Mg/Ca and environmental parameters.

        Uses the formulation of Holland et al (2020). Different parameter sets discussed
        in the paper are available via the 'parameters' argument.

        mgca = mgca_sw**A * B exp((C1 * ca_sw + C2 * carb_sw + D) * temperature)

        Parameters
        ----------
        mgca : array_like
            Foraminiferal Mg/Ca in mmol/mol
        temperature : array_like
            Temperature in celcius.
        carb_sw : array_like
            Seawater carbon parameter - either DIC, CO3 or pH, depending on the species.
        mgca_sw : array_like
            Seawater Mg/Ca, in mol/mol.
        ca_sw : array_like
            Seawater calcium concentration, in mol kg-1. 
        parameters : array_like, str or `params` object
            Either a 'params' object containing parameter values and associated unctertainties,
            a string selecting one of the in-built options, or an array_like of numbers
            to use as parameters. To see a list of pre-defined parameters, use:
               `Holland.available_params()`.
        """
        super().__init__()
        
        self.fn_name = 'Holland et al (2020) Multi-factor Mg/Ca Equation'
        self.fn_text = 'mgca = mgca_sw**A * B exp((C1 * ca_sw + C2 * carb_sw + D) * temperature)'

        # update class attributes for exponential case
        self.variables.update(['mgca', 'temperature', 'carb_sw', 'ca_sw', 'mgca_sw'])
        
        self._var_update(mgca=mgca, temperature=temperature, carb_sw=carb_sw, ca_sw=ca_sw, mgca_sw=mgca_sw)
        
        if parameters is None:
            parameters = 'Multispecies_Anand'
        if isinstance(parameters, str):
            self.parameters = params.load(proxy='mgca', mode='holland_2020', parameters=parameters)
        elif isinstance(parameters, foramgeochem.general.params):
            self.parameters = parameters
        else:
            try:
                self.parameters = params(values=parameters)
            except:
                raise ValueError('`parameters` must be a string, a <foramgeochem.general.params> object or array_like')
    
        self._calc_mgca = holland2020_calc_mgca
        self._calc_temp = holland2020_calc_temp
        self._calc_carb = holland2020_calc_carb_sw
        self._calc_mgca_sw = holland2020_calc_mgca_sw
        self._calc_ca_sw = holland2020_calc_Ca_sw

    def __repr__(self):
        outstr = []
        outstr.append(self.fn_name)
        outstr.append('-' * len(self.fn_name))
        outstr.append(self.fn_text + '\n')
        outstr.append(self.parameters.__repr__())
        outstr.append('\nVariables:')
        outstr.append('  Accepted: {}'.format(self.variables))
        outstr.append('  Provided: {}'.format(self.variables.difference(self.missing)))

        return '\n'.join(outstr)

    @staticmethod
    def available_params():
        """
        List defined parameter sets and associated info.
        """
        params.available_parameters(proxy='mgca', mode='holland_2020')

    def calc_mgca(self, temperature=None, carb_sw=None, ca_sw=None, mgca_sw=None):
        """
        Calculate foram Mg/Ca from Temperature, and seawater Mg/Ca, [Ca] and carbon chemistry.

        Parameters
        ----------
        temperature : array_like
            Temperature in celcius.
        mgca_sw : array_like
            Seawater Mg/Ca, in mol/mol.
        ca_sw : array_like
            Seawater calcium concentration, in mol kg-1.
        carb_sw : array_like
            Seawater carbon parameter - either DIC, CO3 or pH, depending on the species.
        p : array_like
            Parameters for the model, in the order [A, B, C, D, E].

        Returns
        -------
        Foraminiferal Mg/Ca in mmol/mol : array_like
        """
        self._var_update(temperature=temperature, carb_sw=carb_sw, ca_sw=ca_sw, mgca_sw=mgca_sw)
        self._var_check()
        
        req = self.missing.intersection(['temperature', 'carb_sw', 'ca_sw', 'mgca_sw'])
        if len(req) > 0:
            raise ValueError('Please provide {}'.format(req))
        
        self.last_calc = self._calc_mgca(temp=self.temperature, mgca_sw=self.mgca_sw, carb_sw=self.carb_sw, ca_sw=self.ca_sw,
                                           p=self.parameters.values)
        
        return self.last_calc

    def calc_temp(self, mgca=None, carb_sw=None, ca_sw=None, mgca_sw=None):
        """
        Calculate Temperature from foram Mg/Ca, and seawater Mg/Ca, [Ca] and carbon chemistry.

        Parameters
        ----------
        mgca : array_like
            Foraminiferal Mg/Ca in mmol/mol
        mgca_sw : array_like
            Seawater Mg/Ca, in mol/mol.
        ca_sw : array_like
            Seawater calcium concentration, in mol kg-1.
        carb_sw : array_like
            Seawater carbon parameter - either DIC, CO3 or pH, depending on the species.
        p : array_like
            Parameters for the model, in the order [A, B, C, D, E].

        Returns
        -------
        Temperature in celcius. : array_like
        """
        self._var_update(mgca=mgca, carb_sw=carb_sw, ca_sw=ca_sw, mgca_sw=mgca_sw)
        self._var_check()
        
        req = self.missing.intersection(['mgca', 'carb_sw', 'ca_sw', 'mgca_sw'])
        if len(req) > 0:
            raise ValueError('Please provide {}'.format(req))
            
        self.last_calc = self._calc_temp(mgca=self.mgca, carb_sw=self.carb_sw, mgca_sw=self.mgca_sw, ca_sw=self.ca_sw,
                                         p=self.parameters.values)
            
        return self.last_calc
        
    def calc_carb(self, mgca=None, temperature=None, ca_sw=None, mgca_sw=None):
        """
        Calculate seawater carbon chemistry from Temperature, foram Mg/Ca, and seawater Mg/Ca and [Ca].

        Parameters
        ----------
        temperature : array_like
            Temperature in celcius.
        mgca : array_like
            Foraminiferal Mg/Ca in mmol/mol
        mgca_sw : array_like
            Seawater Mg/Ca, in mol/mol.
        ca_sw : array_like
            Seawater calcium concentration, in mol kg-1.        
        p : array_like
            Parameters for the model, in the order [A, B, C, D, E].

        Returns
        -------
        Seawater carbon parameter - either DIC, CO3 or pH, depending on the species. : array_like
        """
        self._var_update(mgca=mgca, temperature=temperature, ca_sw=ca_sw, mgca_sw=mgca_sw)
        self._var_check()
        
        req = self.missing.intersection(['mgca', 'temperature', 'ca_sw', 'mgca_sw'])
        if len(req) > 0:
            raise ValueError('Please provide {}'.format(req))
            
        self.last_calc = self._calc_carb(mgca=self.mgca, temperature=self.temperature, mgca_sw=self.mgca_sw, ca_sw=self.ca_sw,
                                         p=self.parameters.values)
        
        return self.last_calc
    
    def calc_ca(self, mgca=None, temperature=None, carb_sw=None, mgca_sw=None):
        """
        Calculate seawater calcium concentration from Temperature, foram Mg/Ca, and seawater Mg/Ca and carbon chemistry.

        Parameters
        ----------
        temperature : array_like
            Temperature in celcius.
        mgca : array_like
            Foraminiferal Mg/Ca in mmol/mol
        mgca_sw : array_like
            Seawater Mg/Ca, in mol/mol.
        carb_sw : array_like
            Seawater carbon parameter - either DIC, CO3 or pH, depending on the species.
        p : array_like
            Parameters for the model, in the order [A, B, C, D, E].

        Returns
        -------
        Seawater calcium concentration, in mol kg-1. : array_like
        """
        self._var_update(mgca=mgca, temperature=temperature, carb_sw=carb_sw, mgca_sw=mgca_sw)
        self._var_check()
        
        req = self.missing.intersection(['mgca', 'temperature', 'carb_sw', 'mgca_sw'])
        if len(req) > 0:
            raise ValueError('Please provide {}'.format(req))
            
        self.last_calc = self._calc_ca_sw(mgca=self.mgca, temperature=self.temperature, mgca_sw=self.mgca_sw, carb_sw=self.carb_sw,
                                          p=self.parameters.values)
        
        return self.last_calc
    
    def calc_mgca_sw(self, mgca=None, temperature=None, ca_sw=None, carb_sw=None):
        """
        Calculate seawater Mg/Ca from Temperature, foram Mg/Ca, and seawater [Ca] and carbon chemistry.

        Parameters
        ----------
        temperature : array_like
            Temperature in celcius.
        mgca : array_like
            Foraminiferal Mg/Ca in mmol/mol
        ca_sw : array_like
            Seawater calcium concentration, in mol kg-1.
        carb_sw : array_like
            Seawater carbon parameter - either DIC, CO3 or pH, depending on the species.
        p : array_like
            Parameters for the model, in the order [A, B, C, D, E].

        Returns
        -------
        Seawater Mg/Ca, in mol/mol. : array_like
        """
        self._var_update(mgca=mgca, temperature=temperature, ca_sw=ca_sw, carb_sw=carb)
        self._var_check()
        
        req = self.missing.intersection(['mgca', 'temperature', 'ca_sw', 'carb_sw'])
        if len(req) > 0:
            raise ValueError('Please provide {}'.format(req))
            
        self.last_calc = self._calc_mgca_sw(mgca_sw=self.mgca, temperature=self.temperature, carb_sw=self.carb_sw, ca_sw=self.ca_sw,
                                            p=self.parameters.values)
        
        return self.last_calc