"""
Exponential relationship between Mg/Ca and temperature.
"""
import foramgeochem
from foramgeochem.general import proxy, params
from foramgeochem.helpers import load_params
from uncertainties.unumpy import log, exp

# exponential, as in traditional palaeothermometer
def exp_mgca_2_temp(mgca_f, p=None):
    """
    Calculate foram Mg/Ca from temperature using the 'classic' exponential function.

    Parameters
    ----------
    mgca_f : array_like
        Foraminiferal Mg/Ca in mmol/mol
    p : array_like
        Parameters for the model, in the order [A, B].
    
    Returns
    -------
    Temperature in Celcius. : array_like
    """
    A, B = p
    return log(mgca_f / A) / B 

def exp_temp_2_mgca(temp, p=None):
    """
    Calculate temperature from foram Mg/Ca using the 'classic' exponential function.

    Parameters
    ----------
    temp : array_like
        Temperature in Celcius.
    p : array_like
        Parameters for the model, in the order [A, B].
    
    Returns
    -------
    Foraminiferal Mg/Ca in mmol/mol : array_like
    """
    A, B = p
    return A * exp(temp * B)

class exponential(proxy):
    """
    The 'classic' exponential relationship between formainiferal Mg/Ca and temperature.
    """
    def __init__(self, mgca=None, temperature=None, parameters=None):
        """
        The 'classic' exponential relationship between formainiferal Mg/Ca and temperature.

        mgca = A * exp(temperature * B)


        Parameters
        ----------
        mgca : float or array_like
            The Mg/Ca of foraminiferal calcite, in mmol/mol.
        temperature : float or array_like
            The temperature, in degrees celcius.
        parameters : array_like, str or `params` object
            Either a 'params' object containing parameter values and associated unctertainties,
            a string selecting one of the in-built options, or an array_like of numbers
            to use as parameters.To see a list of pre-defined parameters, use:
               `Holland.available_params()`.
        """
        super().__init__()
        
        self.fn_name = 'Exponential Mg/Ca-Temperature Relationship'
        self.fn_text = 'mgca = A * exp(temperature * B)'

        # update class attributes for exponential case
        self.variables.update(['mgca', 'temperature'])
        
        self._var_update(mgca=mgca, temperature=temperature)
        self._var_check()

        if parameters is None:
            parameters = 'Multispecies_Anand'
        if isinstance(parameters, str):
            self.parameters = params.load(proxy='mgca', mode='exponential', parameters=parameters)
        elif isinstance(parameters, foramgeochem.general.params):
            self.parameters = parameters    
        else:
            try:
                self.parameters = params(values=parameters)
            except:
                raise ValueError('`parameters` must be a string, a <foramgeochem.general.params> object or array_like')
    
        self._calc_temp = exp_mgca_2_temp
        self._calc_mgca = exp_temp_2_mgca
    
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
        params.available_parameters(proxy='mgca', mode='exponential')
        
    def calc_temp(self, mgca=None):
        """
        Calculate temperature from Mg/Ca.
        """
        self._var_update(mgca=mgca)
        self._var_check()
        
        if 'mgca' not in self.missing:
            self.last_calc = self._calc_temp(self.mgca, self.parameters.values)
        else:
            raise ValueError('Please provide `mgca`')
            
        return self.last_calc

    def calc_mgca(self, temperature=None):
        """
        Calculate Mg/Ca from temperature
        """
        self._var_update(temperature=temperature)
        self._var_check()
        
        if 'temperature' not in self.missing:
            self.last_calc = self._calc_mgca(self.temperature, self.parameters.values)
        else:
            raise ValueError('Please provide `mgca`')
            
        return self.last_calc
