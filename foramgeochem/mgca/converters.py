import foramgeochem
from foramgeochem.general import proxy, params
from foramgeochem.mgca import tfr
from foramgeochem.helpers import load_params

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
    
        self._calc_temp = tfr.exp_mgca_2_temp
        self._calc_mgca = tfr.exp_temp_2_mgca
    
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


class Holland(proxy):
    """
    The multi-factor modified exponential relationship between Mg/Ca, temperature, carbon, [Ca] and Mg/Casw of Holland et al (sub.)
    """
    def __init__(self, mgca=None, temperature=None, carb_sw=2100e-6, ca_sw=10.2e-3, mgca_sw=5.0, parameters=None):
        """
        Create an object for converting between Mg/Ca and environmental parameters.

        Uses the formulation of Holland et al (sub.). Different parameter sets discussed
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
        
        self.fn_name = 'Holland et al (sub.) Multi-factor Mg/Ca Equation'
        self.fn_text = 'mgca = mgca_sw**A * B exp((C1 * ca_sw + C2 * carb_sw + D) * temperature)'

        # update class attributes for exponential case
        self.variables.update(['mgca', 'temperature', 'carb_sw', 'ca_sw', 'mgca_sw'])
        
        self._var_update(mgca=mgca, temperature=temperature, carb_sw=carb_sw, ca_sw=ca_sw, mgca_sw=mgca_sw)
        
        if parameters is None:
            parameters = 'Multispecies_Anand'
        if isinstance(parameters, str):
            self.parameters = params.load(proxy='mgca', mode='holland_2018', parameters=parameters)
        elif isinstance(parameters, foramgeochem.general.params):
            self.parameters = parameters
        else:
            try:
                self.parameters = params(values=parameters)
            except:
                raise ValueError('`parameters` must be a string, a <foramgeochem.general.params> object or array_like')
    
        self._calc_mgca = tfr.holland2020_calc_mgca
        self._calc_temp = tfr.holland2020_calc_temp
        self._calc_carb = tfr.holland2020_calc_carb_sw
        self._calc_mgca_sw = tfr.holland2020_calc_mgca_sw
        self._calc_ca_sw = tfr.holland2020_calc_Ca_sw

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
        params.available_parameters(proxy='mgca', mode='holland_2018')

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
            Parameters for the model, in the order [A, B, C1, C2, D].

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
            Parameters for the model, in the order [A, B, C1, C2, D].

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
            Parameters for the model, in the order [A, B, C1, C2, D].

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
            Parameters for the model, in the order [A, B, C1, C2, D].

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
            Parameters for the model, in the order [A, B, C1, C2, D].

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