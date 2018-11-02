import foramgeochem
from foramgeochem.general import proxy, params
from foramgeochem.mgca import tfr

class MgCa_exponential(proxy):
    """
    The 'classic' exponential relationship between formainiferal Mg/Ca and temperature.
    """
    def __init__(self, mgca_f=None, temperature=None, parameters=None):
        """
        The 'classic' exponential relationship between formainiferal Mg/Ca and temperature.

        mgca_f = A * exp(temperature * B)


        Parameters
        ----------
        mgca_f : float or array_like
            The Mg/Ca of foraminiferal calcite, in mmol/mol.
        temperature : float or array_like
            The temperature, in degrees celcius.
        parameters : array_like, str or `params` object
            Either a 'params' object containing parameter values and associated unctertainties,
            a string selecting one of the in-built options, or an array_like of numbers
            to use as parameters.
        """
        super().__init__()
        
        # update class attributes for exponential case
        self.variables.update(['mgca_f', 'temperature'])
        
        self._var_update(mgca_f=mgca_f, temperature=temperature)
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
        self._calc_mgca_f = tfr.exp_temp_2_mgca
    
    def __repr__(self):
        outstr = []
        outstr.append('Exponential Mg/Ca-Temperature Relationship')
        outstr.append('------------------------------------------')
        outstr.append('Parameter Info: ' + self.parameters.info)
        outstr.append(' - Variables Accepted: {}'.format(self.variables))
        outstr.append(' - Variables Provided: {}'.format(self.variables.difference(self.missing)))

        return '\n'.join(outstr)
        
    def calc_temp(self, mgca_f=None):
        """
        Calculate temperature from Mg/Ca.
        """
        self._var_update(mgca_f=mgca_f)
        self._var_check()
        
        if 'mgca_f' not in self.missing:
            self.last_calc = self._calc_temp(self.mgca_f, *self.parameters.values)
        else:
            raise ValueError('Please provide `mgca_f`')
            
        return self.last_calc

    def calc_mgca(self, temperature=None):
        """
        Calculate Mg/Ca from temperature
        """
        self._var_update(temperature=temperature)
        self._var_check()
        
        if 'temperature' not in self.missing:
            self.last_calc = self._calc_mgca_f(self.temperature, *self.parameters.values)
        else:
            raise ValueError('Please provide `mgca_f`')
            
        return self.last_calc


class MgCa_Holland(proxy):
    """
    The multi-factor modified exponential relationship between Mg/Ca, temperature, carbon, [Ca] and Mg/Casw of Holland et al (2018)
    """
    def __init__(self, mgca_f=None, temperature=None, carb=2100e-6, ca=10.2e-3, mgca_sw=5.0, parameters=None):
        """

        """
        super().__init__()
        
        # update class attributes for exponential case
        self.variables.update(['mgca_f', 'temperature', 'carb', 'ca', 'mgca_sw'])
        
        self._var_update(mgca_f=mgca_f, temperature=temperature, carb=carb, ca=ca, mgca_sw=mgca_sw)
        
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
    
        self._calc_mgca_f = tfr.holland2018_calc_mgca
        self._calc_temp = tfr.holland2018_calc_temp
        self._calc_carb = tfr.holland2018_calc_dic_sw
        self._calc_mgca_sw = tfr.holland2018_calc_mgca_sw
        self._calc_ca_sw = tfr.holland2018_calc_Ca_sw

    def __repr__(self):
        outstr = []
        outstr.append('Holland et al (2018) Multi-factor Mg/Ca Equation')
        outstr.append('------------------------------------------------')
        outstr.append('Parameter Info: ' + self.parameters.info)
        outstr.append(' - Variables Accepted: {}'.format(self.variables))
        outstr.append(' - Variables Provided: {}'.format(self.variables.difference(self.missing)))

        return '\n'.join(outstr)


    def calc_mgca_f(self, temperature=None, carb=None, ca=None, mgca_sw=None):
        self._var_update(temperature=temperature, carb=carb, ca=ca, mgca_sw=mgca_sw)
        self._var_check()
        
        req = self.missing.intersection(['temperature', 'carb', 'ca', 'mgca_sw'])
        if len(req) > 0:
            raise ValueError('Please provide {}'.format(req))
        
        self.last_calc = self._calc_mgca_f(temp=self.temperature, mgca_sw=self.mgca_sw, dic_sw=self.carb, ca_sw=self.ca,
                                           p=self.parameters.values)
        
        return self.last_calc

    def calc_temp(self, mgca_f=None, carb=None, ca=None, mgca_sw=None):
        self._var_update(mgca_f=mgca_f, carb=carb, ca=ca, mgca_sw=mgca_sw)
        self._var_check()
        
        req = self.missing.intersection(['mgca_f', 'carb', 'ca', 'mgca_sw'])
        if len(req) > 0:
            raise ValueError('Please provide {}'.format(req))
            
        self.last_calc = self._calc_temp(mgca=self.mgca_f, dic_sw=self.carb, mgca_sw=self.mgca_sw, ca_sw=self.ca,
                                         p=self.parameters.values)
            
        return self.last_calc
        
    def calc_carb(self, mgca_f=None, temperature=None, ca=None, mgca_sw=None):
        self._var_update(mgca_f=mgca_f, temperature=temperature, ca=ca, mgca_sw=mgca_sw)
        self._var_check()
        
        req = self.missing.intersection(['mgca_f', 'temperature', 'ca', 'mgca_sw'])
        if len(req) > 0:
            raise ValueError('Please provide {}'.format(req))
            
        self.last_calc = self._calc_carb(mgca=self.mgca_f, temperature=self.temperature, mgca_sw=self.mgca_sw, ca_sw=self.ca,
                                         p=self.parameters.values)
        
        return self.last_calc
    
    def calc_ca(self, mgca_f=None, temperature=None, carb=None, mgca_sw=None):
        self._var_update(mgca_f=mgca_f, temperature=temperature, carb=carb, mgca_sw=mgca_sw)
        self._var_check()
        
        req = self.missing.intersection(['mgca_f', 'temperature', 'carb', 'mgca_sw'])
        if len(req) > 0:
            raise ValueError('Please provide {}'.format(req))
            
        self.last_calc = self._calc_ca_sw(mgca=self.mgca_f, temperature=self.temperature, mgca_sw=self.mgca_sw, dic_sw=self.carb,
                                          p=self.parameters.values)
        
        return self.last_calc
    
    def calc_mgca_sw(self, mgca_f=None, temperature=None, ca=None, carb=None):
        self._var_update(mgca_f=mgca_f, temperature=temperature, ca=ca, dic_sw=carb)
        self._var_check()
        
        req = self.missing.intersection(['mgca_f', 'temperature', 'ca', 'carb'])
        if len(req) > 0:
            raise ValueError('Please provide {}'.format(req))
            
        self.last_calc = self._calc_mgca_sw(mgca=self.mgca_f, temperature=self.temperature, dic_sw=self.carb, ca_sw=self.ca,
                                            p=self.parameters.values)
        
        return self.last_calc