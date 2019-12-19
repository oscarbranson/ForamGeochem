.. _mgca:

=====
Mg/Ca
=====

Functions for converting between foraminiferal Mg/Ca and environmental parameters.

Exponential
-----------
The 'classic' exponential function, where Mg/Ca depends only on temperature, following:

.. math ::

    Mg/Ca_{foram} = A\ e^{Temp\ B}

Functions
^^^^^^^^^

Converter:
 - :class:`foramgeochem.mgca.exponential`

Transfer functions:
 - :func:`~foramgeochem.mgca.tfr.exp_mgca_2_temp`
 - :func:`~foramgeochem.mgca.tfr.exp_temp_2_mgca` 

Holland et al (sub.)
--------------------
The multi-factor formulation of Holland et al (sub.), which accounts for the influence of Temperature, Mg/Ca\ :sub:`SW`, [Ca]\ :sub:`SW` and a species-specific carbon-system parameter (*carb*, i.e. DIC\ :sub:`SW`, [CO\ :sub:`3`]\ :sub:`SW` or pH).

.. math ::

    Mg/Ca_{foram} = Mg/Ca_{SW}\ ^A\ DIC_{SW}\ ^B\ e^{(C\ [Ca]_{SW} + D T + E)}

Functions
^^^^^^^^^

Converter:
 - :class:`foramgeochem.mgca.Holland`

Transfer functions:
 - :func:`~foramgeochem.mgca.tfr.holland2018_calc_mgca`
 - :func:`~foramgeochem.mgca.tfr.holland2018_calc_temp`
 - :func:`~foramgeochem.mgca.tfr.holland2018_calc_mgca_sw`
 - :func:`~foramgeochem.mgca.tfr.holland2018_calc_Ca_sw`
 - :func:`~foramgeochem.mgca.tfr.holland2018_calc_carb_sw`
