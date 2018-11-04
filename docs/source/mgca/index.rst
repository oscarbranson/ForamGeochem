=====
Mg/Ca
=====

Functions for converting between foraminiferal Mg/Ca and environmental parameters.

Exponential
-----------
The 'classic' exponential function, where Mg/Ca depends only on temperature, following:

.. math ::

    Mg/Ca_{foram} = A\ e^{Temp\ B}

This function is accessible via the 

Holland et al (2018)
--------------------
The multi-factor formulation of Holland et al (2018), which accounts for the influence of Temperature, Mg/Ca\ :sub:`SW`, [Ca]\ :sub:`SW` and a species-specific carbon-system parameter (*carb*, i.e. DIC\ :sub:`SW`, [CO\ :sub:`3`]\ :sub:`SW` or pH).

.. math ::

    Mg/Ca_{foram} = Mg/Ca_{SW}\ ^A\ B\ e^{(C_1 [Ca]_{SW} + C_2 [carb] + D) Temp}