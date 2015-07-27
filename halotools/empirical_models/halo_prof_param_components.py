# -*- coding: utf-8 -*-
"""

halotools.halo_prof_param_components contains the classes and functions 
defining relations between dark matter halos and the parameters 
governing their internal structure. The classic example of such a component 
is a relation between halo mass and NFW concentration. More generally, 
a halo profile parameter model is just a mapping between *any* 
halo property and *any* profile parameter. 

"""

import numpy as np
import model_defaults
from ..sim_manager import sim_defaults
import model_defaults

__all__ = ['ConcMass']

class ConcMass(object):
    """ Container class for commonly used concentration-mass 
    relations in the literature. 

    The only currently supported model is `dutton_maccio14`.

    """

    def __init__(self, cosmology=sim_defaults.default_cosmology, 
        redshift = sim_defaults.default_redshift, 
        prim_haloprop_key = model_defaults.prim_haloprop_key, 
        conc_mass_model = model_defaults.conc_mass_model, **kwargs):
        """
        Parameters 
        ----------
        cosmology : object, optional keyword argument
            Astropy cosmology object. Default is set in `~halotools.empirical_models.sim_defaults`.

        redshift : float, optional keyword argument 
            Default is set in `~halotools.empirical_models.sim_defaults`.

        prim_haloprop_key : string, optional keyword argument 
            Specifies the column name of the mass-like halo property, e.g., 'mvir' or 'm200b'. 
            Default is set in `~halotools.empirical_models.sim_defaults`.

        conc_mass_model : string, optional keyword argument 
            Specifies the calibrated fitting function used to model the concentration-mass relation. 
             Default is set in `~halotools.empirical_models.sim_defaults`.

        Examples 
        ---------
        >>> conc_mass_model = ConcMass()
        >>> conc_mass_model = ConcMass(redshift = 2, prim_haloprop_key = 'm500c')

        """
        self.cosmology = cosmology
        self.redshift = redshift
        self.prim_haloprop_key = prim_haloprop_key
        self.conc_mass_model = conc_mass_model

    def __call__(self, **kwargs):
        """ Method used to evaluate the mean NFW concentration as a function of 
        halo mass. 

        This is the primary method seen by the outside world. It has no functionality 
        of its own, it only calls the desired function based on the model keyword. 

        Parameters
        ----------        
        prim_haloprop : array, optional keyword argument 
            Array of mass-like variable upon which occupation statistics are based. 
            If ``prim_haloprop`` is not passed, then ``halo_table`` keyword argument must be passed. 

        halo_table : object, optional keyword argument 
            Data table storing halo catalog. 
            If ``halo_table`` is not passed, then ``prim_haloprop`` keyword argument must be passed. 

        Notes 
        -----
        The testing for this model can be found in 
        `~halotools.empirical_models.test_empirical_models.test_halo_prof_param_components`. 

        """
        # Retrieve the array storing the mass-like variable
        if 'halo_table' in kwargs.keys():
            mass = kwargs['halo_table'][self.prim_haloprop_key]
        elif 'prim_haloprop' in kwargs.keys():
            mass = kwargs['prim_haloprop']
        else:
            raise KeyError("Must pass one of the following keyword arguments to mean_occupation:\n"
                "``halo_table`` or ``prim_haloprop``")

        return getattr(self, self.conc_mass_model)(mass)

    def dutton_maccio14(self, mass):
        """ Power-law fit to the concentration-mass relation from 
        Equations 12 & 13 of Dutton & Maccio 2014, arXiv:1402.7073.

        Parameters 
        ----------
        mass : array_like 

        Returns 
        -------
        c : array_like
            Concentrations of the input halos. 

        Notes 
        -----        
        This model was only calibrated for the Planck 1-year cosmology.

        Model assumes that halo mass definition is Mvir.

        :math:`a = 0.537 + (1.025 - 0.537)\\exp(-0.718z^{1.08})`

        :math:`b = -0.097 + 0.024z`

        :math:`M_{0} = 10^{12}M_{\odot}/h`

        :math:`\\log_{10}c(M) = a + b\\log_{10}(M / M_{0})`

        """

        a = 0.537 + (1.025 - 0.537) * np.exp(-0.718 * self.redshift**1.08)
        b = -0.097 + 0.024 * self.redshift
        m0 = 1.e12

        logc = a + b * np.log10(mass / m0)
        c = 10**logc

        return c
