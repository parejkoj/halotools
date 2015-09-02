# -*- coding: utf-8 -*-
"""
Module composes the behavior of the profile models 
and the velocity models to produce models for the 
full phase space distribution of galaxies within their halos. 
"""

__author__ = ['Andrew Hearin']

from .profile_models import *
from .velocity_models import *
from .monte_carlo_phase_space import *
from . import model_defaults

class NFWPhaseSpace(NFWProfile, NFWJeansVelocity, MonteCarloGalProf):
    """ NFW halo profile, based on Navarro, Frenk and White (1999).

    """

    def __init__(self, **kwargs):
        """
        Parameters 
        ----------
        conc_mass_model : string, optional  
            Specifies the calibrated fitting function used to model the concentration-mass relation. 
             Default is set in `~halotools.empirical_models.sim_defaults`.

        cosmology : object, optional 
            Astropy cosmology object. Default is set in `~halotools.empirical_models.sim_defaults`.

        redshift : float, optional  
            Default is set in `~halotools.empirical_models.sim_defaults`.

        halo_boundary : string, optional  
            String giving the column name of the halo catalog that stores the boundary of the halo. 
            Default is set in the `~halotools.empirical_models.model_defaults` module. 

        """

        super(NFWPhaseSpace, self).__init__(**kwargs)

        cmin, cmax, dc = (
            model_defaults.min_permitted_conc, 
            model_defaults.max_permitted_conc,
            model_defaults.default_dconc
            )
        self._setup_lookup_tables((cmin, cmax, dc))

        self.column_keys_to_allocate = np.dtype([
            ('host_centric_distance', 'f8'), 
            ('x', 'f8'), ('y', 'f8'), ('z', 'f8'), 
            ('halo_x', 'f8'), ('halo_y', 'f8'), ('halo_z', 'f8'), 
            ('vx', 'f8'), ('vy', 'f8'), ('vz', 'f8'), 
            ('halo_vx', 'f8'), ('halo_vy', 'f8'), ('halo_vz', 'f8'), 
            ])

    def assign_phase_space(self, halo_table):
        """
        """
        self.mc_pos(halo_table = halo_table)
        self.mc_vel(halo_table = halo_table)












        