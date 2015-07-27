# -*- coding: utf-8 -*-
"""
This module contains various component features used by 
HOD-style models of the galaxy-halo connection. 

"""

__all__ = (['OccupationComponent','Zheng07Cens','Zheng07Sats', 
    'Leauthaud11Cens', 'Leauthaud11Sats']
    )

from functools import partial
from copy import copy
import numpy as np
import math
from scipy.special import erf 
from scipy.stats import poisson
from scipy.optimize import brentq
from scipy.interpolate import InterpolatedUnivariateSpline as spline

from . import model_defaults
from . import model_helpers as model_helpers
from . import smhm_components

from ..utils.array_utils import array_like_length as custom_len
from ..  import sim_manager
from ..halotools_exceptions import HalotoolsModelInputError

from astropy.extern import six
from abc import ABCMeta, abstractmethod, abstractproperty
import warnings

@six.add_metaclass(ABCMeta)
class OccupationComponent(model_helpers.GalPropModel):
    """ Abstract base class of any occupation model. 
    Functionality is mostly trivial. 
    The sole purpose of the base class is to 
    standardize the attributes and methods 
    required of any HOD-style model for halo occupation statistics. 
    """
    def __init__(self, **kwargs):
        """ 
        Parameters 
        ----------
        gal_type : string, keyword argument 
            Name of the galaxy population whose occupation statistics is being modeled. 

        threshold : float, keyword argument
            Threshold value defining the selection function of the galaxy population 
            being modeled. Typically refers to absolute magnitude or stellar mass. 

        occupation_bound : float, keyword argument
            Upper bound on the number of gal_type galaxies per halo. 
            The only currently supported values are unity or infinity. 

        prim_haloprop_key : string, keyword argument 
            String giving the column name of the primary halo property governing 
            the occupation statistics of gal_type galaxies, e.g., ``halo_mvir``. 

        sec_haloprop_key : string, optional keyword argument
            String giving the column name of the secondary halo property governing 
            the occupation statistics of gal_type galaxies, e.g., ``halo_nfw_conc``.
            Only pertains to galaxy populations with assembly-biased occupations. 
            Default is None. 

        """
        super(OccupationComponent, self).__init__(galprop_key='occupation')

        required_kwargs = ['gal_type',  'threshold', 'occupation_bound', 'prim_haloprop_key']
        model_helpers.bind_required_kwargs(required_kwargs, self, **kwargs)

        if 'sec_haloprop_key' in kwargs.keys():
            self.sec_haloprop_key = kwargs['sec_haloprop_key']

        self.param_dict = {}

        # Enforce the requirement that sub-classes have been configured properly
        required_method_name = 'mean_occupation'
        if not hasattr(self, required_method_name):
            raise SyntaxError("Any sub-class of OccupationComponent must "
                "implement a method named %s " % required_method_name)

        self._additional_methods_to_inherit = ['mean_occupation']

    def mc_occupation(self, seed=None, **kwargs):
        """ Method to generate Monte Carlo realizations of the abundance of galaxies. 

        Parameters
        ----------        
        prim_haloprop : array, optional keyword argument 
            Array of mass-like variable upon which occupation statistics are based. 
            If ``prim_haloprop`` is not passed, then ``halo_table`` keyword argument must be passed. 

        halo_table : object, optional keyword argument 
            Data table storing halo catalog. 
            If ``halo_table`` is not passed, then ``prim_haloprop`` keyword argument must be passed. 

        seed : int, optional keyword argument 
            Random number seed used to generate the Monte Carlo realization. 
            Default is None. 

        Returns
        -------
        mc_abundance : array
            Integer array giving the number of galaxies in each of the input halo_table.     
        """ 
        first_occupation_moment = self.mean_occupation(**kwargs)
        if self.occupation_bound == 1:
            return self._nearest_integer_distribution(first_occupation_moment, seed=seed, **kwargs)
        elif self.occupation_bound == float("inf"):
            return self._poisson_distribution(first_occupation_moment, seed=seed, **kwargs)
        else:
            raise KeyError("The only permissible values of occupation_bound for instances "
                "of OccupationComponent are unity and infinity.")

    def _nearest_integer_distribution(self, first_occupation_moment, seed=None, **kwargs):
        """ Nearest-integer distribution used to draw Monte Carlo occupation statistics 
        for central-like populations with only permissible galaxy per halo.

        Parameters 
        ----------
        first_occupation_moment : array
            Array giving the first moment of the occupation distribution function. 

        seed : int, optional keyword argument 
            Random number seed used to generate the Monte Carlo realization. 
            Default is None. 

        Returns
        -------
        mc_abundance : array
            Integer array giving the number of galaxies in each of the input halo_table. 
        """
        np.random.seed(seed=seed)
        mc_generator = np.random.random(custom_len(first_occupation_moment))
        return np.where(mc_generator < first_occupation_moment, 1, 0)

    def _poisson_distribution(self, first_occupation_moment, seed=None, **kwargs):
        """ Poisson distribution used to draw Monte Carlo occupation statistics 
        for satellite-like populations in which per-halo abundances are unbounded. 

        Parameters 
        ----------
        first_occupation_moment : array
            Array giving the first moment of the occupation distribution function. 

        seed : int, optional keyword argument 
            Random number seed used to generate the Monte Carlo realization. 
            Default is None. 

        Returns
        -------
        mc_abundance : array
            Integer array giving the number of galaxies in each of the input halo_table. 
        """
        np.random.seed(seed=seed)
        # The scipy built-in Poisson number generator raises an exception 
        # if its input is zero, so here we impose a simple workaround
        first_occupation_moment = np.where(first_occupation_moment <=0, 
            model_defaults.default_tiny_poisson_fluctuation, first_occupation_moment)
        return poisson.rvs(first_occupation_moment)

class Zheng07Cens(OccupationComponent):
    """ ``Erf`` function model for the occupation statistics of central galaxies, 
    introduced in Zheng et al. 2005, arXiv:0408564. This implementation uses 
    Zheng et al. 2007, arXiv:0703457, to assign fiducial parameter values. 

    """

    def __init__(self, 
        gal_type='centrals', 
        threshold=model_defaults.default_luminosity_threshold,
        prim_haloprop_key=model_defaults.prim_haloprop_key,
        **kwargs):
        """
        Parameters 
        ----------
        gal_type : string, optional keyword argument
            Name of the galaxy population being modeled. Default is ``centrals``.  

        threshold : float, optional keyword argument
            Luminosity threshold of the mock galaxy sample. If specified, 
            input value must agree with one of the thresholds used in Zheng07 to fit HODs: 
            [-18, -18.5, -19, -19.5, -20, -20.5, -21, -21.5, -22].
            Default value is specified in the `~halotools.empirical_models.model_defaults` module.

        prim_haloprop_key : string, optional keyword argument 
            String giving the column name of the primary halo property governing 
            the occupation statistics of gal_type galaxies. 
            Default value is specified in the `~halotools.empirical_models.model_defaults` module.

        Examples 
        --------
        >>> cen_model = Zheng07Cens()
        >>> cen_model = Zheng07Cens(gal_type='cens', threshold=-19.5)
        >>> cen_model = Zheng07Cens(prim_haloprop_key='halo_m200b')

        Notes 
        -----
        The test suite for this model is documented at 
        `~halotools.empirical_models.test_empirical_models.test_Zheng07Cens`
        """
        occupation_bound = 1.0

        # Call the super class constructor, which binds all the 
        # arguments to the instance.  
        super(Zheng07Cens, self).__init__(gal_type=gal_type, 
            threshold=threshold, occupation_bound=occupation_bound, 
            prim_haloprop_key=prim_haloprop_key, 
            **kwargs)

        self.param_dict = self.get_published_parameters(self.threshold)

        self.publications = ['arXiv:0408564', 'arXiv:0703457']

    def mean_occupation(self, **kwargs):
        """ Expected number of central galaxies in a halo of mass halo_mass.
        See Equation 2 of arXiv:0703457.

        Parameters
        ----------        
        prim_haloprop : array, optional keyword argument 
            Array of mass-like variable upon which occupation statistics are based. 
            If ``prim_haloprop`` is not passed, then ``halo_table`` keyword argument must be passed. 

        halo_table : object, optional keyword argument 
            Data table storing halo catalog. 
            If ``halo_table`` is not passed, then ``prim_haloprop`` keyword argument must be passed. 

        Returns
        -------
        mean_ncen : array
            Mean number of central galaxies in the input halo_table. 

        Examples 
        --------
        >>> cen_model = Zheng07Cens()

        The `mean_occupation` method of all OccupationComponent instances supports 
        two different options for arguments. The first option is to directly 
        pass the array of the primary halo property: 

        >>> testmass = np.logspace(10, 15, num=50)
        >>> mean_ncen = cen_model.mean_occupation(prim_haloprop = testmass)

        The second option is to pass `mean_occupation` a full halo catalog. 
        In this case, the array storing the primary halo property will be selected 
        by accessing the ``cen_model.prim_haloprop_key`` column of the input halo catalog. 
        For illustration purposes, we'll use a fake halo catalog rather than a 
        (much larger) full one:

        >>> fake_sim = sim_manager.FakeSim()
        >>> mean_ncen = cen_model.mean_occupation(halo_table=fake_sim.halo_table)

        Notes 
        -----
        The `mean_occupation` method computes the following function: 

        :math:`\\langle N_{\\mathrm{cen}} \\rangle_{M} = 
        \\frac{1}{2}\\left( 1 + 
        \\mathrm{erf}\\left( \\frac{\\log_{10}M - 
        \\log_{10}M_{min}}{\\sigma_{\\log_{10}M}} \\right) \\right)`

        """
        # Retrieve the array storing the mass-like variable
        if 'halo_table' in kwargs.keys():
            mass = kwargs['halo_table'][self.prim_haloprop_key]
        elif 'prim_haloprop' in kwargs.keys():
            mass = kwargs['prim_haloprop']
        else:
            function_name = "Zheng07Cens.mean_occupation"
            raise HalotoolsModelInputError(function_name)

        logM = np.log10(mass)
        mean_ncen = 0.5*(1.0 + erf(
            (logM - self.param_dict['logMmin'])
            /self.param_dict['sigma_logM']))

        return mean_ncen


    def get_published_parameters(self, threshold, publication='Zheng07'):
        """
        Best-fit HOD parameters from Table 1 of Zheng et al. 2007.

        Parameters 
        ----------

        threshold : float
            Luminosity threshold defining the SDSS sample 
            to which Zheng et al. fit their HOD model. 
            If the ``publication`` keyword argument is set to ``Zheng07``, 
            then ``threshold`` must be agree with one of the published values: 
            [-18, -18.5, -19, -19.5, -20, -20.5, -21, -21.5, -22].

        publication : string, optional keyword argument 
            String specifying the publication that will be used to set  
            the values of ``param_dict``. Default is Zheng et al. (2007). 

        Returns 
        -------
        param_dict : dict
            Dictionary of model parameters whose values have been set to 
            agree with the values taken from Table 1 of Zheng et al. 2007.

        Examples 
        --------
        >>> cen_model = Zheng07Cens()
        >>> cen_model.param_dict = cen_model.get_published_parameters(cen_model.threshold)

        """

        def get_zheng07_params(threshold):
            #Load tabulated data from Zheng et al. 2007, Table 1
            logMmin_array = [11.35,11.46,11.6,11.75,12.02,12.3,12.79,13.38,14.22]
            sigma_logM_array = [0.25,0.24,0.26,0.28,0.26,0.21,0.39,0.51,0.77]
            # define the luminosity thresholds corresponding to the above data
            threshold_array = np.arange(-22,-17.5,0.5)
            threshold_array = threshold_array[::-1]

            threshold_index = np.where(threshold_array==threshold)[0]
            if len(threshold_index)==1:
                param_dict = {
                'logMmin': logMmin_array[threshold_index[0]],
                'sigma_logM' : sigma_logM_array[threshold_index[0]]
                }
            else:
                raise ValueError("Input luminosity threshold "
                    "does not match any of the Table 1 values of "
                    "Zheng et al. 2007 (arXiv:0703457)")

            return param_dict

        if publication in ['zheng07', 'Zheng07', 'Zheng_etal07', 'zheng_etal07','zheng2007','Zheng2007']:
            param_dict = get_zheng07_params(threshold)
            return param_dict
        else:
            raise KeyError("For Zheng07Cens, only supported best-fit models are currently Zheng et al. 2007")


class Leauthaud11Cens(OccupationComponent):
    """ HOD-style model for any central galaxy occupation that derives from 
    a stellar-to-halo-mass relation. 
    """
    def __init__(self, smhm_model=smhm_components.Moster13SmHm, 
        gal_type = 'centrals', 
        threshold = model_defaults.default_stellar_mass_threshold, 
        prim_haloprop_key=model_defaults.prim_haloprop_key,
        **kwargs):
        """
        Parameters 
        ----------
        gal_type : string, optional keyword argument
            Name of the galaxy population being modeled. Default is ``centrals``.  

        threshold : float, optional keyword argument
            Stellar mass threshold of the mock galaxy sample. 
            Default value is specified in the `~halotools.empirical_models.model_defaults` module.

        prim_haloprop_key : string, optional keyword argument 
            String giving the column name of the primary halo property governing 
            the occupation statistics of gal_type galaxies. 
            Default value is specified in the `~halotools.empirical_models.model_defaults` module.

        smhm_model : object, optional keyword argument 
            Sub-class of `~halotools.empirical_models.smhm_components.PrimGalpropModel` governing 
            the stellar-to-halo-mass relation. Default is `Moster13SmHm`. 

        redshift : float, optional keyword argument 
            Redshift of the stellar-to-halo-mass relation. Default is 0. 

        Examples 
        --------
        >>> cen_model = Leauthaud11Cens()
        >>> cen_model = Leauthaud11Cens(threshold = 11.25)
        >>> cen_model = Leauthaud11Cens(prim_haloprop_key = 'halo_m200b')

        """
        occupation_bound = 1.0

        # Call the super class constructor, which binds all the 
        # arguments to the instance.  
        super(Leauthaud11Cens, self).__init__(
            gal_type=gal_type, threshold=threshold, 
            occupation_bound=occupation_bound, 
            prim_haloprop_key = prim_haloprop_key, 
            **kwargs)

        self.smhm_model = smhm_model(
            gal_type=gal_type, prim_haloprop_key = prim_haloprop_key, **kwargs)
        self.param_dict = self.smhm_model.param_dict

        self.publications = ['arXiv:1103.2077', 'arXiv:1104.0928']

    def mean_occupation(self, **kwargs):
        """ Expected number of central galaxies in a halo.
        See Equation 8 of arXiv:1103.2077.

        Parameters
        ----------        
        prim_haloprop : array, optional keyword argument 
            Array of mass-like variable upon which occupation statistics are based. 
            If ``prim_haloprop`` is not passed, then ``halo_table`` keyword argument must be passed. 

        halo_table : object, optional keyword argument 
            Data table storing halo catalog. 
            If ``halo_table`` is not passed, then ``prim_haloprop`` keyword argument must be passed. 

        Returns
        -------
        mean_ncen : array
            Mean number of central galaxies in the halo of the input mass. 

        Notes 
        -----
        Assumes constant scatter in the stellar-to-halo-mass relation. 
        """
        logmstar = np.log10(self.smhm_model.mean_stellar_mass(**kwargs))
        logscatter = math.sqrt(2)*self.smhm_model.mean_scatter(**kwargs)

        mean_ncen = 0.5*(1.0 - 
            erf((self.threshold - logmstar)/logscatter))

        return mean_ncen        


class Zheng07Sats(OccupationComponent):
    """ Power law model for the occupation statistics of satellite galaxies, 
    introduced in Kravtsov et al. 2004, arXiv:0308519.

    :math:`\\langle N_{sat} \\rangle_{M} = \left( \\frac{M - M_{0}}{M_{1}} \\right)^{\\alpha}`

    """

    def __init__(self,
        gal_type='satellites', 
        threshold=model_defaults.default_luminosity_threshold,
        prim_haloprop_key=model_defaults.prim_haloprop_key,
        modulate_with_cenocc = False, 
        **kwargs):
        """
        Parameters 
        ----------
        gal_type : string, optional keyword argument
            Name of the galaxy population being modeled. Default is ``satellites``.  

        threshold : float, optional keyword argument
            Luminosity threshold of the mock galaxy sample. If specified, 
            input value must agree with one of the thresholds used in Zheng07 to fit HODs: 
            [-18, -18.5, -19, -19.5, -20, -20.5, -21, -21.5, -22].
            Default value is specified in the `~halotools.empirical_models.model_defaults` module.

        prim_haloprop_key : string, optional keyword argument 
            String giving the column name of the primary halo property governing 
            the occupation statistics of gal_type galaxies. 
            Default value is specified in the `~halotools.empirical_models.model_defaults` module.

        modulate_with_cenocc : bool, optional keyword argument 
            If True, the first satellite moment will be multiplied by the 
            the first central moment. Default is False. 
            If ``modulate_with_cenocc`` is True, 
            the mean occupation method of `Zheng07Sats` will 
            be multiplied by the the first moment of `Zheng07Cens`, 
            as in Zheng et al. 2007, so that:

            :math:`\\langle N_{\mathrm{sat}}\\rangle_{M}\\Rightarrow\\langle N_{\mathrm{sat}}\\rangle_{M}\\times\\langle N_{\mathrm{cen}}\\rangle_{M}`

        gal_type_centrals : string, optional keyword argument
            Name of the central galaxy population whose first moment 
            modulates the `mean_occupation` method of `Zheng07Sats`. 
            Only required if ``modulate_with_cenocc`` is True. 

        Examples 
        --------
        >>> sat_model = Zheng07Sats()
        >>> sat_model = Zheng07Sats(gal_type='sats')
        >>> sat_model = Zheng07Sats(threshold = -21)

        The ``param_dict`` attribute can be used to build an alternate 
        model from an existing instance. This feature has a variety of uses. For example, 
        suppose you wish to study how the choice of halo mass definition impacts HOD predictions:

        >>> sat_model1 = Zheng07Sats(threshold = -19.5, prim_haloprop_key='m200b')
        >>> sat_model1.param_dict['alpha_satellites'] = 1.05
        >>> sat_model2 = Zheng07Sats(threshold = -19.5, prim_haloprop_key='m500c')
        >>> sat_model2.param_dict = sat_model1.param_dict 

        After executing the above four lines of code, ``sat_model1`` and ``sat_model2`` are 
        identical in every respect, excepting only for the difference in the halo mass definition. 

        A common convention in HOD modeling of satellite populations is for the first 
        occupation moment of the satellites to be multiplied by the first occupation 
        moment of the associated central population. 
        The ``modulate_with_cenocc`` keyword arguments allows you 
        to study the impact of this choice:

        >>> sat_model1 = Zheng07Sats(threshold=-18)
        >>> cen_model_instance = Zheng07Cens(threshold = sat_model1.threshold)
        >>> sat_model2 = Zheng07Sats(threshold = sat_model1.threshold, modulate_with_cenocc=True, gal_type_centrals = 'cens')

        Now ``sat_model1`` and ``sat_model2`` are identical in every respect, 
        excepting only the following difference:

        :math:`\\langle N_{\mathrm{sat}}\\rangle^{\mathrm{model 2}} = \\langle N_{\mathrm{cen}}\\rangle\\times\\langle N_{\mathrm{sat}}\\rangle^{\mathrm{model 1}}`


        Notes 
        -----
        The test suite for this model is documented at 
        `~halotools.empirical_models.test_empirical_models.test_Zheng07Sats`

        """
        occupation_bound = float("inf")

        # Call the super class constructor, which binds all the 
        # arguments to the instance.  
        super(Zheng07Sats, self).__init__(
            gal_type=gal_type, threshold=threshold, 
            occupation_bound=occupation_bound, 
            prim_haloprop_key = prim_haloprop_key, 
            **kwargs)

        self.param_dict = self.get_published_parameters(self.threshold)

        self.modulate_with_cenocc = modulate_with_cenocc
        if self.modulate_with_cenocc is True:
            if 'gal_type_centrals' not in kwargs.keys():
                raise KeyError("If ``modulate_with_cenocc`` is True, must also pass "
                    "the gal_type_centrals keyword.")
            gal_type_centrals = kwargs['gal_type_centrals']
            self.central_occupation_model = Zheng07Cens(
                prim_haloprop_key = prim_haloprop_key, 
                threshold = threshold, gal_type = gal_type_centrals)
            self.ancillary_model_dependencies = ['central_occupation_model']
            

        self.publications = ['arXiv:0308519', 'arXiv:0703457']

    def mean_occupation(self, **kwargs):
        """Expected number of satellite galaxies in a halo of mass logM.
        See Equation 5 of arXiv:0703457.

        Parameters
        ----------        
        prim_haloprop : array, optional keyword argument
            Array storing a mass-like variable that governs the occupation statistics. 
            If ``prim_haloprop`` is not passed, then ``halo_table`` 
            keyword arguments must be passed. 

        halo_table : object, optional keyword argument 
            Data table storing halo catalog. 
            If ``halo_table`` is not passed, then ``prim_haloprop`` 
            keyword arguments must be passed. 

        Returns
        -------
        mean_nsat : float or array
            Mean number of satellite galaxies in a host halo of the specified mass. 

        :math:`\\langle N_{\\mathrm{sat}} \\rangle_{M} = \left( \\frac{M - M_{0}}{M_{1}} \\right)^{\\alpha} \\langle N_{\\mathrm{cen}} \\rangle_{M}`

        or 

        :math:`\\langle N_{\\mathrm{sat}} \\rangle_{M} = \left( \\frac{M - M_{0}}{M_{1}} \\right)^{\\alpha}`, 

        depending on whether a central model was passed to the constructor. 

        Examples 
        --------
        The `mean_occupation` method of all OccupationComponent instances supports 
        two different options for arguments. The first option is to directly 
        pass the array of the primary halo property: 

        >>> sat_model = Zheng07Sats()
        >>> testmass = np.logspace(10, 15, num=50)
        >>> mean_nsat = sat_model.mean_occupation(prim_haloprop =testmass)

        The second option is to pass `mean_occupation` a full halo catalog. 
        In this case, the array storing the primary halo property will be selected 
        by accessing the ``sat_model.prim_haloprop_key`` column of the input halo catalog. 
        For illustration purposes, we'll use a fake halo catalog rather than a 
        (much larger) full one:

        >>> fake_sim = sim_manager.FakeSim()
        >>> mean_nsat = sat_model.mean_occupation(halo_table=fake_sim.halo_table)

        """
        # Retrieve the array storing the mass-like variable
        if 'halo_table' in kwargs.keys():
            mass = kwargs['halo_table'][self.prim_haloprop_key]
        elif 'prim_haloprop' in kwargs.keys():
            mass = kwargs['prim_haloprop']
        else:
            function_name = "Zheng07Sats.mean_occupation"
            raise HalotoolsModelInputError(function_name)
        mass = np.array(mass)
        if np.shape(mass) == ():
            mass = np.array([mass])

        M0 = 10.**self.param_dict['logM0']
        M1 = 10.**self.param_dict['logM1']

        # Call to np.where raises a harmless RuntimeWarning exception if 
        # there are entries of input logM for which mean_nsat = 0
        # Evaluating mean_nsat using the catch_warnings context manager 
        # suppresses this warning
        mean_nsat = np.zeros_like(mass)

        idx_nonzero = np.where(mass - M0 > 0)[0]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)

            mean_nsat[idx_nonzero] = ((mass[idx_nonzero] - M0)/M1)**self.param_dict['alpha']

        # If a central occupation model was passed to the constructor, 
        # multiply mean_nsat by an overall factor of mean_ncen
        if self.modulate_with_cenocc is True:
            mean_ncen = self.central_occupation_model.mean_occupation(**kwargs)
            mean_nsat *= mean_ncen

        return mean_nsat

    def get_published_parameters(self, threshold, publication='Zheng07'):
        """
        Best-fit HOD parameters from Table 1 of Zheng et al. 2007.

        Parameters 
        ----------
        threshold : float
            Luminosity threshold of the mock galaxy sample. 
            Input value must agree with one of the thresholds used in Zheng07 to fit HODs: 
            [-18, -18.5, -19, -19.5, -20, -20.5, -21, -21.5, -22].

        publication : string, optional keyword argument 
            String specifying the publication that will be used to set  
            the values of ``param_dict``. Default is Zheng et al. (2007). 

        Returns 
        -------
        param_dict : dict
            Dictionary of model parameters whose values have been set to 
            the values taken from Table 1 of Zheng et al. 2007.

        Examples 
        --------
        >>> sat_model = Zheng07Sats()
        >>> sat_model.param_dict = sat_model.get_published_parameters(sat_model.threshold)
        """

        def get_zheng07_params(threshold):
            #Load tabulated data from Zheng et al. 2007, Table 1
            logM0_array = [11.2,10.59,11.49,11.69,11.38,11.84,11.92,13.94,14.0]
            logM1_array = [12.4,12.68,12.83,13.01,13.31,13.58,13.94,13.91,14.69]
            alpha_array = [0.83,0.97,1.02,1.06,1.06,1.12,1.15,1.04,0.87]
            # define the luminosity thresholds corresponding to the above data
            threshold_array = np.arange(-22,-17.5,0.5)
            threshold_array = threshold_array[::-1]

            threshold_index = np.where(threshold_array==threshold)[0]
            if len(threshold_index)==1:
                param_dict = {
                'logM0' : logM0_array[threshold_index[0]],
                'logM1' : logM1_array[threshold_index[0]],
                'alpha' : alpha_array[threshold_index[0]]
                }
            else:
                raise ValueError("Input luminosity threshold "
                    "does not match any of the Table 1 values of Zheng et al. 2007 (arXiv:0703457).")
            return param_dict

        if publication in ['zheng07', 'Zheng07', 'Zheng_etal07', 'zheng_etal07','zheng2007','Zheng2007']:
            param_dict = get_zheng07_params(threshold)
            return param_dict
        else:
            raise KeyError("For Zheng07Sats, only supported best-fit models are currently Zheng et al. 2007")


class Leauthaud11Sats(OccupationComponent):
    """ HOD-style model for any satellite galaxy occupation that derives from 
    a stellar-to-halo-mass relation. 
    """
    def __init__(self, smhm_model=smhm_components.Moster13SmHm, 
        gal_type = 'satellites', 
        threshold = model_defaults.default_stellar_mass_threshold, 
        prim_haloprop_key=model_defaults.prim_haloprop_key,
        modulate_with_cenocc = False, 
        **kwargs):
        """
        Parameters 
        ----------
        gal_type : string, optional keyword argument
            Name of the galaxy population being modeled. Default is ``satellites``.  

        threshold : float, optional keyword argument
            Stellar mass threshold of the mock galaxy sample. 
            Default value is specified in the `~halotools.empirical_models.model_defaults` module.

        prim_haloprop_key : string, optional keyword argument 
            String giving the column name of the primary halo property governing 
            the occupation statistics of gal_type galaxies. 
            Default value is specified in the `~halotools.empirical_models.model_defaults` module.

        smhm_model : object, optional keyword argument 
            Sub-class of `~halotools.empirical_models.smhm_components.PrimGalpropModel` governing 
            the stellar-to-halo-mass relation 

        redshift : float, optional keyword argument 
            Redshift of the stellar-to-halo-mass relation. Default is 0. 

        modulate_with_cenocc : bool, optional keyword argument 
            If True, the first satellite moment will be multiplied by the 
            the first central moment. Default is False. 

        gal_type_centrals : string, keyword argument
            Name of the central galaxy population whose first moment 
            modulates the `mean_occupation` method of `Zheng07Sats`. 
            Only required if ``modulate_with_cenocc`` is True. 

        Examples 
        --------
        >>> sat_model = Leauthaud11Sats()
        """

        if modulate_with_cenocc is True:
            if 'gal_type_centrals' not in kwargs.keys():
                raise KeyError("If ``modulate_with_cenocc`` is True, must also pass "
                    "the gal_type_centrals keyword.")
            gal_type_centrals = kwargs['gal_type_centrals']
            self.central_occupation_model = Leauthaud11Cens(
                gal_type=gal_type_centrals, threshold=threshold,
                prim_haloprop_key = prim_haloprop_key, smhm_model = smhm_model, 
                **kwargs)
        else:
            self.central_occupation_model = Leauthaud11Cens(
                threshold=threshold, prim_haloprop_key = prim_haloprop_key, 
                smhm_model = smhm_model, **kwargs)

        self.ancillary_model_dependencies = ['central_occupation_model']

        super(Leauthaud11Sats, self).__init__(
            gal_type=gal_type, threshold=threshold, 
            occupation_bound=float("inf"), 
            prim_haloprop_key = prim_haloprop_key, 
            **kwargs)

        self._initialize_param_dict()

        self.modulate_with_cenocc = modulate_with_cenocc


        self.publications = self.central_occupation_model.publications

    def mean_occupation(self, **kwargs):
        """ Expected number of central galaxies in a halo of mass halo_mass.
        See Equation 12-14 of arXiv:1103.2077.

        Parameters
        ----------        
        prim_haloprop : array, optional keyword argument
            array of masses of halo_table in the catalog

        halo_table : object, optional keyword argument 
            Data table storing halo catalog. 

        Returns
        -------
        mean_nsat : array
            Mean number of central galaxies in the halo of the input mass. 

        Examples 
        --------
        >>> sat_model = Leauthaud11Sats()
        >>> mean_nsat = sat_model.mean_occupation(prim_haloprop = 1.e13)

        Notes 
        -----
        Assumes constant scatter in the stellar-to-halo-mass relation. 
        """
        # Retrieve the array storing the mass-like variable
        if 'halo_table' in kwargs.keys():
            mass = kwargs['halo_table'][self.prim_haloprop_key]
        elif 'prim_haloprop' in kwargs.keys():
            mass = kwargs['prim_haloprop']
        else:
            function_name = "Leauthaud11Sats.mean_occupation"
            raise HalotoolsModelInputError(function_name)

        self._update_satellite_params()

        mean_nsat = (
            np.exp(-self._mcut/mass)*
            (mass/self._msat)**self.param_dict['alphasat']
            )

        if self.modulate_with_cenocc is True:
            mean_nsat *= self.central_occupation_model.mean_occupation(**kwargs)

        return mean_nsat

    def _initialize_param_dict(self):
        """ Set the initial values of ``self.param_dict`` according to 
        the SIG_MOD1 values of Table 5 of arXiv:1104.0928 for the 
        lowest redshift bin. 

        Notes 
        -----
        These values are only for ballpark purposes, and are 
        not self-consistent with arXiv:1104.0928, 
        because a different stellar-to-halo-mass relation is used here. 
        """

        self._msat_mcut_abcissa = np.logspace(9, 15, num=500)

        self.param_dict['alphasat'] = 1.0
        self.param_dict['bsat'] = 10.62
        self.param_dict['bcut'] = 1.47
        self.param_dict['betacut'] = -0.13
        self.param_dict['betasat'] = 0.859

        self._update_satellite_params()


    def _update_satellite_params(self):
        """ Private method to update the model parameters. 

        """

        # Tabulate the inverse stellar-to-halo-mass relation
        ordinates = self.central_occupation_model.smhm_model.mean_stellar_mass(
            prim_haloprop =self._msat_mcut_abcissa)
        spline_function = spline(ordinates, self._msat_mcut_abcissa)

        # Call the interpolater to compute the knee
        knee = spline_function(10.**self.threshold)

        self._msat = (
            1.e12*self.param_dict['bsat']*
            (knee / 1.e12)**self.param_dict['betasat'])

        self._mcut = (
            1.e12*self.param_dict['bcut']*
            (knee / 1.e12)**self.param_dict['betacut'])

























