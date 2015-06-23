#!/usr/bin/env python

import numpy as np 
from astropy.table import Table 
from scipy.stats import spearmanr

from ..abunmatch import ConditionalAbunMatch
from .. import model_defaults
from ...sim_manager import FakeMock

from ..preloaded_subhalo_model_blueprints import Campbell15_blueprint


def test_Campbell15():
	"""
	prim_haloprop_key = 'mpeak'
	prim_galprop_key = 'stellar_mass'
	sec_haloprop_key = 'halo_zhalf'
	sec_galprop_key = 'ssfr'
	fake_data = FakeMock()
	sm_min = fake_data.galaxy_table['stellar_mass'].min()
	sm_max = fake_data.galaxy_table['stellar_mass'].max()
	sm_bins = np.logspace(np.log10(sm_min)-0.01, np.log10(sm_max)+0.01, 50)

	blueprint = Campbell15_blueprint(
		prim_haloprop_key=prim_haloprop_key, 
		prim_galprop_key=prim_galprop_key, 
		sec_haloprop_key=sec_haloprop_key, 
		galprop_key=sec_galprop_key)


	cam_noscatter = ConditionalAbunMatch(
		galprop_key=galprop_key, 
		prim_galprop_key = prim_galprop_key, 
		sec_haloprop_key = sec_haloprop_key, 
		input_galaxy_table = fake_data.galaxy_table, 
		prim_galprop_bins = sm_bins
		)


	"""

