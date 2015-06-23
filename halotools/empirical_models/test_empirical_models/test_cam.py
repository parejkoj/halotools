#!/usr/bin/env python

import numpy as np 
from astropy.table import Table 
from scipy.stats import spearmanr

from ..abunmatch import ConditionalAbunMatch
from .. import model_defaults
from ...sim_manager import FakeMock

def retrieve_prim_galprop_subsample(model, mock, lower_bound, upper_bound):
	mask = np.where(
		(mock.galaxy_table[model.prim_galpropkey] > lower_bound) & 
		(mock.galaxy_table[model.prim_galpropkey] < upper_bound))[0]
	return mock.galaxy_table[mask]

def check_conditional_one_point(model, mock, data, 
	lower_bound, upper_bound):

	mock_subsample = retrieve_prim_galprop_subsample(mock, 
		lower_bound, upper_bound)

	data_subsample = retrieve_prim_galprop_subsample(data, 
		lower_bound, upper_bound)

	data_mean = data_subsample[mock.galprop_key].mean()
	mock_mean= mock_subsample[mock.galprop_key].mean()
	mean_fracdiff = (mock_mean-data_mean)/data_mean
	assert np.allclose(mean_fracdiff, 0, atol=0.3, rtol=0.3)

def check_range(model, mock, data):
	min_mock_galprop = mock.galaxy_table[model.galprop_key].min()
	max_mock_galprop = mock.galaxy_table[model.galprop_key].max()
	min_data_galprop = data.galaxy_table[model.galprop_key].min()
	max_data_galprop = data.galaxy_table[model.galprop_key].max()

	assert min_data_galprop > -10
	assert max_data_galprop < 10

	assert min_mock_galprop > -10

#	min_galprop_fracdiff = (min_mock_galprop-min_data_galprop)/min_data_galprop
#	max_galprop_fracdiff = (max_mock_galprop-max_data_galprop)/max_data_galprop
#	assert np.allclose(min_galprop_fracdiff, 0, atol=0.3, rtol=0.3)
#	assert np.allclose(max_galprop_fracdiff, 0, atol=0.3, rtol=0.3)

def check_spearmanr(model, mock, lower_bound, upper_bound, desired_correlation):

	mock_subsample = retrieve_prim_galprop_subsample(mock, 
		lower_bound, upper_bound)

	corr = spearmanr(mock_subsample[mock.galprop_key], 
		mock_subsample[mock.sec_haloprop_key])[0]
	corr_fracdiff = (corr-desired_correlation)/desired_correlation
	assert np.allclose(corr_fracdiff, 0, atol=0.3, rtol=0.3)

#############################################

def test_cam_no_scatter():

	galprop_key = 'gr_color'
	prim_galprop_key = 'stellar_mass'
	sec_haloprop_key = 'zhalf'

	fake_data = FakeMock()
	sm_min = fake_data.galaxy_table[prim_galprop_key].min()
	sm_max = fake_data.galaxy_table[prim_galprop_key].max()
	sm_bins = np.logspace(np.log10(sm_min)-0.01, np.log10(sm_max)+0.01, 50)

	cam_noscatter = ConditionalAbunMatch(
		galprop_key=galprop_key, 
		prim_galprop_key = prim_galprop_key, 
		sec_haloprop_key = sec_haloprop_key, 
		input_galaxy_table = fake_data.galaxy_table, 
		prim_galprop_bins = sm_bins
		)
	fake_mock_noscatter = FakeMock(approximate_ngals = 1e4)
	fake_mock_noscatter.galaxy_table[galprop_key] = (
		cam_noscatter.mc_gr_color(galaxy_table = fake_mock_noscatter.galaxy_table))

	assert hasattr(cam_noscatter, 'galprop_key')

	# Check no-scatter mock
	check_range(cam_noscatter, fake_mock_noscatter, fake_data)
	#sm_low, sm_high = 1.e10, 5.e10
	#check_conditional_one_point(fake_mock_noscatter, fake_data, 
#		'stellar_mass', 'gr_color', sm_low, sm_high)

#	check_spearmanr(fake_mock_noscatter, fake_data, sm_low, sm_high, 0.99)

#	check_spearmanr(fake_mock_noscatter, 'stellar_mass', 'gr_color', sec_haloprop_key, 
#		lower_bound, upper_bound, desired_correlation)


#	sm_low, sm_high = 5.e10, 1.e11
#	check_conditional_one_point(fake_mock_noscatter, fake_data, 
#		'stellar_mass', 'gr_color', sm_low, sm_high)
#	check_spearmanr(fake_mock_noscatter, fake_data, sm_low, sm_high, 0.99)


def test_cam_gr_color():
	galprop_key = 'gr_color'
	prim_galprop_key = 'stellar_mass'
	sec_haloprop_key = 'zhalf'

	fake_data = FakeMock()
	sm_min = fake_data.galaxy_table[prim_galprop_key].min()
	sm_max = fake_data.galaxy_table[prim_galprop_key].max()
	sm_bins = np.logspace(np.log10(sm_min)-0.01, np.log10(sm_max)+0.01, 50)

	cam_scatter_50 = ConditionalAbunMatch(
		galprop_key=galprop_key, 
		prim_galprop_key = prim_galprop_key, 
		sec_haloprop_key = sec_haloprop_key, 
		input_galaxy_table = fake_data.galaxy_table, 
		prim_galprop_bins = sm_bins, 
		correlation_strength = 0.5
		)
	fake_mock_scatter_50 = FakeMock(approximate_ngals = 1e4)
	fake_mock_scatter_50.galaxy_table[galprop_key] = (
		cam_scatter_50.mc_gr_color(galaxy_table = fake_mock_scatter_50.galaxy_table))

	cam_variable_scatter = ConditionalAbunMatch(
		galprop_key=galprop_key, 
		prim_galprop_key = prim_galprop_key, 
		sec_haloprop_key = sec_haloprop_key, 
		input_galaxy_table = fake_data.galaxy_table, 
		prim_galprop_bins = sm_bins, 
		correlation_strength = [0.25, 0.75], 
		correlation_strength_abcissa = [2.e10, 7.e10]
		)
	fake_mock_variable_scatter = FakeMock(approximate_ngals = 1e4)
	fake_mock_variable_scatter.galaxy_table[galprop_key] = (
		cam_variable_scatter.mc_gr_color(galaxy_table = fake_mock_variable_scatter.galaxy_table))


	# Check mock with 50% correlation strength
#	check_range(fake_mock_scatter_50, fake_data)
#	sm_low, sm_high = 1.e10, 5.e10
#	check_conditional_one_point(fake_mock_scatter_50, fake_data, 
#		'stellar_mass', 'gr_color', sm_low, sm_high)
#	check_spearmanr(fake_mock_scatter_50, fake_data, sm_low, sm_high, 0.5)

#	sm_low, sm_high = 5.e10, 1.e11
#	check_conditional_one_point(fake_mock_scatter_50, fake_data, 
#		'stellar_mass', 'gr_color', sm_low, sm_high)
#	check_spearmanr(fake_mock_scatter_50, fake_data, sm_low, sm_high, 0.5)



	# Check mock with variable correlation strength
#	check_range(fake_mock_variable_scatter, fake_data)
#	sm_low, sm_high = 1.e10, 5.e10
#	check_conditional_one_point(fake_mock_variable_scatter, fake_data, 
#		'stellar_mass', 'gr_color', sm_low, sm_high)
#	check_spearmanr(fake_mock_variable_scatter, fake_data, sm_low, sm_high, 0.34)

#	sm_low, sm_high = 5.e10, 1.e11
#	check_conditional_one_point(fake_mock_variable_scatter, fake_data, 
#		'stellar_mass', 'gr_color', sm_low, sm_high)
#	check_spearmanr(fake_mock_variable_scatter, fake_data, sm_low, sm_high, 0.835)

