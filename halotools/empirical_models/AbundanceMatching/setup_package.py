# Licensed under a 3-clause BSD style license
from __future__ import absolute_import

import os
from distutils.extension import Extension

ROOT = os.path.relpath(os.path.dirname(__file__))

def get_extensions():
    sources = (
    	[os.path.join(ROOT, 'wrapper_fiducial_deconvolute.pyx'), 
    	os.path.join(ROOT, 'src', 'fiducial_deconvolute.c')]
    	)
    ext = Extension(
        name="halotools.empirical_models.AbundanceMatching.wrapped_fiducial_deconvolute",
        sources=sources)
    return [ext]

def requires_2to3():
    return False

def get_package_data():
    return {'halotools.empirical_models.AbundanceMatching': ['include/*.h']}