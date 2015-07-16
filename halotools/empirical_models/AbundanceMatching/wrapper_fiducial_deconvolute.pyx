"""
Cython wrapper for Peter Behroozi's convolved_fit C function, 
the primary workhorse of SHAM deconvolution. 
""" 

import numpy as np

def fiducial_deconvolute(af, smm, mf, scatter, repeat=40, sm_step=0.01):
    """
    Parameters 
    ----------
    af : dict 
        Some description. 

    smm : array_like 
        Numpy array storing the stellar mass function 

    mf : array_like 
        Numpy array storing the subhalo mass function 

    scatter : float 
        Amount of scatter in the abundance matching relation (what units?)

    repeat : int, optional 
        Number of iterations to use in the deconvolution. Default is 40. 

    sm_step : float, optional 
        Some description. Default is 0.01. 

    Returns 
    -------
    smm : array_like 
        Why does this function return a modified input array 
        rather than a new array?

    Notes 
    -----
    Here is the place to put a free-form description of the gotchas. 

    """
    if len(smm) != len(mf):
        raise ValueError('`smf` and `mf` must have the same size!')
    sm_step = np.fabs(float(sm_step))
    sm_min = min(af['key'].min(), smm.min())
    if sm_min <= 0:
        offset = sm_step-sm_min
        af['key'] += offset
        smm += offset
    fiducial_deconvolute_declarations.convolved_fit(af, len(af), smm, mf, 
        len(mf), float(scatter), int(repeat), sm_step)
    if sm_min <= 0:
        smm -= offset
        af['key'] -= offset
    return smm

