__all__ = ['fiducial_deconvolute']

import numpy as np
cimport numpy as np 

cimport fiducial_deconvolute_declarations

def fiducial_deconvolute(af_key, af_val, smm, mf, scatter, repeat=40, sm_step=0.01):
    if len(smm) != len(mf):
        raise ValueError('`smf` and `mf` must have the same size!')
    sm_step = np.fabs(float(sm_step))
    sm_min = min(af_key.min(), smm.min())
    if sm_min <= 0:
        offset = sm_step-sm_min
        af_key += offset
        smm += offset
    fiducial_deconvolute_declarations.convolved_fit(af_key, af_val, len(af_key), 
        smm, mf, len(mf), float(scatter), int(repeat), sm_step)
    if sm_min <= 0:
        smm -= offset
        af_key -= offset
    return smm

