__all__ = ['fiducial_deconvolute']
import os
import numpy as np
import numpy.ctypeslib as C

cimport fiducial_deconvolute_declarations

# define types
_double_ctype = C.ndpointer(np.float64, ndim=1, flags='C')

# load the c library
here = os.path.abspath(os.path.dirname(__file__))
_C_LIB = C.load_library('fiducial_deconvolute', here)
_C_LIB.convolved_fit.restype = None
_C_LIB.convolved_fit.argtypes = [_double_ctype, _double_ctype, C.ctypes.c_int, \
        _double_ctype, _double_ctype, C.ctypes.c_int, C.ctypes.c_double,\
        C.ctypes.c_int, C.ctypes.c_double]

def fiducial_deconvolute(af_key, af_val, smm, mf, scatter, repeat=40, sm_step=0.01):
    if len(smm) != len(mf):
        raise ValueError('`smf` and `mf` must have the same size!')
    sm_step = np.fabs(float(sm_step))
    sm_min = min(af_key.min(), smm.min())
    if sm_min <= 0:
        offset = sm_step-sm_min
        af_key += offset
        smm += offset
    # _C_LIB.convolved_fit(af_key, af_val, len(af_key), 
    #     smm, mf, len(mf), float(scatter), int(repeat), sm_step)
    fiducial_deconvolute_declarations.convolved_fit(af_key, af_val, len(af_key), 
        smm, mf, len(mf), float(scatter), int(repeat), sm_step)
    if sm_min <= 0:
        smm -= offset
        af_key -= offset
    return smm

