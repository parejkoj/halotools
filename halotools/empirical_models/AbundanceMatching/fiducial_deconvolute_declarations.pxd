""" Cython declarations for fiducial_deconvolute. 
"""

cdef extern from "include/fiducial_deconvolute.h":
    void convolved_fit(struct af_point * af_points, int num_af_points, 
        double * smm, double * mf, int MASS_BINS, double scatter, 
        int repeat, double sm_step);
