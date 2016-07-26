"""
Numerical integration using the GSL: https://goo.gl/ZOF5K3
Requires CythonGSL: https://github.com/twiecki/CythonGSL
"""
from __future__ import print_function
from cython_gsl cimport *
import numpy as np
cimport numpy as np

## ctypedef int (* odefunc) (
##     double t, const double y[], double dydt[], void * params) nogil

# simple example: binary Markov model (mk2)
cdef int mk2(double t, double y[], double f[], void *params) nogil:
    cdef double * Q = <double *>params
    cdef double q01 = Q[0]
    cdef double q10 = Q[1]
    cdef double p0 = y[0], p1 = y[1]
    f[0] = -q01*p0 + q01*p1
    f[1] = -q10*p1 + q10*p0
    return GSL_SUCCESS

def integrate():
    cdef double[:] params = np.array([0.3,0.1])
    cdef int ndim = 2
    cdef gsl_odeiv2_system sys
    sys.function = mk2
    sys.dimension = ndim
    sys.params = <void *>&params[0]

    cdef double hstart = 1e-6 # initial step size
    # keep the local error on each step within:
    cdef double epsabs = 1e-8 # absolute error
    cdef double epsrel = 0.0  # relative error

    cdef gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new(
        &sys, gsl_odeiv2_step_rkf45, hstart, epsabs, epsrel)
    
    cdef int i
    cdef double t, t1#, y[ndim]
    cdef double[:] y = np.empty(ndim, dtype=np.double)
    t = 0.0
    t1 = 1.0
    y[0] = 0.0
    y[1] = 1.0

    cdef int status, nsteps = 100
    cdef double ti, x = t1/nsteps
    for i in range(nsteps):
        ti = (i+1) * x
        status = gsl_odeiv2_driver_apply(d, &t, ti, &y[0])

        if (status != GSL_SUCCESS):
            print("error, return value=%d" % status)
            break

        print("%.5e %.5e %.5e" %(t, y[0], y[1]))

    t = 0.0
    y[0] = 0.0
    y[1] = 1.0
    status = gsl_odeiv2_driver_apply(d, &t, ti, &y[0])

    if (status != GSL_SUCCESS):
        print("error, return value=%d" % status)
    print("%.5e %.5e %.5e" %(t, y[0], y[1]))

    gsl_odeiv2_driver_free(d)
    return y

