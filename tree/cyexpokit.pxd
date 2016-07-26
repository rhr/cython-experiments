cimport numpy as np

cdef void lndexpm3(np.double_t[:,:,:] q, np.double_t[:] t, Py_ssize_t[:] qi,
                   np.uint8_t[:] tmask,
                   np.double_t[:,:,:] p, int ideg, np.double_t[:] wsp) nogil

cdef np.double_t logsumexp(np.double_t[:] a) nogil
