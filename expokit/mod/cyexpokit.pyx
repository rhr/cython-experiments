import numpy as np
cimport numpy as np
from scipy.sparse import coo_matrix
from scipy.linalg import expm
from libc.math cimport pow

cdef extern:
    void wrapsingledmexpv(int n, double t, double* v, double * w, double* wsp, int lwsp, int* iwsp, int liwsp, int *ia,  int *ja,  double *a, int nz)

def test():
    cdef np.double_t[:,:] dense
    # rate matrix with rows as ancestors and columns as descendants
    dense = np.array([[-1,1,0,0],
                      [0,-2,2,0],
                      [0,0,-3,3],
                      [0,0,4,-4]], dtype=np.double)

    # !! transpose to get forward probs from dmexpv, given an initial
    # !! prob. vector
    ## Q = coo_matrix(dense.T)

    Q = coo_matrix(dense)

    # (otherwise, the returned vector is the column of conditional
    # probs given a vector of outcome probs; does not sum to 1 in this
    # case)

    cdef int n = Q.shape[0]
    cdef int m = min(n-1, 30)
    cdef int nz = Q.nnz
    ## cdef int[:] ia = Q.indptr.copy()
    ## cdef int[:] ja = Q.indices.copy()
    cdef int[:] ia = Q.row+1
    cdef int[:] ja = Q.col+1
    cdef double[:] a = Q.data
    ## cdef double[:] v = np.zeros(n)
    ## v[0] = 1
    cdef double[:] v = np.array([0.5, 0, 0, 0.5])
    cdef double[:] w = np.empty(n)
    cdef int ideg = 6
    ## cdef int lwsp = n*(m+1)+n+(m+2)**2+4*((m+2)**2)+ideg+1
    cdef int lwsp = n*(m+2)+n+5*(m+2)**2+ideg+1
    cdef double[:] wsp = np.empty(lwsp)
    cdef int liwsp = m+2
    cdef int[:] iwsp = np.empty(liwsp*2, dtype=np.int32) # !! double the length, otherwise crash
    cdef double t = 0.5
    print expm(np.asarray(dense)*t)
    wrapsingledmexpv(n, t, &v[0], &w[0], &wsp[0], lwsp, &iwsp[0], liwsp, &ia[0], &ja[0], &a[0], nz)
    print list(w)
    ## cdef int i
    ## for i in range(4):
    ##     v = np.zeros(n)
    ##     v[i] = 1
    ##     wrapsingledmexpv(&n, &m, &t1, &v[0], &w[0], &tol, &anorm, &wsp[0], &lwsp, &iwsp[0], &liwsp, &itrace, &iflag, &ia[0], &ja[0], &a[0], &nz)
    ##     print list(w)

