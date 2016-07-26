# cython: profile=True
import numpy as np
from libc.math cimport exp, log
from numpy.math cimport INFINITY
cimport cython
from tree cimport Tree

# Indexing arrays - what type to use?
# -----------------------------------
#
# TL;DR: int is probably fine; BUT
# 
# In Cython, for maximum correctness and portability, variables used
# for indexing arrays should be declared Py_ssize_t, not int. This
# raises a question: If you want to create/use a numpy array of
# indices, what type should it be? Py_ssize_t has no defined
# equivalent in numpy. The closest seems to be np.intp: see
# 
# https://github.com/scikit-learn/scikit-learn/wiki/C-integer-types:-the-missing-manual
#
# and
#
# https://github.com/numpy/numpy/issues/1654
#
# So, e.g:
#
# cdef np.ndarray[np.intp_t, ndim=1] my_index_array = np.arange(10, dtype=np.intp)
#
# (Note np.intp_t on the cdef side, np.intp on the Python side.)

cdef extern:
    void f_dexpm(int nstates, double* H, double t, double* expH) nogil
    void f_dexpm_wsp(int nstates, double* H, double t, int i,
                     double* wsp, double* expH) nogil

@cython.boundscheck(False)
cdef void dexpm3(np.double_t[:,:,:] q, np.double_t[:] t, Py_ssize_t[:] qi,
                 np.uint8_t[:] tmask,
                 np.double_t[:,:,:] p, int ideg, np.double_t[:] wsp) nogil:
    """
    Compute transition probabilities exp(q*t) for a 'stack' of q
    matrices, over all values of t from a 1-d array, where q is selected
    from the stack by indices in qi. Uses pre-allocated arrays for
    intermediate calculations and output, to minimize overhead of
    repeated calls (e.g. for ML optimization or MCMC).

    Args:

        q (np.double_t[m,k,k]): stack of m square rate matrices of dimension k

        t (np.double_t[n]): 1-d array of times (branch lengths) of length n

        qi (int[n]): 1-d array indicating assigning q matrices to times

        tmask (bint[n]): 1-d array indicating which t values to process
        
        p (np.double_t[n,k,k]): stack of n square p matrices holding results
          of exponentiation, i.e., p[i] = exp(q[qi[i]]*t[i])

        ideg (int): used in expokit Fortran code; a good default is 6

        wsp (np.double_t[:]): expokit "workspace" array, must have
          min. length = 4*k*k+ideg+1
    """
    cdef Py_ssize_t i
    cdef int nstates = q.shape[1]
    for i in range(t.shape[0]):
        if tmask[i]:
            f_dexpm_wsp(nstates, &q[qi[i],0,0], t[i], ideg, &wsp[0], &p[i,0,0])

@cython.boundscheck(False)
cdef void lndexpm3(np.double_t[:,:,:] q, np.double_t[:] t, Py_ssize_t[:] qi,
                   np.uint8_t[:] tmask,
                   np.double_t[:,:,:] p, int ideg, np.double_t[:] wsp) nogil:
    "same as dmexp3, but log-transforms p"
    cdef Py_ssize_t i, j, k
    cdef int nstates = q.shape[1]
    for i in range(t.shape[0]):
        if tmask[i]:
            f_dexpm_wsp(nstates, &q[qi[i],0,0], t[i], ideg, &wsp[0], &p[i,0,0])
            for j in range(nstates):
                for k in range(nstates):
                    p[i,j,k] = log(p[i,j,k])

def test_dexpm3(int N=1000):
    m = 3  # number of q matrices
    k = 4  # number of states
    q = np.zeros((m,k,k))  # 3-d 'stack' of q matrices
    for i in range(m):
        # fill the off-diagonals of each q matrix, incrementing:
        # 0.1, 0.2, ...
        a = q[i]
        a += 0.1*(i+1)
        # make the rows sum to zero
        np.fill_diagonal(a, 0)
        a[np.diag_indices_from(a)] = -a.sum(axis=1)
    t = np.ones(3)
    tmask = np.ones(3, dtype=np.uint8)
    qi = np.array([0,1,2], dtype=np.intp)
    
    cdef int x, ideg = 6
    wsp = np.empty(4*k*k+ideg+1)
    n = len(t)
    p = np.empty((n,k,k))
    
    # with all arrays allocated, can call dexpm3
    for x in range(N):
        dexpm3(q, t, qi, tmask, p, ideg, wsp)
    return p

def test_lndexpm3(int N=1000):
    m = 3  # number of q matrices
    k = 4  # number of states
    q = np.zeros((m,k,k))  # 3-d 'stack' of q matrices
    for i in range(m):
        # fill the off-diagonals of each q matrix, incrementing:
        # 0.1, 0.2, ...
        a = q[i]
        a += 0.1*(i+1)
        # make the rows sum to zero
        np.fill_diagonal(a, 0)
        a[np.diag_indices_from(a)] = -a.sum(axis=1)
    t = np.ones(3)
    tmask = np.ones(3, dtype=np.uint8)
    qi = np.array([0,1,2], dtype=np.intp)
    
    cdef int x, ideg = 6
    wsp = np.empty(4*k*k+ideg+1)
    n = len(t)
    p = np.empty((n,k,k))
    
    # with all arrays allocated, can call dexpm3
    for x in range(N):
        lndexpm3(q, t, qi, tmask, p, ideg, wsp)
    return p

@cython.boundscheck(False)
cdef np.double_t mklnl(Tree t,
                       np.double_t[:,:] fraclnl,
                       np.double_t[:,:,:] p,
                       int k,
                       np.double_t[:] tmp) nogil:
    """
    Standard Mk log-likelihood calculator.

    Args:

    t (Tree)

    fraclnl (np.double_t[m,k]): array to hold computed fractional log-likelihoods,
      where m = number of nodes, including leaf nodes; k = number of states
      * fraclnl[i,j] = fractional log-likelihood of node i for charstate j
      * leaf node values should be pre-filled, e.g. 0 for observed state,
        -np.inf everywhere else
      * this function calculates the internal node values (where a node
        could be a branch 'knee') and returns the log-likelihood at the root

    p (np.double_t[m,k,k]): p matrix

    k (int): number of states

    tmp (np.double_t[k]): to hold intermediate values
    """

    cdef Py_ssize_t i, parent, j, child, ancstate, childstate

    # For each internal node (in postorder sequence)...
    for i in range(t.nnodes):
        parent = t.postorder[i]
        if t.nchildren[parent]==0:
            continue
        # parent indexes the current internal node
        # For each child of this node...
        child = t.leftchild[parent]
        for j in range(t.nchildren[parent]):
            for ancstate in range(k):
                for childstate in range(k):
                    # Multiply child's likelihood by p-matrix entry
                    # (addition in log space)
                    tmp[childstate] = (p[child,ancstate,childstate] +
                                       fraclnl[child,childstate])
                # Sum of log-likelihoods of children
                if fraclnl[parent,ancstate] == -INFINITY:
                    fraclnl[parent,ancstate] = logsumexp(tmp)
                else:
                    fraclnl[parent,ancstate] += logsumexp(tmp)
            child = t.rightsib[child]

    # logsumexp
    cdef np.double_t result = 0.0, largest = fraclnl[0,0]
    for i in range(1, k):
        if (fraclnl[0,i] > largest):
            largest = fraclnl[0,i]
    for i in range(k):
        result += exp(fraclnl[0,i] - largest)
    return largest + log(result)

@cython.boundscheck(False)
cdef np.double_t logsumexp(np.double_t[:] a) nogil:
    """
    nbviewer.jupyter.org/gist/sebastien-bratieres/285184b4a808dfea7070
    Faster than scipy.misc.logsumexp
    """
    cdef Py_ssize_t i, n = a.shape[0]
    cdef np.double_t result = 0.0, largest_in_a = a[0]
    for i in range(1, n):
        if (a[i] > largest_in_a):
            largest_in_a = a[i]
    for i in range(n):
        result += exp(a[i] - largest_in_a)
    return largest_in_a + log(result)

def make_mklnl_func(Tree tree, data, int k, int nq, Py_ssize_t[:,:] qidx):
    ## cdef np.double_t[:,:,:] p = np.empty((tree.nnodes, k, k), dtype=np.double)
    cdef np.double_t[:] t = tree.length.copy()
    cdef np.double_t[:,:,:] p = np.empty((tree.nnodes, k, k), dtype=np.double)
    cdef int i, nc, ideg = 6
    cdef np.double_t[:] wsp = np.empty(4*k*k+ideg+1)
    cdef Py_ssize_t[:] qi = np.zeros(tree.nnodes, dtype=np.intp)
    cdef np.uint8_t[:] tmask = np.ones(tree.nnodes, dtype=np.uint8)
    tmask[0] = 0
    cdef np.double_t[:,:] fraclnl = np.empty((tree.nnodes, k), dtype=np.double)
    fraclnl[:] = -INFINITY
    for i, nc in enumerate(tree.nchildren):
        if nc == 0:
            fraclnl[i, <Py_ssize_t>data[tree.label[i]]] = 0
    cdef np.double_t[:] tmp = np.empty(k, dtype=np.double)
    cdef np.ndarray q = np.zeros((nq,k,k), dtype=np.double)
    
    def f(np.double_t[:] params):
        """
        params: array of free rate parameters, assigned to q by indices in qidx

        qidx columns:
            0, 1, 2 - index axes of q
            3 - index of params

        This scheme allows flexible specification of models. E.g.:

        Symmetric mk2:
            params = [0.2]; qidx = [[0,0,1,0],[0,1,0,0]]
            
        Asymmetric mk2:
            params = [0.2,0.6]; qidx = [[0,0,1,0],[0,1,0,1]]
        """
        cdef Py_ssize_t i, j, r, a, b, c, d
        cdef np.double_t x = 0
        for r in range(qidx.shape[0]):
            a = qidx[r,0]; b = qidx[r,1]; c = qidx[r,2]; d = qidx[r,3]
            q[a,b,c] = params[d]
        for r in range(nq):
            for i in range(k):
                x = 0
                for j in range(k):
                    if i != j:
                        x -= q[r,i,j]
                q[r,i,i] = x
        lndexpm3(q, t, qi, tmask, p, ideg, wsp)
        ## np.log(p, out=p)
        return mklnl(tree, fraclnl, p, k, tmp)

    # attached allocated arrays to function object
    f.t = t
    f.fraclnl = fraclnl
    f.q = q
    f.p = p
    f.qi = qi
    f.tmask = tmask
    return f
