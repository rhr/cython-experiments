cimport numpy as np

cdef class Node(object):
    cdef public Node parent, leftchild, leftsib, rightsib
    cdef readonly int nchildren
    cdef public double length
    cdef public int ni, left, right
    cdef public str label, treename
    cpdef Node prune(self)
    cpdef add_child(self, Node n)

cdef class Tree(object):
    cdef:
        readonly Node root
        readonly Py_ssize_t[:] parent
        readonly Py_ssize_t[:] leftchild
        readonly Py_ssize_t[:] rightsib
        readonly np.int32_t[:] nchildren
        readonly Py_ssize_t[:] postorder
        readonly Py_ssize_t[:] postith
        readonly np.double_t[:] length
        readonly int nnodes
        readonly list label

    cpdef index(self, root=*)
    cpdef double rootpathlen(self, Py_ssize_t i, Py_ssize_t j=*)
    cdef int subtree_nnodes(self, Py_ssize_t node)
    cdef int subtree_nleaves(self, Py_ssize_t node)
    ## cpdef traverse(self, Py_ssize_t i)
