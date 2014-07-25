from cython.operator cimport dereference as deref, preincrement as inc
from cython import address as addr

cdef class Node

cdef class Node(object):
    "node in a rooted, potentially multifurcating tree"
    ## cdef object label
    cdef public Node parent, leftchild, leftsib, rightsib
    cdef readonly int nchildren
    cdef public float length
    cdef public int ni, pi, li, ii, left, right
    cdef public str label, treename
    ## cdef int _pre, _post, _ipre, _ipost, _epre, _epost

    property isleaf:
        def __get__(self):
            return self.nchildren == 0

    property isroot:
        def __get__(self):
            return self.parent is None

    property children:
        def __get__(self):
            return list(self.iter_children())

    property leaves:
        def __get__(self):
            pass
    
    def __cinit__(self):
        pass

    def __iter__(self):
        for c in self.iter_children():
            for n in c.iter_nodes():
                yield n

    cpdef add_child(self, Node n):
        cdef Node c
        if self.nchildren == 0:
            self.leftchild = n
        else:
            c = self.leftchild
            while c.rightsib is not None: c = c.rightsib
            c.rightsib = n
            n.leftsib = c
        n.parent = self
        self.nchildren += 1

    cpdef prune(self):
        cdef Node p = self.parent, lsib = self.leftsib, rsib = self.rightsib
        if p is not None:
            if lsib is None and rsib is None:
                # self is the only child of parent
                return p
            if lsib is None and rsib is not None:
                # self is the first child of parent
                p.leftchild = rsib
                rsib.leftsib = None
            elif lsib is not None and rsib is None:
                # self is the last child of parent
                lsib.rightsib = None
            elif lsib is not None and rsib is not None:
                # self has both left and right sibs
                lsib.rightsib = rsib
                rsib.leftsib = lsib
            else:
                pass
            p.nchildren -= 1
        self.parent = None
        self.leftsib = None
        self.rightsib = None
        return p
        
    def iter_children(self):
        cdef Node n = self.leftchild
        while n is not None:
            yield n
            n = n.rightsib
    
    def iter_nodes(self):
        """
        iterate over nodes descendant from self - including self
        """
        yield self
        for child in self.iter_children():
            for n in child.iter_nodes():
                yield n

    def bft(self):
        "breadth-first traversal of descendants"
        v = self.children
        while v:
            w = []
            for n in v:
                yield n
                w.extend(n.iter_children())
            v = w

        
