from libc.stdlib cimport malloc, free

ctypedef struct Node

ctypedef struct Node:
    Node * leftchild
    Node * rightsib
    Node * parent
    Py_ssize_t i # postorder index

cdef class Tree(object):
    cdef:
        Node * _root

cdef inline Node * node_alloc(int n) nogil:
    cdef Node * node = <Node *>malloc(n * sizeof(Node))
    node.leftchild = node.rightsib = node.parent = NULL
    node.i = 0
    return node

cdef void add_child(Node * parent, Node * child):
    cdef Node * n = parent.leftchild
    if n == NULL:
        parent.leftchild = child
    else:
        while n.rightsib != NULL:
            n = n.rightsib
        n.rightsib = child
    child.parent = parent

cdef Node * n = node_alloc(1)
n.i = 0
cdef Node * m = node_alloc(1)
m.i = 1
add_child(n, m)
print n.leftchild.i

## n.label = 'foo'
## m.label = 'bar'
## print n.label
## n.leftchild = m
## print n.leftchild.label
## n.leftchild.rightsib = o
#n.leftchild.label = 'bar'
## n.leftchild.rightsib = <Node *>malloc(sizeof(Node))
## n.leftchild.rightsib.label = 'bar'
#add_child(n, <Node *>malloc(sizeof(Node)))
#n.leftchild.label = 'bar'
#print n.leftchild.parent.label
## free(n.leftchild.rightsib)
## free(n.leftchild)
free(n)

