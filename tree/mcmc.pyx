cimport tree, cyexpokit
from numpy.math cimport NAN

def make_switchpoint_tree(root):
    nodes = root.iternodes()
    nodes.next() # skip root
    for n in list(nodes):
        x = tree.Node()
        p = n.prune()
        p.add_child(x)
        x.add_child(n)
        if n.length != NAN:
            n.length /= 2.0
            x.length = n.length
    for i, n in enumerate(root.iternodes()): n.ni = i
