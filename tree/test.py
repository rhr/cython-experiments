import newick, tree, cyexpokit
import numpy as np

s = '((((Homo:0.21,Pongo:0.21)A:0.28,Macaca:0.49)B:0.13,Ateles:0.62)C:0.38,Galago:1.00)root;'
r = newick.parse(s)
t = tree.Tree(r)
data = dict(zip(t.leaf_labels(), [1,0,0,0,0]))

qidx = np.array([[0,0,1,0],
                 [0,1,0,0]])

f = cyexpokit.make_mklnl_func(t, data, 2, 1, qidx)
