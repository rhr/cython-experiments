import newick, tree, mcmc, random, theano
import pymc3 as pm
import theano.tensor as T
import numpy as np

def stem_subtreelen(t):
    a = np.empty(t.nnodes)
    for i in t.postorder:
        nc = t.nchildren[i]
        if nc == 0:
            a[i] = t.length[i]
        else:
            c = t.leftchild[i]
            x = 0.0
            for j in range(nc):
                x += a[c]
                c = t.rightsib[c]
            if i > 0: # not root
                x += t.length[i]
            a[i] = x
    return a

def make_switchpoint_tree(root):
    nodes = root.iternodes()
    nodes.next() # skip root
    for n in list(nodes):
        x = tree.Node()
        p = n.prune()
        p.add_child(x)
        x.add_child(n)
        if n.length != np.nan:
            n.length /= 2.0
            x.length = n.length
    for i, n in enumerate(root.iternodes()): n.ni = i
    return tree.Tree(root)

r = newick.parse(open('/home/rree/src/ivy/examples/primates.newick').read())
t = tree.Tree(r)
#t = make_switchpoint_tree(r)
tlen = np.sum(t.length)
nh = tree.nodeheights(t)
stl = stem_subtreelen(t)

def get_ith_child(t, n, i):
    c = t.leftchild[n]
    for i in range(i):
        c = t.rightsib[c]
    return c

## @theano.compile.ops.as_op(itypes=[T.dscalar], otypes=[T.iscalar])
## def sample_path(u):
##     v = 0.0
##     i = 0
##     while v < u:
##         i = get_ith_child(t, i, np.random.choice(t.nchildren[i]))
##         v += t.length[i]
##     return np.array([i])
        
## #@theano.compile.ops.as_op(itypes=[T.dscalar, T.iscalar], otypes=[T.dscalar])
## def logp(u, c):
##     if T.le(u, 0.5):
##         return T.log(T.switch(T.eq(c,0), 0.1, 0.8))
##     else:
##         return T.log(0.05)

with pm.Model() as model:
    u = pm.Uniform('u')
    dd = {}
    cc = {}
    for i in range(t.nnodes):
        nc = t.nchildren[i]
        if nc > 2:
            dd[i] = pm.Dirichlet('d', np.array([0.1]*nc), shape=nc)
            cc[i] = pm.Categorical('c', d)

    lnl = pm.DensityDist('lnl', logp, observed=dict(u=u, c=c))
    tr = pm.sample(1000, pm.Slice([d]))

    ## g = {}
    ## for i in range(t.nnodes):
    ##     nc = t.nchildren[i]
    ##     if nc > 1:
    ##       g[i] = pm.Dirichlet('dir_{}'.format(i), np.array([0.1]*nc), shape=nc)
    ## x = pm.DensityDist('x', tst)
    ## tr = pm.sample(10)
