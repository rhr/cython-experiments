import newick, mcmc
r = newick.parse(open('/home/rree/src/ivy/examples/primates.newick').read())
mcmc.insert_knees(r)
