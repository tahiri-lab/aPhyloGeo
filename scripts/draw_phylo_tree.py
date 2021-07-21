import pylab
from Bio import Phylo
import sys
import os
import matplotlib

matplotlib.use('Agg')

file_ = sys.argv[1]

tree = Phylo.read(file_, 'newick')
Phylo.draw(tree, do_show=True)

workdir = os.getcwd()

pylab.savefig('Phylo.jpg', dpi=2000)

print("Phylogenetic tree generated and saved as Phylo.jpg")
print("Phylogenetic tree generated and saved in directory:%s" % workdir)

