#sample inds from each species randomly, output their relationships as a newick

cd ~/Desktop
cd ~/Applications
source activate msprime-env
python3
import msprime, pyslim
import numpy as np
import matplotlib.pyplot as plt

ts = pyslim.load("3tree_lazy_t5.trees")
recap_ts = ts.recapitate(recombination_rate=1e-9, Ne=6000)
recap_ts.dump("3tree_lazy_t5_recap.trees")
ts = pyslim.load("3tree_lazy_t5_recap.trees")
sum([t.num_roots == 1 for t in ts.trees()])
ts.num_trees

#Draw one individual from each population

alive = ts.individuals_alive_at(0)
locs = ts.individual_locations[alive, :]

groups = {
   'A' : alive[np.logical_and(locs[:, 0] > 72, locs[:,0] < 92)],
   'B' : alive[np.logical_and(locs[:, 0] > 48, locs[:, 0] < 68)],
   'C' : alive[np.logical_and(locs[:, 0] > 24, locs[:, 0] < 44)],
   'D' : alive[np.logical_and(locs[:, 0] < 20, locs[:,0] > 0)]
    }

for k in groups:
	print(f"We have {len(groups[k])} individuals in the {k} group")

np.random.choice(groups['A'],1)
np.random.choice(groups['B'],1)
np.random.choice(groups['C'],1)
np.random.choice(groups['D'],1)


A = [8000, 10000, 12000, 14000]
sts = ts.simplify(A)
with open('nonspatial.trees', 'w') as f:
   for t in sts.trees():
      print(t.newick(), file=f)

#These numbers are the draws from the random choice functions above; they are the
#individual id's that you'll use from the edge table to search for genome id's

1833
1446
432
2867

#Long table, you'll search for the four id's above
print(ts.tables.nodes)

#The numbers contained in A are the genome id's not the individual id's; each individual
#has two genomes, so you can pick either one at random to put here
A = [22066, 21292, 19264, 24134] 
sts = ts.simplify(A)
with open('random_t45.trees', 'w') as f:
   for t in sts.trees():
      print(t.newick(), file=f)

sts.dump("territorial_t4_sp.trees")

#The newick function will print out all of the trees
#Copy and paste the entire console into a text file; select just the newick trees;
#should be ~5-10,000 or so, and delete everything else;
#do a search and replace: replace 1: with A:, 2: with B:, 3: with C:, and 4: with D:
#save as a .trees file and import to R