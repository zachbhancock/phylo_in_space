#Python code
source activate msprime-env
python3
import msprime, pyslim
import numpy as np
import matplotlib.pyplot as plt

ts = pyslim.load("RF_div_sigma1.trees")
recap_ts = ts.recapitate(recombination_rate=1e-9, Ne=16000)
ts = msprime.mutate(recap_ts, rate=1e-8, model=msprime.InfiniteSites(msprime.NUCLEOTIDES), random_seed=1)
ts.dump("RF_div_sigma1_recap.trees")
ts = pyslim.load("RF_div_sigma25_recap.trees")
sum([t.num_roots == 1 for t in ts.trees()])
ts.num_trees

#Draw one individual from each population

alive = ts.individuals_alive_at(0)
locs = ts.individual_locations[alive, :]

#random
groups_1 = {
   'A' : alive[np.logical_and(locs[:, 0] > 120, locs[:, 0] < 140)],
   'B' : alive[np.logical_and(locs[:, 0] > 96, locs[:, 0] < 116)],
   'C' : alive[np.logical_and(locs[:, 0] > 72, locs[:,0] < 92)],
   'D' : alive[np.logical_and(locs[:, 0] > 48, locs[:,0] < 68)],
   'E' : alive[np.logical_and(locs[:, 0] > 24, locs[:,0] < 44)],
   'F' : alive[np.logical_and(locs[:, 0] > 0, locs[:,0] < 20)],
    }

#close
groups_2 = {
   'A' : alive[np.logical_and(np.logical_and(np.logical_and(locs[:, 0] > 120, locs[:, 0] < 122), locs[:,1] < 17), locs[:,1] > 13)],
   'B' : alive[np.logical_and(np.logical_and(np.logical_and(locs[:, 0] > 114, locs[:, 0] < 116), locs[:,1] < 17), locs[:,1] > 13)],
   'C' : alive[np.logical_and(np.logical_and(np.logical_and(locs[:, 0] > 72, locs[:,0] < 75), locs[:,1] < 17), locs[:,1] > 13)],
   'D' : alive[np.logical_and(np.logical_and(np.logical_and(locs[:, 0] > 66, locs[:,0] < 68), locs[:,1] < 17), locs[:,1] > 13)],
   'E' : alive[np.logical_and(np.logical_and(np.logical_and(locs[:, 0] > 24, locs[:,0] < 26), locs[:,1] < 17), locs[:,1] > 13)],
   'F' : alive[np.logical_and(np.logical_and(np.logical_and(locs[:, 0] > 18, locs[:,0] < 20), locs[:,1] < 17), locs[:,1] > 13)]
    }
    
#far
groups_3 = {
   'A' : alive[np.logical_and(np.logical_and(np.logical_and(locs[:, 0] > 138, locs[:, 0] < 140), locs[:,1] < 17), locs[:,1] > 13)],
   'B' : alive[np.logical_and(np.logical_and(np.logical_and(locs[:, 0] > 96, locs[:, 0] < 98), locs[:,1] < 17), locs[:,1] > 13)],
   'C' : alive[np.logical_and(np.logical_and(np.logical_and(locs[:, 0] > 72, locs[:,0] < 74), locs[:,1] < 17), locs[:,1] > 13)],
   'D' : alive[np.logical_and(np.logical_and(np.logical_and(locs[:, 0] > 48, locs[:,0] < 50), locs[:,1] < 17), locs[:,1] > 13)],
   'E' : alive[np.logical_and(np.logical_and(np.logical_and(locs[:, 0] > 24, locs[:,0] < 26), locs[:,1] < 17), locs[:,1] > 13)],
   'F' : alive[np.logical_and(np.logical_and(np.logical_and(locs[:, 0] > 0, locs[:,0] < 2), locs[:,1] < 17), locs[:,1] > 13)]
    }

#center
groups_4 = {
   'A' : alive[np.logical_and(np.logical_and(np.logical_and(locs[:, 0] > 128, locs[:, 0] < 132), locs[:,1] < 17), locs[:,1] > 13)],
   'B' : alive[np.logical_and(np.logical_and(np.logical_and(locs[:, 0] > 104, locs[:, 0] < 108), locs[:,1] < 17), locs[:,1] > 13)],
   'C' : alive[np.logical_and(np.logical_and(np.logical_and(locs[:, 0] > 80, locs[:,0] < 84), locs[:,1] < 17), locs[:,1] > 13)],
   'D' : alive[np.logical_and(np.logical_and(np.logical_and(locs[:, 0] > 56, locs[:,0] < 60), locs[:,1] < 17), locs[:,1] > 13)],
   'E' : alive[np.logical_and(np.logical_and(np.logical_and(locs[:, 0] > 32, locs[:,0] < 36), locs[:,1] < 17), locs[:,1] > 13)],
   'F' : alive[np.logical_and(np.logical_and(np.logical_and(locs[:, 0] > 8, locs[:,0] < 12), locs[:,1] < 17), locs[:,1] > 13)]
    }

for k in groups_3:
	print(f"We have {len(groups_3[k])} individuals in the {k} group_3")

A = np.random.choice(groups_1['A'],5, replace=False)
B = np.random.choice(groups_1['B'],5, replace=False)
C = np.random.choice(groups_1['C'],5, replace=False)
D = np.random.choice(groups_1['D'],5, replace=False)
E = np.random.choice(groups_1['E'],5, replace=False)
F = np.random.choice(groups_1['F'],5, replace=False)

keep_indivs = np.concatenate((A, B, C, D, E, F), axis=None)

2314, 1896, 2219, 2475, 1701, 3125, 3061, 3127, 3379, 2956, 3574,
       3457, 3401, 3829, 3469, 5028, 4521, 4953, 4601, 4775,  851, 1281,
       1437,  792, 1279,  415,  436,   26,  648,  466

print(ts.tables.nodes)

a = [32628, 31792, 32438, 32950, 31402, 34250, 34122, 34254, 34758, 33912, 35148, 34914, 34802, 35658, 34938, 38056, 37042, 37906, 37202, 37550, 29702, 30562, 30874, 29584, 30558, 28830, 28872, 28052, 29296, 28932]
sts = ts.simplify(a)



group_order = ['A', 'B', 'C', 'D', 'E', 'F']
sampled_nodes = [[] for _ in groups_1]
for j, k in enumerate(group_order):
   for ind in groups_1[k]:
      sampled_nodes[j].extend(ts.individual(ind).nodes)

pairs = [(i, j) for i in range(6) for j in range(6)]
group_div = ts.divergence(sampled_nodes, indexes=pairs).reshape((6, 6))

print("\t" + "\t".join(group_order))
for i, group in enumerate(group_order):
   print(f"{group_order[i]}:\t" + "\t".join(map(str, np.round(group_div[i], 7))))

Fst = ts.Fst(sampled_nodes, indexes = pairs)

haps = []
for i in sts.haplotypes():
    haps.append(i)

sequence_IDs = []
for i in range(len(haps)):
    sequence_IDs.append(f'sample_{sts.samples()[i]}_pop_{sts.node(i).population}')

with open('RF_div_sigma25_random_spdel.fas', 'w') as f:
    for i in range(len(haps)):
        f.write(f'>{sequence_IDs[i]}\n{haps[i]}\n')


with open('RF_sigma1.trees', 'w') as f:
   for t in sts.trees():
      print(t.newick(), file=f)

end()