#Obtain SNPs in R

library(devtools)
install_github('coleoguy/evobir', build_vignettes = F)
library(evobiR)
collection <- "CA_WF_sigma05.trees"
ref <- "topologies.txt"
top.counts <- countTrees(collection, ref)
top.counts

#SNPs

mutated = msprime.mutate(ts, rate=1e-8, model=msprime.InfiniteSites(msprime.NUCLEOTIDES), random_seed=1)

ts_sp = mutated.simplify(A)
haps = []
for i in ts_sp.haplotypes():
    haps.append(i)

sequence_IDs = []
for i in range(len(haps)):
    sequence_IDs.append(f'sample_{ts_sp.samples()[i]}_pop_{ts_sp.node(i).population}')

with open('CA_WF_sigma25.fas', 'w') as f:
    for i in range(len(haps)):
        f.write(f'>{sequence_IDs[i]}\n{haps[i]}\n')
