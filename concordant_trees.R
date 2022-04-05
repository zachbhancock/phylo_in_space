library(evobiR)
library(phytools)

#function that describes expected concordant trees
f1 <- function(t){
  return(1 - (2/3)*exp(-t))
}
#function that describes expected discordant trees
f2 <- function(t){
  return(exp(-t)/3)
}

collection <- "3tree_territorial_t25.trees"
ref <- "3tree_topologies.txt"
top.counts <- countTrees(collection, ref)

tc <- top.counts[[1]]
sum.tc <- sum(top.counts[[1]])

#proportion of trees
pr.tc <- tc / sum.tc

#true tau value - should be close to the expected; if not, may have to adjust 
t <- log(3*pr.tc[1])

#prints expected values of concordant and discordant given the true tau value
f3(t)
f4(t)