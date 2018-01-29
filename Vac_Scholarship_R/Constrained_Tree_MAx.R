
library(tree)
library(ape)
library(phangorn)
library(seqinr)
library("phangorn")
library("ape")
library("Biostrings")
library("ggplot2")
library("ggtree")
library("devtools")
library("strataG")
library("phylobase")

# read in the fasta sequences
S <- read.phyDat("test.fa",format="fasta", type="DNA")
#read in the fsc tree
T <-read.nexus("1PopDNA_1_true_trees.trees")
#pull out an example tree
T <- T$NumGen_tree_1_1_pos_0
plot(T)



n<- length(S)
names<- c(rep(1:2))
for( i in 1:2) { 
  names[i] <- paste( i, 1, sep = ".")
} 

names(S) <- names

fit = pml(T, data=S, rearrangment="none",rate = 0.00002)
fit




library(ape)






tree$edge.length[1]
t1 <- tree$edge.length[1]
t2 <- tree$edge.length[2]
t3 <- tree$edge.length[3]

# Let A denoted by 1, C 2, G 3, T 4 

# transition probabilities function for JC model
transition_prob <- function(i,j,t,u){
  if(i==j){ output <- exp(-u*t) + (1-exp(-u*t))*0.25    # note for future: need == operator
  } else{
    output <- (1-exp(-u*t))*0.25
  }
  return(output)
}



transition_prob(1,1,t1,1)
# all transition equilibria are 0.25 
# mutation rate = 1

0.25*transition_prob(1,4,t1,1)



getAnywhere('pml')
getAnywhere('optim.pml')
getAnywhere('pml.fit')

fit = pml(T, data=S, rearrangment="none",rate = 0.00002)


tree <- reorder(tree, "postorder")
tree












#ggtree

Tr <- ggtree(tree)
Tr
# Add Scale
Tr <- Tr + geom_treescale()

# Internal Node Numbers
Tr <- Tr +geom_text2(aes(subset=!isTip, label=node), hjust=-.3) 
# Tip Labels
Tr <- Tr + geom_tiplab()
# Tip Points
Tr <- Tr + geom_tippoint()
# Tip Points
Tr +geom_text2(aes(subset=!isTip, label=node), hjust=1) 



