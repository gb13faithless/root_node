######################################################################
######## 4 sequence example
######################################################################

# These Packages are Required for Manipulating the Tree
library(tree)
library(ape)
library(phangorn)
library(seqinr)
library(Biostrings)
library(ggplot2)
library(ggtree)

# This Package is Required to Exponentiate the Rate Matrix 
#install.packages('expm')
library(expm)


# read in the fasta sequences
S35 <- read.phyDat("test.fa",format="fasta", type="DNA")
#read in the fsc tree (only 1 )
T35 <-read.nexus("35_1_true_trees.trees")
T35 <-T35$NumGen_tree_1_1_pos_0
T35 <- as.phylo(T35)
T35

plot(T35)
T35$edge


#rename the sequences to match the tree nodes 


n<- length(S35)
names<- c(rep(1:n))
for( i in 1:n) { 
  names[i] <- paste( i, 1, sep = ".")
} 

names(S35) <- names



# reorder tree so it is postorder traversal (for bigger sequences)

T35 <- reorder(T35, "postorder")
#check postorder
T35$edge

#rearranges S35 so that the order of the sequences matches the tree
S35 <- S35[order(T35$tip.label)]


#Convert Tree to ggtree format (more data)

Tr <- ggtree(T35)
# Add Scale and node numbers 
Tr <- Tr + geom_treescale() + geom_text(aes(label=node), hjust=-.3)
Tr

# Internal Node Numbers
# Tr <- Tr +geom_text2(aes(subset=!isTip, label=node), hjust=-.3) 
# Tip Labels
#Tr <- Tr + geom_tiplab() + geom_nodepoint()
#Tr



################################################################
#Calculate the conditional likelihood matrices of the tips
################################################################



T35$edge.length
T35$edge


###############
# Creating the Conditional Likelihood Matrices 
###############

#sequence length
ns <- length(S35[[1]])
#number of edges
ne <- length(T35$edge.length)
#number of nodes
nn <- length(Tr$data$node)

#creates number of blank conditional likelihood matirces for each edge
Ct<- replicate(nn, matrix(nrow = 4, ncol = ns), simplify=FALSE)
Ct





# 6 edge lengths
#t1,t1,t2,t3,t',t''
# t'+t'' depends the maximisation, not each separately

t<- rep(1,ne)
t



Likelihood35<- function(t){
  
  #t<- c(1,1,2,3,1,1)
  
  # to ensure that t1=t1 over the first 2 branch lengths
  #t[1]<- t[2]
  Edge <- T35$edge
  Edge.Length <- t
  
  
  
  
  
  
  
  # edge lengths correspond to the edges
  
  
  Q<- matrix( 
    c(-0.75, 0.25, 0.25, 0.25, 0.25, -0.75,0.25,0.25,0.25,0.25,-0.75,0.25,0.25,0.25,0.25,-0.75), 
    nrow=4, 
    ncol=4)
  
  mu <- 1 
  
  Edges <- cbind(Edge,Edge.Length)
  Edges
  # produces the tip node conditional likelihood matrices
  for(i in 1:ne){
    if(Edges[i,2]<=length(S35)){
      for(j in 1:ns){
        tmp <- matrix( c(0,0,0,0), nrow=4, ncol=1)      # creates empty column matrix 
        tmp[S35[[Edges[i,2]]][j],1] <- 1                             # Inputs 1 into column where there is the nucleotide
        Ct[[Edges[i,2]]][,j] <- tmp                        # Stores the L_j values (Sequence)
      }
    }   
  }
  
  
  # to fill in the internal nodes
  for(i in 1:ne){
    if(Edges[i,2] > length(S35)){ 
      tmp <- Edges[Edges[,1] == Edges[i,2],]
      # This is the step I need to add in the Pij matrix
      M1<- expm(Q*mu*tmp[1,3]) %*% Ct[[tmp[1,2]]]
      M2<- expm(Q*mu*tmp[2,3]) %*% Ct[[tmp[2,2]]]
      Ct[[tmp[1,1]]] <- M1 * M2
    }
  }
  
  
  #to to the root node
  tmp <- tail(Edges, 2)
  # This is the step I need to add in the Pij matrix
  M1<- expm(Q*mu*tmp[1,3]) %*% Ct[[tmp[1,2]]]
  M2<- expm(Q*mu*tmp[2,3]) %*% Ct[[tmp[2,2]]]
  Ct[[tmp[1,1]]] <- M1 * M2
  
  
  ##
  ## Add the columns and multiply by 0.25 to get the full likelihoods at each site
  L <- colSums (Ct[[tmp[1,1]]], na.rm = FALSE, dims = 1)
  L <- L * 0.25
  
  #Log the likelihoods and add
  LL <- log(L)
  LL <- sum(LL)

  
  #regular likelihoods by exponentiating
  #L <- exp(L)
  #L
  
  return(-LL)
}


Likelihood35(t)


#Optimising Pαckage
#install.packages("optimr")
library(optimr)

#is a minimising optimiser unless you add control$maximize=TRUE


out <- optim(par = t, fn=Likelihood35, lower=0.00000001,method="L-BFGS-B")


out





# need to make the lower bound just above 0 or else optimising fails

# Unfortunately Minimises to 0



