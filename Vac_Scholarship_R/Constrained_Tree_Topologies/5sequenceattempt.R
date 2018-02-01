######################################################################
######## 5 sequence example
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
S5 <- read.phyDat("test.fa",format="fasta", type="DNA")
#read in the fsc tree (only 1 )
T5 <-read.nexus("5seq_1_true_trees.trees")
T5 <- as.phylo(T5)
T5

plot(T5)
T5$edge


#rename the sequences to match the tree nodes 


n<- length(S5)
names<- c(rep(1:n))
for( i in 1:n) { 
  names[i] <- paste( i, 1, sep = ".")
} 

names(S5) <- names



# reorder tree so it is postorder traversal (for bigger sequences)

T5 <- reorder(T5, "postorder")
#check postorder
T5$edge[1,]

#rearranges S5 so that the order of the sequences matches the tree
S5 <- S5[order(T5$tip.label)]


#Convert Tree to ggtree format (more data)

Tr <- ggtree(T5)
# Add Scale and node numbers 
Tr <- Tr + geom_treescale() + geom_text(aes(label=node), hjust=-.3)
Tr

# Internal Node Numbers
# Tr <- Tr +geom_text2(aes(subset=!isTip, label=node), hjust=-.3) 
# Tip Labels
#Tr <- Tr + geom_tiplab() + geom_nodepoint()
#Tr





par <- T5$edge.length


  
################################################################
#Calculate the conditional likelihood matrices of the tips
################################################################
  


T5$edge.length
T5$edge


###############
# Creating the Conditional Likelihood Matrices 
###############

#sequence length
ns <- length(S5[[1]])
#number of edges
ne <- length(T5$edge.length)
#number of nodes
nn <- length(Tr$data$node)

#creates number of blank conditional likelihood matirces for each edge
Ct<- replicate(nn, matrix(nrow = 4, ncol = ns), simplify=FALSE)
Ct
Ct[[1]][2,3]
# if it is a lower node, we calculate as in the 2 sequence example





#identical( Edge.Length[1], Edge.Length[2] )

Likelihood5<- function(t){

# Set a bunch of iniailising edge lengths as a vector, length(T5$edge.length) (t1-t8 in this case)
#t<- c(0,2,3,4,2,6,7,9)
Edge <- T5$edge
Edge.Length <- t







# edge lengths correspond to the edges


Q<- matrix( 
  c(-0.75, 0.25, 0.25, 0.25, 0.25, -0.75,0.25,0.25,0.25,0.25,-0.75,0.25,0.25,0.25,0.25,-0.75), 
  nrow=4, 
  ncol=4)

mu <- 1 
Edge.Length
Edge

Edges <- cbind(Edge,Edge.Length)
# produces the tip node conditional likelihood matrices
Edges

for(i in 1:ne){
if(Edges[i,2]<=length(S5)){
 for(j in 1:ns){
      tmp <- matrix( c(0,0,0,0), nrow=4, ncol=1)      # creates empty column matrix 
      tmp[S5[[Edges[i,2]]][j],1] <- 1                             # Inputs 1 into column where there is the nucleotide
      Ct[[Edges[i,2]]][,j] <- tmp                        # Stores the L_j values (Sequence)
    }
  }   
}
Ct

Edge

Edges
Ct
# to fill in the internal nodes
for(i in 1:ne){
  if(Edges[i,2] > length(S5)){ 
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
Ct

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
  
  
Likelihood5(c(0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01))
Edge.Length
Edges
#Optimising Pαckage
#install.packages("optimr")
library(optimr)

#is a minimising optimiser unless you add control$maximize=TRUE
optim(par = c(0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01), fn=Likelihood5, lower=0.00000001,method="L-BFGS-B")

# need to make the lower bound just above 0 or else optimising fails

# Unfortunately Minimises to 0




  #conditional likelihoods over each site (sequence 1) and each nucleotide possibility

  
  
  
  
  ################################################################
  #Repeat for Sequence 2
  ################################################################
  
  
  #Create L vector with number of rows equal to sites, to store likelihoods
  Ct2 <- matrix(nrow = n, ncol = 4)
  
  
  for(i in 1:n){
    for(j in 1:4){
      tmp <- matrix( 
        c(0,0,0,0), 
        nrow=4, 
        ncol=1)
      tmp[S3[[2]][i],1]<-1
      Ct2[i,j] <- Q2[j,]%*%tmp # Stores the P_ij(t)L_j(1) values (Sequence)
    }
  }
  
  # conditional likelihoods at sequence 2 for the 4 nucleotides
  Ct2
  
  
  ################################################################
  #Putting It Together
  ################################################################
  
  #product the two matrices to create conditional likelihoods 4 nucleotides over all sites
  
  
  LM5<- Ct1 *Ct2
  LM5
  
  #Add the rows and multiply by 0.25 to get the full likelihoods at each site
  W <- rowSums (LM5, na.rm = FALSE, dims = 1)
  W <- W * 0.25
  
  #Multiply the likelihoods to get the full likelihood
  L <- prod(W)
  L
  
  #log likelihood at Node 5
  LL <- log(L)
  LL
  
  
  
  ################################################################
  # the big clade now ( i.e. the likelihood at node 4)
  ################################################################
  
  
  ############################
  ### Sequence 3
  ############################
  
  
  
  #Create L vector with number of rows equal to sites, to store likelihoods
  Ct3 <- matrix(nrow = n, ncol = 4)
  
  
  for(i in 1:n){
    for(j in 1:4){
      tmp <- matrix( 
        c(0,0,0,0), 
        nrow=4, 
        ncol=1)
      tmp[S3[[3]][i],1]<-1
      Ct3[i,j] <- Q4[j,]%*%tmp # Stores the P_ij(t)L_j(1) values (Sequence)
    }
  }
  
  # conditional likelihoods (multiplied by P_ij) at sequence 3 for the 4 nucleotides
  Ct3
  
  # Conditional Likelihoods at node 5 for the 4 nucletoides over the 16 sites (from previous)
  LM5
  # to introduce the P_ij
  Ct4<- t(Q3 %*%t(LM5))
  Ct4
  
  
  
  
  ################################################################
  #Putting It Together
  ################################################################
  
  #product the two matrices to create conditional likelihoods 4 nucleotides over all sites
  
  
  LM4<- Ct3 *Ct4
  LM4
  
  #Add the rows and multiply by 0.25 to get the full likelihoods at each site
  W <- rowSums (LM4, na.rm = FALSE, dims = 1)
  W <- W * 0.25
  W
  #Multiply the likelihoods to get the full likelihood
  L <- prod(W)
  L
  
  #log likelihood at Node 4
  LL <- log(L)
  LL
  
  para <- c(LL,par[1],par[2],t3)
  bar <- c("LogLikelihood", "t1", "t2","t3")
  newList <- list("integer" = para, "names" = bar)
  
  
  return(newList)





# Checking it works
Likelihood(c(0.1,0.3),1,S3,T3)


#Optimising Pαckage
#install.packages("optimr")
library(optimr)

optimr(par = c(0.3, 0.4),fn=Likelihood)


# Shoots off to infinity as t1 goes to 0
Likelihood(c(0.00000000000000000000001,0.5))











# comment out junk code
.f = function() {
  
  
  
  
  #sum the rows and multiply by 0.25 to create full likelihoods for the 4 nucleotides
  LL3 <- colSums (Ct1, na.rm = FALSE, dims = 1)
  
  
  
  ## Conditional likelihoods at sequence 3 for the 4 nucleotides (Including P_ij)
  LL3
  
  ## Conditional Likelihoods at Node 5 for the 4 Nucleotides
  L5
  # Multiple L5 By the ( P_ij) to get the 4 conditional likelihoods at node 5
  L5<- Q3%*%L5
  L5
  
  # get the 4 conditional likelihoods at node 4 (multiply elementwise)
  
  LC<- LL3*L5
  LC
  #Total Likelihood is the sum multiplied by 0.25
  L <-sum(LC)*0.25
  L
  #Log Likelihood
  LL <- log(L)
  LL
  
  