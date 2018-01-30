#install.packages('expm')
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

#comment
### JC Rate Matrix

Q<- matrix( 
c(-0.75, 0.25, 0.25, 0.25, 0.25, -0.75,0.25,0.25,0.25,0.25,-0.75,0.25,0.25,0.25,0.25,-0.75), 
nrow=4, 
ncol=4)

## JC equilibrium matrix

R<- matrix( 
  c(0.25, 0.25, 0.25, 0.25, 0.25, 0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25), 
  nrow=4, 
  ncol=4)

### Exponentiate the Rate Matrix

library(expm)
# Q is the rate matrix, mu is the transition rate, t is the time taken
t <- 0.5
mu <- 1



P <- expm(Q*(mu*t))
P 
# P is the matrix of transition probabilities


####################################
# 2 sequences
####################################

#1Pop directory needed

# read in the fasta sequences
S <- read.phyDat("test.fa",format="fasta", type="DNA")
#read in the fsc tree
T <-read.nexus("1PopDNA_1_true_trees.trees")
#pull out an example tree
T <- T$NumGen_tree_1_1_pos_0

#examine the tree
plot(T)
#output a list of the edges
T$edge


#rename the sequences to match the tree nodes

n<- length(S)
names<- c(rep(1:n))
for( i in 1:n) { 
  names[i] <- paste( i, 1, sep = ".")
}

names(S) <- names

# Setting initial parameters

mu <- 1
t<- 0.4 #T$edge.length[1]
Q1<- expm(Q*(mu*t)) #for branch length 1
Q2<- expm(Q*(mu*t)) #for branch length 2 (same in this case)

#Create L vector with number of rows equal to sites, to store likelihoods
Ct1 <- matrix(nrow = length(S$`1.1`), ncol = 4)


for(i in 1:length(S$`1.1`)){
  for(j in 1:4){
  tmp <- matrix( 
    c(0,0,0,0), 
    nrow=4, 
    ncol=1)
  tmp[S$`1.1`[i],1]<-1
  Ct1[i,j] <- Q1[j,]%*%tmp # Stores the P_ij(t)L_j(1) values (Sequence)
}
}

Ct1



################################################################
#Repeat for Sequence 2
################################################################


#Create L vector with number of rows equal to sites, to store likelihoods
Ct2 <- matrix(nrow = length(S$`1.1`), ncol = 4)


for(i in 1:length(S$`1.1`)){
  for(j in 1:4){
    tmp <- matrix( 
      c(0,0,0,0), 
      nrow=4, 
      ncol=1)
    tmp[S$`2.1`[i],1]<-1
    Ct2[i,j] <- Q2[j,]%*%tmp # Stores the P_ij(t)L_j(1) values (Sequence)
  }
}

Ct2



#Conditional Likelihoods at the node just the product of the 2 sequences (at each nucleotide)
LT<- Ct1*Ct2
LT
#Add the rows and multiply by 0.25 to get the full likelihoods at each site
LT <- rowSums (LT, na.rm = FALSE, dims = 1)
LT <- LT * 0.25



#Multiply the likelihoods to get the full likelihood
L <- prod(LT)
L

#log likelihood at Node 5
LL <- log(L)
LL



Likelihood<- function(par){


######################################################################
######## 3 sequence example
######################################################################

# read in the fasta sequences
S3 <- read.phyDat("test.fa",format="fasta", type="DNA")
#read in the fsc tree
T3 <-read.nexus("3seq_1_true_trees.trees")
#pull out an example tree
T3 <- T3$NumGen_tree_1_1_pos_0



plot(T3)
T3$edge

#rename the sequences to match the tree nodes



n<- length(S3)
names<- c(rep(1:3))
for( i in 1:3) { 
  names[i] <- paste( i, 1, sep = ".")
} 

names(S3) <- names

# reorder tree so it is postorder traversal (for bigger sequences)

T3 <- reorder(T3, "postorder")
#check postorder
T3$edge[1,]

#check the plots of postorder Trees

Tr <- ggtree(T3)
# Add Scale
Tr <- Tr + geom_treescale()
# Internal Node Numbers
Tr <- Tr +geom_text2(aes(subset=!isTip, label=node), hjust=-.3) 
# Tip Labels
Tr <- Tr + geom_tiplab()
Tr








# Initialised Parameters

mu <- 1
#Branch Lenghts
#t1<- 0.4 #T$edge.length[1]
#t2<- 0.6
t3<- par[1]+par[2] #constrain the edge lengths

Q1<- expm(Q*(mu*par[1]))
#Same Length clade
Q2<- expm(Q*(mu*par[1]))
Q3<- expm(Q*(mu*par[2]))
Q4<- expm(Q*(mu*t3))


################################################################
# first clade - to Calculate the likelihood at node 5
################################################################

#Create L vector with number of rows equal to sites, to store likelihoods

n<- length(S3[[1]])

Ct1 <- matrix(nrow = n, ncol = 4)


for(i in 1:n){
  for(j in 1:4){
    tmp <- matrix( 
      c(0,0,0,0), 
      nrow=4, 
      ncol=1)
    tmp[S3[[1]][i],1]<-1
    Ct1[i,j] <- Q1[j,]%*%tmp # Stores the P_ij(t)L_j(1) values (Sequence)
  }
}

#conditional likelihoods over each site (sequence 1) and each nucleotide possibility
Ct1




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
return(LL)
}


# works
Likelihood(c(0.1,0.3))


#Optimising PAckage
install.packages("optimr")
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
  
  
  
  
  
  
  
  
  
  
  
  




### JC stationary distribribtuion matrix

Pi<- matrix( 
  c(0.25, 0, 0, 0, 0, 0.25,0,0,0,0,0.25,0,0,0,0,0.25), 
  nrow=4, 
  ncol=4)






#traversal example

library(ape)
set.seed(1)
x = rtree(5,FALSE)
plot(x, show.tip.label=FALSE)
nodelabels(); # shows internal node numbers
tiplabels();
# the original setting just counts from the bottom up
x$edge
#postorder traversal
x = reorder(x, "postorder")
x$edge
# applies the clade at a time style order








# Junk Code

# Let A denoted by 1, C 2, G 3, T 4 

# transition probabilities function for JC model
transition_prob <- function(i,j,t,u){
  if(i==j){ output <- exp(-u*t) + (1-exp(-u*t))*0.25    # note for future: need == operator
                                 } else{
output <- (1-exp(-u*t))*0.25
  }
  return(output)
}

transition_prob(1,3,1,1)
A <- eigen(Q)$vectors 
D <- diag(eigen(Q)$values)
A %*% D %*% solve(A)   # Q= ADA^-1 

A %*% exp(D) %*% solve(A) 
#exp(Qmut) = exp(mu t) A %*%exp(d)%*%A-1

P_mat<- exp(mu*t)*A%*%exp(D)%*%solve(A)


A
D
A%*%exp(D)%*%solve(A)

}

