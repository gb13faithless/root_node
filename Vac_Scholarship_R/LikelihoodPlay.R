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
  Ct1[i,j] <- Q1[i,j]%*%tmp[j,1] # Stores the P_ij(t)L_j(1) values (Sequence)
}
}

Ct1
#sum the columns to create conditional likelihoods for the 4 nucleotides
LL1 <- colSums (Ct1, na.rm = FALSE, dims = 1)
#sums the conditional likelihoods
LL1<- sum(LL1)


################################################################
#Repeat for Sequence 2
################################################################


#Create L vector with number of rows equal to sites, to store likelihoods
Ct1 <- matrix(nrow = length(S$`1.1`), ncol = 4)


for(i in 1:length(S$`1.1`)){
  for(j in 1:4){
    tmp <- matrix( 
      c(0,0,0,0), 
      nrow=4, 
      ncol=1)
    tmp[S$`1.1`[i],1]<-1
    Ct1[i,j] <- Q1[i,j]%*%tmp[j,1] # Stores the P_ij(t)L_j(1) values (Sequence)
  }
}

Ct1
#sum the columns to create conditional likelihoods for the 4 nucleotides
LL1 <- colSums (Ct1, na.rm = FALSE, dims = 1)
#sums the conditional likelihoods
LL1<- sum(LL1)


################################################################
#Repeat for Sequence 2
################################################################

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
    Ct1[i,j] <- Q1[i,j]%*%tmp[j,1] # Stores the P_ij(t)L_j(1) values (Sequence)
  }
}

Ct1
#sum the columns to create conditional likelihoods for the 4 nucleotides
LL1 <- colSums (Ct1, na.rm = FALSE, dims = 1)
#Conditional likelihood at sequence 1 for the 4 nucleotides
LL1



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
    Ct2[i,j] <- Q2[i,j]%*%tmp[j,1] # Stores the P_ij(t)L_j(1) values (Sequence)
  }
}

Ct2
#sum the columns to create conditional likelihoods at sequence 2 for the 4 nucleotides
LL2 <- colSums (Ct2, na.rm = FALSE, dims = 1)
#Conditional likelihoods at sequence 2 for the 4 nucleotides
LL2



#Conditional Likelihoods at the node just the product of the 2 sequences (at each nucleotide)
LT<- LL1*LL2
#Total Likelihood is the sum multiplied by 0.25 (JC)
LT <- sum(LT)*0.25
#log likelihood
LLT<- log(LT)
LLT


###################################
######## 3 sequence example
###################################

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

#reorder tree so it is postorder traversal

T3 <- reorder(T3, "postorder")
#check postorder
T3$edge[1,]

#check the plot of postorder

Tr <- ggtree(T3)
# Add Scale
Tr <- Tr + geom_treescale()
# Internal Node Numbers
Tr <- Tr +geom_text2(aes(subset=!isTip, label=node), hjust=-.3) 
# Tip Labels
Tr <- Tr + geom_tiplab()
Tr


# conditional likelihoods

mu <- 1
t1<- 0.4 #T$edge.length[1]
t2<- 0.6
t2<- 0.8
Q1<- expm(Q*(mu*t1))
#T$edge.length[1]
Q2<- expm(Q*(mu*t1))
Q3<- expm(Q*(mu*t2))
Q4<- expm(Q*(mu*t3))


# first clade - likelihood at node 5

# set up likelihoods for all 16 sites


L <- matrix(, nrow = length(S3$`1.1`), ncol = 4)

#need 4 conditional likelihoods for node 4 (ACGT)


for(i in 1:4){
for(j in 1:length(S3[[1]])){
  #creates conditional likelihood (in this case 0 and 1)
  tmp <- matrix( 
    c(0,0,0,0), 
    nrow=4, 
    ncol=1)
  # edge order
  t<- T3$edge[1,2]
  # runs through the sequence
  tmp[S3[[t]][j],1]<-1
  #
  Ct1<- Q1[i,]%*%tmp
  #creates conditional likelihood (in this case 0 and 1)
  tmp <- matrix( 
    c(0,0,0,0), 
    nrow=4, 
    ncol=1)
  t<- T3$edge[2,2]
  tmp[S3[[t]][j],1]<-1
  Ct2<- Q2[i,]%*%tmp
  L[j,i] <- Ct1%*%Ct2
}
}
L


#4 conditional likelihoods at each site (rename)
L5<-L
#log
LL5 <- log(L5)
LL5
#log likelihood totals for all sites (just sum them up)
total <- colSums (LL5, na.rm = FALSE, dims = 1)
#node 5 conditional log likelihoods
total
#total conditional likelihoods at each nucleotide (ACGT) (unlogged)
exp(total)


# the big clade now ( i.e. the likelihood at node 4)


L <- matrix(, nrow = length(S3$`1.1`), ncol = 4)

for(
  i in 1:4){
  for(j in 1:length(S3[[1]])){
    #creates conditional likelihood for each site 
    tmp <- matrix( 
      c(0,0,0,0), 
      nrow=4, 
      ncol=1)
    # runs through the sequence 
    tmp[S3[[3]][j],1]<-1
    # Calculates the conditional likelihood component by sequence 3
    Ct3<- Q3[i,]%*%tmp
    # calculates the conditional likelihood for this site by nucleotide i (ACGT) at site j (1-length of the sequence)
    tmp2 <- Q4[i,]%*%L5[j,]
    #creates conditional likelihood at site j for nucleotide i 
    L[j,i] <- tmp2%*%Ct3
  }
}

L

#log likelihoods
LL <- log(L)
LL
# sum to get conditional likelihoods
total <- colSums (LL, na.rm = FALSE, dims = 1)
#node 4 conditional log likelihoods
total
#total conditional likelihoods (unlogged)
exp(total)
#sum to get total
L<- sum(exp(total))
#multiply by 0.25 for pi0 (the same in JV)
L<- 0.25*L
#log likelihood
LL <- log(L)
LL





















# comment out junk code
.f = function() {




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

