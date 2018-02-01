# This file will attempt to put together 2 and 3 sequence Likelihoods to Try and do what Ben
# Recommended on 31-1

#####################################################################
######## Packages
#####################################################################


# This Package is Required to Exponentiate the Rate Matrix 
#install.packages('expm')
library(expm)

# These Packages are Required for Manipulating the Tree
library(tree)
library(ape)
library(phangorn)
library(seqinr)
library("Biostrings")
library("ggplot2")
library("ggtree")

#####################################################################
######## Read in the Data
#####################################################################

# Directory is in root_node/ Vac_Scholarship/ Constrainted.../Data 

# read in the 3 fasta sequences
S3 <- read.phyDat("XXY.fa",format="fasta", type="DNA")

## Extract the first 2 for the 2 sequence tree 
S2 <- S3[c(1,2)]

#read in the 3 sequence tree
T3 <-read.nexus("T3.trees")
#pull out an example tree
T3 <- T3$NumGen_tree_1_1_pos_0
#read in the 2 sequence tree
T2 <-read.nexus("T2.trees")
#If there are multiple trees, pull out an example tree
T2 <- T2$NumGen_tree_1_1_pos_0

# example tree plots
plot(T3) 
plot(T2)

#rename the sequences to match the tree node labels
#S2
n<- length(S2)
names<- c(rep(1:n))
for( i in 1:n) { 
  names[i] <- paste( i, 1, sep = ".")
}
names(S2) <- names
#S3
n<- length(S3)
names<- c(rep(1:n))
for( i in 1:n) { 
  names[i] <- paste( i, 1, sep = ".")
}
names(S3) <- names

########################################################################
# 2 sequences
########################################################################

# t edge length

Likelihood2 <- function(t){
  mu <- 1
  #t<- 1
  
  Q<- matrix( 
    c(-0.75, 0.25, 0.25, 0.25, 0.25, -0.75,0.25,0.25,0.25,0.25,-0.75,0.25,0.25,0.25,0.25,-0.75), 
    nrow=4, 
    ncol=4)
  
  
  Q1<- expm(Q*(mu*t)) #Pij for branch length 1
  Q2<- expm(Q*(mu*t)) #Pij for branch length 2 (same in this case)
  
  #Create vector with number of cols equal to number of polymorphic sites, to store likelihoods
  n <- length(S2[[1]])
  Ct1 <- matrix( nrow = 4, ncol = n)
  
  # Loop to Calculate the Conditional Likelihoods for each nucletoide, for each site
  # i loops the different sites (filling across the cols)
  # 4 rows, ACGT
  
  
  for(i in 1:n){
    tmp <- matrix( c(0,0,0,0), nrow=4, ncol=1)     
    tmp[S2[[1]][i],1]<-1                             
    Ct1[,i] <- tmp                        
  }
  
 #conditional likelihoods of seq 1 no. rows correspond to ACGT, no. columns correspond to length of sequence 
  Ct1
  
################################################################
#Repeat for Sequence 2
################################################################
  
#conditional likelihoods of seq 2 no. rows correspond to ACGT, no. columns correspond to length of sequence 
  Ct2 <- matrix( nrow = 4, ncol = n)
  
  
  for(i in 1:n){
    tmp <- matrix( c(0,0,0,0), nrow=4, ncol=1)      
    tmp[S2[[2]][i],1]<-1                            
    Ct2[,i] <- tmp                        
  }
  
  
  Ct2
  
## Matrix Multiply Through the Pij to get conditional likelihoods of ACGT (rows) at each site (cols) with PIJ
  
  Ct1 <- Q1%*%Ct1
  Ct2 <- Q2%*%Ct2
  
  
## Conditional Likelihoods for each ACGT at the node just the product of the 2 sequences (sites down the cols)
# Elementwise multiply
  CL <- Ct1 *Ct2 
  CL
  
## Add the columns and multiply by 0.25 to get the full likelihoods at each site
  L <- colSums (CL, na.rm = FALSE, dims = 1)
  L <- L * 0.25
  
#Multiply the likelihoods over all sites to get the full likelihood
  L <- prod(L)
  L
  
  #log likelihood 
  LL <- log(L)
  LL
  return(-LL)
}

######################################################################
######## 3 sequence Likelihood Function
######################################################################

# t1 is a scalar and what we want to maximise
# t2 is a vector of 2 elements (t' and t'')
# t2's sum matters for likelihood, not the individual elements

Likelihood3<- function(t){
  
  ####################  
  # Initial Parameters
  ####################
  
  Q<- matrix( 
    c(-0.75, 0.25, 0.25, 0.25, 0.25, -0.75,0.25,0.25,0.25,0.25,-0.75,0.25,0.25,0.25,0.25,-0.75), 
    nrow=4, 
    ncol=4)

  # Mutation Rate
  mu <- 1
  # Branch Lenghts
  #t <- c(0,0,0)
  #t1<- 0.4 #T$edge.length[1]
  #t2<- 0.6

  # Exponentiated Rate Matrices for Likelihood Calculations
  # Minor Clade
  library(expm)
  Q1<- expm(Q*(mu*t[1]))
  Q2<- expm(Q*(mu*t[1]))
  # Connect Minor Clade to Root
  Q3<- expm(Q*(mu*t[2]))
  # Outgroup Branch
  Q4<- expm(Q*(mu*t[3]))
  
  
  ################################################################
  # minor clade - to Calculate the likelihood at node 5 (ingorup node)
  ################################################################
  
  #Create vector with number of rows equal to ACGT, cols for number of sites, to store likelihoods
  
  n<- length(S3[[1]])
  Ct1 <- matrix(nrow = 4, ncol = n)
  
  # Loops all the conditional likelihoods for seq 1 along the cols
  
  for(i in 1:n){
    tmp <- matrix( c(0,0,0,0), nrow=4, ncol=1)      # creates empty column matrix 
    tmp[S3[[1]][i],1]<-1                             # Inputs 1 into column where there is the nucleotide
    Ct1[,i] <- tmp                        #L_j values (Sequence)
  }
  
  #conditional likelihoods of seq 1 no. rows correspond to ACGT, no. columns correspond to length of sequence 
  Ct1
  
  
  ################################################################
  #Repeat for Sequence 2
  ################################################################

  Ct2 <- matrix(nrow = 4, ncol = n)
  
  for(i in 1:n){
    tmp <- matrix( c(0,0,0,0), nrow=4, ncol=1)      
    tmp[S3[[2]][i],1]<-1                             
    Ct2[,i] <- tmp                        
  }
  
  # conditional likelihoods at sequence 2 for the 4 nucleotides
  Ct2
  
  
  ################################################################
  #Putting It Together
  ################################################################
  
  ## Matrix Multiply Through the Pij to get conditional likelihoods of ACGT (rows) at each site (cols) with PIJ
  
  Ct1 <- Q1%*%Ct1
  Ct2 <- Q2%*%Ct2
  
  
  ## Conditional Likelihoods for each ACGT at the node just the product of the 2 sequences (sites down the cols)
  # Elementwise multiply to get conditional likelihood at node  
  Ct4 <- Ct1*Ct2 
  ############################################################   
  # These are the Likelihoods at the nodes for Later Calcs
  ############################################################   
  Ct4


  ################################################################
  # the big clade ( i.e. the likelihood at root)
  ################################################################
  
  ############################
  ### Sequence 3
  ############################
  
  
  #Create L vector with number of rows equal to sites, to store likelihoods
  Ct3 <- matrix(nrow = 4, ncol = n)
  
  
  for(i in 1:n){
    tmp <- matrix( c(0,0,0,0), nrow=4, ncol=1)     
    tmp[S3[[1]][i],1]<-1                             
    Ct3[,i] <- tmp                        
  }
  
  # conditional likelihoods at sequence 3 for the 4 nucleotides
  Ct3
  
  ################################################################
  #Putting It Together
  ################################################################
  
  
  ## Matrix Multiply Through the Pij to get conditional likelihoods of ACGT (rows) at each site (cols) with PIJ
  
  Ct3 <- Q3%*%Ct3
  Ct4 <- Q4%*%Ct4
  
  #product the two matrices to create conditional likelihoods 4 nucleotides over all sites
  
  CL<- Ct3*Ct4
  CL
  
  ## Add the columns and multiply by 0.25 to get the full likelihoods at each site
  L <- colSums (CL, na.rm = FALSE, dims = 1)
  L <- L * 0.25
  
  #Multiply the likelihoods over all sites to get the full likelihood
  L <- prod(L)
  L
  
  #log likelihood 
  LL <- log(L)
  LL

  return(-LL)
}


#install.packages("optimr")
library(optimr)

#is a minimising optimiser unless you add control$maximize=TRUE
optim(par =c(1,0.2,1.2),  fn=Likelihood3, lower=0.0000000000000001,method="L-BFGS-B")

optim(par = 1, fn=Likelihood2, lower=0.0000000000000001,method="L-BFGS-B")

Likelihood3(c(1,0.2,1.2))





n <- 1000
Data <- matrix( nrow = n, ncol = 2)
Data[,1] <- c(1:1000)*0.01
for(i in 1:n){
  Data[i,2] <- Likelihood3(c(Data[i,1],0.001,0.001))
}
Data <- as.data.frame(Data)
colnames(Data)<- c("t","Likelihood")

ggplot(data = Data, aes(x = t, y = Likelihood)) + geom_point(shape=1)
# geom_hline(yintercept = 0, colour = "gray65") +
# geom_vline(xintercept = 0, colour = "gray65") +
ggtitle("t against L")





S <- S3[c(1,2)]
S

