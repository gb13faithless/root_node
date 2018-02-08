# This file holds the functions used to calculate the likelihood

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
library(Biostrings)
library(ggplot2)
library(ggtree)


########################################################################
# 2 sequences
########################################################################

# t edge length
# S2 is the list of sequences (phyDat format)
# Tt is the tree (phylo format)

Likelihood2 <- function(t,S2,T2){
  mu <- 1
  #t<- 1
  
  Q<- matrix( 
    c(-0.75, 0.25, 0.25, 0.25, 0.25, -0.75,0.25,0.25,0.25,0.25,-0.75,0.25,0.25,0.25,0.25,-0.75), 
    nrow=4, 
    ncol=4)
  
  
  Q1<- expm(Q*(mu*t)) #Pij for branch length 1
  Q2<- expm(Q*(mu*t)) #Pij for branch length 2 (same in this case)
  
  
  ################################################################
  #Sequence 1
  ################################################################
  
  
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
  
## Matrix Multiply Through the Pij to get conditional likelihoods of ACGT (rows) at each site (cols) with PIJ
  
  Ct1 <- Q1%*%Ct1
  Ct2 <- Q2%*%Ct2
  
  
## Conditional Likelihoods for each ACGT at the node just the product of the 2 sequences (sites down the cols)
# Elementwise multiply
  CL <- Ct1 *Ct2 
  
## Add the columns and multiply by 0.25 to get the full likelihoods at each site
  L <- colSums (CL, na.rm = FALSE, dims = 1)
  L <- L * 0.25
  
#Multiply the likelihoods over all sites to get the full likelihood
  L <- prod(L)
  
  #log likelihood 
  LL <- log(L)
  return(-LL)
}

######################################################################
######## 3 sequence Likelihood Function
######################################################################

# t is a vector with 3 elements
# t[1] is what we want to maximise
# t[2],t[3] are t' and t''
# t2's sum matters for likelihood, not the individual elements
# S3 is the list of sequences (phyDat format)
# Sequence 1 and 2 are the ingorup, 3 is the outgroup
# T3 is the tree (phylo format)

Likelihood3<- function(t,S3,T3){
  
  Q<- matrix( 
    c(-0.75, 0.25, 0.25, 0.25, 0.25, -0.75,0.25,0.25,0.25,0.25,-0.75,0.25,0.25,0.25,0.25,-0.75), 
    nrow=4, 
    ncol=4)

  # Mutation Rate
  mu <- 1

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
  
  
  ################################################################
  # Sequence 1 
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
  #Putting Sequence 1 and 2 Together
  ################################################################
  
  ## Matrix Multiply Through the Pij to get conditional likelihoods of ACGT (rows) at each site (cols) with PIJ
  
  Ct1 <- Q1%*%Ct1
  Ct2 <- Q2%*%Ct2
  
  
  ## Conditional Likelihoods for each ACGT at the node just the product of the 2 sequences (sites down the cols)
  # Elementwise multiply to get conditional likelihood at node  
  Ct4 <- Ct1*Ct2 
  
  ############################################################   
  # These are the Conditional Likelihoods at the node (for Later Calcs)
  ############################################################   


  ################################################################
  # Calculating the big clade ( i.e. the likelihood at root)
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
  
  # conditional likelihoods at sequence 3 for the 4 nucleotides is Ct3
  
  ################################################################
  #Putting It Together
  ################################################################
  
  ## Matrix Multiply Through the Pij to get conditional likelihoods of ACGT (rows) at each site (cols) with PIJ
  
  Ct3 <- Q3%*%Ct3
  Ct4 <- Q4%*%Ct4
  
  #product the two matrices to create conditional likelihoods 4 nucleotides over all sites
  
  CL<- Ct3*Ct4
  
  ## Add the columns and multiply by 0.25 to get the full likelihoods at each site
  L <- colSums (CL, na.rm = FALSE, dims = 1)
  L <- L * 0.25
  
  #Multiply the likelihoods over all sites to get the full likelihood
  L <- prod(L)
  
  #log likelihood 
  LL <- log(L)

  return(-LL)
}


#install.packages("optimr") to do the maximising
library(optimr)



################################################################
# Function to calculate the ML estimate for 2 sequence trees
################################################################

# t is a scalar (t1)
# S2 is the list of sequences (phyDat format)
# T2 is the tree (phylo format)
# l is the lower bound for branch lenghts (above 0)

ML2<- function(t,S2,T2,l){
  #is a minimising optimiser unless you add control$maximize=TRUE
  optim<- optim(par = t, S2=S2, T2=T2,  fn=Likelihood2, lower=l, method="L-BFGS-B")
  names <- c("t1")
  branch.length <- c(optim$par[1])
  df <- data.frame(names,branch.length)
  return(df)
}


################################################################
# Function to calculate the ML estimate for 3 sequence trees
################################################################

# t is a vector with 3 elements
# t[1] is what we want to maximise
# t[2],t[3] are t' and t''
# t2's sum matters for likelihood, not the individual elements
# S3 is the list of sequences (phyDat format)
# T3 is the tree (phylo format)
# l is the lower bound for branch lenghts (above 0)



ML3<- function(t,S3,T3,l){
#is a minimising optimiser unless you add control$maximize=TRUE
optim<- optim(par = t, S3=S3, T3=T3,  fn=Likelihood3, lower=l, method="L-BFGS-B")
names <- c("t1","t2")
branch.length <- c(optim$par[1],optim$par[2]+optim$par[3])
df <- data.frame(names,branch.length)
  return(df)
}

ML3(t3,Sequence[[1]],T1,0.00000001)




#function to get the minimum non negative number
minpositive = function(x) min(x[x > 0])




#junk code
function(n){
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
  
  
}









