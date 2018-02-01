

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

######################################################################
######## 3 sequence Likelihood Function
######################################################################


# par[1] is t1 (between the minor clade)
# par[2] is t2 (connecting the minor clade to the root)
# mu is the mutation rate (usually set as 1)
# S3 is the 3 sequence tree, of class phlyo
# T3 is the 3 sequences phylogenetic tree (read in previously as a nexus file) where sequences 1 and 2 
# are in a subclade, and all 3 sequences are named in the apporpriate format

#,mu,S3,T3 other parameters
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
  #t1<- 0.4 #T$edge.length[1]
 #t2<- 0.6
  
   #constrain the edge lengths
  
  
  # Exponentiated Rate Matrices for Likelihood Calculations
  #Minor Clade
  library(expm)
  Q1<- expm(Q*(mu*t[1]))
  Q2<- expm(Q*(mu*t[1]))
  #Connect Minor Clade to Root
  Q3<- expm(Q*(mu*t[2]))
  # Outgroup Branch
  Q4<- expm(Q*(mu*t[3]))
  
  
  ################################################################
  # minor clade - to Calculate the likelihood at node 5 (ingorup node)
  ################################################################
  
  #Create L vector with number of rows equal to sites, to store likelihoods
  
  n<- length(S3[[1]])
  
  Ct1 <- matrix(nrow = 4, ncol = n)
  

  
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
  
  
  #Create L vector with number of rows equal to sites, to store likelihoods
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
  # Elementwise multiply
  Ct4 <- Ct1*Ct2 
  ############################################################   
  # These are the Likelihoods at the nodes for Later Calcs
  ############################################################   
  Ct4
  
  ############################################################   
  # If we want to calulate the likelihood at 5 (for the subclade)
  ############################################################ 
  
  ##
  ## Add the columns and multiply by 0.25 to get the full likelihoods at each site
  #L <- colSums (Ct4, na.rm = FALSE, dims = 1)
  #L <- L * 0.25
  
  #Multiply the likelihoods over all sites to get the full likelihood
  #L <- prod(L)
  #L
  
  #log likelihood 
  #LL <- log(L)
  #LL
  
  
  
  ################################################################
  # the big clade now ( i.e. the likelihood at node 4)
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
  
 # para <- c(LL,t[1],t[2],t3)
  #bar <- c("LogLikelihood", "t1", "t2","t3")
  #newList <- list("integer" = para, "names" = bar)
  
  return(-LL)
}

Likelihood3(c(1,1,1))


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

#install.packages("optimr")
library(optimr)

#is a minimising optimiser unless you add control$maximize=TRUE
optim(par = c(4.76,1,1), fn=Likelihood3, lower=0.00001,method="L-BFGS-B")
optim(par = c(0.1), fn=Likelihood2, lower=0.00001,method="L-BFGS-B")



S <- S3[c(1,2)]
S

Likelihood2 <- function(t){
  mu <- 1
  #t<- 0.4 
  
  Q<- matrix( 
    c(-0.75, 0.25, 0.25, 0.25, 0.25, -0.75,0.25,0.25,0.25,0.25,-0.75,0.25,0.25,0.25,0.25,-0.75), 
    nrow=4, 
    ncol=4)

  
  Q1<- expm(Q*(mu*t)) #for branch length 1
  Q2<- expm(Q*(mu*t)) #for branch length 2 (same in this case)
  
  #Create L vector with number of rows equal to number of polymorphic sites, to store likelihoods
  n <- length(S[[1]])
  
  Ct1 <- matrix( nrow = 4, ncol = n)
  
  # Loop to Calculate the Conditional Likelihoods for each nucletoide, for each site
  # i loops the different sites (filling across the cols)
  # 4 rows, ACGT
  
  
  for(i in 1:n){
    tmp <- matrix( c(0,0,0,0), nrow=4, ncol=1)     
    tmp[S[[1]][i],1]<-1                             
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
    tmp[S[[2]][i],1]<-1                            
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


Likelihood3(c(1,1.9,0.1))
Likelihood2(2)

#is a minimising optimiser unless you add control$maximize=TRUE
optim(par = 1, fn=Likelihood2, lower=0.01,method="L-BFGS-B")
# Unfortunately Minimises to 0




#install.packages("optimr")
library(optimr)




n <- 1000
Data <- matrix( nrow = n, ncol = 2)
Data[,1] <- c(1:1000)*0.01
for(i in 1:n){
  Data[i,2] <- Likelihood2(Data[i,1])
}
Data <- as.data.frame(Data)
colnames(Data)<- c("t","Likelihood")

ggplot(data = Data, aes(x = t, y = Likelihood)) + geom_point(shape=1)
# geom_hline(yintercept = 0, colour = "gray65") +
# geom_vline(xintercept = 0, colour = "gray65") +
ggtitle("t against L")

#install.packages("R.matlab")
#library(R.matlab)
#Matlab$startServer()
#matlab <- Matlab()
#isOpen <- open(matlab)





