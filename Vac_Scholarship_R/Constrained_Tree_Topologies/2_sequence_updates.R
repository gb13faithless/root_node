
####################################
# 2 sequences
####################################

input.sequence.file <- "data.fa" 
input.tree.file <- "data.trees"

# read in the fasta sequences
# Note, prior to doing this, need to convert fsc arlequin output to fasta format using PGPSpider
S <- read.phyDat("test.fa",format="fasta", type="DNA")
#read in the fast sim coal generated tree
T <-read.nexus("1PopDNA_1_true_trees.trees")
#If there are multiple trees, pull out an example tree
T <- T$NumGen_tree_1_1_pos_0

#examine the tree
plot(T)
#output a list of the edges
T$edge


#rename the sequences to match the tree node labels
# Note this is not really necessary

n<- length(S)
names<- c(rep(1:n))
for( i in 1:n) { 
  names[i] <- paste( i, 1, sep = ".")
}

names(S) <- names


# Setting initial parameters

#mu is mutation rate
# t is the edge length (for a 2 sequence tree, will be the same)
Likelihood2 <- function(t){
#mu <- 1
#t<- 0.4 
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
return(LL)
}

#install.packages("optimr")
library(optimr)

#is a minimising optimiser unless you add control$maximize=TRUE
optim(par = 1, fn=Likelihood2, lower=0,method="L-BFGS-B")
# Unfortunately Minimises to 0


n <- 1000
Data <- matrix( nrow = n, ncol = 2)
Data[,1] <- c(1:1000)*0.01
for(i in 1:n){
Data[i,2] <- Likelihood2(Data[i,1])
}
Data <- as.data.frame(Data)
colnames(Data)<- c("t","Likelihood")

ggplot(data = Data, aes(x = t, y = Likelihood)) + geom_point(shape=1)  +
  ggtitle("t against L")

 # geom_hline(yintercept = 0, colour = "gray65") +
 # geom_vline(xintercept = 0, colour = "gray65")



