# This file will attempt to put together 2 and 3 sequence Likelihoods 

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

# read in the 3 fasta sequence XXX
XXX <- read.phyDat("XXX.fa",format="fasta", type="DNA")
# read in the 3 fasta sequence XXY
XXY <- read.phyDat("XXY.fa",format="fasta", type="DNA")
# read in the 3 fasta sequence XYX
XYX <- read.phyDat("XYX.fa",format="fasta", type="DNA")
# read in the 3 fasta sequence XYX
XYZ <- read.phyDat("XYZ.fa",format="fasta", type="DNA")


## Extract the XX sequence
XX <- XXX[c(1,2)]
## Extract the XY sequence
XY <- XYX[c(1,2)]


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


#seed tvalues
t3 <- c(1,1,1)
t <- 1

# 3 sequence Data

ML3(t3,XYZ,T3,0.001)
ML2(t,XY,T3,0.001)

################################################################
# Complicated Sequence
################################################################

Complex3 <- read.phyDat("Complex.fa",format="fasta", type="DNA")








#Change names
n<- length(Complex3)
names<- c(rep(1:n))
for( i in 1:n) { 
  names[i] <- paste( i, 1, sep = ".")
}
names(Complex3) <- names

## Extract the XX sequence
Complex2 <- Complex3[c(1,2)]


ML3(t3,Complex3,T3,0.001)
ML2(t,Complex2,T2,0.001)

################################################################
#matrix of distances using dist.dna
################################################################

# the issue with using dist.dna is that it only counts the polymorphic sites in determining the distance
# this is because FSC outputs only the polymorphic sites
# however, you specify the number of demes in FSC - that should be the number of overall sites
# therefore, you need to multiply the dna Dist by (# of polymorphic sites/# overall sites)


# polymorphic sites
np<- length(Complex3[[1]]) 
# overall sites
no <- 1000 
Sd <- as.DNAbin(Complex3)


DMraw<- dist.dna(Sd,model="raw", as.matrix = TRUE)
DMraw <- DMraw * (np/no)
DMraw


c#extract submatrix of distances between ingroup and outgroup
sub<- DMraw[3,1:2]

#function to get the minimum non negative number
minpositive = function(x) min(x[x > 0])

av.distance<-mean(sub) # av distance near 0.75 suggests near random evolution
max.distance <-max(sub)
min.distance <-minpositive(sub)

#extract submatrix of distances between ingroup 
sub.in<- DMraw[1:2,1:2]

av.distance.in<-mean(sub.in)
max.distance.in <-max(sub.in)
min.distance.in <-minpositive(sub.in)

#extract submatrix of distances between outgroup 
sub.out<- DMraw[3,3]

av.distance.out<-mean(sub.out)
max.distance.out <-max(sub.out)
min.distance.out <-minpositive(sub.out)















################################################################
# Plotting Graphs
################################################################

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
