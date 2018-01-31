#packages, idk how many are actually needed
install.packages("ape")
install.packages("phangorn")
install.packages("seqinr")
install.packages("tree")
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


# use pgp spider to convert the data to fasta from arp in the fsc output

# read in the fasta sequences
S <- read.phyDat("test.fa",format="fasta", type="DNA")
S2 <- read.phyDat("test2.fa",format="fasta", type="DNA")
#read in the fsc tree
T <-read.nexus("35_1_true_trees.trees")
#pull out an example tree
T <- T$NumGen_tree_1_1_pos_0
plot(T)


#specify ingroup and outgroup sizes
n1 <- 30
n2 <- 5


# some useful code cleans the names of the sequences *like a find and replace function
#names(S) <- gsub(" ", "", names(S))
#names(S) <- gsub("_", ".", names(S))

# loops to match the tip labels of the fasta seqeunces with the tip labels of the tree, to use pml

n<- length(S)
names<- c(rep(1:35))
for( i in 1:30) { 
names[i] <- paste( i, 1, sep = ".")
} 
for( i in 31:35) { 
  names[i] <- paste( i, 2, sep = ".")
} 

names(S) <- names
names(S2) <- names

# calculate the likelihood
# defaults to Jukes Cantor
# Can fit other models

fit = pml(T, data=S, rearrangment="none",rate = 0.00002)
fit
fit2 = pml(T, data=S2, rearrangment="none",rate = 0.00002)
fit2



fit$logLik  #tree log likelihood
fit$data # alignments
fit$tree$edge.length   # branch lengths    
fit$siteLik    # site log likelihoods
plot(fit)

fit2$logLik  #tree log likelihood
fit2$data # alignments
fit2$tree$edge.length   # branch lengths    
fit2$siteLik    # site log likelihoods
plot(fit2)




# optimises the branch lengths

fitJC  <- optim.pml(fit, TRUE)
logLik(fitJC)
fitJC$logLik  #tree log likelihood
fitJC$data # alignments
fitJC$tree$edge.length   # branch lengths    
fitJC$siteLik    # site log likelihoods
plot(fit)
plot(fitJC)



#function to calculate the hamming distance between 2 strings

hamming <- function(X, Y) {
  if ( missing(Y) ) {
    uniqs <- unique(as.vector(X))
    U <- X == uniqs[1]
    H <- t(U) %*% U
    for ( uniq in uniqs[-1] ) {
      U <- X == uniq
      H <- H + t(U) %*% U
    }
  } else {
    uniqs <- union(X, Y)
    H <- t(X == uniqs[1]) %*% (Y == uniqs[1])
    for ( uniq in uniqs[-1] ) {
      H <- H + t(X == uniq) %*% (Y == uniq)
    }
  }
  nrow(X) - H
}


# Convert the character lists into a vectorised format that the hamming format can use

install.packages("stringdist")
install.packages("e1071")
library(stringdist)
library(e1071)

x<- as.character(S$`1.1`)
x<- strsplit(x, "")
l1 <- length(x)
s1 <- array(dim=c(l1,1))
for(i in 1:l1){
  s1[[i,1]]<- x[[i]][1]
}
s1

y<- as.character(S$`2.1`)
y<- strsplit(y, "")
l2 <- length(y)
s2 <- array(dim=c(l2,1))
for(i in 1:l2){
  s2[[i,1]]<- y[[i]][1]
}
s2

#calculate the hamming distance
hamming(s1,s2)




#extract ingroup only

subtrees(T)
# need to find a way to characterise the subtrees better
T.in <- subtrees(T)[[6]]


fit.in = pml(T.in, data=S, rearrangment="none",rate = 0.00002)
fit.in

fit.in$logLik  #tree log likelihood
fit.in$data # alignments
fit.in$tree$edge.length   # branch lengths    
fit.in$siteLik    # site log likelihoods
plot(fit.in)

fitJC.in  <- optim.pml(fit.in, TRUE, rearrangement = "none")
logLik(fitJC.in)
fitJC.in$logLik  #tree log likelihood
fitJC.in$data # alignments
fitJC.in$tree$edge.length   # branch lengths    
fitJC.in$siteLik    # site log likelihoods
plot(fit.in)
plot(fitJC.in)

#matrix of distances using dist.dna

Sd <- as.DNAbin(S)
DM<- dist.dna(Sd,model="JC69", as.matrix = TRUE)

DM

DMraw<- dist.dna(Sd,model="raw", as.matrix = TRUE)
DMraw

#extract submatrix of distances between ingroup and outgroup
sub<- DMraw[31:35,1:30]

#function to get the minimum non negative number
minpositive = function(x) min(x[x > 0])

av.distance<-mean(sub) # av distance near 0.75 suggests near random evolution
max.distance <-max(sub)
min.distance <-minpositive(sub)

#extract submatrix of distances between ingroup 
sub.in<- DMraw[1:30,1:30]

av.distance.in<-mean(sub.in)
max.distance.in <-max(sub.in)
min.distance.in <-minpositive(sub.in)

#extract submatrix of distances between outgroup 
sub.out<- DMraw[31:35,31:35]

av.distance.out<-mean(sub.out)
max.distance.out <-max(sub.out)
min.distance.out <-minpositive(sub.out)



## Note If the sequences are very different, most evolutionary distances are undefined
## and a non-finite value (Inf or NaN) is returned. You may do dist.dna(, model = "raw") 
## to check whether some values are higher than 0.75.


#ggtree

Tr <- ggtree(fit$tree)
Tr
# Add Scale
Tr <- Tr + geom_treescale()

# Internal Node Numbers
Tr +geom_text2(aes(subset=!isTip, label=node), hjust=-.3) 
# Tip Labels
Tr <- Tr + geom_tiplab()
# Tip Points
Tr + geom_tippoint()
# Tip Points
Tr +geom_text2(aes(subset=!isTip, label=node), hjust=1) 






D1  <- read.phyDat( "test.fa",format="fasta", type="DNA")
D1
#dm1  <- dist.ml(D1, "F81")
#dm1
D2 <- read.phyDat( "test2.fa", format = "fasta",type="DNA")
#dm2 <- dist.ml(D2, "F81")

#MyTree$tip.label <- fit$tree$tip.label

# create usual tree
#treeUPGMA  <- upgma(dm) #neighbour joining to get the topology
#treeUPGMA$edge.length 


# create circle tree
#treeNJ  <- NJ(dm) #neighbour joining to get the topology

#read in tree
T <- read.nexus("RootNodeExample_1_true_trees.trees")
T1 <- T$NumGen_tree_1_1_pos_0
T2 <- T$NumGen_tree_10_1_pos_0

T3 <-subtrees(T1)

T<- T3[[2]]

# max likelihood for supplied tree

help(pml)
help(optim.pml)

fit3 = pml(T3[[2]], data=D1,rate=0.0002, rearrangment="none")
fit3
fit3$logLik  #tree log likelihood
fit3$data # alignments
fit3$tree$edge.length   # branch lengths    
fit3$siteLik    # site log likelihoods
plot(fit3)



fit4  <- optim.pml(fit3, TRUE)
fit4
plot(fit4)





names(D1)<-paste0(1:50, ".1")


names(D2)<-c("1.1","2.1","3.1","4.1","5.1","6.1","7.1","8.1","9.1","10.1")

#T1$tip.label <- names(D1)
#T2$tip.label <- names(D2)

names(D1)

#M <- MyTree$NumGen_tree_2_1_pos_0
getAnywhere('pml')
getAnywhere('optim.pml')
getAnywhere('pml.fit')

tree <- reorder(tree, "postorder")
tree




#
layout(matrix(c(1,2), 2, 1), height=c(1,2))
par(mar = c(0,0,2,0)+ 0.1)
plot(T1, main="Tree")
plot.tree(MyTree, "unrooted", main="NJ")

# max likelihood for supplied tree
fit = pml(T1, data=D1,rearrangment="none")
fit
fit$logLik  #tree log likelihood
fit$data # alignments
fit$tree$edge.length   # branch lengths    
fit$siteLik    # site log likelihoods
plot(fit)

pml.control(epsilon = 1e-08, maxit = 10, trace = 1)
fit2 = pml(T2, data=D2,rate=0.0002)
fit2
fit2$logLik  #tree log likelihood
fit2$data # alignments
fit2$tree$edge.length   # branch lengths   

fit2$siteLik    # site log likelihoods
plot(fit2)






# optimises branch length
fitJC  <- optim.pml(fit, TRUE, rearrangement = "none")
logLik(fitJC)
fitJC$logLik  #tree log likelihood
fitJC$data # alignments
fitJC$tree$edge.length   # branch lengths    
fitJC$siteLik    # site log likelihoods
plot(fit)
plot(fitJC)


fitJC2  <- optim.pml(fit2, TRUE, rearrangement = "stochastic")
logLik(fitJC2)
fitJC2$logLik  #tree log likelihood
fitJC2$data # alignments
fitJC2$tree$edge.length   # branch lengths    
fitJC2$siteLik    # site log likelihoods
plot(fitJC2)
plot(fit2)
D2



fit2$tree$edge.length

dev.off()

plot(MyTree)










help(pml)

# lists methods for pml with Jukes Cantor model
methods(class="pml")

help(optim.pml)


# Example of GTR model

fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE, 
                    rearrangement = "NNI", control=pml.control(trace=0))
fitGTR


fit2 <- pml(treeNJ, T_phys)
print(fit2)

#fitJC <- optim.pml(fit2, model = "JC", rearrangement = "stochastic")
#logLik(fitJC)
#bs <- bootstrap.pml(fitJC, bs=100, optNni=TRUE, multicore=TRUE, control = pml.control(trace=0))
#plotBS(midpoint(fitJC$tree), bs, p = 50, type="p")


