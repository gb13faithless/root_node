#### This code is also found in maxlikelihoodtree.R 


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
