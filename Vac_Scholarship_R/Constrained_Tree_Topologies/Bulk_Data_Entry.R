
# These Packages are Required for Manipulating the Tree
library(tree)
library(ape)
library(phangorn)
library(seqinr)
library(Biostrings)
library(ggplot2)
library(ggtree)

#install.packages("optimr") to do the maximising
library(optimr)


# reads in all the fasta files


T <-read.nexus("3seq2_1_true_trees.trees")
plot(T$NumGen_tree_1_1_pos_0)
plot(T$NumGen_tree_100_1_pos_0)

  
ldf3 <- list() # creates a list
listfa <- dir(pattern = "*.fa") # creates the list of all the fa files in the directory
for (k in 1:length(listfa)){
  ldf3[[k]] <- read.phyDat(listfa[k],format="fasta", type="DNA")
}




n <- length(listfa)

Names<- c(rep(1:n))
for( i in 1:n) { 
  Names[i] <- paste( "sim", i, sep = "")
}
Names


#Create names


# Create Vectors for t1
t1.given.f3 <- c(1:n)
t1.given.f2 <- c(1:n)
delta.t1 <- c(1:n)
max.dist <- c(1:n)
min.dist <- c(1:n)
delta.dist <- c(1:n)

df <- DataFrame(Names,t1.given.f2,t1.given.f3,delta.t1,max.dist,min.dist,delta.dist)

##############################################################
#matrix of distances using dist.dna
################################################################

# the issue with using dist.dna is that it only counts the polymorphic sites in determining the distance
# this is because FSC outputs only the polymorphic sites
# however, you specify the number of demes in FSC - that should be the number of overall sites
# therefore, you need to multiply the dna Dist by (# of polymorphic sites/# overall sites)

#function to get the minimum non negative number
minpositive = function(x) min(x[x > 0])



for(i in 1:n){
# polymorphic sites
np<- length(ldf3[[i]][[1]]) 
# overall sites
no <- 10000

Sd <- as.DNAbin(ldf3[[i]])

DMraw<- dist.dna(Sd,model="raw", as.matrix = TRUE)
DMraw <- DMraw * (np/no)
DMraw

df$max.dist[i] <-max(DMraw)
df$min.dist[i] <-minpositive(DMraw)
}

df$delta.dist <- df$max.dist-df$min.dist




t3 <- c(1,1,1)
t <- 1



for(i in 1:n){
Temp <- ML3(t3,ldf3[[i]],T1,0.00001)
df$t1.given.f3[i] <- Temp$branch.length[1]
Temp <- ML2(t,ldf3[[i]],T1,0.00001)
df$t1.given.f2[i] <- Temp$branch.length[1]

}
df$delta.t1 <- df$t1.given.f2-df$t1.given.f3
#df
df <- as.data.frame(df)


log(0.25)

#### Plots

library(ggplot2)

plot <- ggplot(df, aes(delta.dist, delta.t1))+geom_point(size=0.5)
plot

plot2 <- ggplot(df, aes(t1.given.f2, t1.given.f3))+geom_point(size=0.5) + 
  scale_x_continuous(limits = c(0, 2)) + scale_y_continuous(limits = c(0, 1))
plot2 <- plot2 + geom_abline(a=1,b=0)
plot2 <- plot2 + stat_smooth(method="lm", se=TRUE,
                             formula= myformula,colour="red")
plot2
  
# poly function to fit polynomial
  
  
model <- lm(df$t1.given.f3 ~ poly(df$t1.given.f2,2))
myformula <- y ~ poly(x, 2)
summary(model)

  
#plot <- plot + annotate(geom="text",x=0.12, y=0.195, label="Optimal Risky Portfolio (0.17988,0.16304)", color='black')
#plot <- plot + scale_x_continuous(limits = c(0, 0.3))+ scale_y_continuous(limits = c(0, 0.3)) +  
 # annotate("point", x =Standard_Deviation[Row], y =Expected_Return[Row] , colour = "red",size=2)
#plot <- plot+ 
 # ggtitle("Unrestricted Efficient Frontier with CML")
#plot <- plot+ geom_vline(xintercept = 0.1,linetype="dotted")



head(df)




