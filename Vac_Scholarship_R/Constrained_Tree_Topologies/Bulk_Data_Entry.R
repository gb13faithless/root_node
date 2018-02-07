S1 <- read.phyDat("1.fa",format="fasta", type="DNA")
S2 <- read.phyDat("2.fa",format="fasta", type="DNA")
S3 <- read.phyDat("3.fa",format="fasta", type="DNA")
S4 <- read.phyDat("4.fa",format="fasta", type="DNA")
S5 <- read.phyDat("5.fa",format="fasta", type="DNA")
S6 <- read.phyDat("6.fa",format="fasta", type="DNA")
S7 <- read.phyDat("7.fa",format="fasta", type="DNA")
S8 <- read.phyDat("8.fa",format="fasta", type="DNA")
S9 <- read.phyDat("9.fa",format="fasta", type="DNA")
S10 <- read.phyDat("10.fa",format="fasta", type="DNA")
S11 <- read.phyDat("11.fa",format="fasta", type="DNA")
S12 <- read.phyDat("12.fa",format="fasta", type="DNA")
S13 <- read.phyDat("13.fa",format="fasta", type="DNA")
S14 <- read.phyDat("14.fa",format="fasta", type="DNA")
S15 <- read.phyDat("15.fa",format="fasta", type="DNA")
S16 <- read.phyDat("16.fa",format="fasta", type="DNA")
S17 <- read.phyDat("17.fa",format="fasta", type="DNA")
S18 <- read.phyDat("18.fa",format="fasta", type="DNA")
S19 <- read.phyDat("19.fa",format="fasta", type="DNA")
S20 <- read.phyDat("20.fa",format="fasta", type="DNA")




T3 <-read.nexus("T.trees")
T1 <- T3[1]

n <- 20 

#Create names

Names<- c(rep(1:n))
for( i in 1:n) { 
  Names[i] <- paste( "S", i, sep = "")
}
Names



#List of the Sequences

Sequence <- list(S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S16,S17,S18,S19,S20 )



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
np<- length(Sequence[[i]][[1]]) 
# overall sites
no <- 1000 

Sd <- as.DNAbin(Sequence[[i]])

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
Temp <- ML3(t3,Sequence[[i]],T1,0.00000001)
df$t1.given.f3[i] <- Temp$branch.length[1]
Temp <- ML2(t,Sequence[[i]],T1,0.00000001)
df$t1.given.f2[i] <- Temp$branch.length[1]

}
df$delta.t1 <- df$t1.given.f2-df$t1.given.f3
#df
df <- as.data.frame(df)


#### Plots

require(ggplot2)

plot <- ggplot(df, aes(delta.dist, delta.t1))+geom_point(size=0.5)
plot2 <- ggplot(df, aes(t1.given.f2, t1.given.f3))+geom_point(size=0.5) + 
  scale_x_continuous(limits = c(0, 2)) + scale_y_continuous(limits = c(0, 1))
plot2
plot
  
  
  
  
#plot <- plot + annotate(geom="text",x=0.12, y=0.195, label="Optimal Risky Portfolio (0.17988,0.16304)", color='black')
#plot <- plot + scale_x_continuous(limits = c(0, 0.3))+ scale_y_continuous(limits = c(0, 0.3)) +  
 # annotate("point", x =Standard_Deviation[Row], y =Expected_Return[Row] , colour = "red",size=2)
#plot <- plot+ 
 # ggtitle("Unrestricted Efficient Frontier with CML")
#plot <- plot+ geom_vline(xintercept = 0.1,linetype="dotted")








