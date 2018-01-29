#install.packages("readr") 
#install.packages("stringr")
install.packages("seqinr")
library(readr)
library(stringr)
library(seqinr)
#read in data using read.fasta
mystring <- read.fasta("rosalind_gc.txt",as.string=TRUE, forceDNAtolower = FALSE)
Count <- length(mystring)

#extract counts of A C G T
answers <- array(1:Count*4, dim=c(Count,4))
answers[1,]
for (i in 1:Count){
answers[i,] <- c(str_count(mystring[i], "A"),str_count(mystring[i], "C"),str_count(mystring[i], "G"),str_count(mystring[i], "T"))
}
# CG
answers <- cbind(answers, answers[,2]+answers[,3], answers[,1]+answers[,2]+answers[,3]+answers[,4]) 
answers
#Total
answers <- cbind(answers, 100*answers[,5]/answers[,6])
answers[,7]
#create data frame
df <- data.frame(names(mystring),answers)
df

answer <- c(1,2)
answer[2] <- max(answers[,7])
answer[1]<- names(mystring)[which.max(answers[,7])]
answer
