install.packages("ape")
install.packages("Biostrings")
devtools::install_github("GuangchuangYu/ggtree")
install.packages("tidyverse")
install.packages("devtools")
source("http://bioconductor.org/biocLite.R")
biocLite("Biostrings")
install.packages("phangorn")
library("phangorn")
library("ape")
library("Biostrings")
library("ggplot2")
library("ggtree")
library("devtools")
library("strataG")
install.packages("phylobase")
library("phylobase")


# 10 sample tree, 19 nodes
MyTree <- read.nexus("1PopDNA_1_true_trees.trees")



# Make Tree
T <- ggtree(MyTree)
T



# Add Scale
T + geom_treescale()
# Circular Layout
ggtree(MyTree, layout="circular") + ggtitle("(Phylogram) circular layout")
# Internal Node Numbers
T +geom_text2(aes(subset=!isTip, label=node), hjust=-.3) 
# Tip Labels
T + geom_tiplab()
# Tip Points
T + geom_tippoint()
# Viewing Clades
viewClade(T+geom_tiplab(), node=13)


# Group Clades
MyTree <- groupClade(MyTree, node=c(13,19))
# Tree with Grouped Clade in different colours
T2 <- ggtree(MyTree, aes(color=group, linetype=group))
T2 + geom_tiplab(aes(subset=(group==2)))
T2

# Group by Tip Label
MyTree <- groupOTU(MyTree, focus=c("2.1", "5.1", "9.1", "6.1"))
ggtree(MyTree, aes(color=group)) + geom_tiplab()



#Clade Labels
T+geom_cladelabel(node=13, label="Clade 1", color = "red")




#annotations

# number of nodes needed for number of colours

d = data.frame(node=T$data$node, color=sample(c('red', 'blue', 'green'), 19, replace=TRUE),value = 1:19, place = c(rep("GZ", 6), rep("HK", 3), rep("CZ", 10)))


# magrittr
T %<+% d + aes(color=I(color))+ geom_tiplab(aes(color=I(color))) +
  geom_tippoint(aes(size=value, shape=place, color=I(color)), alpha=0.25) +theme(legend.position="right")

# other options
T %<+% d + geom_text(aes(color=I(color), label=node), hjust=1, vjust=-0.4, size=3) +
  geom_text(aes(color=I(color), label=value), hjust=1, vjust=1.4, size=3)



#Highlight Clades
T + geom_hilight(node=13, fill="steelblue", alpha=.6)















