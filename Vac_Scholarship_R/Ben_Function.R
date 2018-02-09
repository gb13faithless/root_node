

species.fas <- read.dna(file='test.fa',format='fas')
species.tree <- makeTree(species.fas)
plot(species.tree)

# Change the tipnames to numbers

# From the sequences, calculate genetic distances based
# on the F81 model. F81 is the MOST complex model possible 
# in the function.


makeTree <- function(obj.phy,excluded='none',method='NJ'){
  # A function for making a ML-tree from a phydat object
  # Created by: AB Rohrlach
  # Requires a sequence object, how to deal with missing data, and a method choice.
  
  # Check to see that the object actually is a DNAbin,
  # if so, convert to a phyDat object
  if(class(obj.phy)=="DNAbin"){
    obj.phy <- as.phyDat(obj.phy)
  }
  
  # Grab the sequence names to save
  obj.names <- names(obj.phy)
  
  # Change the tipnames to numbers
  names(obj.phy) <- 1:length(names(obj.phy))
  
  # From the sequences, calculate genetic distances based
  # on the F81 model. F81 is the MOST complex model possible 
  # in the function.
  dm = dist.ml(obj.phy,exclude=excluded,model = 'F81')
  
  # Can either use Neighbour joining (NJ) or the
  # Unweighted Pair Group Method with Arithmetic Mean (upgma)
  # method (default is NJ)
  
  # Make the tree (built in function)
  if(method=='NJ'){
    tree = NJ(dm)
  }
  if(method=='upgma'){
    tree = upgma(dm)
  }
  
  # Perform some standard tweaking methods to explore
  # treespace a little
  treeRA <- random.addition(obj.phy,method='fitch')
  treeNNI <- optim.parsimony(tree, obj.phy)
  treeRatchet <- pratchet(obj.phy, start=tree)
  treeRatchet$tip.label <- obj.names
  
  # Return the true tip names and convert to tree format
  names(obj.phy) <- obj.names
  return(acctran(treeRatchet, obj.phy))
}
