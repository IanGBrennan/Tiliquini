library(phytools)

############################################################################

# load the basic data including the tree
load("Data/Tiliquini_Data.RData")

# load the processed BayesTraits MODULE data
load("Data/BayesTraits_processed_MODULES.RData")

# load the processed BayesTraits TRAIT data
load("Data/BayesTraits_processed_TRAITS.RData")

# load the focal trees
load("Data/BayesTraits_TREES.RData")

# check to make sure the traits and trees are in the same order
match(colnames(all.LSR.traits), names(trait.BT$all.trees))

############################################################################


# Estimate ancestral states for all traits
anc.traits <- NULL; anc.all.traits <- NULL
for (k in 1:length(trait.trees)){
  base.trait <- all.LSR.traits[,k]; 
  names(base.trait) <- rownames(all.LSR.traits)
  anc.traits[[k]] <- fastAnc(trait.trees[[k]], base.trait)
  anc.all.traits[[k]] <- c(base.trait, anc.traits[[k]])
  #phenogram(til.tree, x=c(base.trait, anc.traits[[k]]), fsize=0.5)
}
names(anc.traits) <- names(trait.trees)
names(anc.all.traits) <- names(trait.trees)

############################################################################


# Extract trait values at 0.1 million year intervals
source("Scripts/trait.at.time.R")
all.anc <- lapply(1:length(anc.all.traits), 
                  function(x) trait.at.time(timeslices=0.1, 
                                            phy=til.tree, 
                                            trait.vector=anc.all.traits[[x]], 
                                            plot=F))
names(all.anc) <- names(anc.all.traits)

############################################################################


# Save these data
save(anc.traits, anc.all.traits, all.anc, file = "Data/Tiliquini_Ancestral_TRAITS.Rdata")



