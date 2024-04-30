source("Scripts/plotting_BayesTraits.R")

# plot the location of inferred rate shifts on trees of interest

# explanation of the circles:
# white is for possible shifts (found in >=50% of sampled posterior)
# grey is for likely shifts (found in >=70% of sampled posterior)
# black is for definite shifts (found in >=95% of sampled posterior)


##############################################################################
# start with the MODULES

# make a plot matrix to catch the plots
layout(matrix(c(1,2,3,4,5,6,7,8), 4, 2, byrow = F))

# load the basic data including the tree
load("Data/Tiliquini_Data.RData")

# load the processed BayesTraits data
load("Data/BayesTraits_processed_MODULES.RData")

plot.VarRates.tree(BT = module.BT$all.res$Head, # res object for focal trait
                   phy = til.tree, # tree for focal trait
                   col.palette = "GnBu", # choose your color palette
                   legend = T, # want a rate legend included?
                   tree.type = "phylogram", # tree shape
                   log.rates = F, # should rates be logged (probably)
                   relative.rates = T, # should rates be relative (more dramatic! log.rates must be = F)
                   trait = "Head", # name of the trait
                   outline = F, # black outline on branches for light colors
                   pos.selection = T, # indicate positive selection (r >= 10)
                   shift.rates = T, # indicate shifted rates (r >= 2)
                   annotations = F) # include info on # of shifts?

plot.VarRates.tree(BT = module.BT$all.res$Limb, # res object for focal trait
                   phy = til.tree, # tree for focal trait
                   col.palette = "YlGnBu", # choose your color palette
                   legend = F, # want a rate legend included?
                   tree.type = "phylogram", # tree shape
                   log.rates = F, # should rates be logged (probably)
                   relative.rates = T, # should rates be relative (more dramatic! log.rates must be = F)
                   trait = "Limb", # name of the trait
                   outline = F, # black outline on branches for light colors
                   pos.selection = T, # indicate positive selection (r >= 10)
                   shift.rates = T, # indicate shifted rates (r >= 2)
                   annotations = F) # include info on # of shifts?

plot.VarRates.tree(BT = module.BT$all.res$Tail, # res object for focal trait
                   phy = til.tree, # tree for focal trait
                   col.palette = "YlGnBu", # choose your color palette
                   legend = F, # want a rate legend included?
                   tree.type = "phylogram", # tree shape
                   log.rates = F, # should rates be logged (probably)
                   relative.rates = T, # should rates be relative (more dramatic! log.rates must be = F)
                   trait = "Tail", # name of the trait
                   outline = F, # black outline on branches for light colors
                   pos.selection = T, # indicate positive selection (r >= 10)
                   shift.rates = T, # indicate shifted rates (r >= 2)
                   annotations = F) # include info on # of shifts?

plot.VarRates.tree(BT = module.BT$all.res$Body, # res object for focal trait
                   phy = til.tree, # tree for focal trait
                   col.palette = "YlGnBu", # choose your color palette
                   legend = F, # want a rate legend included?
                   tree.type = "phylogram", # tree shape
                   log.rates = F, # should rates be logged (probably)
                   relative.rates = T, # should rates be relative (more dramatic! log.rates must be = F)
                   trait = "Body", # name of the trait
                   outline = F, # black outline on branches for light colors
                   pos.selection = T, # indicate positive selection (r >= 10)
                   shift.rates = T, # indicate shifted rates (r >= 2)
                   annotations = F) # include info on # of shifts?





##############################################################################
# now onto all the individual TRAITS

# load the processed BayesTraits data
load("Data/BayesTraits_processed_TRAITS.RData")

# set up a plot space for all trait trees
par(mfrow=c(5,4))

# select the specific trees of interest
trait.trees <- trait.BT$mean.scalar.trees

# sort the order of the trees to match the columns
all.LSR.traits <- dplyr::select(all.LSR.traits, names(trait.trees))

# plot all the trees with edges colored by rates
for(j in 1:length(trait.trees)){plot.VarRates.tree(BT=trait.BT$all.res[[j]], 
                                                   phy=trait.BT$scalar.trees[[j]],
                                                   col.palette="YlGnBu", legend=F, tree.type="phylogram", 
                                                   log.rates=F, relative.rates=T,
                                                   trait=names(trait.BT$all.res)[[j]], outline=F, 
                                                   pos.selection=T, shift.rates=T, annotations=T)}


