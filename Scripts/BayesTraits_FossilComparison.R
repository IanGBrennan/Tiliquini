source("Scripts/plotting_BayesTraits.R")

setwd("/Applications/BayesTraitsV4/PPPostProcess")



fossil.tree <- read.nexus("/Applications/BayesTraitsV4/Tiliquini_BT_Fossil_HL.tre")
extant.tree <- read.nexus("/Applications/BayesTraitsV4/Tiliquini_BT.tre")

extant.BT <- process_PPP(res.files = "HeadLength_noFossil_TraitResults.txt", 
                         phy = extant.tree, col.palette="YlGnBu")

fossil.BT <- process_PPP(res.files = "HeadLength_Fossil_TraitResults.txt", 
                         phy = fossil.tree, col.palette="YlGnBu")

par(mfrow=c(2,1))

plot.VarRates.tree(BT = extant.BT$all.res$Head, # res object for focal trait
                   phy = extant.tree, # tree for focal trait
                   col.palette = "YlGnBu", # choose your color palette
                   legend = T, # want a rate legend included?
                   tree.type = "phylogram", # tree shape
                   log.rates = F, # should rates be logged (probably)
                   relative.rates = F, # should rates be relative (more dramatic! log.rates must be = F)
                   trait = "Head Length", # name of the trait
                   outline = F, # black outline on branches for light colors
                   pos.selection = T, # indicate positive selection (r >= 10)
                   shift.rates = T, # indicate shifted rates (r >= 2)
                   annotations = T) # include info on # of shifts?


plot.VarRates.tree(BT = fossil.BT$all.res$Head, # res object for focal trait
                   phy = fossil.tree, # tree for focal trait
                   col.palette = "YlGnBu", # choose your color palette
                   legend = T, # want a rate legend included?
                   tree.type = "phylogram", # tree shape
                   log.rates = F, # should rates be logged (probably)
                   relative.rates = F, # should rates be relative (more dramatic! log.rates must be = F)
                   trait = "Head Length", # name of the trait
                   outline = F, # black outline on branches for light colors
                   pos.selection = T, # indicate positive selection (r >= 10)
                   shift.rates = T, # indicate shifted rates (r >= 2)
                   annotations = T) # include info on # of shifts?
