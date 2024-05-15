source("Scripts/rate.trajectory.R")

############################################################################

# load the basic data including the tree
load("Data/Tiliquini_Data.RData")

# load the ancestral trait data
load("Data/Tiliquini_Ancestral_TRAITS.Rdata")

# load the processed BayesTraits TRAIT data
load("Data/BayesTraits_processed_MODULES.RData")

############################################################################




HEAD.rtn <- rate.to.node.BT(phy=til.tree, PPP.obj=module.BT$all.res$Head, psize=3, lsize=2, 
                            log.rates=F, col.palette="YlGnBu", relative.rates=T)

TAIL.rtn <- rate.to.node.BT(phy=til.tree, PPP.obj=module.BT$all.res$Tail, psize=3, lsize=2, 
                            log.rates=F, col.palette="YlGnBu", relative.rates=T)

LIMB.rtn <- rate.to.node.BT(phy=til.tree, PPP.obj=module.BT$all.res$Limb, psize=3, lsize=2, 
                            log.rates=F, col.palette="YlGnBu", relative.rates=T)

BODY.rtn <- rate.to.node.BT(phy=til.tree, PPP.obj=module.BT$all.res$Body, psize=3, lsize=2, 
                            log.rates=F, col.palette="Yl", relative.rates=T)


############################################################################

# Wrap plot
(HEAD.rtn + LIMB.rtn) /
  (TAIL.rtn + BODY.rtn)


# Plot it for adding to multi figure panel
(plot_spacer()   + HEAD.rtn + plot_spacer())/
  (plot_spacer() + LIMB.rtn + plot_spacer()) /
  (plot_spacer() + TAIL.rtn + plot_spacer()) /
  (plot_spacer() + BODY.rtn + plot_spacer())






