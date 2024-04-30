source("Scripts/morphotrajectory.R")

############################################################################

# load the basic data including the tree
load("Data/Tiliquini_Data.RData")

# load the scalar trees
load("Data/BayesTraits_TREES.RData")

############################################################################

body.dbn <- distance.btwn.nodes(phy = til.tree,
                                VRphy = module.trees$Body,
                                trait = all.modules$body,
                                tip = "Cyclodomorphus_michaeli",
                                sim.num = 500,
                                stat = "quantile",
                                legend = "right")

head.dbn <- distance.btwn.nodes(phy = til.tree,
                                VRphy = module.trees$Head,
                                trait = all.modules$head,
                                tip = "Tiliqua_rugosa",
                                sim.num = 500,
                                stat = "quantile",
                                legend = "right")

tail.dbn <- distance.btwn.nodes(phy = til.tree,
                                VRphy = module.trees$Tail,
                                trait = all.modules$tail,
                                tip = "Egernia_s.zellingi",
                                sim.num = 500,
                                stat = "quantile",
                                legend = "right")

limb.dbn <- distance.btwn.nodes(phy = til.tree,
                                VRphy = module.trees$Limb,
                                trait = all.modules$limb,
                                tip = "Tiliqua_rugosa",
                                sim.num = 500,
                                stat = "quantile",
                                legend = "right")

############################################################################

# plot as a grid
(limb.dbn + tail.dbn) /
  (body.dbn + head.dbn)

# Plot it for adding to multi figure panel
  ((head.dbn + theme(legend.position="none")) + plot_spacer() + plot_spacer())/
  ((limb.dbn + theme(legend.position="none")) + plot_spacer() + plot_spacer()) /
  ((tail.dbn + theme(legend.position="none")) + plot_spacer() + plot_spacer()) /
  ((body.dbn + theme(legend.position="none")) + plot_spacer() + plot_spacer())

############################################################################

