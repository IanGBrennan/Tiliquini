source("Scripts/morphotrajectory.R")

############################################################################

# load the basic data including the tree
load("Data/Tiliquini_Data.RData")

# load the scalar trees
load("Data/BayesTraits_TREES.RData")

############################################################################

# Generate Sim-to-Node plots

tl.tw.plot <- sim.to.node(phy=til.tree, VRphy=module.trees$Tail, 
                          trait=all.modules$tail[,c("TailLength","TailWidth")],
                          tip="Egernia_s.zellingi", sim.num=500, legend="none")

hl.fw.plot <- sim.to.node(phy=til.tree, VRphy=module.trees$Limb, 
                          trait=all.modules$limb[,c("Hand","Foot")],
                          tip="Tiliqua_rugosa", sim.num=500, legend="none")

il.bw.plot <- sim.to.node(phy=til.tree, VRphy=module.trees$Body, 
                          trait=all.modules$body[,c("Interlimb","BodyWidth")],
                          tip="Cyclodomorphus_michaeli", sim.num=500, legend="none")

hw.hd.plot <- sim.to.node(phy=til.tree, VRphy=module.trees$Head, 
                          trait=all.modules$head[,c("HeadWidth","HeadDepth")],
                          tip="Tiliqua_rugosa", sim.num=500, legend="none")

############################################################################

# plot as a grid
(hl.fw.plot + tl.tw.plot) /
  (il.bw.plot + hw.hd.plot)

# Plot it for adding to multi figure panel
  (hw.hd.plot + plot_spacer() + plot_spacer())/
  (hl.fw.plot + plot_spacer() + plot_spacer()) /
  (tl.tw.plot + plot_spacer() + plot_spacer()) /
  (il.bw.plot + plot_spacer() + plot_spacer())

############################################################################



