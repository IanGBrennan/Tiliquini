source("Scripts/rate.trajectory.R")

############################################################################

# load the basic data including the tree
load("Data/Tiliquini_Data.RData")

# load the ancestral trait data
load("Data/Tiliquini_Ancestral_TRAITS.Rdata")

# load the processed BayesTraits TRAIT data
load("Data/BayesTraits_processed_MODULES.RData")

############################################################################

# Rate Trajectory Plots

# Rate Trajectory HEAD
HEAD.anc.module <- data.frame(HeadWidth = anc.all.traits$HeadWidth, 
                              SnoutEye = anc.all.traits$SnoutEye, 
                              HeadDepth = anc.all.traits$HeadDepth, 
                              PosSkull = anc.all.traits$PosSkull)
HEAD.plot <- rate.trajectory.BT.module(tree = til.tree, module = HEAD.anc.module,
                                       PPP.all.res = module.BT$all.res$Head,
                                       tip.spread = c("Tiliqua_rugosa","Tiliqua_scincoides"),
                                       focus = "clade", psize=3, lsize=2, background.color="grey",
                                       gimme.the.data=F, inset=T, relative.rates=T,
                                       col.palette="YlGnBu")
HEAD.rate <- rate.trajectory.BT.module(tree = til.tree, module = HEAD.anc.module,
                                       PPP.all.res = module.BT$all.res$Head,
                                       tip.spread = c("Cyclodomorphus_michaeli","Cyclodomorphus_praealtus"),
                                       focus = "clade", psize=3, lsize=2, background.color="grey",
                                       gimme.the.data=T, inset=T)


# Rate Trajectory BODY
BODY.anc.module <- data.frame(Interlimb = anc.all.traits$Interlimb, 
                              BodyWidth = anc.all.traits$BodyWidth, 
                              PelvicWidth = anc.all.traits$PelvicWidth, 
                              PelvicHeight = anc.all.traits$PelvicHeight)
BODY.plot <- rate.trajectory.BT.module(tree = til.tree, module = BODY.anc.module,
                                       PPP.all.res = module.BT$all.res$Body,
                                       tip.spread = c("Tiliqua_rugosa","Cyclodomorphus_maximus"),
                                       focus = "clade", psize=3, lsize=2, background.color="grey",
                                       gimme.the.data=F, inset=T, relative.rates=T,
                                       col.palette="YlGnBu")
BODY.rate <- rate.trajectory.BT.module(tree = til.tree, module = BODY.anc.module,
                                       PPP.all.res = module.BT$all.res$Body,
                                       tip.spread = c("Cyclodomorphus_michaeli","Cyclodomorphus_praealtus"),
                                       focus = "clade", psize=3, lsize=2, background.color="grey",
                                       gimme.the.data=T, inset=T)


# Rate Trajectory and Rate-to-Node for TAIL
TAIL.anc.module <- data.frame(TailLength = anc.all.traits$TailLength, 
                              TailWidth = anc.all.traits$TailWidth)
TAIL.plot <- rate.trajectory.BT.module(tree = til.tree, module = TAIL.anc.module,
                                       PPP.all.res = module.BT$all.res$Tail,
                                       tip.spread = c("Egernia_depressa","Egernia_cygnitos"),
                                       focus = "clade", psize=3, lsize=2, background.color="grey",
                                       gimme.the.data=F, inset=T, relative.rates=T,
                                       col.palette="YlGnBu")
TAIL.rate <- rate.trajectory.BT.module(tree = til.tree, module = TAIL.anc.module,
                                       PPP.all.res = module.BT$all.res$Tail,
                                       tip.spread = c("Cyclodomorphus_michaeli","Cyclodomorphus_praealtus"),
                                       focus = "clade", psize=3, lsize=2, background.color="grey",
                                       gimme.the.data=T, inset=T)


# Rate Trajectory and Rate-to-Node for LIMB
LIMB.anc.module <- data.frame(UpperArm = anc.all.traits$UpperArm, 
                              LowerArm = anc.all.traits$LowerArm,
                              Hand = anc.all.traits$Hand,
                              UpperLeg = anc.all.traits$UpperLeg,
                              LowerLeg= anc.all.traits$LowerLeg,
                              Foot = anc.all.traits$Foot)
LIMB.plot <- rate.trajectory.BT.module(tree = til.tree, module = LIMB.anc.module,
                                       PPP.all.res = module.BT$all.res$Limb,
                                       tip.spread = c("Tiliqua_rugosa","Cyclodomorphus_maximus"),
                                       focus = "clade", psize=3, lsize=2, background.color="grey",
                                       gimme.the.data=F, inset=T, relative.rates=T,
                                       col.palette="YlGnBu")
LIMB.rate <- rate.trajectory.BT.module(tree = til.tree, module = LIMB.anc.module,
                                       PPP.all.res = module.BT$all.res$Limb,
                                       tip.spread = c("Cyclodomorphus_michaeli","Cyclodomorphus_praealtus"),
                                       focus = "clade", psize=3, lsize=2, background.color="grey",
                                       gimme.the.data=T, inset=T)

############################################################################


# combine the module rate data frames
rate.list <- list(head=HEAD.rate, body=BODY.rate, tail=TAIL.rate, limb=LIMB.rate)
# and write to file
save(rate.list, file="Data/rate.trajectories.RData")

############################################################################

# Wrap plot
(HEAD.plot + LIMB.plot) /
  (TAIL.plot + BODY.plot)


# Plot it for adding to multi figure panel
(plot_spacer()   + plot_spacer() + HEAD.plot)/
  (plot_spacer() + plot_spacer() + LIMB.plot) /
  (plot_spacer() + plot_spacer() + TAIL.plot) /
  (plot_spacer() + plot_spacer() + BODY.plot)




