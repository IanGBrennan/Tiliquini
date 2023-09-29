source("/Scripts/plotting_BayesTraits.R")
source("/Scripts/trait.at.time.R")
source("/Scripts/rate.trajectory.R")

# Rate Trajectory and Rate-to-Node for HEAD
HEAD.anc.module <- data.frame(Head_Width = anc.all.traits$Head_Width, 
                              Snout_Eye = anc.all.traits$Snout_Eye, 
                              Head_Depth = anc.all.traits$Head_Depth, 
                              Pos_Skull = anc.all.traits$Pos_Skull)
HEAD.plot <- rate.trajectory.BT.module(tree = egernia.tree, module = HEAD.anc.module,
                                       PPP.all.res = all.BT$all.res$Head,
                                       tip.spread = c("Tiliqua_rugosa","Tiliqua_scincoides"),
                                       focus = "clade", psize=3, lsize=2, background.color="grey",
                                       gimme.the.data=F, inset=T)
HEAD.rate <- rate.trajectory.BT.module(tree = egernia.tree, module = HEAD.anc.module,
                                       PPP.all.res = all.BT$all.res$Head,
                                       tip.spread = c("Cyclodomorphus_michaeli","Cyclodomorphus_praealtus"),
                                       focus = "clade", psize=3, lsize=2, background.color="grey",
                                       gimme.the.data=T, inset=T)


# Rate Trajectory and Rate-to-Node for BODY
BODY.anc.module <- data.frame(Interlimb = anc.all.traits$Interlimb, 
                              Body_Width = anc.all.traits$Body_Width, 
                              Pelvic_Width = anc.all.traits$Pelvic_Width, 
                              Pelvic_Height = anc.all.traits$Pelvic_Height)
BODY.plot <- rate.trajectory.BT.module(tree = egernia.tree, module = BODY.anc.module,
                                       PPP.all.res = all.BT$all.res$Body,
                                       tip.spread = c("Tiliqua_rugosa","Cyclodomorphus_maximus"),
                                       focus = "clade", psize=3, lsize=2, background.color="grey",
                                       gimme.the.data=F, inset=T)
BODY.rate <- rate.trajectory.BT.module(tree = egernia.tree, module = BODY.anc.module,
                                       PPP.all.res = all.BT$all.res$Body,
                                       tip.spread = c("Cyclodomorphus_michaeli","Cyclodomorphus_praealtus"),
                                       focus = "clade", psize=3, lsize=2, background.color="grey",
                                       gimme.the.data=T, inset=T)


# Rate Trajectory and Rate-to-Node for TAIL
TAIL.anc.module <- data.frame(Tail_Length = anc.all.traits$Tail_Length, 
                              Tail_Width = anc.all.traits$Tail_Width)
TAIL.plot <- rate.trajectory.BT.module(tree = egernia.tree, module = TAIL.anc.module,
                                       PPP.all.res = all.BT$all.res$Tail,
                                       tip.spread = c("Egernia_depressa","Egernia_cygnitos"),
                                       focus = "clade", psize=3, lsize=2, background.color="grey",
                                       gimme.the.data=F, inset=T)
TAIL.rate <- rate.trajectory.BT.module(tree = egernia.tree, module = TAIL.anc.module,
                                       PPP.all.res = all.BT$all.res$Tail,
                                       tip.spread = c("Cyclodomorphus_michaeli","Cyclodomorphus_praealtus"),
                                       focus = "clade", psize=3, lsize=2, background.color="grey",
                                       gimme.the.data=T, inset=T)


# Rate Trajectory and Rate-to-Node for LIMB
LIMB.anc.module <- data.frame(Upper_Arm = anc.all.traits$Upper_Arm, 
                              Lower_Arm = anc.all.traits$Lower_Arm,
                              Hand = anc.all.traits$Hand,
                              Upper_Leg = anc.all.traits$Upper_Leg,
                              Lower_Leg= anc.all.traits$Lower_Leg,
                              Foot = anc.all.traits$Foot)
LIMB.plot <- rate.trajectory.BT.module(tree = egernia.tree, module = LIMB.anc.module,
                                       PPP.all.res = all.BT$all.res$Limb,
                                       tip.spread = c("Tiliqua_rugosa","Cyclodomorphus_maximus"),
                                       focus = "clade", psize=3, lsize=2, background.color="grey",
                                       gimme.the.data=F, inset=T)
LIMB.rate <- rate.trajectory.BT.module(tree = egernia.tree, module = LIMB.anc.module,
                                       PPP.all.res = all.BT$all.res$Limb,
                                       tip.spread = c("Cyclodomorphus_michaeli","Cyclodomorphus_praealtus"),
                                       focus = "clade", psize=3, lsize=2, background.color="grey",
                                       gimme.the.data=T, inset=T)

rate.list <- list(head=HEAD.rate, body=BODY.rate, tail=TAIL.rate, limb=LIMB.rate)
save(rate.list, file="/Users/ianbrennan/Documents/GitHub/Tiliquini/Data/rate.trajectories.RData")


