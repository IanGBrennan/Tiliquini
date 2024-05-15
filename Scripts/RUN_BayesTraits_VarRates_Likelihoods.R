source("/Scripts/plotting_BayesTraits.R")
source("/Scripts/trait.at.time.R")
source("/Scripts/rate.trajectory.R")
source("/Scripts/Calculate_AICs.R")
library(dplyr)

# summarize the BayesTraits results for each morphological trait
setwd("/Applications/BayesTraitsV4/PPPostProcess")
in.files <- dir(getwd(), pattern="_TraitResults.txt") # will give you the files for individual traits
in.files <- dir(getwd(), pattern="_ModuleResults.txt") # will give you the files for modules

# process the outputs of BayesTraits
all.BT <- process_PPP(res.files = in.files, phy = til.tree, col.palette="YlOrRd")
# the above will give you a series of trees and res objects
# all.BT$all.trees holds a tree for each trait. The tree is your input tree, rescaled by the mean estimated sigma (rate) value per branch (Mean.SigV in BayesTraits terms)
# all.BT$all.sig holds the mean sigma value per branch in a vector (?), this is the value each branch above was scaled by
# all.BT$all.res holds a data frame for each trait including the node and branch information, rate info, etc.
# all.BT$scalar.trees holds a tree for each trait. The tree is your input tree, rescaled by the median scalar r (Median.Scalar in BayesTraits terms)
# all.BT$mean.scalar.trees holds a tree for each trait. The tree is your input tree, rescaled by the mean scalar r (Mean.Scalar in BayesTraits terms)
# all.BT$rate.trees holds a tree for each trait. The tree is your input tree, with each branch rescaled to the mean estimated sigma value (Mean.SigV)
# that's a lot of trees. You can plot each to see what they show you, but ultimately I think you'll want to work with either the mean or median scalar trees (beetleVR$mean.scalar.trees or beetleVR$scalar.trees). Median is slightly more conservative. 
save(all.BT, file="Data/BayesTraits_processed_modules.RData")

# put all the traits into a single dataframe
all.traits <- bind_cols(all.modules$head, all.modules$body, all.modules$tail, all.modules$limb,
                        all.modules$size, all.modules$neck, all.modules$eye)

# get the likelihood by fitting the data to the scaled tree
ScaleTree.AIC <- function(scaled.tree, bt.res, trait){
  # fit a BM model with geiger to the scaled tree
  BT.fit <- fitContinuous(scaled.tree, trait, model="BM")
  # identify shifts as only branches with Median.Scalars >= 2
  no.param <- nrow(dplyr::filter(bt.res, Median.Scalar >= 2)) + 1
  grabAIC(BT.fit$opt$lnL, no.param)
}
BT.AICs <- NULL

BT.AICs[1]  <- ScaleTree.AIC(all.BT$all.trees$BodyWidth, all.BT$all.res$BodyWidth, select(all.traits, Body_Width))
BT.AICs[2]  <- ScaleTree.AIC(all.BT$all.trees$EyeDiameter, all.BT$all.res$EyeDiameter, select(all.traits, Eye_Diameter))
BT.AICs[3]  <- ScaleTree.AIC(all.BT$all.trees$Foot, all.BT$all.res$Foot, select(all.traits, Foot))
BT.AICs[4]  <- ScaleTree.AIC(all.BT$all.trees$Hand, all.BT$all.res$Hand, select(all.traits, Hand))
BT.AICs[5]  <- ScaleTree.AIC(all.BT$all.trees$HeadDepth, all.BT$all.res$HeadDepth, select(all.traits, Head_Depth))
BT.AICs[6]  <- ScaleTree.AIC(all.BT$all.trees$HeadWidth, all.BT$all.res$HeadWidth, select(all.traits, Head_Width))
BT.AICs[7]  <- ScaleTree.AIC(all.BT$all.trees$Interlimb, all.BT$all.res$Interlimb, select(all.traits, Interlimb))
BT.AICs[8]  <- ScaleTree.AIC(all.BT$all.trees$LowerArm, all.BT$all.res$LowerArm, select(all.traits, Lower_Arm))
BT.AICs[9]  <- ScaleTree.AIC(all.BT$all.trees$LowerLeg, all.BT$all.res$LowerLeg, select(all.traits, Lower_Leg))
BT.AICs[10] <- ScaleTree.AIC(all.BT$all.trees$Neck, all.BT$all.res$Neck, select(all.traits, Neck))
BT.AICs[11] <- ScaleTree.AIC(all.BT$all.trees$PelvicHeight, all.BT$all.res$PelvicHeight, select(all.traits, Pelvic_Height))
BT.AICs[12] <- ScaleTree.AIC(all.BT$all.trees$PelvicWidth, all.BT$all.res$PelvicWidth, select(all.traits, Pelvic_Width))
BT.AICs[13] <- ScaleTree.AIC(all.BT$all.trees$PosSkull, all.BT$all.res$PosSkull, select(all.traits, Pos_Skull))
BT.AICs[14] <- ScaleTree.AIC(all.BT$all.trees$Size, all.BT$all.res$Size, select(all.traits, Size))
BT.AICs[15] <- ScaleTree.AIC(all.BT$all.trees$SnoutEye, all.BT$all.res$SnoutEye, select(all.traits, Snout_Eye))
BT.AICs[16] <- ScaleTree.AIC(all.BT$all.trees$TailLength, all.BT$all.res$TailLength, select(all.traits, Tail_Length))
BT.AICs[17] <- ScaleTree.AIC(all.BT$all.trees$TailWidth, all.BT$all.res$TailWidth, select(all.traits, Tail_Width))
BT.AICs[18] <- ScaleTree.AIC(all.BT$all.trees$UpperArm, all.BT$all.res$UpperArm, select(all.traits, Upper_Arm))
BT.AICs[19] <- ScaleTree.AIC(all.BT$all.trees$UpperLeg, all.BT$all.res$UpperLeg, select(all.traits, Upper_Leg))
names(BT.AICs) <- c("Body_Width","Eye_Diameter","Foot","Hand","Head_Depth","Head_Width","Interlimb","Lower_Arm","Lower_Leg","Neck","Pelvic_Height","Pelvic_Width","Pos_Skull","Size","Snout_Eye","Tail_Length","Tail_Width","Upper_Arm","Upper_Leg")
BT.AICs[c(4,5,9,15)] <- 100 # all these traits have VR = BM (no shifts, all Median.Scalar=1), so we shouldn't include the VR

# save the BT results
save(BT.AICs, file="Data/BayesTraits_AICs.RData")

# if you'd like to try the kappa model which is supposedly like punctuated equilibrium
kappa.AIC <- NULL
for (j in 1:length(traits.vec)){kappa.AIC[[j]] <- fitContinuous(egernia.tree, traits.vec[[j]], model="kappa")$opt$aic}
names(kappa.AIC) <- names(traits.vec)








