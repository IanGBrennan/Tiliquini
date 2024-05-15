library("neldermead")
library("fBasics")
library("gsl")
library("stabledist")
library("statmod")

library("pulsR")
library("geiger")
library("parallel")
library("dplyr")
library("ggplot2"); library(patchwork)

##############################################################################

setwd("/Users/ianbrennan/Documents/GitHub/Tiliquini")

source("Scripts/Calculate_AICs.R")

##############################################################################

# load the basic data including the tree
load("Data/Tiliquini_Data.RData")

# load the processed BayesTraits trait data
load("Data/BayesTraits_processed_TRAITS.RData")

##############################################################################


# rename the all.LSR.traits df for ease of use
all.traits <- all.LSR.traits

# establish the list of models we want to fit
model.vec <- c("BM","OU","EB","JN","NIG","BMJN","BMNIG")

# create an empty list to hold the results
trait.fit <- NULL
params.list <- NULL

# make a loop to fit the models to each trait
for (k in 1:ncol(all.traits)){
  # turn the trait into a vector
  trait.vec <- all.traits[,k]
  names(trait.vec) <- rownames(all.traits)
  
  # fit the models in parallel with mclapply
  res.list <- mclapply(1:7, function(x) {
    fit_reml_levy(phy = til.tree, dat = trait.vec, model = model.vec[[x]])}, mc.cores = 8)
  names(res.list) <- model.vec
  trait.fit[[k]] <- unlist(lapply(res.list, function(x) x$AIC))
  params.list[[k]] <- lapply(res.list, function(x) x$params)
}
names(trait.fit) <- colnames(all.traits)
names(params.list) <- colnames(all.traits)

#for(k in 1:length(trait.fit)){trait.fit[[k]] <- trait.fit[[k]][1:7]}

save(trait.fit, file="Data/puLSR_AICs.RData")
save(params.list, file="Data/puLSR_parameters.RData")

##############################################################################

# We need estimate the likelihood of the VR model. We can do this by 
# fitting Brownian Motion to our scaled trees

# get the likelihood by fitting the data to the scaled tree
ScaleTree.AIC <- function(scaled.tree, bt.res, trait){
  # fit a BM model with geiger to the scaled tree
  BT.fit <- fitContinuous(scaled.tree, trait, model="BM")
  # identify shifts as only branches with Median.Scalars >= 2
  no.param <- nrow(dplyr::filter(bt.res, Median.Scalar >= 2)) + 1
  grabAIC(BT.fit$opt$lnL, no.param)
}
BT.AICs <- NULL

BT.AICs[1]  <- ScaleTree.AIC(trait.BT$all.trees$Interlimb, trait.BT$all.res$Interlimb, select(all.traits, Interlimb))
BT.AICs[2]  <- ScaleTree.AIC(trait.BT$all.trees$BodyWidth, trait.BT$all.res$BodyWidth, select(all.traits, BodyWidth))
BT.AICs[3]  <- ScaleTree.AIC(trait.BT$all.trees$PelvicWidth, trait.BT$all.res$PelvicWidth, select(all.traits, PelvicWidth))
BT.AICs[4]  <- ScaleTree.AIC(trait.BT$all.trees$PelvicHeight, trait.BT$all.res$PelvicHeight, select(all.traits, PelvicHeight))
BT.AICs[5]  <- ScaleTree.AIC(trait.BT$all.trees$HeadWidth, trait.BT$all.res$HeadWidth, select(all.traits, HeadWidth))
BT.AICs[6]  <- ScaleTree.AIC(trait.BT$all.trees$SnoutEye, trait.BT$all.res$SnoutEye, select(all.traits, SnoutEye))
BT.AICs[7]  <- ScaleTree.AIC(trait.BT$all.trees$EyeDiameter, trait.BT$all.res$EyeDiameter, select(all.traits, EyeDiameter))
BT.AICs[8]  <- ScaleTree.AIC(trait.BT$all.trees$HeadDepth, trait.BT$all.res$HeadDepth, select(all.traits, HeadDepth))
BT.AICs[9]  <- ScaleTree.AIC(trait.BT$all.trees$TailWidth, trait.BT$all.res$TailWidth, select(all.traits, TailWidth))
BT.AICs[10] <- ScaleTree.AIC(trait.BT$all.trees$UpperArm, trait.BT$all.res$UpperArm, select(all.traits, UpperArm))
BT.AICs[11] <- ScaleTree.AIC(trait.BT$all.trees$LowerArm, trait.BT$all.res$LowerArm, select(all.traits, LowerArm))
BT.AICs[12] <- ScaleTree.AIC(trait.BT$all.trees$Hand, trait.BT$all.res$Hand, select(all.traits, Hand))
BT.AICs[13] <- ScaleTree.AIC(trait.BT$all.trees$UpperLeg, trait.BT$all.res$UpperLeg, select(all.traits, UpperLeg))
BT.AICs[14] <- ScaleTree.AIC(trait.BT$all.trees$LowerLeg, trait.BT$all.res$LowerLeg, select(all.traits, LowerLeg))
BT.AICs[15] <- ScaleTree.AIC(trait.BT$all.trees$Foot, trait.BT$all.res$Foot, select(all.traits, Foot))
BT.AICs[16] <- ScaleTree.AIC(trait.BT$all.trees$Neck, trait.BT$all.res$Neck, select(all.traits, Neck))
BT.AICs[17] <- ScaleTree.AIC(trait.BT$all.trees$PosSkull, trait.BT$all.res$PosSkull, select(all.traits, PosSkull))
BT.AICs[18] <- ScaleTree.AIC(trait.BT$all.trees$PelvicGap, trait.BT$all.res$PelvicGap, select(all.traits, PelvicGap))
BT.AICs[19] <- ScaleTree.AIC(trait.BT$all.trees$TailLength, trait.BT$all.res$TailLength, select(all.traits, TailLength))
BT.AICs[20] <- ScaleTree.AIC(trait.BT$all.trees$Size, trait.BT$all.res$Size, select(all.traits, Size))

names(BT.AICs) <- names(trait.fit)
#names(BT.AICs) <- c("Body_Width","Eye_Diameter","Foot","Hand","Head_Depth","Head_Width","Interlimb","Lower_Arm","Lower_Leg","Neck","Pelvic_Height","Pelvic_Width","Pos_Skull","Size","Snout_Eye","Tail_Length","Tail_Width","Upper_Arm","Upper_Leg")
#BT.AICs[c(4,5,9,15)] <- 100 # all these traits have VR = BM (no shifts, all Median.Scalar=1), so we shouldn't include the VR

# save the BT results
save(BT.AICs, file="Data/BayesTraits_AICs.RData")


##############################################################################

# We will combine the results of the standard and VR model fitting
# summarize and plot them

# if you have results from fitting the VarRates model separately, bring them in here.
for(k in 1:length(trait.fit)){
  curr.trait <- names(trait.fit[k])
  var.rate <- BT.AICs[curr.trait]; names(var.rate) <- "VarRates"
  trait.fit[[k]] <- append(trait.fit[[k]], var.rate)
}

# turn the AIC scores into AIC-weights
trait.aicw.all <- lapply(trait.fit, function(x) aic.w(x))

# then cluster the Jump style models together
trait.fit.types <- NULL
for (j in 1:length(trait.fit)){
  trait.fit.types[[j]] <- trait.fit[[j]][c(1:3, which(trait.fit[[j]]==min(trait.fit[[j]][4:7])),8)]
  names(trait.fit.types[[j]]) <- c("BM","OU","EB","Jump","VarRates")
}
# and turn the AIC scores into AIC-weights again
trait.aicw.types <- lapply(trait.fit.types, function(x) aic.w(x))
names(trait.aicw.types) <- names(trait.fit)

# make a dataframe of the 'type' AICw
aicw.df <- NULL
for (k in 1:length(trait.aicw.types)){
  aicw.df <- rbind(aicw.df, data.frame(trait=names(trait.aicw.types)[[k]], aicw=as.vector(trait.aicw.types[[k]]), model=names(trait.aicw.types[[k]])))
}
# plot the 'type'  results
aicw.types.plot <- ggplot(aicw.df, aes(fill=model, y=aicw, x=trait)) + 
  geom_bar(position="stack", stat="identity") + theme_classic() + 
  scale_fill_brewer(palette="RdYlBu") +
  theme(axis.text.x = element_text(angle = 45, hjust=1), legend.position="bottom")

# make a dataframe of all the model AICw
aicw.all.df <- NULL
for (k in 1:length(trait.aicw.all)){
  aicw.all.df <- rbind(aicw.all.df, data.frame(trait=names(trait.aicw.all)[[k]], aicw=as.vector(trait.aicw.all[[k]]), model=names(trait.aicw.all[[k]])))
}
# plot all the AICw 
aicw.all.plot <- ggplot(aicw.all.df, aes(fill=model, y=aicw, x=trait)) + 
  geom_bar(position="stack", stat="identity") + theme_classic() + 
  scale_fill_brewer(palette="RdYlBu") +
  theme(axis.text.x = element_text(angle = 45, hjust=1), legend.position="bottom")

# save the results externally
save(trait.fit, params.list, trait.aicw.all, trait.fit.types, 
     trait.aicw.all, trait.aicw.types, aicw.df, aicw.all.df, file="Data/ModelFitting_Results.Rdata")

# plot the two figures together
aicw.types.plot / aicw.all.plot

##############################################################################

# report which model is preferred for each trait, where preference is determined
# if the top fitting model has 2x the AICw of the next best model

sig.fit <- NULL
for(k in 1:length(trait.aicw.types)){
  curr.res <- sort(trait.aicw.types[[k]], decreasing=T)
  if(curr.res[[1]]/curr.res[[2]] > 2){sig.fit <- append(sig.fit, names(curr.res)[[1]])}
  if(curr.res[[1]]/curr.res[[2]] <= 2){sig.fit <- append(sig.fit, "none preferred")}
}
names(sig.fit) <- names(trait.aicw.types)
sig.fit[order(names(sig.fit))]

