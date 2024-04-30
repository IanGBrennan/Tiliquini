library(ggplot2)

source("Scripts/trait.at.time.R")

############################################################################

# load the basic data including the tree
load("Data/Tiliquini_Data.RData")

# load the processed BayesTraits TRAIT data
load("Data/BayesTraits_processed_TRAITS.RData")

############################################################################

# annoying, but to get the right order we have to rename the rate objects
sig.BT <- trait.BT$all.sig

############################################################################

# reorder the objects
for(k in 1:length(sig.BT)){names(sig.BT[[k]])[1:Ntip(til.tree)] <- til.tree$tip.label}
for(k in 1:length(sig.BT)){sig.BT[[k]] <- append(sig.BT[[k]], 0, after=Ntip(til.tree))}
for(k in 1:length(sig.BT)){names(sig.BT[[k]])[[Ntip(til.tree)+1]] <- Ntip(til.tree)+1}

############################################################################

# actually extract the rates-at-time
all.RaT <- lapply(1:length(sig.BT), 
                  function(x) rate.at.time(timeslices=0.1, 
                                           phy=til.tree, 
                                           rate.vector=sig.BT[[x]]))
names(all.RaT) <- names(sig.BT)

############################################################################

# extract the 95% quantile around the estimated rate
all.rates <- lapply(1:length(all.RaT), 
                    function(x) extract.stat(all.RaT[[x]], 
                                             stat="mean", 
                                             plot = F, 
                                             range="quantile"))
names(all.rates) <- names(all.RaT)

############################################################################

# combine the mean and quantiles into a single data frame
all.rates <- lapply(1:length(all.rates), 
                    function(x) all.rates[[x]] <- cbind(all.rates[[x]], trait=names(all.rates)[[x]]))
names(all.rates) <- names(all.RaT)

############################################################################

# reshape into a new dataframe to make plotting easier
all.rates.df <- NULL
for(k in 1:length(all.rates)){all.rates.df <- rbind(all.rates.df, all.rates[[k]])}

############################################################################

# plot the data
ggplot(all.rates.df) +
  geom_ribbon(aes(x=time, ymin=`5%`, ymax=`95%`, fill=trait, alpha=0.5)) +
  geom_line(aes(x=time, y=rate, color=trait)) +
  facet_wrap(~trait, scales="free") + 
  theme_classic() + theme(legend.position="none")
