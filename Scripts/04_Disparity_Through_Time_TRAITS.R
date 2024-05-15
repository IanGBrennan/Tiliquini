library(ggplot2)


##############################################################################

# source necessary script
source("Scripts/trait.at.time.R")

# load the ancestral trait data
load("Data/Tiliquini_Ancestral_TRAITS.Rdata")

##############################################################################


# extract trait disparity at each point in time
emp.dis <- lapply(1:length(all.anc), function(x) extract.variance(all.anc[[x]], plot = FALSE, metric = "disparity"))
names(emp.dis) <- names(all.anc)

# add in information about which trait we're looking at
emp.dis <- lapply(1:length(emp.dis), function(x) emp.dis[[x]] <- cbind(emp.dis[[x]], trait=names(emp.dis)[[x]]))
names(emp.dis) <- names(all.anc)

# convert these data into 'long' format for easier plotting
emp.dis.df <- NULL
for(k in 1:length(emp.dis)){emp.dis.df <- rbind(emp.dis.df, emp.dis[[k]])}

# plot the trends for each individual trait
ggplot(emp.dis.df) +
  geom_line(aes(x=time, y=measure, color=trait)) +
  facet_wrap(~trait, scales="free") + 
  theme_classic() + theme(legend.position="none")