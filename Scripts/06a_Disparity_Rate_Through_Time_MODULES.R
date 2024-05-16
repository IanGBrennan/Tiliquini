source("Scripts/disparity.to.BM.R")
source("Scripts/trait.at.time.R")

############################################################################

# load the basic data including the tree
load("Data/Tiliquini_Data.RData")

# load the scalar trees
load("Data/BayesTraits_TREES.RData")

# load the processed BayesTraits MODULE data
load("Data/BayesTraits_processed_MODULES.RData")

############################################################################

# summarize trends in disparity through time for each module
# while simulating data based on estimated parameters

# HEAD
head.corBM <- disparity.to.BM.module(phy=til.tree, scaled.phy=module.trees$Head, 
                                     trait.df=all.modules$head, sim.num=100, metric="variance", 
                                     trait.name="Head", loess.span=0.1, window.size=20, trait.correlation=T)
# BODY
body.corBM <- disparity.to.BM.module(phy=til.tree, scaled.phy=module.trees$Body, 
                                     trait.df=all.modules$body, sim.num=100, metric="variance", 
                                     trait.name="Body", loess.span=0.1, window.size=20, trait.correlation=T)
# TAIL
tail.corBM <- disparity.to.BM.module(phy=til.tree, scaled.phy=module.trees$Tail, 
                                     trait.df=all.modules$tail, sim.num=100, metric="variance", 
                                     trait.name="Tail", loess.span=0.1, window.size=20, trait.correlation=T)
# LIMB
limb.corBM <- disparity.to.BM.module(phy=til.tree, scaled.phy=module.trees$Limb, 
                                     trait.df=all.modules$limb, sim.num=100, metric="variance", 
                                     trait.name="Limb", loess.span=0.1, window.size=20, trait.correlation=T)

# save these data to file
save(head.corBM, body.corBM, tail.corBM, limb.corBM, file="Data/Disparity_processed_MODULES.RData")

# have a quick look at the trends
(head.corBM$disparity | head.corBM$slopes | plot_spacer()) /
  (limb.corBM$disparity | limb.corBM$slopes | plot_spacer()) /
  (tail.corBM$disparity | tail.corBM$slopes | plot_spacer()) /
  (body.corBM$disparity | body.corBM$slopes | plot_spacer())

############################################################################

# Now combine the separate disparity summaries (empirical, simulated, slopes) into workable dataframes.
all.disparity <- rbind(cbind(limb.corBM$empirical.metric, trait="limb"),
                       cbind(head.corBM$empirical.metric, trait="head"),
                       cbind(body.corBM$empirical.metric, trait="body"),
                       cbind(tail.corBM$empirical.metric, trait="tail"))
all.sim.disparity <- rbind(cbind(limb.corBM$simulated.metric, trait="limb"),
                           cbind(head.corBM$simulated.metric, trait="head"),
                           cbind(body.corBM$simulated.metric, trait="body"),
                           cbind(tail.corBM$simulated.metric, trait="tail"))
all.slope <- rbind(cbind(limb.corBM$slopes.df, trait="limb"),
                   cbind(head.corBM$slopes.df, trait="head"),
                   cbind(body.corBM$slopes.df, trait="body"),
                   cbind(tail.corBM$slopes.df, trait="tail"))

############################################################################

# Plot the empirical vs. simulated disparity
disparity.plot <- ggplot() +
  #geom_ribbon(data=all.sim.disparity, aes(x=rev(time), ymin=`5%`, ymax=`95%`, fill=trait), alpha=0.1) +
  #scale_fill_brewer(palette="Spectral") +
  geom_line(data=all.sim.disparity, aes(x=rev(time), y=`50%`, color=trait), alpha=0.75, linetype="dotted") +
  geom_line(data=all.disparity, aes(x=rev(time), y=measure, color=trait)) +
  scale_color_brewer(palette="Spectral") + xlab("Time (Ma)") + ylab("Disparity (variance)") +
  scale_x_reverse() + theme_classic() + theme(legend.position="none")

############################################################################

# Plot the difference in slopes of empirical vs. simulated disparity
library(ggbreak)
slope.plot <- ggplot() +
  geom_abline(slope=0, intercept=0, color="black", linetype="dotted") +
  geom_ribbon(data=all.slope, aes(x=rev(time),ymin=slope.lower,ymax=slope.upper, fill=trait), alpha=0.1) +
  scale_fill_brewer(palette="Spectral") +
  geom_line(data=all.slope, aes(x=rev(time), y=slope.diff, color=trait)) +
  scale_color_brewer(palette="Spectral") + 
  xlab("Time (Ma)") + ylab("Difference in Slopes (observed - BM") +
  scale_x_reverse() + theme_classic() + theme(legend.position="none")

slope.plot <- slope.plot + scale_y_break(c(0.75,3),scales=0.5)


############################################################################

# Do a similar exercise to extract rate-through-time information

# get a rate at time object from a dataframe
RaT <- lapply(module.BT$all.res, function(x) rate.at.time.df(timeslices=0.1, obj=x, plot=F, relative.rates="mean"))
# extract the rates through time
RaT <- lapply(RaT, function(x) extract.stat(x, stat="mean", range="quantile", plot=F))
# add a column with the module name
RaT$Body$module <- "body"; RaT$Head$module <- "head"; RaT$Limb$module <- "limb"; RaT$Tail$module <- "tail"

# reshape the data to make it easier to plot
RaT.df <- NULL
for(k in 1:length(RaT)){RaT.df <- rbind(RaT.df, RaT[[k]])}

# plot the rates through time
rate.plot <- ggplot(RaT.df) +
               geom_hline(yintercept=1, linetype="dotted", color="grey") +
               geom_line(aes(x=-time, y=rate, color=module)) +
               scale_color_brewer(palette="Spectral") +
               xlab("Million Years Ago") +
               ylim(c(0.5,2.5)) +
               scale_x_reverse() + 
               #xlim(c(-40,0)) +
               theme_classic() +
               theme(legend.position="bottom")

############################################################################

# Combine the plots
(disparity.plot + slope.plot + rate.plot) + plot_layout(guides="collect") & theme(legend.position="bottom")
  

############################################################################

# Alternatively plot the same data but split by each module

# Plot the empirical vs. simulated disparity
disparity.plot2 <- ggplot() +
  #geom_ribbon(data=all.sim.disparity, aes(x=rev(time), ymin=`5%`, ymax=`95%`, fill=trait), alpha=0.1) +
  #scale_fill_brewer(palette="Spectral") +
  geom_line(data=all.sim.disparity, aes(x=rev(time), y=`50%`), alpha=0.75, linetype="dotted") +
  geom_line(data=all.disparity, aes(x=rev(time), y=measure)) +
  scale_color_brewer(palette="Spectral") + xlab("Time (Ma)") + ylab("Disparity (variance)") +
  scale_x_reverse() + theme_classic() + theme(legend.position="bottom") +
  facet_grid(rows=vars(trait), scales="free")

# Plot the difference in slopes of empirical vs. simulated disparity
slope.plot2 <- ggplot() +
  geom_abline(slope=0, intercept=0, color="black", linetype="dotted") +
  geom_ribbon(data=all.slope, aes(x=rev(time),ymin=slope.lower,ymax=slope.upper), alpha=0.1) +
  scale_fill_brewer(palette="Spectral") +
  geom_line(data=all.slope, aes(x=rev(time), y=slope.diff)) +
  scale_color_brewer(palette="Spectral") + 
  xlab("Time (Ma)") + ylab("Difference in Slopes (observed - BM") +
  scale_x_reverse() + theme_classic() + theme(legend.position="bottom") +
  facet_grid(rows=vars(trait), scales="free")

# plot the rates through time
rate.plot2 <-   ggplot(RaT.df) +
  geom_hline(yintercept=1, linetype="dotted", color="grey") +
  #geom_ribbon(aes(x=time-max(time), ymin=rate-sd, ymax=rate+sd), alpha=0.25) +
  #geom_ribbon(aes(x=time-max(time), ymin=`5%`, ymax=`95%`, fill=module), alpha=0.25) +
  geom_line(aes(x=-time, y=rate)) +
  scale_fill_brewer(palette="Spectral") +
  scale_color_brewer(palette="Spectral") +
  xlab("Million Years Ago") +
  scale_x_reverse() +
  #ylim(c(0,2.5)) +
  facet_grid(scales="free", rows=vars(module)) +
  theme_classic()

(disparity.plot2 + slope.plot2 + rate.plot2) + plot_layout(guides="collect") & theme(legend.position="bottom")

