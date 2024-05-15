# These functions work to:
      # (1) extract trait values on all branches in a tree at specified time intervals (slices).
      # (2) extract trait values on all branches in a tree at specified time intervals (slices) for multiple traits.
      # (3) calculate the variance in trait values across all branches in a tree at specified time intervals (slices).

# As input the first two functions require the phylogeny, and trait values at terminal (tip) and internal (ancestral) nodes.
## These must be estimated previously, under some Brownian Motion process (constant or variable rates; though maybe not?)
## As long as we're willing to accept that rates are constant along individual branches, we can estimate
## trait values along the branches as a fractional step between the ancestral and descendant node trait values.

# The third function takes the output of the either of the first two functions to calculate the variances.
## These can be presented as "cumulative" in which the variance is the total variance of all taxa living at
## a given point in time. I think this helps to visualize the expectation that variance should increase
## with increasing richness (as a result of adaptive divergence). Alternatively we can present variance as
## "relative" in which we correct for the number of taxa alive at a given period. I think this helps to
## visualize temporal patterns in the accumulation of trait diversity. 

## The scripts here are heavily influenced by this: https://github.com/cichlidx/ronco_et_al.

require(ape)
require(phytools)
require(dplyr)
require(dispRity)
require(Rmisc)

trait.at.time <- function(timeslices, phy, trait.vector, plot=F){
  if(!length(trait.vector) == (Ntip(phy) + Nnode(phy))){stop("your trait vector does not include both tip and ancestral trait values")}
  
  # extract the node heights including the root!
  nH <- nodeHeights(phy)
  root <- max(nH[,2])
  
  # round the root height to avoid issues
  root2 <- floor(root*1000)/1000
  
  # define the timeslices for extracting trait values, and apply it to your tree
  ts <- c(seq(from=0, to=root2, by=timeslices), root2)
  
  # make sure the order of the trait vector matches the tree tip and node order
  trait.vector <- trait.vector[match(names(trait.vector), c(phy$tip, ((Ntip(phy)+1):(Ntip(phy) + Nnode(phy)))))]
  
  # combine the trait values with the node ages (including tips)
  ndel <- node.depth.edgelength(phy) # age of each node (tips are largest numbers)
  trait.ages <-  cbind(trait.vector, ndel) # seting reconstructed phenotype with x coordinate (age) of corresponding node
  # new is just head.anc[,1]
  
  # for each timepoint reconstruct the trait value as a function of the distance between the parent and child node
  Yt <- list()
  for (i in 1:(length(ts)-1)) {      
    # get row number for the nodesHeights within root and timepoint
    # get all edges crossing the timeslcie v:tips, but not ending before (not adding nodes which are not present anymore)
    
    # which edges (row number of nH) overlap timeslice i
    curr.edges      <- intersect(which(nH[,1]<= ts[i] ), which(nH[,2] >= ts[i]))
    
    # get start/end node numbers of all edges between root and timepoint 
    curr.edge.nodes <- phy$edge[curr.edges,]
    
    # extract trait values and ages at start/end nodes
    trait.start <- trait.ages[curr.edge.nodes[,1],]
    trait.end <- trait.ages[curr.edge.nodes[,2],]
    
    # apply the current slice time
    slice <- rep(ts[i],times = length(trait.end[,1]))
    
    # estimate the trait values along the edges at timeslice i
    #  (delta trait / delta time) * (edge fragment length) + (original trait value)
    trait <- (((trait.end[,1]-trait.start[,1])/(trait.end[,2]-trait.start[,2])) * 
                (ts[i] - trait.start[,2])) + trait.start[,1]
    
    # add all the trait values at timeslice i to the list
    Yt[[i]] <- cbind(trait, slice)
  }
  
  # remove the fractional period just before the tips (combine the penultimate and ultimate windows)
  ts2 <- ts[-(length(ts)-1)]
  
  # now make a data frame of the trait values at given times
  shuffled.ages <- NULL
  for (k in 1:length(ts2)){
    curr.slice <- unlist(Yt[[k]])
    curr.slice.df <- data.frame(curr.slice, "time" = ts2[k])
    shuffled.ages <- rbind(shuffled.ages, curr.slice.df)
  }
  
  # extract just the trait and timeslices
  trait.time <- shuffled.ages[,c("trait","time")]
  
  # reorder the columns so 'time' comes first
  trait.time <- relocate(trait.time, time)
  
  # plot the values if you're interested, it should look like the phytools phenogram
  if(plot==T){plot(trait.time$trait ~ trait.time$time, xlab="Time", ylab="Trait")}
  
  # give up the object
  return(trait.time)
}

# example:
# trait.at.time(timeslices = 0.25, phy = egernia.tree, trait.vector = head.anc[,1], plot = T)


trait.at.time.multi <- function(timeslices, phy, trait.matrix, plot=F){
  all.traits <- NULL
  
  # let's just make a loop that runs trait.at.time and renames the output.matrix
  for (k in 1:ncol(trait.matrix)){
    trait.vector <- trait.matrix[,k]
    names(trait.vector) <- rownames(trait.matrix)
    curr.trait <- trait.at.time(timeslices=timeslices, phy=phy, trait.vector=trait.vector, plot=F)
    colnames(curr.trait)[2] <- colnames(trait.matrix)[k]
    if(k == 1){all.traits <- curr.trait}
    else{all.traits <- cbind(all.traits, curr.trait[,2])}
  }

  #all.traits <- relocate(all.traits, time)
  colnames(all.traits)[2:ncol(all.traits)] <- colnames(trait.matrix)
  return(all.traits)
}

# example:
# trait.at.time.multi(timeslices = 0.25, phy = egernia.tree, trait.matrix = head.anc, plot = F)

# helper functions for 'extract.variance'. 
# 'disparity' is measured as average pairwise distance between species following Harmon et al. (2003).
# 'functional divergence' follows Villeger et al. (2008)
# 'variance' is simply difference between maximum and minimum trait values
agg.disparity <-  function(x) geiger::disparity(data=as.matrix(c(x)))
agg.func.div <-  function(x) dispRity::func.div(as.matrix(c(x)))
agg.variance <- function(x) abs(min(x)-max(x))

extract.variance <- function(trait.time.obj, plot=c("cumulative", "relative", "sideXside", FALSE), metric=c("variance", "functional divergence", "disparity")){
  # create a data frame of the individual time slices
  time.variance <- data.frame(time = unique(trait.time.obj$time))
  
  if(metric == "variance"){
    time.variance$measure <- aggregate(trait.time.obj[,2], list(trait.time.obj$time), agg.variance)[,2]
  ## establish the minimum and maximum values at each slice
  #min.variance <- aggregate(trait.time.obj[,2], list(trait.time.obj$time), min)
  #max.variance <- aggregate(trait.time.obj[,2], list(trait.time.obj$time), max)
  #
  ## determine the variance in each time slice abs(trait min - trait max)
  #time.variance$variance <- abs(min.variance[,2] - max.variance[,2])
  #
  ## determine the extant variance
    tip.vals <- filter(trait.time.obj, time == max(trait.time.obj$time))
    extant.var <- abs(min(tip.vals[,2]) - max(tip.vals[,2]))
  }
  if(metric == "functional divergence"){
    time.variance$measure <- aggregate(trait.time.obj[,2], list(trait.time.obj$time), agg.func.div)[,2]
    extant.var <- max(time.variance$measure)
  }
  if(metric == "disparity"){
    time.variance$measure <- aggregate(trait.time.obj[,2], list(trait.time.obj$time), agg.disparity)[,2]
    extant.var <- max(time.variance$measure)
  }

  # establish the number of species living in each time slice
  spp.slice <- aggregate(trait.time.obj[,2], list(trait.time.obj$time), length)
  
  # add the richness information to the data frame
  time.variance$richness <- spp.slice[,2]
  
  # correct variance by the number of species living at that time
  time.variance$measure.rich <- time.variance$measure/time.variance$richness
  
  # get the name of the current variable
  var.name <- colnames(trait.time.obj)[2]
  
  # determine at what time 50% of the morphological diversity had accumulated
  time.50 <- time.variance[which.min(abs(time.variance$measure - (extant.var/2))),"time"]
  
  # plot the results if you'd like
  if(plot=="cumulative"){plot(time.variance$measure ~ time.variance$time, 
                              xlab="Time", ylab=paste("Total", metric), type="l", main=var.name)
                         abline(v=time.50, col="#ABDDA4", lty=3, lwd=3)}
  if(plot=="relative"){plot(time.variance$measure.rich ~ time.variance$time, 
                            xlab="Time", ylab=paste("Relative", metric), type="l", main=var.name)}
  if(plot=="sideXside"){layout(matrix(1:2,ncol=2))
    plot(time.variance$measure ~ time.variance$time, 
         xlab="Time", ylab=paste("Total", metric), type="l", main=var.name)
    abline(v=time.50, col="#ABDDA4", lty=3, lwd=3)
    plot(time.variance$measure.rich ~ time.variance$time, 
         xlab="Time", ylab=paste("Relative", metric), type="l", main=var.name)
    layout(matrix(1))}
  if(plot==F){}
  
  return(time.variance)
}

# example:
# extract.variance(all.anc[,c(1,5)], plot="sideXside")


# get the rates at timeslices from a processed dataframe
rate.at.time.df <- function(timeslices, obj, plot=F, relative.rates=c("F","mean","median")){
  # define the timeslices for extracting trait values, and apply it to your tree
  ts <- c(seq(from=min(c(obj$timestart,obj$timestop)), to=max(c(obj$timestart,obj$timestop)), by=timeslices), max(c(obj$timestart,obj$timestop)))
  
  # for each timepoint reconstruct the trait value as a function of the distance between the parent and child node
  Yt <- list()
  for (i in 1:(length(ts)-1)) {      
    # which edges overlap timeslice i
    curr.edges <- filter(obj, timestart <= ts[[i]] & timestop >= ts[[i]])
    
    # extract the rates for the branches of interest
    curr.rates <- curr.edges$edge.rate
    
    # apply the current slice time
    slice <- rep(ts[i],times = length(curr.rates))
    
    # add all the trait values at timeslice i to the list
    Yt[[i]] <- cbind(curr.rates, slice)
  }
  # remove the fractional period just before the tips (combine the penultimate and ultimate windows)
  ts2 <- ts[-(length(ts)-1)]
  
  # now make a data frame of the trait values at given times
  shuffled.ages <- NULL
  for (k in 1:length(ts2)){
    curr.slice <- unlist(Yt[[k]])
    curr.slice.df <- data.frame(curr.slice, "time" = ts2[k])
    shuffled.ages <- rbind(shuffled.ages, curr.slice.df)
  }
  
  # extract just the trait and timeslices
  rate.time <- shuffled.ages[,c("curr.rates","time")]
  
  # reorder the columns so 'time' comes first
  rate.time <- relocate(rate.time, time)
  colnames(rate.time) <- c("time", "rate")

  # rescale rates relative to the mean if requested
  if(relative.rates=="F"){rate.time$rate <- rate.time$rate}
  if(relative.rates=="mean"){rate.time$rate <- rate.time$rate/mean(rate.time$rate)}
  if(relative.rates=="median"){rate.time$rate <- rate.time$rate/median(rate.time$rate)}
  
  # plot the values if you're interested, it should look like the horizontal branches of the tree
  if(plot==T){plot(rate.time$rate ~ rate.time$time, xlab="Time", ylab="Rate", pch=16)}
  
  # give up the object
  return(rate.time)
}


# a rate.matrix is just the 'multiple.rates' value of an RRphylo output (e.g. RRhead$multiple.rates)
rate.at.time <- function(timeslices, phy, rate.vector, plot=F){
  # RRphylo rates are ridge regression coefficients (so - & +) so we want the absolute value
  #rate.vector <- abs(rate.vector)
  
  if(!length(rate.vector) == (Ntip(phy) + Nnode(phy))){stop("your rate vector does not include both tip and ancestral trait values")}
  
  # make sure the matrix order matches the tip + node order
  rate.vector <- order.by.phylo(rate.vector, phy)
  names(rate.vector) <- 1:length(rate.vector)
  
  # extract the node heights including the root!
  nH <- nodeHeights(phy)
  root <- max(nH[,2])
  
  # round the root height to avoid issues
  root2 <- floor(root*1000)/1000
  
  # define the timeslices for extracting trait values, and apply it to your tree
  ts <- c(seq(from=0, to=root2, by=timeslices), root2)
  
  # combine the trait values with the node ages (including tips)
  ndel <- node.depth.edgelength(phy) # age of each node (tips are largest numbers)
  rate.ages <-  cbind(rate.vector, ndel) # seting reconstructed phenotype with x coordinate (age) of corresponding node
  # new is just head.anc[,1]
  
  # for each timepoint reconstruct the trait value as a function of the distance between the parent and child node
  Yt <- list()
  for (i in 1:(length(ts)-1)) {      
    # get row number for the nodesHeights within root and timepoint
    # get all edges crossing the timeslcie v:tips, but not ending before (not adding nodes which are not present anymore)
    
    # which edges (row number of nH) overlap timeslice i
    curr.edges <- intersect(which(nH[,1]<= ts[i] ), which(nH[,2] >= ts[i]))
    
    # get start/end node numbers of all edges between root and timepoint 
    curr.edge.nodes <- phy$edge[curr.edges,]
    
    # extract the rates for the branches of interest
    curr.rates <- rate.vector[curr.edge.nodes[,2]]
    
    # extract rate values and ages at start/end nodes
    rate.start <- rate.ages[curr.edge.nodes[,1],]
    rate.end <- rate.ages[curr.edge.nodes[,2],]
    
    # apply the current slice time
    slice <- rep(ts[i],times = length(rate.end[,1]))
    
    # add all the trait values at timeslice i to the list
    Yt[[i]] <- cbind(curr.rates, slice)
  }
  
  # remove the fractional period just before the tips (combine the penultimate and ultimate windows)
  ts2 <- ts[-(length(ts)-1)]
  
  # now make a data frame of the trait values at given times
  shuffled.ages <- NULL
  for (k in 1:length(ts2)){
    curr.slice <- unlist(Yt[[k]])
    curr.slice.df <- data.frame(curr.slice, "time" = ts2[k])
    shuffled.ages <- rbind(shuffled.ages, curr.slice.df)
  }
  
  # extract just the trait and timeslices
  rate.time <- shuffled.ages[,c("curr.rates","time")]
  
  # reorder the columns so 'time' comes first
  rate.time <- relocate(rate.time, time)
  colnames(rate.time) <- c("time", "rate")
  
  # plot the values if you're interested, it should look like the horizontal branches of the tree
  if(plot==T){plot(rate.time$rate ~ rate.time$time, xlab="Time", ylab="Rate", pch=16)}
  
  # give up the object
  return(rate.time)
}

# example:
# rate.at.time(timeslices = 0.25, phy = egernia.tree, rate.vector = RRhead$multiple.rates[,4], plot=T)
# we can do an experiment to make sure it's working:
# rvec <- rep(1, 93)
# names(rvec) <- rownames(rate.matrix)
# rvec[which(names(rvec)==77)] <- 2
# rate.at.time(timeslices = 0.25, phy = egernia.tree, rate.vector = rvec, plot = T)




rate.at.time.multi <- function(timeslices, phy, rate.matrix, plot=F){
  all.rates <- NULL
  
  # let's just make a loop that runs trait.at.time and renames the output.matrix
  for (k in 1:ncol(rate.matrix)){
    rate.vector <- rate.matrix[,k]
    curr.rate <- rate.at.time(timeslices=timeslices, phy=phy, rate.vector=rate.vector, plot=F)
    colnames(curr.rate)[2] <- colnames(rate.matrix)[k]
    if(k == 1){all.rates <- curr.rate}
    else{all.rates <- cbind(all.rates, curr.rate[,2])}
  }
  
  #all.traits <- relocate(all.traits, time)
  colnames(all.rates)[2:ncol(all.rates)] <- colnames(rate.matrix)
  return(all.rates)
}

# example:
# rate.at.time.multi(timeslices = 0.25, phy = egernia.tree, rate.matrix = RRhead$multiple.rates, plot=F)

extract.stat <- function(rate.time.obj, stat=c("mean", "median", "scale"), range=c("confidence","quantile"),
                         plot=c("average", "corrected", "sideXside", FALSE)){
  
  # create a data frame of the individual time slices
  time.mean <- data.frame(time = unique(rate.time.obj$time))
  
  # establish the minimum and maximum values at each slice
  if(stat=="mean"){rate.mean <- aggregate(rate.time.obj[,2], list(rate.time.obj$time), mean)}
  if(stat=="median"){rate.mean <- aggregate(rate.time.obj[,2], list(rate.time.obj$time), median)}
  if(stat=="scale"){
    rate.mean <- aggregate(rate.time.obj[,2], list(rate.time.obj$time), mean)
    rate.mean[,2] <- ((rate.mean[,2] - min(rate.mean[,2]))/diff(range(rate.mean[,2])) * 99) + 1
  }
  rate.sd <- aggregate(rate.time.obj[,2], list(rate.time.obj$time), sd)
  

  # if you want to estimate quantiles on the rates:
  if(range=="quantile"){
    rate.qt5 <- aggregate(rate.time.obj[,2], list(rate.time.obj$time), FUN='quantile', probs=0.025)
    rate.qt95 <- aggregate(rate.time.obj[,2], list(rate.time.obj$time), FUN='quantile', probs=0.975)
    rate.mean <- cbind(rate.mean, rate.sd[,2], rate.qt5[,2], rate.qt95[,2])    
  }

  # if you prefer to estimate confidence intervals on the rates:
  if(range=="confidence"){
   rate.all <- aggregate(rate.time.obj[,2], list(rate.time.obj$time), FUN='CI', ci=0.95)
   rate.mean <- data.frame(time = time.mean[,1], rate = rate.mean[,2], 
                           sd = rate.sd[,2], "5%" = rate.all$x[,3], "95%" = rate.all$x[,1])   
  }
  # regardless, rename the columns
  colnames(rate.mean) <- c("time", "rate", "sd", "5%", "95%")

  # establish the number of species living in each time slice
  spp.slice <- aggregate(rate.time.obj[,2], list(rate.time.obj$time), length)
  
  # add the richness information to the data frame
  rate.mean$richness <- spp.slice[,2]
  
  # correct rate by the number of species living at that time
  rate.mean$rate.rich <- rate.mean$rate/rate.mean$richness
  
  # get the name of the current variable
  var.name <- colnames(rate.time.obj)[2]
  
  # plot the results if you'd like
  if(plot=="average"){plot(rate.mean$rate ~ rate.mean$time, 
                              xlab="Time", ylab=paste(stat, "Evolutionary Rate"), type="l")}
  if(plot=="corrected"){plot(rate.mean$rate.rich ~ rate.mean$time, 
                            xlab="Time", ylab=paste(stat, "Evolutionary Rate"), type="l")}
  if(plot=="sideXside"){layout(matrix(1:2,ncol=2))
    plot(rate.mean$rate ~ rate.mean$time, 
         xlab="Time", ylab=paste(stat, "Evolutionary Rate"), type="l")
    plot(rate.mean$rate.rich ~ rate.mean$time, 
         xlab="Time", ylab=paste(stat, "Evolutionary Rate"), type="l")
    layout(matrix(1))}
  if(plot==F){}
  
  return(rate.mean)
  
}

# example:
#extract.mean(rate.time.obj = testo[,c(1,5)], stat = "mean", plot = "sideXside")

# reorder any dataframe or matrix following the phylogenetic order
order.by.phylo <- function(jumbled.order, phy){
  # make sure the matrix order matches the tip + node order
  
  # if the data is in a vector (single trait)
  if(is.vector(jumbled.order)){
    phylo.order <- jumbled.order[order(match(names(jumbled.order), 
                                             c(phy$tip.label, (length(phy$tip.label)+1):(Ntip(phy)+Nnode(phy)))))]
  }
  # if the data is in a matrix or data frame (multiple traits)
  if(!is.vector(jumbled.order)){
  phylo.order <- jumbled.order[order(match(rownames(jumbled.order), 
                                         c(phy$tip.label, (length(phy$tip.label)+1):(Ntip(phy)+Nnode(phy))))),]
  }
  return(phylo.order)
}

# example:
# order.by.phylo(rate.matrix, egernia.tree)

extract.disparity <- function(trait.time.multi.obj, plot=T, 
  metric=c("func.div", "variances", "quantiles", "func.eve", "disparity")){
  
  trait.time.obj <- trait.time.multi.obj
  
  # create a data frame of the individual time slices
  time.mean <- data.frame(time = unique(trait.time.obj$time))
  
  # get the time points
  times <- unique(trait.time.obj$time)
  
  # make a loop to establish functional diversity at each timeslice
  trait.div <- NULL
  if(metric == "func.div"){for  (p in 1:length(times)){trait.div[p] <- dispRity::func.div(filter(trait.time.obj, time == times[[p]]))}}
  if(metric == "func.eve"){for  (p in 1:length(times)){trait.div[p] <- dispRity::func.eve(filter(trait.time.obj, time == times[[p]]))}}
  if(metric == "variances"){for (p in 1:length(times)){trait.div[p] <- variances(filter(trait.time.obj, time == times[[p]]))}}
  if(metric == "quantiles"){for (p in 1:length(times)){trait.div[p] <- quantiles(filter(trait.time.obj, time == times[[p]]))}}
  if(metric == "disparity"){for (p in 1:length(times)){trait.div[p] <- geiger::disparity(data=filter(trait.time.obj, time == times[[p]]))}}
  # geiger::disparity(x) is equal to mean(dispRity::pairwise.dist(x)^2)
  # note: 'variances' and 'quantiles' won't work here because they are applied to each trait individually
  trait.div[which(is.na(trait.div))] <- 1 # correct any NaN to 1
  
  # combine the times and disparity values back together
  trait.div <- data.frame(time = times, measure = trait.div)
  
  # plot the results if you'd like
  ylabel <- paste0("Disparity (", metric, ")")
  if(plot==T){plot(trait.div$measure ~ trait.div$time, 
                           xlab="Time", ylab=ylabel, type="l")}
  if(plot==F){}
  
  return(trait.div)
}

# example:
# all.limb <- as.matrix(rbind(all.modules$limb, allRR$limb$aces))
# testo <- trait.at.time.multi(timeslices = 0.15, phy = egernia.tree, trait.matrix = all.limb, plot = F)
# extract.disparity(testo, plot=T, metric="func.div")
# * previously called 'extract.func.div'


multiplot.module <- function(module.list, phy){
  for(k in 1:length(module.list)){
    par(mfrow=c(1,(length(module.list[[k]]))+1))
    plotTree.barplot(phy, 
                   setNames(module.list[[k]][,1], rownames(module.list[[k]])), 
                   tip.labels = T, add=T,
                   args.barplot=list(xlab=names(module.list[[k]][1])))
    for(j in 2:length(module.list[[k]])){
      plotTree.barplot(phy, 
                     setNames(module.list[[k]][,j], rownames(module.list[[k]])),
                     args.plotTree=list(plot=FALSE),
                     args.barplot=list(xlab=names(module.list[[k]][j])),
                     tip.labels = T, add=T)
    }
  }
}
# example:
# multiplot.module(all.modules, egernia.tree)
# will produce 4 plots (one for each module in all.modules)