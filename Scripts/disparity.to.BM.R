#source("Tiliquini/Scripts/trait.at.time.R")

require("mvMORPH")

disparity.to.BM <- function(phy, scaled.phy, trait, sim.num, 
                            metric=c("disparity","variance"), trait.name,
                            loess.span, window.size){
  # fit BM to get parameter estimates
  sim.fit <- geiger::fitContinuous(phy, trait, model="BM")
  # simulate traits
  sim.sim <- lapply(1:sim.num, function(x) fastBM(phy, a=sim.fit$opt$z0, sig2=sim.fit$opt$sigsq, internal=T))
  # extract trait at time object for each simulation
  #sim.tat <- mclapply(sim.sim, function(x) trait.at.time(timeslices=0.1, phy=egernia.tree, trait=x, plot=F), mc.cores=8)
  sim.tat <- lapply(sim.sim, function(x) trait.at.time(timeslices=0.1, phy=phy, trait=x, plot=F))
  # extract your disparity metric at each time slice
  sim.vat <- lapply(sim.tat, function(x) extract.variance(x, plot=F, metric=metric))
  # pull all our values into a dataframe
  sim.var <- NULL
  for(t in 1:length(sim.vat)){sim.var <- cbind(sim.var, sim.vat[[t]][,2])}
  # get the quantiles at each time
  qts <- data.frame(t(apply(sim.var, 1, function(x) quantile(x, probs=c(0.05,0.5,0.95)))))
  # add the time back in
  sim.qts <- data.frame(time=sim.vat[[1]][,1], qts)
  # correct the column names
  colnames(sim.qts) <- c("time","5%","50%","95%")
  
  # combined empirical trait with estimated ancestral nodes from the scaled tree
  emp.trait <- c(trait, fastAnc(scaled.phy, trait))
  # extract trait at time object for the empirical data
  emp.tat <- trait.at.time(timeslices=0.1, phy=phy, trait=emp.trait, plot=F)
  # extract your disparity metric at each time slice
  emp.vat <- extract.variance(emp.tat, plot=F, metric=metric)
  
  # plot the simulated and empirical accumulation of your disparity metric
  accumulation.plot <- ggplot() +
    geom_ribbon(data=sim.qts, aes(x=max(time)-time, ymin=`5%`, ymax=`95%`), fill="grey", alpha=0.25) +
    geom_line(data=sim.qts, aes(x=max(time)-time, y=`50%`), color="darkGrey", linetype="dotted") +
    geom_line(data=emp.vat, aes(x=max(time)-time, y=measure), color="black", lwd=1.5) +
    geom_segment(data=sim.qts, aes(x=max(time)*0.9, xend=max(time), 
                                   y=`50%`[nrow(sim.qts)], yend=`50%`[nrow(sim.qts)]), 
                 arrow=arrow(length=unit(0.3,"cm")), color="darkGrey") + 
    geom_segment(data=emp.vat, aes(x=max(time)*0.9, xend=max(time), 
                                   y=measure[nrow(emp.vat)], yend=measure[nrow(emp.vat)]), 
                 arrow=arrow(length=unit(0.3,"cm")), color="black") +
    scale_x_reverse() + xlab("Time (Ma)") + ylab(paste("Disparity", paste0("(",metric,")"))) +
    labs(subtitle=trait.name) + theme_classic()
  
  # now we move on to the difference in slopes
  
  # get the slope through time of our empirical data
  curr.slope <- extract.empirical.slopes(emp.vat, window.size=window.size)
  
  # and the slopes through time of our simulated data
  slopes.list <- extract.empirical.slopes(sim.vat, window.size=window.size)
  sim.slopes <- data.frame(matrix(nrow=nrow(slopes.list[[1]]),ncol=length(slopes.list)))
  for(k in 1:length(slopes.list)){sim.slopes[,k] <- slopes.list[[k]]$slope}
  sim.slopes <- data.frame(apply(sim.slopes, 2, function(x) loess(x ~ slopes.list[[1]]$time, span=loess.span)$fitted))
  slope.diffs <- curr.slope$slope - sim.slopes
  slope.qts <- data.frame(time=slopes.list[[1]]$time, 
                           t(data.frame(apply(slope.diffs, 1, function(x) quantile(x, probs=c(0.05,0.5,0.95))))))
  colnames(slope.qts) <- c("time","lower","mean","upper")
  
  # combine the two slopes data sets
  slopes.df <- data.frame(time=curr.slope$time, emp.slope=curr.slope$slope, 
                          slope.diff=slope.qts$mean, slope.lower=slope.qts$lower, 
                          slope.upper=slope.qts$upper)

  # plot the difference in slope of your simulate and empirical disparity metric
  slopes.plot <- ggplot() +
    geom_ribbon(data=slopes.df, aes(x=rev(time), ymin=slope.lower, ymax=slope.upper), fill="#ABDDA4", alpha=0.25) +
    geom_hline(yintercept=0, linetype="dotted") +
    geom_line(data=slopes.df, aes(x=rev(time), y=slope.diff), color="#ABDDA4") +
    xlab("Time (Ma)") + ylab("Difference in Slopes (observed - BM)") + labs(subtitle=trait.name) +
    scale_x_reverse() + theme_classic()
  
  # return the plots and data frames in case they're needed later
  return(list(disparity=accumulation.plot, slopes=slopes.plot, simulated.metric=sim.qts, empirical.metric=emp.vat, slopes.df=slopes.df, sim.vat=sim.vat))
}

# I need the quantiles on the difference in slope between the emp and each simulation

disparity.to.BM.module <- function(phy, scaled.phy, trait.df, sim.num, 
                            metric=c("disparity","variance"), trait.name,
                            loess.span, window.size, trait.correlation=F){
  # get the # of traits
  trait.no <- ncol(trait.df)
  
  # fit BM to get parameter estimates
  if(trait.correlation==F){sim.fit <- mvMORPH::mvBM(phy, trait.df, model="BM1", param=list(constraint="diagonal"), echo=F, diagnostic=F)}
  if(trait.correlation==T){sim.fit <- mvMORPH::mvBM(phy, trait.df, model="BM1", param=list(constraint=FALSE), echo=F, diagnostic=F)}
  # simulate traits
  sim.sim <- lapply(1:sim.num, function(x) mvSIM(tree=phy, nsim=1, model="BM1", param=sim.fit))
  # make sure to get the internal nodes too
  for(j in 1:length(sim.sim)){sim.sim[[j]] <- rbind(sim.sim[[j]], estim(tree=phy, data=sim.sim[[j]], object=sim.fit, asr=T)$estimates)}
  
  # extract the trait at time for each simulated module
  sim.tat <- lapply(sim.sim, function(x) trait.at.time.multi(timeslices=0.1, phy=phy, trait=x, plot=F))

  sim.vat <- NULL
  # loop through each trait and get the VaT
  for(k in 1:(ncol(sim.tat[[1]])-1)){
    curr.tat <- lapply(sim.tat, function(x) x[,c(1,k+1)])
    curr.vat <- lapply(curr.tat, function(x) extract.variance(x, plot=F, metric=metric))
    if(k==1){sim.vat <- lapply(1:sim.num, function(x) curr.vat[[x]][,1:2])}
    if(k>1){
      for(j in 1:sim.num){sim.vat[[j]] <- cbind(sim.vat[[j]], curr.vat[[j]][,2])}
    }
  }
  # get the additive variance
  for (t in 1:length(sim.vat)){sim.vat[[t]]$additive <- apply(sim.vat[[t]], 1, function(x) sum(x[2:length(x)]))}
  # clean it up
  for (t in 1:sim.num){sim.vat[[t]] <- sim.vat[[t]][,c("time","additive")]}
  
  
  # pull all our values into a dataframe
  sim.var <- NULL
  for(t in 1:length(sim.vat)){sim.var <- cbind(sim.var, sim.vat[[t]]$additive)}

  # get the quantiles at each time
  qts <- data.frame(t(apply(sim.var, 1, function(x) quantile(x, probs=c(0.05,0.5,0.95)))))
  # add the time back in
  sim.qts <- data.frame(time=sim.vat[[1]][,1], qts)
  # correct the column names
  colnames(sim.qts) <- c("time","5%","50%","95%")
  
  # combined empirical trait with estimated ancestral nodes from the scaled tree
  emp.fit <- mvBM(tree=egernia.tree, data=trait.df, model="BM1", param=list(constraint=FALSE), echo=F, diagnostic=F)
  emp.trait <- rbind(trait.df, estim(tree=scaled.phy, data=trait.df, object=sim.fit, asr=T)$estimates)
  # extract trait at time object for the empirical data
  emp.tat <- trait.at.time.multi(timeslices=0.1, phy=phy, trait=emp.trait, plot=F)
  # we have to extract the metric for each trait individually
  emp.vat <- data.frame(time=unique(emp.tat$time))
  for (j in 1:(ncol(emp.tat)-1)){
    curr.trait <- emp.tat[,c(1,j+1)]
    temp.vat <- extract.variance(curr.trait, plot=F, metric=metric)
    emp.vat <- cbind(emp.vat, temp.vat$measure)
  }
  # make a new column for the additive measure (additive variance, disparity, et al.)
  emp.vat$additive <- apply(emp.vat, 1, function(x) sum(x[2:length(x)]))
  emp.vat <- data.frame(time=emp.vat$time, measure=emp.vat$additive)

  # plot the simulated and empirical accumulation of your disparity metric
  accumulation.plot <- ggplot() +
    geom_ribbon(data=sim.qts, aes(x=max(time)-time, ymin=`5%`, ymax=`95%`), fill="grey", alpha=0.25) +
    geom_line(data=sim.qts, aes(x=max(time)-time, y=`50%`), color="darkGrey", linetype="dotted") +
    geom_line(data=emp.vat, aes(x=max(time)-time, y=measure), color="black", lwd=1.5) +
    geom_segment(data=sim.qts, aes(x=max(time)*0.9, xend=max(time), 
                                   y=`50%`[nrow(sim.qts)], yend=`50%`[nrow(sim.qts)]), 
                 arrow=arrow(length=unit(0.3,"cm")), color="darkGrey") + 
    geom_segment(data=emp.vat, aes(x=max(time)*0.9, xend=max(time), 
                                   y=measure[nrow(emp.vat)], yend=measure[nrow(emp.vat)]), 
                 arrow=arrow(length=unit(0.3,"cm")), color="black") +
    scale_x_reverse() + xlab("Time (Ma)") + ylab(paste("Disparity", paste0("(",metric,")"))) +
    labs(subtitle=trait.name) + theme_classic()
  
  # now we move on to the difference in slopes
  
  # get the slope through time of our empirical data
  curr.slope <- extract.empirical.slopes(emp.vat, window.size=window.size)
  
  # and the slopes through time of our simulated data
  slopes.list <- extract.empirical.slopes(sim.vat, window.size=window.size)
  sim.slopes <- data.frame(matrix(nrow=nrow(slopes.list[[1]]),ncol=length(slopes.list)))
  for(k in 1:length(slopes.list)){sim.slopes[,k] <- slopes.list[[k]]$slope}
  sim.slopes <- data.frame(apply(sim.slopes, 2, function(x) loess(x ~ slopes.list[[1]]$time, span=loess.span)$fitted))
  slope.diffs <- curr.slope$slope - sim.slopes
  slope.qts <- data.frame(time=slopes.list[[1]]$time, 
                          t(data.frame(apply(slope.diffs, 1, function(x) quantile(x, probs=c(0.05,0.5,0.95))))))
  colnames(slope.qts) <- c("time","lower","mean","upper")
  
  # combine the two slopes data sets
  slopes.df <- data.frame(time=curr.slope$time, emp.slope=curr.slope$slope, 
                          slope.diff=slope.qts$mean, slope.lower=slope.qts$lower, 
                          slope.upper=slope.qts$upper)
  
  # plot the difference in slope of your simulate and empirical disparity metric
  slopes.plot <- ggplot() +
    geom_ribbon(data=slopes.df, aes(x=rev(time), ymin=slope.lower, ymax=slope.upper), fill="#ABDDA4", alpha=0.25) +
    geom_hline(yintercept=0, linetype="dotted") +
    geom_line(data=slopes.df, aes(x=rev(time), y=slope.diff), color="#ABDDA4") +
    xlab("Time (Ma)") + ylab("Difference in Slopes (observed - BM)") + labs(subtitle=trait.name) +
    scale_x_reverse() + theme_classic()
  
  # return the plots and data frames in case they're needed later
  return(list(disparity=accumulation.plot, slopes=slopes.plot, simulated.metric=sim.qts, empirical.metric=emp.vat, slopes.df=slopes.df, sim.vat=sim.vat))
}


# function to plot the diff slopes 
# this function doesn't work properly
diff.slopes.module <- function(disparity.df, window.size, trait.name="trait"){
  
  # set the window size
  window <- round(window.size/2,0)
  
  # get the slope through time of our empirical data
  curr.slope <- extract.empirical.slopes(disparity.df, window.size=window)
  
  # put the simulated values into their own dataframe
  sim.list <- list(`5%`=data.frame(time=disparity.df$time,measure=disparity.df$`5%`),
                   `50%`=data.frame(time=disparity.df$time,measure=disparity.df$`50%`),
                   `95%`=data.frame(time=disparity.df$time,measure=disparity.df$`95%`))
  #sim.list <- list(`5%` =data.frame(time=disparity.df$time, measure=loess(`5%`~time,data=disparity.df, span=1)$fitted),
  #                 `50%`=data.frame(time=disparity.df$time, measure=loess(`50%`~time,data=disparity.df,span=1)$fitted),
  #                 `95%`=data.frame(time=disparity.df$time, measure=loess(`95%`~time,data=disparity.df,span=1)$fitted))
  

  #ggplot() +
  #  geom_line(data=sim.list[[1]],aes(x=time,y=measure),col="lightBlue") +
  #  geom_line(data=sim.list2[[1]],aes(x=time,y=measure),col="darkBlue") +
  #  geom_line(data=sim.list[[3]],aes(x=time,y=measure),col="lightGreen") +
  #  geom_line(data=sim.list2[[3]],aes(x=time,y=measure),col="darkGreen") +
  #  geom_line(data=sim.list[[2]], aes(x=time,y=measure),col="yellow") +
  #  geom_line(data=curr.slope, aes(x=time,y=measure), col="red") +
  #  theme_classic()
    
  
  # and the slopes through time of our simulated data
  slopes.list <- extract.empirical.slopes(sim.list, window.size=window)
  sim.slopes <- NULL
  for(k in 1:length(slopes.list)){sim.slopes <- cbind(sim.slopes, slopes.list[[k]]$slope)}
  sim.slopes <- data.frame(time=slopes.list[[1]]$time, sim.slopes)
  colnames(sim.slopes) <- c("time","5%","50%","95%")
  
  # combine the two slopes data sets
  slopes.df <- data.frame(time=curr.slope$time, emp.slope=curr.slope$slope, sim.slope=sim.slopes$`50%`,
                          sim.lower=sim.slopes$`5%`, sim.upper=sim.slopes$`95%`,
                          slope.diff=curr.slope$slope-sim.slopes$`50%`,
                          slope.diff.lower=curr.slope$slope-sim.slopes$`5%`,
                          slope.diff.upper=curr.slope$slope-sim.slopes$`95%`)
  
  # plot the difference in slope of your simulate and empirical disparity metric
  slopes.plot <- ggplot() +
    #geom_ribbon(data=slopes.df, aes(x=max(time)-time, ymin=slope.diff, ymax=slope.diff), fill="#ABDDA4", alpha=0.25) +
    geom_hline(yintercept=0, linetype="dotted") +
    geom_line(data=slopes.df, aes(x=max(time)-time, y=slope.diff), color="#ABDDA4", lwd=1.5) +
    xlab("Time (Ma)") + ylab("Difference in Slopes (observed - BM)") + labs(subtitle=trait.name) +
    scale_x_reverse() + theme_classic()
  
  return(list(slopes.plot=slopes.plot, slopes.df=slopes.df))
}

# this function estimates diff in slopes between observed and simulated data
# 'diff.slopes.module' does not work properly
plot.module.slope <- function(observed.df, simulated.list, loess.span, window.size){
  # estimate the slopes of the observed data
  trait.obs <- extract.empirical.slopes(observed.df, window.size=window.size)
  
  # estimate the slopes of the simulated data (all of them)
  sim.slope  <- lapply(simulated.list, function(x) extract.empirical.slopes(x, window.size=window.size))
  # make an empty dataframe to hold the slopes across all simulations
  sim.slopes <- data.frame(matrix(nrow=nrow(sim.slope[[1]]),ncol=length(sim.slope)))
  # fill the dataframe with the slope estimates
  for(j in 1:length(sim.slope)){sim.slopes[,j] <- sim.slope[[j]]$slope}
  # loess smooth the data moderately
  sim.slopes <- data.frame(apply(sim.slopes, 2, function(x) loess(x ~ sim.slope[[1]]$time, span=loess.span)$fitted))
  
  # get the difference between observed and simulated slopes (all)
  slope.diffs <- trait.obs$slope - sim.slopes
  # get the quantiles on the diffs
  slope.qts <- data.frame(time=trait.obs$time,
                          emp=trait.obs$slope,
                          t(apply(slope.diffs, 1, function(x) quantile(x, probs=c(0.05,0.5,0.95)))))
  # rename the columns
  colnames(slope.qts) <- c("time","empirical","lower","mean","upper")
  
  # plot the trend
  ggplot() +
    geom_hline(yintercept=0, linetype="dotted") +
    geom_ribbon(data=slope.qts, aes(x=rev(time),ymin=lower,ymax=upper), fill="#ABDDA4", alpha=0.25) +
    geom_line(data=slope.qts, aes(x=rev(time),y=mean), color="#ABDDA4") +
    xlab("Time (Ma)") + ylab("Difference in Slope (observed - BM") +
    scale_x_reverse() + theme_classic()
}



#Make a function to extract the slopes of empirical disparity trends (variance, functional diversity, et al.)
extract.empirical.slopes <- function(disparity.dataframe, window.size){
  window <- round(window.size/2,0)
  
  if(class(disparity.dataframe)=="data.frame"){
    curr.trait <- disparity.dataframe
    curr.slope <- NULL
    
    for(i in 1:(nrow(curr.trait))){
      min.t <- i-window; max.t <- i+window
      if(min.t<0){min.t=0}; if(max.t>nrow(curr.trait)){max.t=nrow(curr.trait)}
      temp.trait <- curr.trait[(min.t):(max.t),]
      temp.trait <- temp.trait[complete.cases(temp.trait),]
      curr.slope <- append(curr.slope, lm(temp.trait[,2] ~ temp.trait$time)$coefficients[2])
    }
    emp.slopes <- cbind(curr.trait, slope = curr.slope)
  }
  
  if(class(disparity.dataframe)=="list"){
      emp.slopes <- list()
      for (k in 1:length(disparity.dataframe)){
        curr.trait <- disparity.dataframe[[k]]

        curr.slope <- NULL
        for(i in 1:(nrow(curr.trait))){
          min.t <- i-window; max.t <- i+window
          if(min.t<0){min.t=0}; if(max.t>nrow(curr.trait)){max.t=nrow(curr.trait)}
          temp.trait <- curr.trait[(min.t):(max.t),]
          temp.trait <- temp.trait[complete.cases(temp.trait),]
          curr.slope <- append(curr.slope, lm(temp.trait[,2] ~ temp.trait$time)$coefficients[2])
        }
        names(curr.slope) <- NULL
        curr.emp <- cbind(curr.trait, slope = curr.slope)
        emp.slopes[[k]] <- curr.emp
      }
    names(emp.slopes) <- names(disparity.dataframe)
  }
return(emp.slopes)
}
# e.g. extract.empirical.slopes(empirical.variances) or extract.empirical.slopes(module.disparity)



