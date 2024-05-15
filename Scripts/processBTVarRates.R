processBT.log <- function(log.file, single.trait=T, skip.lines=43, model=c("VarRates","fabric")){
  # establish how many lines to skip
  if(single.trait==T){
      if(model=="VarRates"){skip.lines=43}
      if(model=="fabric"){skip.lines=45}
  }
  if(single.trait==F){skip.lines <- skip.lines}

  # read in the log file, skipping all the way down to the actual mcmc
  bt.log <- read.delim(log.file, skip=skip.lines)
  # remove the empty column at the end
  bt.log <- bt.log[,1:(ncol(bt.log)-1)]
  #
  ESS <- apply(bt.log, 2, coda::effectiveSize)[4:ncol(bt.log)]
  var.mean <- apply(bt.log, 2, mean)[4:ncol(bt.log)]
  if(single.trait==T){sigma.qts <- quantile(bt.log$Sigma.2.1, probs = c(0.001,0.05,0.5,0.95,0.999))}
  if(single.trait==F){
    sigmas <- dplyr::select(bt.log, starts_with("Sigma"))
    sigma.qts <- apply(sigmas, 2, function(x) quantile(x, probs=c(0.001,0.05,0.5,0.95,0.999)))
  }
  
  return(list(ESS=ESS, Mean=var.mean, Sigma_Quantiles=sigma.qts, Log=bt.log))
}