require(phytools)
require(dplyr)
require(tidyr)

# plot a phylogeny with branches colored by rate
# coloring options are a bit tricky, so stick to 'log.rates=F relative.rates=T'
plot.MCMCRates.tree <- function(mcmc.sum, phy, col.palette = "Blues", legend = F,
                               tree.type = c("phylogram", "fan"),
                               log.rates = F, relative.rates = F,
                               trait=NULL, outline=F,
                               pos.selection = F, shift.rates = F,
                               annotations=F){
  
  #ntraits <- length(RR)
  #opt.layout <- n2mfrow(ntraits)
  #if(!(ntraits %% 2) == 0){mat.len <- ntraits + 1} else{mat.len <- ntraits}
  
  #all.node.rates <- NULL
  
  # set the plot layout
  #layout(matrix(nrow = opt.layout[1], ncol = opt.layout[2], 1:mat.len))
  
  # Set the layout
  if(legend==T){layout(
    matrix(c(1,2,1,3), ncol=2, byrow=T), 
    widths=c(4,1), 
    heights=c(1,1))
  }
  
  # extract the rates and edge/node information
  #if(log.rates==T){BT$Mean.SigV <- log(BT$Mean.SigV)}
  all.edges <- data.frame(edge.rate = mcmc.sum$mean,
                          edge = mcmc.sum$edge,
                          child.node = mcmc.sum$child.node)
  if(log.rates==T && relative.rates==F){all.edges$raw.rate <- all.edges$edge.rate; all.edges$edge.rate <- log(all.edges$edge.rate)}
  #if(log.rates==T && relative.rates==T){all.edges$edge.rate <- mean(all.edges$edge.rate) / all.edges$edge.rate}
  if(log.rates==T && relative.rates==T){all.edges$raw.rate <- all.edges$edge.rate / mean(all.edges$edge.rate); all.edges$edge.rate <- log(all.edges$edge.rate)}
  if(log.rates==F && relative.rates==T){all.edges$edge.rate <- all.edges$edge.rate / mean(all.edges$edge.rate)}
  #all.edges <- all.edges[-which(all.edges$child.node==Ntip(phy)+1),]
  
  # scale the rates by the fastest rate (for ease of plotting)
  all.edges$rounded.rates <- round((all.edges$edge.rate - min(all.edges$edge.rate))/diff(range(all.edges$edge.rate)) * 99) + 1
  
  # sort it by the original order!
  all.edges <- all.edges[order(match(all.edges$child.node, phy$edge[,2])),]  
  
  # make a color ramp that suits the data
  if(col.palette %in% c("magma", "inferno", "plasma", "viridis","cividis", "rocket", "mako")){
    new.cols <- rev(viridis(n=100, option=col.palette))
    #if(col.palette %in% c("mako", "plasma", "inferno")){new.cols <- rev(new.cols)}
  }else if(col.palette == "GyYlRd"){
    colfunc <- colorRampPalette(c("#c7c7c7", "#ffee2e", "#e81d13"))
    new.cols <- colfunc(100)
  }else if(col.palette == "GnBu"){
    colfunc <- colorRampPalette(c("#C7E9B4","#7FCDBB","#41B6C4","#2C7FB8","#253494","#081D58"))
    new.cols <- colfunc(100)
  }else{
    col.ramp <- colorRampPalette(brewer.pal(9, col.palette)[2:9])
    new.cols <- (col.ramp(100))
    if(col.palette %in% c("Spectral","RdYlGn","RdYlBu","RdGy","RdBu","PuOr","PRGn","PiYG","BrBG")){new.cols <- rev(new.cols)}
  }
  
  
  #all.edges$edge.color <- sapply(all.edges$rounded.rates, function(x) new.cols[x]) # this is definitely wrong
  all.edges$edge.color <- new.cols[all.edges$rounded.rates]
  
  # plot the phylogeny with colored branches
  if(outline==T){plot.phylo(phy, edge.width=5, type=tree.type, cex=0.3); par(new=T)}
  plot.phylo(phy, edge.color = unlist(all.edges$edge.color), edge.width=3, type=tree.type, cex=0.3, open.angle=5)
  axisPhylo()
  title(trait)
  
  # plot the location of edges undergoing putative shifts in evolutionary rate
  if(shift.rates==T){
    sr <- dplyr::filter(BT, Pct.time.scaled >= 50 & Mean.Scalar >= 2)
    if(nrow(sr) > 0){
      for (j in 1:nrow(sr)){
        edgelabels(text="", edge=noquote(sr$edge[j]), frame="circle", bg="white", cex=0.4, col=0.5)
      }      
    }
    if(annotations==T){title(xlab = paste(nrow(sr), "instances of rate shifts"))}
  }
  
  # plot the location of edges undergoing shifts in evolutionary rate
  if(shift.rates==T){
    sr <- dplyr::filter(BT, Pct.time.scaled >= 70 & Pct.time.scaled < 95 & Mean.Scalar >= 2)
    if(nrow(sr) > 0){
      for (j in 1:nrow(sr)){
        edgelabels(text="", edge=noquote(sr$edge[j]), frame="circle", bg="grey", cex=0.3, col=0.5)
      }      
    }
    #if(annotations==T){title(xlab = paste(nrow(sr), "instances of rate shifts"))}
  }
  
  # plot the location of edges undergoing positive selection sensu Baker et al. 2016
  if(pos.selection==T){
    #g2 <- dplyr::filter(BT, Mean.deltaVB >= 2 & Pct.time.scaled >= 95)
    g2 <- dplyr::filter(BT, Mean.Scalar >= 2 & Pct.time.scaled >= 95)
    
    if(nrow(g2) > 0){
      for (j in 1:nrow(g2)){
        edgelabels(text="", edge=noquote(g2$edge[j]), frame="circle", bg="black", cex=0.2, col=0.5)
      }   
    }
    if(annotations==T){title(sub = paste(nrow(g2), "instances of positive selection"))}
  }
  
  # if we chose to plot a legend, do that now
  if(legend==T && log.rates==F && relative.rates==T){plot(0,type='n',axes=FALSE,ann=FALSE);
    color.bar(new.cols, 
              min=round(min(all.edges$edge.rate),2),
              max=round(max(all.edges$edge.rate),2))}
  if(legend==T && log.rates==T && relative.rates==F){plot(0,type='n',axes=FALSE,ann=FALSE);
    color.bar(new.cols, 
              min=round(min(all.edges$raw.rate),5),
              max=round(max(all.edges$raw.rate),5))}
  if(legend==T && log.rates==T && relative.rates==T){plot(0,type='n',axes=FALSE,ann=FALSE);
    color.bar(new.cols, 
              min=round(min(log(all.edges$raw.rate)),2),
              max=round(max(log(all.edges$raw.rate)),2))}
  if(legend==T && log.rates==F && relative.rates==F){plot(0,type='n',axes=FALSE,ann=FALSE);
    color.bar(new.cols, 
              min=round(min(all.edges$edge.rate),2),
              max=round(max(all.edges$edge.rate),2))}
}


# this function is just to plot the scale bar, which is dumb, but it works, so hey.
color.bar <- function(lut, min=0, max=100, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  
  #dev.new(width=1.75, height=5)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, ticks, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}

# this is the helper function for processing the MCMC file to extract and summarize branch rate information
process_MCMCrates <- function(mcmc.df) {
  
  # provide adequate warning 
  message(paste("processing may be slow provided", nrow(mcmc.df), "MCMC samples"))
  
  # extract just the rate parameters
  rates <- dplyr::select(mcmc.df, starts_with("r"))
  
  # reshape the data into long format
  full.rates <- rates %>% 
    tidyr::pivot_longer(
      cols = 1:ncol(rates), 
      values_to = "rate",
      names_to = "name"
    )
  
  # add some extra data columns
  full.rates$partition <- sapply(full.rates$name, function(x) strsplit(x, "_")[[1]][2])
  full.rates$nedge <- sapply(full.rates$name, function(x) strsplit(x, "_")[[1]][3])
  full.rates$edge <- sapply(full.rates$nedge, function(x) strsplit(x,"n")[[1]][2])
  
  # summarize mean rates per edge
  mean.rates <- full.rates %>%
    group_by(edge) %>%
    summarise(mean = mean(rate), n = n())
  # translate edge to numeric
  mean.rates$edge <- as.numeric(as.character(mean.rates$edge))
  #mean.rates$edge <- mean.rates$node + Ntip(tree)
  
  # establish the relationship between edge and node numbers
  tree.edge <- data.frame(parent.node = tree$edge[,1],
                          child.node = tree$edge[,2],
                          ape.edge = 1:nrow(tree$edge),
                          mcmc.edge = tree$edge[,2])
  
  # connect the mcmc edge convention with ape
  colnames(mean.rates)[[1]] <- "mcmc.edge"
  mean.rates$edge <- sapply(mean.rates$mcmc.edge, function(x) tree.edge[which(tree.edge$mcmc.edge == x),"ape.edge"])
  mean.rates$child.node <- mean.rates$mcmc.edge
  
  # return the processed dataframe
  return(mean.rates)
}