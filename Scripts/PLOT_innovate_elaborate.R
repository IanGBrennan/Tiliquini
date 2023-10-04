library(dispRity)
library(mcmcmcglmmm)
library(mvMORPH)
library(ggplot2)
library(scico)
library(RColorBrewer)
# https://builtin.com/data-science/step-step-explanation-principal-component-analysis

# Create a function to measure the elaboration and innovation in a given dataset

# this first function is basic and doesn't include phylogenetic information
# innovate = the residuals from the regression (major axis of variation)
# elaborate = the distance along the regression from the geometric mean 
innovate.elaborate <- function(trait.df, plot=F, summary=F){
  require(ggplot2)
  # run PCA on the trait data, we will work with the first two PCs regardless
  pca <- prcomp(trait.df)
  # run a linear model on the first two PC axes
  pca.lm <- lm(pca$x[,2] ~ pca$x[,1])
  # get the geometric mean of the trait data (should lie on the regression line)
  geo.mean <- c(mean(pca$x[,1]), mean(pca$x[,2]))
  # adjust the data to where they intersect the regression
  pca.elab <- data.frame(PC1 = pca$x[,1],
                         PC2 = pca$x[,2] - pca.lm$residuals)
  # get the euclidean distance between data and the geo.mean
  elab.dist <- apply(pca.elab, 1, function(x) euclidean(x, geo.mean))
  # make the dataframe
  output <- data.frame(elaborate = elab.dist, innovate = abs(pca.lm$residuals))
  # plot the results if you'd like
  if(plot==T){print(autoplot(pca, label=T, loadings=T, loadings.label=T) + 
                      geom_smooth(method="lm",se=F) + theme_classic())}
  # show the PC variance scores
  if(summary==T){print(summary(pca))}
  # return the data
  return(output)
}

#trait.df <- data.frame(anc.all.traits)


# this second function incorporates the phylogeny and measures:
# elaboration = the distance along the regression (major axis of variation) from the modelled root value
# innovation = the residual distance from the regression (orthogonal to the major axis)
innovate.elaborate.phy <- function(trait.df, phy, plot=F, summary=F, PCs=c(1,2),
                                   angles=c("equal","conservative","liberal")){
  require(ggplot2)
  require(phytools)
  require(ggfortify)
  require(patchwork)
  # reorder the trait dataframe to match the order of tips on the tree
  trait.df <- trait.df[order(match(rownames(trait.df), phy$tip.label)), , drop=F]
  # run PCA on the trait data, we will work with the first two PCs regardless
  pca <- prcomp(trait.df)
  # run a linear model on the first two PC axes
  pca.lm <- lm(pca$x[,PCs[[2]]] ~ pca$x[,PCs[[1]]])
  # re-center the regression through the root
  root <- pca$x[Ntip(phy)+1,] # get position of the root on PCs 1/2
  pca.lm$coefficients[[1]] <- root[[2]] # change the coefficient of the linear model
  geo.mean <- c(mean(pca$x[,PCs[[1]]]), mean(pca$x[,PCs[[2]]])) # what was the previous center of the data?
  root.mean <- data.frame(PC1=root[1], PC2=root[2]) # set the new data center as the root values
  diff.mean <- geo.mean - root.mean # what's the difference between the geometric (old) and root (new) centers?
  # correct the residuals
  new.resid <- pca.lm$residuals + diff.mean$PC2
  # adjust the data to where they intersect the regression
  pca.elab <- data.frame(PC1 = pca$x[,PCs[[1]]],
                         PC2 = pca$x[,PCs[[2]]] - new.resid)
  # get the origin of the trait data, represented by the root  (should lie on the regression line)
  new.geo.mean <- pca.elab[Ntip(phy)+1,]
  # get the euclidean distance between data and the geo.mean
  elab.dist <- apply(pca.elab, 1, function(x) euclidean(x, new.geo.mean))
  # make the dataframe
  output <- data.frame(elaborate.root = elab.dist, innovate.root = abs(new.resid))
  # get parent nodes for each node/tip
  parent.df <- data.frame(child = 1:(Ntip(phy) + Nnode(phy)),
                          parent = sapply(1:(Ntip(phy) + Nnode(phy)), 
                                          function(x) getParent2(phy,x)))
  # get elaboration and innovation scores from parent to child
  pca.parent <- cbind(pca$x[,c(PCs[[1]],PCs[[2]])], parent.df)
  output$elaborate.increment <- abs(sapply(pca.parent$child, function(x) pca.parent[x,1] - 
                                             pca.parent[pca.parent[x,"parent"],1]))
  output$innovate.increment <- abs(sapply(pca.parent$child, function(x) pca.parent[x,2] - 
                                            pca.parent[pca.parent[x,"parent"],2]))
  # get the angle of change from point to point
  from <- pca.parent[,1:2]
  to <- data.frame(t(sapply(pca.parent$parent, function(x) unlist(pca.parent[x,1:2]))))
  diff.df <- to - from
  diff.angle <- data.frame(useful::cart2pol(diff.df[,1], diff.df[,2], degrees=T))
  output$angle <- diff.angle$theta
  output$distance <- diff.angle$r
  
  # decide how you want to define the angles of "innovation" vs. "elaboration"
  # "equal" option allows allots 180* to each, in 90/90/90/90 increments
  if(angles=="equal"){
    output$change <- unlist(sapply(output$angle, function(x) if(x >= 0 && x <= 45) paste("elaboration")
                                   else if(x > 45 && x < 135) paste("innovation")
                                   else if(x >= 135 && x <= 225) paste("elaboration")
                                   else if(x > 225 && x < 315) paste("innovation")
                                   else if(x >= 315 && x <= 360) paste("elaboration")))
  }
  # "liberal" option favors "innovation" (270*) over "elaboration (90*) in 45/135/45/135 increments
  if(angles=="liberal"){
    output$change <- unlist(sapply(output$angle, function(x) if(x >= 0 && x <= 22.5) paste("elaboration")
                                   else if(x > 22.5 && x < 157.5) paste("innovation")
                                   else if(x >= 157.5 && x <= 202.5) paste("elaboration")
                                   else if(x > 202.5 && x < 337.5) paste("innovation")
                                   else if(x >= 337.5 && x <= 360) paste("elaboration")))
  }
  if(angles=="conservative"){
    output$change <- unlist(sapply(output$angle, function(x) 
                                   if(x >= 0 && x <= 30) paste("elaboration")
                                   else if(x > 30 && x < 60) paste("neither")
                                   else if(x >= 60 && x <= 120) paste("innovation")
                                   else if(x > 120 && x < 150) paste("neither")
                                   else if(x >= 150 && x <= 210) paste("elaboration")
                                   else if(x > 210 && x < 240) paste("neither")
                                   else if(x >= 240 && x<= 300) paste("innovation")
                                   else if(x > 300 && x < 330) paste("neither")
                                   else if(x >= 330 && x <= 360) paste("elaboration")))
  }
  # add the parent node info back in
  output$parent <- pca.parent$parent
  output$child <- pca.parent$child
  # color the edges according to the prevailing direction of change (innovation or elaboration)
  edges <- unlist(sapply(pca.parent$child, function(x) which(phy$edge[,2]==x)))
  output$edges <- append(edges, 0, after=Ntip(phy))
  output$edge.color <- sapply(output$change, function(x) if(x=="elaboration"){paste("lightblue")}
                                                         else if(x=="innovation"){paste("orangered")}
                                                         else if(x=="neither"){paste("lightgrey")})
  #color.df <- output[order(output$edges),]
  #color.df <- color.df[-1,]
  
  # scale and color the amount of change
  output$distance.scaled <- round((output$distance - min(output$distance))/diff(range(output$distance)) * 99) + 1
  # split the dataframes, we'll bring them back together at the end
  out.inn <- dplyr::filter(output, change=="innovation")
  out.ela <- dplyr::filter(output, change=="elaboration")
  if("neither" %in% output$change){out.neither <- dplyr::filter(output, change=="neither")}
  # make appropriate color ramps
  colors.inn <- (colorRampPalette(brewer.pal(9, "Oranges")[2:9])(100))
  colors.ela <- (colorRampPalette(brewer.pal(9, "Blues")[2:9])(100))
  # get scaled colors from our ramps
  out.inn$color.scaled <- colors.inn[out.inn$distance.scaled]
  out.ela$color.scaled <- colors.ela[out.ela$distance.scaled]
  if("neither" %in% output$change){out.neither$color.scaled <- "lightgrey"}
  # combine the dataframe back together
  if("neither" %in% output$change){output <- rbind(out.inn, out.ela, out.neither)}
  else{output <- rbind(out.inn, out.ela)}
  # extract the colors
  color.df <- output[order(output$edges),]
  color.df <- color.df[-1,]
  # plot if you'd like
  if(plot==T){
    layout(
      matrix(c(1,2,1,3), nrow=2, ncol=2, byrow=TRUE), 
      widths=c(5,1), 
      heights=c(1,1)
    )
    plot(phy, edge.col=color.df$color.scaled, cex=0.3, edge.width=2); axisPhylo()
    color.bar(colors.inn, min=0, max=round(max(output$distance),2), title="innovation")
    color.bar(colors.ela, min=0, max=round(max(output$distance),2), title="elaboration")
    }
  
  # plot the results if you'd like
  pca.x <- data.frame(pca$x)
  figure <- (ggplot(data=pca.x, aes(x=pca.x[,PCs[[1]]],y=pca.x[,PCs[[2]]],label=rownames(pca$x))) +
      geom_point() + geom_text(hjust=0, vjust=0, color="lightGrey") +
      geom_abline(intercept=pca.lm$coefficients[[1]], 
                                 slope=pca.lm$coefficients[[2]], color="orangered", lty="dotted") +
      geom_point(aes(x=new.geo.mean$PC1,y=new.geo.mean$PC2), fill="lightblue", size=3, pch=21) + theme_classic() + 
      autoplot(pca, data=trait.df, colour = "white", scale=F, loadings=T, loadings.label=T, loadings.colour="orangered", loadings.label.colour="black") + theme_classic())
  if(plot==T){print(figure)}
  # show the PC variance scores
  if(summary==T){print(summary(pca))}
  # return the data
  #return(output[order(output$innovate.increment, decreasing=T),])
  return(list(ie.df = output, pca = pca, lm = pca.lm, plot = figure))
}

# helper function: adjust the getParent function of phytools
getParent2 <- function (tree, node) 
{
  ind <- which(tree$edge[, 2] == node)
  if (length(ind) > 0) 
    pp <- tree$edge[which(tree$edge[, 2] == node), 1]
  else {
    pp <- node
  }
  pp
}

# make a function to visualize the contribution of traits to the PCs
plot.rotations <- function(pca.obj, plot=T){
  require(reshape2)
  require(ggplot2)
  # make the pca rotation matrix a dataframe
  pca.rot <- data.frame(pca.obj$rotation)
  # change the column names to reflect the relative contributions to the variance  
  colnames(pca.rot) <- sapply(1:ncol(pca.rot), function(x) paste0(colnames(pca.rot)[[x]],"_",round(summary(pca.obj)$importance[2,x]*100,2)))
  # create a variable with the rownames
  pca.rot$trait <- rownames(pca.rot)
  # reshape the data from wide to long format
  pca.rot <- reshape2::melt(pca.rot, id.vars="trait", variable.name="PC")
  # plot the resulting data
  rot.plot <- ggplot(pca.rot,aes(x=trait,y=value)) +
    geom_col() + theme_bw() +
    facet_wrap(~PC) + coord_flip() + scale_x_discrete(limits=rev)
  if(plot==T){print(rot.plot)}
  return(rot.plot)
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

# function to get the euclidean distance between two points
euclidean <- function(a, b) sqrt(sum((a - b)^2))


# this function will plot the regression of a linear model based on PC1 and PC2
# through time from the root of a tree to the tips
inn.elab.tt <- function(phy, trait.df, pca.in, base.color=c("Blues","Greens","Reds","Oranges","Purples","Spectral","YlOrRd")){
  require(ape); require(RColorBrewer); require(phytools)
  #curr.pca <- prcomp(trait.df)
  
  ages.df <- round(data.frame(nodeHeights(phy)),3)
  ages.df$start <- round(max(ages.df) - ages.df[,1],3)
  ages.df$stop  <- round(max(ages.df) - ages.df[,2],3)
  ages.df$parent <- phy$edge[,1]
  ages.df$child  <- phy$edge[,2]
  btimes <- branching.times(phy)
  btimes.reorder <- ages.df[with(ages.df, order(ages.df$start, ages.df$stop, decreasing=T)),]
  btimes.vec <- sort(unique(c(ages.df$start, ages.df$stop)), decreasing=T)
  tree.segments <- data.frame(embed(btimes.vec, 2)[,2:1])
  colnames(tree.segments) <- c("start", "stop")
  tree.segments$mid.point <- apply(tree.segments, 1, mean)
  
  coeffs <- NULL
  for (k in 1:nrow(tree.segments)){
    curr.segments <- dplyr::filter(ages.df, start >= tree.segments$mid[[k]] &
                                     stop  <= tree.segments$mid[[k]])
    curr.data <- data.frame(pca.in[curr.segments$child,c("PC1","PC2")])
    curr.lm <- lm(curr.data$PC2 ~ curr.data$PC1)
    coeffs <- rbind(coeffs, data.frame(a=curr.lm$coefficients[[1]], b=curr.lm$coefficients[[2]]))
  }
  coeffs$color <- colorRampPalette(brewer.pal(9, base.color)[2:9])(nrow(coeffs))
  plot(1, type="n", xlab="", ylab="", xlim=c(-5, 5), ylim=c(-5, 5))
  for(j in 1:nrow(coeffs)){
    abline(a=coeffs[j,"a"], b=coeffs[j,"b"], col=coeffs$color[[j]], lwd=3)
    profvis::pause(0.1)
  }
  
}
