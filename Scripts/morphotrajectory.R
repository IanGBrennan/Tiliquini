require(ggplot2)
require(RColorBrewer)
require(dplyr)
require(tibble)
require(phytools)
require(ggdensity)
require(mvMORPH)

# the 'morphotrajectory.RR' and 'morphotrajectory.VR' functions take empirical traits and plot 
# a phylomorphospace atop a simulated (BM) distribution of traits

morphotrajectory.RR <- function(phy, trait, anc, branch.rates, tip.spread, 
                             focus=c("tip", "clade"), sim.trait=NULL,  
                             brew.pal="YlOrRd", psize=3, lsize=2){
  # check that the tip names are correct
  if(!all(tip.spread %in% phy$tip.label)){stop("WARNING: check tip names")}
  # check that all the tips are in the trait dataframe
  if(nrow(trait) != length(phy$tip))
    stop("your trait dataframe must have the same number of rows as species in your tree")
  
  # make sure your traits are in the same order as your tree
  trait <- trait %>%
    rownames_to_column('taxon') %>%
    arrange(match(taxon, phy$tip.label)) %>%
    column_to_rownames('taxon')
  
  # create a phylomorphospace object using phytools
  obj <- phylomorphospace(tree=phy, X=trait, A=anc, 
                          ftype="off", bty="n", node.size=c(0,1))
  # extract the xy coordinates of the nodes and segments
  phydat <- data.frame(xstart=obj$xx[obj$edge[,1]],
                       ystart=obj$yy[obj$edge[,1]],
                       xstop=obj$xx[obj$edge[,2]],
                       ystop=obj$yy[obj$edge[,2]],
                       nodestart=obj$edge[,1],
                       nodestop=obj$edge[,2])
  
  # log and scale the module rates
  rate.df <- data.frame(raw.rate = branch.rates[,1])
  rate.df$rate <- log(rate.df$raw.rate)
  rate.df$scaled <- round((rate.df$rate - min(rate.df$rate))/diff(range(rate.df$rate)) * 99) + 1
  # choose colors and apply them to the scaled rates
  colors <- (colorRampPalette(brewer.pal(9, brew.pal)[2:9])(100))
  if(brew.pal %in% c("Spectral","RdYlGn","RdYlBu","RdGy","RdBu","PuOr","PRGn","PiYG","BrBG")){colors <- rev(colors)}

  rownames(rate.df)[Ntip(phy):nrow(rate.df)] <- (1:length(phy$tip.label))
  rate.df$nodestart <- as.numeric(as.character(rownames(rate.df)))
  phydat <- left_join(phydat, rate.df, by="nodestart")
  phydat$color <- colors[phydat$scaled]
  
  # if you're working with a clade
  if(focus == "clade"){
    # determine the MRCA node of your target tips
    mrcn <- getMRCA(phy, tip = tip.spread)
    # get all the descendant nodes of your MRCA including internals
    all.desc <- getDescendants(phy, mrcn)    
  }
  # if you're working with a single tip
  if(focus == "tip"){all.desc <- which(phy$tip.label==tip.spread[[1]])}

  # identify which descendant nodes are tips
  tip.desc <- all.desc[which(all.desc <= Ntip(phy))]
  # get the node path from the root to each tip, and make a vector of the unique nodes
  focal.nodes <- unique(unlist(sapply(tip.desc, function(x) nodepath(phy, from=Ntip(phy)+1, to=x))))
  # subset the phydat dataframe to just the nodes of interest
  focal.dat <- dplyr::filter(phydat, nodestart %in% focal.nodes & nodestop %in% focal.nodes)
  
  # plot the branches from root to tips colored by evolutionary rate, with simulated data
  ggplot()+
    {if(!is.null(sim.trait))geom_hdr(data=sim.trait, aes(x=sim.trait[,1], y=sim.trait[,2]), method="mvnorm", fill="lightGrey", probs=c(0.25,0.5,0.8,0.95))} +
    #geom_hdr_lines(data=sim.handfoot, aes(x=Hand, y=Foot), method="mvnorm") +
    geom_segment(data=phydat, aes(x=xstart,y=ystart,xend=xstop,yend=ystop), lwd=lsize-1, color="grey") +
    geom_point(data=trait, aes(x=trait[,1], y=trait[,2]), size=psize, shape=21, fill="lightGrey", color="#757574") +
    geom_point(data=focal.dat[1,], aes(x=xstart,y=ystart), size=psize+1, shape=3, color="black") +
    geom_point(data=focal.dat[1,], aes(x=xstart,y=ystart), size=psize+1, shape=4, color="black") +
    #geom_point(data=anc[which(rownames(anc)==mrcn),], aes(x=anc[,1], y=anc[,2]), size=psize+1, shape=21, fill="white", color="black") +
    geom_segment(data=focal.dat,aes(x=xstart,y=ystart,xend=xstop,yend=ystop),color=focal.dat$color, lwd = lsize) +
    geom_point(data=focal.dat, aes(x=xstart, y=ystart), size=psize, shape=21, fill="white", color="black") +
    geom_point(data=focal.dat, aes(x=xstop, y=ystop),   size=psize, shape=21, fill="white", color="black") +
    labs(x=colnames(trait)[[1]], y=colnames(trait)[[2]]) +
    theme_classic()
  
}

morphotrajectory.VR <- function(phy, scaled.phy, trait, branch.rates, tip.spread, 
                             focus=c("tip", "clade"), sim.num, 
                             trait.correlation=c("no","yes","both"),
                             brew.pal="YlOrRd", psize=3, lsize=2){
  # check that the tip names are correct
  if(!all(tip.spread %in% phy$tip.label)){stop("WARNING: check tip names")}
  # check that all the tips are in the trait dataframe
  if(nrow(trait) != length(phy$tip))
    stop("your trait dataframe must have the same number of rows as species in your tree")
  
  # make sure your traits are in the same order as your tree
  trait <- trait %>%
    rownames_to_column('taxon') %>%
    arrange(match(taxon, phy$tip.label)) %>%
    column_to_rownames('taxon')
  
  # set the traits into vectors
  t1 <- setNames(trait[,1],rownames(trait))
  t2 <- setNames(trait[,2],rownames(trait))
  
  # estimate the sigma of a constant rate BM for each trait
  if(trait.correlation=="no" | trait.correlation=="both"){
    sig.constant.t1 <- geiger::fitContinuous(phy, t1, model="BM")$opt$sigsq
    sig.constant.t2 <- geiger::fitContinuous(phy, t2, model="BM")$opt$sigsq 
  }
  if(trait.correlation=="yes" | trait.correlation=="both"){
    mv.bm.fit <- mvBM(phy, trait, model="BM1")
  }

  
  # estimate ancestral states at nodes from the scaled phylogeny
  anc <- data.frame(t1.anc=matrix(fastAnc(scaled.phy, t1)),
                    t2.anc=matrix(fastAnc(scaled.phy, t2)))
  colnames(anc) <- colnames(trait); 
  rownames(anc) <- (Ntip(phy)+1):(Nnode(phy)+Ntip(phy))
  
  # simulate N trait sets
  if(trait.correlation=="no" | trait.correlation=="both"){
    sim.t1 <- as.data.frame(fastBM(phy, a=anc[1,1], sig2=sig.constant.t1, n=sim.num))
    sim.t2 <- as.data.frame(fastBM(phy, a=anc[1,2], sig2=sig.constant.t2, n=sim.num))
    sim.trait.bm <- data.frame(t1=unlist(sim.t1), t2=unlist(sim.t2))
  }
  sim.trait.mvbm <- NULL
  if(trait.correlation=="yes" | trait.correlation=="both"){
    sim.both <- lapply(1:sim.num, function(x) as.data.frame(mvSIM(phy, nsim=1, model="BM1", param=mv.bm.fit)))
    for (k in 1:length(sim.both)){sim.trait.mvbm <- rbind(sim.trait.mvbm, sim.both[[k]])}
  }

  # create a phylomorphospace object using phytools
  obj <- phylomorphospace(tree=phy, X=trait, A=anc, 
                          ftype="off", bty="n", node.size=c(0,1))
  # extract the xy coordinates of the nodes and segments
  phydat <- data.frame(xstart=obj$xx[obj$edge[,1]],
                       ystart=obj$yy[obj$edge[,1]],
                       xstop=obj$xx[obj$edge[,2]],
                       ystop=obj$yy[obj$edge[,2]],
                       nodestart=obj$edge[,1],
                       nodestop=obj$edge[,2])
  
  # log and scale the module rates
  rate.df <- data.frame(raw.rate = branch.rates$Mean.SigV, 
                        nodestart=branch.rates$parent.node, 
                        nodestop=branch.rates$child.node)
  rate.df$rate <- log(rate.df$raw.rate)
  rate.df$scaled <- round((rate.df$rate - min(rate.df$rate))/diff(range(rate.df$rate)) * 99) + 1
  # choose colors and apply them to the scaled rates
  colors <- (colorRampPalette(brewer.pal(9, brew.pal)[2:9])(100))
  if(brew.pal %in% c("Spectral","RdYlGn","RdYlBu","RdGy","RdBu","PuOr","PRGn","PiYG","BrBG")){colors <- rev(colors)}
  
  #rownames(rate.df)[Ntip(phy):nrow(rate.df)] <- (1:length(phy$tip.label))
  #rate.df$nodestart <- as.numeric(as.character(rownames(rate.df)))
  phydat <- inner_join(phydat, rate.df)
  phydat$color <- colors[phydat$scaled]
  
  # if you're working with a clade
  if(focus == "clade"){
    # determine the MRCA node of your target tips
    mrcn <- getMRCA(phy, tip = tip.spread)
    # get all the descendant nodes of your MRCA including internals
    all.desc <- getDescendants(phy, mrcn)    
  }
  # if you're working with a single tip
  if(focus == "tip"){all.desc <- which(phy$tip.label==tip.spread[[1]])}
  
  # identify which descendant nodes are tips
  tip.desc <- all.desc[which(all.desc <= Ntip(phy))]
  # get the node path from the root to each tip, and make a vector of the unique nodes
  focal.nodes <- unique(unlist(sapply(tip.desc, function(x) nodepath(phy, from=Ntip(phy)+1, to=x))))
  # subset the phydat dataframe to just the nodes of interest
  focal.dat <- dplyr::filter(phydat, nodestart %in% focal.nodes & nodestop %in% focal.nodes)
  
  # plot the branches from root to tips colored by evolutionary rate, with simulated data
  ggplot()+
    {if(!is.null(sim.trait.bm))geom_hdr(data=sim.trait.bm, aes(x=sim.trait.bm[,1], y=sim.trait.bm[,2]), method="mvnorm", fill="lightGrey", probs=c(0.25,0.5,0.8,0.95))} +
    {if(!is.null(sim.trait.mvbm))geom_hdr_lines(data=sim.trait.mvbm, aes(x=sim.trait.mvbm[,1], y=sim.trait.mvbm[,2]), method="mvnorm", color="#3da5d1", probs=c(0.25,0.5,0.8,0.95))} +
    #geom_hdr_lines(data=sim.handfoot, aes(x=Hand, y=Foot), method="mvnorm") +
    geom_segment(data=phydat, aes(x=xstart,y=ystart,xend=xstop,yend=ystop), lwd=lsize-1, color="grey") +
    geom_point(data=trait, aes(x=trait[,1], y=trait[,2]), size=psize, shape=21, fill="lightGrey", color="#757574") +
    geom_point(data=focal.dat[1,], aes(x=xstart,y=ystart), size=psize+1, shape=3, color="black") +
    geom_point(data=focal.dat[1,], aes(x=xstart,y=ystart), size=psize+1, shape=4, color="black") +
    #geom_point(data=anc[which(rownames(anc)==mrcn),], aes(x=anc[,1], y=anc[,2]), size=psize+1, shape=21, fill="white", color="black") +
    geom_segment(data=focal.dat,aes(x=xstart,y=ystart,xend=xstop,yend=ystop),color=focal.dat$color, lwd = lsize) +
    geom_point(data=focal.dat, aes(x=xstart, y=ystart), size=psize, shape=21, fill="white", color="black") +
    geom_point(data=focal.dat, aes(x=xstop, y=ystop),   size=psize, shape=21, fill="white", color="black") +
    labs(x=colnames(trait)[[1]], y=colnames(trait)[[2]]) +
    theme_classic()
  
}
# e.g.:
#morphotrajectory.VR(phy=egernia.tree,
#                    scaled.phy=all.BT$scalar.trees$Foot,
#                    trait=smallLSR[,c("Hand","Foot")],
#                    branch.rates=all.BT$all.res$Foot,
#                    tip.spread=c("Tiliqua_rugosa","Tiliqua_gigas"),
#                    focus="clade", sim.num=100,
#                    trait.correlation="both",
#                    brew.pal="YlGnBu", psize=3, lsize=2)
# 'trait.correlation="yes"' will use mvMORPH to fit multivariate BM
# to the traits, then simulated correlated traits under the fitted
# mvBM parameter estimates
# 'trait.correlation="no"' will assume independent evolution of traits
# under BM
# 'trait.correlations="both"' will do both!


# this function will plot empirical data ontop of a simulated distribution of traits
# under both correlated and uncorrelated BM
# the main difference from 'morphotrajectory.VR' is that it will not plot a phylomorphospace
# and the simulated traits have different root values
sim.morphotrajectory <- function(phy, trait, sim.num){
  # check that the tip names are correct
  if(!all(tip.spread %in% phy$tip.label)){stop("WARNING: check tip names")}
  # check that all the tips are in the trait dataframe
  if(nrow(trait) != length(phy$tip))
    stop("your trait dataframe must have the same number of rows as species in your tree")
  
  # make sure your traits are in the same order as your tree
  trait <- trait %>%
    rownames_to_column('taxon') %>%
    arrange(match(taxon, phy$tip.label)) %>%
    column_to_rownames('taxon')
  
  # set the traits into vectors
  t1 <- setNames(trait[,1],rownames(trait))
  t2 <- setNames(trait[,2],rownames(trait))
  
  # estimate the sigma of a constant rate BM for each trait
  sig.constant.t1 <- geiger::fitContinuous(phy, t1, model="BM")$opt$sigsq
  sig.constant.t2 <- geiger::fitContinuous(phy, t2, model="BM")$opt$sigsq 
  # estimate the parameters of a multivariate correlated BM model
  mv.bm.fit <- mvBM(phy, trait, model="BM1", echo=F)
  
  # simulate N trait sets under multivariate uncorrelated BM
  sim.t1 <- as.data.frame(fastBM(phy, a=mv.bm.fit$theta[[1]], sig2=sig.constant.t1, n=sim.num))
  sim.t2 <- as.data.frame(fastBM(phy, a=mv.bm.fit$theta[[2]], sig2=sig.constant.t2, n=sim.num))
  sim.trait.bm <- data.frame(t1=unlist(sim.t1), t2=unlist(sim.t2))
  # simulate N trait sets under multivariate correlated BM
  sim.both <- lapply(1:sim.num, function(x) as.data.frame(mvSIM(phy, nsim=1, model="BM1", param=mv.bm.fit)))
  sim.trait.mvbm <- NULL
  for (k in 1:length(sim.both)){sim.trait.mvbm <- rbind(sim.trait.mvbm, sim.both[[k]])}

  # plot the branches from root to tips colored by evolutionary rate, with simulated data
  ggplot()+
    geom_hdr(data=sim.trait.bm, aes(x=sim.trait.bm[,1], y=sim.trait.bm[,2]), method="mvnorm", fill="lightGrey", probs=c(0.25,0.5,0.8,0.95)) +
    #geom_hdr(data=sim.trait.mvbm, aes(x=sim.trait.mvbm[,1], y=sim.trait.mvbm[,2]), method="mvnorm", fill="#ABDDA4", probs=c(0.25,0.5,0.8,0.95)) +
    geom_hdr_lines(data=sim.trait.mvbm, aes(x=sim.trait.mvbm[,1], y=sim.trait.mvbm[,2]), method="mvnorm", color="#f26224", probs=c(0.25,0.5,0.8,0.95)) +
    geom_point(data=trait, aes(x=trait[,1], y=trait[,2]), pch=21, fill="white") +
    labs(x=colnames(trait)[[1]], y=colnames(trait)[[2]]) +
    theme_classic() + theme(legend.position = "bottom")
}

# e.g.: sim.morphotrajectory(phy=egernia.tree, trait=smallLSR[,c("Hand","Foot")], sim.num=100)

# the 'sim.to.node' function will plot the trajectory of a single taxon from root to tip
# showing the accumulation of variation at each internal node. This is plotted atop 
# a distribution of traits simulated under BM showing the possible trait values for 
# each node for a direct comparison of simulated to empirical node values.
sim.to.node <- function(phy, VRphy, trait, tip, sim.num){
  # get the path from the root to the focal tip
  tip.node <- which(phy$tip.label==tip)
  tip.path <- nodepath(phy, from=Ntip(phy)+1, to=tip.node)
  tip.path.named <- tip.path; tip.path.named[length(tip.path)] <- tip
  # pad the tip.path numbers with 0s so they plot in the right order
  tip.padded <- gsub("\\s", "0", format(tip.path, width=max(nchar(tip.path))))
  
  # estimate correlated multivariate BM on the empirical data
  mv.bm.fit.cor <- mvBM(VRphy, trait, model="BM1", echo=F, diagnostic=F)
  # estimate uncorrelated multivariate BM on the empirical data
  mv.bm.fit.unc <- mvBM(phy, trait, model="BM1", param=list(constraint="diagonal"), echo=F, diagnostic=F)
  # make a model object with theta values from correlated and sigma values from uncorrelated
  mv.bm.fit.hyb <- mv.bm.fit.unc; mv.bm.fit.hyb$theta <- mv.bm.fit.cor$theta
  
  # extract ancestral states for the input VR model
  emp.anc <- data.frame(estim(tree=VRphy, data=trait, object=mv.bm.fit.cor, asr=T)$estimates); 
  naked.trait.anc <- rbind(emp.anc, trait)
  emp.anc$node <- rownames(emp.anc)
  trait$node <- sapply(rownames(trait), function(x) which(phy$tip.label==x))
  trait.anc <- rbind(trait, emp.anc)
  
  # simulate data under the uncorrelated model with root theta values from the correlated model
  sim.list <- lapply(1:sim.num, function(x) data.frame(mvSIM(tree=phy, nsim=1, model="BM1", 
                                                             param=list(theta=mv.bm.fit.cor$theta, ntraits=2, 
                                                                        sigma=mv.bm.fit.unc$sigma))))
  # estimate ancestral states for each simulated dataset with root theta values from the correlated model
  sim.list <- lapply(1:sim.num, function(x) rbind(estim(tree=phy, data=sim.list[[x]], object=mv.bm.fit.hyb, asr=T)$estimates, sim.list[[x]]))
  
  # make a df of the simulated trait values at each node along the node path
  node.sims <- NULL
  for(k in 2:(length(tip.path)-1)){
    curr.node <- data.frame(extract.node.sim(sim.list, node=tip.path[[k]]))
    curr.node$node <- tip.padded[[k]]
    node.sims <- rbind(node.sims, curr.node)
  }
  tip.sim <- data.frame(extract.node.sim(sim.list, node=tip)); tip.sim$node <- tip
  node.sims <- rbind(node.sims, tip.sim)
  
  # make a df of the empirical trait values at each node along the node path
  node.emps <- NULL
  for(k in 1:(length(tip.path)-1)){
    curr.node <- trait.anc[which(rownames(trait.anc)==tip.path[[k]]),]
    curr.node$node <- tip.padded[[k]]
    node.emps <- rbind(node.emps, curr.node)
  }
  tip.emp <- trait[which(rownames(trait)==tip),]; tip.emp$node <- tip
  node.emps <- rbind(node.emps, tip.emp)
  
  # make two new columns so we can plot the trajectory as line segments
  # these are just the x/y coordinates from the previous node
  x.prev <- NULL
  for(k in 1:nrow(node.emps)){
    if(k==1){x.prev <- append(x.prev,node.emps[1,1])}
    if(k>1){x.prev <- append(x.prev,node.emps[k-1,1])}
  }
  node.emps$x.previous <- x.prev
  # see comment above
  y.prev <- NULL
  for(k in 1:nrow(node.emps)){
    if(k==1){y.prev <- append(y.prev,node.emps[1,2])}
    if(k>1){y.prev <- append(y.prev,node.emps[k-1,2])}
  }
  node.emps$y.previous <- y.prev
  
  # extend your preferred color palette to have more colors
  more.colors <- colorRampPalette(brewer.pal(11, "Spectral"))(length(tip.path))
  
  # plot our output
  ggplot() +
    geom_hdr_lines(data=node.sims, aes(x=X1, y=X2, color=node), method="mvnorm", probs=0.95, alpha=0.75) +
    #geom_hdr(data=node.sims, aes(x=X1, y=X2, fill=node), method="mvnorm", probs=0.95, alpha=0.1) +
    #scale_color_brewer(palette="Spectral") +
    scale_color_manual(values=more.colors) +
    geom_segment(data=node.emps, aes(x=get(names(node.emps)[[1]]), xend=x.previous, 
                                     y=get(names(node.emps)[[2]]), yend=y.previous), 
                 color="black", lwd=1, linetype="dotted") +
    geom_point(data=node.emps, aes(x=get(names(node.emps)[[1]]), y=get(names(node.emps)[[2]]), color=node), size=3) +
    xlab(names(trait.anc)[[1]]) + ylab(names(trait.anc)[[2]]) +
    theme_classic() + theme(legend.position="bottom")
}

# e.g.:
#sim.to.node(phy = egernia.tree,
#            VRphy = all.BT$scalar.trees$Limb,
#            trait = all.modules$limb[,c("Hand","Foot")],
#            tip = "Tiliqua_rugosa",
#            sim.num = 100)


# the 'distance.btwn.nodes' function shows the difference in trait variation accumulating 
# along a path from the root to a specified tip
# along the way, the ancestral states for each node are estimated under BM and VR
# in the multivariate case, distances are measured among nodes using mahalanobis
# distances measured to the root node
# then the plot shows the difference in trait accumulated between nodes compared
# between the observed and uncorrelated BM
distance.btwn.nodes <- function(phy, VRphy, trait, tip, sim.num, stat=c("confidence","quantile")){
  # get the path from the root to the focal tip
  tip.node <- which(phy$tip.label==tip)
  tip.path <- nodepath(phy, from=Ntip(phy)+1, to=tip.node)
  tip.path.named <- tip.path; tip.path.named[length(tip.path)] <- tip
  # pad the tip.path numbers with 0s so they plot in the right order
  tip.padded <- gsub("\\s", "0", format(tip.path, width=max(nchar(tip.path))))
  tip.padded.named <- tip.padded; tip.padded.named[length(tip.padded.named)] <- tip
  
  # estimate correlated multivariate BM on the empirical data
  mv.bm.fit.cor <- mvBM(VRphy, trait, model="BM1")
  # estimate uncorrelated multivariate BM on the empirical data
  mv.bm.fit.unc <- mvBM(phy, trait, model="BM1", param=list(constraint="diagonal"))
  # make a model object with theta values from correlated and sigma values from uncorrelated
  mv.bm.fit.hyb <- mv.bm.fit.unc; mv.bm.fit.hyb$theta <- mv.bm.fit.cor$theta
  
  # extract ancestral states for the input VR model
  emp.anc <- data.frame(estim(tree=VRphy, data=trait, object=mv.bm.fit.cor, asr=T)$estimates); 
  naked.trait.anc <- rbind(emp.anc, trait)

  
  # simulate data under the uncorrelated model with root theta values from the correlated model
  sim.list <- lapply(1:sim.num, function(x) data.frame(mvSIM(tree=phy, nsim=1, model="BM1", 
                                                             param=list(theta=mv.bm.fit.cor$theta, ntraits=ncol(trait), 
                                                                        sigma=mv.bm.fit.unc$sigma))))
  # estimate ancestral states for each simulated dataset with root theta values from the correlated model
  sim.list <- lapply(1:sim.num, function(x) rbind(estim(tree=phy, data=sim.list[[x]], object=mv.bm.fit.hyb, asr=T)$estimates, sim.list[[x]]))
  
  # combine the observed and ancestral trait values
  emp.anc$node <- rownames(emp.anc)
  trait$node <- sapply(rownames(trait), function(x) which(phy$tip.label==x))
  trait.anc <- rbind(trait, emp.anc)
  
  
  # I think there's two ways to do the next bit:
  # (1) use euclidean distances between pairs of points, though this can be biased
  # (2) use mahalanobis distances between each point and the root, then
  # get distances among pairs. But this isn't right I don't think.
  # I'm using (1) but will leave the code for (2) in place.
  # ultimately they give very similar results.
  
  # get multivariate distance among empirical nodes using mahalanobis distances
  #mahal.dist <- mahalanobis(x=naked.trait.anc, center=unlist(naked.trait.anc[1,]), cov=cov(naked.trait.anc))
  # get the mulitvariate distance among nodes using mahalanobis distances
  euc.pair <- NULL
  for(k in 1:(length(tip.path)-1)){
    # euclidean distances for multidimensional space
    euc.dists.sim <- unlist(lapply(sim.list, function(x) euclidean(x[which(rownames(x)==tip.path.named[[k]]),], 
                                                                   x[which(rownames(x)==tip.path.named[[k+1]]),])))
    
    # get maha for each node relative to the root node
    #mahal.sims <- lapply(sim.list, function(y) mahalanobis(x=y, center=unlist(y[1,]), cov=cov(y)))
    # determine the distance between pairs of focal nodes
    #euc.dists.sim <- unlist(lapply(mahal.sims, function(x) x[which(names(x)==tip.path.named[[k+1]])] - 
    #                                                       x[which(names(x)==tip.path.named[[k]])]))
    
    # get the mean and CI on each node
    if(stat=="confidence"){ed.cis <- Rmisc::CI(euc.dists.sim, ci=0.95)}
    if(stat=="quantile"){ed.cis <- quantile(euc.dists.sim, probs=c(0.95,0.5,0.05))}
    
    # make a dataframe of the required information
    euc.pair <- rbind(euc.pair, data.frame(pair=k, 
                                           timestart=nodeheight(phy,tip.path[[k]]),
                                           timestop =nodeheight(phy,tip.path[[k+1]]),
                                           upper=ed.cis[[1]], mean=ed.cis[[2]], lower=ed.cis[[3]],
                                           empirical=euclidean(naked.trait.anc[which(rownames(naked.trait.anc)==tip.path.named[[k+1]]),],
                                                               naked.trait.anc[which(rownames(naked.trait.anc)==tip.path.named[[k]]),]),
                                           #empirical=mahal.dist[which(names(mahal.dist)==tip.path.named[[k+1]])] - 
                                           #  mahal.dist[which(names(mahal.dist)==tip.path.named[[k]])],
                                           nodestart=tip.path[[k]],
                                           nodestop =tip.path[[k+1]]))
  }
  # add in a row representing the root
  euc.pair <- rbind(data.frame(pair=0, timestart=0, timestop=0, upper=0, mean=0, lower=0, empirical=0,
                               nodestart=57, nodestop=57), euc.pair)
  
  # extend your preferred color palette to have more colors
  more.colors <- colorRampPalette(brewer.pal(11, "Spectral"))(length(tip.path))
  euc.pair$node.name <- tip.padded.named
  
  # plot the resulting df
  ggplot() +
    #geom_segment(data=euc.pair, aes(x=timestart,xend=timestop,y=mean,yend=mean)) +
    geom_hline(yintercept=0, color="black", linetype="dotted") +
    geom_ribbon(data=euc.pair, aes(x=timestop, ymin=empirical-lower, ymax=empirical-upper), fill="lightGrey", alpha=0.25) +
    geom_line(data=euc.pair, aes(x=timestop,y=empirical-mean), color="grey") +
    #geom_point(data=euc.pair, aes(x=timestop,y=empirical-mean), pch=21, fill="white") +
    geom_point(data=euc.pair, aes(x=timestop,y=empirical-mean,color=node.name)) +
    scale_color_manual(values=more.colors) +
    xlab("Time (Millions of Years)") + ylab("Trait Difference Along Edge (observed - BM)") +
    theme_classic() + theme(legend.position="bottom")
    
}


# measure euclidean distance of multidimensional data
euclidean <- function(a, b) sqrt(sum((a - b)^2))

# extract node values from lots of simulated data
extract.node.sim <- function(xlist, node){
  curr.list <- lapply(xlist, function(x) x[which(rownames(x)==node),])
  curr.df <- as.data.frame(do.call(rbind, curr.list))
  return(curr.df)
}

  