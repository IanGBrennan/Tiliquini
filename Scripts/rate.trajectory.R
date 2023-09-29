require(phytools)
require(ggplot2)
require(patchwork)
require(dplyr)
require(tibble)
require(RColorBrewer)
source("/Users/ianbrennan/Documents/GitHub/Egernia_Evolution/Scripts/processRRrates.R")
source("/Users/ianbrennan/Documents/GitHub/Egernia_Evolution/Scripts/trait.at.time.R")

## RATE.TRAJECTORY HASN'T BEEN UPDATED WITH THE NEW METHOD THAT'S IN RATE.TRAJECTORY.MODULE!
# this function will plot a phytools-style phenogram for a trait
# with the root-tip path of a clade or tip colored according to 
# the evolutionary rate of those edges estimated by RRphylo.
# if requested, it will also return the data for plotting
rate.trajectory <- function(RR, trait, anc, branch.rates, tip.spread,
                            focus=c("tip","clade"), brew.pal="YlOrRd",
                            psize=3, lsize=2, background.color=c("grey","brew.pal"),
                            gimme.the.data=F){
  # set the tree
  phy <- RR$tree
  
  # check that the tip names are correct
  if(!all(tip.spread %in% phy$tip.label)){stop("WARNING: check tip names")}
  # check that all the tips are in the trait dataframe
  if(nrow(trait) != length(phy$tip))
    stop("your trait dataframe must have the same number of rows as species in your tree")
  
  # make sure your traits are in the same order as your tree
  module <- module %>%
    rownames_to_column('taxon') %>%
    arrange(match(taxon, phy$tip.label)) %>%
    column_to_rownames('taxon')
  
  ### THIS IS WHERE I"M AT, NEED TO MAKE A NEW FUNCTION FOR JUST ONE TRAIT
  # use the processRRrates function to get a dataframe with all the necessary info
  rates.df <- processRRrates.multi(RR = RR, col.palette = brew.pal, log.rates = log.rates)

  
  
  
  
  
  # dissolve your tree into edge segments with discrete start/stop
  edges <- data.frame(parent=phy$edge[,1], child=phy$edge[,2])
  edges$xstart <- sapply(edges$parent, function(x) nodeheight(phy, x))
  edges$xstop <- sapply(edges$child, function(x) nodeheight(phy, x))
  edges$xstart <- edges$xstart - max(nodeHeights(phy))
  edges$xstop <- edges$xstop - max(nodeHeights(phy))
  edges$xstop <- round(edges$xstop, 5)
  
  # match the trait values (y axis) to the nodes
  rownames(trait) <- sapply(rownames(trait), function(x) which(phy$tip.label==x))
  trait <- rbind(trait, anc)
  edges$ystart <- sapply(edges$parent, function(x) trait[which(rownames(trait)==x),1])
  edges$ystop  <- sapply(edges$child,  function(x) trait[which(rownames(trait)==x),1])
  rownames(branch.rates)[Ntip(phy):nrow(branch.rates)] <- (1:length(phy$tip.label))
  branch.rates <- abs(branch.rates)
  edges$ratestart <- sapply(edges$parent, function(x) branch.rates[which(rownames(branch.rates)==x),1])
  edges$ratestop  <- sapply(edges$child, function(x) branch.rates[which(rownames(branch.rates)==x),1])
  
  # log and scale the module rates
  #rate.df <- data.frame(raw.rate = branch.rates[,1])
  rate.df <- abs(branch.rates); colnames(rate.df)[[1]] <- "raw.rate"
  rate.df$rate <- log(rate.df$raw.rate)
  rate.df$scaled <- round((rate.df$rate - min(rate.df$rate))/diff(range(rate.df$rate)) * 99) + 1
  # choose colors and apply them to the scaled rates
  colors <- (colorRampPalette(brewer.pal(9, brew.pal)[2:9])(100))
  palette.direction <- 1
  if(brew.pal %in% c("Spectral","RdYlGn","RdYlBu","RdGy","RdBu","PuOr","PRGn","PiYG","BrBG")){colors <- rev(colors); palette.direction <- -1}
  
  rownames(rate.df)[Ntip(phy):nrow(rate.df)] <- (1:length(phy$tip.label))
  rate.df$parent <- as.numeric(as.character(rownames(rate.df)))
  edges <- left_join(edges, rate.df, by="parent")
  edges$color <- colors[edges$scaled]
  edges$trait <- colnames(trait)[[1]]
  
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
  focal.dat <- dplyr::filter(edges, parent %in% focal.nodes & child %in% focal.nodes)
  
  # make the plot
  rt.plot <- ggplot() +
    {if(background.color=="brew.pal")geom_segment(data=edges, aes(x=xstart,y=ystart,xend=xstop,yend=ystop), lwd=lsize-1, color=edges$color, alpha=0.75)} +
    {if(background.color=="grey")geom_segment(data=edges, aes(x=xstart,y=ystart,xend=xstop,yend=ystop), lwd=lsize-1, color="grey", alpha=0.5)} +
    geom_segment(data=focal.dat, aes(x=xstart,y=ystart,xend=xstop,yend=ystop), lwd=lsize+1, color="#706f6f") +
    geom_segment(data=focal.dat, aes(x=xstart,y=ystart,xend=xstop,yend=ystop), lwd=lsize, color=focal.dat$color) +
    geom_point(data=focal.dat, aes(x=xstart, y=ystart), size=psize, shape=21, fill="white", color="#706f6f") +
    geom_point(data=focal.dat, aes(x=xstop, y=ystop),   size=psize, shape=21, fill="white", color="#706f6f") +
    xlab("Million Years Ago") + ylab(paste(colnames(trait)[[1]])) +
    theme_classic()
  
  # get a color scale
  color.df$scale <- seq(min(rate.df$raw.rate), max(rate.df$raw.rate), length.out=100)
  
  scale <- ggplot() + geom_tile(data=color.df, aes(x=1, y=scale, fill=scale)) +
    scale_fill_distiller(palette = brew.pal, direction=palette.direction) +
    theme_classic() +  ylab("Evolutionary Rate") +
    theme(legend.position="none", 
          axis.ticks = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank())
  
  if(gimme.the.data==F){return(rt.plot + (plot_spacer() / scale) + plot_layout(widths = c(20, 1)))}
  
  if(gimme.the.data==T){return(edges)}
}
# e.g.:
#rate.trajectory(phy = egernia.tree,
#                trait = dplyr::select(all.modules$head, Head_Width),
#                anc = dplyr::select(data.frame(allRR$head$aces), Head_Width),
#                branch.rates = dplyr::select(data.frame(allRR$head$multiple.rates), Head_Width),
#                tip.spread = c("Tiliqua_rugosa", "Tiliqua_occipitalis"),
#                focus = "clade",
#                brew.pal = "YlOrRd",
#                psize = 3, lsize = 2, background.color = "grey")

# This function produces the node data similar to above, but for
# a module comprising several traits. It does not produce a plot
# but will return the extracted data.
# Maybe most importantly, this will estimate the euclidean distance
# between nodes in multidimensional space: column *euc.dist*
rate.trajectory.module <- function(RR, module, tip.spread,
                                   focus=c("tip","clade"), brew.pal="YlOrRd",
                                   psize=3, lsize=2, background.color=c("grey","brew.pal"),
                                   gimme.the.data=F, log.rates=F){
  
  # set the tree
  phy <- RR$tree
  
  # make sure your traits are in the same order as your tree
  module <- module %>%
    rownames_to_column('taxon') %>%
    arrange(match(taxon, phy$tip.label)) %>%
    column_to_rownames('taxon')
  
  # use the processRRrates function to get a dataframe with all the necessary info
  rates.df <- processRRrates.multi(RR = RR, col.palette = brew.pal, log.rates = log.rates)
  
  # combine all traits (contemporary and ancestors)
  all.trait <- rbind(RR$aces, module)
  rownames(all.trait)[Ntip(phy):nrow(all.trait)] <- 1:Ntip(phy)
  
  # resolve distances among multidimensional nodes via Euclidean distance
  tstart <- NULL; tstop <- NULL
  for (k in 1:nrow(rates.df)){
    curr.traits <- all.trait[which(rownames(all.trait) %in% c(unlist(rates.df[k,c("parent.node","child.node")]))),]
    tstart[k] <- round(norm(as.matrix(curr.traits[1,]), "F"),5)
    tstop[k]  <- round(norm(as.matrix(curr.traits[2,]), "F"),5)
  }
  rates.df$traitstart <- tstart
  rates.df$traitstop <- tstop
  rates.df$module <- colnames(RR$rates)
  
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
  focal.dat <- dplyr::filter(rates.df, parent.node %in% focal.nodes & child.node %in% focal.nodes)
  
  # make the plot
rt.plot <- ggplot() +
    {if(background.color=="brew.pal")geom_segment(data=rates.df, aes(x=timestart,y=traitstart,xend=timestop,yend=traitstop), lwd=lsize-1, color=rates.df$edge.color, alpha=0.75)} +
    {if(background.color=="grey"    )geom_segment(data=rates.df, aes(x=timestart,y=traitstart,xend=timestop,yend=traitstop), lwd=lsize-1, color="grey", alpha=0.5)} +
    geom_segment(data=focal.dat, aes(x=timestart,y=traitstart,xend=timestop,yend=traitstop), lwd=lsize+1, color="#706f6f") +
    geom_segment(data=focal.dat, aes(x=timestart,y=traitstart,xend=timestop,yend=traitstop), lwd=lsize, color=focal.dat$edge.color) +
    geom_point(data=focal.dat, aes(x=timestart, y=traitstart), size=psize, shape=21, fill="white", color="#706f6f") +
    geom_point(data=focal.dat, aes(x=timestop,  y=traitstop),  size=psize, shape=21, fill="white", color="#706f6f") +
    xlab("Million Years Ago") + ylab(paste(colnames(RR$rates))) +
    theme_classic()
  
  # get a color scale
  colors <- (colorRampPalette(brewer.pal(9, brew.pal))(100));
  palette.direction = 1
  if(brew.pal %in% c("Spectral","RdYlGn","RdYlBu","RdGy","RdBu","PuOr","PRGn","PiYG","BrBG")){colors <- rev(colors); palette.direction=-1}
  color.df <- data.frame(color = colors)
  color.df$scale <- seq(min(RR$rates), max(RR$rates), length.out=100)
  
  scale <- ggplot() + geom_tile(data=color.df, aes(x=1, y=scale, fill=scale)) +
    scale_fill_distiller(palette = brew.pal, direction=palette.direction) +
    theme_classic() +  ylab("Evolutionary Rate") +
    theme(legend.position="none", 
          axis.ticks = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank())
  
  if(gimme.the.data==F){return(rt.plot + (plot_spacer() / scale) + plot_layout(widths = c(20, 1)))}
  
  if(gimme.the.data==T){return(rates.df)}
  
  
  return(edges)
}
# e.g.: rate.trajectory.module(RR=allRR$head, module=all.modules$head, 
#                              tip.spread=c("Tiliqua_rugosa", "Tiliqua_occipitalis"), focus="clade", brew.pal="YlOrRd",
#                              psize=3, lsize=2, background.color="grey", gimme.the.data=F)


# This function will produce a traitgram/phenogram
# with edges/branches colored according to estimated evolutionary rates.
# You can also use it to export a data frame for use with rate.to.node()
rate.trajectory.BT<- function(tree, module, PPP.all.res, tip.spread,
                              focus=c("tip","clade"),
                              psize=3, lsize=2, background.color=c("grey","brew.pal"),
                              gimme.the.data=F, inset=F){
  
  # set the tree
  phy <- tree
  
  # use the processRRrates function to get a dataframe with all the necessary info
  rates.df <- PPP.all.res
  
  # combine all traits (contemporary and ancestors)
  all.trait <- module
  
  # get the trait values at the parent and child nodes (start and end of edge)
  rates.df$traitstart <- sapply(rates.df$parent.node, function(x) all.trait[which(names(all.trait)==x)])
  rates.df$traitstop  <- sapply(rates.df$Node.No,     function(x) all.trait[which(names(all.trait)==x)])
  
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
  focal.dat <- dplyr::filter(rates.df, parent.node %in% focal.nodes & child.node %in% focal.nodes)
  
  # make the plot
  rt.plot <- ggplot() +
    {if(background.color=="brew.pal")geom_segment(data=rates.df, aes(x=timestart,y=traitstart,xend=timestop,yend=traitstop), lwd=lsize-1, color=rates.df$edge.color, alpha=0.75)} +
    {if(background.color=="grey"    )geom_segment(data=rates.df, aes(x=timestart,y=traitstart,xend=timestop,yend=traitstop), lwd=lsize-1, color="grey", alpha=0.5)} +
    geom_segment(data=focal.dat, aes(x=timestart,y=traitstart,xend=timestop,yend=traitstop), lwd=lsize+1, color="#706f6f") +
    geom_segment(data=focal.dat, aes(x=timestart,y=traitstart,xend=timestop,yend=traitstop), lwd=lsize, color=focal.dat$edge.color) +
    geom_point(data=focal.dat, aes(x=timestart, y=traitstart), size=psize, shape=21, fill="white", color="#706f6f") +
    geom_point(data=focal.dat, aes(x=timestop,  y=traitstop),  size=psize, shape=21, fill="white", color="#706f6f") +
    xlab("Million Years Ago") + ylab("Trait Value") +
    theme_classic()
  
  # make a color ramp that suits the data based on the input colors from rates.df
  palette.direction = 1
  if(rates.df$palette[[1]] %in% c("magma", "inferno", "plasma", "viridis","cividis", "rocket", "mako")){
    colors <- viridis(n=100, option=col.palette)
  }else{
    col.ramp <- colorRampPalette(brewer.pal(9, rates.df$palette[[1]]))
    colors <- (col.ramp(100))
    if(col.palette %in% c("Spectral","RdYlGn","RdYlBu","RdGy","RdBu","PuOr","PRGn","PiYG","BrBG")){colors <- rev(colors); palette.direction=-1}
  }
  
  # plot the scale bar
  scale <- ggplot() + geom_tile(data=color.df, aes(x=1, y=scale, fill=scale)) +
    {if(rates.df$palette[[1]]%in%rownames(brewer.pal.info))scale_fill_distiller(palette = rates.df$palette[[1]], direction=palette.direction)} +
    {if(rates.df$palette[[1]]%in%c("magma", "inferno", "plasma", "viridis","cividis", "rocket", "mako"))scale_fill_viridis_c(option=rates.df$palette[[1]], direction=palette.direction)} +
    theme_classic() +  ylab("Evolutionary Rate") +
    theme(legend.position="none", 
          axis.ticks = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank())
  
  if(gimme.the.data==F & inset==F){return(rt.plot + (plot_spacer() / scale) + plot_layout(widths = c(20, 1)))}
  if(gimme.the.data==F & inset==T){return(rt.plot + inset_element(scale, left=0.2, bottom=0.6, top=1, right=0.4, align_to = "plot"))}
  
  if(gimme.the.data==T){return(rates.df)}
  
  
  return(edges)
}
# e.g.: rate.trajectory.BT(tree=egernia.tree, module=anc.all.traits$Head_Width, PPP.all.res=all.BT$all.res$HeadWidth, 
#                           tip.spread=c("Tiliqua_adelaidensis", "Tiliqua_scincoides"),
#                           focus="clade", psize=3, lsize=2, background.color="grey",
#                           gimme.the.data=F)


# This function will produce a traitgram/phenogram
# with edges/branches colored according to estimated evolutionary rates.
# You can also use it to export a data frame for use with rate.to.node()
rate.trajectory.BT.module<- function(tree, module, PPP.all.res, tip.spread,
                              focus=c("tip","clade"),
                              psize=3, lsize=2, background.color=c("grey","brew.pal"),
                              gimme.the.data=F, inset=F){
  
  # set the tree
  phy <- tree
  
  # make sure your traits are in the same order as your tree
  module <- module %>%
    rownames_to_column('taxon') %>%
    arrange(match(taxon, phy$tip.label)) %>%
    column_to_rownames('taxon')
  
  # use the processRRrates function to get a dataframe with all the necessary info
  rates.df <- PPP.all.res
  
  # combine all traits (contemporary and ancestors)
  all.trait <- module
  
  # resolve distances among multidimensional nodes via Euclidean distance
  tstart <- NULL; tstop <- NULL
  for (k in 1:nrow(rates.df)){
    curr.traits <- all.trait[which(rownames(all.trait) %in% c(unlist(rates.df[k,c("parent.node","Node.No")]))),]
    tstart[k] <- round(norm(as.matrix(curr.traits[which(rownames(curr.traits)==rates.df[k,"parent.node"]),]), "F"),5)
    tstop[k]  <- round(norm(as.matrix(curr.traits[which(rownames(curr.traits)==rates.df[k,"Node.No"]),]), "F"),5)
  }
  rates.df$traitstart <- tstart
  rates.df$traitstop <- tstop
  #rates.df$module <- colnames(RR$rates)
  
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
  focal.dat <- dplyr::filter(rates.df, parent.node %in% focal.nodes & child.node %in% focal.nodes)
  
  # make the plot
  rt.plot <- ggplot() +
    {if(background.color=="brew.pal")geom_segment(data=rates.df, aes(x=timestart,y=traitstart,xend=timestop,yend=traitstop), lwd=lsize-1, color=rates.df$edge.color, alpha=0.75)} +
    {if(background.color=="grey"    )geom_segment(data=rates.df, aes(x=timestart,y=traitstart,xend=timestop,yend=traitstop), lwd=lsize-1, color="grey", alpha=0.5)} +
    geom_segment(data=focal.dat, aes(x=timestart,y=traitstart,xend=timestop,yend=traitstop), lwd=lsize+1, color="#706f6f") +
    geom_segment(data=focal.dat, aes(x=timestart,y=traitstart,xend=timestop,yend=traitstop), lwd=lsize, color=focal.dat$edge.color) +
    geom_point(data=focal.dat, aes(x=timestart, y=traitstart), size=psize, shape=21, fill="white", color="#706f6f") +
    geom_point(data=focal.dat, aes(x=timestop,  y=traitstop),  size=psize, shape=21, fill="white", color="#706f6f") +
    xlab("Million Years Ago") + ylab("Trait Value") +
    theme_classic()
  
  # make a color ramp that suits the data based on the input colors from rates.df
  palette.direction = 1
  if(rates.df$palette[[1]] %in% c("magma", "inferno", "plasma", "viridis","cividis", "rocket", "mako")){
    colors <- viridis(n=100, option=rates.df$palette[[1]])
  }else{
    col.ramp <- colorRampPalette(brewer.pal(9, rates.df$palette[[1]]))
    colors <- (col.ramp(100))
    if(rates.df$palette[[1]]%in% c("Spectral","RdYlGn","RdYlBu","RdGy","RdBu","PuOr","PRGn","PiYG","BrBG")){colors <- rev(colors); palette.direction=-1}
  }
  color.df <- data.frame(color = colors)
  color.df$scale <- seq(min(rates.df$Mean.SigV), max(rates.df$Mean.SigV), length.out=100)
  
  # plot the scale bar
  scale <- ggplot() + geom_tile(data=color.df, aes(x=1, y=scale, fill=scale)) +
    {if(rates.df$palette[[1]]%in%rownames(brewer.pal.info))scale_fill_distiller(palette = rates.df$palette[[1]], direction=palette.direction)} +
    {if(rates.df$palette[[1]]%in%c("magma", "inferno", "plasma", "viridis","cividis", "rocket", "mako"))scale_fill_viridis_c(option=rates.df$palette[[1]], direction=palette.direction)} +
    theme_classic() +  ylab("Evolutionary Rate") +
    theme(legend.position="none", 
          axis.ticks = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank())
  
  if(gimme.the.data==F & inset==F){return(rt.plot + (plot_spacer() / scale) + plot_layout(widths = c(20, 1)))}
  if(gimme.the.data==F & inset==T){return(rt.plot + inset_element(scale, left=0.2, bottom=0.6, top=1, right=0.4, align_to = "plot"))}
  
  if(gimme.the.data==T){return(rates.df)}
  
  
  return(edges)
}
# e.g.: rate.trajectory.BT.module(tree=egernia.tree, module=HEAD.anc.module, PPP.alla.res=all.BT.all.res$Head,
#                           tip.spread=c("Tiliqua_adelaidensis", "Tiliqua_scincoides"),
#                           focus="clade", psize=3, lsize=2, background.color="grey",
#                           gimme.the.data=F)
# HEAD.anc.module <- data.frame(Head_Width = anc.all.traits$Head_Width, 
#                               Snout_Eye = anc.all.traits$Snout_Eye, 
#                               Head_Depth = anc.all.traits$Head_Depth, 
#                               Pos_Skull = anc.all.traits$Pos_Skull)


# this function subsets the data frame produced by *rate.trajectory.module*
# down to just the edges which have phenotypic shifts inferred by *l1ou*.
# It sounds overly simplistic but l1ou operates in edge postorder for the 
# phy object, whereas everyone else just uses the normal order.
shifted.edges.l1ou <- function(phy, l1ou.object, rate.traj.obj){
  # map the edges of your preferred tree
  edge.map <- data.frame(edge = 1:nrow(phy$edge),
                         parent = phy$edge[,1],
                         child = phy$edge[,2])
  # map the edges of the l1ou tree whose edges are in postorder
  edge.map.l1ou <- data.frame(po.edge = 1:nrow(l1ou.object$tree$edge),
                              parent = l1ou.object$tree$edge[,1],
                              child = l1ou.object$tree$edge[,2])
  # combine the two maps to be able to translate
  combo.edge.map <- inner_join(edge.map, edge.map.l1ou)
  
  # identify the edges in our preferred tree which have been shifted
  l1ou.shifts <- dplyr::filter(combo.edge.map, po.edge %in% l1ou.object$shift.configuration)
  
  # filter the rate.trajectory dataframe to just the filtered edges
  shifts <- dplyr::filter(rate.traj.obj, child.node %in% l1ou.shifts$child)
  
  return(shifts)
}


# plot the evolutionary rate along branches leading to a rate shift
# as inferred by BayesTraits
rate.to.node.BT <- function(phy, PPP.obj, psize=3, lsize=2,
                            log.rates=F, col.palette="YlOrRd"){
  # if we want to log the rates
  if(log.rates==T){
    PPP.obj$Mean.SigV <- log(PPP.obj$Mean.SigV)
    PPP.obj$ratestart <- log(PPP.obj$ratestart)
    PPP.obj$ratestop  <- log(PPP.obj$ratestop)
    PPP.obj$edge.rate <- log(PPP.obj$edge.rate)
  }

  # scale the rates by the fastest rate (for ease of plotting)
  PPP.obj$rounded.rates <- round((PPP.obj$edge.rate - min(PPP.obj$edge.rate))/diff(range(PPP.obj$edge.rate)) * 99) + 1
  
  # sort it by the original order!
  PPP.obj <- PPP.obj[order(match(PPP.obj$child.node, phy$edge[,2])),]  
  
  # make a color ramp that suits the data
  if(col.palette %in% c("magma", "inferno", "plasma", "viridis","cividis", "rocket", "mako")){
    new.cols <- viridis(n=100, option=col.palette)
  }else{
    col.ramp <- colorRampPalette(brewer.pal(9, col.palette)[2:9])
    new.cols <- (col.ramp(100))
    if(col.palette %in% c("Spectral","RdYlGn","RdYlBu","RdGy","RdBu","PuOr","PRGn","PiYG","BrBG")){new.cols <- rev(new.cols)}
  }
  
  # correct the edge colors
  PPP.obj$edge.color <- new.cols[PPP.obj$rounded.rates]
  
  # identify shifted edges
  #shifts <- filter(PPP.obj, Median.Scalar > 2 & Pct.time.scaled > 70)
  shifts <- filter(PPP.obj, Mean.Scalar > 2 & Pct.time.scaled > 70)
  if(nrow(shifts)==0){stop("NO SHIFTS IN THIS TRAIT")}
  
  # identify the background edges
  background <- filter(PPP.obj, !edge %in% shifts$edge)
  
  # get the path to the shift(s)
  edge.shifts <- NULL
  for(j in 1:nrow(shifts)){
    # get the path from root to the node of interest
    shift.nodes <- nodepath(phy, Ntip(phy)+1, shifts[j,"child.node"])
    # get the node path from root to shifted edge
    for(k in 1:(length(shift.nodes)-1)){
      curr.edge <- dplyr::filter(PPP.obj, parent.node==shift.nodes[k] & child.node==shift.nodes[k+1])
      curr.edge$shift <- paste0("edge_", shifts[j,"edge"],"_childnode_", shifts[j,"child.node"])
      edge.shifts <- rbind(edge.shifts, curr.edge)
    }  
  }
  edge.shifts.unique <- edge.shifts %>% dplyr::select(!shift) %>% dplyr::distinct()
  
  
  # get a rate at time object from a dataframe
  RaT <- rate.at.time.df(timeslices=0.1, obj = background, plot=F)
  # extract the rates through time
  RaT <- extract.stat(RaT, stat="mean", range="quantile", plot=F)
  
  ggplot()+
    #{if(log.rates==T)geom_ribbon(data=RaT, aes(x=time-max(time), ymin=log(`5%`), ymax=log(`95%`)), fill="grey", alpha=0.5)} +
    #{if(log.rates==F)geom_ribbon(data=RaT, aes(x=time-max(time), ymin=`5%`, ymax=`95%`), fill="grey", alpha=0.5)} +
    #{if(log.rates==T)geom_line(data=RaT, aes(x=time-max(time), y=log(rate)), color="darkGrey", lwd=lsize-1)} +
    #{if(log.rates==F)geom_line(data=RaT, aes(x=time-max(time), y=rate), color="darkGrey", lwd=lsize-1)} +
    geom_ribbon(data=RaT, aes(x=time-max(time), ymin=`5%`, ymax=`95%`), fill="grey", alpha=0.5) +
    geom_line(data=RaT, aes(x=time-max(time), y=rate), color="darkGrey", lwd=lsize) +
    # the two commented lines will plot the trajectory of every branch
    #geom_segment(data=rate.traj.obj, aes(x=timestart, xend=timestart, y=ratestart, yend=ratestop), lwd=lsize, color="lightBlue", alpha=0.5, lineend="round") +
    #geom_segment(data=rate.traj.obj, aes(x=timestart, xend=timestop,  y=edge.rate, yend=edge.rate),lwd=lsize, color="lightBlue", alpha=0.5, lineend="round") +
    geom_segment(data=edge.shifts.unique, aes(x=timestart, xend=timestart, y=ratestart, yend=ratestop),  lwd=lsize, color=edge.shifts.unique$edge.color, lineend="round") +
    geom_segment(data=edge.shifts.unique, aes(x=timestart, xend=timestop,  y=edge.rate, yend=edge.rate), lwd=lsize, color=edge.shifts.unique$edge.color, lineend="round") +
    geom_point(data=edge.shifts.unique,   aes(x=timestart, y=ratestart), size=psize, shape=21, fill="white", color="#706f6f") +
    geom_point(data=edge.shifts.unique,   aes(x=timestop,  y=ratestop),  size=psize, shape=21, fill="white", color="#706f6f") +
    geom_point(data=shifts,        aes(x=timestop,  y=ratestop),  size=psize-1, shape=20, fill="black") +
    {if(log.rates==T)ylab("log Evolutionary Rate")} +
    {if(log.rates==F)ylab("Evolutionary Rate")} +
    xlab("Million Years Ago") +
    theme_classic()
}
# e.g.: rate.to.node.BT(phy=egernia.tree, PPP.obj=all.BT$all.res$Interlimb, psize=3, lsize=2, rates.logged=F)

# One more function to plot the evolutionary rate of branches leading to
# a shift inferred from l1ou
rate.to.node <- function(phy, rate.traj.obj, shift.obj, RaT, 
                         plot.type=c("slanted","square"), lsize=2, psize=3,
                         rates.logged=T, shift.method=c("l1ou","PhyloEM")){
  
  if(shift.method=="l1ou"){
    # identify the edges which have shifted evolutionary rates (l1ou_postorder --> normal order)
    shifts <- shifted.edges.l1ou(phy = phy,
                                 l1ou.object = shift.obj,
                                 rate.traj.obj = rate.traj.obj) 
  }
  if(shift.method=="PhyloEM"){
    shifts <- dplyr::filter(rate.traj.obj, edge %in% shift.obj$shift$edges)
  }

  # filter down to clades, removing single tips
  #shifts <- dplyr::filter(shifts, child.node > Ntip(phy))
  
  # we have to repeat this process for each of the shifts
  edge.shifts <- NULL
  for (j in 1:nrow(shifts)){
    # get the path from root to the node of interest
    shift.nodes <- nodepath(phy, Ntip(phy)+1, shifts[j,"child.node"])
    # get the node path from root to shifted edge
    for(k in 1:(length(shift.nodes)-1)){
      curr.edge <- dplyr::filter(rate.traj.obj, parent.node==shift.nodes[k] & child.node==shift.nodes[k+1])
      curr.edge$shift <- paste0("edge_", shifts[j,"child.node"])
      edge.shifts <- rbind(edge.shifts, curr.edge)
    }
  }
  
  # plot the output
  if(plot.type=="slanted"){
    ggplot() +
      {if(rates.logged==T)geom_ribbon(data=RaT, aes(x=time-max(time), ymin=log(`5%`), ymax=log(`95%`)), fill="grey", alpha=0.5)} +
      {if(rates.logged==F)geom_ribbon(data=RaT, aes(x=time-max(time), ymin=`5%`, ymax=`95%`), fill="grey", alpha=0.5)} +
      {if(rates.logged==T)geom_line(data=RaT, aes(x=time-max(time), y=log(rate)), color="darkGrey", lwd=lsize-1)} +
      {if(rates.logged==F)geom_line(data=RaT, aes(x=time-max(time), y=rate), color="darkGrey", lwd=lsize-1)} +
      geom_segment(data=edge.shifts, aes(x=timestart, y=ratestart, xend=timestop, yend=ratestop), color=edge.shifts$edge.color, lwd=lsize) +
      geom_point(data=edge.shifts,   aes(x=timestart, y=ratestart), size=psize, shape=21, fill="white", color="#706f6f") +
      geom_point(data=edge.shifts,   aes(x=timestop,  y=ratestop),  size=psize, shape=21, fill="white", color="#706f6f") +
      xlab("Million Years Ago") + ylab("Evolutionary Rate") +
      theme_classic()
  }
  if(plot.type=="square"){
    ggplot()+
      {if(rates.logged==T)geom_ribbon(data=RaT, aes(x=time-max(time), ymin=log(`5%`), ymax=log(`95%`)), fill="grey", alpha=0.5)} +
      {if(rates.logged==F)geom_ribbon(data=RaT, aes(x=time-max(time), ymin=`5%`, ymax=`95%`), fill="grey", alpha=0.5)} +
      {if(rates.logged==T)geom_line(data=RaT, aes(x=time-max(time), y=log(rate)), color="darkGrey", lwd=lsize-1)} +
      {if(rates.logged==F)geom_line(data=RaT, aes(x=time-max(time), y=rate), color="darkGrey", lwd=lsize-1)} +
      # the two commented lines will plot the trajectory of every branch
      #geom_segment(data=rate.traj.obj, aes(x=timestart, xend=timestart, y=ratestart, yend=ratestop), lwd=lsize, color="lightBlue", alpha=0.5, lineend="round") +
      #geom_segment(data=rate.traj.obj, aes(x=timestart, xend=timestop,  y=edge.rate, yend=edge.rate),lwd=lsize, color="lightBlue", alpha=0.5, lineend="round") +
      geom_segment(data=edge.shifts, aes(x=timestart, xend=timestart, y=ratestart, yend=ratestop),  lwd=lsize, color=edge.shifts$edge.color, lineend="round") +
      geom_segment(data=edge.shifts, aes(x=timestart, xend=timestop,  y=edge.rate, yend=edge.rate), lwd=lsize, color=edge.shifts$edge.color, lineend="round") +
      geom_point(data=edge.shifts,   aes(x=timestart, y=ratestart), size=psize, shape=21, fill="white", color="#706f6f") +
      geom_point(data=edge.shifts,   aes(x=timestop,  y=ratestop),  size=psize, shape=21, fill="white", color="#706f6f") +
      geom_point(data=shifts,        aes(x=timestop,  y=ratestop),  size=psize-1, shape=20, fill="black") +
      xlab("Million Years Ago") + ylab("Evolutionary Rate") +
      theme_classic()
  }
}
# e.g.: rate.to.node(phy = egernia.tree, rate.traj.obj = mod.limb, l1ou.obj = all.l1ou$limb, RaT = every.trait.rate$limb, plot.type="square")
# applying a gradient color to line is technically possible, but would probably be complicated
# see: https://stackoverflow.com/questions/37241893/is-it-possible-to-apply-color-gradient-to-geom-smooth-with-ggplot-in-r


