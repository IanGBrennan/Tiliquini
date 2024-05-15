# from Crouch & Tobias 2022 Ecology Letters 

visualizeBayesTraitsLog <- function(x, exclude.tree.number = TRUE){
  
  require(gridExtra)
  require(ggplot2)
  require(ggthemes)
  
  colnames(x) <- sapply(colnames(x), .fix.colnames)
  
  if(exclude.tree.number == TRUE){
    g <- grep("Tree.No", colnames(x))
    x <- x[,-g]
  }
  
  
  trace.plots <- lapply(colnames(x)[-1], .plot.single.trace, data = x)
  density.plots <- lapply(colnames(x)[-1], .plot.single.density, data = x)
  
  num.rows <- length(trace.plots)
  total.num.plots <- num.rows * 2
  plot.list <- vector(mode = "list", length = total.num.plots)
  trace.locations <- seq(1, total.num.plots, by=2)
  density.locations <- seq(2, total.num.plots, by=2)
  plot.list[trace.locations] <- trace.plots
  plot.list[density.locations] <- density.plots
  
  do.call("grid.arrange", c(plot.list, ncol=2))
  
}


.fix.colnames <- function(name){
  spl <- strsplit(name, "")[[1]]
  test <- spl=="^"
  back <- spl[test==FALSE]
  back.name <- paste(back, collapse="")
  back.name <- gsub(" ",".", back.name)
  return(back.name)
}



.plot.single.trace <- function(data, var){
  ggplot(data, aes_string(x = "Iteration", y = var)) +
    geom_line() +
    xlab("Iteration") +
    ylab(var) +
    theme_bw() +
    theme(legend.position="right", 
          axis.text=element_text(size=12),
          axis.title=element_text(size=12))
}



.plot.single.density <- function(data, var){
  ggplot(data, aes_string(x = var)) +
    geom_density() +
    xlab(var) +
    theme_bw() +
    theme(legend.position="none", 
          axis.text=element_text(size=12),
          axis.title=element_text(size=12))
}


