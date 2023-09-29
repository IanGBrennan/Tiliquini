library("neldermead")
library("fBasics")
library("gsl")
library("stabledist")
library("statmod")
library("parallel")

library("pulsR")
library("geiger")


# make a data frame of the traits of interest
puLSR <- smallLSR[,1:20]

# establish the list of models we want to fit
model.vec <- c("BM","OU","EB","JN","NIG","BMJN","BMNIG")

# create an empty list to hold the results
trait.fit <- NULL
params.list <- NULL

# make a loop to fit the models to each trait
for (k in 1:ncol(puLSR)){
  # turn the trait into a vector
  trait.vec <- puLSR[,k]
  names(trait.vec) <- rownames(puLSR)
  
  # fit the models in parallel with mclapply
  res.list <- mclapply(1:7, function(x) {
    fit_reml_levy(phy = egernia.tree, dat = trait.vec, model = model.vec[[x]])}, mc.cores = 8)
  names(res.list) <- model.vec
  trait.fit[[k]] <- unlist(lapply(res.list, function(x) x$AIC))
  params.list[[k]] <- lapply(res.list, function(x) x$params)
}
names(trait.fit) <- colnames(puLSR)
names(params.list) <- colnames(puLSR)

#for(k in 1:length(trait.fit)){trait.fit[[k]] <- trait.fit[[k]][1:7]}



# if you have results from fitting the VarRates model separately, bring them in here.
for(k in 1:length(trait.fit)){
  curr.trait <- names(trait.fit[k])
  var.rate <- BT.AICs[curr.trait]; names(var.rate) <- "VarRates"
  trait.fit[[k]] <- append(trait.fit[[k]], var.rate)
}

for(k in 1:length(trait.fit)){
  curr.trait <- names(trait.fit[k])
  kappa.rate <- kappa.AIC[curr.trait][[1]]; names(kappa.rate) <- "Kappa"
  trait.fit[[k]] <- append(trait.fit[[k]], kappa.rate)
}

trait.aicw.all <- lapply(trait.fit, function(x) aic.w(x))
trait.fit.types <- NULL
for (j in 1:length(trait.fit)){
  trait.fit.types[[j]] <- trait.fit[[j]][c(1:3, which(trait.fit[[j]]==min(trait.fit[[j]][4:7])),8,9)]
  names(trait.fit.types[[j]]) <- c("BM","OU","EB","Jump","VarRates","Kappa")
}
trait.aicw.types <- lapply(trait.fit.types, function(x) aic.w(x))
names(trait.aicw.types) <- names(trait.fit)

# make a dataframe of the 'type' AICw
aicw.df <- NULL
for (k in 1:length(trait.aicw.types)){
  aicw.df <- rbind(aicw.df, data.frame(trait=names(trait.aicw.types)[[k]], aicw=as.vector(trait.aicw.types[[k]]), model=names(trait.aicw.types[[k]])))
}
# plot the 'type'  results
aicw.types.plot <- ggplot(aicw.df, aes(fill=model, y=aicw, x=trait)) + 
  geom_bar(position="stack", stat="identity") + theme_classic() + 
  scale_fill_brewer(palette="RdYlBu") +
  theme(axis.text.x = element_text(angle = 45, hjust=1), legend.position="bottom")

# make a dataframe of all the model AICw
aicw.all.df <- NULL
for (k in 1:length(trait.aicw.all)){
  aicw.all.df <- rbind(aicw.all.df, data.frame(trait=names(trait.aicw.all)[[k]], aicw=as.vector(trait.aicw.all[[k]]), model=names(trait.aicw.all[[k]])))
}
# plot all the AICw 
aicw.all.plot <- ggplot(aicw.all.df, aes(fill=model, y=aicw, x=trait)) + 
  geom_bar(position="stack", stat="identity") + theme_classic() + 
  scale_fill_brewer(palette="RdYlBu") +
  theme(axis.text.x = element_text(angle = 45, hjust=1), legend.position="bottom")

# save the results externally
save(trait.fit, params.list, trait.aicw.all, trait.fit.types, 
     trait.aicw.all, trait.aicw.types, aicw.df, aicw.all.df, file="Data/pulsR_Results.Rdata")


aicw.types.plot / aicw.all.plot

sig.fit <- NULL
for(k in 1:length(trait.aicw.types)){
  curr.res <- sort(trait.aicw.types[[k]], decreasing=T)
  if(curr.res[[1]]/curr.res[[2]] > 2){sig.fit <- append(sig.fit, names(curr.res)[[1]])}
  if(curr.res[[1]]/curr.res[[2]] <= 2){sig.fit <- append(sig.fit, "none preferred")}
}
names(sig.fit) <- names(trait.aicw.types)
sig.fit[order(names(sig.fit))]


