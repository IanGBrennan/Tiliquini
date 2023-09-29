
for(j in 1:length(all.modules)){
  curr.mod <- all.modules[[j]]
  for (k in 1:length(curr.mod)){
    trait <- select(curr.mod, k)
    tname <- colnames(trait)
    filename <- paste0("~/Desktop/", tname, ".txt")
    write.table(trait, file=filename, sep=" ", quote=F, col.names=F)
  }
}
