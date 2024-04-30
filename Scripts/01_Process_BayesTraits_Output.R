setwd("/Data/BayesTraits/")


############################################################################

# summarize the BayesTraits results for each morphological trait
tr.files <- dir(getwd(), pattern="_TraitResults.txt") # will give you the files for individual traits
mo.files <- dir(getwd(), pattern="_ModuleResults.txt") # will give you the files for modules

############################################################################

# process the outputs of BayesTraits
trait.BT  <- process_PPP(res.files = tr.files, phy = til.tree, col.palette="YlGnBu")
module.BT <- process_PPP(res.files = mo.files, phy = til.tree, col.palette="YlGnBu")

# the above will give you a series of trees and res objects
# all.BT$all.trees holds a tree for each trait. The tree is your input tree, rescaled by the mean estimated sigma (rate) value per branch (Mean.SigV in BayesTraits terms)
# all.BT$all.sig holds the mean sigma value per branch in a vector (?), this is the value each branch above was scaled by
# all.BT$all.res holds a data frame for each trait including the node and branch information, rate info, etc.
# all.BT$scalar.trees holds a tree for each trait. The tree is your input tree, rescaled by the median scalar r (Median.Scalar in BayesTraits terms)
# all.BT$mean.scalar.trees holds a tree for each trait. The tree is your input tree, rescaled by the mean scalar r (Mean.Scalar in BayesTraits terms)
# all.BT$rate.trees holds a tree for each trait. The tree is your input tree, with each branch rescaled to the mean estimated sigma value (Mean.SigV)

############################################################################

# identify focal trees
trait.trees <- trait.BT$mean.scalar.trees
module.trees <- module.BT$mean.scalar.trees

############################################################################

# write your data and trees to file
save(trait.trees, module.trees, file="Data/BayesTraits_TREES.RData")
save(trait.BT, file="Data/BayesTraits_processed_TRAITS.RData")
save(module.BT, file="Data/BayesTraits_processed_MODULES.RData")