source("Scripts/innovate_elaborate.R")

############################################################################

# load the basic data including the tree
load("Data/Tiliquini_Data.RData")

# load the ancestral trait data
load("Data/Tiliquini_Ancestral_TRAITS.Rdata")

############################################################################


# flatten all the ancestral trait data into a single dataframe
anc.all.df <- data.frame(anc.all.traits)

############################################################################


# start with the full dataset
total.ie.phy <- innovate.elaborate.phy(trait.df = anc.all.df, 
                                       phy=til.tree, plot=T, summary=T, 
                                       angles="equal", PCs=c(1,2))
