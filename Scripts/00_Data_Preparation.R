# set the working directory to the GitHub repo
setwd("/Users/ianbrennan/Documents/GitHub/Tiliquini")

############################################################################

# Load necessary packages
library(dplyr)
library(ggplot2)
library(phytools)
library(RColorBrewer)
library(tidyr)

############################################################################

# Read in the raw morphological data
edata <- read.csv("Data/Tiliquini_Morphology.csv", header=T)

# Select the variables of interest
xdata <- edata %>%
  dplyr::select(Genus_species, Snout_Vent, Snout_Axilla, Interlimb,
                Body_Width, Pelvic_Width, Pelvic_Height, 
                Head_Width, Head_Length, Snout_Eye, Eye_Diameter, Head_Depth, 
                Tail_Width, Tail_Length, Tail_Loss,
                Arm, Upper_Arm, Lower_Arm, Hand,
                Leg, Upper_Leg, Lower_Leg, Foot)
# exclude incomplete samples
xdata <- xdata[complete.cases(xdata),]

############################################################################

# Some measurements are can be split into multiple variables

# determine neck length
xdata <- dplyr::mutate(xdata, Neck = Snout_Axilla - Head_Length)

# split head_length into tripartite model: snout, eye, posterior skull
xdata <- dplyr::mutate(xdata, Pos_Skull = Head_Length - (Snout_Eye + Eye_Diameter))

# redefine existing "snout" measurement to exlude eye diameter
#xdata <- mutate(xdata, Snout = Snout_Eye - Eye_Diameter) 

# identify the diff between SVL and ILL
xdata <- dplyr::mutate(xdata, Pelvic_Gap = Snout_Vent - (Interlimb + Snout_Axilla)) 

# We might want to drop the variables that you split (only if you want)
#xdata <- dplyr::select(xdata, -Snout_Axilla, -Head_Length, -Arm, -Leg)

############################################################################

# Create a tibble to get the species means  
# but, keep Tail_Length out, we'll handle it separately
sp.means <- xdata %>%
  group_by(Genus_species) %>%
  summarise_at(vars(Snout_Vent, Snout_Axilla, Interlimb,
                    Body_Width, Pelvic_Width, Pelvic_Height, 
                    Head_Width, Head_Length, Snout_Eye, Eye_Diameter, Head_Depth, 
                    Tail_Width, #Tail_Length, 
                    Arm, Upper_Arm, Lower_Arm, Hand,
                    Leg, Upper_Leg, Lower_Leg, Foot,
                    Neck, Pos_Skull, Pelvic_Gap), mean)

# We want to summarize tail length only for original tails
tails <- filter(xdata, Tail_Loss == "N")
# make sure there we aren't going to lose any species
#setdiff(xdata$Genus_species, tails$Genus_species)
sp.tails <- tails %>%
  group_by(Genus_species) %>%
  summarise_at(vars(Tail_Length), mean)

# Now add the Tail_Length variable back in
sp.means <- left_join(sp.means, sp.tails)

# Provide a little more info to the data
sp.means <- sp.means %>% 
  separate_wider_delim(Genus_species, delim = " ", names = c("Genus", "Species")) %>%
  mutate(Genus_species = paste0(Genus,"_",Species))

# Save the species mean file
write.csv(sp.means, row.names=FALSE,
          file="Data/Tiliquini_Morphology_spMEANS.csv")

############################################################################

# Drop some variables before we calculate our log-shape ratios (geometric mean of size)
sp.means <- dplyr::select(sp.means, -Snout_Vent, -Snout_Axilla, 
                          -Head_Length, -Arm, -Leg)
all.raw.traits <- sp.means

# Compute the geometric mean for obtaining size
allsize <- apply(sp.means[3:21], 1, prod)^(1/ncol(sp.means[3:21]))

# Compute the log shape ratios
allLSR <- sp.means[3:21]/allsize; allLSR$Size <- log(allsize)
# fix the names of the LSR traits
colnames(allLSR)[1:20] <- sapply(colnames(allLSR)[1:20], function(x) gsub("_", "", x))
# make a custom DF for the lsr traits
all.LSR.traits <- allLSR;
rownames(all.LSR.traits) <- sp.means$Genus_species


# Add the taxon, group, and family, names back onto the data frame
rownames(allLSR) <- sp.means$Genus_species;
allLSR$Genus <- sp.means$Genus
allLSR$Genus_species <- sp.means$Genus_species

# Save the log shape ratio file
write.csv(allLSR, row.names=FALSE,
          file="Data/Tiliquini_AllLSR.csv")

############################################################################

# Trim our preferred tree down to match our data
ttree <- read.tree("Trees/Tiliquini_forR.tre")
til.tree <- drop.tip(ttree, setdiff(ttree$tip.label, sp.means$Genus_species))

############################################################################

# Make our LSR data into a list of module data frames
head <- select(all.LSR.traits, HeadWidth, HeadDepth, PosSkull, SnoutEye)
body <- select(all.LSR.traits, Interlimb, BodyWidth, PelvicWidth, PelvicHeight)
tail <- select(all.LSR.traits, TailLength, TailWidth)
limb <- select(all.LSR.traits, UpperArm, LowerArm, Hand, UpperLeg, LowerLeg, Foot)
size <- select(all.LSR.traits, Size)
neck <- select(all.LSR.traits, Neck)
eye  <- select(all.LSR.traits, EyeDiameter)
pelvis <- select(all.LSR.traits, PelvicGap)
# stick it all together
all.modules <- list(head, body, tail, limb, size, neck, eye, pelvis)
names(all.modules) <- c("head", "body", "tail", "limb", "size", "neck", "eye", "pelvis")

############################################################################

# Make an RData object we can easily load
save(til.tree, all.raw.traits, all.LSR.traits, allLSR, all.modules, file="Data/Tiliquini_Data.RData")

############################################################################







