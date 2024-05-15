# Code below incorporates fossil and extant taxa into the dated molecular phylogeny
# Fossil taxa are added based on combined morph/molecular BEAST analyses
# Additional extant taxa (2x Tribolonotus) added based on Title et al. (2024)
#


# Add four fossil taxa based on placement from BEAST analyses
r1 <- read.tree("/Users/ianbrennan/Desktop/AusARG_MCMCTree/ILN/Tiliquini/mcmctree_FULL_C3_R1234/mcmctree_Tiliquini_FULL_C3_ILN_HKY_R1234.nwk")
r1 <- ladderize(r1)
r1$edge.length <- r1$edge.length*100

r2 <- bind.tip(tree = r1,
                  tip.label = "Proegernia_mikebulli", 
                  where = getMRCA(r1, c("Scincidae_Lissolepis_luctuosa_WAMR90226", "Scincidae_Egernia_hosmeri_SAMAR36707")),
                  position = 5,
                  edge.length = (5 + (max(nodeHeights(r1)) - nodeheight(r1, getMRCA(r1, c("Scincidae_Lissolepis_luctuosa_WAMR90226", "Scincidae_Egernia_hosmeri_SAMAR36707"))))) - 25.6)
plot(r2, show.tip.label = F); axisPhylo()

r3 <- bind.tip(tree = r2,
                  tip.label = "Proegernia_palankarinnensis", 
                  where = which(r2$tip.label=="Scincidae_Corucia_zebrata_ABTC50418"),
                  position = 28,
                  edge.length = 2.4)
plot(r3, show.tip.label = F); axisPhylo()

r4 <- bind.tip(tree = r3,
               tip.label = "Tiliqua_frangens", 
               where = which(r3$tip.label=="Scincidae_Tiliqua_rugosa.rugosa_SAMAR18978"),
               position = 6.5,
               edge.length = 4.2)
plot(r4, show.tip.label = F); axisPhylo()

r5 <- bind.tip(tree = r4,
                  tip.label = "Egernia_gillespieae", 
                  where = getMRCA(r4, c("Scincidae_Bellatorias_obiri_AMSR.100018", "Scincidae_Egernia_hosmeri_SAMAR36707")),
                  position = 2.1,
                  edge.length = (2.1 + (max(nodeHeights(r4)) - nodeheight(r4, getMRCA(r4, c("Scincidae_Bellatorias_obiri_AMSR.100018", "Scincidae_Egernia_hosmeri_SAMAR36707"))))) - 14.73)
plot(r5, show.tip.label = F); axisPhylo()

r6 <- bind.tip(tree = r5,
               tip.label = "Tiliqua_pusilla", 
               where = getMRCA(r5, c("Scincidae_Tiliqua_adelaidensis_SAMAR42426", "Scincidae_Tiliqua_nigrolutea_SAMAR33410")),
               position = 1.5,
               edge.length = (1.5 + (max(nodeHeights(r5)) - nodeheight(r5, getMRCA(r5, c("Scincidae_Tiliqua_adelaidensis_SAMAR42426", "Scincidae_Tiliqua_nigrolutea_SAMAR33410"))))) - 15.16)
plot(r6, show.tip.label = F); axisPhylo()

write.tree(r6, "/Users/ianbrennan/Desktop/AusARG_MCMCTree/ILN/Tiliquini/mcmctree_FULL_C3_R1234/mcmctree_Tiliquini_FULL_C3_ILN_HKY_R1234_wFOSSILS.tre")
##########



# Add two additional Tribolonotus spp. based on dates from Title et al. 2024
t1 <- read.tree("/Users/ianbrennan/Documents/GitHub/Tiliquini/Trees/Tiliquini_forR2.tre")
t1 <- ladderize(t1)
t1 <- drop.tip(t1, c("Sphenodon_punctatus", "Eugongylus_rufescens", "Eutropis_longicaudata"))
plot(t1, cex=0.4); axisPhylo()

t2 <- bind.tip(tree = t1,
               tip.label = "Tribolonotus_blanchardi", 
               where = which(t1$tip.label=="Tribolonotus_pseudoponceleti"),
               position = 20,
               edge.length = 20)
plot(t2, cex=0.3); axisPhylo()

t3 <- bind.tip(tree = t2,
               tip.label = "Tribolonotus_ponceleti", 
               where = which(t2$tip.label=="Tribolonotus_pseudoponceleti"),
               position = 7,
               edge.length = 7)
plot(t3, cex=0.3); axisPhylo()

t3 <- ladderize(t3)
plot(t3, cex=0.3); axisPhylo()
write.tree(t3, "/Users/ianbrennan/Documents/GitHub/Tiliquini/Trees/Tiliquini_forR2.tre")

