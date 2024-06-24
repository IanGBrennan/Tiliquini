# BEAST

This directory holds the XML file necessary to replicate our combined evidence analysis of Tiliquini skinks. The XML file is a modified version of that used in [Thorn et al. 2023](https://royalsocietypublishing.org/doi/full/10.1098/rspb.2023.0704). Major changes to the file concern (1) the exclusion of continuous morphological data, (2) the replacement of existing multilocus mitonuclear molecular data with codon and rate partitioned concatenated multilocus AHE data, (3) the inclusion of *Tiliqua frangens*. 

\  

If you are interested in rerunning our analyses, it's important that you use BEAST v1.8.4, as specific model elements are not compatible with more recent versions of BEAST (e.g. BEAST v1.10.x) or BEAST2. 

+ *Tiliquini_100mil_GTR_frangens_noCont.xml*: input file for BEAST analysis

+ *Tiliquini_100mil_MEDIANtimeCI.pdf*: maximum clade credibility tree with median ages indicated at nodes and 95% confidence intervals noted by blue/purple bars.

+ *Tiliquini_100mil_subst_trees.pdf*: figure comparing rate variability among branches between molecular and morphological datasets.
