There are three scripts used for analyzing 16S microbiome data associated with this publication.
They should be run in the following order:

1. tve16s.R includes sequence processing and calling OTUs using the DADA2 pipeline. The raw sequence files can be found: 
/projectnb/davies-hb/hannah/TVE_Panama/TVE_16S_ITS/tve_prestress_files/lane1and2_16S/fastqs

2. 16S_DiversityAnalyses.Rmd takes the outputs from the tve16s.R script and calculates microbiome diversity metrics.

3. 16s_CommunityComposition.Rmd takes the outputs from the tve16s.R script and calculates community composition by plotting principal coordinate analyses of all, core, and accessory microbiome datasets.

All relevant input files, and output files that are not produced in the scripts, are found in the data_files folder. 