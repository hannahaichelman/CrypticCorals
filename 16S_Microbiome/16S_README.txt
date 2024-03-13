There are three scripts used for analyzing 16S microbiome data associated with this publication that are located in the 'scripts' folder.
They should be run in the following order:

1. tve16s_AllTimepoints.R includes sequence processing and calling OTUs using the DADA2 pipeline. The raw sequence files can be found on NCBI SRA: X

2. 16S_DiversityAnalyses.Rmd takes the outputs from the tve16s_AllTimepoints.R script and calculates microbiome diversity metrics.

3. 16s_CommunityComposition.Rmd takes the outputs from the tve16s_AllTimepoints.R script and calculates community composition by plotting principal coordinate analyses of all, core, and accessory microbiome datasets.

All relevant input files, and output files that are not produced in the scripts, are found in the 'data_files' folder.
