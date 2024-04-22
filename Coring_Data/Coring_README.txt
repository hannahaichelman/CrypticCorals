There are two scripts used for downloading satellite sea surface temperature (SST) data and analyzing historical growth data from coral cores associated with this publication that are located in the 'scripts' folder.
They should be run in the following order:

1. Panama_HadISST.R constructs the .csv file of annual and summer mean temperature data needed in the second script.

2. HistoricalGrowthAnalysis.R takes the output from the Panama_HadISST.R script and the historical growth data from coral cores (Growth_allRegions_JPcoringdata.csv) to analyze historical growth differences across lineage and relate those changes to temperature.

All relevant input files are found in the 'data_files' folder.
