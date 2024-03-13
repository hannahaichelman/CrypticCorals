# This script analyzes the ITS2 data from the Panama Daily Thermal Variability experiment - both T0 and data collected before the thermal challenge
# Author: Hannah E Aichelman
# hannahaichelman@gmail.com


#### Load packages and read in data ####

#packages
#install.packages("decontam")
library(decontam)
packageVersion("decontam") #‘1.10.0’ - NK's version - 1.16.0 HA's version
library(phyloseq)
library(tidyverse)
library(ggpubr)
library(dplyr)
library(viridis)
#library(microbiomeutilities)

# have gone through multiple iterations of the its2 analysis, different input file formats shown below

# SymPortal ITS2 DIV Analysis
# cleaned file up to remove extraneous info in the header in excel, but original file here:
# /Users/hannahaichelman/Dropbox/BU/TVE/16S_ITS2/ITS_All_Timepoints/191_20220121_03_DBV_20220121T084332.profiles.absolute.abund_and_meta.txt
setwd("~/Dropbox/BU/TVE/TVE_Github/DielTempVariability/ITS2_Symbiodiniaceae/")
its2_all = read.csv("data_files/SymPortal_AllTimes_RawDIVs.csv")
head(its2_all)

#### Clean up input data ####

# Remove samples not included in this dataset
its2_all = its2_all %>%
  dplyr::filter(frag != "Alexa2") %>%
  dplyr::filter(frag != "Alexa1") %>%
  dplyr::filter(frag != "MA2") %>%
  dplyr::filter(frag != "ME1") %>%
  dplyr::filter(frag != "negcontrol1") %>%
  dplyr::filter(frag != "ITS-neg-old") %>%
  dplyr::filter(frag != "negcontrol2") %>%
  dplyr::filter(frag != "ITS-neg-new")

# Remove clades that are not symbiotic
its2_all2 = its2_all %>%
  select(-contains("E1c"))

# Look for and remove DIVs with 0 reads after removing samples
colSums(its2_all2[,-1])

# Remove DIVs that have 0 reads
its2_all3 = its2_all2 %>%
  select(-B2.B2e) %>%
  select(-B1)

str(its2_all3)

its2_all3 = its2_all3 %>%
  mutate_if(is.integer,as.numeric)

# Look for individual samples with 0 reads:
rowSums(its2_all3[, 2:length(its2_all3)]) # I3H2 has 0 reads

# remove individuals with 0 reads and keeping lineage 3 for T0 samples
its2_all4 = its2_all3 %>%
  dplyr::filter(frag!= "I3H2") #have to remove sample I3H2 because 0 reads

head(its2_all4)
str(its2_all4)

# keep this df for div analysis
its2_divs = its2_all4

# read in sample data
# DIVs - this sample info includes all T0 and Pre-Stress Individuals that were run through SymPortal together.
samdf = read.csv("data_files/SampleInfo.csv")
dim(samdf)
head(samdf)

# add identifying data that includes lineage 3
phys_metadata = read.csv("data_files/phys_metadata_its2.csv")
head(phys_metadata)
dim(phys_metadata)

# combine with samdf
samdf2 = left_join(samdf, phys_metadata, by = "frag")
head(samdf2)
dim(samdf2)

# remove the rows for samples that have been removed in data
samdf3 = samdf2 %>%
  dplyr::filter(frag != "Alexa2") %>%
  dplyr::filter(frag != "Alexa1") %>%
  dplyr::filter(frag != "MA2") %>%
  dplyr::filter(frag != "ME1") %>%
  dplyr::filter(frag != "negcontrol1") %>%
  dplyr::filter(frag != "ITS-neg-old") %>%
  dplyr::filter(frag != "negcontrol2") %>%
  dplyr::filter(frag != "ITS-neg-new") %>%
  dplyr::filter(frag != "I3H2") %>% #have to remove sample I3H2 'cause 0 reads
  dplyr::select(-dominant_type)

str(samdf3)
samdf3$treat = as.factor(samdf3$treat)
samdf3$timepoint = as.factor(samdf3$timepoint)
samdf3$sitename = as.factor(samdf3$sitename)
samdf3$gen_site = as.factor(samdf3$gen_site)
samdf3$lineage = as.factor(samdf3$lineage)
samdf3$reef = as.factor(samdf3$reef)
samdf3$Sample_or_Control = as.factor(samdf3$Sample_or_Control)

# change levels of treatment factor
samdf3$treat = factor(samdf3$treat,levels = c("Initial","Control", "Control 2","Low Var","Mod Var","High Var"))

# make T0 data frames to separate initial timepoint analyses out first
samdf_t0 = samdf3 %>%
  dplyr::filter(treat == "Initial")

rownames(samdf_t0) = samdf_t0$frag

its2_divs_t0 = its2_divs %>%
  dplyr::filter(frag %in% samdf_t0$frag) %>%
  column_to_rownames("frag") #rownames have to match between counts table & sample data table or else phyloseq will throw a fit

# make pre-stress data frame
# in addition to the filtering we did for the T0 samples, we need to remove Lineage 3 (L3) individuals
samdf_prestress = samdf3 %>%
  dplyr::filter(timepoint == "pre_stress") %>%
  dplyr::filter(is.na(lineage) | lineage!="L3") # remove lineage 3 from dataset but don't remove other NAs for lineage

rownames(samdf_prestress) = samdf_prestress$frag

its2_divs_prestress = its2_divs %>%
  dplyr::filter(frag!="I2D4") %>% #lineage 3
  dplyr::filter(frag!="I2I7") %>% #lineage 3
  dplyr::filter(frag!="I2I4") %>% #lineage 3
  dplyr::filter(frag!="I2I3") %>% #lineage 3
  dplyr::filter(frag!="I2D1") %>% #lineage 3
  dplyr::filter(frag!="I2H13") %>% #lineage 3
  dplyr::filter(frag!="I2H2") %>% #lineage 3
  dplyr::filter(frag!="I2H12") %>% #lineage 3
  dplyr::filter(frag!="I2D3") %>% #lineage 3
  dplyr::filter(frag!="I2I10") %>% #lineage 3
  dplyr::filter(frag!="I2H4") %>% #lineage 3
  dplyr::filter(frag %in% samdf_prestress$frag) %>%
  column_to_rownames("frag")

# make all sample data frame for looking at shifts in syms
samdf_all = samdf3 %>%
  dplyr::filter(treat != "Control 2")

rownames(samdf_all) = samdf_all$frag

its2_divs_all = its2_divs %>%
  dplyr::filter(frag %in% samdf_all$frag) %>%
  column_to_rownames("frag") #rownames have to match between counts table & sample data table or else phyloseq will throw a fit

## making a taxa table for phyloseq - T0
taxa.t0 <- data.frame(colnames(its2_divs_t0)) #extract sym data

colnames(taxa.t0) <- c("DIV") #changing the column name to be more user-friendly

taxa.t0$majority_its2 = c("A4z",	"A4",	"B19",	"B5",	"B5a",	"C1",	"C1",	"C1",	"C1",	"C1",	"C1",	"C1",	"C3",
                       "C1",	"C3",	"C3af",	"C1",	"C1",	"C1",	"C1",	"C1",	"C1",	"C1",	"C1",	"C15",	"C1",	"C1",
                       "C1",	"D1",	"D1",	"D1",	"D1",	"D1",	"D1",	"D1",	"D1",	"D1")

taxa.t0$genus = str_sub(taxa.t0$DIV, 1, 1)

str(taxa.t0)
taxa.t0$DIV = as.factor(taxa.t0$DIV)
taxa.t0$majority_its2 = as.factor(taxa.t0$majority_its2)
taxa.t0$genus = as.factor(taxa.t0$genus)

#rownmaes also have to match between the columns of the counts table & the taxa table or ps freaks out
rownames(taxa.t0) <- taxa.t0$DIV

taxa.t0.m <- as.matrix(taxa.t0) #also has to be a matrix

ps.its2.t0 <- phyloseq(sample_data(samdf_t0),
         otu_table(its2_divs_t0,taxa_are_rows=FALSE),
         tax_table(taxa.t0.m))
ps.its2.t0
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 37 taxa and 54 samples ]
# sample_data() Sample Data:       [ 54 samples by 8 sample variables ]
# tax_table()   Taxonomy Table:    [ 37 taxa by 3 taxonomic ranks ]


## making a taxa table for phyloseq - prestress samples
taxa.ps <- data.frame(colnames(its2_divs_prestress)) #extract sym data

colnames(taxa.ps) <- c("DIV") #changing the column name to be more user-friendly

taxa.ps$majority_its2 = c("A4z",	"A4",	"B19",	"B5",	"B5a",	"C1",	"C1",	"C1",	"C1",	"C1",	"C1",	"C1",	"C3",
                          "C1",	"C3",	"C3af",	"C1",	"C1",	"C1",	"C1",	"C1",	"C1",	"C1",	"C1",	"C15",	"C1",	"C1",
                          "C1",	"D1",	"D1",	"D1",	"D1",	"D1",	"D1",	"D1",	"D1",	"D1")

taxa.ps$genus = str_sub(taxa.ps$DIV, 1, 1)

str(taxa.ps)
taxa.ps$DIV = as.factor(taxa.ps$DIV)
taxa.ps$majority_its2 = as.factor(taxa.ps$majority_its2)
taxa.ps$genus = as.factor(taxa.ps$genus)

#rownmaes also have to match between the columns of the counts table & the taxa table or ps freaks out
rownames(taxa.ps) <- taxa.ps$DIV

taxa.ps.m <- as.matrix(taxa.ps) #also has to be a matrix

ps.its2.ps <- phyloseq(sample_data(samdf_prestress),
                    otu_table(its2_divs_prestress,taxa_are_rows=FALSE),
                    tax_table(taxa.ps.m))
ps.its2.ps
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 37 taxa and 172 samples ]
# sample_data() Sample Data:       [ 172 samples by 8 sample variables ]
# tax_table()   Taxonomy Table:    [ 37 taxa by 3 taxonomic ranks ]

## making a taxa table for phyloseq - all samples
taxa.all <- data.frame(colnames(its2_divs_all)) #extract sym data

colnames(taxa.all) <- c("DIV") #changing the column name to be more user-friendly

taxa.all$majority_its2 = c("A4z",	"A4",	"B19",	"B5",	"B5a",	"C1",	"C1",	"C1",	"C1",	"C1",	"C1",	"C1",	"C3",
                          "C1",	"C3",	"C3af",	"C1",	"C1",	"C1",	"C1",	"C1",	"C1",	"C1",	"C1",	"C15",	"C1",	"C1",
                          "C1",	"D1",	"D1",	"D1",	"D1",	"D1",	"D1",	"D1",	"D1",	"D1")

taxa.all$genus = str_sub(taxa.all$DIV, 1, 1)

str(taxa.all)
taxa.all$DIV = as.factor(taxa.all$DIV)
taxa.all$majority_its2 = as.factor(taxa.all$majority_its2)
taxa.all$genus = as.factor(taxa.all$genus)

#rownmaes also have to match between the columns of the counts table & the taxa table or ps freaks out
rownames(taxa.all) <- taxa.all$DIV

taxa.all.m <- as.matrix(taxa.all) #also has to be a matrix

ps.its2.all <- phyloseq(sample_data(samdf_all),
                       otu_table(its2_divs_all,taxa_are_rows=FALSE),
                       tax_table(taxa.all.m))
ps.its2.all
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 37 taxa and 231 samples ]
# sample_data() Sample Data:       [ 231 samples by 8 sample variables ]
# tax_table()   Taxonomy Table:    [ 37 taxa by 3 taxonomic ranks ]
#


#### DECONTAM ####
# Skipping unless doing pre-med sequence analysis

ps.rel <- transform_sample_counts(ps.its2, function(OTU) OTU/sum(OTU))

#keeps samples with summed counts greater than 0
ps.its2.no0 <- prune_samples(sample_sums(ps.its2)!=0, ps.its2)
ps.its2.no0 # don't lose any for pre-med seqs

# now use package decontam to remove contaminant sequences
#creates a TRUE/FALSE column for whether a sample is a negative control or not
sample_data(ps.its2.no0)$control <- sample_data(ps.its2.no0)$Sample_or_Control == "Control"

contamdf.prev <- isContaminant(ps.its2.no0, neg="control",threshold=0.5)
table(contamdf.prev$contaminant)
# FALSE  TRUE  - pre med seqs
# 440    31

# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(ps.its2.no0, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Control", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_Control != "Control", ps.pa)

# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

#remove contaminants from ps object:
ps.clean1 <- prune_taxa(!contamdf.prev$contaminant,ps.its2.no0)

#also remove negative controls, don't need them anymore I think
ps.clean <- subset_samples(ps.clean1,(Sample_or_Control!="Control"))
ps.clean #31 contaminants removed, 3 samples removed (all from pre stress timepoint)
# POST-MED SEQS:
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 440 taxa and 172 samples ]
# sample_data() Sample Data:       [ 172 samples by 9 sample variables ]
# tax_table()   Taxonomy Table:    [ 440 taxa by 7 taxonomic ranks ]

ps.clean.rel <- transform_sample_counts(ps.clean, function(OTU) OTU/sum(OTU))


#### Bar plot - raw table pre-processing ####
# Skipping decontam for now since we are ignoring negative controls per ben hume
ps.rel.t0 <- transform_sample_counts(ps.its2.t0, function(OTU) OTU/sum(OTU))
plot_bar(ps.rel.t0, x="gen_site",fill="majority_its2")+
  theme_classic()

ps.rel.ps <- transform_sample_counts(ps.its2.ps, function(OTU) OTU/sum(OTU))
plot_bar(ps.rel.ps, x="frag",fill="majority_its2")+
  theme_classic()

ps.rel.all <- transform_sample_counts(ps.its2.all, function(OTU) OTU/sum(OTU))
plot_bar(ps.rel.all, x="frag",fill="majority_its2")+
  theme_classic()

#keeps samples with summed counts greater than 0
ps.its2.t0.no0 <- prune_samples(sample_sums(ps.its2.t0)!=0, ps.its2.t0)
ps.its2.t0.no0 # don't lose any

ps.its2.ps.no0 <- prune_samples(sample_sums(ps.its2.ps)!=0, ps.its2.ps)
ps.its2.ps.no0 # don't lose any

# Remove NA's to create ps object we will plot
# first for T0
ps.cleanest.t0 <- subset_samples(ps.its2.t0.no0,(!is.na(lineage)))
ps.cleanest.rel.t0 <- transform_sample_counts(ps.cleanest.t0, function(OTU) OTU/sum(OTU))

# save phyloseq object
saveRDS(ps.cleanest.t0, "data_files/ps.its2.t0.RDS")
saveRDS(ps.cleanest.rel.t0, "data_files/ps.its2.t0.rel.RDS")

# write this out as dataframe
seqtab.rel.t0 <- data.frame(ps.cleanest.rel.t0@otu_table)
write.csv(seqtab.rel.t0, file="data_files/ITS2.t0.seqtab.rel.csv")
sample.data.t0 = data.frame(sample_data(ps.cleanest.rel.t0))
write.csv(sample.data.t0, file="data_files/ITS2.t0.sample.data.rel.csv")

taxa.t0 <- data.frame(ps.cleanest.rel.t0@tax_table)
mtaxa.t0 <- as.matrix(taxa.t0)
write.csv(taxa.t0, file = "data_files/symportal_taxa.t0.csv")


# then for Pre-Stress - removing low variability and high variability to simplify experimental design
ps.cleaner.ps <- subset_samples(ps.its2.ps.no0,(!is.na(lineage)))
ps.cleaner2.ps <- subset_samples(ps.cleaner.ps,(!is.na(treat)))
ps.cleaner3.ps <- subset_samples(ps.cleaner2.ps,(treat!="Control 2")) #156 samples remain
ps.cleaner4.ps <- subset_samples(ps.cleaner3.ps,(treat!="Low Var")) #122 samples remain
ps.cleanest.ps <- subset_samples(ps.cleaner4.ps,(treat!="High Var")) #84 samples remain

ps.cleanest.rel.ps <- transform_sample_counts(ps.cleanest.ps, function(OTU) OTU/sum(OTU))

# save phyloseq object
saveRDS(ps.cleanest.ps, "data_files/ps.its2.prestress.RDS")
saveRDS(ps.cleanest.rel.ps, "data_files/ps.its2.prestress.rel.RDS")

# write this out as dataframe
seqtab.rel.ps <- data.frame(ps.cleanest.rel.ps@otu_table)
write.csv(seqtab.rel.ps, file="data_files/ITS2.prestress.seqtab.rel.csv")
sample.data.ps = data.frame(sample_data(ps.cleanest.rel.ps))
write.csv(sample.data.ps, file="data_files/ITS2.prestress.sample.data.rel.csv")

taxa.ps <- data.frame(ps.cleanest.rel.ps@tax_table)
mtaxa.ps <- as.matrix(taxa.ps)
write.csv(taxa.ps, file = "data_files/symportal_taxa.prestress.csv")

#### Compare time points ####
p.all.gen_site = plot_bar(ps.rel.all, x="frag", fill="majority_its2") +
  geom_col(position = position_dodge2(width = 0.8, preserve = "single")) +
  scale_fill_manual(name = "Genus", values = its2_cols_greens) +
  theme_bw() +
  facet_wrap(~gen_site, ncol = 6, nrow = 9, scales = "free_x") +
  #theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.title.x = element_blank(), axis.title.y = element_blank())
p.all.gen_site
ggsave(p.all.gen_site, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/16S_ITS2/ITS_All_Timepoints/syms_compare_timepoints.pdf", width=12, height=12, units=c("in"), useDingbats=FALSE)

#### Normalized Results ####
# not using for now, included as a test to see if it changed PCAs much

#BiocManager::install("edgeR", update = FALSE)
library(edgeR)

seqs.types <- as.data.frame(ps.cleanest@otu_table)
seqs.types.t <- t(seqs.types)

its2SeqList = DGEList(counts = seqs.types.t)
head(its2SeqList$samples)

its2SeqNorm =  calcNormFactors(its2SeqList, method = "TMM")
head(its2SeqNorm$samples)

its2TMM = t(cpm(its2SeqNorm, normalized.lib.sizes = TRUE))

#phyloseq
ps.norm <- phyloseq(otu_table(its2TMM, taxa_are_rows=FALSE),
                    sample_data(sample.data),
                    tax_table(mtaxa))
ps.norm
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 440 taxa and 156 samples ]
# sample_data() Sample Data:       [ 176 samples by 9 sample variables ]
# tax_table()   Taxonomy Table:    [ 440 taxa by 7 taxonomic ranks ]
saveRDS(ps.norm, "data_files/ps.its2.norm.RDS")


#### PS Object Versions ####
# T0
ps.cleanest.t0 = readRDS("data_files/ps.its2.t0.RDS")
seqtab.t0 <- data.frame(ps.cleanest.t0@otu_table)
samdf.t0 <- data.frame(ps.cleanest.t0@sam_data)

ps.cleanest.t0.rel = readRDS("data_files/ps.its2.t0.rel.RDS")
seqtab.rel.t0 <- data.frame(ps.cleanest.t0.rel@otu_table)
samdf.rel.t0 <- data.frame(ps.cleanest.t0.rel@sam_data)

taxa.t0 = read.csv(file = "data_files/symportal_taxa.t0.csv", header = TRUE) %>%
  select(-X)
rownames(taxa.t0) <- as.factor(taxa.t0$DIV)
sum(taxa.t0$genus == "B") # 3 div
sum(taxa.t0$genus == "C") # 23 div
sum(taxa.t0$genus == "D") # 9 div

mtaxa.t0 <- as.matrix(taxa.t0)

# Pre-Stress
ps.cleanest.ps = readRDS("data_files/ps.its2.prestress.RDS")
seqtab.ps <- data.frame(ps.cleanest.ps@otu_table)
samdf.ps <- data.frame(ps.cleanest.ps@sam_data)

ps.cleanest.ps.rel = readRDS("data_files/ps.its2.prestress.rel.RDS")
seqtab.rel.ps <- data.frame(ps.cleanest.ps.rel@otu_table)
samdf.rel.ps <- data.frame(ps.cleanest.ps.rel@sam_data)

taxa.ps = read.csv(file = "data_files/symportal_taxa.prestress.csv", header = TRUE) %>%
  select(-X)
rownames(taxa.ps) <- as.factor(taxa.ps$DIV)
sum(taxa.ps$genus == "B") # 3 div
sum(taxa.ps$genus == "C") # 23 div
sum(taxa.ps$genus == "D") # 9 div

mtaxa.ps <- as.matrix(taxa.ps)

#### Bar plot - post-processing T0 ####
its2_cols_greens = c("A4" = "#ffeda0", "A4z" =  "#fd8d3c",
                     "B19" = "#4eb3d3", "B5a" = "#0868ac", "B5" = "#800026",
                     "C1" = "#edf8e9", "C15" = "#feb24c",  "C3" = "#a1d99b",
                     "C3af" = "#238b45", "D1" = "#00441b")

#its2_cols_blues = c("#f6eff7", "#d0d1e6", "#a6bddb", "#67a9cf", "#3690c0", "#02818a", "#016450")

# bar plot of all individuals
p.all = plot_bar(ps.cleanest.t0.rel, x="gen_site", fill="majority_its2") +
  theme_bw() +
  scale_fill_manual(name = "Majority ITS2", values = its2_cols_greens) +
  xlab("Sample") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
p.all
ggsave(p.all, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/16S_ITS2/ITS_All_Timepoints/syms_t0.pdf", width=11, height=4, units=c("in"), useDingbats=FALSE)

# Now facet wrap/aggregate by interesting factors

# Lineage:
ps.lin <- merge_samples(ps.cleanest.t0, "lineage")
ps.rel.lin <- transform_sample_counts(ps.lin, function(x) x / sum(x))

p.lineage = plot_bar(ps.rel.lin, fill="majority_its2") +
  geom_bar(stat = "identity") +
  scale_fill_manual(name = "Majority ITS2", values = its2_cols_greens) +
  theme_bw()
p.lineage
ggsave(p.lineage, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/16S_ITS2/ITS_All_Timepoints/syms_t0_lineage.pdf", width=4, height=3.5, units=c("in"), useDingbats=FALSE)


# Lineage only:
p.lin = plot_bar(ps.cleanest.t0.rel, x="gen_site", fill="majority_its2") +
  #geom_bar(stat = "identity", width = 0.2) +
  geom_col(position = position_dodge2(width = 0.8, preserve = "single")) +
  scale_fill_manual(name = "Genus", values = its2_cols_greens) +
  theme_bw() +
  facet_grid(~lineage, scales = "free", space = "free") +
  #theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.title.x = element_blank(), axis.title.y = element_blank())
p.lin
ggsave(p.lin, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/16S_ITS2/ITS_All_Timepoints/syms_lineagefacet_t0_majorityITS2.pdf", width=10, height=4, units=c("in"), useDingbats=FALSE)

# Site of Origin:
p.site = plot_bar(ps.cleanest.t0.rel, x="gen_site", fill="majority_its2") +
  #geom_bar(stat = "identity") +
  geom_col(position = position_dodge2(width = 0.8, preserve = "single")) +
  scale_fill_manual(values = its2_cols_greens) +
  theme_bw() +
  facet_grid(~sitename, scales = "free", space = "free") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.title.x = element_blank(), axis.title.y = element_blank())
p.site
ggsave(p.site, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/16S_ITS2/ITS_All_Timepoints/syms_sitefacet_t0_majorityITS2.pdf", width=10, height=4, units=c("in"), useDingbats=FALSE)

#### Bar plot - post-processing PreStress ####
#its2_cols_end = c("#D53E4F", "#F46D43", "#3288BD", "#FEE08B", "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD", "#FDAE61")
#its2_cols_purple = c("white", "#efedf5","#bcbddc", "#756bb1")
its2_cols_greens = c("A4" = "#ffeda0", "A4z" =  "#fd8d3c",
                     "B19" = "#4eb3d3", "B5a" = "#0868ac", "B5" = "#800026",
                     "C1" = "#edf8e9", "C15" = "#feb24c",  "C3" = "#a1d99b",
                     "C3af" = "#238b45", "D1" = "#00441b")

#its2_cols_blues = c("#f6eff7", "#d0d1e6", "#a6bddb", "#67a9cf", "#3690c0", "#02818a", "#016450")

# bar plot of all individuals
p.all = plot_bar(ps.cleanest.ps.rel, x="frag", fill="majority_its2") +
  theme_bw() +
  scale_fill_manual(name = "Majority ITS2", values = its2_cols_greens) +
  xlab("Sample") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
p.all
ggsave(p.all, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/16S_ITS2/ITS_All_Timepoints/syms_prestress.pdf", width=11, height=4, units=c("in"), useDingbats=FALSE)

# subset by factors we need:
ps.cleanest.rel.control = subset_samples(ps.cleanest.ps.rel, treat=="Control")
#ps.cleanest.rel.lowvar = subset_samples(ps.cleanest.rel, treat=="Low Var")
ps.cleanest.rel.modvar = subset_samples(ps.cleanest.ps.rel, treat=="Mod Var")
#ps.cleanest.rel.highvar = subset_samples(ps.cleanest.rel, treat=="High Var")


# Now facet wrap/aggregate by interesting factors

# Lineage:
ps.lin <- merge_samples(ps.cleanest.ps, "lineage")
ps.rel.lin <- transform_sample_counts(ps.lin, function(x) x / sum(x))

p.lineage = plot_bar(ps.rel.lin, fill="majority_its2") +
  geom_bar(stat = "identity") +
  scale_fill_manual(name = "Majority ITS2", values = its2_cols_greens) +
  theme_bw()
p.lineage
ggsave(p.lineage, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/16S_ITS2/ITS_All_Timepoints/syms_prestress_lineage.pdf", width=4, height=3.5, units=c("in"), useDingbats=FALSE)

# Treatment:
p.treat = plot_bar(ps.cleanest.ps.rel, x="frag", fill="majority_its2") +
  geom_col(position = position_dodge2(width = 0.8, preserve = "single")) +
  scale_fill_manual(values = its2_cols_greens) +
  theme_bw() +
  facet_grid(~treat, scales = "free", space = "free") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.title.x = element_blank(), axis.title.y = element_blank())
p.treat
ggsave(p.treat, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/16S_ITS2/ITS_All_Timepoints/syms_prestress_treatment.pdf", width=10, height=4, units=c("in"), useDingbats=FALSE)

# Treatment faceted by lineage
p.treat.control = plot_bar(ps.cleanest.rel.control, x="frag", fill="majority_its2") +
  geom_col(position = position_dodge2(width = 0.8, preserve = "single")) +
  scale_fill_manual(name = "Majority ITS2", values = its2_cols_greens) +
  theme_bw() +
  facet_grid(~lineage, scales = "free", space = "free") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.title.x = element_blank())+
  coord_cartesian(expand=FALSE)
p.treat.control

p.treat.mod = plot_bar(ps.cleanest.rel.modvar, x="frag", fill="majority_its2") +
  geom_col(position = position_dodge2(width = 0.8, preserve = "single")) +
  scale_fill_manual(name = "Majority ITS2", values = its2_cols_greens) +
  theme_bw() +
  facet_grid(~lineage, scales = "free", space = "free") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.title.x = element_blank())+
  coord_cartesian(expand=FALSE)
p.treat.mod

p.alltreat = ggarrange(p.treat.control, p.treat.mod,
                       legend = "right", ncol = 2, nrow = 1,
                       label.x = FALSE, common.legend = TRUE)
ggsave(p.alltreat, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/16S_ITS2/ITS_All_Timepoints/syms_treatment_lineagefacet_majorityITS2.pdf", width=8, height=4, units=c("in"), useDingbats=FALSE)

# Lineage only:
p.lin = plot_bar(ps.cleanest.ps.rel, x="frag", fill="majority_its2") +
  geom_col(position = position_dodge2(width = 0.8, preserve = "single")) +
  scale_fill_manual(name = "Genus", values = its2_cols_greens) +
  theme_bw() +
  facet_grid(~lineage, scales = "free", space = "free") +
  #theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.title.x = element_blank(), axis.title.y = element_blank())
p.lin
ggsave(p.lin, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/16S_ITS2/ITS_All_Timepoints/syms_lineagefacet_prestress_majorityITS2.pdf", width=10, height=4, units=c("in"), useDingbats=FALSE)

# Site of Origin:
p.site = plot_bar(ps.cleanest.ps.rel, x="frag", fill="majority_its2") +
  geom_col(position = position_dodge2(width = 0.8, preserve = "single")) +
  scale_fill_manual(values = its2_cols_greens) +
  theme_bw() +
  facet_grid(~sitename, scales = "free", space = "free") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.title.x = element_blank(), axis.title.y = element_blank())
p.site
ggsave(p.site, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/16S_ITS2/ITS_All_Timepoints/syms_sitefacet_prestress_majorityITS2.pdf", width=10, height=4, units=c("in"), useDingbats=FALSE)

#### Bray Curtis Dissimilarity PCoAs - T0 ####
# color schemes
cols_site_diverging <- c("CI" = "#543005", "PD"= "#bf812d",  "SP"= "#dfc27d",  "BN" = "#003c30", "BS"= "#35978f", "CA"= "#80cdc1")

cols_treat_reds <- c("darkgrey", "#FF9966","#CC3300","#7f0000")

cols_lineage <- c("L1" = "#3f007d", "L2" = "#807dba", "L3" = "#bcbddc")

ps.ord.t0 <- ordinate(ps.cleanest.t0.rel,"PCoA",distance="bray")

pcoa.site = plot_ordination(ps.cleanest.t0.rel, ps.ord.t0, color ="sitename", shape = "sitename")+
  geom_point(alpha=0.8)+
  scale_color_manual(name = "Site", breaks = c("BN","BS","CA","CI","PD","SP"), values=cols_site_diverging)+
  scale_shape_manual(name = "Site", breaks = c("BN","BS","CA","CI","PD","SP"), values=c(15,16,17,22,21,24))+
  stat_ellipse() +
  #xlim(-1.3,1) +
  #ylim(-1.2,1) +
  theme_bw()
pcoa.site

pcoa.lineage = plot_ordination(ps.cleanest.rel, ps.ord, color ="lineage")+
  geom_point(alpha=0.8)+
  scale_color_manual(name = "Lineage", values=cols_lineage)+
  stat_ellipse()+
  #xlim(-1.3,1) +
  #ylim(-1.2,1) +
  theme_bw()
pcoa.lineage

gg.pcoa <- ggarrange(pcoa.site,pcoa.lineage,labels=c("A.","B."),nrow=1,common.legend=F,legend="none")
ggsave(gg.pcoa, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/16S_ITS2/ITS_All_Timepoints/bray_pcoas_rel_divs_t0.pdf", width=6, height=3, units=c("in"), useDingbats=FALSE)

#### Bray Curtis Dissimilarity PCoAs - PreStress ####
# color schemes
cols_site_diverging <- c("CI" = "#543005", "PD"= "#bf812d",  "SP"= "#dfc27d",  "BN" = "#003c30", "BS"= "#35978f", "CA"= "#80cdc1")

#cols_treat_reds <- c("darkgrey", "#FF9966","#CC3300","#7f0000")
cols_treat_reds <- c("darkgrey","#CC3300")

cols_lineage <- c("L1" = "#3f007d", "L2" = "#807dba", "L3" = "#bcbddc")

ps.ord.ps <- ordinate(ps.cleanest.ps.rel,"PCoA",distance="bray")

pcoa.site = plot_ordination(ps.cleanest.ps.rel, ps.ord.ps, color ="sitename", shape = "sitename")+
  geom_point(alpha=0.8)+
  scale_color_manual(name = "Site", breaks = c("BN","BS","CA","CI","PD","SP"), values=cols_site_diverging)+
  scale_shape_manual(name = "Site", breaks = c("BN","BS","CA","CI","PD","SP"), values=c(15,16,17,22,21,24))+
  stat_ellipse() +
  #xlim(-1.5,1.0) +
  #ylim(-1.2,1) +
  theme_bw()
pcoa.site

pcoa.treat = plot_ordination(ps.cleanest.ps.rel, ps.ord.ps, color ="treat", shape = "treat")+
  geom_point(alpha=0.8)+
  scale_color_manual(name = "DTV Treatment", values=cols_treat_reds)+
  scale_shape_manual(name = "DTV Treatment", values=c(15,19,17,18))+
  stat_ellipse()+
  #xlim(-1.3,1) +
  #ylim(-1.2,1) +
  theme_bw()
pcoa.treat
#ggsave(pcoa.treat, filename = "/Users/hannahaichelman/Documents/BU/TVE/16S_ITS2/ITS_PreStress_Timepoint/bray_pcoa_treat.pdf", width=4.5, height=3, units=c("in"), useDingbats=FALSE)

pcoa.lineage = plot_ordination(ps.cleanest.ps.rel, ps.ord.ps, color ="lineage")+
  geom_point(alpha=0.8)+
  scale_color_manual(name = "Lineage", values=cols_lineage)+
  stat_ellipse()+
  #xlim(-1.3,1) +
  #ylim(-1.2,1) +
  theme_bw()
pcoa.lineage

gg.pcoa <- ggarrange(pcoa.site,pcoa.treat,pcoa.lineage,labels=c("A.","B.","C."),nrow=1,common.legend=F,legend="none")
ggsave(gg.pcoa, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/16S_ITS2/ITS_All_Timepoints/bray_pcoas_rel_divs_prestress.pdf", width=6, height=3, units=c("in"), useDingbats=FALSE)


## try adding loadings for top 5 sequences
# ps2 <- tax_glom(ps.cleanest.rel, "Species")
# ps2f <- format_to_besthit(ps2)
# ps.ord.2 <- ordinate(ps2f, method = "PCoA",distance="bray")
#
# p <- plot_ordination_utils(ps2f, ps.ord.2,
#                            color = "lineage", plot.arrow = TRUE,
#                            scale.arrow = 1, top.taxa = 5) +
#   scale_color_manual(name = "Lineage", values=cols_lineage_purples)
#   facet_wrap(~sitename)
# p
# ggsave(p, filename = "/Users/hannahaichelman/Documents/BU/TVE/16S_ITS2/ITS_PreStress_Timepoint/bray_pcoas_site_arrows.pdf", width=4, height=3, units=c("in"), useDingbats=FALSE)
#
# pcoa.treat = plot_ordination_utils(ps2f, ps.ord.2, color ="lineage", plot.arrow = TRUE, top.taxa = 5) +
#   geom_point(alpha=0.8)+
#   scale_color_manual(name = "DTV Treatment", values=cols_treat_reds)#+
#   scale_shape_manual(name = "DTV Treatment", values=c(15,19,17,18))+
#   stat_ellipse()+
#   xlim(-1.05,1.1) +
#   ylim(-0.8,1.3) +
#   theme_bw()
# ggsave(pcoa.treat, filename = "/Users/hannahaichelman/Documents/BU/TVE/16S_ITS2/ITS_PreStress_Timepoint/bray_pcoa_treat.pdf", width=4, height=3, units=c("in"), useDingbats=FALSE)

#### PCoA Stats - T0 ####
library(vegan)
#remotes::install_github("Jtrachsel/funfuns")
library(funfuns)
library(dplyr)
library(edgeR)

## Raw (cleaned)
ps = readRDS("data_files/ps.its2.t0.RDS")

seq.ps <- data.frame(ps@otu_table)
samdf.ps <- data.frame(ps@sam_data)
dist.ps <- vegdist(seq.ps)

# dispersion
# by sitename
bet.ps <- betadisper(dist.ps,samdf.ps$sitename)
anova(bet.ps)
# Response: Distances
#           Df  Sum Sq  Mean Sq F value Pr(>F)
# Groups     5 0.60392 0.120783  1.7229 0.1493
# Residuals 44 3.08460 0.070105
permutest(bet.ps,pairwise=TRUE,permutations=999)
plot(bet.ps)

# by lineage
bet.ps <- betadisper(dist.ps,samdf.ps$lineage)
anova(bet.ps)
# Response: Distances
#           Df Sum Sq  Mean Sq F value Pr(>F)
# Groups     2 0.0090 0.004502  0.0652  0.937
# Residuals 47 3.2473 0.069092
permutest(bet.ps,pairwise=TRUE,permutations=999)
plot(bet.ps)

# adonis
adonis2(formula = seq.ps ~ sitename + lineage, data = samdf.ps, permutations = 999)
#           Df SumOfSqs      R2      F Pr(>F)
# sitename  5   5.3338 0.26425 3.2058  0.001 ***
# lineage   2   0.8752 0.04336 1.3150  0.232
# Residual 42  13.9758 0.69239
# Total    49  20.1847 1.00000

## Relative abundance
# Report this since it is the data used for making the PCA's
ps.cleanest.rel = readRDS("data_files/ps.its2.t0.rel.RDS")

seq.ps <- data.frame(ps.cleanest.rel@otu_table)
samdf.ps <- data.frame(ps.cleanest.rel@sam_data)
dist.ps <- vegdist(seq.ps)

# dispersion by sitename
bet.ps <- betadisper(dist.ps,samdf.ps$sitename)
anova(bet.ps)
# Response: Distances
#           Df Sum Sq Mean Sq F value Pr(>F)
# Groups     5 1.0807 0.21613  1.4332 0.2313
# Residuals 44 6.6355 0.15081
permutest(bet.ps,pairwise=TRUE,permutations=999)
# Pairwise comparisons:
#   (Observed p-value below diagonal, permuted p-value above diagonal)
# BN         BS         CA         CI         PD    SP
# BN            7.0200e-01 7.3000e-02 6.4500e-01 9.9200e-01 0.808
# BS 6.9232e-01            1.0000e-03 9.0000e-01 6.7800e-01 0.488
# CA 7.4554e-02 3.0221e-04            1.0000e-03 9.6000e-02 0.143
# CI 6.4652e-01 9.0342e-01 6.4303e-06            6.4500e-01 0.450
# PD 9.9026e-01 6.9351e-01 9.1433e-02 6.5199e-01            0.827
# SP 8.1815e-01 4.9378e-01 1.4203e-01 4.5997e-01 8.3228e-01

plot(bet.ps)

# by lineage
bet.ps <- betadisper(dist.ps,samdf.ps$lineage)
anova(bet.ps)
# Response: Distances
#           Df Sum Sq  Mean Sq F value Pr(>F)
# Groups     2 0.0213 0.010628  0.0604 0.9414
# Residuals 47 8.2646 0.175843
permutest(bet.ps,pairwise=TRUE,permutations=999) # no sig diffs
plot(bet.ps)

# adonis
adonis2(formula = seq.ps ~ sitename + lineage, data = samdf.ps, permutations = 999)
#           Df SumOfSqs      R2      F Pr(>F)
# sitename  5   5.5299 0.28089 3.4964  0.001 ***
# lineage   2   0.8716 0.04427 1.3776  0.185
# Residual 42  13.2856 0.67484
# Total    49  19.6871 1.00000

#pairwise.adonis(seq.ps, factors=samdf.ps$lineage, permutations=999) #no significant comparisons
#pairwise.adonis(seq.ps, factors=samdf.ps$sitename, permutations=999) #all comparisons different except PD vs SP and BN vs CA

#### PCoA Stats - PreStress ####
library(vegan)
#remotes::install_github("Jtrachsel/funfuns")
library(funfuns)
library(dplyr)
library(edgeR)

## Raw (cleaned)
ps = readRDS("data_files/ps.its2.prestress.RDS")

seq.ps <- data.frame(ps@otu_table)
samdf.ps <- data.frame(ps@sam_data)
dist.ps <- vegdist(seq.ps)

# dispersion
# by sitename
bet.ps <- betadisper(dist.ps,samdf.ps$sitename)
anova(bet.ps)
# Response: Distances
#           Df  Sum Sq  Mean Sq F value Pr(>F)
# Groups     5 0.60392 0.120783  1.7229 0.1493
# Residuals 44 3.08460 0.070105
permutest(bet.ps,pairwise=TRUE,permutations=999)
plot(bet.ps)

# by lineage
bet.ps <- betadisper(dist.ps,samdf.ps$lineage)
anova(bet.ps)
# Response: Distances
#           Df Sum Sq  Mean Sq F value Pr(>F)
# Groups     2 0.0090 0.004502  0.0652  0.937
# Residuals 47 3.2473 0.069092
permutest(bet.ps,pairwise=TRUE,permutations=999)
plot(bet.ps)

# adonis
adonis2(formula = seq.ps ~ sitename + lineage, data = samdf.ps, permutations = 999)
#           Df SumOfSqs      R2      F Pr(>F)
# sitename  5   5.3338 0.26425 3.2058  0.001 ***
# lineage   2   0.8752 0.04336 1.3150  0.232
# Residual 42  13.9758 0.69239
# Total    49  20.1847 1.00000

## Relative abundance
# Report this since it is the data used for making the PCA's
ps.cleanest.rel = readRDS("data_files/ps.its2.prestress.rel.RDS")

seq.ps <- data.frame(ps.cleanest.rel@otu_table)
samdf.ps <- data.frame(ps.cleanest.rel@sam_data)
dist.ps <- vegdist(seq.ps)

# dispersion by sitename
bet.ps <- betadisper(dist.ps,samdf.ps$sitename)
anova(bet.ps)
# Response: Distances
#           Df Sum Sq Mean Sq F value Pr(>F)
# Groups     5 1.0807 0.21613  1.4332 0.2313
# Residuals 44 6.6355 0.15081
permutest(bet.ps,pairwise=TRUE,permutations=999)
# Pairwise comparisons:
#   (Observed p-value below diagonal, permuted p-value above diagonal)
# BN         BS         CA         CI         PD    SP
# BN            7.0200e-01 7.3000e-02 6.4500e-01 9.9200e-01 0.808
# BS 6.9232e-01            1.0000e-03 9.0000e-01 6.7800e-01 0.488
# CA 7.4554e-02 3.0221e-04            1.0000e-03 9.6000e-02 0.143
# CI 6.4652e-01 9.0342e-01 6.4303e-06            6.4500e-01 0.450
# PD 9.9026e-01 6.9351e-01 9.1433e-02 6.5199e-01            0.827
# SP 8.1815e-01 4.9378e-01 1.4203e-01 4.5997e-01 8.3228e-01

plot(bet.ps)

# by lineage
bet.ps <- betadisper(dist.ps,samdf.ps$lineage)
anova(bet.ps)
# Response: Distances
#           Df Sum Sq  Mean Sq F value Pr(>F)
# Groups     2 0.0213 0.010628  0.0604 0.9414
# Residuals 47 8.2646 0.175843
permutest(bet.ps,pairwise=TRUE,permutations=999) # no sig diffs
plot(bet.ps)

# adonis
adonis2(formula = seq.ps ~ sitename + lineage, data = samdf.ps, permutations = 999)
#           Df SumOfSqs      R2      F Pr(>F)
# sitename  5   5.5299 0.28089 3.4964  0.001 ***
# lineage   2   0.8716 0.04427 1.3776  0.185
# Residual 42  13.2856 0.67484
# Total    49  19.6871 1.00000

#pairwise.adonis(seq.ps, factors=samdf.ps$lineage, permutations=999) #no significant comparisons
#pairwise.adonis(seq.ps, factors=samdf.ps$sitename, permutations=999) #all comparisons different except PD vs SP and BN vs CA

#### Summary of Dominant ITS2 Majority Types and DIVs ####

ps.cleanest.t0.rel = readRDS("data_files/ps.its2.t0.rel.RDS")
seqtab.rel.t0 <- data.frame(ps.cleanest.t0.rel@otu_table)
samdf.rel.t0 <- data.frame(ps.cleanest.t0.rel@sam_data)

ps.cleanest.ps.rel = readRDS("data_files/ps.its2.prestress.rel.RDS")
seqtab.rel.ps <- data.frame(ps.cleanest.ps.rel@otu_table)
samdf.rel.ps <- data.frame(ps.cleanest.ps.rel@sam_data)

## First T0
# change column names to majority its2 sequence
taxa.t0 = read.csv(file = "data_files/symportal_taxa.t0.csv", header = TRUE) %>%
  select(-X)

taxa.t0$majority_its2
colnames(seqtab.rel.t0) = c("A4z","A4","B19","B5","B5a","C1.1","C1.2","C1.3","C1.4","C1.5","C1.6","C1.7","C3.1","C1.8","C3.2","C3af",
                            "C1.9","C1.10","C1.11","C1.12","C1.13","C1.14","C1.15","C1.16","C15","C1.17","C1.18","C1.19",
                            "D1.1","D1.2","D1.3","D1.4","D1.5","D1.6","D1.7","D1.8","D1.9" )

# make new data frame and sum columns with the same majority its2 sequence
seqtab.rel.t0.new = seqtab.rel.t0 %>%
  mutate(A4z_sum = A4z) %>%
  mutate(A4_sum = A4) %>%
  mutate(B19_sum = B19) %>%
  mutate(B5_sum = B5) %>%
  mutate(B5a_sum = B5a) %>%
  mutate(C1_sum = rowSums(select(., starts_with("C1.")))) %>%
  mutate(C3_sum = rowSums(select(., starts_with("C3.")))) %>%
  mutate(C3af_sum = C3af) %>%
  mutate(C15_sum = C15) %>%
  mutate(D1_sum = rowSums(select(., starts_with("D1.")))) %>%
  rownames_to_column(var = "frag") %>%
  select(frag, contains("_sum"))

its2.rel.t0.combined = left_join(seqtab.rel.t0.new, samdf.rel.t0, by = "frag")

its2.rel.t0.combined = its2.rel.t0.combined %>%
  mutate_at(c(2:11), as.numeric)

# convert factors
its2.rel.t0.combined$treat = as.factor(its2.rel.t0.combined$treat)
its2.rel.t0.combined$treat = factor(its2.rel.t0.combined$treat, levels = c("Control", "Low Var","Mod Var","High Var"))
its2.rel.t0.combined$sitename = as.factor(its2.rel.t0.combined$sitename)
its2.rel.t0.combined$lineage = as.factor(its2.rel.t0.combined$lineage)
its2.rel.t0.combined$gen_site = as.factor(its2.rel.t0.combined$gen_site)
its2.rel.t0.combined$reef = as.factor(its2.rel.t0.combined$reef)

# add in dominant and minor distinctions to use in other plots
its2.rel.t0.combined.2 = its2.rel.t0.combined %>%
  mutate(dominant_type = case_when(A4z_sum >= 0.5 ~ "A4z",
                                   A4_sum >= 0.5 ~ "A4",
                              B19_sum >= 0.5 ~ "B19",
                              B5_sum >= 0.5 ~ "B5",
                              B5a_sum >= 0.5 ~ "B5a",
                              C1_sum >= 0.5 ~ "C1",
                              C3_sum >= 0.5 ~ "C3",
                              C3af_sum >= 0.5 ~ "C3af",
                              C15_sum >= 0.5 ~ "C15",
                              D1_sum >= 0.5 ~ "D1")) %>%
  mutate(minor_type = case_when(B19_sum < 0.5 & B19_sum > 0.0 ~ "B19",
                                B5_sum < 0.5 & B5_sum > 0.0 ~ "B5",
                                C1_sum < 0.5 & C1_sum > 0.0 ~ "C1",
                                C3_sum < 0.5 & C3_sum > 0.0 ~ "C3",
                                C3af_sum < 0.5 & C3af_sum > 0.0 ~ "C3af",
                                D1_sum < 0.5 & D1_sum > 0.0 ~ "D1"))

write.csv(its2.rel.t0.combined.2, file = "~/Dropbox/BU/TVE/TVE_Github/DielTempVariability/Physiology_Data/data_files/ITS2.dominanttype.T0.csv", row.names = FALSE)

## Now Prestress
# change column names to majority its2 sequence
taxa.ps = read.csv(file = "data_files/symportal_taxa.prestress.csv", header = TRUE) %>%
  select(-X)

taxa.ps$majority_its2
colnames(seqtab.rel.ps) = c("A4z","A4","B19","B5","B5a","C1.1","C1.2","C1.3","C1.4","C1.5","C1.6","C1.7","C3.1","C1.8","C3.2","C3af",
                            "C1.9","C1.10","C1.11","C1.12","C1.13","C1.14","C1.15","C1.16","C15","C1.17","C1.18","C1.19",
                            "D1.1","D1.2","D1.3","D1.4","D1.5","D1.6","D1.7","D1.8","D1.9" )

# make new data frame and sum columns with the same majority its2 sequence
seqtab.rel.ps.new = seqtab.rel.ps %>%
  mutate(A4z_sum = A4z) %>%
  mutate(A4_sum = A4) %>%
  mutate(B19_sum = B19) %>%
  mutate(B5_sum = B5) %>%
  mutate(B5a_sum = B5a) %>%
  mutate(C1_sum = rowSums(select(., starts_with("C1.")))) %>%
  mutate(C3_sum = rowSums(select(., starts_with("C3.")))) %>%
  mutate(C3af_sum = C3af) %>%
  mutate(C15_sum = C15) %>%
  mutate(D1_sum = rowSums(select(., starts_with("D1.")))) %>%
  rownames_to_column(var = "frag") %>%
  select(frag, contains("_sum"))

its2.rel.ps.combined = left_join(seqtab.rel.ps.new, samdf.rel.ps, by = "frag")

its2.rel.ps.combined = its2.rel.ps.combined %>%
  mutate_at(c(2:11), as.numeric)

# convert factors
its2.rel.ps.combined$treat = as.factor(its2.rel.ps.combined$treat)
its2.rel.ps.combined$treat = factor(its2.rel.ps.combined$treat, levels = c("Control", "Low Var","Mod Var","High Var"))
its2.rel.ps.combined$sitename = as.factor(its2.rel.ps.combined$sitename)
its2.rel.ps.combined$lineage = as.factor(its2.rel.ps.combined$lineage)
its2.rel.ps.combined$gen_site = as.factor(its2.rel.ps.combined$gen_site)
its2.rel.ps.combined$reef = as.factor(its2.rel.ps.combined$reef)

# add in dominant and minor distinctions to use in other plots
its2.rel.ps.combined.2 = its2.rel.ps.combined %>%
  mutate(dominant_type = case_when(A4z_sum >= 0.5 ~ "A4z",
                                   A4_sum >= 0.5 ~ "A4",
                                   B19_sum >= 0.5 ~ "B19",
                                   B5_sum >= 0.5 ~ "B5",
                                   B5a_sum >= 0.5 ~ "B5a",
                                   C1_sum >= 0.5 ~ "C1",
                                   C3_sum >= 0.5 ~ "C3",
                                   C3af_sum >= 0.5 ~ "C3af",
                                   C15_sum >= 0.5 ~ "C15",
                                   D1_sum >= 0.5 ~ "D1")) %>%
  mutate(minor_type = case_when(B19_sum < 0.5 & B19_sum > 0.0 ~ "B19",
                                B5_sum < 0.5 & B5_sum > 0.0 ~ "B5",
                                C1_sum < 0.5 & C1_sum > 0.0 ~ "C1",
                                C3_sum < 0.5 & C3_sum > 0.0 ~ "C3",
                                C3af_sum < 0.5 & C3af_sum > 0.0 ~ "C3af",
                                D1_sum < 0.5 & D1_sum > 0.0 ~ "D1"))

write.csv(its2.rel.ps.combined.2, file = "~/Dropbox/BU/TVE/TVE_Github/DielTempVariability/Physiology_Data/data_files/ITS2.dominanttype.prestress.csv", row.names = FALSE)

head(its2.rel.t0.combined.2)
head(its2.rel.ps.combined.2)

its2.rel.combined = rbind(its2.rel.t0.combined.2, its2.rel.ps.combined.2)
write.csv(its2.rel.combined, file = "~/Dropbox/BU/TVE/TVE_Github/DielTempVariability/Physiology_Data/data_files/ITS2.dominanttype.alltimepoints.csv", row.names = FALSE)

# summarize proportion of individuals with more than 50% C or D
its2.rel.combined.2 %>%
  #filter(sitename=="CI") %>%
  group_by(lineage) %>%
  dplyr::summarize(n = n(), # total sample size
            n_gt50D = sum(D1_sum>=0.5), # number of individuals with proportion of D1 > 0.5
            p_gt50D = n_gt50D/n) # proportion of individuals with D1 > 0.5
# All sites:
# lineage     n n_gt50D p_gt50D
# 1 L1         95      50   0.526
# 2 L2         61       9   0.148

# CI Only:
# lineage     n n_gt50D p_gt50D
# 1 L1         13       1  0.0769
# 2 L2         14       6  0.429

# Kruskal-Wallis (non-parametric alternative to one-way anova) to see if there is a difference
# in mean proportion D based on lineage
its2.rel.combined.CI = its2.rel.combined.2 %>%
  filter(sitename=="CI")

shapiro.test(its2.rel.combined.2$D1_sum) # not normal
kruskal.test(D1_sum ~ lineage, data = its2.rel.combined.2)
#data:  D1_sum by lineage
#Kruskal-Wallis chi-squared = 29.87, df = 1, p-value = 4.619e-08

kruskal.test(D1_sum ~ lineage, data = its2.rel.combined.CI)
#data:  propD by lineage
#Kruskal-Wallis chi-squared = 2.2793, df = 1, p-value = 0.1311

## Stopped here for the revision ##

# Make a plot of the proportion of durusdinium by lineage
library(Rmisc)

cols_lineage <- c("L1" = "#3f007d", "L2" = "#807dba", "L3" = "#bcbddc")

its2.rel.combined_plot = its2.rel.combined.2 %>%
  select(frag, treat, sitename, lineage, D1_sum)

d1_means = summarySE(its2.rel.combined_plot, measurevar = "D1_sum", groupvars = c("lineage"))

its2_dom_plot <- ggplot(its2.rel.combined_plot, aes(x = lineage, y = D1_sum))+
  theme_bw()+
  geom_jitter(aes(color = lineage, fill = lineage),
              #width = 0.25,
              alpha=0.2, pch = 21,
              color = "black") +
  geom_errorbar(data = d1_means, aes(x = lineage, ymax = D1_sum+se, ymin = D1_sum-se, color = lineage), width = .2, position = position_dodge(width=0.4)) +
  geom_point(data = d1_means, mapping = aes(x=lineage, y=D1_sum, color = lineage, fill = lineage), size = 3.5, pch = 21, color = "black", position = position_dodge(width=0.4))+
  scale_fill_manual(name = "Lineage",
                    breaks = c("L1","L2"),
                    values = cols_lineage)+
  scale_color_manual(name = "Lineage",
                     breaks = c("L1","L2"),
                     values = cols_lineage)+
  xlab("Lineage")+
  ylab("Proportion Durusdinium (+/- SE)")+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))
its2_dom_plot

ggsave(its2_dom_plot, filename = "/Users/hannahaichelman/Documents/BU/TVE/16S_ITS2/ITS_PreStress_Timepoint/ProportionD1_majorityITS2.pdf", width=4, height=4, units=c("in"), useDingbats=FALSE)


## Now make a new data frame of dominant DIVs not 'majority its2 sequences'
ps.cleanest.rel = readRDS("/Users/hannahaichelman/Documents/BU/TVE/16S_ITS2/ITS_PreStress_Timepoint/ps.its2.rel.RDS")
seqtab.rel <- data.frame(ps.cleanest.rel@otu_table)
samdf.rel <- data.frame(ps.cleanest.rel@sam_data)
tax.rel <- data.frame(ps.cleanest.rel@tax_table)

# make new data frame and sum columns with the same majority its2 sequence
seqtab.rel.new = seqtab.rel %>%
  rownames_to_column(var = "frag")

div.rel.combined = left_join(seqtab.rel.new, samdf.rel, by = "frag")

div.rel.combined = div.rel.combined %>%
  mutate_at(c(2:21), as.numeric)

# convert factors
div.rel.combined$treat = as.factor(div.rel.combined$treat)
div.rel.combined$treat = factor(div.rel.combined$treat, levels = c("Control", "Low Var","Mod Var","High Var"))
div.rel.combined$sitename = as.factor(div.rel.combined$sitename)
div.rel.combined$lineage = as.factor(div.rel.combined$lineage)
div.rel.combined$gen_site = as.factor(div.rel.combined$gen_site)
div.rel.combined$reef = as.factor(div.rel.combined$reef)

# add in dominant and minor distinctions to use in other plots
div.rel.combined.2 = div.rel.combined %>%
  mutate(dominant_div = case_when(B19.B19j >= 0.5 ~ "B19.B19j",
                                   B5 >= 0.5 ~ "B5",
                                   C1.C3.C1c.C1b.C72k.C1w >= 0.5 ~ "C1.C3.C1c.C1b.C72k.C1w",
                                   C1.C1c.C3.C1al.C1b >= 0.5 ~ "C1.C1c.C3.C1al.C1b",
                                   C1.C3.C1b.C1c.C42.2.C72k >= 0.5 ~ "C1.C3.C1b.C1c.C42.2.C72k",
                                   C3.C1.C3y.C6a >= 0.5 ~ "C3.C1.C3y.C6a",
                                   C3af.C3hq >= 0.5 ~ "C3af.C3hq",
                                   C3.C3y.C6a >= 0.5 ~ "C3.C3y.C6a",
                                   C1.C1c >= 0.5 ~ "C1.C1c",
                                   C1.C1c.C1al.C72k >= 0.5 ~ "C1.C1c.C1al.C72k",
                                   C1.C1c.C1b >= 0.5 ~ "C1.C1c.C1b",
                                   C1.C1c.C1am.C1ao.C1an.C3cm.C72k >= 0.5 ~ "C1.C1c.C1am.C1ao.C1an.C3cm.C72k",
                                   C1.C1c.C1al >= 0.5 ~ "C1.C1c.C1al",
                                   C3.C1 >= 0.5 ~ "C3.C1",
                                   C15 >= 0.5 ~ "C15",
                                   D1.D4.D1cl.D1ck.D4c.D1cm.D1c.D2 >= 0.5 ~ "D1.D4.D1cl.D1ck.D4c.D1cm.D1c.D2",
                                   D1.D4.D4c.D1cl.D1az.D1c.D2 >= 0.5 ~ "D1.D4.D4c.D1cl.D1az.D1c.D2",
                                   D1.D4.D4c.D1c >= 0.5 ~ "D1.D4.D4c.D1c",
                                   D1.D4.D4c >= 0.5 ~ "D1.D4.D4c",
                                   D1.D4.D4c.D2 >= 0.5 ~ "D1.D4.D4c.D2"))

write.csv(div.rel.combined.2, file = "/Users/hannahaichelman/Documents/BU/TVE/16S_ITS2/ITS_PreStress_Timepoint/ITS2.dominantDIVs.csv", row.names = FALSE)


#### Plasticity and Relative Change in Syms ####
## NOT USING SINCE ONLY LOOKING AT PRE-STRESS TIMEPOINT
library(tidyverse)
library(reshape2)
library(ggpubr)

# plot plasticity in proportion of Durusdinium symbionts relative to control, facet by lineage
cols_lineage <- c("L1" = "#3f007d", "L2" = "#807dba", "L3" = "#bcbddc")
cols_site <- c("CI" = "#543005", "PD"= "#bf812d",  "SP"= "#dfc27d",  "BN" = "#003c30", "BS"= "#35978f", "CA"= "#80cdc1")
cols_treat <- c("darkgrey", "#FF9966","#CC3300","#7f0000")

control_lowvar = its2.rel.combined %>%
  dplyr::filter(treat == "Control" | treat == "Low Var")

control_modvar = its2.rel.combined %>%
  dplyr::filter(treat == "Control" | treat == "Mod Var")

control_highvar = its2.rel.combined %>%
  dplyr::filter(treat == "Control" | treat == "High Var")

p.ds.lowvar = control_lowvar %>%
  drop_na(lineage) %>%
  ggplot(aes(x = treat, y = propD, color = gen_site, group = gen_site))+
  theme_bw()+
  geom_point(size = 4, pch = 1)+
  geom_line(aes(group = gen_site)) +
  xlab("Treatment")+
  ylab("Proportion Durusdinium")+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.title.x = element_blank(), legend.position = "none") +
  facet_wrap(~lineage)
p.ds.lowvar

p.ds.modvar = control_modvar %>%
  drop_na(lineage) %>%
  ggplot(aes(x = treat, y = propD, color = gen_site, group = gen_site))+
  theme_bw()+
  geom_point(size = 4, pch = 1)+
  geom_line(aes(group = gen_site)) +
  xlab("Treatment")+
  ylab("Proportion Durusdinium")+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.title.x = element_blank(), legend.position = "none") +
  facet_wrap(~lineage)
p.ds.modvar

p.ds.highvar = control_highvar %>%
  drop_na(lineage) %>%
  ggplot(aes(x = treat, y = propD, color = gen_site, group = gen_site))+
  theme_bw()+
  geom_point(size = 4, pch = 1)+
  geom_line(aes(group = gen_site)) +
  xlab("Treatment")+
  ylab("Proportion Durusdinium")+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.title.x = element_blank(), legend.position = "none") +
  facet_wrap(~lineage)
p.ds.highvar

all_ds = ggarrange(p.ds.lowvar, p.ds.modvar,p.ds.highvar,
                   labels = c("A", "B", "C"),
                   ncol = 3, nrow = 1)
ggsave(all_ds, filename = "/Users/hannahaichelman/Documents/BU/TVE/16S_ITS2/ITS_PreStress_Timepoint/PropD_alltreats.pdf", width=9, height=4, units=c("in"), useDingbats=FALSE)


# plot average proportion durusdinium by factors of interest
# violin plot - treatment
p.ds.lineage = its2.rel.combined %>%
  drop_na(lineage) %>%
  ggplot(aes(x = lineage, y = propD, fill = lineage))+
  geom_violin(trim = FALSE)+
  geom_point(size = 2, color = "black")+
  scale_fill_manual(values = cols_lineage)+
  theme_minimal()+
  xlab("Treatment")+
  ylab("Proportion Durusdinium")+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), legend.position = "none") +
  facet_wrap(~treat)
p.ds.lineage
ggsave(p.ds.lineage, filename = "/Users/hannahaichelman/Documents/BU/TVE/16S_ITS2/ITS_PreStress_Timepoint/PropD_violinplot_treat.pdf", width=5, height=5, units=c("in"), useDingbats=FALSE)

# violin plot - site
p.ds.site = its2.rel.combined %>%
  drop_na(lineage) %>%
  ggplot(aes(x = lineage, y = propD, fill = lineage))+
  geom_violin(trim = FALSE)+
  geom_point(size = 2, color = "black")+
  scale_fill_manual(values = cols_lineage)+
  theme_minimal()+
  xlab("Treatment")+
  ylab("Proportion Durusdinium")+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), legend.position = "none") +
  facet_wrap(~sitename)
p.ds.site
ggsave(p.ds.site, filename = "/Users/hannahaichelman/Documents/BU/TVE/16S_ITS2/ITS_PreStress_Timepoint/PropD_violinplot_site.pdf", width=6, height=5, units=c("in"), useDingbats=FALSE)

# read in data file and combine all post-med seqs into one proportion for genera A,B,C,D
data = read.delim(file = "/Users/hannahaichelman/Documents/BU/TVE/16S_ITS2/ITS_PreStress_Timepoint/ITS2.seqtab.rel.csv", sep = ",")%>%
  mutate(propA = rowSums(select(., contains("_A")))) %>%
  mutate(propB = rowSums(select(., contains("_B")))) %>%
  mutate(propC = rowSums(select(., contains("_C")))) %>%
  mutate(propD = rowSums(select(., contains("_D")))) %>%
  select(frag, propA, propB, propC, propD)

# read in sample info
samples = read.delim(file = "/Users/hannahaichelman/Documents/BU/TVE/16S_ITS2/ITS_All_Timepoints/ITS2.sample.data.rel.csv", sep = ",") %>%
  select(-X)

cols = c("frag", "timepoint", "Sample_or_Control", "treat", "gen_site", "sitename", "reef", "lineage")

# join together sample info and data
data_joined = left_join(samples,
                        data,
                        by = "frag") %>%
  mutate_at(cols, factor)

# add a column that indicates what the sample is dominated by, so we can make more categorical plots

data_joined = data_joined %>%
  mutate(dominant = case_when(propA >= 0.5 ~ "A",
                              propB >= 0.5 ~ "B",
                              propC >= 0.5 ~ "C",
                              propD >= 0.5 ~ "D")) %>%
  mutate(minor = case_when(propA < 0.5 & propA > 0.0 ~ "A",
                           propB < 0.5 & propB > 0.0 ~ "B",
                           propC < 0.5 & propC > 0.0 ~ "C",
                           propD < 0.5 & propD > 0.0 ~ "D"))
str(data_joined)
data_joined$dominant = as.factor(data_joined$dominant)
data_joined$minor = as.factor(data_joined$minor)
data_joined$treat = factor(data_joined$treat,
                               levels = c("Initial", "Control 1", "Low Var","Mod Var","High Var"))

data_joined_melt = data_joined %>%
  select(treat, timepoint, gen_site, sitename, lineage, dominant)
melt = melt(data_joined_melt, id.vars=c("treat","timepoint","gen_site","sitename","lineage"))

cols = c("#D53E4F", "#F46D43", "#FEE08B", "#66C2A5")

melt %>%
  drop_na(treat) %>%
  drop_na(lineage) %>%
  ggplot(aes(x = treat, y = value, fill = value)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cols,
                    breaks = c("A","B","C","D"),
                    name = "Dominant ITS2 Type") +
  xlab("Treatment") +
  ylab("Relative Abundance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  facet_wrap(~lineage, scale = "free")

ggsave(its2_barplot_site, filename = "/Users/hannahaichelman/Documents/BU/TVE/16S_ITS2/T0_Timepoint/its2_t0_CA.pdf", width=6, height=4, units=c("in"), useDingbats=FALSE)

# plot barplot of dominant type
dominant.barplot = data_joined %>%
  filter(treat != "Control 2") %>%
  drop_na(lineage) %>%
  ggplot(aes(x=treat, fill = dominant)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = cols,
                    breaks = c("A","B","C","D"),
                    name = "Dominant ITS2 Type") +
  xlab("Treatment") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  facet_wrap(~lineage)
dominant.barplot
ggsave(dominant.barplot, filename = "/Users/hannahaichelman/Documents/BU/TVE/16S_ITS2/ITS_All_Timepoints/dominant_type_prop.pdf", width=6, height=4, units=c("in"), useDingbats=FALSE)

# plot barplot of dominant type
minor.barplot = data_joined %>%
  filter(treat != "Control 2") %>%
  drop_na(lineage) %>%
  drop_na(minor) %>%
  ggplot(aes(x=treat, fill = minor)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = cols,
                    breaks = c("A","B","C","D"),
                    name = "Minor ITS2 Type") +
  xlab("Treatment") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  facet_wrap(~lineage)
minor.barplot
ggsave(minor.barplot, filename = "/Users/hannahaichelman/Documents/BU/TVE/16S_ITS2/ITS_All_Timepoints/minor_type_prop.pdf", width=6, height=4, units=c("in"), useDingbats=FALSE)



# Now do a kinda roundabout way to be able to plot individual samples and how they change proportion through the experiment
Initial = data_joined %>%
  filter(treat == "Initial") %>%
  mutate(Initial_D = propD) %>%
  mutate(Initial_C = propC) %>%
  select(gen_site, frag, Initial_D, Initial_C)

Control = data_joined %>%
  filter(treat == "Control 1") %>%
  mutate(Control_D = propD) %>%
  mutate(Control_C = propC) %>%
  select(gen_site, frag, Control_D, Control_C)

High = data_joined %>%
  filter(treat == "High Var") %>%
  mutate(High_D = propD) %>%
  mutate(High_C = propC) %>%
  select(gen_site, frag, High_D, High_C)

Low = data_joined %>%
  filter(treat == "Low Var") %>%
  mutate(Low_D = propD) %>%
  mutate(Low_C = propC) %>%
  select(gen_site, frag, Low_D, Low_C)

Mod = data_joined %>%
  filter(treat == "Mod Var") %>%
  mutate(Mod_D = propD) %>%
  mutate(Mod_C = propC) %>%
  select(gen_site, frag, Mod_D, Mod_C)

Combined = Initial %>%
  full_join(Control, by = "gen_site") %>%
  full_join(Low, by = "gen_site") %>%
  full_join(High, by = "gen_site") %>%
  full_join(Mod, by = "gen_site") %>%
#  na.omit() %>%
  mutate(Control_D_delta = Control_D - Initial_D,
         Low_D_delta = Low_D - Initial_D,
         High_D_delta = High_D - Initial_D,
         Mod_D_delta = Mod_D - Initial_D,
         Control_C_delta = Control_C - Initial_C,
         Low_C_delta = Low_C - Initial_C,
         High_C_delta = High_C - Initial_C,
         Mod_C_delta = Mod_C - Initial_C) %>%
  select(gen_site, Control_D_delta, Low_D_delta, High_D_delta, Mod_D_delta,
         Control_C_delta, Low_C_delta, High_C_delta, Mod_C_delta) %>%
  melt(id.vars = c("gen_site"))

Combined2 = Combined %>%
  left_join(samdf_t0, by = "gen_site")

Combined2_Cs = Combined2 %>%
  filter(str_detect(variable, '_C_'))

Combined2_Ds = Combined2 %>%
  filter(str_detect(variable, '_D_'))

str(Combined2_Cs)
Combined2_Cs$variable=as.factor(Combined2_Cs$variable)
Combined2_Cs$variable = factor(Combined2_Cs$variable,
                                       levels = c("Control_C_delta", "Low_C_delta", "Mod_C_delta","High_C_delta"))

Combined2_Ds$variable=as.factor(Combined2_Ds$variable)
Combined2_Ds$variable = factor(Combined2_Ds$variable,
                                       levels = c("Control_D_delta", "Low_D_delta", "Mod_D_delta","High_D_delta"))

# make some plots
cols_lineage <- c("L1" = "#fc8d59", "L2" = "#99d594")
cols_treatment <- c("cornflowerblue","darkolivegreen3","darkgoldenrod1","coral2")
cols_sites <- c("I-Cristobal" = "red4", "I-Punta Donato"= "indianred3",  "I-STRI Point"= "mistyrose3",  "O-Bastimentos N" = "royalblue4", "O-Bastimentos S"= "cornflowerblue", "O-Cayo de Agua"= "lightblue")

delta.Cs = Combined2_Cs %>%
  drop_na(lineage) %>%
  ggplot(aes(x = variable, y = value, color = gen_site, group = gen_site)) +
  theme_bw() +
  geom_point(size = 4) +
  #geom_line(aes(group = gen_site)) +
  ylab("Change in Cladocopium Proportion") +
  geom_hline(aes(yintercept = 0), color = "grey", linetype = "dashed", lwd = 1.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.title.x = element_blank(), legend.position = "none") +
  facet_wrap(~lineage)
delta.Cs

delta.Cs.box = Combined2_Cs %>%
  drop_na(lineage) %>%
  ggplot(aes(x=variable, y=value)) +
  geom_jitter(shape=16,
              size = 3,
              position=position_jitter(0.2),
              alpha=0.8,
              aes(color = variable)) +
  scale_color_manual(values = cols_treatment) + # for jittered points
  geom_boxplot(outlier.shape = NA,
               alpha = 0.85,
               aes(fill = variable))+
  scale_fill_manual(values = cols_treatment) + # for boxplot
  ylab("Change in Cladocopium Proportion") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.title.x = element_blank(), legend.position = "none") +
  facet_wrap(~lineage)
delta.Cs.box

delta.Ds = Combined2_Ds %>%
  drop_na(lineage) %>%
  ggplot(aes(x = variable, y = value, color = gen_site, group = gen_site)) +
  theme_bw() +
  geom_point(size = 4) +
  #geom_line(aes(group = gen_site)) +
  ylab("Change in Durusdinium Proportion") +
  geom_hline(aes(yintercept = 0), color = "grey", linetype = "dashed", lwd = 1.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.title.x = element_blank(), legend.position = "none") +
  facet_wrap(~lineage)
delta.Ds

delta.Ds.box = Combined2_Ds %>%
  drop_na(lineage) %>%
  ggplot(aes(x=variable, y=value)) +
  geom_jitter(shape=16,
              size = 3,
              position=position_jitter(0.2),
              alpha=0.8,
              aes(color = variable)) +
  scale_color_manual(values = cols_treatment) + # for jittered points
  geom_boxplot(outlier.shape = NA,
               alpha = 0.85,
               aes(fill = variable))+
  scale_fill_manual(values = cols_treatment) + # for boxplot
  ylab("Change in Durusdinium Proportion") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.title.x = element_blank(), legend.position = "none") +
  facet_wrap(~lineage)
delta.Ds.box

all_delta = ggarrange(delta.Cs, delta.Ds,
                   labels = c("A", "B"),
                   ncol = 2, nrow = 1)
ggsave(all_delta, filename = "/Users/hannahaichelman/Documents/BU/TVE/16S_ITS2/ITS_All_Timepoints/Delta_CandD_alltreats.pdf", width=12, height=8, units=c("in"), useDingbats=FALSE)

all_delta_box = ggarrange(delta.Cs.box, delta.Ds.box,
                      labels = c("A", "B"),
                      ncol = 2, nrow = 1)
ggsave(all_delta_box, filename = "/Users/hannahaichelman/Documents/BU/TVE/16S_ITS2/ITS_All_Timepoints/Delta_CandD_alltreats_box.pdf", width=12, height=8, units=c("in"), useDingbats=FALSE)
