# This script analyzes the ITS2 data from the Panama Daily Thermal Variability script.
# Author: Hannah E Aichelman
# Last Updated: July 1, 2022


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
# /Users/hannahaichelman/Documents/BU/TVE/16S_ITS2/PreStress_Timepoint/20210421_aichelman_PreStress/its2_type_profiles/154_20210426_DBV_20210427T024417.profiles.absolute.abund_and_meta.txt
setwd("/Users/hannahaichelman/Documents/BU/TVE/16S_ITS2/ITS_PreStress_Timepoint/20210421_aichelman_PreStress/its2_type_profiles")
its2_prestress = read.csv("SymPortal_PreStress_RawDIVs.csv")
head(its2_prestress)

# SymPortal ITS2 Type Analysis -
# consulted Ben Hume about aggregating data across the 'majority its2 type' - he agrees it is appropriate. also need to report how many DIVs are contained within each majority ITS2 sequence
# this is just the number of columns that have the same 'type' in the header
# i convert the div's to majority its2 data in the script below

# SymPortal Post Med Sequence Analysis - not using this one now
#setwd("/Users/hannahaichelman/Documents/BU/TVE/16S_ITS2/ITS_PreStress_Timepoint/20210421_aichelman_PreStress/post_med_seqs")
#its2_prestress_postmed = read.csv("/Users/hannahaichelman/Documents/BU/TVE/16S_ITS2/ITS_PreStress_Timepoint/20210421_aichelman_PreStress/post_med_seqs/SymPortal_PreStress_PostMed.csv")

#### Clean up input data ####

# Remove samples not included in this dataset
its2_prestress = its2_prestress %>%
  filter(frag != "Alexa2") %>%
  filter(frag != "Alexa1") %>%
  filter(frag != "MA2") %>%
  filter(frag != "ME1") %>%
  filter(frag != "negcontrol1") %>%
  filter(frag != "ITS-neg-old") %>%
  filter(frag != "negcontrol2") %>%
  filter(frag != "ITS-neg-new")

# Remove clades that are not symbiotic
its2_prestress2 = its2_prestress %>%
  select(-contains("E1c"))

# Look for and remove DIVs with 0 reads after removing samples
colSums(its2_prestress2[,-1])

# Remove DIVs that have 0 reads
its2_prestress3 = its2_prestress2 %>%
  select(-B2) %>%
  select(-B1)

str(its2_prestress3)

its2_prestress3 = its2_prestress3 %>%
  mutate_if(is.integer,as.numeric)

# Look for individual samples with 0 reads:
rowSums(its2_prestress3[, 2:21]) # I3H2 has 0 reads

# remove individuals with 0 reads and those in under-represented lineage 3
its2_prestress4 <- its2_prestress3 %>%
  filter(frag!= "I3H2") %>% #have to remove sample I3H2 because 0 reads
  filter(frag!="I2D4") %>% #lineage 3
  filter(frag!="I2I7") %>% #lineage 3
  filter(frag!="I2I4") %>% #lineage 3
  filter(frag!="I2I3") %>% #lineage 3
  filter(frag!="I2D1") %>% #lineage 3
  filter(frag!="I2H13") %>% #lineage 3
  filter(frag!="I2H2") %>% #lineage 3
  filter(frag!="I2H12") %>% #lineage 3
  filter(frag!="I2D3") %>% #lineage 3
  filter(frag!="I2I10") %>% #lineage 3
  filter(frag!="I2H4") %>% #lineage 3
  column_to_rownames("frag")

head(its2_prestress4)
str(its2_prestress4)

# keep this df for div analysis
its2_divs = its2_prestress4

# read in sample data
# DIVs:
samdf = read.csv("SampleInfo.csv")

# add identifying data
phys_metadata = read.csv("/Users/hannahaichelman/Documents/BU/TVE/phys_metadata.csv") %>%
  select(-dominant_type)

head(phys_metadata)

# combine with samdf
samdf = left_join(samdf, phys_metadata, by = "frag")

# remove the rows for samples that have been removed in data
samdf = samdf %>%
  filter(frag != "negcontrol1") %>%
  filter(frag != "ITS-neg-old") %>%
  filter(frag != "negcontrol2") %>%
  filter(frag != "ITS-neg-new") %>%
  filter(frag != "I3H2") %>% #have to remove sample I3H2 'cause 0 reads
  filter(is.na(lineage) | lineage!="L3") # remove lineage 3 from dataset but don't remove other NAs for lineage

str(samdf)
samdf$treat = as.factor(samdf$treat)
samdf$sitename = as.factor(samdf$sitename)
samdf$lineage = as.factor(samdf$lineage)

# change levels of treatment factor
samdf$treat = factor(samdf$treat,levels = c("Control","Control 2","Low Var","Mod Var","High Var"))

#rownames have to match between counts table & sample data table or else phyloseq will throw a fit
rownames(samdf) <- samdf$frag

#making a taxa table for phyloseq
taxa <- data.frame(colnames(its2_divs)) #extract sym data

colnames(taxa) <- c("DIV") #changing the column name to be more user-friendly

taxa$majority_its2 = c("B19",	"B5",	"C1",	"C1",	"C1",	"C3",	"C3af",	"C3",	"C1",	"C1",
                       "C1",	"C1",	"C1",	"C1",	"C15",	"D1",	"D1",	"D1",	"D1",	"D1")

taxa$genus = str_sub(taxa$DIV, 1, 1)

str(taxa)
taxa$DIV = as.factor(taxa$DIV)
taxa$majority_its2 = as.factor(taxa$majority_its2)
taxa$genus = as.factor(taxa$genus)

#rownmaes also have to match between the columns of the counts table & the taxa table or ps freaks out
rownames(taxa) <- taxa$DIV

taxa.m <- as.matrix(taxa) #also has to be a matrix

ps.its2 <- phyloseq(sample_data(samdf),
                    otu_table(its2_divs,taxa_are_rows=FALSE),
                    tax_table(taxa.m))
ps.its2
# DIVS:
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 20 taxa and 172 samples ]
# sample_data() Sample Data:       [ 172 samples by 7 sample variables ]
# tax_table()   Taxonomy Table:    [ 20 taxa by 3 taxonomic ranks ]


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
ps.rel <- transform_sample_counts(ps.its2, function(OTU) OTU/sum(OTU))
plot_bar(ps.rel, x="frag",fill="majority_its2")+
  theme_classic()

#keeps samples with summed counts greater than 0
ps.its2.no0 <- prune_samples(sample_sums(ps.its2)!=0, ps.its2)
ps.its2.no0 # don't lose any

# Remove NA's to create ps object we will plot
ps.cleaner <- subset_samples(ps.its2.no0,(!is.na(lineage)))
ps.cleaner2 <- subset_samples(ps.cleaner,(!is.na(treat)))
ps.cleanest <- subset_samples(ps.cleaner2,(treat!="Control 2")) #156 samples remain
ps.cleanest.rel <- transform_sample_counts(ps.cleanest, function(OTU) OTU/sum(OTU))

# save phyloseq object
saveRDS(ps.cleanest, "/Users/hannahaichelman/Documents/BU/TVE/16S_ITS2/ITS_PreStress_Timepoint/ps.its2.RDS")
saveRDS(ps.cleanest.rel, "/Users/hannahaichelman/Documents/BU/TVE/16S_ITS2/ITS_PreStress_Timepoint/ps.its2.rel.RDS")

# write this out as dataframe
seqtab.rel <- data.frame(ps.cleanest.rel@otu_table)
#write.csv(seqtab.rel, file="/Users/hannahaichelman/Documents/BU/TVE/16S_ITS2/ITS_PreStress_Timepoint/ITS2.seqtab.rel.csv")
sample.data = data.frame(sample_data(ps.cleanest.rel))
#write.csv(sample.data, file="/Users/hannahaichelman/Documents/BU/TVE/16S_ITS2/ITS_PreStress_Timepoint/ITS2.sample.data.rel.csv")

taxa <- data.frame(ps.cleanest.rel@tax_table)
mtaxa <- as.matrix(taxa)
#write.csv(taxa, file = "/Users/hannahaichelman/Documents/BU/TVE/16S_ITS2/ITS_PreStress_Timepoint/symportal_taxa.csv")

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
saveRDS(ps.norm, "/Users/hannahaichelman/Documents/BU/TVE/16S_ITS2/ITS_PreStress_Timepoint/ps.its2.norm.RDS")


#### PS Object Versions ####
ps.cleanest = readRDS("/Users/hannahaichelman/Documents/BU/TVE/16S_ITS2/ITS_PreStress_Timepoint/ps.its2.RDS")
seqtab <- data.frame(ps.cleanest@otu_table)
samdf <- data.frame(ps.cleanest@sam_data)

ps.cleanest.rel = readRDS("/Users/hannahaichelman/Documents/BU/TVE/16S_ITS2/ITS_PreStress_Timepoint/ps.its2.rel.RDS")
seqtab.rel <- data.frame(ps.cleanest.rel@otu_table)
samdf.rel <- data.frame(ps.cleanest.rel@sam_data)

taxa = read.csv(file = "/Users/hannahaichelman/Documents/BU/TVE/16S_ITS2/ITS_PreStress_Timepoint/symportal_taxa.csv", header = TRUE) %>%
  select(-X)
rownames(taxa) <- as.factor(taxa$DIV)
sum(taxa$genus == "B") # 2 div
sum(taxa$genus == "C") # 13 div
sum(taxa$genus == "D") # 5 div

mtaxa <- as.matrix(taxa)

#### Bar plot - post-processing ####
#its2_cols_end = c("#D53E4F", "#F46D43", "#3288BD", "#FEE08B", "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD", "#FDAE61")
#its2_cols_purple = c("white", "#efedf5","#bcbddc", "#756bb1")
its2_cols_greens = c("B19" = "#ffeda0", "B5" = "#fd8d3c",
                     "C1" = "#edf8e9", "C15" = "#feb24c",  "C3" = "#a1d99b",
                     "C3af" = "#238b45", "D1" = "#00441b")

its2_cols_blues = c("#f6eff7", "#d0d1e6", "#a6bddb", "#67a9cf", "#3690c0", "#02818a", "#016450")

plot_bar(ps.cleanest.rel, x="frag", fill="majority_its2") +
  theme_bw()

# subset by factors we need:
ps.cleanest.rel.control = subset_samples(ps.cleanest.rel, treat=="Control")
ps.cleanest.rel.lowvar = subset_samples(ps.cleanest.rel, treat=="Low Var")
ps.cleanest.rel.modvar = subset_samples(ps.cleanest.rel, treat=="Mod Var")
ps.cleanest.rel.highvar = subset_samples(ps.cleanest.rel, treat=="High Var")

ps.cleanest.rel.CI = subset_samples(ps.cleanest.rel, sitename=="CI")
ps.cleanest.rel.CI.control = subset_samples(ps.cleanest.rel.CI, treat=="Control")
ps.cleanest.rel.CI.lowvar = subset_samples(ps.cleanest.rel.CI, treat=="Low Var")
ps.cleanest.rel.CI.modvar = subset_samples(ps.cleanest.rel.CI, treat=="Mod Var")
ps.cleanest.rel.CI.highvar = subset_samples(ps.cleanest.rel.CI, treat=="High Var")


# Now facet wrap/aggregate by interesting factors

# Lineage:
ps.lin <- merge_samples(ps.cleanest, "lineage")
ps.rel.lin <- transform_sample_counts(ps.lin, function(x) x / sum(x))

p.lineage = plot_bar(ps.rel.lin, fill="majority_its2") +
  geom_bar(stat = "identity") +
  scale_fill_manual(name = "Majority ITS2", values = its2_cols_greens) +
  theme_bw()
p.lineage
ggsave(p.lineage, filename = "/Users/hannahaichelman/Documents/BU/TVE/16S_ITS2/ITS_PreStress_Timepoint/syms_lineage.pdf", width=4, height=3.5, units=c("in"), useDingbats=FALSE)

# Treatment:
p.treat = plot_bar(ps.cleanest.rel, x="frag", fill="majority_its2") +
  geom_bar(stat = "identity") +
  scale_fill_manual(name = "Majority ITS2", values = its2_cols_greens) +
  theme_bw() +
  facet_wrap(~treat, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
p.treat
ggsave(p.treat, filename = "/Users/hannahaichelman/Documents/BU/TVE/16S_ITS2/ITS_PreStress_Timepoint/syms_treatment.pdf", width=10, height=8, units=c("in"), useDingbats=FALSE)

# Treatment faceted by lineage
p.treat.control = plot_bar(ps.cleanest.rel.control, x="frag", fill="majority_its2") +
  geom_bar(stat = "identity") +
  scale_fill_manual(name = "Majority ITS2", values = its2_cols_greens) +
  theme_bw() +
  facet_wrap(~lineage, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.title.x = element_blank())+
  coord_cartesian(expand=FALSE)
p.treat.control

p.treat.low = plot_bar(ps.cleanest.rel.lowvar, x="frag", fill="majority_its2") +
  geom_bar(stat = "identity") +
  scale_fill_manual(name = "Majority ITS2", values = its2_cols_greens) +
  theme_bw() +
  facet_wrap(~lineage, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.title.x = element_blank(), axis.title.y = element_blank())+
  coord_cartesian(expand=FALSE)
p.treat.low

p.treat.mod = plot_bar(ps.cleanest.rel.modvar, x="frag", fill="majority_its2") +
  geom_bar(stat = "identity") +
  scale_fill_manual(name = "Majority ITS2", values = its2_cols_greens) +
  theme_bw() +
  facet_wrap(~lineage, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.title.x = element_blank())+
  coord_cartesian(expand=FALSE)
p.treat.mod

p.treat.high = plot_bar(ps.cleanest.rel.highvar, x="frag", fill="majority_its2") +
  geom_bar(stat = "identity") +
  scale_fill_manual(name = "Majority ITS2", values = its2_cols_greens) +
  theme_bw() +
  facet_wrap(~lineage, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.title.x = element_blank(), axis.title.y = element_blank())+
  coord_cartesian(expand=FALSE)
p.treat.high

p.alltreat = ggarrange(p.treat.control, p.treat.low, p.treat.mod, p.treat.high,
                       labels = "(a)", legend = "right", ncol = 2, nrow = 2,
                       label.x = FALSE, common.legend = TRUE)
ggsave(p.alltreat, filename = "/Users/hannahaichelman/Documents/BU/TVE/16S_ITS2/ITS_PreStress_Timepoint/syms_treatment_lineagefacet_majorityITS2.pdf", width=8, height=4, units=c("in"), useDingbats=FALSE)


# Make the same figure but for CI only to compare trends
p.treat.control.CI = plot_bar(ps.cleanest.rel.CI.control, x="frag", fill="majority_its2") +
  geom_bar(stat = "identity") +
  scale_fill_manual(name = "Genus", values = its2_cols_greens) +
  theme_bw() +
  facet_wrap(~lineage, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.title.x = element_blank())+
  coord_cartesian(expand=FALSE)
p.treat.control.CI

p.treat.low.CI = plot_bar(ps.cleanest.rel.CI.lowvar, x="frag", fill="majority_its2") +
  geom_bar(stat = "identity") +
  scale_fill_manual(name = "Genus", values = its2_cols_greens) +
  theme_bw() +
  facet_wrap(~lineage, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.title.x = element_blank(), axis.title.y = element_blank())+
  coord_cartesian(expand=FALSE)
p.treat.low.CI

p.treat.mod.CI = plot_bar(ps.cleanest.rel.CI.modvar, x="frag", fill="majority_its2") +
  geom_bar(stat = "identity") +
  scale_fill_manual(name = "Genus", values = its2_cols_greens) +
  theme_bw() +
  facet_wrap(~lineage, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.title.x = element_blank())+
  coord_cartesian(expand=FALSE)
p.treat.mod.CI

p.treat.high.CI = plot_bar(ps.cleanest.rel.CI.highvar, x="frag", fill="majority_its2") +
  geom_bar(stat = "identity") +
  scale_fill_manual(name = "Genus", values = its2_cols_greens) +
  theme_bw() +
  facet_wrap(~lineage, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.title.x = element_blank(), axis.title.y = element_blank())+
  coord_cartesian(expand=FALSE)
p.treat.high.CI

p.alltreat.CI = ggarrange(p.treat.control.CI, p.treat.low.CI, p.treat.mod.CI, p.treat.high.CI,
                          labels = "(a)", legend = "right", ncol = 2, nrow = 2,
                          label.x = FALSE, common.legend = TRUE)
ggsave(p.alltreat.CI, filename = "/Users/hannahaichelman/Documents/BU/TVE/16S_ITS2/ITS_PreStress_Timepoint/syms_treatment_lineagefacet_CIonly_majorityITS2.pdf", width=12, height=8, units=c("in"), useDingbats=FALSE)

# Lineage only:
p.lin = plot_bar(ps.cleanest.rel, x="frag", fill="majority_its2") +
  geom_bar(stat = "identity") +
  scale_fill_manual(name = "Genus", values = its2_cols_greens) +
  theme_bw() +
  facet_wrap(~lineage, scales = "free_x") +
  #theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.title.x = element_blank(), axis.title.y = element_blank())
p.lin
ggsave(p.lin, filename = "/Users/hannahaichelman/Documents/BU/TVE/16S_ITS2/ITS_PreStress_Timepoint/syms_lineage_majorityITS2.pdf", width=10, height=6, units=c("in"), useDingbats=FALSE)

# Site of Origin:
p.site = plot_bar(ps.cleanest.rel, x="frag", fill="majority_its2") +
  geom_bar(stat = "identity") +
  #scale_fill_manual(values = its2_cols_greens) +
  theme_bw() +
  facet_wrap(~sitename, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.title.x = element_blank(), axis.title.y = element_blank())
p.site
ggsave(p.site, filename = "/Users/hannahaichelman/Documents/BU/TVE/16S_ITS2/ITS_PreStress_Timepoint/syms_site_majorityITS2.pdf", width=10, height=6, units=c("in"), useDingbats=FALSE)

#### Bray Curtis Dissimilarity PCoAs ####
# color schemes
#cols_site <- c("I-Cristobal" = "red4", "I-Punta Donato"= "indianred3",  "I-STRI Point"= "mistyrose3",  "O-Bastimentos N" = "royalblue4", "O-Bastimentos S"= "cornflowerblue", "O-Cayo de Agua"= "lightblue")
cols_site_diverging <- c("CI" = "#543005", "PD"= "#bf812d",  "SP"= "#dfc27d",  "BN" = "#003c30", "BS"= "#35978f", "CA"= "#80cdc1")

#cols_treatment <- c("cornflowerblue","darkolivegreen3","darkgoldenrod1","coral2")
cols_treat_reds <- c("darkgrey", "#FF9966","#CC3300","#7f0000")

cols_lineage_purples <- c("L1" = "#3f007d", "L2" = "#807dba")

ps.ord <- ordinate(ps.cleanest.rel,"PCoA",distance="bray")

pcoa.site = plot_ordination(ps.cleanest.rel, ps.ord, color ="sitename", shape = "sitename")+
  geom_point(alpha=0.8)+
  scale_color_manual(name = "Site", breaks = c("BN","BS","CA","CI","PD","SP"), values=cols_site_diverging)+
  scale_shape_manual(name = "Site", breaks = c("BN","BS","CA","CI","PD","SP"), values=c(15,16,17,22,21,24))+
  stat_ellipse() +
  xlim(-1.3,1) +
  ylim(-1.2,1) +
  theme_bw()

pcoa.treat = plot_ordination(ps.cleanest.rel, ps.ord, color ="treat", shape = "treat")+
  geom_point(alpha=0.8)+
  scale_color_manual(name = "DTV Treatment", values=cols_treat_reds)+
  scale_shape_manual(name = "DTV Treatment", values=c(15,19,17,18))+
  stat_ellipse()+
  xlim(-1.3,1) +
  ylim(-1.2,1) +
  theme_bw()
#ggsave(pcoa.treat, filename = "/Users/hannahaichelman/Documents/BU/TVE/16S_ITS2/ITS_PreStress_Timepoint/bray_pcoa_treat.pdf", width=4.5, height=3, units=c("in"), useDingbats=FALSE)

pcoa.lineage = plot_ordination(ps.cleanest.rel, ps.ord, color ="lineage")+
  geom_point(alpha=0.8)+
  scale_color_manual(name = "Lineage", values=cols_lineage_purples)+
  stat_ellipse()+
  xlim(-1.3,1) +
  ylim(-1.2,1) +
  theme_bw()
#ggsave(pcoa.lineage, filename = "/Users/hannahaichelman/Documents/BU/TVE/16S_ITS2/ITS_PreStress_Timepoint/bray_pcoas_lineage.pdf", width=4, height=3, units=c("in"), useDingbats=FALSE)

gg.pcoa <- ggarrange(pcoa.site,pcoa.treat,pcoa.lineage,labels=c("A.","B.","C."),nrow=1,common.legend=F,legend="none")
ggsave(gg.pcoa, filename = "/Users/hannahaichelman/Documents/BU/TVE/16S_ITS2/ITS_PreStress_Timepoint/bray_pcoas_rel_divs.pdf", width=8, height=3, units=c("in"), useDingbats=FALSE)


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

#### PCoA Stats ####
library(vegan)
#remotes::install_github("Jtrachsel/funfuns")
library(funfuns)
library(dplyr)
library(edgeR)

## Raw (cleaned)
ps = readRDS("/Users/hannahaichelman/Documents/BU/TVE/16S_ITS2/ITS_PreStress_Timepoint/ps.its2.RDS")

seq.ps <- data.frame(ps@otu_table)
samdf.ps <- data.frame(ps@sam_data)
dist.ps <- vegdist(seq.ps)

# dispersion
# by treatment
bet.ps <- betadisper(dist.ps,samdf.ps$treat)
anova(bet.ps) #ns, p=0.697
# Response: Distances
# Df  Sum Sq   Mean Sq F value Pr(>F)
# Groups      3 0.00634 0.0021124  0.4791 0.6973
# Residuals 152 0.67011 0.0044086
permutest(bet.ps,pairwise=TRUE,permutations=999) # all pairwise comparisons ns
plot(bet.ps)

# by sitename
bet.ps <- betadisper(dist.ps,samdf.ps$sitename)
anova(bet.ps) #p=0.3244
# Response: Distances
#             Df Sum Sq  Mean Sq      F N.Perm Pr(>F)
# Groups      5 0.2839 0.056789 1.1745    999  0.313
# Residuals 150 7.2526 0.048351
permutest(bet.ps,pairwise=TRUE,permutations=999)
plot(bet.ps)

# by lineage
bet.ps <- betadisper(dist.ps,samdf.ps$lineage)
anova(bet.ps) #p=0.7464
# Response: Distances
#             Df  Sum Sq   Mean Sq F value Pr(>F)
# Groups      1 0.00102 0.0010190   0.105 0.7464
# Residuals 154 1.49501 0.0097078
permutest(bet.ps,pairwise=TRUE,permutations=999) # p=0.756
plot(bet.ps)

# adonis
adonis2(formula = seq.ps ~ treat + sitename + lineage, data = samdf.ps, permutations = 999)
#           Df SumOfSqs      R2      F Pr(>F)
# treat      3    1.496 0.02310 1.4312  0.073 .
# sitename   5   11.616 0.17938 6.6691  0.001 ***
# lineage    1    0.784 0.01211 2.2507  0.023 *
# Residual 146   50.859 0.78541
# Total    155   64.754 1.00000

## Relative abundance
# Report this since it is the data used for making the PCA's
ps.cleanest.rel = readRDS("/Users/hannahaichelman/Documents/BU/TVE/16S_ITS2/ITS_PreStress_Timepoint/ps.its2.rel.RDS")

# dispersion by treatment
seq.ps <- data.frame(ps.cleanest.rel@otu_table)
samdf.ps <- data.frame(ps.cleanest.rel@sam_data)
dist.ps <- vegdist(seq.ps)
bet.ps <- betadisper(dist.ps,samdf.ps$treat)
anova(bet.ps) #p=0.8568
#             Df  Sum Sq   Mean Sq F value Pr(>F)
# Groups      3 0.00413 0.0013761  0.2562 0.8568
# Residuals 152 0.81649 0.0053716
permutest(bet.ps,pairwise=TRUE,permutations=999) # all ns
plot(bet.ps)

# by sitename
bet.ps <- betadisper(dist.ps,samdf.ps$sitename)
anova(bet.ps) #p=0.1152
#             Df  Sum Sq Mean Sq F value Pr(>F)
# Groups      5  1.3691 0.27382  1.8048 0.1152
# Residuals 150 22.7571 0.15171
permutest(bet.ps,pairwise=TRUE,permutations=999) # CI-BS, CI-CA, PD-CA significant differences
plot(bet.ps)

# by lineage
bet.ps <- betadisper(dist.ps,samdf.ps$lineage)
anova(bet.ps) #p=0.3339
#             Df Sum Sq  Mean Sq F value Pr(>F)
# Groups      1 0.0356 0.035633  0.9397 0.3339
# Residuals 154 5.8397 0.037920
plot(bet.ps)

# adonis
head(adonis2(seq.ps ~ treat+sitename+lineage, data=samdf.ps, permutations=999))
adonis2(formula = seq.ps ~ treat + sitename + lineage, data = samdf.ps, permutations = 999)
#           Df SumOfSqs      R2      F Pr(>F)
# treat      3    0.973 0.01574 1.0206  0.402
# sitename   5   13.595 0.21989 8.5538  0.001 ***
# lineage    1    0.849 0.01373 2.6712  0.025 *
# Residual 146   46.411 0.75063
# Total    155   61.828 1.00000

pairwise.adonis(seq.ps, factors=samdf.ps$lineage, permutations=999) #no significant comparisons
pairwise.adonis(seq.ps, factors=samdf.ps$sitename, permutations=999) #all comparisons different except PD vs SP and BN vs CA


#### Summary of ITS2 DIVs ####
ps.cleanest.rel = readRDS("/Users/hannahaichelman/Documents/BU/TVE/16S_ITS2/ITS_PreStress_Timepoint/ps.its2.rel.RDS")
seqtab.rel <- data.frame(ps.cleanest.rel@otu_table)
samdf.rel <- data.frame(ps.cleanest.rel@sam_data)
tax.rel <- data.frame(ps.cleanest.rel@tax_table)

# change column names to majority its2 sequence
colnames(seqtab.rel) = c("B19","B5","C1.1","C1.2","C1.3","C3.1","C3af",
                         "C3.2","C1.4","C1.5","C1.6","C1.7","C1.8","C1.9","C15",
                         "D1.1","D1.2","D1.3","D1.4","D1.5")

# make new data frame and sum columns with the same majority its2 sequence
seqtab.rel.new = seqtab.rel %>%
  mutate(B19_sum = B19) %>%
  mutate(B5_sum = B5) %>%
  mutate(C1_sum = rowSums(select(., starts_with("C1.")))) %>%
  mutate(C3_sum = rowSums(select(., starts_with("C3.")))) %>%
  mutate(C3af_sum = C3af) %>%
  mutate(C15_sum = C15) %>%
  mutate(D1_sum = rowSums(select(., starts_with("D1.")))) %>%
  rownames_to_column(var = "frag") %>%
  select(frag, contains("_sum"))

its2.rel.combined = left_join(seqtab.rel.new, samdf.rel, by = "frag")

its2.rel.combined = its2.rel.combined %>%
  mutate_at(c(2:8), as.numeric)

# convert factors
its2.rel.combined$treat = as.factor(its2.rel.combined$treat)
its2.rel.combined$treat = factor(its2.rel.combined$treat, levels = c("Control", "Low Var","Mod Var","High Var"))
its2.rel.combined$sitename = as.factor(its2.rel.combined$sitename)
its2.rel.combined$lineage = as.factor(its2.rel.combined$lineage)
its2.rel.combined$gen_site = as.factor(its2.rel.combined$gen_site)
its2.rel.combined$reef = as.factor(its2.rel.combined$reef)

# add in dominant and minor distinctions to use in other plots
its2.rel.combined.2 = its2.rel.combined %>%
  mutate(dominant_type = case_when(B19_sum >= 0.5 ~ "B19",
                                   B5_sum >= 0.5 ~ "B5",
                                   C1_sum >= 0.5 ~ "C1",
                                   C15_sum >= 0.5 ~ "C15",
                                   C3_sum >= 0.5 ~ "C3",
                                   C3af_sum >= 0.5 ~ "C3af",
                                   D1_sum >= 0.5 ~ "D1")) %>%
  mutate(minor_type = case_when(B19_sum < 0.5 & B19_sum > 0.0 ~ "B19",
                                B5_sum < 0.5 & B5_sum > 0.0 ~ "B5",
                                C1_sum < 0.5 & C1_sum > 0.0 ~ "C1",
                                C15_sum < 0.5 & C15_sum > 0.0 ~ "C15",
                                C3_sum < 0.5 & C3_sum > 0.0 ~ "C3",
                                C3af_sum < 0.5 & C3af_sum > 0.0 ~ "C3af",
                                D1_sum < 0.5 & D1_sum > 0.0 ~ "D1"))

write.csv(its2.rel.combined.2, file = "/Users/hannahaichelman/Documents/BU/TVE/16S_ITS2/ITS_PreStress_Timepoint/ITS2.dominanttype.csv", row.names = FALSE)


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
