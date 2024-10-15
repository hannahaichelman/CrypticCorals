#### Set Up ####
library(plyr)
library(dplyr)
library(tidyverse)
library(vegan)

#setwd('~/Dropbox/BU/TVE/TVE_Github/DielTempVariability/2bRAD/')
bams=data.frame(read.table("2bRAD/data_files/bams.txt", header=FALSE)) # list of bam files
colnames(bams)<- "bam"
#goods=c(1:length(bams))

# reading table of pairs of replicates (tab-delimited) - skip if there are no clones
# clonepairs=read.table("clonepairs.tab",sep="\t")
# repsa= clonepairs[,1]
# repsb= clonepairs[,2]
# # removing "b" replicates
# goods=which(!(bams %in% repsb))

#### Metadata ####

# loading individual to population correspondences
i2p=read.table("2bRAD/data_files/bam_barcode_names_tuftscustom.csv",sep=",",header=TRUE) # 2-column tab-delimited table of individual assignments to populations; must be in the same order as samples in the bam list or vcf file.
row.names(i2p)=i2p[,1]
#i2p=i2p[goods,]
site=i2p[,2]

# add in site name
i2p$sitename <- ifelse(i2p$pop == 'I2', 'SP',
                       ifelse(i2p$pop == 'I3', 'CI',
                       ifelse(i2p$pop == 'I4', 'PD',
                       ifelse(i2p$pop == 'O2', 'BS',
                       ifelse(i2p$pop == 'O3', 'CA',
                       'BN')))))

# create new data frame using i2p without clones (duplicated preps), but still with two actual clones (I4G + I4F)
# this is removing files that have the lower coverage/fewer reads
i2p_noclones_allsamps = i2p %>%
  dplyr::filter(bam != "4-MullenDavies_S4_TCAC.nosymbio.fastq.bam") %>% # O4E
  dplyr::filter(bam != "1-MullenDavies_S1_GCTT.nosymbio.fastq.bam") %>% # O2A
  dplyr::filter(prepped_id != "I4A_CLONE") %>% # I4A
  dplyr::filter(bam != "2-MullenDavies_S2_AGTG.nosymbio.fastq.bam") %>% # O4A
  dplyr::filter(bam != "3-MullenDavies_S3_TGGT.nosymbio.fastq.bam") %>% # I3H
  dplyr::filter(prepped_id != "O3B_CLONE") %>% # O3B
  dplyr::filter(bam != "4-MullenDavies_S4_TCAG.nosymbio.fastq.bam") %>% # O2E
  dplyr::filter(bam != "1-MullenDavies_S1_TGGT.nosymbio.fastq.bam") %>% # I2E sample 1
  dplyr::filter(bam != "5-MullenDavies_S5_TGGT.nosymbio.fastq.bam") %>% # I2E sample 2
  dplyr::filter(bam != "1-MullenDavies_S1_ACCA.nosymbio.fastq.bam") # I3C

# make a new bams file with this filtered i2p file and write out a csv file
bams_noclones_allsamps = i2p_noclones_allsamps %>%
  select(bam)
#write.csv(bams_noclones_allsamps, file = "/Users/hannahaichelman/Documents/BU/TVE/2bRAD/Analysis/tuftscustompipeline_denovo_nosyms/bams_noclones_allsamps.txt", row.names = FALSE)

# Remove I4G so this is actually no clones here - I4G was a smaller bam file
i2p_noclones = i2p_noclones_allsamps %>%
  dplyr::filter(bam != "1-MullenDavies_S1_GTGA.nosymbio.fastq.bam") # I4G

# make new bams file for this filtered i2p file with no clones at all
bams_noclones = i2p_noclones %>%
  select(bam)
#write.csv(bams_noclones, file = "/Users/hannahaichelman/Documents/BU/TVE/2bRAD/Analysis/tuftscustompipeline_denovo_nosyms/bams_noclones.csv", row.names = FALSE)

site=i2p_noclones[,8]

# setting up colors for plotting
palette(rainbow(length(unique(site))))
colors=as.numeric(as.factor(site))
colpops=as.numeric(as.factor(sort(unique(site))))

cols_site <- c("CI" = "#543005", "PD"= "#bf812d",  "SP"= "#dfc27d",  "BN" = "#003c30", "BS"= "#35978f", "CA"= "#80cdc1")
cols_treat <- c("darkgrey", "#FF9966","#CC3300","#7f0000")
cols_lineage <- c("#bcbddc","#756bb1")

# load in symbiont type data to create combined name
#sym=read.table("/Users/hannahaichelman/Documents/BU/TVE/16S_ITS2/Laura_SymPortal_InitialData_DominantTypes.csv", sep=",", header=TRUE)

# merge the i2p and sym type dataframes together, keep the order of the i2p because this is important! matches with matrix rows and columns below
# i2p_sym = join(i2p, sym, by="sample_id")
# i2p_sym$id_symtype = paste(i2p_sym$sample_id, "_", i2p_sym$Dominant_Type)

#### PCoA based on IBS ####
# clustering / PCoA based on identity by state (IBS) based on single read resampling
# (for low and/or uneven coverage)
# all clones removed
ma = as.matrix(read.table("2bRAD/data_files/myresult2.noclone.ibsMat"))
colnames(ma)=i2p_noclones$sample_id
rownames(ma)=i2p_noclones$sample_id

# to make cluster dendrogram with all samples included:
ma = as.matrix(read.table("2bRAD/data_files/myresult2.ibsMat"))
colnames(ma)=i2p$sample_id
rownames(ma)=i2p$sample_id

hc=hclust(as.dist(ma),"ave")
plot(hc,cex=0.5)  # this shows how similar clones are
abline(h=0.265, lwd = 2, lty = 2, col = "grey")
# exported this dendrogram as a pdf

# performing PCoA and CAP
conds=data.frame(cbind(site))
pp0=capscale(ma~1)
pp=capscale(ma~site,conds)

# significance of by-site divergence
adonis2(ma~site,conds)
#         Df SumOfSqs      R2      F Pr(>F)
# site      5  0.40473 0.21067 2.3487  0.001 ***
# Residual 44  1.51640 0.78933
# Total    49  1.92113 1.00000

# eigenvectors
plot(pp0$CA$eig)

# find % variance explained - eigenvector for MDS1 / sum of remaining eigenvectors
eigs = as.data.frame(pp0$CA$eig)
eigs$MDS = rownames(eigs)
head(eigs)
# % variance explained by MDS1 = 26.7% variance
eigs[1,1]/sum(eigs[2:49, 1])

# % variance explained by MDS2 = 3.9% variance explained
eigs[2,1]/sum(eigs[3:49, 1])


axes2plot=c(1,2)
quartz()
library(adegenet) # for transp()
cmd=pp0
plot(cmd,choices=axes2plot,display="sites",type="n") # choices - axes to display
points(cmd,choices=axes2plot,pch=19,col=transp(colors,alpha=0.7))
#ordihull(cmd,choices= axes2plot,groups= conds$grp,draw="polygon",col=1+as.numeric(unique(as.factor(conds$grp))),label=T)
ordispider(cmd,choices= axes2plot,groups=conds$site,col="grey80")
ordiellipse(cmd,choices= axes2plot,groups= conds$site,draw="polygon",col=colpops,label=T)

# unscaled, to identify outliers
plot(cmd$CA$u[,axes2plot],pch=19,col=colors)
ordispider(cmd$CA$u[,axes2plot],groups=conds$site,col="grey80")
ordiellipse(cmd$CA$u[,axes2plot],groups= conds$site,draw="polygon",col=colpops,label=T)
identify(cmd$CA$u[,axes2plot],labels=colnames(ma),n=3,cex=0.7)
# I2I, I2D, I2H are the three outliers from I2 here - don't correspond to lowest or highest coverage in quality.txt file

# prettier ggplot option
axes2plot=c(1,2)
pca_s <- as.data.frame(cmd$CA$u[,axes2plot])
#colors=c('royalblue4','cornflowerblue','lightblue','red4','indianred3','mistyrose3')
cols_site <- c("CI" = "#543005", "PD"= "#bf812d",  "SP"= "#dfc27d",  "BN" = "#003c30", "BS"= "#35978f", "CA"= "#80cdc1")

MDS1_2 = ggplot(pca_s, aes(MDS1, MDS2)) +
  theme_bw() +
  geom_point(aes(colour=conds$site, shape =conds$site), size=2, stroke = 1)  +
  stat_ellipse(type = "t",
               aes(color = conds$site), show.legend = NA, lwd = 1) +
  scale_color_manual(values=cols_site,
                     breaks = c("BN", "BS", "CA", "CI", "PD", "SP"),
                     name = "Site") +
  scale_shape_manual(values = c(15,16,17,22,21,24),
                     breaks = c("BN", "BS", "CA", "CI", "PD", "SP"),
                     name = "Site") +
  xlab("MDS1 (26.7% explained variance)") +
  ylab("MDS2 (3.9% explained variance)")
MDS1_2

ggsave(MDS1_2, filename = "/Users/hannahaichelman/Documents/BU/TVE/2bRAD/Analysis/tuftscustompipeline_denovo_nosyms/MDS1_2_sitename.pdf", width=6, height=5, units=c("in"), useDingbats=FALSE)


#plotting 2 and 3
axes2plot=c(3,4)
pca_s <- as.data.frame(cmd$CA$u[,axes2plot])

MDS3_4 = ggplot(pca_s, aes(MDS3, MDS4)) +
  theme_bw() +
  geom_point(aes(colour=conds$site), size=2, stroke = 1)  +
  scale_color_manual(values = cols_site,
                     breaks = c("BN", "BS", "CA", "CI", "PD", "SP"),
                     name = "Site") +
  stat_ellipse(geom="polygon", alpha = .2,
               aes(fill = conds$site), show.legend = NA) + scale_fill_manual(values=colors, name = "Site") +
  xlab(paste0("MDS3")) +
  ylab(paste0("MDS4"))
MDS3_4

ggsave(MDS3_4, filename = "/Users/hannahaichelman/Documents/BU/TVE/2bRAD/Analysis/tuftscustompipeline_denovo_nosyms/MDS3_4_sitename.pdf", width=5, height=5, units=c("in"), useDingbats=FALSE)


#### Analyze relatedness ####
# read in relatedness matrix - this is the dataset with 50 individuals, no clones or technical replicates included.
rel=read.table("2bRAD/data_files/ngsrelate.noclone.res", header = TRUE)
head(rel)
# i2p file associated with these data (in same order as bams list so matches output of ngsrelate):
dim(i2p_noclones)

i2p_noclones_relatedness = i2p_noclones %>%
  mutate(merge_id = as.factor(seq(0,49,by=1)))

ggplot(rel, aes(x = a, y = rab)) +
  geom_point()

# read in relatedness matrix - this is the dataset with 51 individuals, no technical replicates included but one true clone pair left in.
rel.allsamps=read.table("2bRAD/data_files/ngsrelate.noclone.allsamps.res", header = TRUE)
str(rel.allsamps)
rel.allsamps$a = as.factor(rel.allsamps$a)
rel.allsamps$b = as.factor(rel.allsamps$b)

# i2p file associated with these data:
dim(i2p_noclones_allsamps)

# add column to allow us to merge with relatedness file
i2p_noclones_allsamps_relatedness = i2p_noclones_allsamps %>%
  mutate(merge_id = as.factor(seq(0,50,by=1)))

# merge metadata and relatedness file
rel_i2p_merged = left_join(rel.allsamps, i2p_noclones_allsamps_relatedness, by = c("a" = "merge_id"))
head(rel_i2p_merged)

# plot
ggplot(rel_i2p_merged, aes(x = sample_id, y = rab)) +
  geom_point() +
  labs(title = "Pairwise Relatedness")

# this way of plotting is a little muddy, will clean things up with a heatmap.
# subset data to only include sample names and the rab metric that we care about
rel_i2p_filt = rel_i2p_merged %>%
  select(a, b, rab, sample_id) %>%
  rename(sample_id_a = sample_id) %>%
  mutate(sample_id_a = as.factor(sample_id_a))

# merge with metadata again to get the sample_id_b
rel_i2p_filt2 = left_join(rel_i2p_filt, i2p_noclones_allsamps_relatedness, by = c("b" = "merge_id"))
head(rel_i2p_filt2)

# now final dataframe with sample_id_a and sample_id_b and relatedness
rel_i2p_filt_final = rel_i2p_filt2 %>%
  select(sample_id_a, sample_id, rab) %>%
  rename(sample_id_b = sample_id) %>%
  mutate(sample_id_b = as.factor(sample_id_b))

str(rel_i2p_filt_final)

# cut "_CLONE" off sample ID's to make reading and merging easier later on
rel_i2p_filt_final$sample_id_a = gsub("_CLONE", "", rel_i2p_filt_final$sample_id_a)
rel_i2p_filt_final$sample_id_b = gsub("_CLONE", "", rel_i2p_filt_final$sample_id_b)

# plot, again not the best way but just to peek.
rel.plot = ggplot(rel_i2p_filt_final, aes(x = sample_id_a, y = rab, text = sample_id_b)) +
  geom_point() +
  labs(title = "Pairwise Relatedness")
ggplotly(rel.plot)

# order this dataframe so that sample_id_a column is alphabetical
rel_i2p_filt_final2 = rel_i2p_filt_final[order(rel_i2p_filt_final$sample_id_a),]

# this looks good so far, but want to organize based on lineage so read that in
lineages = read.csv("Physiology_Data/data_files/tve_lineages_noclones.csv")
head(lineages)

# merge lineages with relatedness data
rel_i2p_filt_final_lineage = left_join(rel_i2p_filt_final, lineages, by = c("sample_id_a" = "gen_site")) %>%
  rename(lineage_a = lineage)
head(rel_i2p_filt_final_lineage)

# merge again to get lineage associated with sample_id_b
rel_i2p_filt_final_lineage2 = left_join(rel_i2p_filt_final_lineage, lineages, by = c("sample_id_b" = "gen_site")) %>%
  rename(lineage_b = lineage) %>%
  # remove any lineage na's
  filter(!is.na(lineage_a)) %>%
  filter(!is.na(lineage_b))

head(rel_i2p_filt_final_lineage2)

# create a new column that combines lineage_a and lineage_b with an '_'
rel_i2p_filt_final_lineage2$lineage_comparison = paste(rel_i2p_filt_final_lineage2$lineage_a, rel_i2p_filt_final_lineage2$lineage_b, sep = "_")
rel_i2p_filt_final_lineage2$lineage_comparison = as.factor(rel_i2p_filt_final_lineage2$lineage_comparison)
levels(rel_i2p_filt_final_lineage2$lineage_comparison)

# re-order lineage_comparison factor
rel_i2p_filt_final_lineage2$lineage_comparison = factor(rel_i2p_filt_final_lineage2$lineage_comparison, levels = c("L1_L1", "L1_L2", "L1_L3", "L2_L2", "L2_L1", "L2_L3", "L3_L3", "L3_L1", "L3_L2"))

# summarize relatedness by lineage comparison
rel_lin_summary = rel_i2p_filt_final_lineage2 %>%
  group_by(lineage_comparison) %>%
  summarize(mean_rab = mean(rab, na.rm = TRUE),
            sd_rab = sd(rab, na.rm = TRUE),
            n = n())

# make boxplot of relatedness by lineage comparison
comparison_of_interest = c("L1_L1", "L2_L2", "L3_L3")
# made sure I4G clone is not included in this plot - it isn't!

rel.boxplot = rel_i2p_filt_final_lineage2 %>%
  # select only certain lineage comparisons for plotting - uncomment this to look at all comparisons
  filter(lineage_comparison %in% comparison_of_interest) %>%
  ggplot(aes(x = lineage_comparison, y = rab,fill = lineage_comparison)) +
  geom_boxplot() +
  #scale_color_manual(values = c("#3f007d", "#807dba", "#bcbddc")) +
  scale_fill_manual(values = c("#3f007d", "#807dba", "#bcbddc")) +
  labs(y = "Pairwise Relatedness (rab)") +
  theme_bw()
rel.boxplot
ggsave(rel.boxplot, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/2bRAD/Analysis/tuftscustompipeline_denovo_nosyms/relatedness_lin.pdf", width=5, height=4, units=c("in"), useDingbats=FALSE)

# create a matrix for the heatmap
rel.mat = pivot_wider(data = rel_i2p_filt_final, names_from = sample_id_a, values_from = rab) %>%
  column_to_rownames(var = "sample_id_b")

rel.mat.num = as.matrix(rel.mat)
str(rel.mat.num)
corrplot(rel.mat.num,
         type="lower", # show only lower diagonal of corellogram
         tl.col="black",  # text color of data label
         #order = 'alphabet',
         tl.srt=45)

# summarize pairwise relatedness by lineage
rel_lin_summary = rel_i2p_filt_final_lineage %>%
  group_by(lineage) %>%
  summarize(mean_rab = mean(rab, na.rm = TRUE),
            sd_rab = sd(rab, na.rm = TRUE),
            n = n())


rel.plot = ggplot(rel_i2p_filt_final_lineage, aes(x = sample_id_a, y = rab, text = sample_id_b, color = lineage)) +
  geom_point() +
  labs(title = "Pairwise Relatedness by Lineage")
rel.plot

#--------------------
# covariance / PCA (not really needed, Misha prefers IBS)

# performing PCoA and CAP
library(vegan)
# load in covariance matrix
co = as.matrix(read.table("tuftscustompipeline_denovo_nosyms/myresult2.noclone.covMat")) # covariance based on single-read sampling

conds=data.frame(cbind(site))
pp0=capscale(as.dist(1-cov2cor(co))~1) # PCoA
pp=capscale(as.dist(1-cov2cor(co))~site,conds) # CAP

# significance of by-site divergence, based on 1-correlation as distance
adonis(as.dist(1-cov2cor(co))~site,conds)

# eigenvectors: how many are interesting?
plot(pp0$CA$eig)

axes2plot=c(1,2)
quartz()
library(adegenet) # for transp()
cmd=pp0
plot(cmd,choices=axes2plot,display="sites",type="n") # choices - axes to display
points(cmd,choices=axes2plot,pch=19,col=transp(colors,alpha=0.7))
#ordihull(cmd,choices= axes2plot,groups= conds$grp,draw="polygon",col=1+as.numeric(unique(as.factor(conds$grp))),label=T)
ordispider(cmd,choices= axes2plot,groups=conds$site,col="grey80")
ordiellipse(cmd,choices= axes2plot,groups= conds$site,draw="polygon",col=colpops,label=T)

# unscaled, to identify outliers
plot(cmd$CA$u[,axes2plot],pch=19,col=colors)
ordispider(cmd$CA$u[,axes2plot],groups=conds$site,col="grey80")
ordiellipse(cmd$CA$u[,axes2plot],groups= conds$site,draw="polygon",col=colpops,label=T)
identify(cmd$CA$u[,axes2plot],labels=colnames(ma),n=3,cex=0.7)

# I2I, I2D, I2H are the three outliers from I2 here - don't correspond to lowest or highest coverage in quality.txt file

# prettier ggplot option
axes2plot=c(1,2)
pca_s <- as.data.frame(cmd$CA$u[,axes2plot])
#colors=c('royalblue4','cornflowerblue','lightblue','red4','indianred3','mistyrose3')
cols_site <- c("CI" = "#543005", "PD"= "#bf812d",  "SP"= "#dfc27d",  "BN" = "#003c30", "BS"= "#35978f", "CA"= "#80cdc1")

MDS1_2 = ggplot(pca_s, aes(MDS1, MDS2)) +
  theme_classic() +
  geom_point(aes(colour=conds$site), size=2, stroke = 1)  +
  scale_color_manual(values = cols_site,
                     breaks = c("BN", "BS", "CA", "CI", "PD", "SP"),
                     name = "Site") +
  stat_ellipse(geom="polygon", alpha = .3,
               aes(fill = conds$site), show.legend = NA) +
  scale_fill_manual(values=cols_site,
                    breaks = c("BN", "BS", "CA", "CI", "PD", "SP"),
                    name = "Site") +
  xlab(paste0("MDS1")) +
  ylab(paste0("MDS2"))
MDS1_2

ggsave(MDS1_2, filename = "/Users/hannahaichelman/Documents/BU/TVE/2bRAD/Analysis/tuftscustompipeline_denovo_nosyms/MDS1_2_sitename.pdf", width=5, height=5, units=c("in"), useDingbats=FALSE)


#plotting 2 and 3
axes2plot=c(3,4)
pca_s <- as.data.frame(cmd$CA$u[,axes2plot])

MDS3_4 = ggplot(pca_s, aes(MDS3, MDS4)) +
  theme_classic() +
  geom_point(aes(colour=conds$site), size=2, stroke = 1)  +
  scale_color_manual(values = colors,
                     breaks = c("BN", "BS", "CA", "CI", "PD", "SP"),
                     name = "Site") +
  stat_ellipse(geom="polygon", alpha = .2,
               aes(fill = conds$site), show.legend = NA) + scale_fill_manual(values=colors, name = "Site") +
  xlab(paste0("MDS3")) +
  ylab(paste0("MDS4"))
MDS3_4

ggsave(MDS3_4, filename = "/Users/hannahaichelman/Documents/BU/TVE/2bRAD/Analysis/tuftscustompipeline_denovo_nosyms/MDS3_4_sitename.pdf", width=5, height=5, units=c("in"), useDingbats=FALSE)



#-------------
# t-SNE:  machine learning to identify groups of samples
# based on genotypes' correlations
# (only makes sense if you have hundreds of samples)

library(Rtsne)
library(vegan)
library(adegenet)
quartz()

# perplexity:  expected number fo neighbors. Set to 0.5x N(samples per pop)
perp=15
rt = Rtsne(as.dist(1-cov2cor(co)), perplexity=perp,max_iter=2,is_distance=T)
for (i in 1:250){
  rt = Rtsne(as.dist(1-cov2cor(co)), perplexity=perp,max_iter=10,Y_init=rt$Y,is_distance=T)
  plot(rt$Y,col=colors,pch=16,cex=0.8,main=i*10)
}
ordispider(rt$Y,groups=site,col="grey80",alpha=0.01)
ordiellipse(rt$Y,groups= site,draw="polygon",col=colpops,label=T)





