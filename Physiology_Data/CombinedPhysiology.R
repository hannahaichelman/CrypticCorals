# Written by Hannah Aichelman (hannahaichelman@gmail.com)
# Physiology Analysis

library(ggplot2)
library(lme4)
library(plotly)
library(ggridges)
library(tidyverse)
library(arsenal) #easily compare two data frames
library(Rmisc)
library(lmerTest)
library(emmeans)
library(reshape)
library(readxl)
library(wesanderson)
library(ggpubr)
library(car)
library(sjPlot)
library(effects)
library(glmmTMB)
library(performance)
library(patchwork)
library(magrittr)
library(SciViews)

##### Read in and format data #####
# set wd
setwd("/Users/hannahaichelman/Documents/BU/TVE")
# read in the data
post_phys <- read.csv('dtvmaster.csv') # physiology data taken at the end of the experiment
init_phys <- read.csv('initial-phys-mod.csv') # physiology data taken at the start of the experiment

# did it load correctly?
head(post_phys)
str(post_phys)

# not reading in calcification or fv/fm data here because we don't have that information for the initial physiology.
# so i am only reading in the data that matches across the initial and final physiology data frames so we can combine them effectively
post_phys = post_phys %>%
  dplyr::select(frag, survivedtoend, treat, t3sastan1, t3sastan2, t3sastan3, t3sarec1, t3sarec2, t3sarec3, blastvol, blaster,
                hprotplate, hprot1, hprot2, hprot3, chl630_1, chl630_2, chl630_3, chl663_1, chl663_2, chl663_3, symcount1, symcount2, symcount3) %>%
  dplyr::rename(sastan1 = t3sastan1, sastan2 = t3sastan2, sastan3 = t3sastan3, sarec1 = t3sarec1, sarec2 = t3sarec2, sarec3 = t3sarec3) %>%
  mutate(treat = as.factor(treat), sarec3 = as.numeric(sarec3), blastvol = as.numeric(blastvol)) %>%
  dplyr::filter(frag!="I2A12", frag!="I2A7", frag!="O4D8", frag!="I3E10", frag!="I3I10") #these frags are being removed because they were duplicate genotypes within treatment and had the least complete information of the two

init_phys = init_phys %>%
  dplyr::select(frag, survivedtoend, treat, pasastan1, pasastan2, pasastan3, pasarec1, pasarec2, pasarec3, blastvol, blaster,
                hprotplate, hprot1, hprot2, hprot3, chl630_1, chl630_2, chl630_3, chl663_1, chl663_2, chl663_3, symcount1, symcount2, symcount3) %>%
  dplyr::rename(sastan1 = pasastan1, sastan2 = pasastan2, sastan3 = pasastan3, sarec1 = pasarec1, sarec2 = pasarec2, sarec3 = pasarec3) %>%
  mutate(treat = as.factor(treat), survivedtoend = as.character(survivedtoend))

dim(post_phys)
dim(init_phys)

combined_phys <- dplyr::bind_rows(post_phys,init_phys)
dim(combined_phys)
head(combined_phys)

# add in the carbohydrate data we are confident in that Olivia Nieves worked on
carbs <- read.csv('/Users/hannahaichelman/Documents/BU/TVE/Carbohydrates/host_sym_carbs_ON.csv')
dim(carbs)

# merge data frames to keep all of the individuals in the combined physiology spreadsheet and add NA's for any that we don't have carb info for
phys <- left_join(combined_phys, carbs, by = "frag")
dim(phys)
head(phys)

# extract the original two letter site code (inshore/offshore #1-4) from the nubbin ID
phys$origsitecode <- substr(phys$frag, 1, 2)
phys$origsitecode <- as.factor(phys$origsitecode)

# add in site name with inshore/offshore indicator
phys$sitename <- ifelse(phys$origsitecode == 'I2', 'SP',
                        ifelse(phys$origsitecode == 'I3', 'CI',
                               ifelse(phys$origsitecode == 'I4', 'PD',
                                      ifelse(phys$origsitecode == 'O2', 'BS',
                                             ifelse(phys$origsitecode == 'O3', 'CA',
                                                    'BN')))))
phys$sitename <- as.factor(phys$sitename)

# make new nubbin IDs based on the new, more informative site codes
# extract the genotype and frag number
phys$fragid <- substr(phys$frag,3,5)

# add inshore/offshore designation
phys$reef <- substr(phys$frag,1,1)
phys$reef <- ifelse(phys$reef == 'O', 'Outer Reef', 'Inner Reef')
phys$reef <- as.factor(phys$reef)

phys$genet <- substr(phys$fragid,1,1)

#create a new column of combined genotype and site for stats later
phys = phys %>%
  unite(gen_site, c(origsitecode,genet), sep = "", remove = FALSE) %>%
  mutate(gen_site = as.factor(gen_site))

# re-level and re-name treatment
phys$treat <- factor(phys$treat, levels = c("init", "1", "2", "3","4","5"))
levels(phys$treat) <- c("Initial","Control","Low Var","Mod Var","High Var","Control 2")
phys$treat <- as.factor(phys$treat)

str(phys)

## combine physiology data with dominant symbiont type data
head(phys)

its2_types = read.csv("/Users/hannahaichelman/Documents/BU/TVE/16S_ITS2/ITS_PreStress_Timepoint/ITS2.dominanttype.csv")
head(its2_types)

its2_types = its2_types %>%
  select(-Sample_or_Control, -treat, -gen_site, -sitename, -reef, -lineage)

phys <- left_join(phys, its2_types, by = "frag")
str(phys)
phys$dominant_type = as.factor(phys$dominant_type)
phys$minor_type = as.factor(phys$minor_type)
head(phys)

# merge with its2 divs
its2_divs = read.csv("/Users/hannahaichelman/Documents/BU/TVE/16S_ITS2/ITS_PreStress_Timepoint/ITS2.dominantDIVs.csv") %>%
  select(frag, dominant_div)

phys <- left_join(phys, its2_divs, by = "frag")

# merge with lineage
# this includes data frame without duplicated preps and with no clones (I4G removed because it was smaller bam file)
lineages = read.csv("/Users/hannahaichelman/Documents/BU/TVE/2bRAD/Analysis/tuftscustompipeline_denovo_nosyms/tve_lineages_noclones.csv")
head(lineages)

phys <- left_join(phys, lineages, by = "gen_site")
phys$lineage = as.factor(phys$lineage)
head(phys)

## make metadata to use throughout the code

phys_metadata = phys %>%
  select(frag, treat, gen_site, sitename, reef, dominant_div, dominant_type, lineage)

head(phys_metadata)
#write.csv(phys_metadata_all, "/Users/hannahaichelman/Documents/BU/TVE/phys_metadata_updatedDIVs.csv",row.names=FALSE)

# set color palettes used throughout
cols_site <- c("CI" = "#543005", "PD"= "#bf812d",  "SP"= "#dfc27d",  "BN" = "#003c30", "BS"= "#35978f", "CA"= "#80cdc1")
cols_treat <- c("darkgrey", "#FF9966","#CC3300","#7f0000")
cols_lineage <- c("L1" = "#3f007d", "L2" = "#807dba", "L3" = "#bcbddc")
its2_cols_greens = c("C1" = "#edf8e9", "C3af" = "#238b45","C3" = "#a1d99b","D1" = "#00441b")

#### Test for Lineage Distribution ####
# Want to test whether distribution of lineages is significantly different across inshore and offshore sites.

str(phys_metadata_all) # includes all 3 lineages

# first try logit regression: https://stats.oarc.ucla.edu/r/dae/logit-regression/
library(aod)

mylogit = glm(formula = lineage ~ reef, data = phys_metadata_all, family = "binomial")
summary(mylogit)

# now try chi-square test:
chisq_dat = phys_metadata_all %>%
  distinct(gen_site, .keep_all = TRUE)

table(chisq_dat$lineage, chisq_dat$reef) # make contingency table

test = chisq.test(table(chisq_dat$lineage, chisq_dat$reef))
test
# Warning message:
#   In chisq.test(table(chisq_dat$lineage, chisq_dat$reef)) :
#   Chi-squared approximation may be incorrect

# Because chisq.test requires an assumption of a minimum of 5 expected counts in the contingency table, we will use a Fisher's exact test.

fisher.test = fisher.test(table(chisq_dat$lineage, chisq_dat$reef))
fisher.test

# Fisher's Exact Test for Count Data
#
# data:  table(chisq_dat$lineage, chisq_dat$reef)
# p-value = 8.707e-08
# alternative hypothesis: two.sided

# Or can remove low count L3 from dataset and still use Chi-Square test - THIS IS REPORTED IN MANUSCRIPT
chisq_dat_noL3 = phys_metadata_all %>%
  distinct(gen_site, .keep_all = TRUE) %>%
  filter(lineage != "L3") %>%
  droplevels()

test = chisq.test(table(chisq_dat_noL3$lineage, chisq_dat_noL3$reef))
test

# Pearson's Chi-squared test with Yates' continuity correction
#
# data:  table(chisq_dat_noL3$lineage, chisq_dat_noL3$reef)
# X-squared = 23.577, df = 1, p-value = 1.2e-06

#### Surface Area Measurements ####
# get surface area measurements in cm^2
# pa (pre-acclimation) measurements for surface area are suitable only for initial phys nubbins. physiology for the experimental nubbins must use T3 surface area
# I used the appropriate surface area measurement for each time point, but combined into one column for simplicity
# look above for more details - read in initial phys and post-variability phys separately, each data frame with its own SA measurements included

phys$SAcm2 <- (4*((phys$sarec1+phys$sarec2+phys$sarec3)/3))/((phys$sastan1+phys$sastan2+phys$sastan3)/3)

# check to see if surface areas look normal
# histogram of nubbin sizes
ggplot(phys, aes(x=SAcm2, color = sitename, fill = sitename))+
  theme_bw()+
  geom_histogram(binwidth = 1)+
  scale_fill_manual(values = cols_site)+
  scale_color_manual(values = cols_site)+
  labs(x=expression(paste("Surface Area (cm"^2*')')))+
  facet_wrap(~sitename, 6,1)

# boxplot of nubbin sizes
plot <- ggplot(phys, aes(x=sitename, y=SAcm2, fill = sitename, text = paste("Nubbin:",frag)))+
  theme_bw()+
  geom_boxplot()+
  geom_jitter(aes(color = sitename),width = 0.3)+
  scale_fill_manual(values = cols_site)+
  scale_color_manual(values = cols_site)+
  theme(axis.title.x = element_blank())+
  theme(axis.text.x=element_text(angle = 45, vjust = 0.5))
#labs(y=expression(paste("Surface Area (cm"^2*')')))
ggplotly(subplot(list(plot),nrows=1,titleY=F) %>% layout(showlegend=T))

##### Tissue Thickness (T0) #####

# tissue thickness is measured in mm and was only taken from nubbins from the start of the experiment (T0)

# read in initial physiology data file again
init_phys_tiss <- read.csv('initial-phys-mod.csv')

tiss_phys = init_phys_tiss %>%
  select(frag,treat,tisthick1,tisthick2,tisthick3,tisthick4,tisthick5,tisthick6,
         origsitecode,sitename,fragid,nubbin,reef,SAcm2) %>%
  mutate(treat = as.factor(treat), reef = as.factor(reef), sitename = as.factor(sitename))

# re-name treatment
levels(tiss_phys$treat) <- c("Initial")

#re-name site names
tiss_phys$sitename <- factor(tiss_phys$sitename, levels = c("I-Cristobal","I-Punta Donato","I-STRI Point",
                                                            "O-Bastimentos N","O-Bastimentos S","O-Cayo de Agua"))
levels(tiss_phys$sitename) <- c("CI","PD","SP","BN","BS","CA")

# add identifying data
# make new nubbin IDs based on the new, more informative site codes
tiss_phys$genet <- substr(tiss_phys$fragid,1,1)

#create a new column of combined genotype and site for stats later
# calculate average tissue thickness for each nubbin
tiss_phys = tiss_phys %>%
  unite(gen_site, c(origsitecode,genet), sep = "", remove = FALSE) %>%
  mutate(gen_site = as.factor(gen_site)) %>%
  mutate(avgtiss=rowMeans(.[ , c("tisthick1","tisthick2","tisthick3","tisthick4","tisthick5","tisthick6")], na.rm=TRUE)) %>%
  dplyr::filter(frag!= "I3E6_v2") %>% # redundant info, same as other I3E sample
  dplyr::filter(gen_site != "I4G") # clone pair, removed this sample from 2brad so removing from phys too


# merge with lineage info for later plotting
# this includes data frame without duplicated preps, and with one of the clone pair (I4G) removed
lineages = read.csv("/Users/hannahaichelman/Documents/BU/TVE/2bRAD/Analysis/tuftscustompipeline_denovo_nosyms/tve_lineages_noclones.csv")

tiss_phys_all_lin <- left_join(tiss_phys, lineages, by = "gen_site")

str(tiss_phys_all_lin)
tiss_phys_all_lin$gen_site = as.factor(tiss_phys_all_lin$gen_site)
tiss_phys_all_lin$lineage = as.factor(tiss_phys_all_lin$lineage)

tiss_phys_2_lin = tiss_phys_all_lin %>%
  dplyr::filter(is.na(lineage) | lineage!="L3") # want to keep NA values for lineage here since they still have other info, will remove na's for lineage specific plots

# can't combine with dominant sym type here because tissue thickness is from initial phys and sym types were from end of variability

# STATS
m1 <- lm(avgtiss ~ lineage+sitename, data = tiss_phys_all_lin)
summary(m1)
anova(m1)
# Analysis of Variance Table
#
# Response: avgtiss
#           Df  Sum Sq Mean Sq F value    Pr(>F)
# lineage    2 10.4316  5.2158  9.2356 0.0005385 ***
# sitename   5  6.5677  1.3135  2.3259 0.0614993 .
# Residuals 38 21.4605  0.5647

lsmeans(m1, pairwise~lineage, adjust="tukey")
# $contrasts
# contrast estimate    SE df t.ratio p.value
# L1 - L2   0.00459 0.404 38   0.011  0.9999
# L1 - L3   1.20236 0.681 38   1.764  0.1952
# L2 - L3   1.19778 0.549 38   2.182  0.0872

# PLOTS

#SummarySE to format data for plotting
tiss_means_site_all_lin <- summarySE(tiss_phys_all_lin, measurevar="avgtiss", groupvars=c("treat","sitename"))
tiss_means_site_2_lin <- summarySE(tiss_phys_2_lin, measurevar="avgtiss", groupvars=c("treat","sitename"))

# plot, treatment x axis colored by site data figure
tiss_plot_site <- ggplot(tiss_means_site_2_lin,aes(x = sitename, y = avgtiss, fill = sitename, color = sitename))+
  theme_bw()+
  geom_errorbar(aes(x = sitename, ymax = avgtiss+se, ymin = avgtiss-se, color = sitename), width = .2, position = position_dodge(width=0.3)) +
  geom_point(size = 3, pch = 21, position = position_dodge(width=0.3), color = "black")+
  scale_fill_manual(name = "Site",
                     values = cols_site)+
  scale_color_manual(name = "Site",
                     values = cols_site)+
  xlab("Site Name")+
  ylab("Tissue Thickness (mm)")+
  ylim(3.75,6) +
  #geom_vline(xintercept = 1.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  theme(legend.position = "none")
tiss_plot_site

ggsave(tiss_plot_site, filename = "/Users/hannahaichelman/Documents/BU/TVE/TissueThickness/tissuethickness_site_2_lin.pdf", width=4, height=4, units=c("in"), useDingbats=FALSE)


#SummarySE to format data for plotting with lineage
tiss_phys_all_lin_nona = tiss_phys_all_lin %>%
  drop_na(lineage)

tiss_phys_2_lin_nona = tiss_phys_2_lin %>%
  drop_na(lineage)

tiss_means_all_lin <- summarySE(tiss_phys_all_lin_nona, measurevar="avgtiss", groupvars=c("lineage"))
tiss_means_2_lin <- summarySE(tiss_phys_2_lin_nona, measurevar="avgtiss", groupvars=c("lineage"))

# plot, lineage on x axis
tiss_plot_lineage <- ggplot(tiss_means_all_lin,aes(x = lineage, y = avgtiss, fill = lineage, color = lineage))+
  theme_bw()+
  geom_errorbar(aes(x = lineage, ymax = avgtiss+se, ymin = avgtiss-se, color = lineage), width = .2, position = position_dodge(width=0.3)) +
  geom_point(size = 3, position = position_dodge(width=0.3), shape = 21, color = "black")+
  scale_fill_manual(name = "Lineage",
                    values = cols_lineage)+
  scale_color_manual(name = "Lineage",
                     values = cols_lineage)+
  xlab("Lineage")+
  ylab("Tissue Thickness (mm)")+
  #ylim(3.75,6) +
  #geom_vline(xintercept = 1.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  theme(legend.position = "none")
tiss_plot_lineage

ggsave(tiss_plot_lineage, filename = "/Users/hannahaichelman/Documents/BU/TVE/TissueThickness/tissuethickness_lineage_alllin.pdf", width=4, height=4, units=c("in"), useDingbats=FALSE)

# plot just CI corals to see if lineage difference holds
tiss_phys_lineage_CI = tiss_phys_2_lin_nona %>%
  subset(sitename == "CI")

tiss_means_lineage_CI <- summarySE(tiss_phys_lineage_CI, measurevar="avgtiss", groupvars=c("lineage"))

m1 <- lm(avgtiss ~ lineage, data = tiss_phys_lineage_CI)
summary(m1)
anova(m1)
# Response: avgtiss
# Df Sum Sq Mean Sq F value Pr(>F)
# lineage    1 0.0764 0.07640  0.2372 0.6469
# Residuals  5 1.6107 0.32215

lsmeans(m1, pairwise~lineage, adjust="tukey")

# plot, lineage on x axis
tiss_plot_lineage_CI <- ggplot(tiss_means_lineage_CI,aes(x = lineage, y = avgtiss, color = lineage))+
  theme_bw()+
  geom_point(size = 3, position = position_dodge(width=0.3))+
  geom_errorbar(aes(x = lineage, ymax = avgtiss+se, ymin = avgtiss-se), width = .2, position = position_dodge(width=0.3)) +
  scale_color_manual(name = "Lineage",
                     labels = c("L1","L2"),
                     values = cols_lineage)+
  xlab("Lineage")+
  ylab("Tissue Thickness (mm)")+
  #ylim(3.75,6) +
  #geom_vline(xintercept = 1.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  theme(legend.position = "none")
tiss_plot_lineage_CI

ggsave(tiss_plot_lineage_CI, filename = "/Users/hannahaichelman/Documents/BU/TVE/TissueThickness/tissuethickness_lineage_CIonly.pdf", width=4, height=4, units=c("in"), useDingbats=FALSE)


# plot just SP corals to see if lineage difference holds
tiss_phys_lineage_SP = tiss_phys_all_lin_nona %>%
  subset(sitename == "SP")

tiss_means_lineage_SP <- summarySE(tiss_phys_lineage_SP, measurevar="avgtiss", groupvars=c("lineage"))

# plot, reef zone x axis
tiss_plot_lineage_SP <- ggplot(tiss_means_lineage_SP,aes(x = lineage, y = avgtiss, color = lineage))+
  theme_bw()+
  geom_point(size = 3, position = position_dodge(width=0.3))+
  geom_errorbar(aes(x = lineage, ymax = avgtiss+se, ymin = avgtiss-se), width = .2, position = position_dodge(width=0.3)) +
  scale_color_manual(name = "Lineage",
                     labels = c("L2","L3"),
                     values = cols_lineage)+
  xlab("Lineage")+
  ylab("Tissue Thickness (mm)")+
  #ylim(3.75,6) +
  #geom_vline(xintercept = 1.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  theme(legend.position = "none")
tiss_plot_lineage_SP

ggsave(tiss_plot_lineage_SP, filename = "/Users/hannahaichelman/Documents/BU/TVE/TissueThickness/tissuethickness_lineage_SPonly.pdf", width=4, height=4, units=c("in"), useDingbats=FALSE)

##### Protein Concentrations #####
# first, subset physiology data add in standard curve data needed to correct for separate protein plates
hprot_phys = phys %>%
  select(frag,treat,hprotplate,hprot1,hprot2,hprot3,
         survivedtoend,blastvol,gen_site,origsitecode,sitename,fragid,
         reef,genet,SAcm2, dominant_type, lineage)

head(hprot_phys)

# we need to adjust the mean absorbance values by the intercept value of each plate
hprot_phys$protx2coef <- as.numeric(ifelse(hprot_phys$hprotplate == 'D1', 0.472799,
                                           ifelse(hprot_phys$hprotplate == 'D2', 0.385059,
                                                  ifelse(hprot_phys$hprotplate == 'D3', 0.618482,
                                                         ifelse(hprot_phys$hprotplate == 'D4', 0.539057,
                                                                ifelse(hprot_phys$hprotplate == 'D5', 0.637505,
                                                                       ifelse(hprot_phys$hprotplate == 'D6', 0.4492627,
                                                                              ifelse(hprot_phys$hprotplate == 'D7', 0.456354,
                                                                                     ifelse(hprot_phys$hprotplate == 'D8', 0.526895,
                                                                                            ifelse(hprot_phys$hprotplate == 'O1', 0.0650601,
                                                                                                   ifelse(hprot_phys$hprotplate == 'O2', 0.753409,
                                                                                                          ifelse(hprot_phys$hprotplate == 'O3', 0.4325502,
                                                                                                                 ifelse(hprot_phys$hprotplate == 'O4', 0.554534,
                                                                                                                        ifelse(hprot_phys$hprotplate == 'O5', 0.680815,
                                                                                                                               ifelse(hprot_phys$hprotplate == 'O6', 0.772513,
                                                                                                                                      'NA')))))))))))))))

hprot_phys$protxcoef <- as.numeric(ifelse(hprot_phys$hprotplate == 'D1', -0.061659,
                                          ifelse(hprot_phys$hprotplate == 'D2', 0.038351,
                                                 ifelse(hprot_phys$hprotplate == 'D3', -0.073079,
                                                        ifelse(hprot_phys$hprotplate == 'D4', -0.048588,
                                                               ifelse(hprot_phys$hprotplate == 'D5', -0.029648,
                                                                      ifelse(hprot_phys$hprotplate == 'D6', 0.0452952,
                                                                             ifelse(hprot_phys$hprotplate == 'D7', -0.014998,
                                                                                    ifelse(hprot_phys$hprotplate == 'D8', -0.031275,
                                                                                           ifelse(hprot_phys$hprotplate == 'O1', 0.0777724,
                                                                                                  ifelse(hprot_phys$hprotplate == 'O2', -0.021612,
                                                                                                         ifelse(hprot_phys$hprotplate == 'O3', 0.0287098,
                                                                                                                ifelse(hprot_phys$hprotplate == 'O4', -0.029577,
                                                                                                                       ifelse(hprot_phys$hprotplate == 'O5', -0.027620,
                                                                                                                              ifelse(hprot_phys$hprotplate == 'O6', -0.048054,
                                                                                                                                     'NA')))))))))))))))

hprot_phys$protint <- as.numeric(ifelse(hprot_phys$hprotplate == 'D1', 0.004528,
                                        ifelse(hprot_phys$hprotplate == 'D2',0.002527,
                                               ifelse(hprot_phys$hprotplate == 'D3', 0.004711,
                                                      ifelse(hprot_phys$hprotplate == 'D4', 0.003552,
                                                             ifelse(hprot_phys$hprotplate == 'D5', 0.003320,
                                                                    ifelse(hprot_phys$hprotplate == 'D6', 0.0018613,
                                                                           ifelse(hprot_phys$hprotplate == 'D7', 0.002955,
                                                                                  ifelse(hprot_phys$hprotplate == 'D8', 0.003667,
                                                                                         ifelse(hprot_phys$hprotplate == 'O1', 0.0005224,
                                                                                                ifelse(hprot_phys$hprotplate == 'O2', 0.002474,
                                                                                                       ifelse(hprot_phys$hprotplate == 'O3', 0.0007926,
                                                                                                              ifelse(hprot_phys$hprotplate == 'O4', 0.002101,
                                                                                                                     ifelse(hprot_phys$hprotplate == 'O5', 0.003050,
                                                                                                                            ifelse(hprot_phys$hprotplate == 'O6', 0.002921,
                                                                                                                                   'NA')))))))))))))))

head(hprot_phys)

# now we can use the quadratic equation to solve for absorbance
# y = ax^2 + bx + c where
# y = protein concentration
# x = raw absorbance from Bradford Assay

# these equations get the protein values in ug/ul
hprot_phys$calcprot1 <- ((hprot_phys$protx2coef*(hprot_phys$hprot1^2)) + (hprot_phys$protxcoef*hprot_phys$hprot1) + hprot_phys$protint)
hprot_phys$calcprot2 <- ((hprot_phys$protx2coef*(hprot_phys$hprot2^2)) + (hprot_phys$protxcoef*hprot_phys$hprot2) + hprot_phys$protint)
hprot_phys$calcprot3 <- ((hprot_phys$protx2coef*(hprot_phys$hprot3^2)) + (hprot_phys$protxcoef*hprot_phys$hprot3) + hprot_phys$protint)

# calculate average protein concentration
hprot_phys = hprot_phys %>%
  mutate(avgprot=rowMeans(.[ , c("calcprot1","calcprot2","calcprot3")], na.rm=TRUE))

# [prot concentration in ug/ul] = ([sample protein value] * [1000 ul in 1 ml / 80 ul of sample volume (where sample volume includes dilution with seawater)] * number of mls in slurry / [surface area of fragment])
# 80/6.4 is the dilution factor - 80 uL total volume with 6.4 uL of sample
hprot_phys$prot_ugcm2 <- (hprot_phys$avgprot*(80/6.4)*1000*hprot_phys$blastvol)/hprot_phys$SAcm2
hprot_phys$prot_mgcm2 <- hprot_phys$prot_ugcm2/1000

# How many protein NA's?
sum(is.na(hprot_phys$prot_mgcm2))
# 110

# But how many of these are missing and how many are just nubbins that died?
missing_prot_inds = hprot_phys %>%
  filter(survivedtoend=="yes") %>%
  filter(is.na(hprot1))
# 5 have no protein record from the absorbance stage: I2C3, I3A7, O2C8, I3B4, I4E7, but all of these have no blaster or blastvol info

missing_blastinfo = hprot_phys %>%
  filter(survivedtoend=="yes") %>%
  filter(is.na(blastvol))
# "O3F2" "O3H2" "I2C3" "I3A7" "O2C8" "I4G3" "O2G7" "I3B4" "I3D8" "I3F4" "I4E7" "O4F4" "I4G4" "I3E5" "I3F5" "I3G5"

# make data subsets for stats and plotting
hprot_phys_all_lin = hprot_phys %>%
  drop_na(prot_mgcm2) %>%
  dplyr::filter(treat!="Control 2") %>%
  dplyr::filter(frag != "I3E6_v2") %>% # no consensus on this I3E6 individual phenotype for this measure - so removing both
  dplyr::filter(frag != "I3E6") %>%
  dplyr::filter(gen_site != "I4G") %>% # clone with I4F, remove from dataset
  select(frag, gen_site, treat, sitename, dominant_type, lineage, prot_mgcm2)

hprot_phys_2_lin = hprot_phys_all_lin %>%
  dplyr::filter(is.na(lineage) | lineage!="L3") %>%
  drop_na(sitename, lineage) # doing this for plotting

hprot_phys_2_lin_symtype = hprot_phys_2_lin %>%
  drop_na(treat, dominant_type)  # we only have symbiont types from the end of the experiment, so do separate filtering for those plots

initial_hprot_phys_all_lin = hprot_phys_all_lin %>%
  dplyr::filter(treat=="Initial")

initial_hprot_phys_2_lin = initial_hprot_phys_all_lin %>%
  dplyr::filter(is.na(lineage) | lineage!="L3") # want to keep NA values for lineage here since they still have other info, will remove na's for lineage specific plots

post_hprot_phys_all_lin = hprot_phys_all_lin %>%
  dplyr::filter(treat!="Initial")

post_hprot_phys_2_lin = post_hprot_phys_all_lin %>%
  dplyr::filter(is.na(lineage) | lineage!="L3")

# the values we have here seem to be similar to protein values from Wright et al. GCB, good sanity check.

# STATS
## NEED TO RUN STATS SEPARATELY FOR INITIAL AND FINAL TIME POINTS

model.prot.initial = lm(prot_mgcm2 ~ lineage+sitename, data = initial_hprot_phys_2_lin)
summary(model.prot.initial)
anova(model.prot.initial)
#           Df   Sum Sq  Mean Sq F value   Pr(>F)
# lineage    1 0.078953 0.078953 12.6405 0.001105 **
# sitename   5 0.044532 0.008906  1.4259 0.239232
# Residuals 35 0.218611 0.006246

check_model(model.prot.initial)


model.prot.final <- lmer(prot_mgcm2 ~ treat+lineage+sitename+dominant_type + (1|gen_site), data = post_hprot_phys_2_lin)
summary(model.prot.final)
anova(model.prot.final)
# Type III Analysis of Variance Table with Satterthwaite's method
#                 Sum Sq   Mean Sq NumDF  DenDF F value Pr(>F)
# treat         0.018050 0.0060165     3 90.664  0.6453 0.5879
# lineage       0.028737 0.0287374     1 23.451  3.0824 0.0922 .
# sitename      0.093492 0.0186985     5 24.719  2.0056 0.1130
# dominant_type 0.020250 0.0067499     3 70.411  0.7240 0.5410

check_model(model.prot.final)


#specify model (because we are interested in pairwise, have to include the interaction)
m.emm<- lmer(prot_mgcm2 ~ treat*lineage + (1|gen_site), data = post_hprot_phys_2_lin, REML=FALSE)

emms<-emmeans(m.emm, ~lineage|treat) #, adjust="Bonferoni"
pairs(emms, interaction = "pairwise") %>% rbind(adjust="fdr")
# treat    lineage_pairwise estimate     SE  df t.ratio p.value
# Control  L1 - L2            0.1277 0.0364 151   3.510  0.0024
# Low Var  L1 - L2            0.0678 0.0380 152   1.785  0.0762
# Mod Var  L1 - L2            0.0978 0.0341 149   2.873  0.0093
# High Var L1 - L2            0.0932 0.0359 151   2.600  0.0137


# PLOTS
#SummarySE to format data for plotting
prot_means_site_2_lin <- summarySE(hprot_phys_2_lin, measurevar="prot_mgcm2", groupvars=c("treat","sitename"))

# plot, treatment x axis colored by site data figure
prot_plot_site <- ggplot(prot_means_site_2_lin,aes(x = treat, y = prot_mgcm2, fill = sitename))+
  theme_bw()+
  geom_errorbar(aes(x = treat, ymax = prot_mgcm2+se, ymin = prot_mgcm2-se), width = .2, position = position_dodge(width=0.3)) +
  geom_point(size = 3, pch = 21, position = position_dodge(width=0.3))+
  scale_fill_manual(name = "Site",
                     values = cols_site)+
  xlab("Treatment")+
  ylab(bquote("Total Protein (mg" ~cm^-2~')'))+
  ylim(0,0.4) +
  geom_vline(xintercept = 1.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        legend.position = c(.95, .99),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(5, 5, 5, 5)
        #legend.key = element_rect(fill = "none")
)
prot_plot_site

ggsave(prot_plot_site, filename = "/Users/hannahaichelman/Documents/BU/TVE/Protein/plots/protein_site_2_lin.pdf", width=5, height=4, units=c("in"), useDingbats=FALSE)


# plot, treatment x axis and lineage color
#SummarySE to format data for plotting
prot_means_lineage_2_lin <- summarySE(hprot_phys_2_lin, measurevar="prot_mgcm2", groupvars=c("treat","lineage"))

# plot, treatment x axis colored by site data figure
prot_plot_lineage <- ggplot(prot_means_lineage_2_lin,aes(x = treat, y = prot_mgcm2, fill = lineage))+
  theme_bw()+
  geom_errorbar(aes(x = treat, ymax = prot_mgcm2+se, ymin = prot_mgcm2-se), width = .2, position = position_dodge(width=0.3)) +
  geom_point(size = 3, pch = 21, position = position_dodge(width=0.3))+
  scale_fill_manual(name = "Lineage",
                     values = cols_lineage)+
  xlab("Treatment")+
  ylab(bquote("Total Protein (mg" ~cm^-2~')'))+
  ylim(0,0.4) +
  geom_vline(xintercept = 1.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        legend.position = c(.95, .99),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(5, 5, 5, 5)
        #legend.key = element_rect(fill = "none")
)
prot_plot_lineage

ggsave(prot_plot_lineage, filename = "/Users/hannahaichelman/Documents/BU/TVE/Protein/plots/protein_lineage_2_lin.pdf", width=5, height=4, units=c("in"), useDingbats=FALSE)


# plot, treatment x axis and symtype color
#SummarySE to format data for plotting
prot_means_sym <- summarySE(hprot_phys_2_lin_symtype, measurevar="prot_mgcm2", groupvars=c("treat","dominant_type"))

# plot, treatment x axis colored by site data figure
prot_plot_sym <- ggplot(prot_means_sym,aes(x = treat, y = prot_mgcm2, fill = dominant_type))+
  theme_bw()+
  geom_errorbar(aes(x = treat, ymax = prot_mgcm2+se, ymin = prot_mgcm2-se), width = .2, position = position_dodge(width=0.3)) +
  geom_point(size = 3, pch = 21, position = position_dodge(width=0.3))+
  scale_fill_manual(name = "Dominant Symbiont",
                     values = its2_cols_greens)+
  xlab("Treatment")+
  ylab(bquote("Total Protein (mg" ~cm^-2~')'))+
  ylim(0,0.4) +
  #geom_vline(xintercept = 1.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        legend.position = c(.95, .99),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(5, 5, 5, 5)
        #legend.key = element_rect(fill = "none")
)
prot_plot_sym

ggsave(prot_plot_sym, filename = "/Users/hannahaichelman/Documents/BU/TVE/Protein/plots/protein_lineage_sym.pdf", width=6, height=4, units=c("in"), useDingbats=FALSE)



# Combine these protein plots
prot_plots_all = ggarrange(prot_plot_site, prot_plot_lineage, prot_plot_sym,
                           ncol = 3, nrow = 1)
ggsave(prot_plots_all, filename = "/Users/hannahaichelman/Documents/BU/TVE/Protein/plots/protein_all_plots.pdf", width=12, height=4, units=c("in"), useDingbats=FALSE)



# plot just CI corals to see if lineage difference holds
hprot_phys_2_lin_CI = hprot_phys_2_lin %>%
  subset(sitename == "CI")

hprot_means_lineage_CI <- summarySE(hprot_phys_2_lin_CI, measurevar="prot_mgcm2", groupvars=c("treat","lineage"))

prot_plot_lineage_CI <- ggplot(hprot_means_lineage_CI,aes(x = treat, y = prot_mgcm2, color = lineage, pch = lineage))+
  theme_bw()+
  geom_point(size = 3, position = position_dodge(width=0.3))+
  geom_errorbar(aes(x = treat, ymax = prot_mgcm2+se, ymin = prot_mgcm2-se), width = .2, position = position_dodge(width=0.3)) +
  scale_color_manual(name = "Lineage",
                     breaks = c("L1","L2"),
                     values = cols_lineage)+
  scale_shape_manual(name = "Lineage",
                     breaks = c("L1","L2"),
                     values=c(19,17))+
  xlab("Treatment")+
  ylab(bquote("Total Protein (mg" ~cm^-2~')'))+
  #ylim(0,0.5) +
  geom_vline(xintercept = 1.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))
prot_plot_lineage_CI

ggsave(prot_plot_lineage_CI, filename = "/Users/hannahaichelman/Documents/BU/TVE/Protein/plots/protein_lineage_CIonly.pdf", width=5, height=4, units=c("in"), useDingbats=FALSE)

##### Host Carbohydrate Concentrations #####
# first, subset physiology data add in standard curve data needed to correct for separate carbohydrate plates
hcarb_phys = phys %>%
  select(frag,treat,hcarbplate_ha,hcarb1_ha,hcarb2_ha,hcarb3_ha,
         survivedtoend,blastvol,gen_site,origsitecode,sitename,fragid,reef,genet,SAcm2,dominant_type,lineage)

head(hcarb_phys)

# used this to sort out who we are still missing carbs from
missing_hcarbs = hcarb_phys %>%
  dplyr::filter(survivedtoend=="yes") %>%
  dplyr::filter(is.na(hcarb1_ha))
# O2C8, I3B4, I4E7

# For both host carbs and symbiont carbs, for all time points,
# we are missing 3, but both of these do not have blaster or blastvol info so not much we can do

# Add in standard curve info for each plate
hcarb_phys$hcarbxcoef <- as.numeric(ifelse(hcarb_phys$hcarbplate_ha == 'O1', 0.1678088,
                                           ifelse(hcarb_phys$hcarbplate_ha == 'O2', 0.187208,
                                                  ifelse(hcarb_phys$hcarbplate_ha == 'O3', 0.229014,
                                                         ifelse(hcarb_phys$hcarbplate_ha == 'O4', 0.231081,
                                                                ifelse(hcarb_phys$hcarbplate_ha == 'O5', 0.196008,
                                                                       ifelse(hcarb_phys$hcarbplate_ha == 'O6', 0.162171,
                                                                              ifelse(hcarb_phys$hcarbplate_ha == 'O7', 0.1765377,
                                                                                     ifelse(hcarb_phys$hcarbplate_ha == 'O8', 0.141797,
                                                                                            ifelse(hcarb_phys$hcarbplate_ha == 'O9', 0.152009,
                                                                                                   ifelse(hcarb_phys$hcarbplate_ha == 'O10', 0.1547928,
                                                                                                          ifelse(hcarb_phys$hcarbplate_ha == 'O11', 0.170403,
                                                                                                                 ifelse(hcarb_phys$hcarbplate_ha == 'O12', 0.142706,
                                                                                                                        ifelse(hcarb_phys$hcarbplate_ha == 'O13', 0.148280,
                                                                                                                               ifelse(hcarb_phys$hcarbplate_ha == 'O14', 0.164326,
                                                                                                                                      ifelse(hcarb_phys$hcarbplate_ha == 'O15', 0.146469,
                                                                                                                                             ifelse(hcarb_phys$hcarbplate_ha == 'O16', 0.175768,
                                                                                                                                                    ifelse(hcarb_phys$hcarbplate_ha == 'O17', 0.163126,
                                                                                                                                                           ifelse(hcarb_phys$hcarbplate_ha == 'O35', 0.1536346,
                                                                                                                                                                  ifelse(hcarb_phys$hcarbplate_ha == 'O36', 0.1717878,
                                                                                                                                                                         'NA'))))))))))))))))))))

hcarb_phys$hcarbint <- as.numeric(ifelse(hcarb_phys$hcarbplate_ha == 'O1', 0.0005952,
                                         ifelse(hcarb_phys$hcarbplate_ha == 'O2', -0.011863,
                                                ifelse(hcarb_phys$hcarbplate_ha == 'O3', -0.005677,
                                                       ifelse(hcarb_phys$hcarbplate_ha == 'O4', -0.001425,
                                                              ifelse(hcarb_phys$hcarbplate_ha == 'O5', -0.001076,
                                                                     ifelse(hcarb_phys$hcarbplate_ha == 'O6', -0.007433,
                                                                            ifelse(hcarb_phys$hcarbplate_ha == 'O7', 0.0005785,
                                                                                   ifelse(hcarb_phys$hcarbplate_ha == 'O8', -0.003569,
                                                                                          ifelse(hcarb_phys$hcarbplate_ha == 'O9', -0.001455,
                                                                                                 ifelse(hcarb_phys$hcarbplate_ha == 'O10', -0.0003914,
                                                                                                        ifelse(hcarb_phys$hcarbplate_ha == 'O11', -0.001060,
                                                                                                               ifelse(hcarb_phys$hcarbplate_ha == 'O12', 0.002935,
                                                                                                                      ifelse(hcarb_phys$hcarbplate_ha == 'O13', 0.002679,
                                                                                                                             ifelse(hcarb_phys$hcarbplate_ha == 'O14', -0.004464,
                                                                                                                                    ifelse(hcarb_phys$hcarbplate_ha == 'O15', 0.008575,
                                                                                                                                           ifelse(hcarb_phys$hcarbplate_ha == 'O16', 0.002615,
                                                                                                                                                  ifelse(hcarb_phys$hcarbplate_ha == 'O17', -0.006809,
                                                                                                                                                         ifelse(hcarb_phys$hcarbplate_ha == 'O35', -0.0004174,
                                                                                                                                                                ifelse(hcarb_phys$hcarbplate_ha == 'O36', 0.0003559,
                                                                                                                                                                       'NA'))))))))))))))))))))

head(hcarb_phys)

# get the carbohydrate values in mg/mL (concentration of D-glucose standard)
# the standard curve equation for carbs is linear:  y = m*x + b
# where y is the calculated concentration in mg/mL, m is the slope, or x coefficient, and b is the intercept
hcarb_phys$calchcarb1 <- ((hcarb_phys$hcarbxcoef*hcarb_phys$hcarb1_ha) + hcarb_phys$hcarbint)
hcarb_phys$calchcarb2 <- ((hcarb_phys$hcarbxcoef*hcarb_phys$hcarb2_ha) + hcarb_phys$hcarbint)
hcarb_phys$calchcarb3 <- ((hcarb_phys$hcarbxcoef*hcarb_phys$hcarb3_ha) + hcarb_phys$hcarbint)

# calculate average carbohydrate concentration, mg/mL
# this code allows us to take into account when there was an outlier replaced with an NA
hcarb_phys = hcarb_phys %>%
  mutate(avghcarb=rowMeans(.[ , c("calchcarb1","calchcarb2","calchcarb3")], na.rm=TRUE))

# [carb concentration in mg/cm2] = ([mean concentration mg/mL] * [dilution factor] * number of mls in slurry / [surface area of fragment]
# dilution factor is 5 here because in the protocol, each reaction has 10 uL of sample and 40 uL of seawater
hcarb_phys$hcarb_mgcm2 <- (hcarb_phys$avghcarb*5*hcarb_phys$blastvol)/hcarb_phys$SAcm2

# make data subsets for stats and plotting
hcarb_phys_all_lin = hcarb_phys %>%
  drop_na(hcarb_mgcm2) %>%
  dplyr::filter(treat!="Control 2") %>%
  dplyr::filter(gen_site!="I4G") %>% # remove clone
  select(frag, gen_site, treat, sitename, dominant_type, lineage, hcarb_mgcm2)

hcarb_phys_2_lin = hcarb_phys_all_lin %>%
  dplyr::filter(is.na(lineage) | lineage!="L3") %>%
  drop_na(sitename, lineage) # doing this for plotting

hcarb_phys_2_lin_symtype = hcarb_phys_2_lin %>%
  drop_na(treat, dominant_type)  # we only have symbiont types from the end of the experiment, so do separate filtering for those plots

initial_hcarb_phys_all_lin = hcarb_phys_all_lin %>%
  dplyr::filter(treat=="Initial")

initial_hcarb_phys_2_lin = initial_hcarb_phys_all_lin %>%
  dplyr::filter(is.na(lineage) | lineage!="L3") # want to keep NA values for lineage here since they still have other info, will remove na's for lineage specific plots

post_hcarb_phys_all_lin = hcarb_phys_all_lin %>%
  dplyr::filter(treat!="Initial")

post_hcarb_phys_2_lin = post_hcarb_phys_all_lin %>%
  dplyr::filter(is.na(lineage) | lineage!="L3")


# STATS
## NEED TO RUN STATS SEPARATELY FOR INITIAL AND FINAL TIME POINTS

model.hcarb.initial = lm(hcarb_mgcm2 ~ lineage+sitename, data = initial_hcarb_phys_2_lin)
summary(model.hcarb.initial)
anova(model.hcarb.initial)
#           Df   Sum Sq  Mean Sq F value  Pr(>F)
# lineage    1 0.052428 0.052428  6.9917 0.01206 *
# sitename   5 0.031164 0.006233  0.8312 0.53612
# Residuals 36 0.269951 0.007499

check_model(model.hcarb.initial)


model.hcarb.final <- lmer(hcarb_mgcm2 ~ lineage+treat*sitename*dominant_type + (1|gen_site), data = post_hcarb_phys_2_lin)
summary(model.hcarb.final)
anova(model.hcarb.final)
# Type III Analysis of Variance Table with Satterthwaite's method
#                               Sum Sq  Mean Sq NumDF  DenDF F value   Pr(>F)
# lineage                      0.08113 0.081127     1 37.872  4.5590 0.039275 *
# treat                        0.19583 0.065277     3 82.374  3.6683 0.015519 *
# sitename                     0.32543 0.065085     5 51.873  3.6575 0.006559 **
# dominant_type                0.06725 0.022418     3 66.851  1.2598 0.295288
# treat:sitename               0.39150 0.026100    15 82.668  1.4667 0.137272
# treat:dominant_type          0.10599 0.015141     7 84.463  0.8509 0.548717
# sitename:dominant_type       0.27350 0.068375     4 72.972  3.8424 0.006878 **
# treat:sitename:dominant_type 0.26459 0.037798     7 89.322  2.1241 0.048912 *

check_model(model.hcarb.final)


#specify model (because we are interested in pairwise, have to include the interaction)
m.emm<- lmer(hcarb_mgcm2 ~ treat*lineage + (1|gen_site), data = post_hcarb_phys_2_lin, REML=FALSE)

emms<-emmeans(m.emm, ~lineage|treat) #, adjust="Bonferoni"
pairs(emms, interaction = "pairwise") %>% rbind(adjust="fdr")
# treat    lineage_pairwise estimate     SE  df t.ratio p.value
# Control  L1 - L2            0.2079 0.0487 154   4.271  0.0001
# Low Var  L1 - L2            0.0433 0.0509 154   0.851  0.3962
# Mod Var  L1 - L2            0.0751 0.0455 154   1.652  0.1341
# High Var L1 - L2            0.0901 0.0479 154   1.879  0.1244



## PLOTS
#SummarySE to format data for plotting
hcarb_means_site_2_lin <- summarySE(hcarb_phys_2_lin, measurevar="hcarb_mgcm2", groupvars=c("treat","sitename"))

# plot, treatment x axis colored by site data figure
hcarb_plot_site <- ggplot(hcarb_means_site_2_lin,aes(x = treat, y = hcarb_mgcm2, fill = sitename))+
  theme_bw()+
  geom_errorbar(aes(x = treat, ymax = hcarb_mgcm2+se, ymin = hcarb_mgcm2-se), width = .2, position = position_dodge(width=0.3)) +
  geom_point(size = 3, pch = 21, position = position_dodge(width=0.3))+
  scale_fill_manual(name = "Site",
                     values = cols_site)+
  xlab("Treatment")+
  ylab(bquote("Total Host Carbohydrate (mg" ~cm^-2~')'))+
  ylim(0,0.8) +
  geom_vline(xintercept = 1.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        legend.position = "none"
)
hcarb_plot_site

ggsave(hcarb_plot_site, filename = "/Users/hannahaichelman/Documents/BU/TVE/Carbohydrates/Plots/hcarb_site_2_lin.pdf", width=5, height=4, units=c("in"), useDingbats=FALSE)

#SummarySE to format data for plotting lineage
hcarb_means_lineage_2_lin <- summarySE(hcarb_phys_2_lin, measurevar="hcarb_mgcm2", groupvars=c("treat","lineage"))

# plot, treatment x axis colored by lineage data figure
hcarb_plot_lineage <- ggplot(hcarb_means_lineage_2_lin,aes(x = treat, y = hcarb_mgcm2, fill = lineage))+
  theme_bw()+
  geom_errorbar(aes(x = treat, ymax = hcarb_mgcm2+se, ymin = hcarb_mgcm2-se), width = .2, position = position_dodge(width=0.3)) +
  geom_point(size = 3, pch = 21, position = position_dodge(width=0.3))+
  scale_fill_manual(name = "Lineage",
                     values = cols_lineage)+
  xlab("Treatment")+
  ylab(bquote("Total Host Carbohydrate (mg" ~cm^-2~')'))+
  ylim(0,0.8) +
  geom_vline(xintercept = 1.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        legend.position = "none"
)
hcarb_plot_lineage

ggsave(hcarb_plot_lineage, filename = "/Users/hannahaichelman/Documents/BU/TVE/Carbohydrates/Plots/hcarb_lineage_2lin.pdf", width=5, height=4, units=c("in"), useDingbats=FALSE)


# plot, treatment x axis and symtype color
#SummarySE to format data for plotting
hcarb_means_sym <- summarySE(hcarb_phys_2_lin_symtype, measurevar="hcarb_mgcm2", groupvars=c("treat","dominant_type"))

# plot, treatment x axis colored by site data figure
hcarb_plot_sym <- ggplot(hcarb_means_sym,aes(x = treat, y = hcarb_mgcm2, fill = dominant_type))+
  theme_bw()+
  geom_errorbar(aes(x = treat, ymax = hcarb_mgcm2+se, ymin = hcarb_mgcm2-se), width = .2, position = position_dodge(width=0.3)) +
  geom_point(size = 3, pch = 21, position = position_dodge(width=0.3))+
  scale_fill_manual(name = "Dominant Symbiont",
                    values = its2_cols_greens)+
  xlab("Treatment")+
  ylab(bquote("Total Host Carbohydrate (mg" ~cm^-2~')'))+
  ylim(0,0.8) +
  #geom_vline(xintercept = 1.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        legend.position = "none"
)
hcarb_plot_sym

ggsave(hcarb_plot_sym, filename = "/Users/hannahaichelman/Documents/BU/TVE/Carbohydrates/Plots/protein_lineage_sym.pdf", width=6, height=4, units=c("in"), useDingbats=FALSE)



# Combine these host carbohydrate plots
hcarb_plots_all = ggarrange(hcarb_plot_site, hcarb_plot_lineage, hcarb_plot_sym,
                           ncol = 3, nrow = 1)
ggsave(hcarb_plots_all, filename = "/Users/hannahaichelman/Documents/BU/TVE/Carbohydrates/Plots/hcarb_all_plots.pdf", width=12, height=4, units=c("in"), useDingbats=FALSE)



# plot just CI corals to see if lineage difference holds
hcarb_phys_2_lin_nona_CI = hcarb_phys_2_lin_nona %>%
  subset(sitename == "CI")

hcarb_means_lineage_CI <- summarySE(hcarb_phys_2_lin_nona_CI, measurevar="hcarb_mgcm2", groupvars=c("treat","lineage"))

hcarb_plot_lineage_CI <- ggplot(hcarb_means_lineage_CI,aes(x = treat, y = hcarb_mgcm2, color = lineage, pch = lineage))+
  theme_bw()+
  geom_point(size = 3, position = position_dodge(width=0.3))+
  geom_errorbar(aes(x = treat, ymax = hcarb_mgcm2+se, ymin = hcarb_mgcm2-se), width = .2, position = position_dodge(width=0.3)) +
  scale_color_manual(name = "Lineage",
                     breaks = c("L1","L2"),
                     values = cols_lineage)+
  scale_shape_manual(name = "Lineage",
                     breaks = c("L1","L2"),
                     values=c(19,17))+
  xlab("Treatment")+
  ylab(bquote("Total Host Carbohydrate (mg" ~cm^-2~')'))+
  #ylim(0,0.4) +
  geom_vline(xintercept = 1.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))
hcarb_plot_lineage_CI

ggsave(hcarb_plot_lineage_CI, filename = "/Users/hannahaichelman/Documents/BU/TVE/Carbohydrates/Plots/hcarb_lineage_CIonly.pdf", width=5, height=4, units=c("in"), useDingbats=FALSE)


##### Sym Carb Concentrations #####
scarb_phys = phys %>%
  select(frag,treat,scarbplate_ha,scarb1_ha,scarb2_ha,scarb3_ha,
         survivedtoend,blastvol,gen_site,origsitecode,sitename,fragid,reef,genet,SAcm2,dominant_type,lineage)

head(scarb_phys)

scarb_phys$scarbxcoef <- as.numeric(ifelse(scarb_phys$scarbplate_ha == 'O18', 0.149980,
                                           ifelse(scarb_phys$scarbplate_ha == 'O19', 0.168115,
                                                  ifelse(scarb_phys$scarbplate_ha == 'O20', 0.1412873,
                                                         ifelse(scarb_phys$scarbplate_ha == 'O21', 0.142275,
                                                                ifelse(scarb_phys$scarbplate_ha == 'O22', 0.123317,
                                                                       ifelse(scarb_phys$scarbplate_ha == 'O23', 0.1539610,
                                                                              ifelse(scarb_phys$scarbplate_ha == 'O24', 0.1672786,
                                                                                     ifelse(scarb_phys$scarbplate_ha == 'O25', 0.2442510,
                                                                                            ifelse(scarb_phys$scarbplate_ha == 'O26', 0.123947,
                                                                                                   ifelse(scarb_phys$scarbplate_ha == 'O27', 0.133494,
                                                                                                          ifelse(scarb_phys$scarbplate_ha == 'O28', 0.136688,
                                                                                                                 ifelse(scarb_phys$scarbplate_ha == 'O29', 0.166897,
                                                                                                                        ifelse(scarb_phys$scarbplate_ha == 'O30', 0.1466,
                                                                                                                               ifelse(scarb_phys$scarbplate_ha == 'O31', 0.209305,
                                                                                                                                      ifelse(scarb_phys$scarbplate_ha == 'O32', 0.1308968,
                                                                                                                                             ifelse(scarb_phys$scarbplate_ha == 'O33', 0.153948,
                                                                                                                                                    ifelse(scarb_phys$scarbplate_ha == 'O34', 0.1407179,
                                                                                                                                                           ifelse(scarb_phys$scarbplate_ha == 'O35', 0.1536346,
                                                                                                                                                                  ifelse(scarb_phys$scarbplate_ha == 'O36', 0.1717878,
                                                                                                                                                                         'NA'))))))))))))))))))))

scarb_phys$scarbint <- as.numeric(ifelse(scarb_phys$scarbplate_ha == 'O18', -0.004060,
                                         ifelse(scarb_phys$scarbplate_ha == 'O19', 0.001898,
                                                ifelse(scarb_phys$scarbplate_ha == 'O20', -0.0003309,
                                                       ifelse(scarb_phys$scarbplate_ha == 'O21', 0.001065,
                                                              ifelse(scarb_phys$scarbplate_ha == 'O22', -0.001284,
                                                                     ifelse(scarb_phys$scarbplate_ha == 'O23', -0.0005789,
                                                                            ifelse(scarb_phys$scarbplate_ha == 'O24', 0.0012109,
                                                                                   ifelse(scarb_phys$scarbplate_ha == 'O25', -0.0003185,
                                                                                          ifelse(scarb_phys$scarbplate_ha == 'O26', 0.002024,
                                                                                                 ifelse(scarb_phys$scarbplate_ha == 'O27', -0.000288,
                                                                                                        ifelse(scarb_phys$scarbplate_ha == 'O28', -0.001351,
                                                                                                               ifelse(scarb_phys$scarbplate_ha == 'O29', -0.004318,
                                                                                                                      ifelse(scarb_phys$scarbplate_ha == 'O30', -0.00004093,
                                                                                                                             ifelse(scarb_phys$scarbplate_ha == 'O31', -0.001261,
                                                                                                                                    ifelse(scarb_phys$scarbplate_ha == 'O32', -0.0001588,
                                                                                                                                           ifelse(scarb_phys$scarbplate_ha == 'O33', -0.005824,
                                                                                                                                                  ifelse(scarb_phys$scarbplate_ha == 'O34', 0.000427,
                                                                                                                                                         ifelse(scarb_phys$scarbplate_ha == 'O35', -0.0004174,
                                                                                                                                                                ifelse(scarb_phys$scarbplate_ha == 'O36', 0.0003559,
                                                                                                                                                                       'NA'))))))))))))))))))))

# get the carbohydrate values in mg/mL (concentration of D-glucose standard)
# the standard curve equation for carbs is linear - y = m*x + b -
# where y is the calculated concentration in mg/mL, m is the slope, or x coefficient, and b is the intercept
scarb_phys$calcscarb1 <- ((scarb_phys$scarbxcoef*scarb_phys$scarb1_ha) + scarb_phys$scarbint)
scarb_phys$calcscarb2 <- ((scarb_phys$scarbxcoef*scarb_phys$scarb2_ha) + scarb_phys$scarbint)
scarb_phys$calcscarb3 <- ((scarb_phys$scarbxcoef*scarb_phys$scarb3_ha) + scarb_phys$scarbint)

# calculate average protein concentration, mg/mL
scarb_phys = scarb_phys %>%
  mutate(avgscarb=rowMeans(.[ , c("calcscarb1","calcscarb2","calcscarb3")], na.rm=TRUE))

# [carb concentration in in mg/cm2] = ([mean concentration mg/mL] * [dilution factor] * number of mls in slurry / [surface area of fragment]
# dilution factor is 5 here because in the protocol, each reaction has 10 uL of sample and 40 uL of seawater
scarb_phys$scarb_mgcm2 <- (scarb_phys$avgscarb*5*scarb_phys$blastvol)/scarb_phys$SAcm2

head(scarb_phys)

# make data subsets for stats and plotting
scarb_phys_all_lin = scarb_phys %>%
  drop_na(scarb_mgcm2) %>%
  dplyr::filter(treat!="Control 2") %>%
  dplyr::filter(gen_site!="I4G") %>% # remove clone
  select(frag, gen_site, treat, sitename, dominant_type, lineage, scarb_mgcm2)

scarb_phys_2_lin = scarb_phys_all_lin %>%
  dplyr::filter(is.na(lineage) | lineage!="L3") %>%
  drop_na(sitename, lineage) # doing this for plotting

scarb_phys_2_lin_symtype = scarb_phys_2_lin %>%
  drop_na(treat, dominant_type)  # we only have symbiont types from the end of the experiment, so do separate filtering for those plots

initial_scarb_phys_all_lin = scarb_phys_all_lin %>%
  dplyr::filter(treat=="Initial")

initial_scarb_phys_2_lin = initial_scarb_phys_all_lin %>%
  dplyr::filter(is.na(lineage) | lineage!="L3") # want to keep NA values for lineage here since they still have other info, will remove na's for lineage specific plots

post_scarb_phys_all_lin = scarb_phys_all_lin %>%
  dplyr::filter(treat!="Initial")

post_scarb_phys_2_lin = post_scarb_phys_all_lin %>%
  dplyr::filter(is.na(lineage) | lineage!="L3")


# STATS
## NEED TO RUN STATS SEPARATELY FOR INITIAL AND FINAL TIME POINTS
# use lmer()
model.scarb.initial = lm(scarb_mgcm2 ~ lineage+sitename, data = initial_scarb_phys_2_lin)
summary(model.scarb.initial)
anova(model.scarb.initial)
#           Df   Sum Sq   Mean Sq F value   Pr(>F)
# lineage    1 0.021801 0.0218014  9.1326 0.004604 **
# sitename   5 0.013268 0.0026536  1.1116 0.371517
# Residuals 36 0.085940 0.0023872

check_model(model.scarb.initial)


model.scarb.final <- lmer(scarb_mgcm2 ~ lineage+treat+sitename+dominant_type + (1|gen_site), data = post_scarb_phys_2_lin)
summary(model.scarb.final)
anova(model.scarb.final)
# Type III Analysis of Variance Table with Satterthwaite's method
#                 Sum Sq  Mean Sq NumDF DenDF F value   Pr(>F)
# lineage       0.034092 0.034092     1   131  7.4980 0.007036 **
# treat         0.052405 0.017468     3   131  3.8419 0.011258 *
# sitename      0.054374 0.010875     5   131  2.3918 0.041123 *
# dominant_type 0.012618 0.004206     3   131  0.9250 0.430712

check_model(model.scarb.final)


#specify model (because we are interested in pairwise, have to include the interaction)
m.emm<- lmer(scarb_mgcm2 ~ treat*lineage + (1|gen_site), data = post_scarb_phys_2_lin, REML=FALSE)

emms<-emmeans(m.emm, ~lineage|treat) #, adjust="Bonferoni"
pairs(emms, interaction = "pairwise") %>% rbind(adjust="fdr")
# treat    lineage_pairwise estimate     SE  df t.ratio p.value
# Control  L1 - L2           0.04923 0.0237 154   2.080  0.0784
# Low Var  L1 - L2           0.04359 0.0248 154   1.760  0.1071
# Mod Var  L1 - L2           0.00979 0.0221 154   0.442  0.6588
# High Var L1 - L2           0.07896 0.0233 154   3.386  0.0036


# PLOTS
#SummarySE to format data for plotting
scarb_means_site_2_lin <- summarySE(scarb_phys_2_lin, measurevar="scarb_mgcm2", groupvars=c("treat","sitename"))

# plot, treatment x axis colored by site data figure
scarb_plot_site <- ggplot(scarb_means_site_2_lin,aes(x = treat, y = scarb_mgcm2, fill = sitename))+
  theme_bw()+
  geom_errorbar(aes(x = treat, ymax = scarb_mgcm2+se, ymin = scarb_mgcm2-se), width = .2, position = position_dodge(width=0.3)) +
  geom_point(size = 3, pch = 21, position = position_dodge(width=0.3))+
  scale_fill_manual(name = "Site",
                    values = cols_site)+
  xlab("Treatment")+
  ylab(bquote("Total Symbiont Carbohydrate (mg" ~cm^-2~')'))+
  ylim(0,0.35) +
  geom_vline(xintercept = 1.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        legend.position = "none"
)
scarb_plot_site

ggsave(scarb_plot_site, filename = "/Users/hannahaichelman/Documents/BU/TVE/Carbohydrates/Plots/scarb_site_2_lin.pdf", width=5, height=4, units=c("in"), useDingbats=FALSE)

#SummarySE to format data for plotting lineage
scarb_means_lineage_2_lin <- summarySE(scarb_phys_2_lin, measurevar="scarb_mgcm2", groupvars=c("treat","lineage"))

# plot, treatment x axis colored by lineage data figure
scarb_plot_lineage <- ggplot(scarb_means_lineage_2_lin,aes(x = treat, y = scarb_mgcm2, fill = lineage))+
  theme_bw()+
  geom_errorbar(aes(x = treat, ymax = scarb_mgcm2+se, ymin = scarb_mgcm2-se), width = .2, position = position_dodge(width=0.3)) +
  geom_point(size = 3, pch = 21, position = position_dodge(width=0.3))+
  scale_fill_manual(name = "Lineage",
                     values = cols_lineage)+
  xlab("Treatment")+
  ylab(bquote("Total Symbiont Carbohydrate (mg" ~cm^-2~')'))+
  ylim(0,0.35) +
  geom_vline(xintercept = 1.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        legend.position = "none"
)

scarb_plot_lineage

ggsave(scarb_plot_lineage, filename = "/Users/hannahaichelman/Documents/BU/TVE/Carbohydrates/Plots/scarb_lineage_2lin.pdf", width=5, height=4, units=c("in"), useDingbats=FALSE)


# plot, treatment x axis and symtype color
#SummarySE to format data for plotting
scarb_means_sym <- summarySE(scarb_phys_2_lin_symtype, measurevar="scarb_mgcm2", groupvars=c("treat","dominant_type"))

# plot, treatment x axis colored by site data figure
scarb_plot_sym <- ggplot(scarb_means_sym,aes(x = treat, y = scarb_mgcm2, fill = dominant_type))+
  theme_bw()+
  geom_errorbar(aes(x = treat, ymax = scarb_mgcm2+se, ymin = scarb_mgcm2-se), width = .2, position = position_dodge(width=0.3)) +
  geom_point(size = 3, pch = 21, position = position_dodge(width=0.3))+
  scale_fill_manual(name = "Dominant Symbiont",
                    values = its2_cols_greens)+
  xlab("Treatment")+
  ylab(bquote("Total Symbiont Carbohydrate (mg" ~cm^-2~')'))+
  ylim(0,0.35) +
  #geom_vline(xintercept = 1.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        legend.position = "none"
)
scarb_plot_sym

ggsave(scarb_plot_sym, filename = "/Users/hannahaichelman/Documents/BU/TVE/Carbohydrates/Plots/scarb_lineage_sym.pdf", width=6, height=4, units=c("in"), useDingbats=FALSE)



# Combine these symbiont carbohydrate plots
scarb_plots_all = ggarrange(scarb_plot_site, scarb_plot_lineage, scarb_plot_sym,
                            ncol = 3, nrow = 1)
ggsave(scarb_plots_all, filename = "/Users/hannahaichelman/Documents/BU/TVE/Carbohydrates/Plots/scarb_all_plots.pdf", width=12, height=4, units=c("in"), useDingbats=FALSE)


##### Symbiont Density #####
sym_phys = phys %>%
  select(frag,treat,symcount1,symcount2,symcount3,
         survivedtoend,blastvol,gen_site,origsitecode,sitename,fragid,reef,genet,SAcm2,dominant_type,lineage)

head(sym_phys)

# Calculate average per surface area
# make data subsets for stats and plotting
sym_phys_all_lin = sym_phys %>%
  mutate(avgsymcount=rowMeans(.[ , c("symcount1","symcount2","symcount3")], na.rm=TRUE)) %>%
  mutate(sym_cellsmL=avgsymcount*10000) %>% # this is based on the following website http://insilico.ehu.eus/counting_chamber/neubauer_improved.php
  mutate(sym_cm2=(sym_cellsmL*blastvol)/SAcm2) %>%
  mutate(sym_cm2_div = sym_cm2/1000000) %>%
  dplyr::filter(treat!="Control 2") %>%
  dplyr::filter(frag!="I3E6_v2") %>% # same data for sym count for original and v2 of I3E6
  dplyr::filter(frag!="I3A4") %>%
  dplyr::filter(gen_site != "I4G") %>% # remove clone
  drop_na(sym_cm2) %>%
  select(frag, gen_site, treat, sitename, dominant_type, lineage, sym_cm2, sym_cm2_div)


sym_phys_2_lin = sym_phys_all_lin %>%
  dplyr::filter(is.na(lineage) | lineage!="L3") %>%
  drop_na(sitename, lineage) # doing this for plotting

sym_phys_2_lin_symtype = sym_phys_2_lin %>%
  drop_na(treat, dominant_type)  # we only have symbiont types from the end of the experiment, so do separate filtering for those plots

initial_sym_phys_all_lin = sym_phys_all_lin %>%
  dplyr::filter(treat=="Initial")

initial_sym_phys_2_lin = initial_sym_phys_all_lin %>%
  dplyr::filter(is.na(lineage) | lineage!="L3") # want to keep NA values for lineage here since they still have other info, will remove na's for lineage specific plots

post_sym_phys_all_lin = sym_phys_all_lin %>%
  dplyr::filter(treat!="Initial")

post_sym_phys_2_lin = post_sym_phys_all_lin %>%
  dplyr::filter(is.na(lineage) | lineage!="L3")


# STATS
## NEED TO RUN STATS SEPARATELY FOR INITIAL AND FINAL TIME POINTS
# use lmer()
model.sym.initial = lm(sym_cm2 ~ lineage+sitename, data = initial_sym_phys_2_lin)
summary(model.sym.initial)
anova(model.sym.initial)
#           Df     Sum Sq    Mean Sq F value    Pr(>F)
# lineage    1 9.3732e+12 9.3732e+12 19.0628 0.0001022 ***
# sitename   5 3.2179e+12 6.4358e+11  1.3089 0.2820781
# Residuals 36 1.7701e+13 4.9170e+11

check_model(model.sym.initial)


model.sym.final <- lmer(sym_cm2 ~ lineage+treat+sitename+dominant_type + (1|gen_site), data = post_sym_phys_2_lin)
summary(model.sym.final)
anova(model.sym.final)
# Type III Analysis of Variance Table with Satterthwaite's method
#                   Sum Sq    Mean Sq NumDF  DenDF F value    Pr(>F)
# lineage       8.3103e+11 8.3103e+11     1   35.8  4.6251   0.03833 *
# treat         5.2411e+12 1.7470e+12     3 7942.4  9.7232 2.115e-06 ***
# sitename      1.7616e+12 3.5232e+11     5   36.7  1.9608   0.10776
# dominant_type 2.2922e+11 7.6407e+10     3   87.8  0.4252   0.73538

check_model(model.sym.final)


#specify model (because we are interested in pairwise, have to include the interaction)
m.emm<- lmer(sym_cm2 ~ treat*lineage + (1|gen_site), data = post_sym_phys_2_lin, REML=FALSE)

emms<-emmeans(m.emm, ~lineage|treat) #, adjust="Bonferoni"
pairs(emms, interaction = "pairwise") %>% rbind(adjust="fdr")
# treat    lineage_pairwise estimate     SE  df t.ratio p.value
# Control  L1 - L2            492769 152557 151   3.230  0.0061
# Low Var  L1 - L2            188787 159361 152   1.185  0.2380
# Mod Var  L1 - L2            429179 142727 150   3.007  0.0062
# High Var L1 - L2            243636 152229 151   1.600  0.1488


# PLOTS
#SummarySE to format data for plotting
sym_means_site_2_lin <- summarySE(sym_phys_2_lin, measurevar="sym_cm2_div", groupvars=c("treat","sitename"))

# plot, treatment x axis colored by site data figure
sym_plot_site <- ggplot(sym_means_site_2_lin,aes(x = treat, y = sym_cm2_div, fill = sitename))+
  theme_bw()+
  geom_errorbar(aes(x = treat, ymax = sym_cm2_div+se, ymin = sym_cm2_div-se), width = .2, position = position_dodge(width=0.3)) +
  geom_point(size = 3, pch = 21, position = position_dodge(width=0.3))+
  scale_fill_manual(name = "Site",
                    values = cols_site)+
  xlab("Treatment")+
  ylab(bquote("Symbiont Density ("~x10^6~ 'cells' ~cm^-2~')'))+
  ylim(0,2.55) +
  geom_vline(xintercept = 1.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        legend.position = "none"
)
sym_plot_site

ggsave(sym_plot_site, filename = "/Users/hannahaichelman/Documents/BU/TVE/SymbiontDensity/plots/sym_site_2_lin.pdf", width=5, height=4, units=c("in"), useDingbats=FALSE)

#SummarySE to format data for plotting
sym_means_lineage_2_lin <- summarySE(sym_phys_2_lin, measurevar="sym_cm2_div", groupvars=c("treat","lineage"))

# plot, treatment x axis colored by lineage data figure
sym_plot_lineage <- ggplot(sym_means_lineage_2_lin,aes(x = treat, y = sym_cm2_div, fill = lineage))+
  theme_bw()+
  geom_errorbar(aes(x = treat, ymax = sym_cm2_div+se, ymin = sym_cm2_div-se), width = .2, position = position_dodge(width=0.3)) +
  geom_point(size = 3, pch = 21, position = position_dodge(width=0.3))+
  scale_fill_manual(name = "Lineage",
                    values = cols_lineage)+
  xlab("Treatment")+
  ylab(bquote("Symbiont Density ("~x10^6~ 'cells' ~cm^-2~')'))+
  ylim(0,2.55) +
  geom_vline(xintercept = 1.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        legend.position = "none"
  )
sym_plot_lineage

ggsave(sym_plot_lineage, filename = "/Users/hannahaichelman/Documents/BU/TVE/SymbiontDensity/plots/sym_lineage_2_lin.pdf", width=5, height=4, units=c("in"), useDingbats=FALSE)


# plot, treatment x axis and symtype color
#SummarySE to format data for plotting
sym_means_sym <- summarySE(sym_phys_2_lin_symtype, measurevar="sym_cm2_div", groupvars=c("treat","dominant_type"))

# plot, treatment x axis colored by site data figure
sym_plot_sym <- ggplot(sym_means_sym,aes(x = treat, y = sym_cm2_div, fill = dominant_type))+
  theme_bw()+
  geom_errorbar(aes(x = treat, ymax = sym_cm2_div+se, ymin = sym_cm2_div-se), width = .2, position = position_dodge(width=0.3)) +
  geom_point(size = 3, pch = 21, position = position_dodge(width=0.3))+
  scale_fill_manual(name = "Dominant Symbiont",
                    values = its2_cols_greens)+
  xlab("Treatment")+
  ylab(bquote("Symbiont Density ("~x10^6~ 'cells' ~cm^-2~')'))+
  ylim(0,2.55) +
  #geom_vline(xintercept = 1.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        legend.position = "none"
  )
sym_plot_sym

ggsave(sym_plot_sym, filename = "/Users/hannahaichelman/Documents/BU/TVE/SymbiontDensity/plots/sym_lineage_2_lin_symtype.pdf", width=6, height=4, units=c("in"), useDingbats=FALSE)



# Combine these symbiont density plots
sym_plots_all = ggarrange(sym_plot_site, sym_plot_lineage, sym_plot_sym,
                            ncol = 3, nrow = 1)
ggsave(sym_plots_all, filename = "/Users/hannahaichelman/Documents/BU/TVE/SymbiontDensity/plots/sym_all_plots.pdf", width=12, height=4, units=c("in"), useDingbats=FALSE)




## Checking for outliers - used this to remove I3A4
#plot, facet by treatment
plot <- ggplot(sym_phys_nona,aes(x = sitename, y = sym_cm2, color = sitename, pch = sitename, text = paste("Nubbin:",frag)))+
  theme_bw()+
  geom_boxplot()+
  ggtitle("Symbiont Density (cells/cm2)")+
  scale_fill_manual(values = cols_site)+
  geom_jitter(aes(color = sitename),width = 0.3)+
  scale_color_manual(values = cols_site)+
  theme(axis.title.x = element_blank())+
  theme(axis.text.x=element_text(angle = 45, vjust = 0.5))+
  facet_wrap(~treat, nrow=5)
#  labs(y=expression(paste("Protein (mg cm"^-2*' )')))
ggplotly(subplot(list(plot),titleY=F) %>% layout(showlegend=F))

# facet by genotype
plot <- ggplot(phys,aes(x = treat, y = sym_cm2, color = sitename, pch = sitename, text = paste("Nubbin:",frag)))+
  theme_bw()+
  geom_jitter(width = 0.3, aes(fill = sitename))+
  scale_fill_manual(values = palsite)+
  scale_color_manual(values = palsite)+
  theme(axis.title.x = element_blank())+
  theme(axis.text.x=element_text(angle = 45, vjust = 0.5))+
  facet_wrap(~genet+origsitecode, nrow=4)+
  ggtitle("Symbiont Density (cells/cm2)")
#  labs(y=expression(paste("Protein (mg cm"^-2*' )')))
ggplotly(subplot(list(plot),titleY=T) %>% layout(showlegend=T))


##### Chlorophyll Concentration #####

## read in additional chlorophyll data that Olivia added
chl <- read.csv('/Users/hannahaichelman/Documents/BU/TVE/Chlorophyll/OliviaSamples/TVEChlorphyll_ON.csv')

phys_for_chl = phys %>%
  select(frag, survivedtoend, treat, blastvol, blaster, SAcm2, reef, sitename, gen_site, dominant_type, lineage) %>%
  mutate(treat = as.factor(treat), blastvol = as.numeric(blastvol))

#want to keep all of the individuals in the original dtv spreadsheet and add NA's for any that don't have chl info for
chl_phys <- left_join(phys_for_chl, chl, by = "frag")

# find mean absorbance values for each wavelength
chl_phys = chl_phys %>%
  mutate(chl630_avg=rowMeans(.[ , c("blk_630_1","blk_630_2","blk_630_3")], na.rm=TRUE)) %>%
  mutate(chl663_avg=rowMeans(.[ , c("blk_663_1","blk_663_2","blk_663_3")], na.rm=TRUE))

# calculate chlorophyll A and C2 from values in ug/mL (equations come from Carly Kenkel's method, originally a Jeffrey and Haxo 1968 equation)
# 20 cancels out mL so you end up with ug/cm2
chl_phys$chlA <- ((((13.31*chl_phys$chl663_avg) - (0.27*chl_phys$chl630_avg))*20)/chl_phys$SAcm2)
chl_phys$chlC2 <- ((((-8.37*chl_phys$chl663_avg) + (51.72*chl_phys$chl630_avg))*20)/chl_phys$SAcm2)

head(chl_phys)

# make data subsets for stats and plotting
chl_phys_all_lin = chl_phys %>%
  drop_na(chlA) %>%
  dplyr::filter(treat!="Control 2") %>%
  dplyr::filter(frag!="I2A12", frag!="I2A7", frag!="O4D8", frag!="I3E10", frag!="I3I10", frag!="I3E6_v2") %>% #these frags are being removed because they were duplicate genotypes within treatment and had the most complete information of the two
  dplyr::filter(frag!="I3A4") %>% # unexplained high outlier
  dplyr::filter(gen_site != "I4G") %>% # clone with I4F, remove from dataset
  select(frag, gen_site, treat, sitename, dominant_type, lineage, chlA, chlC2)

chl_phys_2_lin = chl_phys_all_lin %>%
  dplyr::filter(is.na(lineage) | lineage!="L3") %>%
  drop_na(sitename, lineage) # doing this for plotting

chl_phys_2_lin_symtype = chl_phys_2_lin %>%
  drop_na(treat, dominant_type)  # we only have symbiont types from the end of the experiment, so do separate filtering for those plots

initial_chl_phys_all_lin = chl_phys_all_lin %>%
  dplyr::filter(treat=="Initial")

initial_chl_phys_2_lin = initial_chl_phys_all_lin %>%
  dplyr::filter(is.na(lineage) | lineage!="L3") # want to keep NA values for lineage here since they still have other info, will remove na's for lineage specific plots

post_chl_phys_all_lin = chl_phys_all_lin %>%
  dplyr::filter(treat!="Initial")

post_chl_phys_2_lin = post_chl_phys_all_lin %>%
  dplyr::filter(is.na(lineage) | lineage!="L3")


# STATS
## NEED TO RUN STATS SEPARATELY FOR INITIAL AND FINAL TIME POINTS
# use lmer()
model.chl.initial = lm(chlA ~ lineage+sitename, data = initial_chl_phys_2_lin)
summary(model.chl.initial)
anova(model.chl.initial)
#           Df Sum Sq Mean Sq F value    Pr(>F)
# lineage    1 297.51 297.505 23.5187 2.377e-05 ***
# sitename   5 141.08  28.215  2.2305   0.07232 .
# Residuals 36 455.39  12.650

check_model(model.chl.initial)


model.chl.final <- lmer(chlA ~ treat+lineage+sitename+dominant_type + (1|gen_site), data = post_chl_phys_2_lin)
summary(model.chl.final)
anova(model.chl.final)
# Type III Analysis of Variance Table with Satterthwaite's method
#               Sum Sq Mean Sq NumDF   DenDF F value  Pr(>F)
# treat         25.757  8.5855     3 107.320  3.2167 0.02575 *
# lineage        9.504  9.5043     1  33.741  3.5610 0.06778 .
# sitename      40.302  8.0604     5  35.670  3.0200 0.02249 *
# dominant_type 13.384  4.4614     3  97.284  1.6716 0.17816

check_model(model.chl.final)


#specify model (because we are interested in pairwise, have to include the interaction)
m.emm<- lmer(chlA ~ treat*lineage + (1|gen_site), data = post_chl_phys_2_lin, REML=FALSE)

emms<-emmeans(m.emm, ~lineage|treat) #, adjust="Bonferoni"
pairs(emms, interaction = "pairwise") %>% rbind(adjust="fdr")
# treat    lineage_pairwise estimate    SE  df t.ratio p.value
# Control  L1 - L2             1.687 0.658 151   2.561  0.0373
# Low Var  L1 - L2             0.624 0.669 152   0.932  0.3527
# Mod Var  L1 - L2             0.972 0.598 141   1.627  0.1414
# High Var L1 - L2             1.497 0.629 146   2.380  0.0373

# PLOTS
#SummarySE to format data for plotting
chla_means_site_2_lin <- summarySE(chl_phys_2_lin, measurevar="chlA", groupvars=c("treat","sitename"))

# plot, treatment x axis colored by site data figure
chla_plot_site <- ggplot(chla_means_site_2_lin,aes(x = treat, y = chlA, fill = sitename))+
  theme_bw()+
  geom_errorbar(aes(x = treat, ymax = chlA+se, ymin = chlA-se), width = .2, position = position_dodge(width=0.3)) +
  geom_point(size = 3, pch = 21, position = position_dodge(width=0.3))+
  scale_fill_manual(name = "Site",
                    values = cols_site)+
  xlab("Treatment")+
  ylab(bquote("Chl a ("*mu*"g" ~cm^-2~')'))+
  ylim(0,15) +
  geom_vline(xintercept = 1.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        legend.position = "none")
chla_plot_site

ggsave(chla_plot_site, filename = "/Users/hannahaichelman/Documents/BU/TVE/Chlorophyll/plots/chla_site_2_lin.pdf", width=5, height=4, units=c("in"), useDingbats=FALSE)

#SummarySE to format data for plotting
chla_means_lineage_2_lin <- summarySE(chl_phys_2_lin, measurevar="chlA", groupvars=c("treat","lineage"))

# plot, treatment x axis colored by lineage data figure
chla_plot_lineage <- ggplot(chla_means_lineage_2_lin,aes(x = treat, y = chlA, fill = lineage))+
  theme_bw()+
  geom_errorbar(aes(x = treat, ymax = chlA+se, ymin = chlA-se), width = .2, position = position_dodge(width=0.3)) +
  geom_point(size = 3, pch = 21, position = position_dodge(width=0.3))+
  scale_fill_manual(name = "Lineage",
                    values = cols_lineage)+
  xlab("Treatment")+
  ylab(bquote("Chl a ("*mu*"g" ~cm^-2~')'))+
  ylim(0,15) +
  geom_vline(xintercept = 1.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        legend.position = "none")
chla_plot_lineage

ggsave(chla_plot_lineage, filename = "/Users/hannahaichelman/Documents/BU/TVE/Chlorophyll/plots/chla_lineage_2_lin.pdf", width=5, height=4, units=c("in"), useDingbats=FALSE)

# plot, treatment x axis and symtype color
#SummarySE to format data for plotting
chla_means_sym <- summarySE(chl_phys_2_lin_symtype, measurevar="chlA", groupvars=c("treat","dominant_type"))

# plot, treatment x axis colored by site data figure
chla_plot_sym <- ggplot(chla_means_sym,aes(x = treat, y = chlA, fill = dominant_type))+
  theme_bw()+
  geom_errorbar(aes(x = treat, ymax = chlA+se, ymin = chlA-se), width = .2, position = position_dodge(width=0.3)) +
  geom_point(size = 3, pch = 21, position = position_dodge(width=0.3))+
  scale_fill_manual(name = "Dominant Symbiont",
                    values = its2_cols_greens)+
  xlab("Treatment")+
  ylab(bquote("Chl a ("*mu*"g" ~cm^-2~')'))+
  ylim(0,15) +
  #geom_vline(xintercept = 1.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        legend.position = "none")
chla_plot_sym

ggsave(chla_plot_sym, filename = "/Users/hannahaichelman/Documents/BU/TVE/Chlorophyll/plots/chla_symtype_2_lin.pdf", width=6, height=4, units=c("in"), useDingbats=FALSE)



# Combine these chlorophyll A plots
chla_plots_all = ggarrange(chla_plot_site, chla_plot_lineage, chla_plot_sym,
                           ncol = 3, nrow = 1)
ggsave(chla_plots_all, filename = "/Users/hannahaichelman/Documents/BU/TVE/Chlorophyll/plots/chla_all_plots.pdf", width=12, height=4, units=c("in"), useDingbats=FALSE)




# plot just CI corals to see if lineage difference holds
chl_phys_lineage_nona_CI = chl_phys_2_lin_nona %>%
  subset(sitename == "CI")

chla_means_lineage_CI <- summarySE(chl_phys_lineage_nona_CI, measurevar="chlA", groupvars=c("treat","lineage"))

chla_plot_lineage_CI <- ggplot(chla_means_lineage_CI,aes(x = treat, y = chlA, color = lineage, pch = lineage))+
  theme_bw()+
  geom_point(size = 3, position = position_dodge(width=0.3))+
  geom_errorbar(aes(x = treat, ymax = chlA+se, ymin = chlA-se), width = .2, position = position_dodge(width=0.3)) +
  scale_color_manual(name = "Lineage",
                     breaks = c("L1","L2"),
                     values = cols_lineage)+
  scale_shape_manual(name = "Lineage",
                     breaks = c("L1","L2"),
                     values=c(19,17))+
  xlab("Treatment")+
  ylab("Chlorophyll A")+
  #ylim(0,0.4) +
  geom_vline(xintercept = 1.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))
chla_plot_lineage_CI

ggsave(chla_plot_lineage_CI, filename = "/Users/hannahaichelman/Documents/BU/TVE/Chlorophyll/plots/chla_lineage_CIonly.pdf", width=5, height=4, units=c("in"), useDingbats=FALSE)



# look for outliers with ggplotly
# chl A
# violin plot with individual data point overlay
plot <- ggplot(chl_phys,aes(x = sitename, y = chlA, color = sitename, pch = sitename,text = paste("Frag:",frag)))+
  theme_bw()+
  geom_boxplot(alpha = 0.2, aes(fill = sitename))+
  scale_fill_manual(values = cols_site)+
  geom_jitter(aes(color = sitename),width = 0.3)+
  ggtitle("Chl A")+
  scale_color_manual(values = cols_site)+
  theme(axis.title.x = element_blank())+
  theme(axis.text.x=element_text(angle = 45, vjust = 0.5))
#  labs(y=expression(paste("Chl A (ug cm"^-2*' )')))
ggplotly(subplot(list(plot),nrows=1,titleY=T) %>% layout(showlegend=T))

# chl C2
plot <- ggplot(chl_phys,aes(x = sitename, y = chlC2, color = sitename, pch = sitename,text = paste("Frag:",frag)))+
  theme_bw()+
  geom_boxplot(alpha = 0.2, aes(fill = sitename))+
  scale_fill_manual(values = palsite)+
  geom_jitter(aes(color = sitename),width = 0.3)+
  ggtitle("Chl C2")+
  scale_color_manual(values = palsite)+
  theme(axis.title.x = element_blank())+
  theme(axis.text.x=element_text(angle = 45, vjust = 0.5))
#  labs(y=expression(paste("Chl C2 (ug cm"^-2*' )')))
ggplotly(subplot(list(plot),nrows=1,titleY=T) %>% layout(showlegend=T))

# there are lots of negative chlorophyll c2 values...seems like something was off with these readings so omitting from the manuscript.

##### Growth #####
# need to re-read in the data sheet for calcification only
post_phys_forcalc <- read.csv('dtvmaster.csv')

# remove unexplained outliers, and subset by only the data we need here
calc_phys = post_phys_forcalc %>%
  select(frag, survivedtoend, treat, blastvol, blaster,
         t2sastan1, t2sastan2, t2sastan3, t2sarec1, t2sarec2, t2sarec3,
         t0sastan1, t0sastan2, t0sastan3, t0sarec1, t0sarec2, t0sarec3,
         t3sastan1, t3sastan2, t3sastan3, t3sarec1, t3sarec2, t3sarec3,
         pabwrec1,pabwrec2,pabwrec3,t0bwrec1, t0bwrec2,t0bwrec3,t1bwrec1, t1bwrec2,t1bwrec3,t2bwrec1, t2bwrec2,t2bwrec3,t3bwrec1, t3bwrec2,t3bwrec3) %>%
  mutate(treat = as.factor(treat), blastvol = as.numeric(blastvol)) %>%
  mutate_at(c(6:23), as.numeric) %>%
  dplyr::filter(frag!="I2A12", frag!="I2A7", frag!="O4D8", frag!="I3E10", frag!="I3I10") %>% #these frags are being removed because they were duplicate genotypes within treatment and had the most complete information of the two
  dplyr::filter(frag!="O2F1", frag!="O4G7", frag!="O2I3", frag!="O3H3", frag!="I3D8", frag!="O3F4", frag!="O3G10", frag!="O3F3", frag!="I4F8") #these frags are being removed because they were unexplainable outliers with weird weight values

#calculate surface area of the coral fragments
calc_phys$SAcm2_t0 <- (4*((calc_phys$t0sarec1+calc_phys$t0sarec2+calc_phys$t0sarec3)/3))/((calc_phys$t0sastan1+calc_phys$t0sastan2+calc_phys$t0sastan3)/3)
calc_phys$SAcm2_t2 <- (4*((calc_phys$t2sarec1+calc_phys$t2sarec2+calc_phys$t2sarec3)/3))/((calc_phys$t2sastan1+calc_phys$t2sastan2+calc_phys$t2sastan3)/3)
calc_phys$SAcm2_t3 <- (4*((calc_phys$t3sarec1+calc_phys$t3sarec2+calc_phys$t3sarec3)/3))/((calc_phys$t3sastan1+calc_phys$t3sastan2+calc_phys$t3sastan3)/3)

#add in descriptive information for coral samples
calc_phys$origsitecode <- substr(calc_phys$frag, 1, 2)
calc_phys$origsitecode <- as.factor(calc_phys$origsitecode)

# add in descriptive site name with inshore/offshore indicator
calc_phys$sitename <- ifelse(calc_phys$origsitecode == 'I2', 'SP',
                             ifelse(calc_phys$origsitecode == 'I3', 'CI',
                                    ifelse(calc_phys$origsitecode == 'I4', 'PD',
                                           ifelse(calc_phys$origsitecode == 'O2', 'BS',
                                                  ifelse(calc_phys$origsitecode == 'O3', 'CA',
                                                         'BN')))))
calc_phys$sitename <- as.factor(calc_phys$sitename)

# extract the genotype and fragment number
calc_phys$fragid <- substr(calc_phys$frag,3,5)

# add inshore/offshore designation
calc_phys$reef <- substr(calc_phys$frag,1,1)
calc_phys$reef <- ifelse(calc_phys$reef == 'O', 'Outer Reef', 'Inner Reef')
calc_phys$reef <- as.factor(calc_phys$reef)

# add in genotype
calc_phys$genet <- substr(calc_phys$fragid,1,1)

#create a new column of combined genotype that will be used as a random effect in stats later
calc_phys = calc_phys %>%
  unite(gen_site, c(origsitecode,genet), sep = "", remove = FALSE) %>%
  mutate(gen_site = as.factor(gen_site))

# re-level and re-name treatment
calc_phys$treat <- factor(calc_phys$treat, levels = c("init", "1", "2", "3","4","5"))
levels(calc_phys$treat) <- c("Initial","Control","Low Var","Mod Var","High Var","Control 2")

# calculate average weight for each time point
# want mg/cm2/day (weight was recorded in grams)
# PA buoyant weight on 9.6.2016 & 9.7.2016
# T0 buoyant weight on 9.22.2016 & 9.23.2016
# T1 buoyant weight on 10.17.2016 & 10.18.2016
# T2 buoyant weight on 11.8.2016 & 11.9.2016 (taken just before heat stress)
# T3 buoyant weight on 12.10.2016 & 12.11.2016 (taken after stress, at end of recovery)
# T0-PA = 16 days, T1-T0 = 25 days, T2-T1 = 22 days, T3-T2 = 32 days
# T3-T0 =

# calculate average weight at each time point
calc_phys = calc_phys %>%
  mutate(paavgbw=rowMeans(.[ , c("pabwrec1","pabwrec2","pabwrec3")], na.rm=TRUE)) %>%
  mutate(t0avgbw=rowMeans(.[ , c("t0bwrec1","t0bwrec2","t0bwrec3")], na.rm=TRUE)) %>%
  mutate(t1avgbw=rowMeans(.[ , c("t1bwrec1","t1bwrec2","t1bwrec3")], na.rm=TRUE)) %>%
  mutate(t2avgbw=rowMeans(.[ , c("t2bwrec1","t2bwrec2","t2bwrec3")], na.rm=TRUE)) %>%
  mutate(t3avgbw=rowMeans(.[ , c("t3bwrec1","t3bwrec2","t3bwrec3")], na.rm=TRUE))

# calculate growth rate at the time point intervals - can't use this though because we do not have dry weights
calc_phys = calc_phys %>%
  mutate(T0_PA_g_cm2_day=((t0avgbw-paavgbw)/16)/SAcm2_t0) %>%
  mutate(T2_T0_g_cm2_day=((t2avgbw-t0avgbw)/57)/SAcm2_t2) %>%
  mutate(T2_T1_g_cm2_day=((t2avgbw-t1avgbw)/22)/SAcm2_t2) %>%
  mutate(T3_T2_g_cm_day2=((t3avgbw-t2avgbw)/32)/SAcm2_t3) %>%
  mutate(T3_T0_g_cm_day2=((t3avgbw-t0avgbw)/79)/SAcm2_t3)

# calculate percent change in growth
calc_phys = calc_phys %>%
  mutate(T0_PA_perc=((t0avgbw-paavgbw)/paavgbw)*100) %>%
  mutate(T1_T0_perc=((t1avgbw-t0avgbw)/t0avgbw)*100) %>%
  mutate(T2_T1_perc=((t2avgbw-t1avgbw)/t1avgbw)*100) %>%
  mutate(T3_T2_perc=((t3avgbw-t2avgbw)/t2avgbw)*100) %>%  # this is from before to after stress
  mutate(T2_T0_perc=((t2avgbw-t0avgbw)/t0avgbw)*100) %>%  # this is time 0 til before the stress
  mutate(T3_T0_perc=((t3avgbw-t0avgbw)/t0avgbw)*100)      # this is time 0 til after the stress

# calculate relative growth rate (RGR)
calc_phys = calc_phys %>%
  mutate(T0_PA_rgr=((ln(t0avgbw)-ln(paavgbw))/16)) %>%
  mutate(T2_T0_rgr=((ln(t2avgbw)-ln(t0avgbw))/57)) %>%
  mutate(T2_T1_rgr=((ln(t2avgbw)-ln(t1avgbw))/22)) %>%
  mutate(T3_T2_rgr=((ln(t3avgbw)-ln(t2avgbw))/32)) %>%
  mutate(T3_T0_rgr=((ln(t3avgbw)-ln(t0avgbw))/79))

# take a look at the dataset
str(calc_phys)

# un-comment the drop_na() corresponding to the time point of data you want to look at.
calc_phys2 = calc_phys %>%
  #drop_na(T3_T2_perc) %>%
  #drop_na(T2_T0_perc) %>%
  #drop_na(T2_T0_g_cm2_day) %>%
  #drop_na(T2_T0_rgr) %>%
  drop_na(T3_T2_rgr) %>%
  dplyr::filter(treat!="Control 2") %>%
  dplyr::filter(gen_site != "I4G") # clone with I4F, remove from dataset

# merge with lineage info for later plotting
lineages = read.csv("/Users/hannahaichelman/Documents/BU/TVE/2bRAD/Analysis/tuftscustompipeline_denovo_nosyms/tve_lineages_noclones.csv")

calc_phys_all_lin <- left_join(calc_phys2, lineages, by = "gen_site")
calc_phys_all_lin$lineage = as.factor(calc_phys_all_lin$lineage)

## combine with dominant symbiont type info
# merge with its2 divs
its2_divs = read.csv("/Users/hannahaichelman/Documents/BU/TVE/16S_ITS2/ITS_PreStress_Timepoint/ITS2.dominantDIVs.csv") %>%
  select(frag, dominant_div)

calc_phys_all_lin <- left_join(calc_phys_all_lin, its2_divs, by = "frag")

# merge with majority its type info
its2_types = read.csv("/Users/hannahaichelman/Documents/BU/TVE/16S_ITS2/ITS_PreStress_Timepoint/ITS2.dominanttype.csv") %>%
  select(frag, dominant_type)

calc_phys_all_lin <- left_join(calc_phys_all_lin, its2_types, by = "frag")
head(calc_phys_all_lin)
calc_phys_all_lin$dominant_type = as.factor(calc_phys_all_lin$dominant_type)

calc_phys_2_lin = calc_phys_all_lin %>%
  dplyr::filter(is.na(lineage) | lineage!="L3") # want to keep NA values for lineage here since they still have other info, will remove na's for lineage specific plots

# exploratory figure
#scatter plot with linear regression and confidence interval
ggplot(calc_phys, aes(treat, T2_T0_rgr, color = sitename))+
  geom_point()+
  geom_smooth(aes(group=sitename), method=lm)+
  theme_classic()


## Mixed Model
# interested in the effects of dtv treatment and lineage

m1 <- lmer(T3_T2_rgr ~ treat+lineage + (1|gen_site), data = calc_phys_2_lin, REML=TRUE)
summary(m1)
# T2-T0 percent change:
# Fixed effects:
#                   Estimate Std. Error       df t value Pr(>|t|)
#   (Intercept)     1.3703     0.2375  85.1334   5.770 1.25e-07 ***
#   treatLow Var    0.2880     0.2180 112.2567   1.321  0.18912
#   treatMod Var    0.5803     0.2073 111.5431   2.799  0.00605 **
#   treatHigh Var   0.5806     0.2135 112.3344   2.719  0.00758 **
#   lineageL2      -0.6754     0.3202  45.2749  -2.109  0.04050 *

# T3-T2 percent change:
# Fixed effects:
#                  Estimate Std. Error        df t value Pr(>|t|)
#   (Intercept)     2.30217    0.24044  94.95830   9.575 1.34e-15 ***
#   treatLow Var    0.08456    0.24913 108.96658   0.339 0.734936
#   treatMod Var    0.28463    0.23236 107.33329   1.225 0.223266
#   treatHigh Var   0.24759    0.23721 107.49086   1.044 0.298935
#   lineageL2      -1.18362    0.30972  45.22985  -3.822 0.000403 ***

# T2-T0 RGR:
#                   Estimate Std. Error         df t value Pr(>|t|)
# (Intercept)    2.373e-04  4.102e-05  8.510e+01   5.786 1.17e-07 ***
# treatLow Var   5.005e-05  3.764e-05  1.123e+02   1.329  0.18638
# treatMod Var   9.971e-05  3.580e-05  1.115e+02   2.785  0.00629 **
# treatHigh Var  9.946e-05  3.686e-05  1.123e+02   2.698  0.00804 **
# lineageL2     -1.161e-04  5.531e-05  4.528e+01  -2.099  0.04140 *

# T3-T2 RGR:
# Fixed effects:
#                 Estimate Std. Error         df t value Pr(>|t|)
# (Intercept)    7.091e-04  7.355e-05  9.495e+01   9.642 9.67e-16 ***
# treatLow Var   2.575e-05  7.622e-05  1.089e+02   0.338 0.736157
# treatMod Var   8.626e-05  7.109e-05  1.073e+02   1.213 0.227660
# treatHigh Var  7.325e-05  7.258e-05  1.075e+02   1.009 0.315112
# lineageL2     -3.620e-04  9.473e-05  4.520e+01  -3.821 0.000403 ***

Anova(m1)
# Response: T2_T0_perc
#           Chisq Df Pr(>Chisq)
# treat   10.3300  3    0.01596 *
# lineage  4.4484  1    0.03493 *

# Response: T2_T0_rgr
# Chisq Df Pr(>Chisq)
# treat   10.1798  3     0.0171 *
# lineage  4.4067  1     0.0358 *

# Response: T3_T2_rgr
# Chisq Df Pr(>Chisq)
# treat    1.8694  3  0.5999525
# lineage 14.6009  1  0.0001329 ***

# Response: T3_T2_perc
#           Chisq Df Pr(>Chisq)
# treat    1.9426  3  0.5844053
# lineage 14.6046  1  0.0001326 ***

# now let's do some more looking into the model
plot_model(m1, "eff", terms="treat")
plot_model(m1, "eff", terms="lineage")
plot_model(m1, type="re")

#Check model fit
r2(m1) #get r2 and adjusted r2
check_model(m1) #check assumptions and model fit

# Now look at custom contrasts with emmmeans

#specify model (because we are interested in pairwise, have to include the interaction)
m.emm<- lmer(T3_T2_rgr ~ treat*lineage + (1|gen_site), data = calc_phys_2_lin, REML=FALSE)

emms<-emmeans(m.emm, ~lineage|treat) #, adjust="Bonferoni"
pairs(emms, interaction = "pairwise") %>% rbind(adjust="fdr")
# for T2_T0_perc
# treat    lineage_pairwise estimate    SE  df t.ratio p.value
# Control  L1 - L2             0.468 0.420 116   1.115  0.2673
# Low Var  L1 - L2             0.831 0.431 122   1.927  0.1125
# Mod Var  L1 - L2             0.787 0.403 106   1.954  0.1125
# High Var L1 - L2             0.613 0.415 113   1.478  0.1896

# for T2_T0_RGR
# treat    lineage_pairwise estimate       SE  df t.ratio p.value
# Control  L1 - L2          0.000081 7.26e-05 116   1.117  0.2664
# Low Var  L1 - L2          0.000143 7.45e-05 122   1.926  0.1130
# Mod Var  L1 - L2          0.000135 6.96e-05 106   1.942  0.1130
# High Var L1 - L2          0.000105 7.17e-05 113   1.461  0.1958

# for T3_T2_perc:
# treat    lineage_pairwise estimate    SE  df t.ratio p.value
# Control  L1 - L2              1.09 0.444 133   2.455  0.0205
# Low Var  L1 - L2              1.06 0.462 140   2.297  0.0231
# Mod Var  L1 - L2              1.19 0.414 120   2.886  0.0092
# High Var L1 - L2              1.35 0.423 124   3.194  0.0071

# for T3_T2_RGR:
# treat    lineage_pairwise estimate       SE  df t.ratio p.value
# Control  L1 - L2          0.000336 0.000136 133   2.473  0.0195
# Low Var  L1 - L2          0.000325 0.000141 140   2.298  0.0230
# Mod Var  L1 - L2          0.000364 0.000127 120   2.876  0.0095
# High Var L1 - L2          0.000412 0.000129 124   3.185  0.0073


# Make growth figures in manuscript:
#SummarySE to format data for plotting by treatment and lineage

calc_phys_2_lin_nona = calc_phys_2_lin %>%
  drop_na(lineage)

#growth_means_2_lin <- summarySE(calc_phys_2_lin_nona, measurevar="T2_T0_rgr", groupvars=c("treat","lineage"))
growth_means_2_lin <- summarySE(calc_phys_2_lin_nona, measurevar="T3_T2_rgr", groupvars=c("treat","lineage"))

# plot, treatment x axis colored by lineage data figure
calc_plot_lineage <- ggplot(calc_phys_2_lin_nona,aes(x = treat, y = T3_T2_rgr))+
  theme_bw()+
  geom_jitter(aes(color = lineage, fill = lineage),
              position=position_dodge(width=0.3),
              alpha=0.2, pch = 21,
              color = "black") +
  geom_errorbar(data = growth_means_2_lin, aes(x = treat, ymax = T3_T2_rgr+se, ymin = T3_T2_rgr-se, color = lineage), width = .2, position = position_dodge(width=0.4)) +
  geom_point(data = growth_means_2_lin, mapping = aes(x=treat, y=T3_T2_rgr, color = lineage, fill = lineage), size = 3.5, pch = 21, color = "black", position = position_dodge(width=0.4))+
  scale_fill_manual(name = "Lineage",
                    breaks = c("L1","L2"),
                    values = cols_lineage)+
  scale_color_manual(name = "Lineage",
                     breaks = c("L1","L2"),
                     values = cols_lineage)+
  xlab("Treatment")+
  ylab(bquote("Specific Growth Rate ("~day^-1~')'))+
  geom_hline(yintercept=0, linetype='dotted', color = 'gray')+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))
calc_plot_lineage

ggsave(calc_plot_lineage, filename = "/Users/hannahaichelman/Documents/BU/TVE/Growth/RGR_T3T2_lineage_2_lin.pdf", width=5, height=4, units=c("in"), useDingbats=FALSE)


## More exploratory figures below
#SummarySE to format data for plotting by treatment and lineage_symtype
calc_phys_2_lin_nona = calc_phys_2_lin %>%
  drop_na(lineage,dominant_type)

# Plot based on dominant symbiont type
growth_means_lineage_sym <- summarySE(calc_phys_2_lin_nona, measurevar="T2_T0_perc", groupvars=c("treat","lineage","dominant_type"))

# plot, treatment x axis colored by site data figure
calc_plot_lineage_sym <- ggplot(growth_means_lineage_sym,aes(x = treat, y = T2_T0_perc, color = dominant_type))+
  theme_bw()+
  geom_point(size = 3, position = position_dodge(width=0.3))+
  geom_smooth(aes(group=dominant_type), position = position_dodge(width=0.3), size = 0.5, method='lm', se = FALSE)+
  #geom_line(aes(group = sitename), size = 0.5, linetype="dashed", position = position_dodge(width=0.3))+
  geom_errorbar(aes(x = treat, ymax = T2_T0_perc+se, ymin = T2_T0_perc-se), width = .2, position = position_dodge(width=0.3)) +
  #scale_color_manual(name = "Lineage-Sym",
  #                   values = c("#3f007d","#3f007d","#807dba","#807dba"))+
  #scale_shape_manual(name = "Lineage-Sym",
  #                   values=c(19,17,19,17))+
  xlab("Treatment")+
  ylab("Percent Change in Weight")+
  #ylab(bquote("Calcification (g" ~cm^-2~ ~day^-1~')'))+
  #ggtitle("Growth during variability")+
  #ylim(-0.5,4) +
  geom_hline(yintercept=0, linetype='dotted', color = 'gray')+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  facet_wrap(~lineage)
calc_plot_lineage_sym

ggsave(calc_plot_lineage_sym, filename = "/Users/hannahaichelman/Documents/BU/TVE/Growth/calcification_T2T0_lineage_dominantType.pdf", width=8, height=5, units=c("in"), useDingbats=FALSE)

# plot just CI corals to see if lineage difference holds
calc_phys_2_lin_nona_CI = calc_phys_2_lin_nona %>%
  subset(sitename == "CI")

growth_means_lineage_CI <- summarySE(calc_phys_2_lin_nona_CI, measurevar="T2_T0_g_cm2_day", groupvars=c("treat","lineage"))

calc_plot_lineage_CI <- ggplot(growth_means_lineage_CI,aes(x = treat, y = T2_T0_g_cm2_day, color = lineage, pch = lineage))+
  theme_bw()+
  geom_point(size = 3, position = position_dodge(width=0.3))+
  geom_smooth(aes(group=lineage), position = position_dodge(width=0.3), size = 0.5, method='lm', se = FALSE)+
  #geom_line(aes(group = sitename), size = 0.5, linetype="dashed", position = position_dodge(width=0.3))+
  geom_errorbar(aes(x = treat, ymax = T2_T0_g_cm2_day+se, ymin = T2_T0_g_cm2_day-se), width = .2, position = position_dodge(width=0.3)) +
  scale_color_manual(name = "Lineage",
                     breaks = c("L1","L2"),
                     values = cols_lineage)+
  scale_shape_manual(name = "Lineage",
                     breaks = c("L1","L2"),
                     values=c(19,17))+
  xlab("Treatment")+
  #ylab("Percent Change in Weight")+
  ylab(bquote("Calcification (g" ~cm^-2~ ~day^-1~')'))+
  #ggtitle("Growth during variability")+
  #ylim(-0.5,4) +
  geom_hline(yintercept=0, linetype='dotted', color = 'gray')+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))
calc_plot_lineage_CI

ggsave(calc_plot_lineage_CI, filename = "/Users/hannahaichelman/Documents/BU/TVE/Growth/calcificationo_T2T0_lineage_CIonly.pdf", width=5, height=4, units=c("in"), useDingbats=FALSE)

#SummarySE to format data for plotting by treatment and site
growth_means_site_all_lin <- summarySE(calc_phys_all_lin, measurevar="T2_T0_perc", groupvars=c("treat","sitename"))
growth_means_site_2_lin <- summarySE(calc_phys_2_lin, measurevar="T2_T0_perc", groupvars=c("treat","sitename"))

# plot, treatment x axis colored by site data figure
calc_plot_site <- ggplot(growth_means_site_2_lin,aes(x = treat, y = T2_T0_perc, color = sitename, pch = sitename))+
  theme_bw()+
  geom_point(size = 3, position = position_dodge(width=0.3))+
  geom_smooth(aes(group=sitename), position = position_dodge(width=0.3), size = 0.5, method='lm', se = FALSE)+
  #geom_line(aes(group = sitename), size = 0.5, linetype="dashed", position = position_dodge(width=0.3))+
  geom_errorbar(aes(x = treat, ymax = T2_T0_perc+se, ymin = T2_T0_perc-se), width = .2, position = position_dodge(width=0.3)) +
  scale_color_manual(name = "Site",
                     labels = c("CI","PD","SP","BN","BS","CA"),
                     values = cols_site)+
  scale_shape_manual(name = "Site",
                     labels = c("CI","PD","SP","BN","BS","CA"),
                     values=c(19,19,19,17,17,17))+
  xlab("Treatment")+
  ylab("Percent Change in Weight")+
  #ylab(bquote("Calcification (g" ~cm^-2~ ~day^-1~')'))+
  #ggtitle("Growth during variability")+
  #ylim(-0.5,4) +
  geom_hline(yintercept=0, linetype='dotted', color = 'gray')+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))
calc_plot_site

ggsave(calc_plot_site, filename = "/Users/hannahaichelman/Documents/BU/TVE/Growth/calcification_T2T0_site_2_lin.pdf", width=5, height=4, units=c("in"), useDingbats=FALSE)

#plot, faceted by treatment to looking for outliers
plot <- ggplot(calc_phys,aes(x = sitename, y = T2_T0_perc, color = sitename, pch = sitename, text = paste("Nubbin:",frag)))+
  theme_bw()+
  geom_boxplot()+
  ggtitle("Percent Change in Weight T0-T2")+
  scale_fill_manual(values = palsite)+
  geom_jitter(aes(color = sitename),width = 0.3)+
  scale_color_manual(values = palsite)+
  theme(axis.title.x = element_blank())+
  theme(axis.text.x=element_text(angle = 45, vjust = 0.5))+
  facet_wrap(~treat, nrow=5)
#  labs(y=expression(paste("Protein (mg cm"^-2*' )')))
ggplotly(subplot(list(plot),titleY=F) %>% layout(showlegend=F))


# manipulate the data to create a data frame we can use to look at the time course of percent change in growth
calc_phys_timecourse = calc_phys %>%
  select(frag, treat, gen_site, sitename, reef,
         T0_PA_perc,T1_T0_perc,T2_T1_perc,T3_T2_perc)

head(calc_phys_timecourse)

calc_phys_timecourse.melt <- melt(calc_phys_timecourse, id.vars = c('frag','treat','gen_site','reef','sitename'), measure.vars = c('T0_PA_perc','T1_T0_perc','T2_T1_perc','T3_T2_perc'))
head(calc_phys_timecourse.melt)
str(calc_phys_timecourse.melt)

colnames(calc_phys_timecourse.melt)[6] <- "time"
levels(calc_phys_timecourse.melt$time)
levels(calc_phys_timecourse.melt$time) <- c("PA-T0","T0-T1","T1-T2","T2-T3")

# Stats for time course changes:
# first try with lmer()
phys_calc_time_nona = calc_phys_timecourse.melt %>%
  drop_na(value) %>%
  filter(treat!="Control 2")

model.calc.time <- lmer(value ~ treat + sitename + time + (1|gen_site), data = phys_calc_time_nona)
summary(model.calc.time)

# SummarySE and plot for time course
growth_time_means <- summarySE(phys_calc_time_nona, measurevar="value", groupvars=c("treat","sitename","time"))

# plot, treatment x axis colored by site data figure
plot.time <- ggplot(growth_time_means,aes(x = time, y = value, color = sitename, pch = treat))+
  theme_bw()+
  geom_point(size = 3, position = position_dodge(width=0.3))+
  #geom_smooth(aes(group=sitename), position = position_dodge(width=0.3), size = 0.5, method='lm', se = FALSE)+
  #geom_line(aes(group = sitename), size = 0.5, linetype="dashed", position = position_dodge(width=0.3))+
  geom_errorbar(aes(x = time, ymax = value+se, ymin = value-se), width = .2, position = position_dodge(width=0.3)) +
  scale_color_manual(name = "Site",
                     labels = c("I-Cristobal","I-Punta Donato","I-STRI Point","O-Bastimentos N","O-Bastimentos S","O-Cayo de Agua"),
                     values = cols_site)+
  scale_shape_manual(name = "Treatment",
                     labels = c("Control 1","Low Var","Mod Var","High Var"),
                     values=c(15,16,17,18))+
  xlab("Time")+
  ylab("Percent Change in Weight")+
  #ggtitle("Growth during variability")+
  #ylim(-0.5,4) +
  geom_hline(yintercept=0, linetype='dotted', color = 'gray')+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))
plot.time
ggsave(file="/Users/hannahaichelman/Documents/BU/TVE/Growth/PercentChange_timecourse.pdf", plot.time, width=8, height=5, units=c("in"), useDingbats=FALSE)

# plot, treatment x axis colored by site data figure
plot.time.2 <- ggplot(growth_time_means,aes(x = time, y = value, color = treat, pch = sitename))+
  theme_bw()+
  geom_point(size = 3, position = position_dodge(width=0.3))+
  geom_errorbar(aes(x = time, ymax = value+se, ymin = value-se), width = .2, position = position_dodge(width=0.3)) +
  scale_color_manual(name = "Treatment",
                     labels = c("Control","Low Var","Mod Var","High Var"),
                     values = cols_treat)+
  scale_shape_manual(name = "Site",
                     labels = c("CI","PD","SP","BN","BS","CA"),
                     values=c(15,16,17,18,4,8))+
  xlab("Time")+
  ylab("Percent Change in Weight")+
  #ggtitle("Growth during variability")+
  #ylim(-0.5,4) +
  geom_hline(yintercept=0, linetype='dotted', color = 'gray')+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  facet_wrap(~sitename, nrow=2, ncol=3)
plot.time.2
ggsave(file="/Users/hannahaichelman/Documents/BU/TVE/Growth/PercentChange_timecourse_treatfacet.pdf", plot.time.2, width=8, height=5, units=c("in"), useDingbats=FALSE)


## Plot 90% daily temperature range by % change in growth
# an estimate of the magnitude of high-frequency temp fluctuations

# logger          percent90
#  Punta.IR1 (PD)           0.999
#  STRI.IR2  (SP)           0.833
#  Cristo.IR3.arr1 (CI)     1.92
#  Cayo.OR3.arr3 (CA)      1.49
#  Drago.OR4 (DM)          1.24

# add in column for the 90% DTR
phys_calc_nona_lineage$dtr_90 <- ifelse(phys_calc_nona_lineage$sitename == 'PD', '0.999',
                                        ifelse(phys_calc_nona_lineage$sitename == 'SP', '0.833',
                                               ifelse(phys_calc_nona_lineage$sitename == 'CI', '1.92',
                                                      ifelse(phys_calc_nona_lineage$sitename == 'CA', '1.49',
                                                             'NA'))))
phys_calc_nona_lineage$dtr_90 <- as.numeric(phys_calc_nona_lineage$dtr_90)
phys_calc_nona_lineage_dtr = phys_calc_nona_lineage %>%
  drop_na(dtr_90)

#SummarySE to format data for plotting by treatment and site
growth_means_dtr <- summarySE(phys_calc_nona_lineage_dtr, measurevar="T2_T0_perc", groupvars=c("dtr_90","treat"))

# plot, treatment x axis colored by site data figure

calc_plot_dtr <- ggplot(growth_means_dtr,aes(x = dtr_90, y = T2_T0_perc, color = treat, pch = treat))+
  theme_bw()+
  geom_point(size = 3, position = position_dodge(width=0.1))+
  geom_smooth(aes(group=treat), position = position_dodge(width=0.1), size = 0.5, method='lm', se = FALSE)+
  #geom_line(aes(group = sitename), size = 0.5, linetype="dashed", position = position_dodge(width=0.3))+
  geom_errorbar(aes(x = dtr_90, ymax = T2_T0_perc+se, ymin = T2_T0_perc-se), width = .2, position = position_dodge(width=0.1)) +
  scale_color_manual(name = "Treatment",
                     labels = c("Control","Low Var","Mod Var","High Var"),
                     values = cols_treat)+
  scale_shape_manual(name = "Treatment",
                     labels = c("Control","Low Var","Mod Var","High Var"),
                     values=c(15,16,17,18))+
  xlab("90% Daily Temperature Range (°C)")+
  ylab("Percent Change in Weight")+
  #ggtitle("Growth during variability")+
  #ylim(-0.5,4) +
  geom_hline(yintercept=0, linetype='dotted', color = 'gray')+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))
calc_plot_dtr

ggsave(calc_plot_dtr, filename = "/Users/hannahaichelman/Documents/BU/TVE/Growth/percentchange_T2T0_dtr.pdf", width=5, height=4, units=c("in"), useDingbats=FALSE)

# seems to not be an effect of reef zone so let's check that out
library(jtools)
m.dtr <- lmer(T2_T0_perc ~ dtr_90 + treat + (1|gen_site), data = phys_calc_nona_lineage_dtr, REML=TRUE)
summ(m.dtv)


##### PAM #####
# need to re-read in the data sheet for pam only until I can organize better
post_phys_forpam <- read.csv('dtvmaster.csv')

# analyzing separately because not sure of timing yet
pam_phys = post_phys_forpam %>%
  dplyr::select(frag, survivedtoend, treat, t3sastan1, t3sastan2, t3sastan3, t3sarec1, t3sarec2, t3sarec3, blastvol, blaster,
                papamrec1,papamrec2,papamrec3,t0pamrec1,t0pamrec2,t0pamrec3,t1pamrec1,t1pamrec2,t1pamrec3,t2pamrec1,t2pamrec2,t2pamrec3,t3pamrec1,t3pamrec2,t3pamrec3,
                t4pamrec1,t4pamrec2,t4pamrec3,t5pamrec1,t5pamrec2,t5pamrec3,t6pamrec1,t6pamrec2,t6pamrec3,t7pamrec1,t7pamrec2,t7pamrec3,t8pamrec1,t8pamrec2,t8pamrec3,t9pamrec1,t9pamrec2,t9pamrec3) %>%
  dplyr::rename(sastan1 = t3sastan1, sastan2 = t3sastan2, sastan3 = t3sastan3, sarec1 = t3sarec1, sarec2 = t3sarec2, sarec3 = t3sarec3) %>%
  mutate(treat = as.factor(treat), sarec3 = as.numeric(sarec3), blastvol = as.numeric(blastvol)) %>%
  dplyr::filter(frag!="I2A12", frag!="I2A7", frag!="O4D8", frag!="I3E10", frag!="I3I10") #these frags are being removed because they were duplicate genotypes within treatment and had the most complete information of the two

# add identifying data
pam_phys$origsitecode <- substr(pam_phys$frag, 1, 2)

# add in site name with inshore/offshore indicator
pam_phys$sitename <- ifelse(pam_phys$origsitecode == 'I2', 'SP',
                            ifelse(pam_phys$origsitecode == 'I3', 'CI',
                                   ifelse(pam_phys$origsitecode == 'I4', 'PD',
                                          ifelse(pam_phys$origsitecode == 'O2', 'BS',
                                                 ifelse(pam_phys$origsitecode == 'O3', 'CA',
                                                        'BN')))))
pam_phys$sitename <- as.factor(pam_phys$sitename)

# make new nubbin IDs based on the new, more informative site codes
# extract the genotype and frag number
pam_phys$fragid <- substr(pam_phys$frag,3,5)

# add inshore/offshore designation
pam_phys$reef <- substr(pam_phys$frag,1,1)
pam_phys$reef <- ifelse(pam_phys$reef == 'O', 'Outer Reef', 'Inner Reef')
pam_phys$reef <- as.factor(pam_phys$reef)
pam_phys$genet <- substr(pam_phys$fragid,1,1)

#create a new column of combined genotype and site for stats later
pam_phys = pam_phys %>%
  unite(gen_site, c(origsitecode,genet), sep = "", remove = FALSE) %>%
  mutate(gen_site = as.factor(gen_site))

# re-level and re-name treatment
pam_phys$treat <- factor(pam_phys$treat, levels = c("init", "1", "2", "3","4","5"))
levels(pam_phys$treat) <- c("Initial","Control","Low Var","Mod Var","High Var","Control 2")
pam_phys$treat <- as.factor(pam_phys$treat)

# calculate average fv/fm for each time point
# PA PAM on 9.2.2016 & 9.3.2016 (day -18)
# T0 PAM on 9.20.2016 (day 0)
# T1 PAM on 10.5.2016 (day 15)
# T2 PAM on 10.25.2016 (day 35)
# T3 PAM on 11.5.2016 (taken 4 days before start of heat stress ramp) (day 45)
# T4 PAM on 11.14.2016 (first day corals were at 32°C) (day 54) - temp increase started day 51
# T5 PAM on 11.21.2016 (end of 32 hold, which lasted 8 days) (day 61)
# T6 PAM on 11.25.2016 (after cool down ramp, back at 29.5°C) (day 65)
# T7 PAM on 11.30.2016 (on day 5 of recovery at 29.5°C) (day 70)
# T8 PAM on 12.4.2016 (on day 9 of recovery at 29.5°C) (day 74)
# T9 PAM on 12.9.2016 (on day 14 of recovery at 29.5°C) (day 79)
# the experiment ended on 12.13.2016

# calculate average PAM at each time point
pam_phys = pam_phys %>%
  mutate(paavgpam=rowMeans(.[ , c("papamrec1","papamrec2","papamrec3")], na.rm=TRUE)) %>%
  mutate(t0avgpam=rowMeans(.[ , c("t0pamrec1","t0pamrec2","t0pamrec3")], na.rm=TRUE)) %>%
  mutate(t1avgpam=rowMeans(.[ , c("t1pamrec1","t1pamrec2","t1pamrec3")], na.rm=TRUE)) %>%
  mutate(t2avgpam=rowMeans(.[ , c("t2pamrec1","t2pamrec2","t2pamrec3")], na.rm=TRUE)) %>%
  mutate(t3avgpam=rowMeans(.[ , c("t3pamrec1","t3pamrec2","t3pamrec3")], na.rm=TRUE)) %>%
  mutate(t4avgpam=rowMeans(.[ , c("t4pamrec1","t4pamrec2","t4pamrec3")], na.rm=TRUE)) %>%
  mutate(t5avgpam=rowMeans(.[ , c("t5pamrec1","t5pamrec2","t5pamrec3")], na.rm=TRUE)) %>%
  mutate(t6avgpam=rowMeans(.[ , c("t6pamrec1","t6pamrec2","t6pamrec3")], na.rm=TRUE)) %>%
  mutate(t7avgpam=rowMeans(.[ , c("t7pamrec1","t7pamrec2","t7pamrec3")], na.rm=TRUE)) %>%
  mutate(t8avgpam=rowMeans(.[ , c("t8pamrec1","t8pamrec2","t8pamrec3")], na.rm=TRUE)) %>%
  mutate(t9avgpam=rowMeans(.[ , c("t9pamrec1","t9pamrec2","t9pamrec3")], na.rm=TRUE)) %>%
  mutate(pamdiff=t9avgpam-t3avgpam)

# explore difference in PAM values from the end of recovery period to start of the heat stress ramp
str(pam_phys)
plot(x=pam_phys$sitename,y=pam_phys$pamdiff, data = pam_phys)

lm02 = aov(pamdiff ~ sitename + treat, data = pam_phys)
summary(lm02)
par(mfrow=c(2,2))
plot(lm02)
TukeyHSD(lm02)
# not a lot interesting happening here, on to stats and plotting

# summarySE doesn't work with NA's, so use this!!
phys_pam_wide = pam_phys %>%
  dplyr::filter(complete.cases(paavgpam,t0avgpam,t1avgpam,t2avgpam,t3avgpam,t4avgpam,t5avgpam,t6avgpam,t7avgpam,t8avgpam,t9avgpam)) %>% #drop any row that has an NA for any time point
  dplyr::filter(treat!="Control 2") %>%
  select(frag,treat,reef,gen_site,sitename,paavgpam,t0avgpam,t1avgpam,t2avgpam,t3avgpam,t4avgpam,t5avgpam,t6avgpam,t7avgpam,t8avgpam,t9avgpam) %>%
  dplyr::filter(gen_site != "I4G") # remove clone

# transform the data to long format so time point is its own column
phys_pam_long = phys_pam_wide %>%
  gather(time, pam, paavgpam:t9avgpam)

# re-level and re-name treatment
phys_pam_long$time <- as.factor(phys_pam_long$time)
levels(phys_pam_long$time) <- c("-18","0", "15","35","45","54","61","65","70","74","79")

# merge with lineage info for later plotting
lineages = read.csv("/Users/hannahaichelman/Documents/BU/TVE/2bRAD/Analysis/tuftscustompipeline_denovo_nosyms/tve_lineages_noclones.csv")

pam_phys_all_lin <- left_join(phys_pam_long, lineages, by = "gen_site")
pam_phys_all_lin$lineage = as.factor(pam_phys_all_lin$lineage)

# merge with its2 types for plotting
its2_types = read.csv("/Users/hannahaichelman/Documents/BU/TVE/16S_ITS2/ITS_PreStress_Timepoint/ITS2.dominanttype.csv") %>%
  select(frag, dominant_type)

pam_phys_all_lin <- left_join(pam_phys_all_lin, its2_types, by = "frag")
pam_phys_all_lin$dominant_type = as.factor(pam_phys_all_lin$dominant_type)

# merge with its2 divs
its2_divs = read.csv("/Users/hannahaichelman/Documents/BU/TVE/16S_ITS2/ITS_PreStress_Timepoint/ITS2.dominantDIVs.csv") %>%
  select(frag, dominant_div)

pam_phys_all_lin <- left_join(pam_phys_all_lin, its2_divs, by = "frag")

# filter to only include T3-T9 for plotting and stats
phys_pam_all_lin_plots = pam_phys_all_lin %>%
  dplyr::filter(time != "-18" & time != "0" & time != "15" & time != "35")

phys_pam_2_lin_plots = phys_pam_all_lin_plots %>%
  dplyr::filter(is.na(lineage) | lineage!="L3") # want to keep NA values for lineage here since they still have other info, will remove na's for lineage specific plots

phys_pam_lin1 = phys_pam_2_lin_plots %>%
  dplyr::filter(lineage == "L1")

phys_pam_lin2 = phys_pam_2_lin_plots %>%
  dplyr::filter(lineage == "L2")

# Stats
# Using the full model with interactions of all our parameters of interest
m.full <- lmer(pam ~ time*lineage*dominant_type*treat + (1|gen_site), data = phys_pam_2_lin_plots, REML=TRUE)
summary(m.full)
Anova(m.full)
# report chi square and p-values to designate significance of factors
# Response: pam
#                                     Chisq Df Pr(>Chisq)
# time                             89.6885  6  < 2.2e-16 ***
# lineage                           7.5172  1   0.006111 **
# dominant_type                    39.0563  3  1.689e-08 ***
# treat                             7.8660  3   0.048864 *
# time:lineage                      8.1242  6   0.229142
# time:dominant_type               52.1631 18  3.532e-05 ***
# lineage:dominant_type             0.0432  1   0.835349
# time:treat                       34.5186 18   0.010861 *
# lineage:treat                     3.9443  3   0.267545
# dominant_type:treat              64.9150  7  1.564e-11 ***
# time:lineage:dominant_type       17.6541  6   0.007157 **
# time:lineage:treat                8.7538 18   0.965175
# time:dominant_type:treat         36.3252 42   0.717630
# lineage:dominant_type:treat       2.7817  3   0.426519
# time:lineage:dominant_type:treat 16.3648 18   0.567104
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# check interactions
m.emm<- lmer(pam ~ time*lineage + (1|gen_site), data = phys_pam_2_lin_plots, REML=FALSE)

emms<-emmeans(m.emm, ~lineage|time) #, adjust="Bonferoni"
pairs(emms, interaction = "pairwise") %>% rbind(adjust="fdr")
# time lineage_pairwise estimate     SE   df t.ratio p.value
# 45   L1 - L2            0.0196 0.0161 70.9   1.219  0.2270
# 54   L1 - L2            0.0368 0.0161 70.9   2.286  0.0590
# 61   L1 - L2            0.0500 0.0161 70.9   3.107  0.0095
# 65   L1 - L2            0.0537 0.0161 70.9   3.339  0.0094
# 70   L1 - L2            0.0296 0.0161 70.9   1.839  0.1185
# 74   L1 - L2            0.0281 0.0161 70.9   1.749  0.1185
# 79   L1 - L2            0.0218 0.0161 70.9   1.357  0.2089


#Check model fit
r2(m.full) #get r2 and adjusted r2
check_model(m.full) #check assumptions and model fit

# Now plot Fv/Fm data

#SummarySE to format data for plotting - dtv treatment
pam_means_treat_all_lin <- summarySE(phys_pam_all_lin_plots, measurevar="pam", groupvars=c("treat","time"))
pam_means_treat_all_lin$time <- as.numeric(as.character(pam_means_treat_all_lin$time))
pam_means_treat_2_lin <- summarySE(phys_pam_2_lin_plots, measurevar="pam", groupvars=c("treat","time"))
pam_means_treat_2_lin$time <- as.numeric(as.character(pam_means_treat_2_lin$time))

pam_plot_treatment <- ggplot(pam_means_treat_2_lin,aes(x = time, y = pam, color = treat, fill = treat))+
  theme_bw()+
  annotate("rect", xmin = 51, xmax = 64, ymin = - Inf, ymax = Inf, fill = "red4", alpha = 0.15)+
  annotate("rect", xmin = 64, xmax = 80, ymin = - Inf, ymax = Inf, fill = "royalblue4", alpha = 0.15)+
  geom_errorbar(aes(x = time, ymax = pam+se, ymin = pam-se), width = .2, position = position_dodge(width=0.3)) +
  geom_point(size = 3, pch = 21, color = "black", position = position_dodge(width=0.3))+
  geom_line(aes(group = treat), size = 0.5, linetype="dashed", position = position_dodge(width=0.3))+
  scale_fill_manual(name = "Treatment",
                    labels = c("Control","Low Var","Mod Var","High Var"),
                    values = cols_treat)+
  scale_color_manual(name = "Treatment",
                     labels = c("Control","Low Var","Mod Var","High Var"),
                     values = cols_treat)+
  ylab("Fv/Fm")+
  scale_x_continuous(name = "Time Point", breaks = c(45,54,61,65,70,74,79))
#ylim(0.4,0.7) +
pam_plot_treatment

ggsave(pam_plot_treatment, file="/Users/hannahaichelman/Documents/BU/TVE/PAM/PAM_treat_2_lin.pdf", width=6, height=4, units=c("in"), useDingbats=FALSE)

#SummarySE to format data for plotting - dtv treatment faceted by lineage
pam_means_treat_all_lin <- summarySE(phys_pam_all_lin_plots, measurevar="pam", groupvars=c("treat","time","lineage"))
pam_means_treat_all_lin$time <- as.numeric(as.character(pam_means_treat_all_lin$time))
pam_means_treat_2_lin <- summarySE(phys_pam_2_lin_plots, measurevar="pam", groupvars=c("treat","time","lineage"))
pam_means_treat_2_lin$time <- as.numeric(as.character(pam_means_treat_2_lin$time))

pam_means_treat_2_lin = pam_means_treat_2_lin %>%
  drop_na(lineage)

pam_plot_treatment_lin <- ggplot(pam_means_treat_2_lin,aes(x = time, y = pam, color = treat, fill = treat))+
  theme_bw()+
  annotate("rect", xmin = 51, xmax = 64, ymin = - Inf, ymax = Inf, fill = "red4", alpha = 0.15)+
  annotate("rect", xmin = 64, xmax = 80, ymin = - Inf, ymax = Inf, fill = "royalblue4", alpha = 0.15)+
  geom_errorbar(aes(x = time, ymax = pam+se, ymin = pam-se), width = .2, position = position_dodge(width=0.3)) +
  geom_point(size = 3, pch = 21, color = "black", position = position_dodge(width=0.3))+
  geom_line(aes(group = treat), size = 0.5, linetype="dashed", position = position_dodge(width=0.3))+
  scale_fill_manual(name = "Treatment",
                    labels = c("Control","Low Var","Mod Var","High Var"),
                    values = cols_treat)+
  scale_color_manual(name = "Treatment",
                     labels = c("Control","Low Var","Mod Var","High Var"),
                     values = cols_treat)+
  ylab("Fv/Fm")+
  scale_x_continuous(name = "Time Point", breaks = c(45,54,61,65,70,74,79))+
  facet_wrap(~lineage)
#ylim(0.4,0.7) +
pam_plot_treatment_lin

ggsave(pam_plot_treatment_lin, file="/Users/hannahaichelman/Documents/BU/TVE/PAM/PAM_treat_lin_facet.pdf", width=6, height=4, units=c("in"), useDingbats=FALSE)


#SummarySE to format data for plotting - lineage
phys_pam_all_lin_plots_nona = phys_pam_all_lin_plots %>%
  drop_na(lineage)

phys_pam_2_lin_plots_nona = phys_pam_2_lin_plots %>%
  drop_na(lineage)

pam_means_all_lin <- summarySE(phys_pam_all_lin_plots_nona, measurevar="pam", groupvars=c("lineage","time"))
pam_means_all_lin$time <- as.numeric(as.character(pam_means_all_lin$time))
pam_means_2_lin <- summarySE(phys_pam_2_lin_plots_nona, measurevar="pam", groupvars=c("lineage","time"))
pam_means_2_lin$time <- as.numeric(as.character(pam_means_2_lin$time))

# plot, treatment x axis colored by site data figure
pam_plot_lineage <- ggplot(pam_means_2_lin,aes(x = time, y = pam, color = lineage, fill = lineage))+
  theme_bw()+
  annotate("rect", xmin = 51, xmax = 64, ymin = - Inf, ymax = Inf, fill = "red4", alpha = 0.15)+
  annotate("rect", xmin = 64, xmax = Inf, ymin = - Inf, ymax = Inf, fill = "royalblue4", alpha = 0.15)+
  geom_errorbar(aes(x = time, ymax = pam+se, ymin = pam-se), width = .2, color = "black", position = position_dodge(width=0.3)) +
  geom_line(aes(group = lineage), size = 1, linetype="dashed", position = position_dodge(width=0.3))+
  geom_point(size = 3.5, color = "black", pch = 21, position = position_dodge(width=0.3))+
  #geom_smooth(aes(group=sitename), position = position_dodge(width=0.3), size = 0.5, method='lm', se = FALSE)+
  scale_fill_manual(name = "Lineage",
                    breaks = c("L1","L2"),
                    #breaks = c("L1","L2","L3"),
                    values = cols_lineage)+
  scale_color_manual(name = "Lineage",
                     breaks = c("L1","L2"),
                     #breaks = c("L1","L2","L3"),
                     values = cols_lineage)+
  ylab("Fv/Fm")+
  scale_x_continuous(name = "Time Point", breaks = c(45,54,61,65,70,74,79)) +
  #ylim(-0.5,4) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))
pam_plot_lineage

ggsave(pam_plot_lineage, filename = "/Users/hannahaichelman/Documents/BU/TVE/PAM/pam_lineage_2.pdf", width=5, height=4, units=c("in"), useDingbats=FALSE)

# plot just CI corals to see if lineage difference holds
phys_pam_lineage_plots_CI = phys_pam_2_lin_plots_nona %>%
  subset(sitename == "CI")

pam_lineagemeans_CI = summarySE(phys_pam_lineage_plots_CI, measurevar="pam", groupvars=c("lineage","time"))
pam_lineagemeans_CI$time = as.numeric(as.character(pam_lineagemeans_CI$time))

# plot, time x axis, shape = symbiont type, color = lineage
pam_plot_lineage_CI <- ggplot(pam_lineagemeans_CI,aes(x = time, y = pam, color = lineage, fill = lineage))+
  theme_bw()+
  annotate("rect", xmin = 51, xmax = 64, ymin = - Inf, ymax = Inf, fill = "red4", alpha = 0.15)+
  annotate("rect", xmin = 64, xmax = Inf, ymin = - Inf, ymax = Inf, fill = "royalblue4", alpha = 0.15)+
  geom_point(size = 3, color = "black", pch = 21, position = position_dodge(width=0.3))+
  geom_line(aes(color = lineage), size = 0.5, linetype="dashed", position = position_dodge(width=0.3))+
  geom_errorbar(aes(x = time, ymax = pam+se, ymin = pam-se), width = .2, position = position_dodge(width=0.3)) +
  scale_fill_manual(name = "Lineage",
                    breaks = c("L1","L2"),
                    values = cols_lineage)+
  scale_color_manual(name = "Lineage",
                     breaks = c("L1","L2"),
                     values = cols_lineage)+
  ylab("Fv/Fm")+
  scale_x_continuous(name = "Time Point", breaks = c(45,54,61,65,70,74,79)) +
  #ylim(-0.5,4) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))
pam_plot_lineage_CI

ggsave(pam_plot_lineage_CI, filename = "/Users/hannahaichelman/Documents/BU/TVE/PAM/pam_lineage_CIonly.pdf", width=6, height=4, units=c("in"), useDingbats=FALSE)


#SummarySE to format data for plotting - dominant symtype
phys_pam_2_linsym_nona = phys_pam_2_lin_plots %>%
  drop_na(lineage, dominant_type)

phys_pam_all_linsym_nona = phys_pam_all_lin_plots %>%
  drop_na(lineage, dominant_type)

pam_2_lineagesymtype_means <- summarySE(phys_pam_2_linsym_nona, measurevar="pam", groupvars=c("time","dominant_type","lineage"))
pam_2_lineagesymtype_means$time <- as.numeric(as.character(pam_2_lineagesymtype_means$time))

pam_all_lineagesymtype_means <- summarySE(phys_pam_all_linsym_nona, measurevar="pam", groupvars=c("time","dominant_type","lineage"))
pam_all_lineagesymtype_means$time <- as.numeric(as.character(pam_all_lineagesymtype_means$time))

its2_cols_greens = c("#edf8e9","#a1d99b",  "#238b45", "#00441b")

# plot, time x axis colored by dominant sym type data figure
pam_plot_lineage_symtype = pam_2_lineagesymtype_means %>%
  #filter(lineage == "L2") %>%
  ggplot(aes(x = time, y = pam, fill = dominant_type, color = dominant_type, shape = dominant_type))+
  theme_bw()+
  annotate("rect", xmin = 51, xmax = 64, ymin = - Inf, ymax = Inf, fill = "red4", alpha = 0.15)+
  annotate("rect", xmin = 64, xmax = Inf, ymin = - Inf, ymax = Inf, fill = "royalblue4", alpha = 0.15)+
  geom_errorbar(aes(x = time, ymax = pam+se, ymin = pam-se), width = .2, color = "black", position = position_dodge(width=1.5)) +
  geom_line(aes(group = dominant_type), size = 1, linetype="dashed", position = position_dodge(width=1.5))+
  geom_point(size = 3.5, color = "black", position = position_dodge(width=1.5))+
  scale_color_manual(name = "Dominant ITS2",
                     values = its2_cols_greens)+
  scale_fill_manual(name = "Dominant ITS2",
                    values = its2_cols_greens)+
  scale_shape_manual(name = "Dominant ITS2",
                     values = c(21,22,23,24))+
  ylab("Fv/Fm")+
  scale_x_continuous(name = "Time Point", breaks = c(45,54,61,65,70,74,79)) +
  #ylim(-0.5,4) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  facet_wrap(~lineage)
pam_plot_lineage_symtype

ggsave(pam_plot_lineage_symtype, filename = "/Users/hannahaichelman/Documents/BU/TVE/PAM/pam_lineage_sym_2_lin_facet_dominanttype.pdf", width=8.5, height=4, units=c("in"), useDingbats=FALSE)

#SummarySE to format data for plotting - CI only dominant symtype
phys_pam_2_linsym_CI = phys_pam_2_linsym_nona %>%
  filter(sitename == "CI")

phys_pam_2_linsym_CI_means <- summarySE(phys_pam_2_linsym_CI, measurevar="pam", groupvars=c("time","dominant_type","lineage"))
phys_pam_2_linsym_CI_means$time <- as.numeric(as.character(phys_pam_2_linsym_CI_means$time))

its2_cols_greens_CI = c("#edf8e9", "#00441b")

# plot, time x axis colored by dominant sym type data figure
pam_plot_lineage_symtype = phys_pam_2_linsym_CI_means %>%
  #filter(lineage == "L2") %>%
  ggplot(aes(x = time, y = pam, fill = dominant_type, color = dominant_type, shape = dominant_type))+
  theme_bw()+
  annotate("rect", xmin = 51, xmax = 64, ymin = - Inf, ymax = Inf, fill = "red4", alpha = 0.15)+
  annotate("rect", xmin = 64, xmax = Inf, ymin = - Inf, ymax = Inf, fill = "royalblue4", alpha = 0.15)+
  geom_errorbar(aes(x = time, ymax = pam+se, ymin = pam-se), width = .2, color = "black", position = position_dodge(width=1.5)) +
  geom_line(aes(group = dominant_type), size = 1, linetype="dashed", position = position_dodge(width=1.5))+
  geom_point(size = 3.5, color = "black", position = position_dodge(width=1.5))+
  scale_color_manual(name = "Dominant ITS2",
                     values = its2_cols_greens_CI)+
  scale_fill_manual(name = "Dominant ITS2",
                    values = its2_cols_greens_CI)+
  scale_shape_manual(name = "Dominant ITS2",
                     values = c(21,24))+
  ylab("Fv/Fm")+
  scale_x_continuous(name = "Time Point", breaks = c(45,54,61,65,70,74,79)) +
  #ylim(-0.5,4) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  facet_wrap(~lineage)
pam_plot_lineage_symtype

ggsave(pam_plot_lineage_symtype, filename = "/Users/hannahaichelman/Documents/BU/TVE/PAM/pam_lineage_sym_2_lin_facet_dominanttype_CIonly.pdf", width=8.5, height=4, units=c("in"), useDingbats=FALSE)

#SummarySE to format data for plotting - dominant DIV
phys_pam_2_linsym_nona = phys_pam_2_lin_plots %>%
  drop_na(lineage, dominant_div)

phys_pam_all_linsym_nona = phys_pam_all_lin_plots %>%
  drop_na(lineage, dominant_div)

pam_2_lineagediv_means <- summarySE(phys_pam_2_linsym_nona, measurevar="pam", groupvars=c("time","dominant_div","lineage"))
pam_2_lineagediv_means$time <- as.numeric(as.character(pam_2_lineagediv_means$time))

pam_all_lineagediv_means <- summarySE(phys_pam_all_linsym_nona, measurevar="pam", groupvars=c("time","dominant_type","lineage"))
pam_all_lineagediv_means$time <- as.numeric(as.character(pam_all_lineagediv_means$time))

# plot, treatment x axis colored by site data figure
pam_plot_lineage_symtype <- ggplot(pam_2_lineagediv_means,aes(x = time, y = pam, color = dominant_div))+
  theme_bw()+
  annotate("rect", xmin = 46, xmax = 63, ymin = - Inf, ymax = Inf, fill = "red4", alpha = 0.15)+
  annotate("rect", xmin = 63, xmax = 80, ymin = - Inf, ymax = Inf, fill = "royalblue4", alpha = 0.15)+
  geom_point(size = 3, position = position_dodge(width=0.3))+
  geom_line(aes(group = dominant_div), size = 0.5, linetype="dashed", position = position_dodge(width=0.3))+
  geom_errorbar(aes(x = time, ymax = pam+se, ymin = pam-se), width = .2, position = position_dodge(width=0.3)) +
  #scale_color_manual(name = "Lineage-Sym",
  #                   values = c("#D53E4F", "#F46D43", "#3288BD", "#FEE08B"))+
  #values = c("#3f007d","#3f007d","#807dba","#807dba","#bcbddc","#bcbddc"))+
  #scale_fill_manual(name = "Lineage-Sym",
  #                  values = c("#D53E4F", "#F46D43", "#3288BD", "#FEE08B"))+
  #values = c("#3f007d","#3f007d","#807dba","#807dba","#bcbddc","#bcbddc"))+
  #scale_shape_manual(name = "Lineage-Sym",
  #                   values=c(21,24,21,24,21,24))+
  ylab("Fv/Fm")+
  scale_x_continuous(name = "Time Point", breaks = c(46,55,62,66,71,75,80)) +
  #ylim(-0.5,4) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  facet_wrap(~lineage)
pam_plot_lineage_symtype

ggsave(pam_plot_lineage_symtype, filename = "/Users/hannahaichelman/Documents/BU/TVE/PAM/pam_lineage_sym_2_lin_DIVs.pdf", width=8, height=4, units=c("in"), useDingbats=FALSE)

#SummarySE to format data for plotting - combine lineage and dominant sym type
phys_pam_2_lin_plots2 = phys_pam_2_lin_plots %>%
  drop_na(lineage, dominant_type) %>%
  unite(lin_sym, c(lineage,dominant_type), sep = "_", remove = FALSE) %>%
  mutate(lin_sym = as.factor(lin_sym))


pam_2_lineagediv_means <- summarySE(phys_pam_2_lin_plots2, measurevar="pam", groupvars=c("time","lin_sym","sitename"))
pam_2_lineagediv_means$time <- as.numeric(as.character(pam_2_lineagediv_means$time))

# plot, treatment x axis colored by site data figure
pam_plot_lineage_symtype <- ggplot(pam_2_lineagediv_means,aes(x = time, y = pam, color = lin_sym))+
  theme_bw()+
  annotate("rect", xmin = 46, xmax = 63, ymin = - Inf, ymax = Inf, fill = "red4", alpha = 0.15)+
  annotate("rect", xmin = 63, xmax = 80, ymin = - Inf, ymax = Inf, fill = "royalblue4", alpha = 0.15)+
  geom_point(size = 3, position = position_dodge(width=0.3))+
  geom_line(aes(group = lin_sym), size = 0.5, linetype="dashed", position = position_dodge(width=0.3))+
  geom_errorbar(aes(x = time, ymax = pam+se, ymin = pam-se), width = .2, position = position_dodge(width=0.3)) +
  #scale_color_manual(name = "Lineage-Sym",
  #                   values = c("#D53E4F", "#F46D43", "#3288BD", "#FEE08B"))+
  #values = c("#3f007d","#3f007d","#807dba","#807dba","#bcbddc","#bcbddc"))+
  #scale_fill_manual(name = "Lineage-Sym",
  #                  values = c("#D53E4F", "#F46D43", "#3288BD", "#FEE08B"))+
  #values = c("#3f007d","#3f007d","#807dba","#807dba","#bcbddc","#bcbddc"))+
  #scale_shape_manual(name = "Lineage-Sym",
  #                   values=c(21,24,21,24,21,24))+
  ylab("Fv/Fm")+
  scale_x_continuous(name = "Time Point", breaks = c(46,55,62,66,71,75,80)) +
  #ylim(-0.5,4) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  facet_wrap(~sitename)
pam_plot_lineage_symtype

ggsave(pam_plot_lineage_symtype, filename = "/Users/hannahaichelman/Documents/BU/TVE/PAM/pam_lineage_sym_2_lin_DIVs.pdf", width=8, height=4, units=c("in"), useDingbats=FALSE)

#SummarySE to format data for plotting
pam_means_site_all_lin <- summarySE(phys_pam_all_lin_plots, measurevar="pam", groupvars=c("sitename","time"))
pam_means_site_all_lin$time <- as.numeric(as.character(pam_means_site_all_lin$time))
pam_means_site_2_lin <- summarySE(phys_pam_2_lin_plots, measurevar="pam", groupvars=c("sitename","time"))
pam_means_site_2_lin$time <- as.numeric(as.character(pam_means_site_2_lin$time))

# plot, treatment x axis colored by site data figure
pam_plot_site <- ggplot(pam_means_site_2_lin,aes(x = time, y = pam, color = sitename, fill = sitename))+
  theme_bw()+
  annotate("rect", xmin = 51, xmax = 64, ymin = - Inf, ymax = Inf, fill = "red4", alpha = 0.15)+
  annotate("rect", xmin = 64, xmax = 80, ymin = - Inf, ymax = Inf, fill = "royalblue4", alpha = 0.15)+
  geom_point(size = 3, pch = 21, color = "black", position = position_dodge(width=0.3))+
  #geom_smooth(aes(group=sitename), position = position_dodge(width=0.3), size = 0.5, method='lm', se = FALSE)+
  geom_line(aes(group = sitename), size = 0.5, linetype="dashed", position = position_dodge(width=0.3))+
  geom_errorbar(aes(x = time, ymax = pam+se, ymin = pam-se), width = .2, position = position_dodge(width=0.3)) +
  scale_fill_manual(name = "Site",
                    labels = c("CI","PD","SP","BN","BS","CA"),
                    values = cols_site)+
  scale_color_manual(name = "Site",
                     labels = c("CI","PD","SP","BN","BS","CA"),
                     values = cols_site)+
  ylab("Fv/Fm")+
  scale_x_continuous(name = "Time Point", breaks = c(45,54,61,65,70,74,79)) +
  #ylim(-0.5,4) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))
pam_plot_site
ggsave(pam_plot_site, file="/Users/hannahaichelman/Documents/BU/TVE/PAM/PAM_site_2_lin.pdf", width=6, height=4, units=c("in"), useDingbats=FALSE)

# Plot separated by site and treatment to compare with ITS2 types
phys_pam_sites = phys_pam_lineage_plots %>%
  subset(sitename == "O-Cayo de Agua") %>%
  subset(treat == "High Var")

str(phys_pam_sites)
phys_pam_sites$frag = as.factor(phys_pam_sites$frag)
phys_pam_sites$time = as.numeric(as.character(phys_pam_sites$time))
phys_pam_sites$lineage = as.factor(phys_pam_sites$lineage)

# plot, time x axis, shape = symbiont type, color = lineage
pam_plot_lineage <- ggplot(phys_pam_sites,aes(x = time, y = pam, shape = frag, color = lineage))+
  theme_bw()+
  geom_point(size = 3, position = position_dodge(width=0.3))+
  geom_line(aes(color = lineage), size = 0.5, linetype="dashed", position = position_dodge(width=0.3))+
  #geom_errorbar(aes(x = time, ymax = pam+se, ymin = pam-se), width = .2, position = position_dodge(width=0.3)) +
  scale_color_manual(name = "Lineage",
                     labels = c("L1","L2"),
                     values = c("#fc8d59", "#99d594"))+
  scale_shape_manual(name = "Frag",
                     values=c(7,8,9,10, 11, 12, 13, 14))+
  ylab("Fv/Fm")+
  scale_x_continuous(name = "Time Point", breaks = c(46,55,62,66,71,75,80)) +
  #ylim(-0.5,4) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))
pam_plot_lineage

ggsave(pam_plot_lineage, filename = "/Users/hannahaichelman/Documents/BU/TVE/PAM/pam_CA_HighVar.pdf", width=6, height=4, units=c("in"), useDingbats=FALSE)


# include temperatures during PAM measurements
pam_temps <- read.csv("/Users/hannahaichelman/Documents/BU/TVE/PAM/pam_temp_timeline.csv")
str(pam_temps)

pam_temps$day <- as.numeric(as.character(pam_temps$day))

# plot
tempplot <-
  ggplot(pam_temps, aes(x = day, y = temp)) +
  theme_bw() +
  geom_point(size = 3) +
  geom_line() +
  ylim(27,33) +
  ylab("Temperature (°C)") +
  scale_x_continuous(name = "Time Point", breaks = c(46,55,62,66,71,75,80))
tempplot

quartz() # to be able to squish and save temp plot figure:
tempplot


#### Corallite Surface Area ####

# corallite surface area has a few extra samples from other initial phys dataframes, but
# it is still all taken from initial phys data

corrsa_phys <- read_excel("/Users/hannahaichelman/Documents/BU/TVE/Corallite_SA/Corallite_SA_Measurements.xlsx", sheet = "data")

head(corrsa_phys)

# remove test fragments and column of notes to look at data
corrsa_phys = corrsa_phys %>%
  dplyr::filter(sample.ID != "I2A6 (test)") %>%
  dplyr::filter(sample.ID != "I2A6 (test2)") %>%
  select(sample.ID, std.area, corallite.avg.area, corallite.avg.poly) %>%
  dplyr::rename(frag = sample.ID)

# re-name treatment
corrsa_phys$treat <- c("Initial")

# add identifying data
corrsa_phys$origsitecode <- substr(corrsa_phys$frag, 1, 2)

# add in site name with inshore/offshore indicator
corrsa_phys$sitename <- ifelse(corrsa_phys$origsitecode == 'I2', 'SP',
                               ifelse(corrsa_phys$origsitecode == 'I3', 'CI',
                                      ifelse(corrsa_phys$origsitecode == 'I4', 'PD',
                                             ifelse(corrsa_phys$origsitecode == 'O2', 'BS',
                                                    ifelse(corrsa_phys$origsitecode == 'O3', 'CA',
                                                           'BN')))))
corrsa_phys$sitename <- as.factor(corrsa_phys$sitename)

# make new nubbin IDs based on the new, more informative site codes
# add inshore/offshore designation
corrsa_phys$reef <- substr(corrsa_phys$frag,1,1)
corrsa_phys$reef <- ifelse(corrsa_phys$reef == 'O', 'Outer Reef', 'Inner Reef')
corrsa_phys$reef <- as.factor(corrsa_phys$reef)

corrsa_phys$genet <- substr(corrsa_phys$frag,3,3)

#create a new column of combined genotype and site for stats later
corrsa_phys = corrsa_phys %>%
  unite(gen_site, c(origsitecode,genet), sep = "", remove = FALSE) %>%
  mutate(gen_site = as.factor(gen_site)) %>%
  dplyr::filter(gen_site != "I4G") # clone with I4F, remove from dataset

# corallite.avg.area = Area of the corallite calculated by using area of ellipse formula (corallite.min/2*corallite.max/2*pi)
# corallite.avg.poly = Area of the corallite measured in ImageJ using the polygon tool

# look at relationship between the two ways that Laura measured corallite surface area
# correlate the two measurements

plot(corrsa_phys$corallite.avg.area, corrsa_phys$corallite.avg.poly, pch = 16)
abline(lm(corrsa_phys$corallite.avg.poly ~ corrsa_phys$corallite.avg.area), col = "red", lwd = 3) # add regression line
text(paste("Correlation:", round(cor(corrsa_phys$corallite.avg.area, corrsa_phys$corallite.avg.poly), 2)), x = 10000, y = 4500) # add Pearson correlation

# they are highly correlated, so will stick with the measurement taken using the polygon tool in ImageJ: corallite.avg.poly

# now need to convert between pixel # and mm^2 for the corallite.avg.poly measurement
# cross multiply to solve for unknown number of mm^2
corrsa_phys = corrsa_phys %>%
  mutate(corallite.avg.poly.mm2 = ((corallite.avg.poly*400)/std.area))

# merge with lineage info for later plotting
lineages = read.csv("/Users/hannahaichelman/Documents/BU/TVE/2bRAD/Analysis/tuftscustompipeline_denovo_nosyms/tve_lineages_noclones.csv")

corrsa_phys_all_lin <- left_join(corrsa_phys, lineages, by = "gen_site")
corrsa_phys_all_lin$lineage = as.factor(corrsa_phys_all_lin$lineage)

corrsa_phys_2_lin = corrsa_phys_all_lin %>%
  filter(is.na(lineage) | lineage!="L3") # want to keep NA values for lineage here since they still have other info, will remove na's for lineage specific plots

# Stats by lineage
str(corrsa_phys_all_lin)
corrsa_phys_all_lin$gen_site = as.factor(corrsa_phys_all_lin$gen_site)
m1 <- lmer(corallite.avg.poly.mm2 ~ lineage + (1|gen_site), data = corrsa_phys_all_lin, REML=TRUE)
summary(m1)
# Fixed effects:
#                 Estimate Std. Error      df t value Pr(>|t|)
#   (Intercept)   8.7252     0.3435 44.7154  25.404  < 2e-16 ***
#   lineageL2     5.1702     0.5718 44.8748   9.042 1.13e-11 ***
#   lineageL3     8.5199     1.1102 44.9221   7.674 1.04e-09 ***

anova(m1)
# Type III Analysis of Variance Table with Satterthwaite's method
#         Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)
# lineage 206.02  103.01     2 44.904  59.483 2.384e-13 ***

lsmeans(m1, pairwise~lineage, adjust="tukey")
# $contrasts
#             contrast estimate    SE   df t.ratio p.value
# L1 - L2     -5.17 0.574 44.2  -9.009  <.0001
# L1 - L3     -8.52 1.111 44.5  -7.666  <.0001
# L2 - L3     -3.35 1.151 44.6  -2.911  0.0152

#SummarySE to format data for plotting site name
corrsa_means_site_all_lin <- summarySE(corrsa_phys_all_lin, measurevar="corallite.avg.poly.mm2", groupvars=c("sitename"))
corrsa_means_site_2_lin <- summarySE(corrsa_phys_2_lin, measurevar="corallite.avg.poly.mm2", groupvars=c("sitename"))

# plot, treatment x axis colored by site data figure
corrsa_plot_site <- ggplot(corrsa_means_site_2_lin,aes(x = sitename, y = corallite.avg.poly.mm2, color = sitename, pch = sitename))+
  theme_bw()+
  geom_point(size = 3, position = position_dodge(width=0.3))+
  geom_errorbar(aes(x = sitename, ymax = corallite.avg.poly.mm2+se, ymin = corallite.avg.poly.mm2-se), width = .2, position = position_dodge(width=0.3)) +
  scale_color_manual(name = "Site",
                     labels = c("CI","PD","SP","BN","BS","CA"),
                     values = cols_site)+
  scale_shape_manual(name = "Site",
                     labels = c("CI","PD","SP","BN","BS","CA"),
                     values=c(19,19,19,17,17,17))+
  xlab("Site Name")+
  ylab(bquote("Corallite Area (mm" ^2~')'))+
  #ylim(3.75,6) +
  #geom_vline(xintercept = 1.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), legend.position = "none")
corrsa_plot_site

ggsave(corrsa_plot_site, filename = "/Users/hannahaichelman/Documents/BU/TVE/Corallite_SA/corrsa_site_2_lin.pdf", width=4, height=4, units=c("in"), useDingbats=FALSE)

# SummarySE for reef zone
corrsa_means_reef_all_lin <- summarySE(corrsa_phys_all_lin, measurevar="corallite.avg.poly.mm2", groupvars=c("reef"))
corrsa_means_reef_2_lin <- summarySE(corrsa_phys_2_lin, measurevar="corallite.avg.poly.mm2", groupvars=c("reef"))

# plot, reef zone x axis
corrsa_plot_reef <- ggplot(corrsa_means_reef_2_lin,aes(x = reef, y = corallite.avg.poly.mm2, color = reef, pch = reef))+
  theme_bw()+
  geom_point(size = 3, position = position_dodge(width=0.3))+
  geom_errorbar(aes(x = reef, ymax = corallite.avg.poly.mm2+se, ymin = corallite.avg.poly.mm2-se), width = .2, position = position_dodge(width=0.3)) +
  scale_color_manual(name = "Reef Zone",
                     labels = c("Inshore","Offshore"),
                     values = c("red4","royalblue4"))+
  scale_shape_manual(name = "Reef Zone",
                     labels = c("Inshore","Offshore"),
                     values=c(19,17))+
  xlab("Reef Zone")+
  ylab(bquote("Corallite Area (mm" ^2~')'))+
  #ylim(3.75,6) +
  #geom_vline(xintercept = 1.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  theme(legend.position = "none")
corrsa_plot_reef

ggsave(corrsa_plot_reef, filename = "/Users/hannahaichelman/Documents/BU/TVE/Corallite_SA/corrsa_reef_2_lin.pdf", width=4, height=4, units=c("in"), useDingbats=FALSE)


#SummarySE to format data for plotting with lineage
corrsa_phys_all_lin_nona = corrsa_phys_all_lin %>%
  drop_na(lineage)

corrsa_phys_2_lin_nona = corrsa_phys_2_lin %>%
  drop_na(lineage)

corrsa_means_all_lin <- summarySE(corrsa_phys_all_lin_nona, measurevar="corallite.avg.poly.mm2", groupvars=c("lineage"))
corrsa_means_2_lin <- summarySE(corrsa_phys_2_lin_nona, measurevar="corallite.avg.poly.mm2", groupvars=c("lineage"))

# plot, lineage x axis
corrsa_plot_lineage <- ggplot(corrsa_means_all_lin,aes(x = lineage, y = corallite.avg.poly.mm2, color = lineage, fill = lineage))+
  theme_bw()+
  geom_errorbar(aes(x = lineage, ymax = corallite.avg.poly.mm2+se, ymin = corallite.avg.poly.mm2-se, color = lineage), width = .2, position = position_dodge(width=0.3)) +
  geom_point(size = 3, position = position_dodge(width=0.3), shape = 21, color = "black")+
  scale_color_manual(name = "Lineage",
                     #breaks = c("L1","L2"),
                     breaks = c("L1","L2","L3"),
                     values = cols_lineage)+
  scale_fill_manual(name = "Lineage",
                    #breaks = c("L1","L2"),
                    breaks = c("L1","L2","L3"),
                    values=cols_lineage)+
  xlab("Lineage")+
  ylab(bquote("Corallite Area (mm" ^2~')'))+
  scale_y_continuous(limits = c(6,19), breaks = seq(6,19, by = 3))+
  #geom_vline(xintercept = 1.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  theme(legend.position = "none")
corrsa_plot_lineage

ggsave(corrsa_plot_lineage, filename = "/Users/hannahaichelman/Documents/BU/TVE/Corallite_SA/corrsa_lineage_all.pdf", width=3, height=3, units=c("in"), useDingbats=FALSE)

# plot just CI corals to see if lineage difference holds
corrsa_phys_lineage_CI = corrsa_phys_2_lin_nona %>%
  subset(sitename == "CI")

# stats for CI only
str(corrsa_phys_lineage_CI)
m1 <- lm(corallite.avg.poly.mm2 ~ lineage, data = corrsa_phys_lineage_CI, REML=TRUE)
summary(m1)
anova(m1)
lsmeans(m1, pairwise~lineage, adjust="tukey")

corrsa_means_lineage_CI <- summarySE(corrsa_phys_lineage_CI, measurevar="corallite.avg.poly.mm2", groupvars=c("lineage"))

corrsa_plot_lineage_CI <- ggplot(corrsa_means_lineage_CI,aes(x = lineage, y = corallite.avg.poly.mm2, color = lineage, fill = lineage))+
  theme_bw()+
  geom_errorbar(aes(x = lineage, ymax = corallite.avg.poly.mm2+se, ymin = corallite.avg.poly.mm2-se), width = .2, position = position_dodge(width=0.3)) +
  geom_point(size = 3, position = position_dodge(width=0.3), shape = 21, color = "black")+
  scale_color_manual(name = "Lineage",
                     breaks = c("L1","L2"),
                     values = cols_lineage)+
  scale_fill_manual(name = "Lineage",
                    breaks = c("L1","L2"),
                    values=cols_lineage)+
  xlab("Lineage")+
  ylab(bquote("Corallite Area (mm" ^2~')'))+
  scale_y_continuous(limits = c(6,19), breaks = seq(6,19, by = 3))+
  #geom_vline(xintercept = 1.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  theme(legend.position = "none")
corrsa_plot_lineage_CI

ggsave(corrsa_plot_lineage_CI, filename = "/Users/hannahaichelman/Documents/BU/TVE/Corallite_SA/corrsa_lineage_CIonly.pdf", width=3, height=3, units=c("in"), useDingbats=FALSE)


#### Prep data for PCA ####

# separating data frames by metric and by time point to consider
calc_pca_end = calc_phys_all_lin %>%
  select(frag,treat,sitename,lineage,T3_T0_rgr) # want to consider growth over the whole experiment

hcarb_pca_t0 = hcarb_phys_all_lin %>%
  dplyr::filter(treat=="Initial") %>%
  select(frag,treat,sitename,lineage,hcarb_mgcm2)

hcarb_pca_end = hcarb_phys_all_lin %>%
  dplyr::filter(treat!="Initial") %>%
  select(frag,treat,sitename,lineage,hcarb_mgcm2)

scarb_pca_t0 = scarb_phys_all_lin %>%
  dplyr::filter(treat=="Initial") %>%
  select(frag,treat,sitename,lineage,scarb_mgcm2)

scarb_pca_end = scarb_phys_all_lin %>%
  dplyr::filter(treat!="Initial") %>%
  select(frag,treat,sitename,lineage,scarb_mgcm2)

hprot_pca_t0 = hprot_phys_all_lin %>%
  dplyr::filter(treat=="Initial") %>%
  select(frag,treat,sitename,lineage,prot_mgcm2)

hprot_pca_end = hprot_phys_all_lin %>%
  dplyr::filter(treat!="Initial") %>%
  select(frag,treat,sitename,lineage,prot_mgcm2)

tissthick_pca_t0 = tiss_phys_all_lin %>%
  select(frag,treat,sitename,lineage,avgtiss)

chla_pca_t0 = chl_phys_all_lin %>%
  dplyr::filter(treat=="Initial") %>%
  select(frag,treat,sitename,lineage,chlA)

chla_pca_end = chl_phys_all_lin %>%
  dplyr::filter(treat!="Initial") %>%
  select(frag,treat,sitename,lineage,chlA)

corrsa_pca_t0 = corrsa_phys_all_lin %>%
  select(frag,treat,sitename,lineage,corallite.avg.poly.mm2)

# Initial PAM measurements were not taken, so no t0 PAM
pam_pca_end = pam_phys_all_lin %>%
  dplyr::filter(time=="79") %>% # corresponds to the closest timepoint from when phenomic metrics were sampled
  select(frag,treat,sitename,lineage,pam)

sym_pca_t0 = sym_phys_all_lin %>%
  dplyr::filter(treat=="Initial") %>%
  select(frag,treat,sitename,lineage,sym_cm2)

sym_pca_end = sym_phys_all_lin %>%
  dplyr::filter(treat!="Initial") %>%
  select(frag,treat,sitename,lineage,sym_cm2)

# Make T0 combined data frame
t0_full = join_all(list(hcarb_pca_t0,scarb_pca_t0,hprot_pca_t0,sym_pca_t0,tissthick_pca_t0,chla_pca_t0,corrsa_pca_t0), by = c("frag","treat","sitename","lineage"), type = "full")

str(t0_full)
#write.csv(t0_full, "/Users/hannahaichelman/Documents/BU/TVE/PCAs/t0_full.csv",row.names=FALSE)

t0_full_log = t0_full %>%
  dplyr::filter(complete.cases(.)) %>% #drop any row that has an NA for any time point
  mutate_if(is.numeric, log)

str(t0_full_log)

t0_full_log$sitename <- as.factor(t0_full_log$sitename)
t0_full_log$treat <- as.factor(t0_full_log$treat)
#write.csv(t0_full_log, "/Users/hannahaichelman/Documents/BU/TVE/PCAs/t0_full_log.csv", row.names=FALSE)

# Make end of variability combined data frame
end_full = join_all(list(calc_pca_end,hcarb_pca_end,scarb_pca_end,hprot_pca_end,pam_pca_end,sym_pca_end,chla_pca_end), by = c("frag","treat","sitename","lineage"), type = "full")

str(end_full)
#write.csv(end_full, "/Users/hannahaichelman/Documents/BU/TVE/PCAs/end_full.csv", row.names = FALSE)

end_full_log = end_full %>%
  dplyr::filter(complete.cases(.)) %>% #drop any row that has an NA for any time point
  mutate(T3_T0_rgr_2 = T3_T0_rgr + 2) %>%
  select(-T3_T0_rgr) %>% #get rid of the column with negative growth values before log transforming
  mutate_if(is.numeric, log) # log transform the numeric columns

str(end_full_log)

end_full_log$sitename <- as.factor(end_full_log$sitename)
end_full_log$treat <- as.factor(end_full_log$treat)
end_full_log$lineage <- as.factor(end_full_log$lineage)
#write.csv(end_full_log, "/Users/hannahaichelman/Documents/BU/TVE/PCAs/end_full_log.csv", row.names = FALSE)

##### Make our PCAs #####
# once you have made the output files in the section above you can skip right to here to finalize filtering and make pca's

library(ggpubr)
library(ggfortify)
library(ggplot2)
library(cluster)
library(FactoMineR) # lots of options for pca visuals and summary stats
library(factoextra)
library(corrplot)
library(dplyr)
#library(ggbiplot)
library(cowplot)
library(vegan)

# try exporting data files as csv and see what we get
t0_pca = read.csv("/Users/hannahaichelman/Documents/BU/TVE/PCAs/t0_full_log.csv")
end_pca = read.csv("/Users/hannahaichelman/Documents/BU/TVE/PCAs/end_full_log.csv")

# format t0 physiology data:
str(t0_pca)

t0_pca$sitename <- as.factor(t0_pca$sitename)
t0_pca$treat <- as.factor(t0_pca$treat)
t0_pca$lineage <- as.factor(t0_pca$lineage)
# all of the samples that have NA's for lineage are ones that didn't prep so we don't have lineage data for those genets

colnames(t0_pca)[colnames(t0_pca)=="hcarb_mgcm2"] <-"hcarb"
colnames(t0_pca)[colnames(t0_pca)=="scarb_mgcm2"] <-"scarb"
colnames(t0_pca)[colnames(t0_pca)=="prot_mgcm2"] <-"prot"
colnames(t0_pca)[colnames(t0_pca)=="sym_cm2"] <-"syms"
colnames(t0_pca)[colnames(t0_pca)=="avgtiss"] <-"tiss"
colnames(t0_pca)[colnames(t0_pca)=="corallite.avg.poly.mm2"] <-"corr_sa"

t0_pca_all_lin = t0_pca %>%
  drop_na(lineage)

t0_pca_2_lin = t0_pca_all_lin %>%
  dplyr::filter(is.na(lineage) | lineage!="L3") # want to keep NA values for lineage here since they still have other info, will remove na's for lineage specific plots

t0_pca_lineage_CI = t0_pca_2_lin %>%
  dplyr::filter(sitename == "CI")


# now end of variability pca data:
str(end_pca)

end_pca$sitename <- as.factor(end_pca$sitename)
end_pca$treat <- as.factor(end_pca$treat)
end_pca$lineage <- as.factor(end_pca$lineage)

# re-level and re-name treatment
end_pca$treat <- factor(end_pca$treat, levels = c("Control","Low Var","Mod Var","High Var"))

# re-name columns for legibility on plots
colnames(end_pca)[colnames(end_pca)=="hcarb_mgcm2"] <-"hcarb"
colnames(end_pca)[colnames(end_pca)=="scarb_mgcm2"] <-"scarb"
colnames(end_pca)[colnames(end_pca)=="prot_mgcm2"] <-"prot"
colnames(end_pca)[colnames(end_pca)=="pam"] <-"fv_fm"
colnames(end_pca)[colnames(end_pca)=="sym_cm2"] <-"syms"
colnames(end_pca)[colnames(end_pca)=="T3_T0_rgr_2"] <-"growth"

# add in dominant symbiont type dataframe
its2_types = read.csv("/Users/hannahaichelman/Documents/BU/TVE/16S_ITS2/ITS_PreStress_Timepoint/ITS2.dominanttype.csv") %>%
  select(frag, dominant_type)

its2_types$dominant_type = as.factor(its2_types$dominant_type)

end_pca_sym <- left_join(end_pca, its2_types, by = "frag")
str(end_pca_sym)

end_pca_all_lin = end_pca_sym %>%
  drop_na(lineage) %>%
  drop_na(dominant_type) %>%
#  unite(lin_sym, c(lineage,dominant_type), sep = "_", remove = FALSE) %>% # make new lineage_dominant type combined factor
#  mutate(lin_sym = as.factor(lin_sym)) %>%
  select(frag, treat, sitename, lineage, dominant_type, hcarb, scarb, prot, fv_fm, syms, chlA, growth)

end_pca_2_lin = end_pca_all_lin %>%
  dplyr::filter(lineage!="L3")

end_pca_lineage_CI = end_pca_2_lin %>%
  dplyr::filter(sitename == "CI")

###
# this part of the code uses the PCA function and stats from the FactoMineR package
# T0 both lineage df's
facto_t0_all_lin <- PCA(t0_pca_all_lin[,5:11], scale.unit = TRUE, ncp = 10, graph = TRUE)
facto_t0_2_lin <- PCA(t0_pca_2_lin[,5:11], scale.unit = TRUE, ncp = 10, graph = TRUE)

facto_end_all_lin <- PCA(end_pca_all_lin[,6:12], scale.unit = TRUE, ncp = 10, graph = TRUE) # no sym type info
facto_end_2_lin <- PCA(end_pca_2_lin[,6:12], scale.unit = TRUE, ncp = 10, graph = TRUE)

#facto_t0_CI <- PCA(t0_pca_lineage_CI[,5:11], scale.unit = TRUE, ncp = 10, graph = TRUE)
#facto_end_CI <- PCA(end_pca_lineage_CI[,5:11], scale.unit = TRUE, ncp = 10, graph = TRUE)

# correlations variables - dimensions
facto_end_2_lin$var$cor
# mean of the variables
facto_end_2_lin$call$centre
# standard error of variables
facto_end_2_lin$call$ecart.type
eig.val1 <- get_eigenvalue(facto_end_2_lin)
eig.val1
fviz_eig(facto_end_2_lin, addlabels = TRUE) # get percentage of explained variances per dimension
var1 <- get_pca_ind(facto_end_2_lin)
var1
var2 <- get_pca_var(facto_end_2_lin)
var2

# Contributions to the principal components
head(var2$contrib, 10)

# make a corrplot on contribution to PC
corrplot(var2$contrib, is.corr=FALSE)

## red line = expected average contribution... "In this case, the expected average contribution (cutoff) is calculated as follow: As mentioned above, if the contributions of the 10 variables were uniform, the expected average contribution on a given PC would be 1/10 = 10%. The expected average contribution of a variable for PC1 and PC2 is : [(10* Eig1) + (10 * Eig2)]/(Eig1 + Eig2)"
# bar graphs of contributions for all of the dimensions... for variables
c <- fviz_contrib(facto_end_2_lin, choice = "var", axes = 1, top = 10)
save_plot("pc1_contribution_of_var.png",c, base_aspect_ratio = 1.6)
d <- fviz_contrib(facto_end_2_lin, choice = "var", axes = 2, top = 10)
save_plot("pc2_contribution_of_var.png",d, base_aspect_ratio = 1.6)
fviz_contrib(facto_end_2_lin, choice = "var", axes = 3, top = 10)
fviz_pca_var(facto_end_2_lin, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")
)

### describe dimensions
dim_descript <- dimdesc(facto_t0_2_lin, axes = c(1,2,3,4), proba = 0.05)
# dim1
dim_descript$Dim.1
# dim2
dim_descript$Dim.2

### make the biplot
# this fuses the facto package and ggplot features (uses base facto biplot and overlays with geom_point so can have shape demark one factor and color another factor)
fviz_pca_biplot(facto_t0_all_lin)
fviz_pca_biplot(facto_t0_2_lin)

fviz_pca_biplot(facto_end_all_lin)
fviz_pca_biplot(facto_end_2_lin)

fviz_pca_biplot(facto_end_CI)
fviz_pca_biplot(facto_t0_CI)

# these are the color schemes, defined in first section of code
cols_treat
cols_site
cols_lineage
its2_cols_greens

# PCA for T0 Physiology - color = sitename
fviz_pca_biplot(facto_t0_2_lin)
pca_t0_site <- fviz_pca_biplot(facto_t0_2_lin,
                               label = "var",
                               col.var = "black", labelsize = 4,
                               alpha.ind = 0) + # makes individs transparent so they can be overwritten by geom_point()
  theme_bw()+
  geom_point(aes(colour=t0_pca_2_lin$sitename), size = 2, stroke = 1) +
  scale_color_manual(values = cols_site,
                     #breaks=c("CI","I-Punta Donato","I-STRI Point","O-Bastimentos N","O-Bastimentos S","O-Cayo de Agua"),
                     labels=c("CI","PD","SP","BN","BS","CA"),
                     name = "Site") +
  #  stat_ellipse(geom = "polygon", type = "t", alpha = 0.2,
  #               aes(fill= t0_pca_lineage$sitename), show.legend = FALSE) + scale_fill_manual(values=cols_site) + # ellipses assumes multivariate distribution using default confidence level (0.95)
  stat_ellipse(aes(color=t0_pca_2_lin$sitename), type = "t", lwd = 1)+
  labs(x = "PC1 (51.1% explained variance)",
       y = "PC2 (13.3% explained variance)") +
  theme(plot.title = element_blank(),
        axis.title = element_text(face = "bold"),
        legend.title.align =  0.5, legend.text.align = 0,
        legend.title = element_text(face = "bold"))
pca_t0_site

ggsave(pca_t0_site, filename = "/Users/hannahaichelman/Documents/BU/TVE/PCAs/pca_t0_site_2_lin.pdf", width=6, height=5, units=c("in"), useDingbats=FALSE)


# PCA for T0 Physiology - color = lineage, shape = sitename
pca_t0_lineage <- fviz_pca_biplot(facto_t0_2_lin,
                                  label = "var",
                                  col.var = "black", labelsize = 4,
                                  alpha.ind = 0) + # makes individs transparent so they can be overwritten by geom_point()
  theme_bw()+
  geom_point(aes(colour=t0_pca_2_lin$lineage, shape = t0_pca_2_lin$sitename), size = 2, stroke = 1) +
  scale_color_manual(values = cols_lineage, breaks=c("L1","L2"), name = "Lineage") +
  scale_shape_manual(values = c(15,16,17,22,21,24),
                     breaks=c("BN", "BS", "CA", "CI", "PD", "SP"),
                     labels=c("BN", "BS", "CA", "CI", "PD", "SP"),
                     name = "Site") +
  #  stat_ellipse(geom = "polygon", type = "t", alpha = 0.4,
  #               aes(fill= t0_pca_lineage$lineage), show.legend = FALSE) + scale_fill_manual(values=cols_lineage) + # ellipses assumes multivariate distribution using default confidence level (0.95)
  stat_ellipse(aes(color=t0_pca_2_lin$lineage), type = "t", lwd = 1)+
  labs(x = "PC1 (51.1% explained variance)",
       y = "PC2 (13.3% explained variance)") +
  theme(plot.title = element_blank(),
        axis.title = element_text(face = "bold"),
        legend.title.align =  0.5, legend.text.align = 0,
        legend.title = element_text(face = "bold"))
pca_t0_lineage

ggsave(pca_t0_lineage, filename = "/Users/hannahaichelman/Documents/BU/TVE/PCAs/pca_t0_lineage_2.pdf", width=6, height=5, units=c("in"), useDingbats=FALSE)

# PCA for T0 Physiology CI ONLY- color = lineage
pca_t0_lineage_CI <- fviz_pca_biplot(facto_t0_CI,
                                     label = "var",
                                     col.var = "black", labelsize = 4,
                                     alpha.ind = 0) + # makes individs transparent so they can be overwritten by geom_point()
  geom_point(aes(colour=t0_pca_lineage_CI$lineage), size = 2, stroke = 1) +
  scale_color_manual(values = cols_lineage, breaks=c("L1","L2"), name = "Lineage", labels = c("SSID-A","SSID-B")) +
  #  stat_ellipse(geom = "polygon", type = "t", alpha = 0.4,
  #               aes(fill= t0_pca_lineage_CI$lineage), show.legend = FALSE) + scale_fill_manual(values=cols_lineage) + # ellipses assumes multivariate distribution using default confidence level (0.95)
  stat_ellipse(aes(color=t0_pca_lineage_CI$lineage), type = "t", lwd = 1)+
  labs(x = "PC1 (57.98% explained variance)",
       y = "PC2 (19.99% explained variance)") +
  theme(plot.title = element_blank(),
        axis.title = element_text(face = "bold"),
        legend.title.align =  0.5, legend.text.align = 0,
        legend.title = element_text(face = "bold"))
pca_t0_lineage_CI

ggsave(pca_t0_lineage_CI, filename = "/Users/hannahaichelman/Documents/BU/TVE/PCAs/pca_t0_lineage_CIonly.pdf", width=5, height=5, units=c("in"), useDingbats=FALSE)


# PCA for End of Variability Physiology - color = site, shape = site
fviz_pca_biplot(facto_end_2_lin)
pca_end_site <- fviz_pca_biplot(facto_end_2_lin,
                                label = "var",
                                col.var = "black", labelsize = 4,
                                alpha.ind = 0) + # makes individs transparent so they can be overwritten by geom_point()
  theme_bw()+
  geom_point(aes(colour=end_pca_2_lin$sitename, shape =end_pca_2_lin$sitename), size = 2, stroke = 1) +
  scale_color_manual(values = cols_site,
                     breaks=c("BN", "BS", "CA", "CI", "PD", "SP"),
                     labels=c("BN", "BS", "CA", "CI", "PD", "SP"),
                     name = "Site") +
  scale_shape_manual(values = c(15,16,17,22,21,24),
                     breaks = c("BN", "BS", "CA", "CI", "PD", "SP"),
                     labels=c("BN", "BS", "CA", "CI", "PD", "SP"),
                     name = "Site") +
  stat_ellipse(aes(color=end_pca_2_lin$sitename), type = "t", lwd = 1, show.legend = FALSE)+
  labs(x = "PC1 (45.4% explained variance)",
       y = "PC2 (14.9% explained variance)") +
  theme(plot.title = element_blank(),
        axis.title = element_text(face = "bold"),
        legend.title.align =  0.5, legend.text.align = 0,
        legend.title = element_text(face = "bold"))
pca_end_site

ggsave(pca_end_site, filename = "/Users/hannahaichelman/Documents/BU/TVE/PCAs/pca_end_site_2.pdf", width=5, height=4, units=c("in"), useDingbats=FALSE)

# PCA for End of Variability Physiology - color = treatment, shape = site
fviz_pca_biplot(facto_end_2_lin)
pca_end_treat <- fviz_pca_biplot(facto_end_2_lin,
                                 label = "var",
                                 col.var = "black", labelsize = 4,
                                 alpha.ind = 0) + # makes individs transparent so they can be overwritten by geom_point()
  theme_bw()+
  geom_point(aes(colour=end_pca_2_lin$treat, shape = end_pca_2_lin$sitename), size = 2, stroke = 1) +
  scale_color_manual(values = cols_treat,
                     #breaks=c("Control", "Low Var", "Mod Var", "High Var"),
                     #labels=c("Control", "Low Var", "Mod Var", "High Var"),
                     name = "Treatment") +
  scale_shape_manual(values = c(15,16,17,22,21,24),
                     breaks = c("BN", "BS", "CA", "CI", "PD", "SP"),
                     labels=c("BN", "BS", "CA", "CI", "PD", "SP"),
                     name = "Site") +
  #  stat_ellipse(geom = "polygon", type = "t", alpha = 0.3,
  #               aes(fill= end_pca_lineage$treat), show.legend = FALSE) + scale_fill_manual(values=cols_treat) + # ellipses assumes multivariate distribution using default confidence level (0.95)
  stat_ellipse(aes(color=end_pca_2_lin$treat), type = "t", lwd = 1)+
  labs(x = "PC1 (45.4% explained variance)",
       y = "PC2 (14.9% explained variance)") +
  theme(plot.title = element_blank(),
        axis.title = element_text(face = "bold"),
        legend.title.align =  0.5, legend.text.align = 0,
        legend.title = element_text(face = "bold"))
pca_end_treat

ggsave(pca_end_treat, filename = "/Users/hannahaichelman/Documents/BU/TVE/PCAs/pca_end_treat_2_lin.pdf", width=5, height=4, units=c("in"), useDingbats=FALSE)

# PCA for End Variability Physiology - color = lineage, shape = sitename
fviz_pca_biplot(facto_end_2_lin)
pca_end_lineage <- fviz_pca_biplot(facto_end_2_lin,
                                   label = "var",
                                   col.var = "black", labelsize = 4,
                                   alpha.ind = 0) + # makes individs transparent so they can be overwritten by geom_point()
  theme_bw()+
  geom_point(aes(colour=end_pca_2_lin$lineage, shape = end_pca_2_lin$sitename), size = 2, stroke = 1) +
  scale_color_manual(values = cols_lineage, breaks=c("L1","L2"), name = "Lineage") +
  scale_shape_manual(values = c(15,16,17,22,21,24),
                     breaks = c("BN", "BS", "CA", "CI", "PD", "SP"),
                     labels=c("BN", "BS", "CA", "CI", "PD", "SP"),
                     name = "Site") +
  #  stat_ellipse(geom = "polygon", type = "t", alpha = 0.4,
  #               aes(fill= end_pca_lineage$lineage), show.legend = FALSE) + scale_fill_manual(values=cols_lineage) + # ellipses assumes multivariate distribution using default confidence level (0.95)
  stat_ellipse(aes(color= end_pca_2_lin$lineage), type = "t", lwd = 1)+
  labs(x = "PC1 (45.4% explained variance)",
       y = "PC2 (14.9% explained variance)") +
  theme(plot.title = element_blank(),
        axis.title = element_text(face = "bold"),
        legend.title.align =  0.5, legend.text.align = 0,
        legend.title = element_text(face = "bold"))
pca_end_lineage

ggsave(pca_end_lineage, filename = "/Users/hannahaichelman/Documents/BU/TVE/PCAs/pca_end_lineage_2.pdf", width=5, height=4, units=c("in"), useDingbats=FALSE)

# PCA for End Variability Physiology - color = lineage+dominant sym type, shape = sitename
its2_cols_greens = c("C1" = "#edf8e9", "C3af" = "#238b45","C3" = "#a1d99b","D1" = "#00441b")

fviz_pca_biplot(facto_end_2_lin)
pca_end_lineagesym <- fviz_pca_biplot(facto_end_2_lin,
                                      label = "var",
                                      col.var = "black", labelsize = 4,
                                      alpha.ind = 0) + # makes individs transparent so they can be overwritten by geom_point()
  theme_bw()+
  geom_point(aes(fill=end_pca_2_lin$dominant_type, shape = end_pca_2_lin$lineage), size = 2, stroke = 1) +
  scale_fill_manual(values = its2_cols_greens,
                    name = "Majority ITS2") +
  scale_shape_manual(values = c(21,24),
                     name = "Lineage") +
  stat_ellipse(aes(color= end_pca_2_lin$dominant_type), type = "t", lwd = 1, show.legend = FALSE)+
  scale_color_manual(values = its2_cols_greens,
                     name = "Majority ITS2") +
  labs(x = "PC1 (45.4% explained variance)",
       y = "PC2 (14.9% explained variance)") +
  theme(plot.title = element_blank(),
        axis.title = element_text(face = "bold"),
        legend.title.align =  0.5, legend.text.align = 0,
        legend.title = element_text(face = "bold"))+
  guides(color = guide_legend(override.aes = list(color = its2_cols_greens)))
pca_end_lineagesym

ggsave(pca_end_lineagesym, filename = "/Users/hannahaichelman/Documents/BU/TVE/PCAs/pca_end_dominantsym_2lin.pdf", width=5, height=4, units=c("in"), useDingbats=FALSE)

# PCA for End Variability Physiology CI ONLY- color = lineage
pca_end_lineage_CI <- fviz_pca_biplot(facto_end_CI,
                                      label = "var",
                                      col.var = "black", labelsize = 4,
                                      alpha.ind = 0) + # makes individs transparent so they can be overwritten by geom_point()
  geom_point(aes(colour=end_pca_lineage_CI$lineage), size = 2, stroke = 1) +
  scale_color_manual(values = cols_lineage, breaks=c("L1","L2"), name = "Lineage", labels = c("SSID-A","SSID-B")) +
  #  stat_ellipse(geom = "polygon", type = "t", alpha = 0.4,
  #               aes(fill= end_pca_lineage_CI$lineage), show.legend = FALSE) + scale_fill_manual(values=cols_lineage) + # ellipses assumes multivariate distribution using default confidence level (0.95)
  stat_ellipse(aes(color= end_pca_lineage_CI$lineage), type = "t", lwd = 1)+
  labs(x = "PC1 (36.3% explained variance)",
       y = "PC2 (21.1% explained variance)") +
  theme(plot.title = element_blank(),
        axis.title = element_text(face = "bold"),
        legend.title.align =  0.5, legend.text.align = 0,
        legend.title = element_text(face = "bold"))
pca_end_lineage_CI

ggsave(pca_end_lineage_CI, filename = "/Users/hannahaichelman/Documents/BU/TVE/PCAs/pca_end_lineage_CIonly.pdf", width=5, height=5, units=c("in"), useDingbats=FALSE)


# try faceting by sitename using ggplot
t0_df <- prcomp(t0_pca[,5:9], scale. = TRUE)
end_df <- prcomp(end_pca_2_lin[,9:15], scale. = TRUE)

facto_end_all_lin <- PCA(end_pca_all_lin[,9:15], scale.unit = TRUE, ncp = 10, graph = TRUE) # no sym type info
facto_end_2_lin <- PCA(end_pca_2_lin[,9:15], scale.unit = TRUE, ncp = 10, graph = TRUE)

# T0 autoplot faceted by sitename
autoplot(t0_df, data = t0_pca,
         colour = "sitename",
         loadings = TRUE, loadings.colour = "black", loadings.label = TRUE, loadings.label.size = 3, loadings.label.colour = "black",
         frame = TRUE, frame.type = "norm", frame.color = "sitename") +
  scale_color_manual(name = "Site",
                     labels = c("I-Cristobal","I-Punta Donato","I-STRI Point","O-Bastimentos N","O-Bastimentos S","O-Cayo de Agua"),
                     values = palsite)+
  scale_fill_manual(name = "Site",
                    labels = c("I-Cristobal","I-Punta Donato","I-STRI Point","O-Bastimentos N","O-Bastimentos S","O-Cayo de Agua"),
                    values = palsite)+
  facet_wrap(~ sitename, ncol = 2) +
  ggtitle("T0 Physiology")

# End of variability autoplot faceted by sitename and colored by treatment
autoplot(end_df, data = end_pca_2_lin,
         colour = "dominant_type",
         loadings = TRUE, loadings.colour = "black", loadings.label = TRUE, loadings.label.size = 3, loadings.label.colour = "black",
         frame = TRUE, frame.type = "norm", frame.color = "lineage") +
  scale_color_manual(name = "Dominant Type",
                     #labels = c("Control", "Low Var", "Mod Var", "High Var"),
                     values = c("#b8e186","#d73027","#fdae61","#74add1"))+
  scale_fill_manual(name = "Dominant Type",
                    #labels = c("Control", "Low Var", "Mod Var", "High Var"),
                    values = c("#b8e186","#d73027","#fdae61","#74add1"))+
  facet_wrap(~ lineage, ncol = 2) +
  ggtitle("End of Variability Physiology") +
  theme_bw()


##### PCA Adonis Tests #####
#Use an Adonis test to get significance of factors on holobiont physiology
library(vegan)
library(MCMC.OTU)
library(MicEco)

t0_full = read.csv("/Users/hannahaichelman/Documents/BU/TVE/PCAs/t0_full.csv")
end_full = read.csv("/Users/hannahaichelman/Documents/BU/TVE/PCAs/end_full.csv")

str(t0_full)

t0_full$sitename <- as.factor(t0_full$sitename)
t0_full$treat <- as.factor(t0_full$treat)
t0_full$lineage <- as.factor(t0_full$lineage)

str(end_full)

end_full$sitename <- as.factor(end_full$sitename)
end_full$treat <- as.factor(end_full$treat)
end_full$lineage <- as.factor(end_full$lineage)

# add in dominant symbiont type dataframe
its2_types = read.csv("/Users/hannahaichelman/Documents/BU/TVE/16S_ITS2/ITS_PreStress_Timepoint/ITS2.dominanttype.csv") %>%
  select(frag, dominant_type) %>%
  mutate(dominant_type = as.factor(dominant_type))

end_full_its2 <- left_join(end_full, its2_types, by = "frag")
str(end_full_its2)

# Dont want to use log transformed data here
# but need to remove NAs
t0_full_adonis = t0_full %>%
  dplyr::filter(complete.cases(.)) %>% #drop any row that has an NA for any time point
  dplyr::filter(lineage != "L3") # remove L3 individuals since these aren't in our PCAs

end_full_adonis = end_full_its2 %>%
  dplyr::filter(complete.cases(.)) %>% #drop any row that has an NA for any time point
  mutate(T3_T0_rgr_2 = T3_T0_rgr + 2) %>%
  select(-T3_T0_rgr) %>% #get rid of the column with negative growth values before log transforming, just making sure to use the same data as we do in the PCAs
  dplyr::filter(lineage != "L3") # remove L3 individuals since these aren't in our PCAs

end_full_adonis <- end_full_adonis[, c(1,2,3,4,11,5,6,7,8,9,10,12)]

end_full_adonis_L1 = end_full_adonis %>%
  filter(lineage == "L1")

end_full_adonis_L2 = end_full_adonis %>%
  filter(lineage == "L2")

# Change dataframe here based on the comparison you are interested in
nl=startedLog(data=end_full_adonis,count.columns=6:12, logstart=1)
#nl=startedLog(data=t0_full_adonis,count.columns=5:11, logstart=1)

goods.dist=vegdist(nl, method="bray", na.rm = TRUE)
goods.pcoa=pcoa(goods.dist)

# PCA:
pcp=prcomp(nl, retx=TRUE, center=TRUE)
scores=goods.pcoa$vectors
summary(goods.pcoa)
conditions=end_full_adonis[, c("frag","treat","sitename","lineage","dominant_type")] #make sure to change dataframe here
#conditions=t0_full_adonis[, c("frag","treat","sitename","lineage")] #make sure to change dataframe here

# PERMANOVA
head(scores)
head(conditions)

t0_model = adonis(scores~lineage+sitename, data=conditions, method="euclidean", permutations = 10000)
end_model = adonis(scores~lineage+dominant_type+treat+sitename, data=conditions, method="euclidean", permutations = 10000)

adonis_OmegaSq(t0_model, partial = TRUE)
adonis_OmegaSq(end_model, partial = TRUE)

# T0:
#           Df SumsOfSqs   MeanSqs F.Model      R2 parOmegaSq    Pr(>F)
# lineage    1  0.031258 0.0312576  20.837 0.33256    0.32080 9.999e-05 ***
# sitename   5  0.010230 0.0020461   1.364 0.10884    0.04153    0.1923
# Residuals 35  0.052503 0.0015001         0.55860
# Total     41  0.093991                   1.00000

# End:
# Terms added sequentially (first to last)
#
#               Df SumsOfSqs  MeanSqs F.Model      R2 parOmegaSq Pr(>F)
# lineage        1   0.04145 0.041446 14.5635 0.11268    0.12050 0.0002 ***
# dominant_type  3   0.00253 0.000845  0.2969 0.00689   -0.02177 0.9736
# treat          3   0.04381 0.014605  5.1319 0.11912    0.11128 0.0002 ***
# sitename       5   0.03527 0.007055  2.4789 0.09590    0.06950 0.0115 *
# Residuals     86   0.24475 0.002846         0.66541
# Total         98   0.36781                  1.00000

# End - Prop D1 NOT Included - justification for not including any interaction terms in the output below
# adonis2(formula = scores ~ lineage * dominant_type * treat * sitename, data = conditions, permutations = 10000, method = "euclidean")
#                                 Df SumsOfSqs  MeanSqs F.Model      R2 parOmegaSq    Pr(>F)
# lineage                         1   0.04145 0.041446 12.9383 0.11268   0.107612 9.999e-05 ***
# dominant_type                   3   0.00253 0.000845  0.2638 0.00689  -0.022819 0.9844016
# treat                           3   0.04381 0.014605  4.5592 0.11912   0.097355 0.0008999 ***
# sitename                        5   0.03527 0.007055  2.2023 0.09590   0.057244 0.0298970 *
# lineage:dominant_type           1   0.00518 0.005185  1.6186 0.01410   0.006209 0.1792821
# lineage:treat                   3   0.01632 0.005440  1.6981 0.04437   0.020716 0.1240876
# dominant_type:treat             6   0.01241 0.002068  0.6457 0.03374  -0.021943 0.7858214
# lineage:sitename                2   0.00338 0.001690  0.5277 0.00919  -0.009633 0.7169283
# dominant_type:sitename          4   0.00805 0.002014  0.6286 0.02190  -0.015234 0.7232277
# treat:sitename                 15   0.03039 0.002026  0.6325 0.08263  -0.058959 0.9148085
# lineage:dominant_type:treat     2   0.00254 0.001269  0.3963 0.00690  -0.012347 0.8391161
# lineage:dominant_type:sitename  1   0.00357 0.003574  1.1156 0.00972   0.001166 0.2891711
# lineage:treat:sitename          3   0.00568 0.001894  0.5912 0.01545  -0.012542 0.7184282
# dominant_type:treat:sitename    1   0.00345 0.003446  1.0757 0.00937   0.000764 0.2975702
# Residuals                      48   0.15376 0.003203         0.41804
# Total                          98   0.36781                  1.00000



#### Correlation Matrices ####
# these data are not included in the manuscript #
#source the rquery.cormat function, will need it below
source("http://www.sthda.com/upload/rquery_cormat.r")

#Function to calculate p-value of correlations for corrplot()
# mat : is a matrix of data
# ... : further arguments to pass to the native R cor.test function
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}


# use same files as for adonis test - not log transformed
t0_full = read.csv("/Users/hannahaichelman/Documents/BU/TVE/PCAs/t0_full.csv")
end_full = read.csv("/Users/hannahaichelman/Documents/BU/TVE/PCAs/end_full.csv")

# format t0 physiology data:
str(t0_full)

t0_full$sitename <- as.factor(t0_full$sitename)
t0_full$treat <- as.factor(t0_full$treat)
t0_full$reef <- as.factor(t0_full$reef)

# add in genet id and combine with lineage dataframe to include 2bRAD population data in our T0 PCAs
t0_full$gen_site <- substr(t0_full$frag,1,3)

lineages = read.csv("/Users/hannahaichelman/Documents/BU/TVE/2bRAD/Analysis/tuftscustompipeline_denovo_nosyms/tve_lineages_noclones.csv")

t0_full_lineage <- left_join(t0_full, lineages, by = "gen_site")
# all of the samples that have NA's for lineage are ones that didn't prep so we don't have lineage data for those genets

colnames(t0_full_lineage)[colnames(t0_full_lineage)=="hcarb_mgcm2"] <-"hcarb"
colnames(t0_full_lineage)[colnames(t0_full_lineage)=="scarb_mgcm2"] <-"scarb"
colnames(t0_full_lineage)[colnames(t0_full_lineage)=="prot_mgcm2"] <-"prot"
colnames(t0_full_lineage)[colnames(t0_full_lineage)=="sym_cm2"] <-"syms"
colnames(t0_full_lineage)[colnames(t0_full_lineage)=="avgtiss"] <-"tiss"
colnames(t0_full_lineage)[colnames(t0_full_lineage)=="corallite.avg.poly.mm2"] <-"corr_sa"

t0_full_all_lin = t0_full_lineage %>%
  drop_na()

t0_full_2_lin = t0_full_all_lin %>%
  filter(lineage!="L3")

# now end of variability data:
str(end_full)

end_full$sitename <- as.factor(end_full$sitename)
end_full$treat <- as.factor(end_full$treat)
end_full$reef <- as.factor(end_full$reef)

# re-level and re-name treatment
end_full$treat <- factor(end_full$treat, levels = c("Control","Low Var","Mod Var","High Var"))

# add in genet id and combine with lineage dataframe to include 2bRAD population data in our T0 PCAs
end_full$gen_site <- substr(end_full$frag,1,3)
end_full_lineage <- left_join(end_full, lineages, by = "gen_site")

# re-name columns for legibility on plots
colnames(end_full_lineage)[colnames(end_full_lineage)=="T2_T0_perc"] <-"growth"
colnames(end_full_lineage)[colnames(end_full_lineage)=="hcarb_mgcm2"] <-"hcarb"
colnames(end_full_lineage)[colnames(end_full_lineage)=="scarb_mgcm2"] <-"scarb"
colnames(end_full_lineage)[colnames(end_full_lineage)=="prot_mgcm2"] <-"prot"
colnames(end_full_lineage)[colnames(end_full_lineage)=="sym_cm2"] <-"syms"

end_full_all_lin = end_full_lineage %>%
  drop_na()

end_full_2_lin = end_full_all_lin %>%
  filter(lineage!="L3")

# add in dominant symbiont type dataframe
its2_types = read.csv("/Users/hannahaichelman/Documents/BU/TVE/16S_ITS2/ITS_PreStress_Timepoint/ITS2.dominanttype.csv") %>%
  select(frag, dominant_type) #%>%
#mutate(D1_1 = D1_sum+1) %>%
#select(-D1_sum) %>%
#mutate(propD1 = log(D1_1)) %>%
#select(-D1_1)
its2_types$dominant_type = as.factor(its2_types$dominant_type)

end_full_all_lin_symtype <- left_join(end_full_all_lin, its2_types, by = "frag")
str(end_full_all_lin_symtype)

end_full_2_lin_symtype <- left_join(end_full_2_lin, its2_types, by = "frag")
str(end_full_2_lin_symtype)

end_full_all_lin_symtype$dominant_type = as.factor(end_full_all_lin_symtype$dominant_type)
end_full_2_lin_symtype$dominant_type = as.factor(end_full_2_lin_symtype$dominant_type)


# df's
#t0_full_all_lin
#t0_full_2_lin
#end_full_all_lin_symtype
#end_full_2_lin_symtype

# make subsets for making correlation plots
t0_full_L1 = t0_full_2_lin %>%
  filter(lineage == "L1") %>%
  select(hcarb, scarb, prot, syms, tiss, chlA, corr_sa)

t0_full_L2 = t0_full_2_lin %>%
  filter(lineage == "L2") %>%
  select(hcarb, scarb, prot, syms, tiss, chlA, corr_sa)

end_full_L1_D1 = end_full_2_lin_symtype %>%
  filter(lineage == "L1") %>%
  filter(dominant_type == "D1") %>%
  select(growth, hcarb, scarb, prot, pam, syms, chlA) %>%
  drop_na()

end_full_L1_C1 = end_full_2_lin_symtype %>%
  filter(lineage == "L1") %>%
  filter(dominant_type == "C1") %>%
  select(growth, hcarb, scarb, prot, pam, syms, chlA) %>%
  drop_na()

end_full_L2 = end_full_2_lin_symtype %>%
  filter(lineage == "L2") %>%
  select(growth, hcarb, scarb, prot, pam, syms, chlA) %>%
  drop_na()

## Lineage 1 correlations
#compute correlation matrix
end_L1_m <- cor(end_full_L1) #create a correlation matrix of data
end_full_L1_p.mat <- cor.mtest(end_full_L1) #calculate significance of correlations, confidence level = 0.95

end_L1_D1_m <- cor(end_full_L1_D1) #create a correlation matrix of data
end_full_L1_D1_p.mat <- cor.mtest(end_full_L1_D1) #calculate significance of correlations, confidence level = 0.95

end_L1_C1_m <- cor(end_full_L1_C1) #create a correlation matrix of data
end_full_L1_C1_p.mat <- cor.mtest(end_full_L1_C1) #calculate significance of correlations, confidence level = 0.95

#Use the rquery.cormat() function to get a table of correlation coefficients and p-values
end_full_L1_table <- rquery.cormat(end_full_L1, type="flatten", graph=FALSE)
end_full_L1_D1_table <- rquery.cormat(end_full_L1_D1, type="flatten", graph=FALSE)
end_full_L1_C1_table <- rquery.cormat(end_full_L1_C1, type="flatten", graph=FALSE)

#plot correlation matrix
#positive correlations in blue and negative correlations are in red. Color intensity and size of the circle are proportional to the correlation coefficients
pdf(file = "/Users/hannahaichelman/Documents/BU/TVE/CorrelationMatrix/L1_end_corrmat.pdf", width=4, height=4, useDingbats=FALSE)
par(mfrow=c(1,3))
corrplot(end_L1_m, method = "circle",
         type = "upper",
         order="alphabet", #hierarchical clustering order of factors
         tl.col="black", tl.srt=45, #text label color and rotation
         p.mat = end_full_L1_p.mat, sig.level = 0.05, insig = "blank", #combine with significance
         title = "L1-all" #can add title, but adds it really high for some reason
)

corrplot(end_L1_D1_m, method = "circle",
         type = "upper",
         order="alphabet", #hierarchical clustering order of factors
         tl.col="black", tl.srt=45, #text label color and rotation
         p.mat = end_full_L1_D1_p.mat, sig.level = 0.05, insig = "blank", #combine with significance
         title = "L1-D1" #can add title, but adds it really high for some reason
)

corrplot(end_L1_C1_m, method = "circle",
         type = "upper",
         order="alphabet", #hierarchical clustering order of factors
         tl.col="black", tl.srt=45, #text label color and rotation
         p.mat = end_full_L1_C1_p.mat, sig.level = 0.05, insig = "blank", #combine with significance
         title = "L1-C1" #can add title, but adds it really high for some reason
)

dev.off()

# T0
#compute correlation matrix
t0_L1_m <- cor(t0_full_L1) #create a correlation matrix of data
t0_full_L1_p.mat <- cor.mtest(t0_full_L1) #calculate significance of correlations, confidence level = 0.95

#Use the rquery.cormat() function to get a table of correlation coefficients and p-values
t0_full_L1_table <- rquery.cormat(t0_full_L1, type="flatten", graph=FALSE)

#plot correlation matrix
#positive correlations in blue and negative correlations are in red. Color intensity and size of the circle are proportional to the correlation coefficients
pdf(file = "/Users/hannahaichelman/Documents/BU/TVE/CorrelationMatrix/L1_t0_corrmat.pdf", width=4, height=4, useDingbats=FALSE)
corrplot(t0_L1_m, method = "circle",
         type = "upper",
         order="alphabet", #hierarchical clustering order of factors
         tl.col="black", tl.srt=45, #text label color and rotation
         p.mat = t0_full_L1_p.mat, sig.level = 0.05, insig = "blank" #combine with significance
         #title = "T30" #can add title, but adds it really high for some reason
)
dev.off()

## Lineage 2 correlations
#compute correlation matrix
end_L2_m <- cor(end_full_L2) #create a correlation matrix of data
end_full_L2_p.mat <- cor.mtest(end_full_L2) #calculate significance of correlations, confidence level = 0.95

#Use the rquery.cormat() function to get a table of correlation coefficients and p-values
end_full_L2_table <- rquery.cormat(end_full_L2, type="flatten", graph=FALSE)

#plot correlation matrix
#positive correlations in blue and negative correlations are in red. Color intensity and size of the circle are proportional to the correlation coefficients
pdf(file = "/Users/hannahaichelman/Documents/BU/TVE/CorrelationMatrix/L2_end_corrmat.pdf", width=4, height=4, useDingbats=FALSE)
corrplot(end_L2_m, method = "circle",
         type = "upper",
         order="alphabet", #hierarchical clustering order of factors
         tl.col="black", tl.srt=45, #text label color and rotation
         p.mat = end_full_L2_p.mat, sig.level = 0.05, insig = "blank" #combine with significance
         #title = "T30" #can add title, but adds it really high for some reason
)
dev.off()

#t0
#compute correlation matrix
t0_L2_m <- cor(t0_full_L2) #create a correlation matrix of data
t0_full_L2_p.mat <- cor.mtest(t0_full_L2) #calculate significance of correlations, confidence level = 0.95

#Use the rquery.cormat() function to get a table of correlation coefficients and p-values
t0_full_L2_table <- rquery.cormat(t0_full_L2, type="flatten", graph=FALSE)

#plot correlation matrix
#positive correlations in blue and negative correlations are in red. Color intensity and size of the circle are proportional to the correlation coefficients
pdf(file = "/Users/hannahaichelman/Documents/BU/TVE/CorrelationMatrix/L2_t0_corrmat.pdf", width=4, height=4, useDingbats=FALSE)
corrplot(t0_L2_m, method = "circle",
         type = "upper",
         order="alphabet", #hierarchical clustering order of factors
         tl.col="black", tl.srt=45, #text label color and rotation
         p.mat = t0_full_L2_p.mat, sig.level = 0.05, insig = "blank" #combine with significance
         #title = "T30" #can add title, but adds it really high for some reason
)
dev.off()

#### JP Coring Data ####
# these data are not included in the manuscript #
cores = read.csv("/Users/hannahaichelman/Documents/BU/TVE/CoringData/Growth_allRegions_JPcoringdata.csv")

head(cores)

pan_cores = cores %>%
  filter(region == "w") %>%
  filter(spp == "s")

str(pan_cores)

pan_cores$rz = as.factor(pan_cores$rz)
pan_cores$site = as.factor(pan_cores$site)
pan_cores$spp = as.factor(pan_cores$spp)

# Re-name factor levels for sites
pan_cores$site <- recode(pan_cores$site,
                         wirci = 'CI',
                         wirpd = 'PD',
                         wirpl = 'PL',
                         wirsp = 'SP',
                         worbn = 'BN',
                         worbs = 'BS',
                         worca = 'CA',
                         wordm = 'DM')

# Remove sites that we don't have corals from for this experiment
pan_cores_tve = pan_cores %>%
  subset(site != "PL") %>%
  subset(site != "DM")

# Subset to only include recent years (more confident in these data) as in JP's analysis
pan_cores_tve_recent = pan_cores_tve %>%
  dplyr::filter(Year >= "2005" & Year <= "2014")

# Average density, calcification, and linear extension by site and plot

sitemeans_linext = summarySE(data = pan_cores_tve_recent, measurevar = c("linext"), groupvars = c("site","rz"))
sitemeans_den = summarySE(data = pan_cores_tve_recent, measurevar = c("density"), groupvars = c("site","rz"))
sitemeans_calc = summarySE(data = pan_cores_tve_recent, measurevar = c("calc"), groupvars = c("site","rz"))

# plot linear extension, treatment x axis colored by site data figure
p.linext <- ggplot(sitemeans_linext,aes(x = site, y = linext, color = site, pch = site))+
  theme_bw()+
  geom_point(size = 3, position = position_dodge(width=0.3))+
  geom_errorbar(aes(x = site, ymax = linext+se, ymin = linext-se), width = .2, position = position_dodge(width=0.3)) +
  scale_color_manual(name = "Site",
                     labels = c("CI","PD","SP","BN","BS","CA"),
                     values = cols_site)+
  scale_shape_manual(name = "Site",
                     labels = c("CI","PD","SP","BN","BS","CA"),
                     values=c(19,19,19,17,17,17))+
  xlab("Site Name")+
  ylab("Linear Extension")+
  #ylim(3.75,6) +
  #geom_vline(xintercept = 1.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), legend.position = "none")
p.linext

p.den <- ggplot(sitemeans_den,aes(x = site, y = density, color = site, pch = site))+
  theme_bw()+
  geom_point(size = 3, position = position_dodge(width=0.3))+
  geom_errorbar(aes(x = site, ymax = density+se, ymin = density-se), width = .2, position = position_dodge(width=0.3)) +
  scale_color_manual(name = "Site",
                     labels = c("CI","PD","SP","BN","BS","CA"),
                     values = cols_site)+
  scale_shape_manual(name = "Site",
                     labels = c("CI","PD","SP","BN","BS","CA"),
                     values=c(19,19,19,17,17,17))+
  xlab("Site Name")+
  ylab("Density")+
  #ylim(3.75,6) +
  #geom_vline(xintercept = 1.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), legend.position = "none")
p.den

p.calc <- ggplot(sitemeans_calc,aes(x = site, y = calc, color = site, pch = site))+
  theme_bw()+
  geom_point(size = 3, position = position_dodge(width=0.3))+
  geom_errorbar(aes(x = site, ymax = calc+se, ymin = calc-se), width = .2, position = position_dodge(width=0.3)) +
  scale_color_manual(name = "Site",
                     labels = c("CI","PD","SP","BN","BS","CA"),
                     values = cols_site)+
  scale_shape_manual(name = "Site",
                     labels = c("CI","PD","SP","BN","BS","CA"),
                     values=c(19,19,19,17,17,17))+
  xlab("Site Name")+
  ylab("Calcification")+
  #ylim(3.75,6) +
  #geom_vline(xintercept = 1.5) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), legend.position = "none")
p.calc


core.plots = ggarrange(p.linext, p.den, p.calc,
                       labels = c("A", "B", "C"),
                       ncol = 3, nrow = 1)

ggsave(core.plots, filename = "/Users/hannahaichelman/Documents/BU/TVE/CoringData/coring_plots_2005-2014.pdf", width=8, height=3, units=c("in"), useDingbats=FALSE)
