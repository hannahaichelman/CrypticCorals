#### Coring Data ####
cores = read.csv("Physiology_Data/data_files/Growth_allRegions_JPcoringdata.csv")
head(cores)

pan_cores = cores %>%
  filter(region == "w") %>%
  filter(spp == "s")

str(pan_cores)

pan_cores$rz = as.factor(pan_cores$rz)
pan_cores$site = as.factor(pan_cores$site)
pan_cores$spp = as.factor(pan_cores$spp)

# Re-name factor levels for sites
pan_cores$site <- dplyr::recode(pan_cores$site,
                                wirci = 'CI',
                                wirpd = 'PD',
                                wirpl = 'PL',
                                wirsp = 'SP',
                                worbn = 'BN',
                                worbs = 'BS',
                                worca = 'CA',
                                wordm = 'DM')
# Re-name core IDs
pan_cores$coreID <- dplyr::recode(pan_cores$coreID,
                                  wordm37s='W_OR_DM_37_S',
                                  worbn15s='W_OR_BN_15_S',
                                  worbs17s='W_OR_BS_17_S',
                                  worca35s='W_OR_CA_35_S',
                                  wirsp05s='W_IR_SP_5_S',
                                  worbn09s='W_OR_BN_9_S',
                                  wirsp02s='W_IR_SP_2_S',
                                  worca27s='W_OR_CA_27_S',
                                  wirpd59s='W_IR_PD_59_S',
                                  worbn14s='W_OR_BN_14_S',
                                  wordm42s='W_OR_DM_42_S',
                                  worbn16s='W_OR_BN_16_S',
                                  worbs23s='W_OR_BS_23_S',
                                  worbs22s='W_OR_BS_22_S',
                                  wirsp03s='W_IR_SP_3_S',
                                  worca30s='W_OR_CA_30_S',
                                  worbs25s='W_OR_BS_25_S',
                                  wordm40s='W_OR_DM_40_S',
                                  worbn07s='W_OR_BN_7_S',
                                  wordm41s='W_OR_DM_41_S',
                                  wirpd63s='W_IR_PD_63_S',
                                  wirci57s='W_IR_CI_57_S',
                                  wirpl51s='W_IR_PL_51_S',
                                  wirpd60s='W_IR_PD_60_S',
                                  wirpd61s='W_IR_PD_61_S',
                                  wirpl47s='W_IR_PL_47_S',
                                  wirsp04s='W_IR_SP_4_S',
                                  worbs20s='W_OR_BS_20_S',
                                  wordm38s='W_OR_DM_38_S',
                                  wirpl48s='W_IR_PL_48_S')

# W_IR_PL_47_S and W_IR_PD_63_S don't have growth chronologies (bad data) even though we have lineage assignments
# W_OR_BS_17_S had low coverage, no lineage assignment
# W_OR_BS_20_S, W_OR_DM_38_S, W_IR_PL_48_S all didn't sequence, no lineage assignment
# N=24 total cores with lineage assignment and growth chronology

# take a look at the metadata that Annabel gave us
#loading in RAD data with core ID and bam file name
rad_meta <- read.csv("Physiology_Data/data_files/all_SSID_updated.csv")
rad_meta <- filter(rad_meta, Project=="cores")
head(rad_meta)
dim(rad_meta)
#[1] 28 12

# loading in and cleaning up lineage assignment data to only include core info
lineages <- read.csv("Physiology_Data/data_files/samples_on_plate_somelineages.csv")
colnames(lineages)[1] ="SampleID"
colnames(lineages)[3] ="lineage"
lineages <- dplyr::filter(lineages, PROJECT=="Panama cores")
head(lineages)
dim(lineages)
#[1] 30  4


lineages$SampleID %in% rad_meta$SampleID
# the 3 samples that didn't sequence aren't included in the rad_meta file because they didn't demultiplex

# looks like there isn't a reason to use the rad_meta info, so will move forward with just lineage info

# loading in and cleaning up lineage assignment data to only include core info
lineages <- read.csv("Physiology_Data/data_files/samples_on_plate_somelineages.csv")
lineages2 = lineages %>%
  dplyr::rename(coreID = SAMPLE_ID) %>%
  dplyr::rename(lineage = X) %>%
  dplyr::filter(PROJECT=="Panama cores") %>%
  select(coreID, lineage)
head(lineages2)
dim(lineages2)

#combining JP metadata with lineage associations
lineages_meta <- left_join(pan_cores, lineages2, by="coreID")
head(lineages_meta)
str(lineages_meta)
lineages_meta$coreID = as.factor(lineages_meta$coreID)
lineages_meta$lineage = as.factor(lineages_meta$lineage)

## Data Analysis
# Subset to only include recent years (more confident in these core data) as in JP's analysis
pan_cores_recent = lineages_meta %>%
  dplyr::filter(Year >= "1980" & Year <= "2014") %>%
  drop_na(lineage) %>%
  dplyr::filter(lineage != "3")
head(pan_cores_recent)
count(unique(pan_cores_recent$coreID)) # lines up with 24 cores with lineage and growth data - phew.

pan_cores_allyears = lineages_meta %>%
  drop_na(lineage) %>%
  dplyr::filter(lineage != "3")
count(unique(pan_cores_allyears$coreID)) # lines up with 24 cores with lineage and growth data - phew.

# make a histogram of how many cores we have through time
cols_lineage <- c("#3f007d","#807dba")

str(pan_cores_allyears)
yearmeans = summarySE(data = pan_cores_allyears, measurevar = c("calc"), groupvars = c("Year","lineage"))
yearmeans = yearmeans %>%
  select(Year, lineage, N)

samplesize_hist =
  ggplot(yearmeans) +
  theme_bw()+
  geom_bar(aes(x=Year, y=N, fill = lineage), stat="identity") +
  scale_fill_manual(name = "Lineage",
                    breaks = c("1","2"),
                    values = cols_lineage)+
  scale_x_continuous(name="Year",breaks = seq(1880,2014,by=10))+
  facet_wrap(~lineage)
samplesize_hist
ggsave(samplesize_hist, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/CoringData/coring_samplesize.pdf", width=10, height=3, units=c("in"), useDingbats=FALSE)


# Average density, calcification, and linear extension by lineage only - recent data only
sitemeans_linext_all_recent = summarySE(data = pan_cores_recent, measurevar = c("linext"), groupvars = c("lineage"))
sitemeans_den_all_recent = summarySE(data = pan_cores_recent, measurevar = c("density"), groupvars = c("lineage"))
sitemeans_calc_all_recent = summarySE(data = pan_cores_recent, measurevar = c("calc"), groupvars = c("lineage"))

sitemeans_linext_cores_recent = summarySE(data = pan_cores_recent, measurevar = c("linext"), groupvars = c("lineage","coreID","site"))
sitemeans_den_cores_recent = summarySE(data = pan_cores_recent, measurevar = c("density"), groupvars = c("lineage","coreID","site"))
sitemeans_calc_cores_recent = summarySE(data = pan_cores_recent, measurevar = c("calc"), groupvars = c("lineage","coreID","site"))

shapes_lineage = c("CI" = 22,"PD" = 21,"PL" = 23,"SP" = 24,"BN" = 15,"BS" = 16,"CA" = 17,"DM" = 25)

p.linext_recent <- ggplot(sitemeans_linext_cores_recent, aes(x = lineage, y = linext, group = lineage))+
  theme_bw()+
  geom_jitter(aes(color = lineage, shape = site),
              width=0.1,
              alpha=0.6,
              size = 1.5)+
  geom_errorbar(data = sitemeans_linext_all_recent, aes(x = lineage, ymax = linext+se, ymin = linext-se, color = lineage), width = .2, position = position_dodge(width=0.3)) +
  geom_point(data = sitemeans_linext_all_recent, mapping = aes(x = lineage, y = linext, color = lineage, fill = lineage), position = position_dodge(width=0.3), size = 3.5, pch = 21, color = 'black')+
  xlab("Lineage")+
  ylab("Linear Extension")+
  scale_fill_manual(name = "Lineage",
                    breaks = c("1","2"),
                    values = cols_lineage)+
  scale_color_manual(name = "Lineage",
                     breaks = c("1","2"),
                     values = cols_lineage)+
  scale_shape_manual(name = "Site",
                     breaks = c("CI", "PD", "PL", "SP", "BN", "BS", "CA", "DM"),
                     values = shapes_lineage) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))
#theme(legend.position = "none")
p.linext_recent

p.den_recent <- ggplot(sitemeans_den_cores_recent, aes(x = lineage, y = density, group = lineage))+
  theme_bw()+
  geom_jitter(aes(color = lineage, shape = site),
              width=0.1,
              alpha=0.6,
              size = 1.5)+
  geom_errorbar(data = sitemeans_den_all_recent, aes(x = lineage, ymax = density+se, ymin = density-se, color = lineage), width = .2, position = position_dodge(width=0.4)) +
  geom_point(data = sitemeans_den_all_recent, mapping = aes(x = lineage, y = density, color = lineage, fill = lineage), position = position_dodge(width=0.4), size = 3.5, pch = 21, color = 'black')+
  xlab("Lineage")+
  ylab("Density")+
  scale_fill_manual(name = "Lineage",
                    breaks = c("1","2"),
                    values = cols_lineage)+
  scale_color_manual(name = "Lineage",
                     breaks = c("1","2"),
                     values = cols_lineage)+
  scale_shape_manual(name = "Site",
                     breaks = c("CI", "PD", "PL", "SP", "BN", "BS", "CA", "DM"),
                     values = shapes_lineage) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))
#theme(legend.position = "none")
p.den_recent

p.calc_recent <- ggplot(sitemeans_calc_cores_recent, aes(x = lineage, y = calc, group = lineage))+
  theme_bw()+
  geom_jitter(aes(color = lineage, shape = site),
              width=0.1,
              alpha=0.6,
              size = 1.5)+
  geom_errorbar(data = sitemeans_calc_all_recent, aes(x = lineage, ymax = calc+se, ymin = calc-se, color = lineage), width = .2, position = position_dodge(width=0.4)) +
  geom_point(data = sitemeans_calc_all_recent, mapping = aes(x = lineage, y = calc, color = lineage, fill = lineage), position = position_dodge(width=0.4), size = 3.5, pch = 21, color = 'black')+
  xlab("Lineage")+
  ylab("Calcification")+
  scale_fill_manual(name = "Lineage",
                    breaks = c("1","2"),
                    values = cols_lineage)+
  scale_color_manual(name = "Lineage",
                     breaks = c("1","2"),
                     values = cols_lineage)+
  scale_shape_manual(name = "Site",
                     breaks = c("CI", "PD", "PL", "SP", "BN", "BS", "CA", "DM"),
                     values = shapes_lineage) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))
#theme(legend.position = "none")
p.calc_recent

core.plots_recent = ggarrange(p.linext_recent, p.den_recent, p.calc_recent,
                              labels = c("A", "B", "C"),
                              ncol = 3, nrow = 1,
                              common.legend = TRUE)
core.plots_recent
ggsave(core.plots_recent, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/CoringData/coring_plots_1980-2014.pdf", width=8, height=4, units=c("in"), useDingbats=FALSE)

# stats for these recent data
calc_recent_lm = lm(calc~lineage, data = sitemeans_calc_cores_recent)
summary(calc_recent_lm)

linext_recent_lm = lm(linext~lineage, data = sitemeans_linext_cores_recent)
summary(linext_recent_lm)

density_recent_lm = lm(density~lineage, data = sitemeans_den_cores_recent)
summary(density_recent_lm)

# try stats with random effect of core:
calc_recent_lmer = lmer(calc~lineage + (1|coreID), data = pan_cores_recent)
summary(calc_recent_lmer)

linext_recent_lmer = lmer(linext~lineage + (1|coreID), data = pan_cores_recent)
summary(linext_recent_lmer)

density_recent_lmer = lmer(density~lineage + (1|coreID), data = pan_cores_recent)
summary(density_recent_lmer)

# Average density, calcification, and linear extension by lineage only - All years

sitemeans_linext_all = summarySE(data = pan_cores_allyears, measurevar = c("linext"), groupvars = c("lineage"))
sitemeans_den_all = summarySE(data = pan_cores_allyears, measurevar = c("density"), groupvars = c("lineage"))
sitemeans_calc_all = summarySE(data = pan_cores_allyears, measurevar = c("calc"), groupvars = c("lineage"))

sitemeans_linext_cores = summarySE(data = pan_cores_allyears, measurevar = c("linext"), groupvars = c("lineage","coreID","site"))
sitemeans_den_cores = summarySE(data = pan_cores_allyears, measurevar = c("density"), groupvars = c("lineage","coreID","site"))
sitemeans_calc_cores = summarySE(data = pan_cores_allyears, measurevar = c("calc"), groupvars = c("lineage","coreID","site"))

p.linext_all <- ggplot(sitemeans_linext_cores, aes(x = lineage, y = linext, group = lineage))+
  theme_bw()+
  geom_jitter(aes(color = lineage, shape = site),
              width=0.1,
              alpha=0.6,
              size = 1.5)+
  geom_errorbar(data = sitemeans_linext_all, aes(x = lineage, ymax = linext+se, ymin = linext-se, color = lineage), width = .2, position = position_dodge(width=0.3)) +
  geom_point(data = sitemeans_linext_all, mapping = aes(x = lineage, y = linext, color = lineage, fill = lineage), position = position_dodge(width=0.3), size = 3.5, pch = 21, color = 'black')+
  xlab("Lineage")+
  ylab("Linear Extension")+
  scale_fill_manual(name = "Lineage",
                    breaks = c("1","2"),
                    values = cols_lineage)+
  scale_color_manual(name = "Lineage",
                     breaks = c("1","2"),
                     values = cols_lineage)+
  scale_shape_manual(name = "Site",
                     breaks = c("CI", "PD", "PL", "SP", "BN", "BS", "CA", "DM"),
                     values = shapes_lineage) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))
#theme(legend.position = "none")
p.linext_all

p.den_all <- ggplot(sitemeans_den_cores, aes(x = lineage, y = density, group = lineage))+
  theme_bw()+
  geom_jitter(aes(color = lineage, shape = site),
              width=0.1,
              alpha=0.6,
              size = 1.5)+
  geom_errorbar(data = sitemeans_den_all, aes(x = lineage, ymax = density+se, ymin = density-se, color = lineage), width = .2, position = position_dodge(width=0.4)) +
  geom_point(data = sitemeans_den_all, mapping = aes(x = lineage, y = density, color = lineage, fill = lineage), position = position_dodge(width=0.4), size = 3.5, pch = 21, color = 'black')+
  xlab("Lineage")+
  ylab("Density")+
  scale_fill_manual(name = "Lineage",
                    breaks = c("1","2"),
                    values = cols_lineage)+
  scale_color_manual(name = "Lineage",
                     breaks = c("1","2"),
                     values = cols_lineage)+
  scale_shape_manual(name = "Site",
                     breaks = c("CI", "PD", "PL", "SP", "BN", "BS", "CA", "DM"),
                     values = shapes_lineage) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))
#theme(legend.position = "none")
p.den_all

p.calc_all <- ggplot(sitemeans_calc_cores, aes(x = lineage, y = calc, group = lineage))+
  theme_bw()+
  geom_jitter(aes(color = lineage, shape = site),
              width=0.1,
              alpha=0.6,
              size = 1.5)+
  geom_errorbar(data = sitemeans_calc_all, aes(x = lineage, ymax = calc+se, ymin = calc-se, color = lineage), width = .2, position = position_dodge(width=0.4)) +
  geom_point(data = sitemeans_calc_all, mapping = aes(x = lineage, y = calc, color = lineage, fill = lineage), position = position_dodge(width=0.4), size = 3.5, pch = 21, color = 'black')+
  xlab("Lineage")+
  ylab("Calcification")+
  scale_fill_manual(name = "Lineage",
                    breaks = c("1","2"),
                    values = cols_lineage)+
  scale_color_manual(name = "Lineage",
                     breaks = c("1","2"),
                     values = cols_lineage)+
  scale_shape_manual(name = "Site",
                     breaks = c("CI", "PD", "PL", "SP", "BN", "BS", "CA", "DM"),
                     values = shapes_lineage) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))
#theme(legend.position = "none")
p.calc_all
core.plots_all = ggarrange(p.linext_all, p.den_all, p.calc_all,
                           labels = c("A", "B", "C"),
                           ncol = 3, nrow = 1, common.legend = TRUE)
core.plots_all
ggsave(core.plots_all, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/CoringData/coring_plots_allyears.pdf", width=8, height=4, units=c("in"), useDingbats=FALSE)

# stats for all data
calc_all_lm = lm(calc~lineage, data = sitemeans_calc_cores)
summary(calc_all_lm)

linext_all_lm = lm(linext~lineage, data = sitemeans_linext_cores)
summary(linext_all_lm)

density_all_lm = lm(density~lineage, data = sitemeans_den_cores)
summary(density_all_lm)

## Visualizing growth trends across lineages - all years

linext_trends = summarySE(data = pan_cores_allyears, measurevar = c("linext"), groupvars = c("lineage","Year"))
den_trends = summarySE(data = pan_cores_allyears, measurevar = c("density"), groupvars = c("lineage","Year"))
calc_trends = summarySE(data = pan_cores_allyears, measurevar = c("calc"), groupvars = c("lineage","Year"))

# Linear extension
linext = ggplot(data=linext_trends, aes(x=Year, y=linext, color=lineage)) +
  theme_bw()+
  geom_point() +
  geom_smooth(aes(fill=lineage)) +
  scale_color_manual(name = "Lineage",
                     breaks = c("1","2"),
                     values = cols_lineage)+
  scale_fill_manual(name = "Lineage",
                    breaks = c("1","2"),
                    values = cols_lineage)+
  ylab("Linear Extension")+
  theme(legend.position="none")
linext

# Skeletal density
den = ggplot(data=den_trends, aes(x=Year, y=density, color=lineage)) +
  theme_bw()+
  geom_point() +
  geom_smooth(aes(fill=lineage)) +
  scale_color_manual(name = "Lineage",
                     breaks = c("1","2"),
                     values = cols_lineage)+
  scale_fill_manual(name = "Lineage",
                    breaks = c("1","2"),
                    values = cols_lineage)+
  ylab("Density") +
  theme(legend.position="none")
den

# Calcification
calc = ggplot(data=calc_trends, aes(x=Year, y=calc, color=lineage)) +
  theme_bw()+
  geom_point() +
  geom_smooth(aes(fill=lineage)) +
  scale_color_manual(name = "Lineage",
                     breaks = c("1","2"),
                     values = cols_lineage) +
  scale_fill_manual(name = "Lineage",
                    breaks = c("1","2"),
                    values = cols_lineage)+
  ylab("Calcification") +
  theme(legend.position="none")
calc

trends = ggarrange(linext, den, calc,
                   labels = c("A", "B", "C"),
                   ncol = 3, nrow = 1)
trends

ggsave(trends, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/CoringData/coring_trends_alltime.pdf", width=8, height=4, units=c("in"), useDingbats=FALSE)

# now looking at growth trends only recent years (1980-2014)
linext_recent = summarySE(data = pan_cores_recent, measurevar = c("linext"), groupvars = c("lineage","Year"))
den_recent = summarySE(data = pan_cores_recent, measurevar = c("density"), groupvars = c("lineage","Year"))
calc_recent = summarySE(data = pan_cores_recent, measurevar = c("calc"), groupvars = c("lineage","Year"))

# Linear extension
linext = ggplot(data=linext_recent, aes(x=Year, y=linext, color=lineage)) +
  theme_bw()+
  geom_point() +
  geom_smooth(aes(fill=lineage)) +
  scale_color_manual(name = "Lineage",
                     breaks = c("1","2"),
                     values = cols_lineage)+
  scale_fill_manual(name = "Lineage",
                    breaks = c("1","2"),
                    values = cols_lineage)+
  ylab("Linear Extension")+
  scale_x_continuous(name="Year",breaks = seq(1980,2014,by=10))+
  theme(legend.position="none")
linext

# Skeletal density
den = ggplot(data=den_recent, aes(x=Year, y=density, color=lineage)) +
  theme_bw()+
  geom_point() +
  geom_smooth(aes(fill=lineage)) +
  scale_color_manual(name = "Lineage",
                     breaks = c("1","2"),
                     values = cols_lineage)+
  scale_fill_manual(name = "Lineage",
                    breaks = c("1","2"),
                    values = cols_lineage)+
  ylab("Density") +
  scale_x_continuous(name="Year",breaks = seq(1980,2014,by=10))+
  theme(legend.position="none")
den

# Calcification
calc = ggplot(data=calc_recent, aes(x=Year, y=calc, color=lineage)) +
  theme_bw()+
  geom_point() +
  geom_smooth(aes(fill=lineage)) +
  scale_color_manual(name = "Lineage",
                     breaks = c("1","2"),
                     values = cols_lineage) +
  scale_fill_manual(name = "Lineage",
                    breaks = c("1","2"),
                    values = cols_lineage)+
  ylab("Calcification") +
  scale_x_continuous(name="Year",breaks = seq(1980,2014,by=10))+
  theme(legend.position="none")
calc

trends_recent = ggarrange(linext, den, calc,
                          labels = c("A", "B", "C"),
                          ncol = 3, nrow = 1)
trends_recent

ggsave(trends_recent, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/CoringData/coring_trends_1980_2014.pdf", width=8, height=4, units=c("in"), useDingbats=FALSE)



# Linear extension by site with jitter points by lineage

linext_allyears = summarySE(data = pan_cores_allyears, measurevar = c("linext"), groupvars = c("lineage","site","coreID"))
str(linext_allyears)

linext_plot_lineage <- ggplot(linext_allyears, aes(x = site, y = linext, group=coreID))+
  theme_bw()+
  geom_jitter(aes(color = lineage, fill = lineage),
              position=position_dodge(width=0.3),
              alpha=0.2, pch = 21,
              color = "black") +
  geom_errorbar(data = linext_allyears, aes(x = site, ymax = linext+se, ymin = linext-se, color = lineage), width = .2, position = position_dodge(width=0.4)) +
  geom_point(data = linext_allyears, mapping = aes(x=site, y=linext, color = lineage, fill = lineage, group=coreID), size = 3.5, pch = 21, color = "black", position = position_dodge(width=0.4))+
  scale_fill_manual(name = "Lineage",
                    breaks = c("1","2"),
                    values = cols_lineage)+
  scale_color_manual(name = "Lineage",
                     breaks = c("1","2"),
                     values = cols_lineage)+
  xlab("Site")+
  ylab("Linear Extension")+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))+
  theme(legend.position = "blank")
linext_plot_lineage

# Skeletal density by site with jitter points by lineage

den_allyears = summarySE(data = pan_cores_allyears, measurevar = c("density"), groupvars = c("lineage","site","coreID"))

den_plot_lineage <- ggplot(den_allyears,aes(x = site, y = density, group=coreID))+
  theme_bw()+
  geom_jitter(aes(color = lineage, fill = lineage),
              position=position_dodge(width=0.3),
              alpha=0.2, pch = 21,
              color = "black") +
  geom_errorbar(data = den_allyears, aes(x = site, ymax = density+se, ymin = density-se, color = lineage), width = .2, position = position_dodge(width=0.4)) +
  geom_point(data = den_allyears, mapping = aes(x=site, y=density, color = lineage, fill = lineage, group=coreID), size = 3.5, pch = 21, color = "black", position = position_dodge(width=0.4))+
  scale_fill_manual(name = "Lineage",
                    breaks = c("1","2"),
                    values = cols_lineage)+
  scale_color_manual(name = "Lineage",
                     breaks = c("1","2"),
                     values = cols_lineage)+
  xlab("Site")+
  ylab("Skeletal Density")+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  theme(legend.position = "blank")

den_plot_lineage

# Calcification by site with jitter points by lineage

calc_allyears = summarySE(data = pan_cores_allyears, measurevar = c("calc"), groupvars = c("lineage","site","coreID"))

calc_plot_lineage <- ggplot(calc_allyears,aes(x = site, y = calc, group=coreID))+
  theme_bw()+
  geom_jitter(aes(color = lineage, fill = lineage),
              position=position_dodge(width=0.3),
              alpha=0.2, pch = 21,
              color = "black") +
  geom_errorbar(data = calc_allyears, aes(x = site, ymax = calc+se, ymin = calc-se, color = lineage), width = .2, position = position_dodge(width=0.4)) +
  geom_point(data = calc_allyears, mapping = aes(x=site, y=calc, color = lineage, fill = lineage, group=coreID), size = 3.5, pch = 21, color = "black", position = position_dodge(width=0.4))+
  scale_fill_manual(name = "Lineage",
                    breaks = c("1","2"),
                    values = cols_lineage)+
  scale_color_manual(name = "Lineage",
                     breaks = c("1","2"),
                     values = cols_lineage)+
  xlab("Site")+
  ylab("Calcification")+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))+
  theme(legend.position = "blank")
calc_plot_lineage

plots_allyears = ggarrange(linext_plot_lineage, den_plot_lineage, calc_plot_lineage,
                           labels = c("A", "B", "C"),
                           ncol = 3, nrow = 1)
plots_allyears

ggsave(plots_allyears, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/CoringData/coring_plots_individualcores_alltime.pdf", width=9, height=4, units=c("in"), useDingbats=FALSE)


# same 3 panels, but now with only recent data
# Linear extension by site with jitter points by lineage

linext_recent = summarySE(data = pan_cores_recent, measurevar = c("linext"), groupvars = c("lineage","site","coreID"))
str(linext_recent)

linext_plot_lineage <- ggplot(linext_recent, aes(x = site, y = linext, group=coreID))+
  theme_bw()+
  geom_jitter(aes(color = lineage, fill = lineage),
              position=position_dodge(width=0.3),
              alpha=0.2, pch = 21,
              color = "black") +
  geom_errorbar(data = linext_recent, aes(x = site, ymax = linext+se, ymin = linext-se, color = lineage), width = .2, position = position_dodge(width=0.4)) +
  geom_point(data = linext_recent, mapping = aes(x=site, y=linext, color = lineage, fill = lineage, group=coreID), size = 3.5, pch = 21, color = "black", position = position_dodge(width=0.4))+
  scale_fill_manual(name = "Lineage",
                    breaks = c("1","2"),
                    values = cols_lineage)+
  scale_color_manual(name = "Lineage",
                     breaks = c("1","2"),
                     values = cols_lineage)+
  xlab("Site")+
  ylab("Linear Extension")+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))+
  theme(legend.position = "blank")
linext_plot_lineage

# Skeletal density by site with jitter points by lineage

den_recent = summarySE(data = pan_cores_recent, measurevar = c("density"), groupvars = c("lineage","site","coreID"))

den_plot_lineage <- ggplot(den_recent,aes(x = site, y = density, group=coreID))+
  theme_bw()+
  geom_jitter(aes(color = lineage, fill = lineage),
              position=position_dodge(width=0.3),
              alpha=0.2, pch = 21,
              color = "black") +
  geom_errorbar(data = den_recent, aes(x = site, ymax = density+se, ymin = density-se, color = lineage), width = .2, position = position_dodge(width=0.4)) +
  geom_point(data = den_recent, mapping = aes(x=site, y=density, color = lineage, fill = lineage, group=coreID), size = 3.5, pch = 21, color = "black", position = position_dodge(width=0.4))+
  scale_fill_manual(name = "Lineage",
                    breaks = c("1","2"),
                    values = cols_lineage)+
  scale_color_manual(name = "Lineage",
                     breaks = c("1","2"),
                     values = cols_lineage)+
  xlab("Site")+
  ylab("Skeletal Density")+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  theme(legend.position = "blank")
den_plot_lineage

# Calcification by site with jitter points by lineage

calc_recent = summarySE(data = pan_cores_recent, measurevar = c("calc"), groupvars = c("lineage","site","coreID"))

calc_plot_lineage <- ggplot(calc_recent,aes(x = site, y = calc, group=coreID))+
  theme_bw()+
  geom_jitter(aes(color = lineage, fill = lineage),
              position=position_dodge(width=0.3),
              alpha=0.2, pch = 21,
              color = "black") +
  geom_errorbar(data = calc_recent, aes(x = site, ymax = calc+se, ymin = calc-se, color = lineage), width = .2, position = position_dodge(width=0.4)) +
  geom_point(data = calc_recent, mapping = aes(x=site, y=calc, color = lineage, fill = lineage, group=coreID), size = 3.5, pch = 21, color = "black", position = position_dodge(width=0.4))+
  scale_fill_manual(name = "Lineage",
                    breaks = c("1","2"),
                    values = cols_lineage)+
  scale_color_manual(name = "Lineage",
                     breaks = c("1","2"),
                     values = cols_lineage)+
  xlab("Site")+
  ylab("Calcification")+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))+
  theme(legend.position = "blank")
calc_plot_lineage

plots_recent = ggarrange(linext_plot_lineage, den_plot_lineage, calc_plot_lineage,
                         labels = c("A", "B", "C"),
                         ncol = 3, nrow = 1)
plots_recent

ggsave(plots_recent, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/CoringData/coring_plots_individualcores_2005_2014.pdf", width=9, height=4, units=c("in"), useDingbats=FALSE)




## Correlate coring data with temperature
# read in temperature data
panama_temp = read.csv("Coring_Data/data_files/Panama_SST_31Jan24.csv")
str(panama_temp)
panama_temp$param = as.factor(panama_temp$param)

mean_sst = panama_temp %>%
  filter(param == "Annual Mean SST") %>%
  dplyr::rename("AnnualMean" = "SST") %>%
  select(-param)

summer_mean_sst = panama_temp %>%
  filter(param == "Annual Summer Mean") %>%
  dplyr::rename("SummerMean" = "SST") %>%
  select(-param)

# combine with coring data
str(pan_cores_allyears)
str(pan_cores_recent)

pan_cores_allyears$lineage = as.factor(pan_cores_allyears$lineage)
pan_cores_recent$lineage = as.factor(pan_cores_recent$lineage)

pan_cores_allyears_trim = pan_cores_allyears %>%
  select(Year, coreID, linext, density, calc, lineage)

pan_cores_recent_trim = pan_cores_recent %>%
  select(Year, coreID, linext, density, calc, lineage)

pan_cores_allyears_trim2 = left_join(pan_cores_allyears_trim, mean_sst, by = "Year")
pan_cores_allyears_trim3 = left_join(pan_cores_allyears_trim2, summer_mean_sst, by = "Year")

head(pan_cores_allyears_trim3)
View(pan_cores_allyears_trim3)

pan_cores_recent_trim2 = left_join(pan_cores_recent_trim, mean_sst, by = "Year")
pan_cores_recent_trim3 = left_join(pan_cores_recent_trim2, summer_mean_sst, by = "Year")

head(pan_cores_recent_trim3)
View(pan_cores_recent_trim3)

corrplot_annualmean = ggplot(pan_cores_allyears_trim3, aes(x = AnnualMean, y = calc, color = lineage, group = lineage)) +
  theme_bw() +
  geom_point() +
  geom_smooth(aes(fill=lineage)) +
  scale_color_manual(name = "Lineage",
                     breaks = c("1","2"),
                     values = cols_lineage) +
  scale_fill_manual(name = "Lineage",
                    breaks = c("1","2"),
                    values = cols_lineage) +
  xlab("Annual Mean Temperature (°C)") +
  ylab("Calcification")
corrplot_annualmean

corrplot_summermean = ggplot(pan_cores_allyears_trim3, aes(x = SummerMean, y = calc, color = lineage, group = lineage)) +
  theme_bw() +
  geom_point() +
  geom_smooth(aes(fill=lineage)) +
  scale_color_manual(name = "Lineage",
                     breaks = c("1","2"),
                     values = cols_lineage) +
  scale_fill_manual(name = "Lineage",
                    breaks = c("1","2"),
                    values = cols_lineage) +
  xlab("Summer Mean Temperature (°C)") +
  ylab("Calcification")
corrplot_summermean

corrplot_annualmean_recent = ggplot(pan_cores_recent_trim3, aes(x = AnnualMean, y = calc, color = lineage, group = lineage)) +
  theme_bw() +
  geom_point() +
  geom_smooth(aes(fill=lineage)) +
  scale_color_manual(name = "Lineage",
                     breaks = c("1","2"),
                     values = cols_lineage) +
  scale_fill_manual(name = "Lineage",
                    breaks = c("1","2"),
                    values = cols_lineage) +
  xlab("Annual Mean Temperature (°C)") +
  ylab("Calcification")
corrplot_annualmean_recent

corrplot_summermean_recent = ggplot(pan_cores_recent_trim3, aes(x = SummerMean, y = calc, color = lineage, group = lineage)) +
  theme_bw() +
  geom_point() +
  geom_smooth(aes(fill=lineage)) +
  scale_color_manual(name = "Lineage",
                     breaks = c("1","2"),
                     values = cols_lineage) +
  scale_fill_manual(name = "Lineage",
                    breaks = c("1","2"),
                    values = cols_lineage) +
  xlab("Summer Mean Temperature (°C)") +
  ylab("Calcification")
corrplot_summermean_recent


corrplots_temperature = ggarrange(corrplot_annualmean, corrplot_annualmean_recent, corrplot_summermean, corrplot_summermean_recent,
                                  ncol = 2, nrow = 2,
                                  common.legend = TRUE, legend = 'right',
                                  labels = c("A. All Years", "B. 1980-2014", "C. All Years", "D. 1980-2014"))
ggsave(corrplots_temperature, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/CoringData/coring_temperature_correlations.pdf", width=8, height=6, units=c("in"), useDingbats=FALSE)



# plot panama temp increase over time
panama_temp_plot = ggplot(panama_temp, aes(x = Year, y = SST, color = param, group = param)) +
  theme_bw() +
  geom_point() +
  stat_smooth(aes(fill=param),
              method = "lm",
              formula = y ~ x,
              geom = "smooth") +
  scale_color_manual(name = "Parameter",
                     values = c("#998ec3","#f1a340")) +
  scale_fill_manual(name = "Parameter",
                    values = c("#998ec3","#f1a340")) +
  scale_x_continuous(breaks = seq(1870,2020,20))+
  ylab("Temperature (°C)") +
  xlab("Year")
panama_temp_plot
ggsave(panama_temp_plot, filename = "/Users/hannahaichelman/Dropbox/BU/TVE/CoringData/temperature_increase.pdf", width=6, height=4, units=c("in"), useDingbats=FALSE)


# write out file to use in fluorish visualization
pan_cores_allyears_fluorish = pan_cores_allyears_trim3 %>%
  mutate(sitename = str_split_i(coreID, "_", 3))
write.csv(pan_cores_allyears_fluorish, file = "/Users/hannahaichelman/Dropbox/BU_Postdoc/Fluorish_Workshop/coring_temp_data.csv", row.names = FALSE)

head(pan_cores_allyears)

# try again
data_allyears_fluorish = pan_cores_allyears %>%
  select(Year, coreID, calc, site, lineage)
data_allyears_fluorish2 = summarySE(data_allyears_fluorish, measurevar = "calc", groupvars = c("site","Year"))
data_allyears_fluorish3 = summarySE(data_allyears_fluorish, measurevar = "calc", groupvars = c("lineage","Year"))

data_allyears_fluorish_site = data_allyears_fluorish2 %>%
  select(Year, site, calc)

pivoted_data_site = pivot_wider(data_allyears_fluorish_site, names_from = Year, values_from = calc)
write.csv(pivoted_data_site, file = "/Users/hannahaichelman/Dropbox/BU_Postdoc/Fluorish_Workshop/pivoted_data_fluorish.csv", row.names = FALSE)

data_allyears_fluorish_lineage = data_allyears_fluorish3 %>%
  select(Year, lineage, calc)

pivoted_data_lineage = pivot_wider(data_allyears_fluorish_lineage, names_from = Year, values_from = calc)

write.csv(pivoted_data_lineage, file = "/Users/hannahaichelman/Dropbox/BU_Postdoc/Fluorish_Workshop/pivoted_data_lineage_fluorish.csv", row.names = FALSE)

