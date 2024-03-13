#Hannah's TVE 16S analysis
#Almost entirely based on DADA2 Pipeline 1.8 Walkthrough:
#https://benjjneb.github.io/dada2/tutorial.html
#with edits by Carly D. Kenkel and modifications by Nicola Kriefall

# Skip to line 310 to avoid stuff that happened on the BU cluster.
#~########################~#
##### PRE-PROCESSING #######
#~########################~#

#fastq files should have R1 & R2 designations for PE reads
#Also - some pre-trimming. Retain only PE reads that match amplicon primer. Remove reads containing Illumina sequencing adapters

#in Terminal home directory:
#following instructions of installing BBtools from https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/installation-guide/
#1. download BBMap package, sftp to installation directory
#2. untar:
#tar -xvzf BBMap_(version).tar.gz
#3. test package:
#cd bbmap
#~/bin/bbmap/stats.sh in=~/bin/bbmap/resources/phix174_ill.ref.fa.gz

# my adaptors for 16S, which I saved as "adaptors.fasta"
# >forward
# AATGATACGGCGACCAC
# >forwardrc
# GTGGTCGCCGTATCATT
# >reverse
# CAAGCAGAAGACGGCATAC
# >reverserc
# GTATGCCGTCTTCTGCTTG

#primers for 16S:
# >forward
# GTGYCAGCMGCCGCGGTA
# >reverse
# GGACTACHVGGGTWTCTAAT

##Still in terminal - making a sample list based on the first phrase before the underscore in the .fastq name
#ls *R1_001.fastq | cut -d '_' -f 1 > samples.list

##cuts off the extra words in the .fastq files
#for file in $(cat samples.list); do  mv ${file}_*R1*.fastq ${file}_R1.fastq; mv ${file}_*R2*.fastq ${file}_R2.fastq; done

##gets rid of reads that still have the adaptor sequence, shouldn't be there, I didn't have any
#for file in $(cat samples.list); do ~/bin/bbmap/bbduk.sh in1=${file}_R1.fastq in2=${file}_R2.fastq ref=adaptors.fasta out1=${file}_R1_NoIll.fastq out2=${file}_R2_NoIll.fastq; done &>bbduk_NoIll.log

##getting rid of first 4 bases (degenerate primers created them)
#for file in $(cat samples.list); do ~/bin/bbmap/bbduk.sh in1=${file}_R1_NoIll.fastq in2=${file}_R2_NoIll.fastq ftl=4 out1=${file}_R1_NoIll_No4N.fastq out2=${file}_R2_NoIll_No4N.fastq; done &>bbduk_No4N.log

##only keeping reads that start with the 16S primer
#for file in $(cat samples.list); do ~/bin/bbmap/bbduk.sh in1=${file}_R1_NoIll_No4N.fastq in2=${file}_R2_NoIll_No4N.fastq restrictleft=20 k=10 literal=GTGYCAGCMGCCGCGGTA,GGACTACHVGGGTWTCTAAT copyundefined=t outm1=${file}_R1_NoIll_No4N_16S.fastq outu1=${file}_R1_check.fastq outm2=${file}_R2_NoIll_No4N_16S.fastq outu2=${file}_R2_check.fastq; done &>bbduk_16S.log
##higher k = more reads removed, but can't surpass k=20 or 21

##using cutadapt to remove primer
# for file in $(cat samples.list)
# do
# cutadapt -g GTGYCAGCMGCCGCGGTA -a ATTAGAWACCCVHGTAGTCC -G GGACTACHVGGGTWTCTAAT -A TACCGCGGCKGCTGRCAC -n 2 --discard-untrimmed -o ${file}_R1.fastq -p ${file}_R2.fastq ${file}_R1_NoIll_No4N_16S.fastq ${file}_R2_NoIll_No4N_16S.fastq
# done &> clip.log
##-g regular 5' forward primer
##-G regular 5' reverse primer
##-o forward out
##-p reverse out
##-max-n 0 means 0 Ns allowed
##this overwrote my original renamed files

# did sftp of *_R1.fastq & *_R2.fastq files to the folder to be used in dada2

# the 16s files for the pre-stress timepoint are here on scc: /projectnb/davies-hb/hannah/TVE_16S_ITS/tve_prestress_files/lane1and2_16S

#~########################~#
##### DADA2 BEGINS #########
#~########################~#

# ran into issues downstream, need to install & load more recent version of dada2:
# library(devtools)
# devtools::install_github("benjjneb/dada2")

library(dada2); packageVersion("dada2")
#Version 1.29.0
library(ShortRead); packageVersion("ShortRead")
#1.56.0
library(Biostrings); packageVersion("Biostrings")
#2.66.0

path <- "/Users/hannahaichelman/Dropbox/BU/TVE/16S_ITS2/16S_All_Timepoints/fastqs_t0" # CHANGE ME to the directory containing the fastq files after unzipping.

fnFs <- sort(list.files(path, pattern = "R1_16S_final.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "R2_16S_final.fastq", full.names = TRUE))

get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(fnFs, get.sample.name))
head(sample.names)
sample.names

#### check for primers ####
FWD <- "GTGYCAGCMGCCGCGGTA"  ## CHANGE ME to your forward primer sequence
REV <- "GGACTACHVGGGTWTCTAAT"  ## CHANGE ME...

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna),
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients
REV.orients

fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[9]]),
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[9]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[9]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[9]]))
#spot checked a few files - don't see any primers

#### Visualizing raw data ####

#First, lets look at quality profile of R1 reads
plotQualityProfile(fnFs.filtN[c(1,2,3,4)])
plotQualityProfile(fnFs.filtN[c(50,51,52,53)])
#looks mostly good up to 180 I think

#Then look at quality profile of R2 reads
plotQualityProfile(fnRs.filtN[c(1,2,3,4)])
plotQualityProfile(fnRs.filtN[c(50,51,52,53)])
#180 again

# Make directory and filenames for the filtered fastqs
filt_path <- file.path(path, "trimmed")
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

#changing a bit from default settings - maxEE=1 (1 max expected error, more conservative), truncating length at 200 bp for both forward & reverse [leaves ~50bp overlap], added "trimleft" to cut off primers [18 for forward, 20 for reverse]
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     truncLen=c(175,175), #leaves ~50bp overlap
                     maxN=0, #DADA does not allow Ns
                     maxEE=c(1,1), #allow 1 expected errors, where EE = sum(10^(-Q/10)); more conservative, model converges
                     truncQ=2,
                     #trimLeft=c(18,20), #N nucleotides to remove from the start of each read
                     rm.phix=TRUE, #remove reads matching phiX genome
                     matchIDs=TRUE, #enforce matching between id-line sequence identifiers of F and R reads
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE

head(out)
tail(out)

# Had to make the trimmed files on the cluster when i ran it with all files (prestress and T0) because of lack of space on my local machine
# Did this here:
# "/projectnb/davies-hb/hannah/TVE_16S_ITS/tve_prestress_files/lane1and2_16S/fastqs"
# by qsub'ing the trimmed_files_R file
# was able to do this on my local machine for the T0 files

#~############################~#
##### Learn Error Rates ########
#~############################~#

#setDadaOpt(MAX_CONSIST=30) #increase number of cycles to allow convergence
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

#sanity check: visualize estimated error rates
#error rates should decline with increasing qual score
#red line is based on definition of quality score alone
#black line is estimated error rate after convergence
#dots are observed error rate for each quality score

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

#~############################~#
##### Dereplicate reads ########
#~############################~#
#Dereplication combines all identical sequencing reads into “unique sequences” with a corresponding “abundance”: the number of reads with that unique sequence.
#Dereplication substantially reduces computation time by eliminating redundant comparisons.
#DADA2 retains a summary of the quality information associated with each unique sequence. The consensus quality profile of a unique sequence is the average of the positional qualities from the dereplicated reads. These quality profiles inform the error model of the subsequent denoising step, significantly increasing DADA2’s accuracy.
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#~###############################~#
##### Infer Sequence Variants #####
#~###############################~#

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

#now, look at the dada class objects by sample
#will tell how many 'real' variants in unique input seqs
#By default, the dada function processes each sample independently, but pooled processing is available with pool=TRUE and that may give better results for low sampling depths at the cost of increased computation time. See our discussion about pooling samples for sample inference.
dadaFs[[1]]
dadaRs[[1]]

#~############################~#
##### Merge paired reads #######
#~############################~#

#To further cull spurious sequence variants
#Merge the denoised forward and reverse reads
#Paired reads that do not exactly overlap are removed

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

summary((mergers[[1]]))

#We now have a data.frame for each sample with the merged $sequence, its $abundance, and the indices of the merged $forward and $reverse denoised sequences. Paired reads that did not exactly overlap were removed by mergePairs.

#~##################################~#
##### Construct sequence table #######
#~##################################~#
#a higher-resolution version of the “OTU table” produced by classical methods

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
#[1]    54 10735

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

plot(table(nchar(getSequences(seqtab)))) #real variants appear to be right in that 244-264 window

#The sequence table is a matrix with rows corresponding to (and named by) the samples, and
#columns corresponding to (and named by) the sequence variants.
#Sequences that are much longer or shorter than expected may be the result of non-specific priming, and may be worth removing

# trying to figure out what these two peaks are, make seq tables of both peaks
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(240,260)] #again, being fairly conservative wrt length

#~############################~#
##### Remove chimeras ##########
#~############################~#
#The core dada method removes substitution and indel errors, but chimeras remain.
#Fortunately, the accuracy of the sequences after denoising makes identifying chimeras easier
#than it is when dealing with fuzzy OTUs: all sequences which can be exactly reconstructed as
#a bimera (two-parent chimera) from more abundant sequences.

seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
#Identified 1675 bimeras out of 10619 input sequences.

#The fraction of chimeras varies based on factors including experimental procedures and sample complexity,
#but can be substantial.
sum(seqtab.nochim)/sum(seqtab2)
#0.952486

saveRDS(seqtab.nochim, file="/Users/hannahaichelman/Dropbox/BU/TVE/16S_ITS2/16S_All_Timepoints/tve16s_all_seqtab.nochim.rds")
write.csv(seqtab.nochim, file="/Users/hannahaichelman/Dropbox/BU/TVE/16S_ITS2/16S_All_Timepoints/tve16s_all_seqtab.nochim.csv")

#~############################~#
##### Track Read Stats #########
#~############################~#

# note that because I created the trimmed files on the cluster, I can't make this file here. But the older version is still relevant (just has the lineage 3 individuals included)
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab2), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)
tail(track)

write.csv(track,file="/Users/hannahaichelman/Dropbox/BU/TVE/16S_ITS2/16S_All_Timepoints/tve16s_t0_readstats.csv",row.names=TRUE,quote=FALSE)

#~############################~#
##### Assign Taxonomy ##########
#~############################~#

#Assign Taxonomy
# downloaded most recent silva files from here: https://zenodo.org/record/4587955#.Yd8ZRxPMJmo

# scp seqtab.nochim.rds file to scc

# this step had to happen on the SCC because it took too much memory to run locally.
# on scc: /projectnb/davies-hb/hannah/TVE_16S_ITS/tve_prestress_files/lane1and2_16S/assign_tax
# or for T0 analysis: /projectnb/davies-hb/hannah/TVE_Panama/TVE_16S_ITS/tve_16s_allfiles/assign_tax

# first, module load R and type R to open interactive window.
# next, install dada2 using biocmanager
# then exit the R interactive environment and submit the job "assign_tax_R" using qsub.

# scp the .rds files back to local machine when finished running

# here is how I would do it locally if it didn't take too much memory:
taxa <- assignTaxonomy(seqtab.nochim, "/Users/hannahaichelman/Documents/BU/TVE/16S_ITS2/16S_PreStress_Timepoint/silva_nr_v132_train_set.fa.gz",tryRC=TRUE)
unname(head(taxa))
taxa.plus <- addSpecies(taxa, "/Users/hannahaichelman/Documents/BU/TVE/16S_ITS2/16S_PreStress_Timepoint/silva_species_assignment_v132.fa.gz",tryRC=TRUE,verbose=TRUE)
# 247 out of 2608 were assigned to the species level.
# Of which 223 had genera consistent with the input table.

saveRDS(taxa.plus, file="/Users/hannahaichelman/Documents/BU/TVE/16S_ITS2/16S_PreStress_Timepoint/tve16s_t0_taxaplus.rds")
saveRDS(taxa, file="/Users/hannahaichelman/Documents/BU/TVE/16S_ITS2/16S_PreStress_Timepoint/tve16s_t0_taxa.rds")
#write.csv(taxa.plus, file="mr16s_taxaplus.csv")
#write.csv(taxa, file="mr16s_taxa.csv")

#### Read in previously saved datafiles  - START HERE ####
# Trying now with all samples included, which was previously done following the same steps above just with all samples included
setwd("~/Dropbox/BU/TVE/TVE_Github/DielTempVariability")

seqtab.nochim <- readRDS("16S_Microbiome/data_files/tve16s_all_seqtab.nochim_newids.rds")
taxa <- readRDS("16S_Microbiome/data_files/tve16s_all_seq_taxa_newids.rds")
#taxa.plus <- readRDS("mr16s_revised_taxaplus.rds")

# check for samples with 0 reads - all looks good for T0, but I2C4 has 0 reads
rowSums(seqtab.nochim)

row_names_to_remove<-c("I2C4") # I2C4 being removed because of 0 reads, ITS because they were negative controls specific to ITS
seqtab.nochim <- seqtab.nochim[!(row.names(seqtab.nochim) %in% row_names_to_remove),]

#~############################~#
##### handoff 2 phyloseq #######
#~############################~#

#BiocManager::install("phyloseq")
library('phyloseq')
library('ggplot2')
library('Rmisc')
library('cowplot')
library('ShortRead')
library('tidyverse')

#import dataframe holding sample information
samdf = read.csv("16S_Microbiome/data_files/SampleInfo.csv")
dim(samdf)
head(samdf)

# add identifying data that includes lineage 3 - best to include from ITS2 metadata
phys_metadata = read.csv("ITS2_Symbiodiniaceae/data_files/phys_metadata_its2.csv")
head(phys_metadata)
dim(phys_metadata)

# combine with samdf
samdf2 = left_join(samdf, phys_metadata, by = "frag")
head(samdf2)
dim(samdf2)

samdf_t0 = samdf2 %>%
  dplyr::filter(treat == "Initial") %>%
  select(-dominant_type)
head(samdf_t0)
dim(samdf_t0)

samdf_all = samdf2 %>%
  select(-dominant_type) %>%
  filter(frag!="I2C4")
head(samdf_all)
dim(samdf_all)

rownames(samdf_all) <- samdf_all$frag

# Construct phyloseq object (straightforward from dada2 outputs)
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),
               sample_data(samdf_all),
               tax_table(taxa))

ps
# 8944 taxa and 54 samples for T0
# for all samples:
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 17903 taxa and 239 samples ]
# sample_data() Sample Data:       [ 239 samples by 8 sample variables ]
# tax_table()   Taxonomy Table:    [ 17903 taxa by 8 taxonomic ranks ]


#### first look at data ####
ps_glom <- tax_glom(ps, "Phylum")
plot_bar(ps_glom, x="lineage", fill="Phylum")+
  theme(legend.position="none")

#phyloseq object with shorter names - doing this one instead of one above
ids <- paste0("sq", seq(1, length(colnames(seqtab.nochim))))

#making output fasta file for lulu step & maybe other things
library(dada2)
path='16S_Microbiome/data_files/tve16s_all.fasta'
uniquesToFasta(seqtab.nochim, path, mode = "w", width = 20000)

colnames(seqtab.nochim)<-ids
taxa2 <- cbind(taxa, rownames(taxa)) #retaining raw sequence info before renaming
rownames(taxa2)<-ids

#phyloseq object with new taxa ids
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),
               sample_data(samdf_all),
               tax_table(taxa2))

ps #17903 taxa and 239 samples

#save(taxa2,file="taxa2.Rdata")
#save(taxa,file="taxa.Rdata")

#### remove mitochondria, chloroplasts, non-bacteria ####
ps.mito <- subset_taxa(ps, (Family=="Mitochondria"))
ps.mito #463 taxa to remove
ps.chlor <- subset_taxa(ps, (Order=="Chloroplast"))
ps.chlor #816 taxa to remove
ps.notbact <- subset_taxa(ps, (Kingdom!="Bacteria") | is.na(Kingdom))
ps.notbact #208 taxa to remove

ps.nomito <- subset_taxa(ps, (Family!="Mitochondria") | is.na(Family))
ps.nomito #17440 taxa
ps.nochlor <- subset_taxa(ps.nomito, (Order!="Chloroplast") | is.na(Order))
ps.nochlor #16624 taxa
ps.clean <- subset_taxa(ps.nochlor, (Kingdom=="Bacteria"))
ps.clean #16416 taxa

#just archaea
ps.arch <- subset_taxa(ps.nomito, (Kingdom=="Archaea"))
ps.arch #204 taxa

#### identifying contamination ####
# using negative control from prestress sequencing to remove contamination in all samples
#install.packages("decontam")
library(decontam)

df <- as.data.frame(sample_data(ps.clean)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps.clean)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=sitename)) + geom_point()

sample_data(ps.clean)$is.neg <- sample_data(ps.clean)$Sample_or_Control == "Control"
contamdf.prev <- isContaminant(ps.clean, neg="is.neg",threshold=0.5)
table(contamdf.prev$contaminant)
# FALSE  TRUE
# 16256   160

# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(ps.clean, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Control", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Control", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

#remove from ps.clean:
ps.clean1 <- prune_taxa(!contamdf.prev$contaminant,ps.clean)
#also remove negative controls, don't need them anymore I think
ps.cleaner <- subset_samples(ps.clean1,(Sample_or_Control!="Control"))

#### blast asvs to NCBI to see if any eukaryotes got through ####
##Running blast on BU SCC to make organism match files for my 16s data
##used 'tve16s_t0.fasta' made above
## Working here on cluster:
# /projectnb/davies-hb/hannah/TVE_Panama/TVE_16S_ITS/tve_16s_allfiles/euk_blast/blast_attempt

# have to split fasta file up to make it run faster
#awk -v size=5000 -v pre=tve16s_t0_split -v pad=5 '/^>/{n++;if(n%size==1){close(f);f=sprintf("%s.%0"pad"d",pre,n)}}{print>>f}' tve16s_t0.fasta

#module load blast+
##submitted the following array job: qsub split_blast_array.qsub
#made that qsub file by doing this, make sure to add -pe omp 10:
#scc6_qsub_launcher.py -N split_blast_t0 -P coral -M haich@bu.edu -j y -h_rt 64:00:00 -jobsfile split_blast_t0

# [haich@scc1 euk_blast]$ cat split_blast_t0
# blastn -query tve16s_t0_split.00001 -db nt -outfmt "6 std staxids sskingdoms" -evalue 1e-5 -max_target_seqs 5 -out tve16s_t0_split.00001_taxids.out -remote
# blastn -query tve16s_t0_split.05001 -db nt -outfmt "6 std staxids sskingdoms" -evalue 1e-5 -max_target_seqs 5 -out tve16s_t0_split.05001_taxids.out -remote


##takes a very long time (Nicola had ~2600 ASVs, took 11 hours)
# when job finishes, concatenate the split output files into one
# [haich@scc1 euk_blast]$ cat tve16s_prestress_split.00001_taxids.out tve16s_prestress_split.05001_taxids.out > tve16s_prestress_split_taxids.out

##now getting taxonomy info:

# #download/install taxonkit things, more instructions here:
# #https://bioinf.shenwei.me/taxonkit/usage/

# cd /net/scc-pa2/scratch/
# wget -c ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
# tar -zxvf taxdump.tar.gz
# cd
# [haich@scc1 taxa]$ cp *.dmp /usr3/graduate/haich/.taxonkit

# module load miniconda
# conda install -c bioconda taxonkit -p .
# cd /net/scc-pa2/scratch/taxa/
# cp *.dmp ~/.taxonkit
# #command taxonkit should work now

##extracting taxa ids from blast output for taxonkit:
# awk -F " " '{print $13}' tve16s_prestress_split_taxids.out > ids
# taxonkit lineage ids > ids.tax
# cut -f1 tve16s_prestress_split_taxids.out > ids.seq; paste ids.seq ids.tax > ids.seq.tax
# grep "Eukaryota" ids.seq.tax | cut -f1 | sort | uniq > euk.contam.asvs

##transferring euk.contam.asvs to back here
##remove from ps.cleaner
##should be 151 to remove
euks <- read.csv("16S_Microbiome/data_files/euk.contam.asvs.csv",header=FALSE)
euks_names <- euks$V1
alltaxa <- taxa_names(ps.cleaner)
keepers <- alltaxa[(!alltaxa %in% euks_names)] #doesn't look like any were removed
ps.cleanest <- prune_taxa(keepers, ps.cleaner)
#16256 taxa and 237 samples

seqtab.cleanest <- data.frame(otu_table(ps.cleanest))
#write.csv(seqtab.cleanest,file="tve16s_seqtab.rev.cleanest.csv")

##save cleaned phyloseq object
saveRDS(ps.cleanest,file="phyloseq.cleanest.all.rds")

#### Decontaminated (Euk contamination removed) files ####
setwd("~/Dropbox/BU/TVE/TVE_Github/DielTempVariability")
ps.cleanest = readRDS("16S_Microbiome/data_files/phyloseq.cleanest.all.rds")
seqtab.cleanest <- data.frame(ps.cleanest@otu_table)
samdf.cleanest <- data.frame(ps.cleanest@sam_data)

# subset T0 and pre-stress files here
ps.cleanest.t0 = subset_samples(ps.cleanest, timepoint=="t0")
seqtab.t0.cleanest <- data.frame(ps.cleanest.t0@otu_table)
samdf.t0.cleanest <- data.frame(ps.cleanest.t0@sam_data)
saveRDS(ps.cleanest.t0,file="phyloseq.cleanest.t0.rds")

ps.cleanest.prestress = subset_samples(ps.cleanest, timepoint=="pre_stress")
seqtab.prestress.cleanest <- data.frame(ps.cleanest.prestress@otu_table)
samdf.prestress.cleanest <- data.frame(ps.cleanest.prestress@sam_data)
saveRDS(ps.cleanest.prestress,file="phyloseq.cleanest.prestress.rds")

#### rarefy decontaminated data #####
library(vegan)
load("16S_Microbiome/data_files/taxa2.Rdata")

# more info on rarefying: https://micca.readthedocs.io/en/latest/phyloseq.html
# plot rarefaction curve
rarecurve(seqtab.t0.cleanest,step=100,label=TRUE) #after removing contaminants
rarecurve(seqtab.prestress.cleanest,step=100,label=TRUE) #after removing contaminants

# Plot reads per sample - all separately for the two timepoints
df.t0 = data.frame(ASVs=rowSums(otu_table(ps.cleanest.t0)>0), reads=sample_sums(ps.cleanest.t0), sample_data(ps.cleanest.t0))
df.prestress = data.frame(ASVs=rowSums(otu_table(ps.cleanest.prestress)>0), reads=sample_sums(ps.cleanest.prestress), sample_data(ps.cleanest.prestress))

ggplot(df.prestress, aes(x=reads)) +
  geom_histogram(bins=50, color='black', fill='grey') +
  theme_bw() +
  geom_vline(xintercept=1000, color= "red", linetype='dashed') +
  labs(title="Histogram: Reads per Sample") + xlab("Read Count") + ylab("Sample Count")

total.t0 <- rowSums(seqtab.t0.cleanest)
total.prestress <- rowSums(seqtab.prestress.cleanest)

min(total.t0)
min(total.prestress)
subset(total.t0, total.t0 <1000)
subset(total.prestress, total.prestress <1000)

# lose 35 samples at 1000 seqs for prestress
# lose 0 samples at 1000 seqs for T0
# Justification for 1000 seq cut-off for rarefying: https://www.nature.com/articles/s41467-018-07275-x

row.names.remove.t0 <- names(subset(total.t0, total.t0 <1000)) # no samples to remove here
row.names.remove.prestress <- names(subset(total.prestress, total.prestress <1000))

seqtab.less.t0 <- seqtab.t0.cleanest[!(row.names(seqtab.t0.cleanest) %in% row.names.remove.t0),]
seqtab.less.prestress <- seqtab.prestress.cleanest[!(row.names(seqtab.prestress.cleanest) %in% row.names.remove.prestress),]

samdf.rare.t0 <- samdf.t0.cleanest[!(row.names(samdf.t0.cleanest) %in% row.names.remove.t0), ]
samdf.rare.prestress <- samdf.prestress.cleanest[!(row.names(samdf.prestress.cleanest) %in% row.names.remove.prestress), ]

# rarefy to 1000 reads per sample
seqtab.rare.t0 <- rrarefy(seqtab.less.t0,sample=1000)
rarecurve(seqtab.rare.t0,step=100,label=TRUE)

seqtab.rare.prestress <- rrarefy(seqtab.less.prestress,sample=1000)
rarecurve(seqtab.rare.prestress,step=100,label=TRUE)

#phyloseq object but rarefied
ps.rare.t0 <- phyloseq(otu_table(seqtab.rare.t0, taxa_are_rows=FALSE),
               sample_data(samdf.rare.t0),
               tax_table(taxa))
ps.rare.t0 #16256 taxa and 54 samples

ps.rare.prestress <- phyloseq(otu_table(seqtab.rare.prestress, taxa_are_rows=FALSE),
                    sample_data(samdf.rare.prestress),
                    tax_table(taxa))
ps.rare.prestress #16256 taxa and 148 samples

#removing missing taxa - lost after rarefying
ps.rare.t0 <- prune_taxa(taxa_sums(ps.rare.t0) > 0, ps.rare.t0)
ps.rare.t0 #6111 taxa and 54 samples with 1000 rarefy

ps.rare.prestress <- prune_taxa(taxa_sums(ps.rare.prestress) > 0, ps.rare.prestress)
ps.rare.prestress #5852 taxa and 148 samples with 1000 rarefy

seqtab.rare.t0 <- data.frame(otu_table(ps.rare.t0))
seqtab.rare.prestress <- data.frame(otu_table(ps.rare.prestress))

#saving
#saveRDS(ps.rare.t0,file="phyloseq.t0.rarefied.1k.rds")
#saveRDS(ps.rare.prestress,file="phyloseq.prestress.rarefied.1k.rds")

#write.csv(seqtab.rare.t0, file="tve16s_seqtab.t0.rev.cleanest.rare_1k")
#write.csv(seqtab.rare.prestress, file="tve16s_seqtab.prestress.rev.cleanest.rare_1k")

### data files - decontaminated, rarefied ####

setwd("~/Dropbox/BU/TVE/TVE_Github/DielTempVariability")
ps.rare.1k.t0 = readRDS("16S_Microbiome/data_files/phyloseq.t0.rarefied.1k.rds")
seqtab.rare.1k.t0 <- data.frame(ps.rare.1k.t0@otu_table)
samdf.rare.1k.t0 <- data.frame(ps.rare.1k.t0@sam_data)

ps.rare.1k.prestress = readRDS("16S_Microbiome/data_files/phyloseq.prestress.rarefied.1k.rds")
seqtab.rare.1k.prestress <- data.frame(ps.rare.1k.prestress@otu_table)
samdf.rare.1k.prestress <- data.frame(ps.rare.1k.prestress@sam_data)

load("16S_Microbiome/data_files/taxa2.Rdata")

ps.rare <- phyloseq(otu_table(seqtab.rare.1k, taxa_are_rows=FALSE),
                    sample_data(samdf.rare.1k),
                    tax_table(taxa2))
ps.rare # 11369 taxa and 202 samples

#### trim underrepresented otus ####
# don't use rarefied data for this
library(MCMC.OTU)

#formatting the table for mcmc.otu - requires one first column that's 1 through whatever
#& has "X" as column name
nums.t0 <- 1:nrow(seqtab.t0.cleanest)
samples.t0 <- rownames(seqtab.t0.cleanest)

nums.prestress <- 1:nrow(seqtab.prestress.cleanest)
samples.prestress <- rownames(seqtab.prestress.cleanest)

int.t0 <- cbind(sample = 0, seqtab.t0.cleanest)
seq.formcmc.t0 <- cbind(X = 0, int.t0)

int.prestress <- cbind(sample = 0, seqtab.prestress.cleanest)
seq.formcmc.prestress <- cbind(X = 0, int.prestress)

seq.formcmc.t0$X <- nums.t0
seq.formcmc.t0$sample <- samples.t0

seq.formcmc.prestress$X <- nums.prestress
seq.formcmc.prestress$sample <- samples.prestress

seq.trim.allinfo.t0 <- purgeOutliers(seq.formcmc.t0,count.columns=3:16258,sampleZcut=-2.5,otu.cut=0.0001,zero.cut=0.02)
# [1] "samples with counts below z-score -2.5 :"
# [1] "I2A" "I3A"
# [1] "zscores:"
# I2A       I3A
# -3.113507 -3.854946
# [1] "OTUs passing frequency cutoff  1e-04 : 1505"
# [1] "OTUs with counts in 0.02 of samples:"
#
# FALSE  TRUE
# 338  1167

seq.trim.allinfo.prestress <- purgeOutliers(seq.formcmc.prestress,count.columns=3:16258,sampleZcut=-2.5,otu.cut=0.0001,zero.cut=0.02)
# [1] "samples with counts below z-score -2.5 :"
# [1] "I2F1"
# [1] "zscores:"
# I2F1
# -3.941567
# [1] "OTUs passing frequency cutoff  1e-04 : 1038"
# [1] "OTUs with counts in 0.02 of samples:"
#
# FALSE  TRUE
# 295   743

#remove sample info
seq.trim.t0 <- seq.trim.allinfo.t0[,3:1169]
seq.trim.prestress <- seq.trim.allinfo.prestress[,3:745]

#write.csv(seq.trim.t0,file="tve16s_seqtab.rev.cleanest.trim.t0.csv")
#write.csv(seq.trim.prestress,file="tve16s_seqtab.rev.cleanest.trim.prestress.csv")

#remake phyloseq objects
ps.trim.t0 <- phyloseq(otu_table(seq.trim.t0, taxa_are_rows=FALSE),
                         sample_data(samdf.t0.cleanest),
                         tax_table(taxa2))
ps.trim.t0 #1167 taxa and 52 samples

#saveRDS(ps.trim,file="phyloseq.cleanest.trim.rds")

ps.trim.prestress <- phyloseq(otu_table(seq.trim.prestress, taxa_are_rows=FALSE),
                    sample_data(samdf.prestress.cleanest),
                    tax_table(taxa2))
ps.trim.prestress #743 taxa and 182 samples

#saveRDS(ps.trim.t0,file="phyloseq.cleanest.trim.t0.rds")
#saveRDS(ps.trim.prestress,file="phyloseq.cleanest.trim.prestress.rds")


#### rarefy trimmed data #####
library(vegan)

ps.trim.t0 = readRDS("16S_Microbiome/data_files/phyloseq.cleanest.trim.t0.rds")
seqtab.trim.t0 <- data.frame(ps.trim.t0@otu_table)
samdf.trim.t0 <- data.frame(ps.trim.t0@sam_data)

ps.trim.prestress = readRDS("16S_Microbiome/data_files/phyloseq.cleanest.trim.prestress.rds")
seqtab.trim.prestress <- data.frame(ps.trim.prestress@otu_table)
samdf.trim.prestress <- data.frame(ps.trim.prestress@sam_data)

# plot rarefaction curve
rarecurve(seqtab.trim.t0,step=100,label=TRUE)
rarecurve(seqtab.trim.prestress,step=100,label=TRUE)

# establish samples to remove
total.t0 <- rowSums(seqtab.trim.t0)
total.prestress <- rowSums(seqtab.trim.prestress)

min(total.t0)
min(total.prestress)
subset(total.t0, total.t0 <1000)
subset(total.prestress, total.prestress <1000)

# lose 42 samples at 1000 seqs for prestress
# lose 0 samples at 1000 seqs for T0
# Justification for 1000 seq cut-off for rarefying: https://www.nature.com/articles/s41467-018-07275-x

row.names.remove.t0 <- names(subset(total.t0, total.t0 <1000)) # no samples to remove here
row.names.remove.prestress <- names(subset(total.prestress, total.prestress <1000))

seqtab.less.t0 <- seqtab.trim.t0[!(row.names(seqtab.trim.t0) %in% row.names.remove.t0),]
seqtab.less.prestress <- seqtab.trim.prestress[!(row.names(seqtab.trim.prestress) %in% row.names.remove.prestress),]

samdf.rare.t0 <- samdf.trim.t0[!(row.names(samdf.trim.t0) %in% row.names.remove.t0), ]
samdf.rare.prestress <- samdf.trim.prestress[!(row.names(samdf.trim.prestress) %in% row.names.remove.prestress), ]

# rarefy to 1000 reads per sample
seqtab.rare.t0 <- rrarefy(seqtab.less.t0,sample=1000)
rarecurve(seqtab.rare.t0,step=100,label=TRUE)

seqtab.rare.prestress <- rrarefy(seqtab.less.prestress,sample=1000)
rarecurve(seqtab.rare.prestress,step=100,label=TRUE)

#phyloseq object but rarefied
ps.trim.rare.t0 <- phyloseq(otu_table(seqtab.rare.t0, taxa_are_rows=FALSE),
                       sample_data(samdf.rare.t0),
                       tax_table(taxa))
ps.trim.rare.t0 #1167 taxa and 52 samples

ps.trim.rare.prestress <- phyloseq(otu_table(seqtab.rare.prestress, taxa_are_rows=FALSE),
                              sample_data(samdf.rare.prestress),
                              tax_table(taxa))
ps.trim.rare.prestress #743 taxa and 140 samples

#removing missing taxa - lost after rarefying
ps.trim.rare.t0 <- prune_taxa(taxa_sums(ps.trim.rare.t0) > 0, ps.trim.rare.t0)
ps.trim.rare.t0 #1167 taxa and 52 samples

ps.trim.rare.prestress <- prune_taxa(taxa_sums(ps.trim.rare.prestress) > 0, ps.trim.rare.prestress)
ps.trim.rare.prestress #743 taxa and 140 samples

seqtab.trim.rare.t0 <- data.frame(otu_table(ps.trim.rare.t0))
seqtab.trim.rare.prestress <- data.frame(otu_table(ps.trim.rare.prestress))

#saving
saveRDS(ps.trim.rare.t0,file="phyloseq.t0.trim.rarefied.1k.rds")
saveRDS(ps.trim.rare.prestress,file="phyloseq.prestress.trim.rarefied.1k.rds")

#write.csv(seqtab.rare.t0, file="tve16s_seqtab.t0.rev.cleanest.rare_1k")
#write.csv(seqtab.rare.prestress, file="tve16s_seqtab.prestress.rev.cleanest.rare_1k")

### data files - decontaminated, trimmed, rarefied ####

setwd("~/Dropbox/BU/TVE/TVE_Github/DielTempVariability")
ps.trim.rare.1k.t0 = readRDS("16S_Microbiome/data_files/phyloseq.t0.trim.rarefied.1k.rds")
seqtab.trim.rare.1k.t0 <- data.frame(ps.trim.rare.1k.t0@otu_table)
samdf.trim.rare.1k.t0 <- data.frame(ps.trim.rare.1k.t0@sam_data)

ps.trim.rare.1k.prestress = readRDS("16S_Microbiome/data_files/phyloseq.prestress.trim.rarefied.1k.rds")
seqtab.trim.rare.1k.prestress <- data.frame(ps.trim.rare.1k.prestress@otu_table)
samdf.trim.rare.1k.prestress <- data.frame(ps.trim.rare.1k.prestress@sam_data)

load("16S_Microbiome/data_files/taxa2.Rdata")

ps.trim.rare.t0 <- phyloseq(otu_table(seqtab.trim.rare.1k.t0, taxa_are_rows=FALSE),
                    sample_data(samdf.trim.rare.1k.t0),
                    tax_table(taxa2))
ps.trim.rare.t0 # 11369 taxa and 202 samples

ps.trim.rare.prestress <- phyloseq(otu_table(seqtab.trim.rare.1k.prestress, taxa_are_rows=FALSE),
                    sample_data(samdf.trim.rare.1k.prestress),
                    tax_table(taxa2))
ps.trim.rare.prestress # 11369 taxa and 202 samples

#### making fasta file for picrust2 - trimmed not rarefied - SKIPPING BC NOT RUNNING PICRUST ####
library(phyloseq)
library(dada2)

#if needed:
setwd("~/Dropbox/BU/TVE/TVE_Github/DielTempVariability")

ps.trim = readRDS("phyloseq.cleanest.trim.rds")
seqtab.trim <- data.frame(ps.trim@otu_table)
samdf.trim <- data.frame(ps.trim@sam_data)
load("taxa2.Rdata")

ps.trim <- phyloseq(otu_table(seqtab.trim, taxa_are_rows=FALSE),
                         sample_data(samdf.trim),
                         tax_table(taxa2))
ps.trim #641 taxa and 171 samples

trim.otu <- as.matrix(ps.trim@otu_table)
trim.taxa <- data.frame(ps.trim@tax_table)
rownames(trim.taxa)==colnames(trim.otu)

colnames(trim.otu) <- trim.taxa$V8
ids <- rownames(trim.taxa)

path="/Users/hannahaichelman/Documents/BU/TVE/16S_ITS2/16S_PreStress/tve16s_rev.cleanest.trimmed.fasta"
uniquesToFasta(trim.otu, path, ids = ids, mode = "w", width = 20000)

#re-formatting seq table so picrust likes it:
#a tab-delimited table with ASV ids as the first column and sample abundances as all subsequent columns
seqtab.trim.t <- t(seqtab.trim)
write.table(seqtab.trim.t,file="tve16s_seqtab.cleanest.trim.t.txt")
#manually removed the quotation marks that appeared in the file, and converted to tab delimited file from Excel

#### moving on to tve16s_diversity_analysis.R script in other folder ####

#### Try Batch Correction ####
## This step not necessary when including only the Pre-Stress data.
## This was useful when trying to analyze T0 and Pre-Stress data together.

cran.packages <- c('knitr', 'xtable', 'ggplot2', 'vegan', 'cluster',
                  'gridExtra', 'pheatmap', 'ruv', 'lmerTest', 'bapred')
install.packages(cran.packages)
bioconductor.packages <- c('sva', 'limma', 'AgiMicroRna',
                          'variancePartition', 'pvca')
if (!requireNamespace('BiocManager', quietly = TRUE))
    install.packages('BiocManager')
BiocManager::install(bioconductor.packages)

BiocManager::install("mixOmics", force = TRUE)

library(knitr)
library(xtable) # table
library(mixOmics)
library(sva) # ComBat
library(ggplot2) # PCA sample plot with density
library(gridExtra) # PCA sample plot with density
library(limma) # removeBatchEffect (LIMMA)
library(vegan) # RDA
library(AgiMicroRna) # RLE plot
library(cluster) # silhouette coefficient
library(variancePartition) # variance calculation
library(pvca) # PVCA
library(pheatmap) # heatmap
library(ruv) # RUVIII
library(lmerTest) # lmer
library(bapred) # FAbatch

setwd("/Users/hannahaichelman/Documents/BU/TVE/16S_ITS2/16S_All_Timepoints")

# data that is cleaned of euk contamination but nothing else
ps.cleanest = readRDS("phyloseq.cleanest.rds")
seqtab.cleanest <- data.frame(ps.cleanest@otu_table)
samdf.cleanest <- data.frame(ps.cleanest@sam_data)

# data that is cleaned and rarefied
ps.rare.1k = readRDS("phyloseq.rarefied.1k.rds")
seqtab.cleanest <- data.frame(ps.rare.1k@otu_table)
samdf.cleanest <- data.frame(ps.rare.1k@sam_data)

# data that is cleaned, trimmed, rarefied
ps.trim.rare = readRDS("phyloseq.trim.rare.rds")
seqtab.cleanest <- data.frame(ps.trim.rare@otu_table)
samdf.cleanest <- data.frame(ps.trim.rare@sam_data)


## Prefiltering
# ad data
seqtab.index.keep <- which(colSums(seqtab.cleanest)*100/(sum(colSums(seqtab.cleanest))) > 0.01)
seqtab.cleanest.keep <- seqtab.cleanest[, seqtab.index.keep]
dim(seqtab.cleanest.keep)
# [1]  237 1125

# Add offset to handle zeros
seqtab.cleanest.keep <- seqtab.cleanest.keep + 1

# Centered log-ratio transformation
seqtab.clr <- logratio.transfo(seqtab.cleanest.keep, logratio = 'CLR')
class(seqtab.clr) <- 'matrix'

## Batch effect detection
# pca
seqtab.pca.before <- pca(seqtab.clr, ncomp = 3)

data = as.data.frame(seqtab.pca.before$variates$X)
batch = as.factor(samdf.cleanest$time)
trt = as.factor(samdf.cleanest$treat)
expl.var = seqtab.pca.before$explained_variance

pMain <- ggplot(data = data, aes(x = data[ ,1], y = data[ ,2], colour = batch, shape = trt)) +
  geom_point()+
  xlab(paste0('PC1: ', round(as.numeric(expl.var[1])*100), '% expl.var')) +
  ylab(paste0('PC2: ', round(as.numeric(expl.var[2])*100), '% expl.var')) +
  scale_color_manual(values = color.mixo(1:10)) +
  theme_bw()+
  stat_ellipse()
#xlim(xlim[1], xlim[2]) +
#ylim(ylim[1], ylim[2]) +
#labs(colour = batch.legend.title, shape = trt.legend.title)
pMain
pTop <- ggplot(data,aes(x = data[ ,1], fill = batch, linetype = trt)) +
  geom_density(alpha = 0.5) +
  ylab('Density')
pRight <- ggplot(data, aes(x=data[ ,2], fill = batch, linetype=trt)) +
  geom_density(alpha = 0.5) +  coord_flip() + ylab('Density')

grid.arrange(pTop, pMain, pRight)



