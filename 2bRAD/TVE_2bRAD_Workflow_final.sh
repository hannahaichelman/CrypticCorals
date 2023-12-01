#This document details the data acquisition and analysis of the SSID 2bRAD data from the TVE experiment. 
#The goal of this analysis is to look for cryptic lineages of Siderastrea siderea from Bocas del Toro, Panama

#Generally following Misha Matz's pipeline, found here: https://github.com/z0on/2bRAD_denovo

#Download scripts from https://github.com/z0on/2bRAD_denovo
[haich@scc1 TVE_2bRAD]$ pwd
/projectnb/davies-hb/hannah/TVE_2bRAD

git clone https://github.com/z0on/2bRAD_denovo.git

# Make everything executable
chmod +x *.pl
chmod +x *.py
chmod +x *.R

# Install modules - we will need bowtie2, samtools, and picard. They are pre-installed as modules on TACC; you will have to install them if you don't have these modules on your cluster

module load perl
module load bowtie2
module load samtools
module load picard


# Working directory for this project is:
pwd
/projectnb/davies-hb/hannah/TVE_2bRAD


#------------------------------GET DATA FROM ILLUMINA BASE SPACE

# followed instructions here: https://docs.google.com/document/d/14S89ps4NB34CsFcis0XdWhObgBoLbM3t6-PJe5OHYcU/edit

# ftp in after you are already logged in to the SCC

ftp -i 130.64.74.72

#you will be prompted for username and password:

# change to working directory where our files are located
ftp> cd 210714-0477H_Sarah_Davies_6470_custom_demultiplex

# move the files
ftp> mget *.gz /projectnb/davies-hb/hannah/TVE_2bRAD/tufts_custom_pipeline_files/  

## The directory above is for the most recent (on 10/4/21) version of the customized files from Tufts
## This includes an analysis allowing for 2 mis-matches when assigning barcodes to sequences
## So, pulled the .gz files from the "standard_pipeline_no_mismatch" folder and the "customized_pipeline_2_mismatches" folder.

## These files are used for the first analysis listed here, the others below used the files from:
/projectnb/davies-hb/hannah/TVE_2bRAD/original_fastq_files
# zipp'd this folder so it looks like this: original_fastq_files.tar.gz

#------------------------------CHECK THAT DATA DOWNLOADED CORRECTLY

# Use the software md5 Master to compare the file id's with what we get from the sequencing facility, located in the md5sum.txt file, located here:

pwd  
/projectnb/davies-hb/hannah/TVE_2bRAD/tufts_custom_pipeline_files/customized_pipeline_2_mismatches

#and here:
pwd
/projectnb/davies-hb/hannah/TVE_2bRAD/tufts_custom_pipeline_files/standard_pipeline_no_mismatch

# everything matched except for the Undetermined.fastq file, which we won't use anyways.

#------------------------------UNZIP FILES - AND COMBINE CUSTOM FILES

# unzip files
# submitted job gunzip_files

#cat gunzip_files 
#!/bin/bash
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N gunzip_2 # job name, anything you want
#$ -l h_rt=24:00:00 #maximum run time
#$ -M hannahaichelman@gmail.com #your email
#$ -m be

gunzip *.gz


# renamed the standard files and the two mismatch files in each of their own directories and then cat'd them with this job:
pwd
/projectnb/davies-hb/hannah/TVE_2bRAD/tufts_custom_pipeline_files

cat cat_files
#!/bin/bash
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N cat_files # job name, anything you want
#$ -l h_rt=24:00:00 #maximum run time
#$ -M hannahaichelman@gmail.com #your email
#$ -m be

cat ./standard_pipeline_no_mismatch/1-MullenDavies_S1_std.fastq ./customized_pipeline_2_mismatches/1-MullenDavies_S1_msmtch.fastq > 1-MullenDavies_S1.fastq
cat ./standard_pipeline_no_mismatch/2-MullenDavies_S2_std.fastq ./customized_pipeline_2_mismatches/2-MullenDavies_S2_msmtch.fastq > 2-MullenDavies_S2.fastq
cat ./standard_pipeline_no_mismatch/3-MullenDavies_S3_std.fastq ./customized_pipeline_2_mismatches/3-MullenDavies_S3_msmtch.fastq > 3-MullenDavies_S3.fastq
cat ./standard_pipeline_no_mismatch/4-MullenDavies_S4_std.fastq ./customized_pipeline_2_mismatches/4-MullenDavies_S4_msmtch.fastq > 4-MullenDavies_S4.fastq
cat ./standard_pipeline_no_mismatch/5-MullenDavies_S5_std.fastq ./customized_pipeline_2_mismatches/5-MullenDavies_S5_msmtch.fastq > 5-MullenDavies_S5.fastq
cat ./standard_pipeline_no_mismatch/6-MullenDavies_S6_std.fastq ./customized_pipeline_2_mismatches/6-MullenDavies_S6_msmtch.fastq > 6-MullenDavies_S6.fastq
cat ./standard_pipeline_no_mismatch/7-MullenDavies_S7_std.fastq ./customized_pipeline_2_mismatches/7-MullenDavies_S7_msmtch.fastq > 7-MullenDavies_S7.fastq
cat ./standard_pipeline_no_mismatch/8-MullenDavies_S8_std.fastq ./customized_pipeline_2_mismatches/8-MullenDavies_S8_msmtch.fastq > 8-MullenDavies_S8.fastq

# Because a report was not generated for the custom_pipeline_2_mismatches files, TUCF provided the read count information here:
# Order 6470, 2 mismatches
# Sample Read Count
# 1-MullenDavies_S1_001.fastq.gz 13632065
# 2-MullenDavies_S2_001.fastq.gz 11908010
# 3-MullenDavies_S3_001.fastq.gz 13388678
# 4-MullenDavies_S4_001.fastq.gz 9662566
# 5-MullenDavies_S5_001.fastq.gz 10082082
# 6-MullenDavies_S6_001.fastq.gz 8504871
# 7-MullenDavies_S7_001.fastq.gz 12029357
# 8-MullenDavies_S8_001.fastq.gz 9263545
# Undetermined_S0_001.fastq.gz 129861144  (59.48%)


#------------------------------BACKING UP DATA
#tarball'd the back-up files in case something goes wrong 
#these are the original files without the demultiplexing with the Tufts custom script:

pwd
/projectnb/davies-hb/hannah/TVE_2bRAD

tar -zcvf original_fastq_files.tar.gz original_fastq_files/

#and for the custom demultiplexed files from Tufts' custom script (same directory):

tar -zcvf tufts_custom_pipeline_fastqs.tar.gz tufts_custom_pipeline_fastqs/


#------------------------------DEMULTIPLEXING, TRIMMING AND MAPPING 

#Download scripts from https://github.com/z0on/2bRAD_denovo
pwd
/projectnb/davies-hb/hannah/TVE_2bRAD

git clone https://github.com/z0on/2bRAD_denovo.git

# Make everything executable
chmod +x *.pl
chmod +x *.py
chmod +x *.R

# Install modules - we will need bowtie2, samtools, and picard. 

module load perl
module load bowtie2
module load samtools
module load picard


# Using the denovo method here to create a reference since we do not have a Siderastrea siderea genome

#==================
# Step 1: Splitting by in-read barcode, deduplicating and quality-filtering the reads

pwd
/projectnb/davies-hb/hannah/TVE_2bRAD/tufts_custom_pipeline_files

# creating a file of commands to run (assuming reads are in fastq files, one file per sample.)
 
../2bRAD_denovo/2bRAD_trim_launch_dedup.pl fastq > trims

# Modify this trims file to be able to submit it as a job on the SCC, designate where the perl script is, and add '&'s to the end of the line
# It looks like this:

cat trims
#!/bin/bash
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N trims # job name, anything you want
#$ -l h_rt=24:00:00 #maximum run time
#$ -M hannahaichelman@gmail.com #your email
#$ -m be

module load perl

../2bRAD_denovo/trim2bRAD_2barcodes_dedup.pl input=7-MullenDavies_S7.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../2bRAD_denovo/trim2bRAD_2barcodes_dedup.pl input=1-MullenDavies_S1.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../2bRAD_denovo/trim2bRAD_2barcodes_dedup.pl input=8-MullenDavies_S8.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../2bRAD_denovo/trim2bRAD_2barcodes_dedup.pl input=3-MullenDavies_S3.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../2bRAD_denovo/trim2bRAD_2barcodes_dedup.pl input=5-MullenDavies_S5.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../2bRAD_denovo/trim2bRAD_2barcodes_dedup.pl input=4-MullenDavies_S4.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../2bRAD_denovo/trim2bRAD_2barcodes_dedup.pl input=6-MullenDavies_S6.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../2bRAD_denovo/trim2bRAD_2barcodes_dedup.pl input=2-MullenDavies_S2.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
wait


# do we have expected number of *.tr0 files created?
ls -l *.tr0 | wc -l

# We have 93 samples vs the 86 that we had before Tufts' custom pipeline!


# for denovo analysis: removing reads with qualities at ends less than Q15

>trimse
for file in *.tr0; do
echo "cutadapt -q 15,15 -m 36 -o ${file/.tr0/}.trim $file > ${file}_trimlog.txt" >> trimse;
done

module load cutadapt

# add job header to the top of the file and qsub
# looks like this:

head -20 trimse 
#!/bin/bash
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N trimse # job name, anything you want
#$ -l h_rt=24:00:00 #maximum run time
#$ -M hannahaichelman@gmail.com #your email
#$ -m be

cutadapt -q 15,15 -m 36 -o 1-MullenDavies_S1_L002_R1_001_ACCA.trim 1-MullenDavies_S1_L002_R1_001_ACCA.tr0 > 1-MullenDavies_S1_L002_R1_001_ACCA.tr0_trimlog.txt
cutadapt -q 15,15 -m 36 -o 1-MullenDavies_S1_L002_R1_001_AGAC.trim 1-MullenDavies_S1_L002_R1_001_AGAC.tr0 > 1-MullenDavies_S1_L002_R1_001_AGAC.tr0_trimlog.txt
cutadapt -q 15,15 -m 36 -o 1-MullenDavies_S1_L002_R1_001_AGTG.trim 1-MullenDavies_S1_L002_R1_001_AGTG.tr0 > 1-MullenDavies_S1_L002_R1_001_AGTG.tr0_trimlog.txt
cutadapt -q 15,15 -m 36 -o 1-MullenDavies_S1_L002_R1_001_CTAC.trim 1-MullenDavies_S1_L002_R1_001_CTAC.tr0 > 1-MullenDavies_S1_L002_R1_001_CTAC.tr0_trimlog.txt
cutadapt -q 15,15 -m 36 -o 1-MullenDavies_S1_L002_R1_001_GACT.trim 1-MullenDavies_S1_L002_R1_001_GACT.tr0 > 1-MullenDavies_S1_L002_R1_001_GACT.tr0_trimlog.txt
cutadapt -q 15,15 -m 36 -o 1-MullenDavies_S1_L002_R1_001_GCTT.trim 1-MullenDavies_S1_L002_R1_001_GCTT.tr0 > 1-MullenDavies_S1_L002_R1_001_GCTT.tr0_trimlog.txt
cutadapt -q 15,15 -m 36 -o 1-MullenDavies_S1_L002_R1_001_TCAC.trim 1-MullenDavies_S1_L002_R1_001_TCAC.tr0 > 1-MullenDavies_S1_L002_R1_001_TCAC.tr0_trimlog.txt
cutadapt -q 15,15 -m 36 -o 1-MullenDavies_S1_L002_R1_001_TCAG.trim 1-MullenDavies_S1_L002_R1_001_TCAG.tr0 > 1-MullenDavies_S1_L002_R1_001_TCAG.tr0_trimlog.txt
cutadapt -q 15,15 -m 36 -o 1-MullenDavies_S1_L002_R1_001_TGGT.trim 1-MullenDavies_S1_L002_R1_001_TGGT.tr0 > 1-MullenDavies_S1_L002_R1_001_TGGT.tr0_trimlog.txt
cutadapt -q 15,15 -m 36 -o 1-MullenDavies_S1_L002_R1_001_TGTC.trim 1-MullenDavies_S1_L002_R1_001_TGTC.tr0 > 1-MullenDavies_S1_L002_R1_001_TGTC.tr0_trimlog.txt
cutadapt -q 15,15 -m 36 -o 2-MullenDavies_S2_L002_R1_001_ACCA.trim 2-MullenDavies_S2_L002_R1_001_ACCA.tr0 > 2-MullenDavies_S2_L002_R1_001_ACCA.tr0_trimlog.txt
cutadapt -q 15,15 -m 36 -o 2-MullenDavies_S2_L002_R1_001_AGAC.trim 2-MullenDavies_S2_L002_R1_001_AGAC.tr0 > 2-MullenDavies_S2_L002_R1_001_AGAC.tr0_trimlog.txt

# Cutadapt info:
# -q is a quality filter, used to trim low-quality ends from reads. 
# The comma separated argument here trims the 5' end with a cutoff of 15 and the 3' end with a cutoff of 15.
# -m is a minimum length filter, it discards processed reads that are shorter than LENGTH provided (here, 36).


# do we have expected number of *.trim files created?
ls -l *.trim | wc -l

# still 93, makes sense!


#------------------------------SEPARATING SEAN AND I'S SAMPLES

# We sequenced two projects together, so I separated files out based on SampleMap.xlsx, and FinalPlateMap.xlsx. The other project's samples are now in their own folder
# The original files without doing the custom script are here:
pwd
/projectnb/davies-hb/hannah/TVE_2bRAD/original_files_analyses/sean_files

# More info on this google drive link: https://docs.google.com/spreadsheets/d/1IETRukTXwhb5Yd9co5DSqDst1m6YbrSppdpGZNItGS4/edit?usp=sharing

# Sean now has 32 .trim files
# I now have 61 samples, which is all of my samples!
 
#------------------------------REMOVING SYM READS BEFORE STARTING DE NOVO ANALYSIS

# Following JP Rippe's github here:
# https://github.com/jprippe/AdultJuv_Depth_FL/blob/master/Analysis_walkthroughs/ssid_processing.txt 

#------------ setting up Symbiodinium genomes

export SYM_FASTA=/projectnb/davies-hb/hannah/TVE_2bRAD/Sym_Genomes/symABCD_genome.fasta
export SYM_DICT=/projectnb/davies-hb/hannah/TVE_2bRAD/Sym_Genomes/symABCD_genome.dict 

# index genome for bowtie2 mapper
# submitted this as a job (name: btb) because it got killed the first time - ran too long

pwd
/projectnb/davies-hb/hannah/TVE_2bRAD/Sym_Genomes

cat btb
#!/bin/bash
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N btb # job name, anything you want
#$ -l h_rt=24:00:00 #maximum run time
#$ -M hannahaichelman@gmail.com #your email
#$ -m be

module load bowtie2

bowtie2-build symABCD_genome.fasta symABCD_genome.fasta

# this took ~2 hours to run as a job

module load samtools
samtools faidx symABCD_genome.fasta 
picard CreateSequenceDictionary R=symABCD_genome.fasta O=symABCD_genome.dict

#------------ Mapping to Symbiodiniaceae genomes, discarding reads that stick


SYM_REF=/projectnb/davies-hb/hannah/TVE_2bRAD/Sym_Genomes/symABCD_genome.fasta

# we are starting with the files that were demultiplexed and trimmed for denovo above (tufts custom files)

pwd
/projectnb/davies-hb/hannah/TVE_2bRAD/tufts_custom_pipeline_files

>mapsym
for F in `ls *.trim`; do
echo "bowtie2 -x $SYM_REF -U $F -S $F.sam">>mapsym
done

qsub mapsym


# saving fastq reads that do NOT map to chr11-14 (Symbiodinium genomes)

pwd
/projectnb/davies-hb/hannah/TVE_2bRAD/tufts_custom_pipeline_files

>sam2fq
for S in `ls *.sam`; do
F=`echo $S | perl -pe 's/\..+//'`;
echo "cat $S | awk '\$3==\"*\"' | cut -f1,10,11 | sed 's/^/@/' | sed 's/\t/\n/' | sed 's/\t/\n+\n/' > $F.nosymbio.fastq">>sam2fq;
done

qsub sam2fq

# the command above converts SAM files back to FASTQ
# Subset all lines where the third element is "*" (i.e., the sequence did not align to the Symbiodiniaceae genome), extract elements 1, 10 and 11, add @ to beginning of each line, change first tab delimiter (\t) to new line (\n) delimiter, change second tab delimiter to new line delimiter with "+" on a new line in between


#------------------------------DE NOVO RAD BUSINESS

# mv all *nosymbio* files to a new direcotry

mkdir denovo_nosym
mv *nosymbio* denovo_nosym/

pwd
/projectnb/davies-hb/hannah/TVE_2bRAD/tufts_custom_pipeline_files/denovo_nosym


# 'uniquing' ('stacking') individual fastq reads:

ls *.nosymbio.fastq | perl -pe 's/^(.+)$/uniquerOne.pl $1 >$1\.uni/' >unii

head -20 unii
#!/bin/bash
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N unii # job name, anything you want
#$ -l h_rt=24:00:00 #maximum run time
#$ -M hannahaichelman@gmail.com #your email
#$ -m be

../../2bRAD_denovo/uniquerOne.pl 1-MullenDavies_S1_ACCA.nosymbio.fastq >1-MullenDavies_S1_ACCA.nosymbio.fastq.uni
../../2bRAD_denovo/uniquerOne.pl 1-MullenDavies_S1_AGAC.nosymbio.fastq >1-MullenDavies_S1_AGAC.nosymbio.fastq.uni
../../2bRAD_denovo/uniquerOne.pl 1-MullenDavies_S1_AGTG.nosymbio.fastq >1-MullenDavies_S1_AGTG.nosymbio.fastq.uni
../../2bRAD_denovo/uniquerOne.pl 1-MullenDavies_S1_CATC.nosymbio.fastq >1-MullenDavies_S1_CATC.nosymbio.fastq.uni


# uniquerOne.pl :
# Makes uniqued 2bRAD read file for a single fastq file. This is analogous to making 'stacks' in STACKS. The script records unique tag sequences, the total number of their appearances, and how many of those were in reverse-complement orientation. (STACKS would consider reverse-complements as separate loci)

# merging uniqued files (set minInd to >10, or >10% of total number of samples, whatever is greater)
# since I have 61 samples I am keeping the minInd at 10

../../2bRAD_denovo/mergeUniq.pl uni minInd=10>all.uniq

# 313722 tags before tufts custom code
# 178617 tags after


# discarding tags that have more than 7 observations without reverse-complement
awk '!($3>7 && $4==0) && $2!="seq"' all.uniq >all.tab

wc -l all.tab  

# 147539 all.tab # before tufts custom code
# 177670 all.tab # after tufts custom code


# creating fasta file out of merged and filtered tags:
awk '{print ">"$1"\n"$2}' all.tab > all.fasta

# clustering reads into loci using cd-hit
# clustering allowing for up to 3 mismatches (-c 0.91); the most abundant sequence becomes reference
# -aL : alignment coverage for the longer sequence (setting at 1 means the alignment must cover 100% of the sequence)
# -aS : alignment coverage for the shorter sequence (setting to 1 means the alignment must cover 100% of the sequence
# -c : sequence identity threshold, calculated as # of identical bases in alignment / full length of the shorter sequence
# -g : the program will cluster it into the most similar cluster that meets the threshold, rather than the first
# -M -T : memory limit and # of threads to use (0 sets no limits)

module load blast
module load cdhit/4.6.8

cd-hit-est -i all.fasta -o cdh_alltags.fas -aL 1 -aS 1 -g 1 -c 0.91 -M 0 -T 0  

#------------
# making fake reference genome (of 30 chromosomes) out of cd-hit cluster representatives
# need bowtie2, samtools and picard_tools for indexing

module load bowtie2
module load samtools
module load picard
module load perl
../../2bRAD_denovo/concatFasta.pl fasta=cdh_alltags.fas num=30

concatenating 99753 records into 30 pseudo-chromosomes 
3326 records per chromosome

# formatting fake genome

export GENOME_FASTA=cdh_alltags_cc.fasta

bowtie2-build $GENOME_FASTA $GENOME_FASTA

samtools faidx $GENOME_FASTA


#==============
# Mapping reads to denovo reference and formatting bam files 

# for denovo: map with bowtie2 with default parameters
 
>maps
for F in `ls *.nosymbio.fastq`; do
echo "bowtie2 --no-unal -x $GENOME_FASTA -U $F -S $F.sam">>maps
done

qsub maps

ls *nosymbio.fastq.sam > sams
cat sams | wc -l  # number should match number of.fastq files
# 61 - matches

# next stage is compressing, sorting and indexing the SAM files, so they become BAM files:
cat sams | perl -pe 's/(\S+)\.sam/samtools view -bS $1\.sam >$1\.unsorted\.bam && samtools sort $1\.unsorted\.bam -o $1\.bam && samtools index $1\.bam/' >s2b

qsub s2b

rm -f *unsorted*
ls *bam | wc -l  # should be the same number as number of .fastq files
# 61! 

# BAM files are the input into various genotype calling / popgen programs, this is the main interim result of the analysis. Archive them.

#------------ quality assessment

>alignmentRates
>align
for F in `ls *nosymbio.fastq`; do 
echo "grep -E '^[ATGCN]+$' $F | wc -l | grep -f - maps.e* -A 4 | tail -1 | perl -pe 's/maps\.e\d+-|% overall alignment rate//g' | xargs echo $F.sam >> alignmentRates" >>align; 
done

qsub align

# this is still not working for me...alignmentRates is just a list of the .sam file names after the job finishes running


#------------------------------ANGSD (fuzzy genotyping)
# this first run-through of ANGSD is just to take a look at base qualities and coverage depth, 
# we will run angsd again with filters informed by this first run-through

pwd
/projectnb/davies-hb/hannah/TVE_2bRAD/tufts_custom_pipeline_files/denovo_nosym

# "FUZZY genotyping" with ANGSD - without calling actual genotypes but working with genotype likelihoods at each SNP. Optimal for low-coverage data (<10x).

# listing all bam filenames 
ls *bam >bams

module load angsd

#----------- assessing base qualities and coverage depth

# angsd settings:
# -minMapQ 20 : only highly unique mappings (prob of erroneous mapping =< 1%)
# -baq 1 : realign around indels (not terribly relevant for 2bRAD reads mapped with --local option) 
# -maxDepth : highest total depth (sum over all samples) to assess; set to 10x number of samples
# -minInd : the minimal number of individuals the site must be genotyped in. Reset to 50% of total N at this stage.

FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -maxDepth 610 -minInd 30"

# T O   D O : 
TODO="-doQsDist 1 -doDepth 1 -doCounts 1 -dumpCounts 2"

# in the following line, -r argument is one chromosome or contig to work with (no need to do this for whole genome as long as the chosen chromosome or contig is long enough, ~1 Mb. When mapping to a real genome, consider chr1:1-1000000 )
# (look up lengths of your contigs in the header of *.sam files)
# i'm not using here and just running on entire reference

angsd -b bams -GL 1 $FILTERS $TODO -P 1 -out dd


# summarizing results (using modified script by Matteo Fumagalli)
module load R
Rscript ../../2bRAD_denovo/plotQC.R prefix=dd

# proportion of sites covered at >5x:
cat quality.txt
# 1-MullenDavies_S1_GTGA.nosymbio.fastq.bam 0.084345091093719
# 4-MullenDavies_S4_AGAC.nosymbio.fastq.bam 0.0968114006980149
# 3-MullenDavies_S3_GACT.nosymbio.fastq.bam 0.107732290119309
# 5-MullenDavies_S5_AGTG.nosymbio.fastq.bam 0.134728985500961
# 1-MullenDavies_S1_CATC.nosymbio.fastq.bam 0.135823141736968
# 4-MullenDavies_S4_CATC.nosymbio.fastq.bam 0.143300634937243
# 1-MullenDavies_S1_ACCA.nosymbio.fastq.bam 0.147309278575941
# 1-MullenDavies_S1_AGAC.nosymbio.fastq.bam 0.149964447724493
# 4-MullenDavies_S4_TGTC.nosymbio.fastq.bam 0.157553578264309
# 5-MullenDavies_S5_ACCA.nosymbio.fastq.bam 0.158757360917282
# 2-MullenDavies_S2_CATC.nosymbio.fastq.bam 0.161220716419526
# 3-MullenDavies_S3_TGGT.nosymbio.fastq.bam 0.164771381355608
# 2-MullenDavies_S2_TCAG.nosymbio.fastq.bam 0.167408587585371
# 4-MullenDavies_S4_GCTT.nosymbio.fastq.bam 0.1685672958044
# 3-MullenDavies_S3_CTAC.nosymbio.fastq.bam 0.171335854702653
# 5-MullenDavies_S5_CATC.nosymbio.fastq.bam 0.17507698641693
# 4-MullenDavies_S4_GACT.nosymbio.fastq.bam 0.178608608009673
# 4-MullenDavies_S4_CTAC.nosymbio.fastq.bam 0.186536106441305
# 1-MullenDavies_S1_TGGT.nosymbio.fastq.bam 0.187915980778309
# 2-MullenDavies_S2_CTAC.nosymbio.fastq.bam 0.191263804784384
# 2-MullenDavies_S2_ACCA.nosymbio.fastq.bam 0.201229072650576
# 1-MullenDavies_S1_GCTT.nosymbio.fastq.bam 0.203782032107712
# 1-MullenDavies_S1_TCAG.nosymbio.fastq.bam 0.211321045549797
# 4-MullenDavies_S4_ACCA.nosymbio.fastq.bam 0.213664646269978
# 3-MullenDavies_S3_CATC.nosymbio.fastq.bam 0.213936015238985
# 5-MullenDavies_S5_TGGT.nosymbio.fastq.bam 0.21572050055017
# 4-MullenDavies_S4_AGTG.nosymbio.fastq.bam 0.241255367938926
# 4-MullenDavies_S4_TGGT.nosymbio.fastq.bam 0.249611357664273
# 4-MullenDavies_S4_TCAC.nosymbio.fastq.bam 0.253727085795271
# 4-MullenDavies_S4_GTGA.nosymbio.fastq.bam 0.256016582413145
# 2-MullenDavies_S2_GTGA.nosymbio.fastq.bam 0.257896683348385
# 2-MullenDavies_S2_GACT.nosymbio.fastq.bam 0.261223191930802
# 1-MullenDavies_S1_GACT.nosymbio.fastq.bam 0.295926730439897
# 3-MullenDavies_S3_AGTG.nosymbio.fastq.bam 0.32260284915647
# 2-MullenDavies_S2_AGAC.nosymbio.fastq.bam 0.323732946972281
# 8-MullenDavies_S8_GTGA.nosymbio.fastq.bam 0.325772648513097
# 5-MullenDavies_S5_GCTT.nosymbio.fastq.bam 0.330937966873282
# 2-MullenDavies_S2_TCAC.nosymbio.fastq.bam 0.337015007128901
# 2-MullenDavies_S2_AGTG.nosymbio.fastq.bam 0.34767540167728
# 5-MullenDavies_S5_GTGA.nosymbio.fastq.bam 0.351132286858971
# 5-MullenDavies_S5_GACT.nosymbio.fastq.bam 0.360611860249588
# 5-MullenDavies_S5_TCAC.nosymbio.fastq.bam 0.37669904363539
# 3-MullenDavies_S3_TCAC.nosymbio.fastq.bam 0.385333375345226
# 2-MullenDavies_S2_GCTT.nosymbio.fastq.bam 0.388500871930495
# 3-MullenDavies_S3_AGAC.nosymbio.fastq.bam 0.397788901368651
# 1-MullenDavies_S1_AGTG.nosymbio.fastq.bam 0.400808378560385
# 5-MullenDavies_S5_AGAC.nosymbio.fastq.bam 0.410006551256408
# 3-MullenDavies_S3_ACCA.nosymbio.fastq.bam 0.413834531060923
# 4-MullenDavies_S4_TCAG.nosymbio.fastq.bam 0.417344477729921
# 6-MullenDavies_S6_TGGT.nosymbio.fastq.bam 0.484359772831813
# 3-MullenDavies_S3_TGTC.nosymbio.fastq.bam 0.506624705504068
# 8-MullenDavies_S8_AGAC.nosymbio.fastq.bam 0.513818396977772
# 3-MullenDavies_S3_TCAG.nosymbio.fastq.bam 0.536820451974425
# 3-MullenDavies_S3_GCTT.nosymbio.fastq.bam 0.580769465856865
# 8-MullenDavies_S8_TGGT.nosymbio.fastq.bam 0.582074588588805
# 2-MullenDavies_S2_TGGT.nosymbio.fastq.bam 0.601116685293984
# 1-MullenDavies_S1_TCAC.nosymbio.fastq.bam 0.626660892952064
# 1-MullenDavies_S1_TGTC.nosymbio.fastq.bam 0.628979622484748
# 2-MullenDavies_S2_TGTC.nosymbio.fastq.bam 0.629282005629554
# 3-MullenDavies_S3_GTGA.nosymbio.fastq.bam 0.731616565034842
# 1-MullenDavies_S1_CTAC.nosymbio.fastq.bam 0.837796101341018


# scp dd.pdf to laptop to look at distribution of base quality scores, fraction of sites in each sample passing coverage thresholds, and fraction of sites passing genotyping rates cutoffs. Use these to guide choices of -minQ,  -minIndDepth and -minInd filters in subsequent ANGSD runs

# try again without the -baq, since this was a new filter for me - JP Rippe used it but Misha doesn't have it in his run through

FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -maxDepth 590 -minInd 30"

# T O   D O : 
TODO="-doQsDist 1 -doDepth 1 -doCounts 1 -dumpCounts 2"

# in the following line, -r argument is one chromosome or contig to work with (no need to do this for whole genome as long as the chosen chromosome or contig is long enough, ~1 Mb. When mapping to a real genome, consider chr1:1-1000000 )
# (look up lengths of your contigs in the header of *.sam files)
angsd -b bams -GL 1 $FILTERS $TODO -P 1 -out dd

# looks like we get the same result with or without -baq filter, just slightly fewer SNPs with the -baq filter


#----------- clones detection (**minInd ~80% of samples)

# the -doVcf flag produces a vcf file at this stage.

pwd
/projectnb/davies-hb/hannah/TVE_2bRAD/tufts_custom_pipeline_files/denovo_nosym

FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 49 -snp_pval 1e-5 -minMaf 0.05"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doVcf 1 -doPost 1 -doGlf 2"

angsd -b bams -GL 1 $FILTERS $TODO -P 1 -out myresult

NSITES=`zcat myresult.mafs.gz | wc -l`
echo $NSITES
# 11635 before tufts pipeline
# 14319 after tufts custom pipeline


# Add the -minIndDepth filter to see what happens, this filter was really important in Nicola's project but not included in Misha's pipeline. 
# Starting with -minIndDepth 2 and see what happens when we increase it

FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 47 -minIndDepth 2 -snp_pval 1e-5 -minMaf 0.05"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doVcf 1 -doPost 1 -doGlf 2"

angsd -b bams -GL 1 $FILTERS $TODO -P 1 -out myresult2

NSITES=`zcat myresult2.mafs.gz | wc -l`
echo $NSITES
# 3867 before tufts pipeline
# 7859 after tufts custom pipeline

# Now the -minIndDepth 3

FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 47 -minIndDepth 3 -snp_pval 1e-5 -minMaf 0.05"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doVcf 1 -doPost 1 -doGlf 2"

angsd -b bams -GL 1 $FILTERS $TODO -P 1 -out myresult3

NSITES=`zcat myresult3.mafs.gz | wc -l`
echo $NSITES
# 1172 before tufts pipeline
# 3653 after tufts custom pipeline

# scp *.ibsMat and bams to laptop, use ibs.R to analyze

# On my local machine, the denovo without sym contamination files are located here:

pwd
/Users/hannahaichelman/Documents/BU/TVE/2bRAD/Analysis/old_denovo_nosyms

# the denovo without sym contamination files plus the tufts custom pipeline files are here:

pwd
/Users/hannahaichelman/Documents/BU/TVE/2bRAD/Analysis/tuftscustompipeline_denovo_nosyms


#------------------------------RE-RUN WITHOUT CLONES - Lineage assignments

pwd
TVE_2bRAD/tufts_custom_pipeline_files/denovo_nosym

# first with all duplicate preps removed, but with the two clones detected still included
# note the -minInd changes, to ~80% of our now 51 samples

FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 41 -minIndDepth 2 -snp_pval 1e-5 -minMaf 0.05"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doVcf 1 -doPost 1 -doGlf 2"


angsd -b bams_noclones_allsamps -GL 1 $FILTERS $TODO -P 1 -out myresult2.noclone.allsamps

NSITES=`zcat myresult2.noclone.allsamps.mafs.gz | wc -l`
echo $NSITES
# 7176

# now with all clones removed, left with 50 samples

FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 40 -minIndDepth 2 -snp_pval 1e-5 -minMaf 0.05"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doVcf 1 -doPost 1 -doGlf 2"


angsd -b bams_noclones -GL 1 $FILTERS $TODO -P 1 -out myresult2.noclone

NSITES=`zcat myresult2.noclone.mafs.gz | wc -l`
echo $NSITES
# 8105

#------------------------------NgsAdmix

# first, install:
pwd
/projectnb/davies-hb/hannah/TVE_2bRAD/bin

wget popgen.dk/software/download/NGSadmix/ngsadmix32.cpp 

# all set, downloaded in bin.

# NgsAdmix for K from 2 to 5 : do not run if the dataset contains clones or genotyping replicates!

for K in `seq 2 5` ;
do
../../bin/NGSadmix -likes myresult2.noclone.beagle.gz -K $K -P 10 -o mydata.noclone_k${K};
done

# scp *qopt files to local machine, use admixturePlotting_v5.R to plot


## alternatively, to use real ADMIXTURE on called SNPs (requires plink and ADMIXTURE):

pwd
/projectnb/davies-hb/hannah/TVE_2bRAD/tufts_custom_pipeline_files/denovo_nosym

module load admixture/1.3.0
module load plink/1.90b6.4

# unzip vcf file:
gunzip myresult2.noclone.vcf.gz

# remove the "(angsd version)" from the head of the file, this is what was causing errors before

plink --vcf myresult2.noclone.vcf --make-bed --allow-extra-chr 0 --out myresult2.noclone
for K in `seq 1 5`; \
do admixture --cv myresult2.noclone.bed $K | tee myresult2.noclone_${K}.out; done

## This isn't working for me, I get the following error:
Error opening .bim file!


# if it was could use this to check K of least CV error:
#grep -h CV myresult.noLD_*.out

#CV error (K=1): 0.53296
#CV error (K=2): 0.45103
#CV error (K=3): 0.43497
#CV error (K=4): 0.43349
#CV error (K=5): 0.44908


#------------------------------Analysis of Genetic Divergence
pwd
/projectnb/davies-hb/hannah/TVE_2bRAD/tufts_custom_pipeline_files/denovo_nosym

# with all clones removed, we have 50 samples

module load angsd/0.923

# filtering sites to work on - use only filters that do not distort allele frequency
# set minInd to 75-90% of the total number fo individuals in the project
# if you are doing any other RAD than 2bRAD or GBS, remove '-sb_pval 1e-5' from FILTERS
# -maxHetFreq 0.5 is the lumped paralog filter - important for doing SFS analyses
FILTERS="-uniqueOnly 1 -remove_bads 1  -skipTriallelic 1 -minMapQ 25 -minQ 30 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -minInd 40"
TODO="-doMajorMinor 1 -doMaf 1 -dosnpstat 1 -doPost 2 -doGeno 11 -doGlf 2 -doVcf 1"


angsd -b bams_noclones -GL 1 $FILTERS $TODO -out sfilt


# -maxHetFreq 0.5 is the lumped paralog filter - important for doing SFS analyses. seeing if it makes a difference for Fst values
FILTERS="-uniqueOnly 1 -remove_bads 1  -skipTriallelic 1 -minMapQ 25 -minQ 30 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -minInd 40 -maxHetFreq 0.5"
TODO="-doMajorMinor 1 -doMaf 1 -dosnpstat 1 -doPost 2 -doGeno 11 -doGlf 2 -doBcf 1"


angsd -b bams_noclones -GL 1 $FILTERS $TODO -out sfilt_maxhet

NSITES=`zcat sfilt_maxhet.mafs.gz | wc -l`
echo $NSITES
# 770398

# make vcf only for running bayescan:
FILTERS="-uniqueOnly 1 -remove_bads 1  -skipTriallelic 1 -minMapQ 25 -minQ 30 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -minInd 40 -snp_pval 1e-5 -minMaf 0.05"
TODO="-doMajorMinor 1 -doMaf 1 -dosnpstat 1 -doPost 2 -doGeno 11 -doGlf 2 -doVcf 1"

angsd -b bams_noclones -GL 1 $FILTERS $TODO -out sfilt_forbayescan


# collecting and indexing filter-passing sites
zcat sfilt.mafs.gz | cut -f 1,2 | tail -n +2 >allSites
angsd sites index allSites

# max het filter
zcat sfilt_maxhet.mafs.gz | cut -f 1,2 | tail -n +2 >allSites_maxhet
angsd sites index allSites_maxhet

export GENOME_REF=/projectnb/davies-hb/hannah/TVE_2bRAD/tufts_custom_pipeline_files/denovo_nosym/cdh_alltags_cc.fasta
TODO="-doSaf 1 -doMajorMinor 1 -doMaf 1 -doPost 1 -anc $GENOME_REF -ref $GENOME_REF"


angsd -sites allSites -b bams_L1 -GL 1 -P 1 $TODO -out L1
angsd -sites allSites -b bams_L2 -GL 1 -P 1 $TODO -out L2
angsd -sites allSites -b bams_L3 -GL 1 -P 1 $TODO -out L3

# max het filter
angsd -sites allSites_maxhet -b bams_L1 -GL 1 -P 1 $TODO -out L1_maxhet
angsd -sites allSites_maxhet -b bams_L2 -GL 1 -P 1 $TODO -out L2_maxhet
angsd -sites allSites_maxhet -b bams_L3 -GL 1 -P 1 $TODO -out L3_maxhet

# generating per-population SFS
realSFS L1.saf.idx >L1.sfs
realSFS L2.saf.idx >L2.sfs
realSFS L3.saf.idx >L3.sfs


# generating per-population SFS - max het filter
realSFS L1_maxhet.saf.idx >L1_maxhet.sfs
realSFS L2_maxhet.saf.idx >L2_maxhet.sfs
realSFS L3_maxhet.saf.idx >L3_maxhet.sfs

# writing down 2d-SFS priors - L1 vs L2
realSFS L1_maxhet.saf.idx L2_maxhet.saf.idx -P 24 > p12_maxhet.sfs ; realSFS fst index L1_maxhet.saf.idx L2_maxhet.saf.idx -sfs p12_maxhet.sfs -fstout p12_maxhet

# global Fst between populations
realSFS fst stats p12_maxhet.fst.idx

#output:
	-> Assuming idxname:p12.fst.idx
	-> Assuming .fst.gz file: p12.fst.gz
	-> FST.Unweight[nObs:771675]:0.023862 Fst.Weight:0.159775
0.023862 0.159775

	-> Assuming idxname:p12_maxhet.fst.idx
	-> Assuming .fst.gz file: p12_maxhet.fst.gz
	-> FST.Unweight[nObs:770220]:0.023784 Fst.Weight:0.172208
0.023784	0.172208

# The difference between unweighted and weighted values is averaging versus ratio of sums method. 
# The weighted value is derived from the ratio of the separate per-locus sums of numerator and denominator values, 
# while the unweighted value is the average of per-locus values. [If that is unclear: weighted is sum(a)/sum(a+b), while unweighted is average(a/(a+b))].


# writing down 2d-SFS priors - L1 vs L3
realSFS L1_maxhet.saf.idx L3_maxhet.saf.idx -P 24 > p13_maxhet.sfs ; realSFS fst index L1_maxhet.saf.idx L3_maxhet.saf.idx -sfs p13_maxhet.sfs -fstout p13_maxhet

# global Fst between populations
realSFS fst stats p13_maxhet.fst.idx

#output:
	-> Assuming idxname:p13.fst.idx
	-> Assuming .fst.gz file: p13.fst.gz
	-> FST.Unweight[nObs:761904]:0.326915 Fst.Weight:0.162287
0.326915 0.162287

	-> Assuming idxname:p13_maxhet.fst.idx
	-> Assuming .fst.gz file: p13_maxhet.fst.gz
	-> FST.Unweight[nObs:760462]:0.327435 Fst.Weight:0.178868
0.327435	0.178868

# writing down 2d-SFS priors - L2 vs L3
realSFS L2_maxhet.saf.idx L3_maxhet.saf.idx -P 24 > p23_maxhet.sfs ; realSFS fst index L2_maxhet.saf.idx L3_maxhet.saf.idx -sfs p23_maxhet.sfs -fstout p23_maxhet

# global Fst between populations
realSFS fst stats p23_maxhet.fst.idx

#output:
	-> Assuming idxname:p23.fst.idx
	-> Assuming .fst.gz file: p23.fst.gz
	-> FST.Unweight[nObs:761878]:0.178944 Fst.Weight:0.108837
0.178944 0.108837

	-> Assuming idxname:p23_maxhet.fst.idx
	-> Assuming .fst.gz file: p23_maxhet.fst.gz
	-> FST.Unweight[nObs:760437]:0.178861 Fst.Weight:0.117511
0.178861	0.117511

# using the max het filter here (removing lumped paralogs) since it is most appropriate for SFS analyses

#------------------------------Identifying Outlier Loci
# To do this, will need to download PGDspider and Bayescan

#----- PGDspider :

cd ~/bin
pwd
/projectnb/davies-hb/hannah/TVE_2bRAD/bin

wget http://www.cmpg.unibe.ch/software/PGDSpider/PGDSpider_2.0.7.1.zip
unzip PGDSpider_2.0.7.1.zip

#----- Bayescan :

cd ~/bin
pwd
/projectnb/davies-hb/hannah/TVE_2bRAD/bin

wget http://cmpg.unibe.ch/software/BayeScan/files/BayeScan2.1.zip
unzip BayeScan2.1.zip
cp BayeScan2.1/binaries/BayeScan2.1_linux64bits bayescan
chmod +x bayescan
rm -r -f BayeScan*


#============= Bayescan: looking for Fst outliers
# Converting vcf (using PGDspider) to Bayescan format: 
pwd
/projectnb/davies-hb/hannah/TVE_2bRAD/tufts_custom_pipeline_files/denovo_nosym

# make tab-delimited file called bspops LISTING assignments of individuals (as they are named in the vcf file) to populations, for example:
ind1	pop0
ind2	pop0
ind3	pop1
ind4	pop1

# create a file called vcf2bayescan.spid containing this text:
nano vcf2bayescan.spid

echo "############
# VCF Parser questions
PARSER_FORMAT=VCF
# Do you want to include a file with population definitions?
VCF_PARSER_POP_QUESTION=true
# Only input following regions (refSeqName:start:end, multiple regions: whitespace separated):
VCF_PARSER_REGION_QUESTION=
# What is the ploidy of the data?
VCF_PARSER_PLOIDY_QUESTION=DIPLOID
# Only output following individuals (ind1, ind2, ind4, ...):
VCF_PARSER_IND_QUESTION=
# Output genotypes as missing if the read depth of a position for the sample is below:
VCF_PARSER_READ_QUESTION=
# Take most likely genotype if "PL" or "GL" is given in the genotype field?
VCF_PARSER_PL_QUESTION=true
# Do you want to exclude loci with only missing data?
VCF_PARSER_EXC_MISSING_LOCI_QUESTION=false
# Select population definition file:
VCF_PARSER_POP_FILE_QUESTION=./bspops
# Only output SNPs with a phred-scaled quality of at least:
VCF_PARSER_QUAL_QUESTION=
# Do you want to include non-polymorphic SNPs?
VCF_PARSER_MONOMORPHIC_QUESTION=false
# Output genotypes as missing if the phred-scale genotype quality is below:
VCF_PARSER_GTQUAL_QUESTION=
# GESTE / BayeScan Writer questions
WRITER_FORMAT=GESTE_BAYE_SCAN
# Specify which data type should be included in the GESTE / BayeScan file  (GESTE / BayeScan can only analyze one data type per file):
GESTE_BAYE_SCAN_WRITER_DATA_TYPE_QUESTION=SNP
############" >vcf2bayescan.spid

# launching bayescan (Misha says this might take 12-24 hours so submitting as a job)
head bayescan

#!/bin/bash
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N bayescan # job name, anything you want
#$ -l h_rt=48:00:00 #maximum run time
#$ -M hannahaichelman@gmail.com #your email
#$ -m be

java -Xmx49152m -Xms512m -jar /projectnb/davies-hb/hannah/TVE_2bRAD/bin/PGDSpider_2.0.7.1/PGDSpider2-cli.jar -inputfile sfilt_forbayescan.vcf -outputfile Best.bayescan -spid vcf2bayescan.spid 

/projectnb/davies-hb/hannah/TVE_2bRAD/bin/bayescan Best.bayescan -threads=20

# then submit job:
qsub -pe omp 28 bayescan


# use bayescan_plots.R to examine results

# removing outliers from VCF file in order to re-calculate Fst without outliers (neutral sites only)
# still have to use sfilt here!

#Best.baye_fst.txt is your bayescan output
cut -d" " -f2- Best.baye_fst.txt | tail -n +2> aaa
#find line number in vcf file that contains header
grep -n "#CHROM" sfilt_forbayescan.vcf
#13, add 1
tail --lines=+14 sfilt_forbayescan.vcf | cut -f 1,2 | paste --delimiters "\t" - aaa > baye_fst_pos.txt

awk '$5<0.5 {print $1"\t"$2}' ./baye_fst_pos.txt > bayesOut4FST

# bayeOuts4FST is the list of outlier loci
grep -Fwvf bayesOut4FST sfilt.vcf > NeutralLoci.vcf
grep -Fwf bayesOut4FST sfilt.vcf > OutlierLoci.vcf


zgrep "^[^#]" NeutralLoci.vcf | awk '{print $1,$2}' >NeutralLoci.sites
zgrep "^[^#]" OutlierLoci.vcf | awk '{print $1,$2}' >OutlierLoci.sites


#------------------------------Analysis of Genetic Divergence Again - Neutral Loci Only
pwd
/projectnb/davies-hb/hannah/TVE_2bRAD/tufts_custom_pipeline_files/denovo_nosym

# with all clones removed, we have 50 samples

module load angsd/0.923

# collecting and indexing filter-passing sites but without outliers
# cut -f 1,2 sfilt.nooutlier.vcf > Sites_noOutliers
# remove the header lines from the Sites_noOutliers file before indexing
angsd sites index NeutralLoci.sites

export GENOME_REF=/projectnb/davies-hb/hannah/TVE_2bRAD/tufts_custom_pipeline_files/denovo_nosym/cdh_alltags_cc.fasta
TODO="-doSaf 1 -doMajorMinor 1 -doMaf 1 -doPost 1 -anc $GENOME_REF -ref $GENOME_REF"


angsd -sites NeutralLoci.sites -b bams_L1 -GL 1 -P 1 $TODO -out L1_neutral
angsd -sites NeutralLoci.sites -b bams_L2 -GL 1 -P 1 $TODO -out L2_neutral
angsd -sites NeutralLoci.sites -b bams_L3 -GL 1 -P 1 $TODO -out L3_neutral

# generating per-population SFS
realSFS L1_neutral.saf.idx >L1_neutral.sfs
realSFS L2_neutral.saf.idx >L2_neutral.sfs
realSFS L3_neutral.saf.idx >L3_neutral.sfs

# writing down 2d-SFS priors - L1 vs L2
realSFS L1_neutral.saf.idx L2_neutral.saf.idx -P 1 > p12_neutral.sfs ; realSFS fst index L1_neutral.saf.idx L2_neutral.saf.idx -sfs p12_neutral.sfs -fstout p12_neutral

# global Fst between populations
realSFS fst stats p12_neutral.fst.idx

#output:
	-> Assuming idxname:p12_neutral.fst.idx
	-> Assuming .fst.gz file: p12_neutral.fst.gz
	-> FST.Unweight[nObs:771324]:0.023573 Fst.Weight:0.136383
0.023573 0.136383

# The difference between unweighted and weighted values is averaging versus ratio of sums method. 
# The weighted value is derived from the ratio of the separate per-locus sums of numerator and denominator values, 
# while the unweighted value is the average of per-locus values. [If that is unclear: weighted is sum(a)/sum(a+b), while unweighted is average(a/(a+b))].


# writing down 2d-SFS priors - L1 vs L3
realSFS L1_neutral.saf.idx L3_neutral.saf.idx -P 1 > p13_neutral.sfs ; realSFS fst index L1_neutral.saf.idx L3_neutral.saf.idx -sfs p13_neutral.sfs -fstout p13_neutral

# global Fst between populations
realSFS fst stats p13_neutral.fst.idx

#output:
	-> Assuming idxname:p13_neutral.fst.idx
	-> Assuming .fst.gz file: p13_neutral.fst.gz
	-> FST.Unweight[nObs:761562]:0.325370 Fst.Weight:0.143440
0.325370 0.143440

# writing down 2d-SFS priors - L2 vs L3
realSFS L2_neutral.saf.idx L3_neutral.saf.idx -P 1 > p23_neutral.sfs ; realSFS fst index L2_neutral.saf.idx L3_neutral.saf.idx -sfs p23_neutral.sfs -fstout p23_neutral

# global Fst between populations
realSFS fst stats p23_neutral.fst.idx

#output:
	-> Assuming idxname:p23_neutral.fst.idx
	-> Assuming .fst.gz file: p23_neutral.fst.gz
	-> FST.Unweight[nObs:761536]:0.178447 Fst.Weight:0.100000
0.178447 0.100000


#------------------------------Analysis of Genetic Divergence Again - Outlier Loci Only
pwd
/projectnb/davies-hb/hannah/TVE_2bRAD/tufts_custom_pipeline_files/denovo_nosym

# with all clones removed, we have 50 samples

module load angsd/0.923

# remove the header lines from the Sites_noOutliers file before indexing
angsd sites index OutlierLoci.sites

export GENOME_REF=/projectnb/davies-hb/hannah/TVE_2bRAD/tufts_custom_pipeline_files/denovo_nosym/cdh_alltags_cc.fasta
TODO="-doSaf 1 -doMajorMinor 1 -doMaf 1 -doPost 1 -anc $GENOME_REF -ref $GENOME_REF"


angsd -sites OutlierLoci.sites -b bams_L1 -GL 1 -P 1 $TODO -out L1_outliers
angsd -sites OutlierLoci.sites -b bams_L2 -GL 1 -P 1 $TODO -out L2_outliers
angsd -sites OutlierLoci.sites -b bams_L3 -GL 1 -P 1 $TODO -out L3_outliers

# generating per-population SFS
realSFS L1_outliers.saf.idx >L1_outliers.sfs
realSFS L2_outliers.saf.idx >L2_outliers.sfs
realSFS L3_outliers.saf.idx >L3_outliers.sfs

# writing down 2d-SFS priors - L1 vs L2
realSFS L1_outliers.saf.idx L2_outliers.saf.idx -P 1 > p12_outliers.sfs ; realSFS fst index L1_outliers.saf.idx L2_outliers.saf.idx -sfs p12_outliers.sfs -fstout p12_outliers

# global Fst between populations
realSFS fst stats p12_outliers.fst.idx

#output:
	-> Assuming idxname:p12_outliers.fst.idx
	-> Assuming .fst.gz file: p12_outliers.fst.gz
	-> FST.Unweight[nObs:351]:0.645093 Fst.Weight:0.732279
0.645093 0.732279

# The difference between unweighted and weighted values is averaging versus ratio of sums method. 
# The weighted value is derived from the ratio of the separate per-locus sums of numerator and denominator values, 
# while the unweighted value is the average of per-locus values. [If that is unclear: weighted is sum(a)/sum(a+b), while unweighted is average(a/(a+b))].


# writing down 2d-SFS priors - L1 vs L3
realSFS L1_outliers.saf.idx L3_outliers.saf.idx -P 1 > p13_outliers.sfs ; realSFS fst index L1_outliers.saf.idx L3_outliers.saf.idx -sfs p13_outliers.sfs -fstout p13_outliers

# global Fst between populations
realSFS fst stats p13_outliers.fst.idx

#output:
	-> Assuming idxname:p13_outliers.fst.idx
	-> Assuming .fst.gz file: p13_outliers.fst.gz
	-> FST.Unweight[nObs:342]:0.659837 Fst.Weight:0.713172
0.659837 0.713172

# writing down 2d-SFS priors - L2 vs L3
realSFS L2_outliers.saf.idx L3_outliers.saf.idx -P 1 > p23_outliers.sfs ; realSFS fst index L2_outliers.saf.idx L3_outliers.saf.idx -sfs p23_outliers.sfs -fstout p23_outliers

# global Fst between populations
realSFS fst stats p23_outliers.fst.idx

#output:
	-> Assuming idxname:p23_outliers.fst.idx
	-> Assuming .fst.gz file: p23_outliers.fst.gz
	-> FST.Unweight[nObs:342]:0.466759 Fst.Weight:0.603811
0.466759 0.603811


#------------------------------DOWNSTREAM ANALYSES

# Used angsd_ibs_pca.R to plot hierarchical clustering of samples based on identity by state matrix (file = myresult2.noclone.ibsMat)
# in addition to principal coordinate analysis (PCoA) also using the same ibsMat. 
# Used admixturePlotting_v5.R to plot results of NgsAdmix (file = mydata.noclone_k3.qopt)


