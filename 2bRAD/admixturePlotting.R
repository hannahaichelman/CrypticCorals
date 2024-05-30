#setwd('~/Dropbox/BU/TVE/TVE_Github/DielTempVariability/')

# assembling the input table
dir="2bRAD/data_files/" # path to input files
inName="mydata.noclone_k3.qopt" # name of the input file to plot, output of ngsAdmix or ADMIXTURE run
#pops="inds2pops" # 2-column tab-delimited table of individual assignments to populations; must be in the same order as samples in the bam list or vcf file.

#------------

npops=as.numeric(sub(".+(\\d+)\\..+","\\1",inName))
tbl=read.table("2bRAD/data_files/mydata.noclone_k3.qopt")
i2p=read.table("2bRAD/data_files/bam_barcode_names_tuftscustom_noclone.csv", header = TRUE, sep = ",")
i2p = i2p[,1:2]

names(i2p)=c("ind","pop")
tbl=cbind(tbl,i2p)
row.names(tbl)=tbl$ind
head(tbl)
# putting populaitons in desired order (edit pop names as needed or skip to plot them alphabetically)
tbl$pop=factor(tbl$pop,levels=c("I2","I3","I4","O2","O3","O4"))

# re-name levels of factor to two letter site code
levels(tbl$pop) <- c("SP","CI","PD","BS","CA","BN")

cols_lineage_k2 <- c("#3f007d", "#807dba")
cols_lineage_k3 <- c("#bcbddc", "#807dba", "#3f007d")

#Download below from https://github.com/z0on/2bRAD_denovo/blob/master/plot_admixture_v5_function.R
source("2bRAD/data_files/plot_admixture_v5_function.R")
quartz()
ords=plotAdmixture(data=tbl,npops=npops,grouping.method="distance",vshift=0.1,colors=cols_lineage_k3)

# recording cluster affiliations
cluster.admix=apply(tbl[,1:npops],1,function(x) {return(paste(which(x>0.25),collapse=".")) })
#save(cluster.admix,file=paste(inName,"_clusters.RData",sep=""))

cluster.admix.df = as.data.frame(cluster.admix)
