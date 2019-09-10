# script to download and process sequence data
#==============================================================
# Author Maxime Sweetlove 2019
# lisence CC 4.0

#--------------------------------------------------------------
# initiation
#--------------------------------------------------------------
#install.packages("BiocManager")
#BiocManager::install("Biostrings")

# available taxonomy databases:
  #silva_nr_v128_train_set.fasta | gg_13_8_train_set_97.fasta | 
  #rdp_train_set_16.fasta | silva_nr_v132_train_set.fasta


#libraries to process the sequences
library(vegan)
library(ggplot2)
library(dada2)
library(gridExtra)
library(phyloseq)
library(DECIPHER)
library(phangorn)
library(SRAdb)

# load the table with the selected datasets
metaDataDB_sub <- read.csv("/Users/msweetlove/OneDrive_RBINS/mARS_NewSeqData/R_wDir/selected_Data_Bacteria.csv")
good_datasets <- unique(metaDataDB_sub$datasetID)
PRJ_IDs <- unique(metaDataDB_sub$BioprojectNumber)

# load additional info (for the dataset literature reference)
require(xlsx)
metaData_meta <-read.xlsx("/Users/msweetlove/OneDrive_RBINS/mARS_NewSeqData/MARS_Overview_PotentialDatasets.xlsx", sheetName = "Antarctic_microbial_studies")

# prepare metadata file to save new important informtion
metaDataDB_sub$maxee <- metaDataDB_sub$trunclen <- metaDataDB_sub$truncQ <- NA
metaDataDB_sub$trimleft <- metaDataDB_sub$maxN <- NA
metaDataDB_sub$numSeq_raw <- metaDataDB_sub$numSeq_pocessed <- metaDataDB_sub$numASV <- NA

seqdatadir="/Users/msweetlove/OneDrive_RBINS/mARS_NewSeqData/R_wDir/Data"
#dir.create(seqdatadir)
setwd(seqdatadir)

# download the latest SRA database ##library(SRAdb)
# last download was on: 30-07-2019
#SRAdb::getSRAdbFile(destdir=seqdatadir, destfile = "SRAmetadb.sqlite.gz")

##### for each dataset:
# length(PRJ_IDs)
## Done: PRJEB23732 (9); PRJNA355879 (10)
#### Failed: PRJNA208431 (not demultiplexed?), PRJEB1602 (not demultiplexed?)

### , PRJEB16346, PRJEB22851, PRJNA433310, PRJEB15586
### PRJEB11496, PRJEB27415, PRJNA344476, PRJNA415906, PRJNA401941, PRJNA305344

### failed on server: PRJNA395930, PRJNA244335, PRJNA421293
# remove samples:
currentSampleSet <- metaDataDB_sub$BioprojectNumber== as.character(PRJ_IDs[prj])
metaDataDB_sub<-metaDataDB_sub[!currentSampleSet,]

### succedeon server: PRJNA359740
prj <- 27
print(PRJ_IDs[prj])

which(as.character(droplevels(PRJ_IDs)) %in% "PRJEB16346")
#--------------------------------------------------------------
# 1 connect to INSDC and download the sequence data
#--------------------------------------------------------------
sample_numbers <- as.character(metaDataDB_sub[metaDataDB_sub$BioprojectNumber==as.character(PRJ_IDs[prj]),]$geneticAccessionNumber)
mth <- as.character(unique(metaDataDB_sub[metaDataDB_sub$BioprojectNumber==as.character(PRJ_IDs[prj]),]$seqMethod))
grp <- as.character(unique(metaDataDB_sub[metaDataDB_sub$BioprojectNumber==as.character(PRJ_IDs[prj]),]$targetGroup))
mkr <- as.character(unique(metaDataDB_sub[metaDataDB_sub$BioprojectNumber==as.character(PRJ_IDs[prj]),]$markerSubfragment))
num_s <- nrow(metaDataDB_sub[metaDataDB_sub$BioprojectNumber==as.character(PRJ_IDs[prj]),])

print(paste("sequencing method:", mth)); print(paste("target group:", grp));print(paste("v-region:", mkr)); print(paste(num_s, "samples in total"))
# the article:
as.character(metaData_meta[metaData_meta$INSDC_study %in% as.character(PRJ_IDs[prj]),]$publication)

  
# downloading the data
sra_dbname <- file.path(paste(getwd(), "SRAmetadb.sqlite", sep="/"))
sra_con <- dbConnect( dbDriver("SQLite"), sra_dbname )
getSRAfile(sample_numbers , sra_con = sra_con, destDir = getwd(), fileType = 'fastq', srcType = 'ftp')
dbDisconnect( sra_con )

#--------------------------------------------------------------
# 2 unzip files and remove compressed files
#--------------------------------------------------------------
fileList<- list.files(path=getwd(), pattern='*.fastq.gz')
for(file in fileList){
  fileName= strsplit(file, ".gz")[[1]]
  system2(command="gunzip", args=paste(" -c ", file," > ", fileName, sep=''), stdout=TRUE)
  system2(command="rm", args=file)
}

#--------------------------------------------------------------
# 3 quality controll
#--------------------------------------------------------------
print(mth)
FullPath <- getwd()

###############################################################
# A) for 454 data or single end Illumina
###############################################################
    fileList<- sort(list.files(path=getwd(), pattern='*.fastq'));fileList
    
    # extract file names 
    sampleNames = unlist(lapply(X=fileList, FUN=function(x){strsplit(x,".fastq")[[1]]}))
    
    # Specify the full path to the forward and reverse files
    rawseqs <- file.path(FullPath, fileList)
    
    # have a quick look at the data quality
    plotQualityProfile(rawseqs[2])
    
    # Specify the full path to the Quality Checked files
    QCs <- file.path(FullPath, paste(sampleNames, "_QC.fastq", sep=''))
    
    # quality filtering
    currentSampleSet <- metaDataDB_sub$BioprojectNumber== as.character(PRJ_IDs[prj])
    metaDataDB_sub[currentSampleSet,]$maxee <- m <- 23 #23
    metaDataDB_sub[currentSampleSet,]$trunclen <- t <- 250
    metaDataDB_sub[currentSampleSet,]$truncQ <- q <- 0
    metaDataDB_sub[currentSampleSet,]$trimleft <- l <- 0
    metaDataDB_sub[currentSampleSet,]$maxN <- n <- 0
    
    out <- filterAndTrim(rawseqs, QCs, truncLen= t, verbose=TRUE,
                         maxN= n, maxEE= m, truncQ= q, rm.phix=TRUE,
                         compress=FALSE, multithread=TRUE, trimLeft= l)
    
    out; plotQualityProfile(QCs[1])

    rawSeqs <- out[paste(metaDataDB_sub[currentSampleSet,]$geneticAccessionNumber , ".fastq", sep=""),1]
    metaDataDB_sub[currentSampleSet,]$numSeq_raw <- rawSeqs
      
    # dereplication: getting (temporarily) ride of redundant sequences
    derep <- derepFastq(QCs, verbose=TRUE)
    
    # Name the derep-class objects by the sample names
    names(derep) <- sampleNames
    
    # unsupervised learning step to estimate run-specific error parameters
    # first step of generating AVS
    err <- learnErrors(QCs, multithread=TRUE)
    plotErrors(err)

    # detect ASVs
    # OMEGA_A : p-value at which an ASV is accepted (the lower the more conservative, default 1e-40)
    # BAND_SIZE : Banding restricts the net cumulative number of insertion of one sequence relative to the other. 
    #             The default value of BAND_SIZE is 16. If DADA is applied to marker genes with high rates of indels, 
    #             such as the ITS region in fungi, the BAND_SIZE parameter should be increased. 
    dadas <- dada(derep, err=err, multithread=TRUE, pool=TRUE,
                   OMEGA_A=1e-70, BAND_SIZE=20, GREEDY=TRUE, MIN_ABUNDANCE=2)
    dadas[[1]]
    
    # make sequence table
    seqtabAll <- makeSequenceTable(dadas)
    
    # number of ASVs:
    ncol(seqtabAll)

###############################################################
# B) for Paired-end Illumina data
###############################################################
    # get list of forward and reverse files in the same order
    fileList_1 <- sort(list.files(path=getwd(), pattern='*_1.fastq'));fileList_1
    fileList_2 <- sort(list.files(path=getwd(), pattern='*_2.fastq'))
    
    # extract file names 
    sampleNames = unlist(lapply(X=fileList_1, FUN=function(x){strsplit(x,"_1.fastq")[[1]]}))
    
    # Specify the full path to the forward and reverse files
    rawseqFs <- file.path(FullPath, fileList_1)
    rawseqRs <- file.path(FullPath, fileList_2)
    
    # have a quick look at the data quality
    plotQualityProfile(c(rawseqFs[1], rawseqRs[1]))
    
    # Specify the full path to the forward and reverse Quality Checked files
    QCFs <- file.path(FullPath, paste(sampleNames, "_1_QC.fastq", sep=''))
    QCRs <- file.path(FullPath, paste(sampleNames, "_2_QC.fastq", sep=''))
    
    # quality filtering
    currentSampleSet <- metaDataDB_sub$BioprojectNumber== as.character(PRJ_IDs[prj])
    m <- c(23,23); metaDataDB_sub[currentSampleSet,]$maxee <- m[1]
    t <- c(240,240); metaDataDB_sub[currentSampleSet,]$trunclen <- t[1]
    metaDataDB_sub[currentSampleSet,]$truncQ <- q <- 0
    l <- c(15,15); metaDataDB_sub[currentSampleSet,]$trimleft <- l[1]
    metaDataDB_sub[currentSampleSet,]$maxN <- n <- 0
    
    out <- filterAndTrim(rawseqFs, QCFs, rawseqRs, QCRs, truncLen= t,
                         maxN= n, maxEE= m, truncQ= q, rm.phix=TRUE,
                         compress=FALSE, multithread=TRUE, trimLeft= l)
    
    out; plotQualityProfile(c(QCFs[1], QCRs[1]))
    
    samp_names <- paste(metaDataDB_sub[currentSampleSet,]$geneticAccessionNumber , "_1.fastq", sep="")
    rawSeqs <- out[samp_names,1]
    metaDataDB_sub[currentSampleSet,]$numSeq_raw <- rawSeqs
    
    # dereplication: getting (temporarily) ride of redundant sequences
    derepFs <- derepFastq(QCFs, verbose=TRUE)
    derepRs <- derepFastq(QCRs, verbose=TRUE)
    
    # Name the derep-class objects by the sample names
    names(derepFs) <- sampleNames
    names(derepRs) <- sampleNames
    
    # unsupervised learning step to estimate run-specific error parameters
    # first step of generating AVS
    errF <- learnErrors(QCFs, multithread=TRUE)
    plotErrors(errF)
    errR <- learnErrors(QCRs, multithread=TRUE)
    plotErrors(errR)
    
    # detect ASV
    # OMEGA_A : p-value at which an ASV is accepted (the lower the more conservative, default 1e-40)
    # BAND_SIZE : Banding restricts the net cumulative number of insertion of one sequence relative to the other. 
    #             The default value of BAND_SIZE is 16. If DADA is applied to marker genes with high rates of indels, 
    #             such as the ITS region in fungi, the BAND_SIZE parameter should be increased. 
    dadaFs <- dada(derepFs, err=errF, multithread=TRUE, pool=TRUE,
                   OMEGA_A=1e-70, BAND_SIZE=20, GREEDY=TRUE, MIN_ABUNDANCE=2)
    dadaRs <- dada(derepRs, err=errR, multithread=TRUE, pool=TRUE,
                   OMEGA_A=1e-70, BAND_SIZE=20, GREEDY=TRUE, MIN_ABUNDANCE=2)
    dadaFs[[1]]
    dadaRs[[1]]
    
    # merge sequence (paired-end assembly)
    mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)
    
    # make sequence table
    seqtabAll <- makeSequenceTable(mergers)

###############################################################
# Continuation, for all data
###############################################################
# number of ASVs:
ncol(seqtabAll)

#filter on lengths
seqlens <- nchar(getSequences(seqtabAll)); table(seqlens)
  # option 1: the lengths were trimed in the QC: do noting
  seqtabAllF <- seqtabAll
  # option 2: there is variation in lengths: remove the largest and smales sequences
  seqtabAllF <- seqtabAll[,seqlens >= 268 & seqlens <= 270]
table(nchar(getSequences(seqtabAllF)))

#remove chimeras
seqtabNoC <- removeBimeraDenovo(seqtabAllF)

# number of ASVs:
ncol(seqtabNoC)

#--------------------------------------------------------------
# 4 Taxonomic assignemnt
#--------------------------------------------------------------
### case of Bacteria
  # option 1 Greengenes
    #fastaRef <- "/Users/msweetlove/Documents/TaxClassificationDBs/gg_13_8_train_set_97.fasta"
    #taxTab <- assignTaxonomy(seqtabNoC, refFasta = fastaRef, multithread=TRUE, verbose=TRUE,
    #                     minBoot = 80)
  #option 2 Silva
    fastaRef <- "/Users/msweetlove/Documents/TaxClassificationDBs/silva_nr_v132_train_set.fasta"
    taxTab <- assignTaxonomy(seqtabNoC, refFasta = fastaRef, multithread=TRUE, verbose=TRUE,
                             minBoot = 80)
    
   
### case of Eukaryotes
    fastaRef <- "/Users/msweetlove/Documents/TaxClassificationDBs/pr2_v4_11_1.fasta"
    
    taxTab <- assignTaxonomy(seqtabNoC, refFasta = fastaRef, multithread=TRUE, verbose=TRUE, minBoot = 80, 
                         taxLevels = c("kingdom","superphylum","phylum","class","order","family",
                                       "genus","species")) #PR2 has different tax levels than default in dada


### continuation for all data
unname(head(taxTab))
#--------------------------------------------------------------
# 5 Finalize and save OTU table
#--------------------------------------------------------------
seqtabNoC_df<-data.frame(t(seqtabNoC))
taxTab_df<-data.frame(taxTab)

### case of Silva
colnames(taxTab_df)<-c("kingdom", "phylum", "class", "order", "family", "genus")

### case greengenes
    colnames(taxTab_df)<-c("kingdom", "phylum", "class", "order", "family", "genus", "species")
    taxTab_df$kingdom<-gsub("k__", "", taxTab_df$kingdom)
    taxTab_df$phylum<-gsub("p__", "", taxTab_df$phylum)
    taxTab_df$class<-gsub("c__", "", taxTab_df$class)
    taxTab_df$order<-gsub("o__", "", taxTab_df$order)
    taxTab_df$family<-gsub("f__", "", taxTab_df$family)
    taxTab_df$genus<-gsub("g__", "", taxTab_df$genus)
    taxTab_df$species<-gsub("s__", "", taxTab_df$species)
  

### case of PR2
    colnames(taxTab_df)<-c("kingdom", "superphylum", "phylum", "class", "order", "family", "genus", "species")

### continuation for all data
for(i in 1:ncol(taxTab_df)){
  if(is.factor(taxTab_df[,i])) levels(taxTab_df[,i]) <- c(levels(taxTab_df[,i]),"unclassified")
} 
taxTab_df[taxTab_df==""] <- c("unclassified")
taxTab_df[taxTab_df=="NA"] <- c("unclassified")
taxTab_df[is.na(taxTab_df)] <- c("unclassified")

OTUtab<-cbind(taxTab_df[match(rownames(seqtabNoC_df), rownames(taxTab_df)),], seqtabNoC_df)
OTUtab$RefSeq<-rownames(OTUtab)
OTUtab<-droplevels(OTUtab)

# rename the ASVs
prefix <- paste(rep("0", 2-nchar(prj)), prj, "_", sep='')
for(i in 1:nrow(OTUtab)){
  n_zero <- nchar(nrow(OTUtab))-nchar(i) #number of digits
  rownames(OTUtab)[[i]]<-paste("ASV_", prefix, paste(rep("0", n_zero), collapse=''), i, sep='')
}

# save data before further clean-up
write.csv(OTUtab, "Data_OTUtab.intermediate.csv")
table(OTUtab$kingdom)

### continuation for all data
  OTUtab<-OTUtab[rowSums(is.na(OTUtab))<ncol(OTUtab),]

# write reference sequences to fasta
OTUtab_fasta<-OTUtab[,colnames(OTUtab) %in% "RefSeq"]
names(OTUtab_fasta)<-rownames(OTUtab)

# rename OTUs and write fasta
prefix <- paste(rep("0", 2-nchar(prj)), prj, "_", sep='')
fastafile<-''
for(i in 1:nrow(OTUtab)){
  n_zero <- nchar(nrow(OTUtab))-nchar(i) #number of digits
  old_name<-rownames(OTUtab)[[i]]
  new_name<-paste("ASV_", prefix, paste(rep("0", n_zero), collapse=''), i, sep='')
  fastaseq<-OTUtab_fasta[names(OTUtab_fasta)==old_name][[1]]
  fastafile<-paste(fastafile, '>', new_name, '\n', fastaseq, '\n', sep='')
  rownames(OTUtab)[[i]]<-new_name
}

# final clean up
OTUtab <- OTUtab[,!colnames(OTUtab) %in% c("RefSeq")]

#save cleaned dataset and fasta file
write.table(fastafile, paste(as.character(PRJ_IDs[prj]), ".fasta", sep=""), quote=FALSE, row.names=FALSE, col.names = FALSE)
write.csv(OTUtab, paste(as.character(PRJ_IDs[prj]), ".csv", sep=""))

# save any changes made to the metadata
write.csv(metaDataDB_sub, 
          "/Users/msweetlove/OneDrive_RBINS/mARS_NewSeqData/R_wDir/selected_Data_Bacteria.csv",
          row.names=FALSE)

# delete the redundant files
list.files()
file.remove(list.files(path=getwd(), pattern='*.fastq')); list.files()


