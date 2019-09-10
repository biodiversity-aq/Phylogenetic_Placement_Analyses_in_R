# script to pool processed sequence data into one dataset
#==============================================================
# Author Maxime Sweetlove 2019
# lisence CC 4.0

#--------------------------------------------------------------
# initiation
#--------------------------------------------------------------
#libraries to process the sequences
library(vegan)
library(ggplot2)

source("/Users/msweetlove/OneDrive_RBINS/mARS_NewSeqData/R_wDir/Utils_SelectDataMars.R")

# load the table with the selected datasets
metaDataDB_sub <- read.csv("/Users/msweetlove/OneDrive_RBINS/mARS_NewSeqData/R_wDir/selected_Data_Bacteria.csv")

# set working directory
seqdatadir="/Users/msweetlove/OneDrive_RBINS/mARS_NewSeqData/R_wDir/Data"
setwd(seqdatadir)



#--------------------------------------------------------------
# 1 gathering the project data
#--------------------------------------------------------------
# get the names of BioProjects with data
BioPrj <- list.files(seqdatadir)
BioPrj <- BioPrj[grep(".fasta", BioPrj)]
BioPrj <- unname(sapply(BioPrj, FUN=function(x){gsub(".fasta", "",x)}))

BioPrj_OTUtab <- list.files(seqdatadir)
BioPrj_OTUtab <- BioPrj_OTUtab[grep(".csv", BioPrj_OTUtab)]
BioPrj_OTUtab <- BioPrj_OTUtab[grep("PRJ", BioPrj_OTUtab)]


#--------------------------------------------------------------
# 2 pooling ASV tables
#--------------------------------------------------------------
OTUtab<-data.frame()
OTUtax<-data.frame()
for(p_csv in BioPrj_OTUtab){
  print(p_csv)
  data_csv <- read.csv(p_csv, header=TRUE, row.names=1)
  data_tax <- data_csv[,colnames(data_csv) %in% c("kingdom", "phylum", "class", "order", "family", "genus")]
  data_tab <- data_csv[,!colnames(data_csv) %in% c("kingdom", "phylum", "class", "order", "family", "genus"),drop=FALSE]
  OTUtax <- rbind(OTUtax, data_tax)
  OTUtab <- combine.data.frame(OTUtab, data_tab, fill=0)
}

#--------------------------------------------------------------
# 3 collecting metadata
#--------------------------------------------------------------
metaDataDB <- metaDataDB_sub[metaDataDB_sub$geneticAccessionNumber %in% colnames(OTUtab),]

# only Antarctic samples
metaDataDB <- metaDataDB[metaDataDB$decimalLatitude <= -60,]

#table(metaDataDB$env_biome)
mar <- c("Ocean", "ocean", "seawater")
aqu <- c("polar lake")
ter <- c("Polar desert", "polar desert")

N_mar <- as.character(metaDataDB[metaDataDB$env_biome %in% mar,]$geneticAccessionNumber)
N_aqu <- as.character(metaDataDB[metaDataDB$env_biome %in% aqu,]$geneticAccessionNumber)
N_ter <- as.character(metaDataDB[metaDataDB$env_biome %in% ter,]$geneticAccessionNumber)

metaDataDB <- metaDataDB[metaDataDB$geneticAccessionNumber %in% c(N_mar, N_aqu, N_ter),]

# restructure metadata
metaDataDB$ASVnumID <- NA
BioPrj_fasta <- list.files(seqdatadir)
BioPrj_fasta <- BioPrj_fasta[grep(".fasta", list.files(seqdatadir))]

for(f in BioPrj_fasta){
  BioPrj <- gsub(".fasta", "", f) 
  if(BioPrj %in% names(table(droplevels(metaDataDB$BioprojectNumber)))){
    con <- file(f,"r")
    first_line <- readLines(con,n=1)
    close(con)
    dataset_num <-as.numeric(strsplit(first_line, "_")[[1]][2])
    metaDataDB[metaDataDB$BioprojectNumber == BioPrj, c("ASVnumID")] <- dataset_num
  }
}

# make a new table with more logical names and structure
meta <- data.frame(ASVnum_ID = as.numeric(metaDataDB$ASVnumID), 
                   mARS_ID = as.character(metaDataDB$datasetID),
                   sample_ID = as.character(metaDataDB$geneticAccessionNumber),
                   fileName = as.character(metaDataDB$datasetName),
                   sampleName = as.character(metaDataDB$eventName),
                   locality = as.character(metaDataDB$locality),
                   date = as.character(metaDataDB$eventDate),
                   targetGroup = as.character(metaDataDB$targetGroup),
                   marker = as.character(metaDataDB$marker),
                   markerSubfragment = as.character(metaDataDB$markerSubfragment),
                   seqMethod = as.character(metaDataDB$seqMethod),
                   studyNumber = as.character(metaDataDB$studyNumber),
                   BioprojectNumber = as.character(metaDataDB$BioprojectNumber),
                   primerSequenceForward = as.character(metaDataDB$primerSequenceForward),
                   primerSequenceReverse = as.character(metaDataDB$primerSequenceReverse),
                   decimalLatitude = as.numeric(metaDataDB$decimalLatitude),
                   decimalLongitude = as.numeric(metaDataDB$decimalLongitude),
                   env_biome = as.character(metaDataDB$env_biome),
                   env_feature = as.character(metaDataDB$env_feature),
                   env_material = as.character(metaDataDB$env_material),
                   depth = as.numeric(metaDataDB$depth),
                   elev = as.numeric(metaDataDB$elev),
                   conduc = as.numeric(metaDataDB$conduc),
                   salinity = as.numeric(metaDataDB$salinity),
                   temp = as.numeric(metaDataDB$temp),
                   truncQ = metaDataDB$truncQ,
                   trunclen = metaDataDB$trunclen,
                   maxee = metaDataDB$maxee,
                   maxN = metaDataDB$maxN,
                   trimleft = metaDataDB$trimleft,
                   numASV = metaDataDB$numASV,
                   lib_reads_seqd = metaDataDB$lib_reads_seqd,
                   numSeq_pocessed = metaDataDB$numSeq_pocessed,
                   stringsAsFactors=FALSE)


# remove samples from the tax and tab tables that did not make the cut  
OTUtab <- OTUtab[,colnames(OTUtab) %in% meta$sample_ID]
OTUtab <- OTUtab[rowSums(OTUtab) > 0,] #remove empty OTUs

OTUtax <- OTUtax[rownames(OTUtax) %in% rownames(OTUtab),]

# remove Eukaryote OTUs
OTUtax <- OTUtax[!(OTUtax$kingdom %in% "Eukaryota"),]
OTUtab <- OTUtab[rownames(OTUtab) %in% rownames(OTUtax),]


#--------------------------------------------------------------
# 4 the sequences
#--------------------------------------------------------------
#BiocManager::install("Rsamtools")
library(Rsamtools)

OTUfasta <-  readDNAStringSet(BioPrj_fasta[1])
for (f in BioPrj_fasta[-1]){
  OTUfasta <- c(OTUfasta, readDNAStringSet(f)) 
} 
OTUfasta <- OTUfasta[rownames(OTUtab)]


#--------------------------------------------------------------
# 5 writing data
#--------------------------------------------------------------
# MMA short for Mars Meta Analysis
write.csv(OTUtab, paste(seqdatadir, "MMA_OTUtab.csv", sep="/"))
write.csv(OTUtax, paste(seqdatadir, "MMA_OTUtax.csv", sep="/"))
write.csv(meta, paste(seqdatadir, "MMA_metaData.csv", sep="/"))

writeXStringSet(OTUfasta, paste(seqdatadir, "MMA_OTUseq.fasta", sep="/"), format="fasta")

  