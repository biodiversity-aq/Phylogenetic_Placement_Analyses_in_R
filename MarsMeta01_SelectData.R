# script to select appropriate datasets from the mARS database
#==============================================================
# Author Maxime Sweetlove 2019
# lisence CC 4.0


#--------------------------------------------------------------
# initiation
#--------------------------------------------------------------
#install.packages("BiocManager")
#BiocManager::install("Biostrings")

#libraries select the data
library(Biostrings)
library(seqinr)
library(seqRFLP)
library(stringr)
library(DescTools)
library(SRAdb)

setwd("/Users/msweetlove/OneDrive_RBINS/mARS_NewSeqData/R_wDir")
#setwd("/Users/msweetlove/OneDrive - Royal Belgian Institute of Natural Sciences/mARS_NewSeqData/R_wDir")

#--------------------------------------------------------------
# create working dataset containing all data from mARS
#--------------------------------------------------------------
mARS_data_dir="/Users/msweetlove/OneDrive_RBINS/mARS_NewSeqData/done"
source("/Users/msweetlove/OneDrive_RBINS/mARS_NewSeqData/R_wDir/Utils_SelectDataMars.R")
source("/Users/msweetlove/OneDrive_RBINS/mARS_NewSeqData/R_wDir/Utils_Alignments.R")
metaDataDB <- create.dataset(SeqSetPath = mARS_data_dir)

#--------------------------------------------------------------
# select datasets that agree with some proposed criteria
#--------------------------------------------------------------
# Bacteria 16S-v4
#--------------------------------------------------------------
#selecting criteria:
#  target region has to be 16S
metaDataDB_sub<-metaDataDB[metaDataDB$marker %in% c("16S", "16S rRNA", "16S SSU rRNA", "16S ssu rRNA"),]

# the sequences are available
metaDataDB_sub<-metaDataDB_sub[metaDataDB_sub$submitted_to_insdc == TRUE,]
metaDataDB_sub<-metaDataDB_sub[metaDataDB_sub$lib_reads_seqd > 10,]
metaDataDB_sub<-metaDataDB_sub[metaDataDB_sub$investigation_type == "mimarks-survey",]


#  metabarcode must contain the v4 region (based on primer sequence)
RefSeq <- read.fasta(file='/Users/msweetlove/OneDrive_RBINS/mARS_NewSeqData/R_wDir/eColiRef.fasta')
RefSeq <- c(toupper(paste(unlist(RefSeq)[[2]], sep='', collapse = '')))
primers <- data.frame(fprimers=metaDataDB_sub$primerSequenceForward,
                      rprimers=metaDataDB_sub$primerSequenceReverse)
primers<-unique(primers)
primers$events <- primers$p_start <- primers$p_end <- primers$v4 <- NA
for(rw in 1:nrow(primers)){
  eventIDs<-metaDataDB_sub[metaDataDB_sub$primerSequenceForward == primers[rw,]$fprimers &
                    metaDataDB_sub$primerSequenceReverse == primers[rw,]$rprimers,]$eventName
  primers[rw,]$events<-list(eventIDs)
  
  fprimers<-split_DNAstringList_unknownSep(primers[rw,]$fprimer)
  rprimers<-split_DNAstringList_unknownSep(primers[rw,]$rprimer)
  startpos <- endpos <- c()
  for(fprimer in fprimers){
    startpos <- c(startpos, primerPosition(fprimer, RefSeq, start=TRUE))
  }
  for(rprimer in rprimers){
    endpos <- c(endpos, primerPosition(revComp(rprimer), RefSeq, start=FALSE))
  }

  S <- primers[rw,]$p_start <- min(startpos)
  E <- primers[rw,]$p_end <- max(endpos) 
  
  if((S < 488 & E > 746) | #v4 is part of amplicon
     (S > 488 & E < 746 & E-S > 150) | #amplicon is part of v4
     (S < 488 & E < 746 & E > 640) | #amplicon overlaps v4 at start
     (S > 488 & E > 746 & S < 600)){ #amplicon overlaps v4 at end
    primers$v4 <- TRUE
    } else{primers$v4 <- FALSE}
}

v4_samples <- unlist(primers$events)
metaDataDB_sub<-metaDataDB_sub[metaDataDB_sub$eventName %in% v4_samples,]

## Finalize
write.csv(metaDataDB_sub, 
          "/Users/msweetlove/OneDrive_RBINS/mARS_NewSeqData/R_wDir/selected_Data_Bacteria.csv",
          row.names=FALSE)


#--------------------------------------------------------------
# Eukaryotes 18S-v4
#--------------------------------------------------------------