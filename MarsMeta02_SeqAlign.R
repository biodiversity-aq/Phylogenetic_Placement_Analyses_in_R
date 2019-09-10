# script that aligns ASVs to a reference alignment
# Author Maxime Sweetlove 2019
# lisence CC 4.0

#--------------------------------------------------------------
# initiation
#--------------------------------------------------------------
#if(!requireNamespace("BiocManager", quietly = TRUE)){
#  install.packages("BiocManager")
#}
#BiocManager::install("ggtree")
#BiocManager::install("ips")

#to install treeio
#source("https://bioconductor.org/biocLite.R")
## biocLite("BiocUpgrade") ## you may need this
#biocLite("treeio")
#biocLite("annotate")


# libraries
library(ggplot2)
library(ape) 
library(treeio)
library(ggtree)
library(ips) #phylogenetic analyses
library(phytools)
library(seqinr)

# set working directory
seqdatadir="/Users/msweetlove/OneDrive_RBINS/mARS_NewSeqData/R_wDir/Data"
CorePath="/Users/msweetlove/OneDrive_RBINS/mARS_NewSeqData"
setwd(seqdatadir)

# load subsetting- and support functions
source("/Users/msweetlove/OneDrive_RBINS/mARS_NewSeqData/R_wDir/Utils_SelectDataMars.R")
source("/Users/msweetlove/OneDrive_RBINS/mARS_NewSeqData/R_wdir/Utils_Phylogenies.R")

# load data
OTUtab <- read.csv(paste(seqdatadir, "MMA_OTUtab.csv", sep="/"), row.names=1, header=TRUE)
OTUtax <- read.csv(paste(seqdatadir, "MMA_OTUtax.csv", sep="/"), row.names=1, header=TRUE)
meta <- read.csv(paste(seqdatadir, "MMA_metaData.csv", sep="/"), row.names=1, header=TRUE)
OTUfasta <- readDNAStringSet(paste(seqdatadir, "MMA_OTUseq.fasta", sep="/"))

#load the mothur recreated SILVA SEED reference alignment, release 132 
SILVA_ref<-readDNAStringSet(file.path(CorePath, "R_wdir/silva.seed_v132.align"))

#--------------------------------------------------------------
# align sequences with Silva SEED as reference
#--------------------------------------------------------------
# used taxon and outgroup pairs:
sort(table(OTUtax$phylum), decreasing=TRUE)
taxon <- "Cyanobacteria"; taxon_root <- "Firmicutes"
taxon <- "Actinobacteria"; taxon_root <- "Chloroflexi"
taxon <- "Proteobacteria"; taxon_root <- "Chloroflexi"
#--------------------------------------------------------------
# make a reference tree of a specific phylum
reftree <- make.referenceTree(SILVA_ref, taxon=taxon, outgroup=taxon_root)
reftree_info <- reftree[[1]]
reftree_file <- reftree[[2]]
reftree_align <- reftree[[3]]

#--------------------------------------------------------------
# place the amplicon sequences in the reference tree
#-------------------------------------------------------------- 
  # in R console
  # get the taxon sequences based on Wang classifier (step to reduce computing time)
  Seqdata_sub <- get.taxon(fasta = OTUfasta, taxtab = OTUtax, 
                            taxon = taxon, taxrank=c("phylum", "class"))
  
  writeXStringSet(Seqdata_sub, file.path(seqdatadir, "seqdata.temp.fasta"),
                  format="fasta")
  
  # align the sequences with references using HMMR h
  # execute in Terminal: alt+cmd+enter 
  ### instalation 
    #cd Applications
    #conda install -c bioconda hmmer
  
  export workDir="/Users/msweetlove/OneDrive_RBINS/mARS_NewSeqData/R_wDir/Data"
  #refalign = the fasta alignment of the reference tree
  export refalign="/Users/msweetlove/OneDrive_RBINS/mARS_NewSeqData/R_wDir/Data/refTree.temp.align.fasta"
  hmmbuild $workDir/hmmRef.temp.hmm $refalign
  
  # seqdata = the query sequences (unaligned)
  export seqdata="/Users/msweetlove/OneDrive_RBINS/mARS_NewSeqData/R_wdir/Data/seqdata.temp.fasta"
  hmmalign -o $workDir/alignment.temp.sto --dna --mapali $refalign $workDir/hmmRef.temp.hmm $seqdata
  
  esl-reformat -o $workDir/alignment.temp.align.fasta afa $workDir/alignment.temp.sto
  rm $workDir/alignment.temp.sto
  
  # in console: split alignment.temp.fasta into the reference an querry sequences
  seq_align_comp<-readDNAStringSet(file.path(seqdatadir, "alignment.temp.align.fasta"))
  # including the dots as gaps seem to cause an error in EPA-NG, replace them with dashes
  seq_align_comp <- DNAStringSet(sapply(seq_align_comp, function(x){gsub(".", "-", x, fixed=T)}))
  #subset the aligned querry sequences out of the alignment
  seq_align_sub <- seq_align_comp[names(Seqdata_sub)]
  writeXStringSet(seq_align_sub, file.path(seqdatadir, "querySeq.temp.align.fasta"),
                  format="fasta")
  #subset the aligned reference sequences out of the alignment
  seq_align_sub2 <- seq_align_comp[setdiff(names(seq_align_comp), names(Seqdata_sub))]
  writeXStringSet(seq_align_sub2, file.path(seqdatadir, "referenceSeq.temp.align.fasta"),
                  format="fasta")
  
  ### downloading and setting up EPA-ng via Terminal: alt+cmd+enter
    # see https://github.com/Pbdas/epa-ng
    # cd /Applications
    # conda install -c bioconda epa-ng
    # epa-ng -h
  
  reftree_info #show the fasttree model parameters 
  # running EPA-ng (Terminal)
  # reftree = the reference tree in nwk format, same as refTree.tree
  export reftree="/Users/msweetlove/OneDrive_RBINS/mARS_NewSeqData/R_wDir/Data/refTree.temp.tree.nwk"
  # refalign2 = the alignment of the reference sequences
  export refalign2="/Users/msweetlove/OneDrive_RBINS/mARS_NewSeqData/R_wDir/Data/referenceSeq.temp.align.fasta"
  # the alignment of the query sequences, length-compatible and aligned to the references (without the references)
  export queryalign="/Users/msweetlove/OneDrive_RBINS/mARS_NewSeqData/R_wDir/Data/querySeq.temp.align.fasta"
  epa-ng --tree $reftree --ref-msa $refalign2 --query $queryalign --outdir $workDir --redo --model "GTR{0.7294/1.7603/1.1684/0.7102/3.2658/1.000}" --baseball-heur
  
  # output is called epa_result.jplace
  # rename the output to prevent it from being overwritten #in console
  jplaceout <- file.path(seqdatadir, "epa_result.jplace")
  new_japlace <- file.path(CorePath, paste("R_wdir/EPA_", taxon,".jplace", sep=""))
  paste("EPA_", taxon, ".jplace", sep="")
  file.copy(from = jplaceout, to = new_japlace, overwrite = TRUE)
  
#####################################################
## test: check how the placement looks
#####################################################
  # read the placement data back to R
  jplaceData <- read.jplace2(jplace_file=new_japlace)
  placements <- jplaceData@placements
  placements <- get.placements.best(placements)
  
  # add the placements in the reference tree, based on their "best" location
  source("/Users/msweetlove/OneDrive_RBINS/mARS_NewSeqData/R_wdir/Utils_Phylogenies.R")
  placement_tree<-get.placements.phylo(jplaceData, verbose=TRUE)
  
  
  placements2 <- jplaceData@placements
  placements2 <- get.placements.best(placements2)
  placements2<- placements2[placements2$edge_num>119 & placements2$edge_num<126,]
  #pp<- placements2[placements2$edge_num ==121,][1,]
  #placements2<- placements2[placements2$edge_num !=121,]
  #placements2<-rbind(placements2, pp)
  #placements2<- placements2[placements2$edge_num ==121,][1,]
  jplaceData2 <- jplaceData
  jplaceData2@placements<-placements2
  source("/Users/msweetlove/OneDrive_RBINS/mARS_NewSeqData/R_wdir/Utils_Phylogenies.R")
  placement_tree2<-get.placements.phylo(jplaceData2, verbose=TRUE)
  checkValidPhylo(placement_tree2)
  
  # just to be sure, check if the resuling tree is in the correct format
  checkValidPhylo(placement_tree)
  
  # subset tree
  otu_pnum <- 13
  p_otuName <- as.character(placements[otu_pnum,]$name)
  bi_subset <- subset.tree(placement_tree, p_otuName, ancestorNodes_back=3)

  # check with Wang classifier taxonomic annotation
  OTUtax[row.names(OTUtax)==p_otuName,c("phylum", "class", "order", "family", "genus")]
  
  # plot subset tree
  colvec<-attributes(bi_subset)$group
  colvec<-gsub(0, "black", colvec); colvec<-gsub(1, "red", colvec)
  par(mar=c(0,0,0,0))
  plot(bi_subset, edge.color=colvec,tip.color=colvec,
       edge.width=as.numeric(attributes(bi_subset)$group)+1,
       cex=0.7)


  


