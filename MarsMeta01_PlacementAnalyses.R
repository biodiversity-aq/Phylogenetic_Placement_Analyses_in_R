# script that aligns ASVs to a reference alignment
# Author Maxime Sweetlove 2019
# lisence CC 4.0

#--------------------------------------------------------------
# initiation
#--------------------------------------------------------------
# libraries
library(ggplot2)
library(ape) 
library(treeio)
library(ggtree)
library(ips) #phylogenetic analyses
library(phytools)
library(seqinr)
library(tidyverse)
library(adephylo)


# set working directory
seqdatadir="/Users/msweetlove/OneDrive_RBINS/mARS_NewSeqData/R_wDir/Data"
CorePath="/Users/msweetlove/OneDrive_RBINS/mARS_NewSeqData"
wdir="/Users/msweetlove/OneDrive_RBINS/mARS_NewSeqData/R_wDir"
setwd(wdir)

# load subsetting- and support functions
source("/Users/msweetlove/OneDrive_RBINS/mARS_NewSeqData/R_wDir/Utils_SelectDataMars.R")
source("/Users/msweetlove/OneDrive_RBINS/mARS_NewSeqData/R_wdir/Utils_Phylogenies.R")

# load data
OTUtab <- read.csv(paste(seqdatadir, "MMA_OTUtab.csv", sep="/"), row.names=1, header=TRUE)
OTUtax <- read.csv(paste(seqdatadir, "MMA_OTUtax.csv", sep="/"), row.names=1, header=TRUE)
meta <- read.csv(paste(seqdatadir, "MMA_metaData.csv", sep="/"), row.names=1, header=TRUE)

#--------------------------------------------------------------
# load jplace data
taxon <- "Cyanobacteria"
taxon <- "Actinobacteria"
taxon <- "Proteobacteria"

jplace <- file.path(CorePath, paste("R_wdir/EPA_", taxon,".jplace", sep=""))
jplaceData <- read.jplace2(jplace_file=jplace)

#--------------------------------------------------------------
# preparing the jplace file
#-------------------------------------------------------------- 
# get the placements
placements <- jplaceData@placements %>%
  get.placements.best()
# add the placements in the reference tree, based on their "best" location
placement_tree<-get.placements.phylo(jplaceData, verbose=TRUE)
# just to be sure, check if the resuling tree is in the correct format
checkValidPhylo(placement_tree)

#--------------------------------------------------------------
# visualizing the complete tree
#-------------------------------------------------------------- 
bi_subset <- rownames(OTUtax[OTUtax$phylum==taxon,]) %>%
  make.groupvec(placement_tree, .)
colvec <- attributes(bi_subset)$group %>%
  gsub(0, "black", .) %>%
  gsub(1, "red", .)
par(mar=c(0,0,0,0))
plot(bi_subset, edge.color=colvec,
     edge.width=((as.numeric(attributes(bi_subset)$group)-2)*-1)+0.3,
     show.tip.label = FALSE)

#--------------------------------------------------------------
# visualizing the position of specific taxonomic groups in the tree
#-------------------------------------------------------------- 
groupNames<-names(table(droplevels(OTUtax[OTUtax$phylum==taxon,]$order)));groupNames
target_taxon<-groupNames[24]
p_otuName<-rownames(OTUtax[OTUtax$phylum==taxon & OTUtax$order==target_taxon,])
ref_taxon <- placement_tree$tip.label %>%
  grep(gsub("ales", "", target_taxon), .) %>%
  placement_tree$tip.label[.]
if(identical(ref_taxon, character(0))){
  bi_subset<-make.groupvec(placement_tree, p_otuName)
} else{
  bi_subset<-make.groupvec(placement_tree, p_otuName, ref_taxon)
}
colvec <- attributes(bi_subset)$group %>%
  gsub(0, "black", .) %>%
  gsub(1, "red", .) %>%
  gsub(2, "green", .)
par(mar=c(0,0,0,0))
plot(bi_subset, edge.color=colvec,
     edge.width=((as.numeric(attributes(bi_subset)$group)-1)*2)+0.2,
     show.tip.label = FALSE)

#### of subgroup: look at those ASVs most distantly related from the references
tree_dist <- adephylo::distTips(placement_tree)
tree_dist2 <- tree_dist %>%
  as.matrix() %>%
  .[p_otuName, ref_taxon] %>%
  apply(., 1, FUN=min)
dev.off()
hist(tree_dist2, breaks=100)

s_names <- names(tree_dist2[tree_dist2>0.25])
bi_subset<-make.groupvec(placement_tree, s_names, ref_taxon)
colvec <- attributes(bi_subset)$group %>%
  gsub(0, "black", .) %>%
  gsub(1, "red", .) %>%
  gsub(2, "green", .)
par(mar=c(0,0,0,0))
plot(bi_subset, edge.color=colvec,
     edge.width=((as.numeric(attributes(bi_subset)$group)-1)*2)+0.2,
     show.tip.label = FALSE)
# get clostest reference
OTUtax[s_names[1],]
phytools::getSisters(placement_tree, s_names[1], mode=c("label"))

#--------------------------------------------------------------
# OTUs versus to habitat
#--------------------------------------------------------------
## distribution marine non marine in the phylogeny
meta$biom <- meta$env_biome %>%
  gsub("Ocean", "mar", .) %>%
  gsub("ocean", "mar", .) %>%
  gsub("seawater", "mar", .) %>%
  gsub("Polar desert", "ter", .) %>%
  gsub("polar desert", "ter", .) %>%
  gsub("polar lake", "aqu", .)

bi_subset <- rownames(OTUtax[OTUtax$phylum==taxon,])
marOTUs <- c()
terOTUs <- c()
aquOTUs <- c()
marterOTUs <- c()
maraquOTUs <- c()
teraquOTUs <- c()
cosmOTUs <- c()
for(tax in bi_subset){
  sampNms<- OTUtab %>%
    dplyr::filter(row.names(.)==tax)
  sampNms<-sampNms[,colSums(sampNms)>0] %>%
    colnames()
  if(is.null(sampNms)){next}
  biom <- meta %>% 
    dplyr::filter(sample_ID %in% sampNms) %>%
    dplyr::select(biom)
  if("ter" %in% biom & !"mar" %in% biom & !"aqu" %in% biom){biom <- "ter"
  } else if(!"ter" %in% biom & "mar" %in% biom & !"aqu" %in% biom){biom <- "mar"
  } else if(!"ter" %in% biom & !"mar" %in% biom & "aqu" %in% biom){biom <- "aqu"
  } else if("ter" %in% biom & "mar" %in% biom & !"aqu" %in% biom){biom <- "marter"
  } else if(!"ter" %in% biom & "mar" %in% biom & "aqu" %in% biom){biom <- "maraqu"
  } else if("ter" %in% biom & !"mar" %in% biom & "aqu" %in% biom){biom <- "teraqu"
  } else if("ter" %in% biom & "mar" %in% biom & "aqu" %in% biom){biom <- "cosm"
  } 
  if(biom=="ter"){
    terOTUs <- c(tax, terOTUs)
  } else if(biom=="mar"){
    marOTUs <- c(tax, marOTUs)
  } else if(biom=="aqu"){
    aquOTUs <- c(tax, aquOTUs)
  } else if(biom=="marter"){
    marterOTUs <- c(tax, marterOTUs)
  } else if(biom=="teraqu"){
    teraquOTUs <- c(tax, teraquOTUs)
  } else if(biom=="maraqu"){
    maraquOTUs <- c(tax, maraquOTUs)
  } else if(biom=="cosm"){
    cosmOTUs <- c(tax, cosmOTUs)
  }
}

tree_biom <- make.groupvec(placement_tree, marOTUs, terOTUs, aquOTUs)
colvec <- attributes(tree_biom)$group %>%
  gsub(0, "black", .) %>% ## no group
  gsub(1, "blue", .) %>% ## marine
  gsub(2, "green", .) %>% ## terrestrial
  gsub(3, "red", .) ## inland waters
par(mar=c(0,0,0,0))
plot(tree_biom, edge.color=colvec,
     edge.width=0.3,
     show.tip.label = FALSE)

# next part: compare clade grouping to habitat grouping??


# multiple linear regression models
# number of ASVs = primer (f, r and interaction? number of v-regions?) + habitat + number of seqs



#--------------------------------------------------------------
# phylogeography
#-------------------------------------------------------------- 
## geographic distance vs. phylogenetic decay
library(geosphere)
## phylo distance
phylo_dist <- adephylo::distTips(placement_tree) %>%
  as.matrix()
phyNms <- rownames(OTUtab)

phylo_dist_sub <- phylo_dist[rownames(phylo_dist) %in% phyNms, colnames(phylo_dist) %in% phyNms]

geo_dist <- matrix(ncol=ncol(phylo_dist_sub), nrow=nrow(phylo_dist_sub), data=0)
colnames(geo_dist) <- colnames(phylo_dist_sub) 
rownames(geo_dist) <- rownames(phylo_dist_sub) 

# calculate geo distance for each OTU pair in the phylo_dist_sub matrix
# when occuring in multiple samples, the minimal distance is chosen
### need to look for a shorter method here, takes a loooong time to compute
for(i in 1:nrow(geo_dist)){
  require(geosphere)
  tax_i <- rownames(geo_dist)[i]
  sampNms_i <- OTUtab %>%
    dplyr::filter(row.names(.)==tax_i)
  sampNms_i<-sampNms_i[,colSums(sampNms_i)>0] %>%
    colnames()
  if(is.null(sampNms_i)){next}
  for(j in 1:ncol(geo_dist)){
    tax_j <- colnames(geo_dist)[j]
    sampNms_j<- OTUtab %>%
      dplyr::filter(row.names(.)==tax_j)
    sampNms_j<-sampNms_j[,colSums(sampNms_j)>0] %>%
      colnames()
    if(is.null(sampNms_j)){next}
    dist_vec<-c()
    for(d_i in sampNms_i){
      d1 <- unlist(meta[meta$sample_ID==d_i,c("decimalLongitude", "decimalLatitude")][1,])
      for(d_j in sampNms_j){
        d2 <- unlist(meta[meta$sample_ID==d_j,c("decimalLongitude", "decimalLatitude")][1,])
        dv <- geosphere::distm(d1, d2, fun = distGeo)
        dist_vec <- c(dist_vec, dv)
      }
    }
    geo_dist[i,j] <- min(dist_vec) 
  }
}


# plot the distance matrices against each other
plot(as.vector(geo_dist), as.vector(phylo_dist_sub), pch=19, cex=0.5)


phylo_dist_sub <- phylo_dist[rownames(phylo_dist) %in% terOTUs, 
                             colnames(phylo_dist) %in% terOTUs]
geo_dist_sub <- geo_dist[rownames(geo_dist) %in% terOTUs, 
                             colnames(geo_dist) %in% terOTUs]
plot(as.vector(geo_dist_sub), as.vector(phylo_dist_sub), pch=19, cex=0.5)


cor.test(log(upper.tri(geo_dist)), log(upper.tri(phylo_dist_sub)))




#--------------------------------------------------------------
# visualizing the position of specific ASVs in the tree
#-------------------------------------------------------------- 
# subset tree
otu_pnum <- 43
p_otuName <- as.character(placements[otu_pnum,]$name)
#p_otuName<-"ASV_18_01"
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

ggplot(as.dendrogram(bi_subset))







#--------------------------------------------------------------
# Analysis: phylogenetic clades
#--------------------------------------------------------------
## are there special clades for antarctica that are poorly represented in the references?
## what are the dominant terrestrial polar/tundra, marine and freshwater clades?

SILVA_ref<-readDNAStringSet(file.path(CorePath, "R_wdir/silva.seed_v132.align"))
ref_taxa <- names(SILVA_ref); SILVA_ref<-NULL
ref_taxa2 <- sapply(ref_taxa, function(x){strsplit(x, '\t')[[1]][3]})
ref_taxa2 <- sapply(ref_taxa2, function(x){paste(tail(strsplit(x, ';')[[1]], n=3), collapse='_')})
ref_taxa2 <- sapply(ref_taxa2, function(x){strsplit(x, ';')[[1]]})
ref_taxa <- data.frame(ncol=6)
colnames(ref_taxa)<-c("kingdom", "phylum", "class", "order", "family", "genus")
for(i in 1:length(ref_taxa2)){
  for(j in 1:6){
    ref_taxa[i,j]<-ref_taxa2[[i]][j]
  }
}


## what is the percentage of sequences with 100% id with references?
##  graph: no difference, 1%, 5%,...
#plot branchlength distribution
groupNames<-row.names(OTUtax[OTUtax$phylum==taxon,])
bi_subset<-make.groupvec(placement_tree, groupNames)
par(mar=c(5,3,2,3))
otu_edgelengths <- bi_subset$edge.length[attributes(bi_subset)$group==1]
hist(otu_edgelengths)
















#--------------------------------------------------------------
# other ideas
#--------------------------------------------------------------
## unifrac versis biome
## unifrac versus region
## diversity latitude gradient?
## is there turnover in dominant groups wit increasing latitude








