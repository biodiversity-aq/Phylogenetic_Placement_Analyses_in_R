# script that aligns ASVs to a reference alignment

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

# libraries
library(ggplot2)
library(reshape)
library(muscle)
library(DECIPHER)
library(ape) 
library(treeio)
library(ggtree)
library(ips) #phylogenetic analyses
library(phytools)

# load data, important libraries and working variables
suppressWarnings(source("/Users/msweetlove/OneDrive_RBINS/mARS_NewSeqData/R_wdir/MarsMeta01_LoadData.R"))

#load the mothur recreated SILVA SEED reference alignment, release 132 
SILVA_ref<-readDNAStringSet(file.path(CorePath, "R_wdir/silva.seed_v132.align"))

#--------------------------------------------------------------
# align sequences with Silva SEED as reference
#--------------------------------------------------------------
taxon <- "Alphaproteobacteria"
taxon_root <- "Deltaproteobacteria"

# get taxon out of reference (+2 seq that serve to root the tree)
  phy_subsetNames <- names(SILVA_ref)
  phy_subsetNames_root <- phy_subsetNames[grepl(paste("^.*;",taxon_root,";.*$", sep=""), phy_subsetNames)]
  phy_subsetNames_root <- phy_subsetNames_root[1]
  phy_subsetNames <- phy_subsetNames[grepl(paste("^.*;",taxon,";.*$", sep=""), phy_subsetNames)]
  phy_subsetNames <- c(phy_subsetNames, phy_subsetNames_root)
  SILVA_ref_sub <- SILVA_ref[phy_subsetNames]
  # change names to family_genus_species
    s.names <- c(names(SILVA_ref_sub))
    s.names <- sapply(s.names, function(x){strsplit(x, '\t')[[1]][3]})
    s.names <- sapply(s.names, function(x){paste(tail(strsplit(x, ';')[[1]], n=3), collapse='_')})
    s.names <- sapply(s.names, function(x){gsub("_    ", "_sp", x, fixed=T)})
    s.names2 <- make.names(s.names, unique=TRUE) #make each name unique
    names(s.names2)<-names(s.names)
    names(SILVA_ref_sub)<-s.names2[phy_subsetNames]
    
  # remove references of same genus if >5 occurences
    taxToKeep <- c(); taxList <- c()
    for(i in 1:length(names(SILVA_ref_sub))){
      taxName <- gsub("\\.[0-9]*$", "", names(SILVA_ref_sub)[i])
      if(!taxName %in% taxList){
        taxToKeep <- c(taxToKeep, names(SILVA_ref_sub)[i])
        taxList <- c(taxList, taxName)
      } else if(table(taxList[taxList==taxName])<6){
        taxToKeep <- c(taxToKeep, names(SILVA_ref_sub)[i])
        taxList <- c(taxList, taxName)
      }
    }
    SILVA_ref_sub <- SILVA_ref_sub[taxToKeep]

  # delete columns with only gaps for all sequences
  SILVA_ref_sub<-RemoveGaps(SILVA_ref_sub, "common")
  
  writeXStringSet(SILVA_ref_sub, file.path(CorePath, "R_wdir/silva_ref.temp.fasta"),
                  format="fasta")
  
  
# get the taxon sequences based on Wang classifier
  phy_subsetNames_OTU <- row.names(OTUdata_vi[OTUdata_vi$phylum==taxon | 
                                            OTUdata_vi$class==taxon |
                                            OTUdata_vi$order==taxon,])
  Seqdata_sub <- Seqdata_vi[phy_subsetNames_OTU]
  writeXStringSet(Seqdata_sub, file.path(CorePath, "R_wdir/seqdata.temp.fasta"),
                  format="fasta")
  
# align the sequences with reference (= profile alignment)
  seq_align<-muscle::muscle(c(SILVA_ref_sub, Seqdata_sub), profile=T,
                       in1=file.path(CorePath, "R_wdir/silva_ref.temp.fasta"),
                       in2=file.path(CorePath, "R_wdir/seqdata.temp.fasta"))
  
  writeXStringSet(as(seq_align,"DNAStringSet"), file.path(CorePath, "R_wdir/alignment.temp.fasta"),
                  format="fasta")
  
  
  
#--------------------------------------------------------------
# make a tree of the alignment
#-------------------------------------------------------------- 
  
# instalaion of FastTree, only run once:
  # system(paste("gcc -O3 -finline-functions -funroll-loops -Wall -o /Applications/FastTree /Users/msweetlove/Desktop/FastTree.c -lm"))

# making a reference tree with FastTree
  system(paste("OMP_NUM_THREADS=1 /Applications/FastTree -n 10 -nt -gtr -gamma",
               "<", file.path(CorePath, "R_wdir/silva_ref.temp.fasta"), 
               ">", file.path(CorePath, "R_wdir/RefTree.temp.nwk")))

#--------------------------------------------------------------
# root the tree
#-------------------------------------------------------------- 
  # read the reference phylogenetic tree  
  phyloTreeRef <- ape::read.tree(file.path(CorePath, "R_wdir/RefTree.temp.nwk"))

  # root phylogenetic tree with outgroup = taxon_root
  # get correct tip.labels of outgroup
  phy_subsetNames_root2 <- strsplit(phy_subsetNames_root, '\t')
  phy_subsetNames_root2 <- paste(tail(strsplit(phy_subsetNames_root2[[1]], ';')[[3]], n=3), collapse='_')
  phy_subsetNames_root2 <- as.vector(gsub("_    ", "_sp", phy_subsetNames_root2, fixed=T))
  
  # root tree
  phyloTreeRef_root <- ape::root(phyloTreeRef, phy_subsetNames_root2, resolve.root=TRUE)
  is.rooted(phyloTreeRef_root)
  ape::write.tree(phyloTreeRef_root, file = file.path(CorePath, "R_wdir/RefTree.root.nwk"),
             tree.names = TRUE)
  
#--------------------------------------------------------------
# place the amplicon sequences in the reference tree
#-------------------------------------------------------------- 
# downloading and setting up EPA-ng => just do this once
  # see https://github.com/Pbdas/epa-ng
  # run in Terminal tab (NOT R Console); run with alt+cmd+enter
  # cd /Applications
  # conda install -c bioconda epa-ng
  ## or with Homebrew: brew install brewsci/bio/epa-ng
  # epa-ng -h
  
  # in R console
  seq_align_comp<-readDNAStringSet(file.path(CorePath, "R_wdir/alignment.temp.fasta"))
  Seqdata_subsetNames <- names(Seqdata_sub)
  seq_align_sub <- seq_align_comp[Seqdata_subsetNames]
  writeXStringSet(seq_align_sub, file.path(CorePath, "R_wdir/alignment_sub.temp.fasta"),
                  format="fasta")
  
  # running EPA-ng (Terminal)
  export Rwd = "/Users/msweetlove/OneDrive_RBINS/mARS_NewSeqData/R_wdir"
  epa-ng --tree $Rwd/RefTree.root.nwk --ref-msa $Rwd/silva_ref.temp.fasta --query $Rwd/alignment_sub.temp.fasta --outdir $Rwd --redo --model GTR+G --baseball-heur
  # output called epa_result.jplace
  
#--------------------------------------------------------------
# working with the epa_result.jplace file
#-------------------------------------------------------------- 
  tree_place <- treeio::read.jplace(file.path(CorePath, "R_wdir/epa_result.jplace"))
  tree.placements <- as.data.frame(treeio::get.placements(tree_place, by = "best"))
  
  tree.phylo <- get.tree(tree_place)

    for(i in 1:length(tree.placements$name)){
      a<-tree.phylo$edge[tree.placements[i,]$edge_num,][1]
      b<-tree.placements[i,]$node
      c <- tree.placements[i,]$name
      d <- paste(OTUdata_vi[as.character(tree.placements[i,]$name),c("class", "family", "genus")], collapse="; ")
      
      print(paste(i, c, a, b, d))
    }
    
    otu.pnum<-10
    btree<-bind.tip(tree.phylo,
                    as.character(tree.placements[otu.pnum,]$name),
                    where=which(base::grepl(paste("^", tree.placements[otu.pnum,]$node, "$", sep=""),
                                            tree.phylo$edge[,2])))
    
    btree<-bind.tip(tree.phylo,
                    as.character(tree.placements[otu.pnum,]$name),
                    where=949)
    

      # subset tree with treeio
  bi_subset <- treeio::tree_subset(btree, 
                                   as.character(tree.placements[otu.pnum,]$name), 
                                   levels_back = 1)
  
  #plot data
  ggtree(bi_subset, aes(color = group)) + 
    geom_tiplab() + 
    theme_tree2() + 
    xlim(0, 100)

  
  
  
  
  tree.phylo$node.label
  

  #temporary files code:
  #tmpFile <- tempfile()
  #cat(tree_place@treetext, file=tmpFile)



  
  # collapse poorly suported nodes (optional)
  phyloTreeRef_root2<-di2multi(phyloTreeRef_root, tol=0.01)
  
  # subset tree with treeio
  seqToKeep <- names(Seqdata_sub)
  bi_subset <- treeio::tree_subset(phyloTree_root2, seqToKeep[3], levels_back = -3)
  
  paste(OTUdata_vi[seqToKeep[3],c("class", "order", "family", "genus", "species")], 
        collapse="; ")
  
  #plot data
  ggtree(bi_subset, aes(color = group)) + 
    geom_tiplab() + 
    theme_tree2() + 
    scale_color_manual(values = c(`1` = "red", `0` = "black")) +
    xlim(0, 4)
  
  #plot data
  bi_subset2<-treeio::tree_subset(phyloTreeRef_root2, phyloTreeRef_root2$tip.label[7], 
                                  levels_back = -4)
  ggtree(bi_subset2, aes(color = group)) + 
    geom_tiplab() + 
    theme_tree2() + 
    scale_color_manual(values = c(`1` = "red", `0` = "black")) +
    xlim(0, 4)
  
  
  
  
  
  library(BoSSA)
  tree_place <- read_jplace(file.path(CorePath, "R_wdir/epa_result.jplace"), full = TRUE)

  plot(tree_place,type="number",main="number",cex.number=1.5)
  
  plot(tree_place,type="fattree",main="fattree")
  plot(tree_place,type="precise",transfo=function(X){X*20})
  pplace_table <- pplace_to_table(tree_place,type="best")

  test<-pplace_to_matrix(tree_place)


  
  # delete the temp files
  #file.remove(file.path(CorePath, "R_wdir/silva_ref.temp.fasta"))
  #file.remove(file.path(CorePath, "R_wdir/seqdata.temp.fasta"))
  
  
#--------------------------------------------------------------
# Analysis: deliminating phylogenetic clades
#--------------------------------------------------------------

  #library(abind)
  #library(caper) # to get the Branch * sites matrices
  #library(betapart) # to compute Beta-Diversity matrices
  #library(ecodist) # to compute correlations profiles (function 'MRM')
  #library(bigmemory)
  #library(adephylo)



#--------------------------------------------------------------
# Analysis: geographic distribution
#--------------------------------------------------------------



#--------------------------------------------------------------
# Analysis: distribution over environmental biomes
#--------------------------------------------------------------



#--------------------------------------------------------------
# alpha diversity according to habitat and primer used
#--------------------------------------------------------------
# multiple linear regression models
# number of ASVs = primer (f, r and interaction? number of v-regions?) + habitat + number of seqs


  
  
#--------------------------------------------------------------
# Bayesian phylogenetic tree
#--------------------------------------------------------------

  Sys.setenv(PATH = paste(Sys.getenv("PATH"), "<path-to-mrbayes", sep = ":"))
  ips::mrbayes(as.DNAbin(seq_align), file=file.path(CorePath, "R_wdir/bayes_tree.nex"), 
              run=FALSE, nruns=10, ngen=10000, rates="gamma")
  
  phyloTree <- ape::read.nexus(file.path(CorePath, "R_wdir/bayes_tree.nex"))
  phyloTree <- treeio::read.mrbayes(file.path(CorePath, "R_wdir/bayes_tree.nex"))
  


