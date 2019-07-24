# utility function
# July 2019
# author: Maxime Sweetlove
# NOTE: based on the packages ape and treeio

# functions:
# extract.placement()
# get.descendants()
# get.ancestors()
# get.placements.best()
# read.jplace2()
# get.placements.phylo()


library(ggplot2)
library(reshape)
library(muscle)
library(DECIPHER)
library(ape) 
library(treeio)
library(ggtree)
library(ips) #phylogenetic analyses
library(phytools)

jplace_file <- file.path(CorePath, "R_wdir/epa_result.jplace")

### test JPLACE
tree.text<-"(((((((A:4{0},B:4{1}):6{2},C:5{3}):8{4},D:6{5}):3{6},E:21{7}):10{8},((F:4{9},G:12{10}):14{11},H:8{12}):13{13}):13{14},((I:5{15},J:2{16}):30{17},(K:11{18},L:11{19}):2{20}):17{21}):4{22},M:56{23});"
Ref.Phylo<- read.tree(text=tree.text)
placements<-data.frame(name=c("AAA", "EEE", "XXX", "XXX", "XXX", "YYY", "outgroup", 
                              "AA2", "EE2", "XX2", "out2", "out3",
                              "AA3", "XX3"), 
                       edge_num=c(1, 5, 18, 8, 9, 25, 14,
                                  1, 5, 18, 14,
                                  1, 18, 14), 
                       likelihood=c(-4987.366, -4977.344, -4487.390, -4523.130, -4589.930, -4977.344, -4977.344,
                                    -4987.366, -4977.344, -4977.344, -4977.344,
                                    -4977.344, -4977.344, -4977.344),
                       like_weight_ratio=c(0.08358164, 0.04068357, 0.02859273, 0.02859273, 0.02859273, 0.02859273, 0.02859273,
                                           0.08358164, 0.04068357, 0.02859273, 0.02859273,
                                           0.02859273, 0.04068357, 0.04068357),
                       distal_length=c(0.5, 1, 1, 1, 1, 1, 1,
                                       1, 1, 1, 1,
                                       1, 1,1),
                       pendant_length=c(4, 5, 3, 3, 3, 5, 24,
                                        8, 12, 3, 30,
                                        7, 23,2))


info <- list("test tree")
test_jplace <- new("jplace2",
                treetext   = tree.text,
                phylo      = Ref.Phylo,
                placements = placements,
                info       = info,
                file       = normalizePath(jplace_file)
                )

#############################################################################
###    Classes and methodes
#############################################################################

#----------------------------------------------------------------------------
# create the S4 class "jplace2" and associated methods
#----------------------------------------------------------------------------

setClass("jplace2", slots=list(treetext="character", 
                                 phylo="phylo", 
                                 placements="data.frame",
                                 info="list",
                                 file="character")
  )

setMethod("show",
            "jplace2",
            function(object) {
              cat("a jplace2 class object\n")
              warning("show method still in devellopment")
            }
  )

#############################################################################
###    generic functions
#############################################################################
extract.placement <- function(object, phylo) {
  placements <- object$placements
  placements <- json_data$placements
  
  
  place <- placements[,1]
  
  ids <- NULL
  if (length(placements) == 2) {
    ids <- vapply(placements[,2], function(x) x[1], character(1))
    names(place) <- ids
  }
  
  place.df <- do.call("rbind", place)
  row.names(place.df) <- NULL
  if (!is.null(ids)) {
    nn <- rep(ids, vapply(place, function(x) {
      nr <- nrow(x)
      if (is.null(nr))
        return(1)
      return(nr)
    }, numeric(1)))
    place.df <- data.frame(name=nn, place.df)
    colnames(place.df) <- c("name", object$fields)
  } else {
    colnames(place.df) <- object$fields
  }
  edgeNum.df <- attr(phylo, "edgeNum")
  place.df <- merge(place.df, edgeNum.df, by.x = "edge_num", by.y = "edgeNum")
  as_tibble(place.df)
}

get.placements.best<-function(placements=NULL){
  if(!is.null(placements)){
    if (!'likelihood' %in% names(placements)){
      stop('no column called \"likelihood\"')
    }
    placements_best <- placements[0,]
    names<-as.character(unique(placements$name))
    for(name in names){
      placements_sub <- placements[placements$name==name,]
      placements_sub <- placements_sub[order(placements_sub$likelihood, decreasing = TRUE),]
      placements_best <- rbind(placements_best, placements_sub[1,])
      ## ToDo_LowPriority: dealing with possibilities that have equal likelihood
    }
  }else{
    stop('no placement data provided to get.placements.best()')
  }
  return(placements_best)
}

get.descendants <- function(target_node, edgeList, type=c("all", "tips", "nodes", "direct")){
  # function will get the descending nodes (in the direction away from the root) from a given node
  # input requires an phylo-style edgeList instead of a phylo object
  # type:all = all descendants
  # type:tips = only the tips
  # type:nodes = only internal nodes
  # type:direct = only the direct descendants
  descList<-edgeList[edgeList[,1]==target_node,2] #the direct descendants
  # get the root, if not present: stop function
  treeRoot <- tryCatch({
    setdiff(edgeList[,1], edgeList[,2])
  }, error=function(e){
    stop("there is no root in the edgeList. Are you sure the tree was rooted?")
  })

  if(is.null(type) | !(type %in% c("all", "tips", "nodes", "direct"))){
    stop("incorrect \"type\" argument. 
         argument must be either \"all\", \"tips\", \"nodes\", or \"direct\"")
  }
  
  if(type != "direct"){
    doneList<-c()
    n_added <- length(descList)
    while(! n_added == 0){ #loop through the tree to add all the other descendants
      n_added = 0
      for(d in descList){
        if(!d %in% doneList){
          d2 <- edgeList[edgeList[,1]==d,2]
          doneList <- c(doneList,d)
          if(length(d2)>0){
            descList<- c(descList, d2)
            n_added <- n_added +1
          }
        }
      }
    }
    if(type=="tips"){
      descList <- descList[descList < treeRoot]
    } else if(type=="nodes"){
      descList <- descList[descList > treeRoot]
    }
  }

  return(descList)
}

get.ancestors <- function(target_node, edgeList, type=c("all", "parent")){
  # function will get the ancestor nodes (in the direction to the root) from a given node
  # input requires an phylo-style edgeList instead of a phylo object
  # type:all = all ancestors
  # type:parent = only the direct ancestor
  ancsList<-edgeList[edgeList[,2]==target_node,1] #the parent
  
  if(is.null(type) | !(type %in% c("all", "parent"))){
    stop("incorrect \"type\" argument. 
         argument must be either \"all\", \"parent\"")
  } 
  if(is.na(target_node) | !(target_node %in% c(edgeList))){
    stop("incorrect \"target_node\" argument.")
  }
  
  if(type == "all"){
    doneList<-c()
    n_added <- length(ancsList)
    while(! n_added == 0){ #loop through the tree to add all the other descendants
      n_added = 0
      for(a in ancsList){
        if(!a %in% doneList){
          a2 <- edgeList[edgeList[,2]==a,1]
          doneList <- c(doneList,a)
          if(length(a2)>0){
            ancsList<- c(ancsList, a2)
            n_added <- n_added +1
          }
        }
      }
    }
  }
  return(ancsList)
}


#############################################################################
###    Reading jplace data
#############################################################################
require(jsonlite)
require(treeio)
read.jplace2 <- function(jplace_file=""){
  ## function to read a jplace file
  ## output will be a jplace2 class object
  
  # 1) read in the data
  json_data <- fromJSON(jplace_file)
  RefTree.jplace <- json_data$tree
    
  # check c.1: is the tree realy present?
  if (identical(RefTree.jplace, character(0))) {
    stop("the tree in the jplace file is an empty character string")
  }
    
  # 2) process the extended newick file
  ### ToDo_LowPriority check if {} is present, otherwise quit or use ape read.tree
  # 2.1) change the {} to another delimitor: @@@
  RefTree.jplace <- gsub("}","",gsub("{","@@@",RefTree.jplace,fixed=TRUE),fixed=TRUE)
  RefTree.jplace <- gsub("[ \t]", "", RefTree.jplace)
  # 2.2) split the tree sting into different parts
  RefTree.parts <- unlist(strsplit(RefTree.jplace, "(?=[\\(\\),;])", perl = TRUE))

  TreeData <- data.frame(label=NA, edge_num=NA, childNode=NA, parentNode=NA, 
                           branchLength=NA, is_tip=NA)
  # 2.3) examine the parts and extract the phylo-like data, keeping the original edge numbering
  # this part differs from the ape::read.tree() function, which will re-number the edges and ignore the original edge numbering
  N_runNode <- 1
  N_internalNode <- sum(RefTree.parts == ",") + 2
  for(i in 1:length(RefTree.parts)){
      #first identify the nodes
      #assume edge number = child node
      if(!RefTree.parts[i] %in% c("(", ")", ",", ";")){
        NodeInfo <- unlist(strsplit(RefTree.parts[i], ":"))
        tipName <- NodeInfo[1]
        NodeInfo2 <- unlist(strsplit(NodeInfo[2],"@@@"))
        edgeNum <- as.numeric(as.character(NodeInfo2[2]))
        
        TreeData[N_runNode,]$edge_num <- edgeNum
        TreeData[N_runNode,]$childNode <- edgeNum
        TreeData[N_runNode,]$branchLength <- as.numeric(as.character(NodeInfo2[1]))
        
        if(!tipName==""){ # it's a tip
          TreeData[N_runNode,]$label <- tipName
          TreeData[N_runNode,]$is_tip <- "a_tip"   
        } else{ # it's an internal branch
          TreeData[N_runNode,]$label <- sprintf("n_%03d", N_internalNode)
          TreeData[N_runNode,]$is_tip <- "c_internal"
          N_internalNode <- N_internalNode + 1
        }
        
        # find the parent node
        j <- i # j = position within RefTree.parts during this while-loop
        k <- 0 # k = number of parantheses opened along the way
        while(TRUE){
          if(RefTree.parts[j]=="("){
            j <- j + 1
            k <- k + 1
          } else if(RefTree.parts[j]==")" && k==0){
            break
            print(j)
          } else if(RefTree.parts[j]==")" && k!=0){
            j <- j + 1
            k <- k - 1
          } else{
            j <- j + 1
          }
        }
        ParentNode <- unlist(strsplit(RefTree.parts[j+1], ":"))
        ParentNode <- unlist(strsplit(ParentNode[2],"@@@"))
        TreeData[N_runNode,]$parentNode <- as.numeric(as.character(ParentNode[2]))

        N_runNode <- N_runNode + 1
    } 
  }
  row.names(TreeData) <- c(1:nrow(TreeData))
    
  if(0 %in% TreeData$edge_num){
    nodes.shifted <- 1
    TreeData$edge_num <- TreeData$edge_num + nodes.shifted
    TreeData$childNode <- TreeData$childNode + nodes.shifted
    TreeData$parentNode <- TreeData$parentNode + nodes.shifted
  } else{
    nodes.shifted <- 0
  }
  
  # 2.4) dealing with the root
  if(NA %in% TreeData$parentNode){
    # give higest edge number to rootEdge
    TreeData[is.na(TreeData$parentNode),]$parentNode <- rootEdge <- max(TreeData$childNode)+1
    } else{rootEdge=NULL}
  
  # 3) finish making the RefTree into a phylo object
  # 3.1) reorder the branches conform the ape package
  #   - the first n nodes correspond to the tips
  #   - node n+1 corresponds to the root
  #   - n+2 to m are the internal nodes
  TreeData_tp<-TreeData[order(TreeData$is_tip, decreasing=FALSE),]
  TreeData_tp$childNode_orig <- TreeData_tp$childNode
  TreeData_tp$parentNode_orig <- TreeData_tp$parentNode
  TreeData_tp$edge_num_orig <- TreeData_tp$edge_num
  
  newOrder_nodes<- c(1:nrow(TreeData_tp))
  row.names(TreeData_tp) <- TreeData_tp$childNode <- TreeData_tp$edge_num <- newOrder_nodes
  
  for(i in 1:nrow(TreeData_tp)){
    if(!TreeData_tp[i,]$parentNode_orig==max(TreeData_tp$parentNode_orig)){
      TreeData_tp[i,]$parentNode <- TreeData_tp[TreeData_tp$childNode_orig==TreeData_tp[i,]$parentNode_orig,]$childNode
    } else{
      TreeData_tp[i,]$parentNode <- TreeData_tp[i,]$parentNode_orig
    }
  }
  
  if(!is.null(rootEdge)){
    # in 2.4, rootEdge was given the higest node number, i.e m+1
    # needs to change to node n+1
    newrootplace <- length(which(TreeData_tp$is_tip=="a_tip"))+1
    # all child_nodes from newrootplace onward need to be +1
    TreeData_tp[TreeData_tp$childNode %in% newrootplace:rootEdge,]$childNode <- 1 + TreeData_tp[TreeData_tp$childNode %in% newrootplace:rootEdge,]$childNode
    # root needs to become newrootplace
    # all parent_nodes from newrootplace onward need to be +1
    for(i in 1:nrow(TreeData_tp)){
      if(TreeData_tp[i,]$parentNode >= newrootplace & 
         TreeData_tp[i,]$parentNode != rootEdge){
        TreeData_tp[i,]$parentNode <- TreeData_tp[i,]$parentNode + 1
      } else if(TreeData_tp[i,]$parentNode == rootEdge){
        TreeData_tp[i,]$parentNode <- newrootplace
      }
    }
  }
  
  
  # 3.2) use the information of the TreeData object to make a phylo object of the reference tree
  edge <- cbind(sapply(TreeData_tp[,c("parentNode")],function(x){x<-as.integer(x)}),
                sapply(TreeData_tp[,c("childNode")],function(x){x<-as.integer(x)}))
  
  edge.length <- as.vector(TreeData_tp$branchLength)
  tip.label <- TreeData_tp[TreeData_tp$is_tip=="a_tip",]$label
  Nnode <- sum(RefTree.parts == ")")
  Ref.Phylo <- list(edge = edge, edge.length = edge.length, Nnode = Nnode, 
                    tip.label = tip.label)
  class(Ref.Phylo) <- "phylo"
  attr(Ref.Phylo, "order") <- "cladewise"
  
    ## ToDo_HighPriority: output tree has issue with ape::ladderize()
    
  # 4) getting the phylogenetic placement data
  # 4.1) extracting the placements from JSON
    placements_full <- json_data$placements
    placements <- placements_full[,1]
    
  # 4.2) converting to a single dataframe and add the names
    for(i in 1:length(placements)){
      placements[[i]] <- cbind(rep(i,nrow(placements[[i]])),placements[[i]])
    }
    placements <- as.data.frame(do.call(rbind,placements),stringsAsFactors=FALSE)
    colnames(placements) <- c("name",json_data$fields)
    placements <- data.frame(apply(placements,2,as.numeric))

    if (length(placements_full) == 2) {
      PlaceNames<-unlist(placements_full[,2])
      placements$name <- sapply(placements$name, function(x){x<-PlaceNames[x]})
     } else{
      warning("check the placements: did not get the expected 2 lists, so no names are provided")
     }
    
    # 4.3) re-order columns to fit with the treeio placement format
    placements$edge_num <- placements$edge_num + nodes.shifted
    placements <- placements[,c("edge_num", "name", "likelihood", "like_weight_ratio",
                                "distal_length", "pendant_length")]
    placements$node<-NA
    
    # 4.4) convert the edge numbers to the new edge numbers in the reference tree (Ref.Phylo)
    for(i in 1:nrow(placements)){
      placements[i,]$edge_num <- TreeData_tp[TreeData_tp$childNode_orig==placements[i,]$edge_num,]$childNode
      placements[i,]$node <- TreeData_tp[TreeData_tp$childNode_orig==placements[i,]$edge_num,]$childNode
    }

    ##ToDo_LowPriority: convert to tibble and check comlete compatibility with treeio
    
    ##ToDo_LowPriority: put all the data into a treeio-style jplace object:
      #@placements  #placements
      #@file #normalizePath(jplace_file)
      #@treetext #json_data$tree
      #@phylo  #Ref.Phylo
      #@data
      #@extraInfo #matrix(nrow=0, ncol=0)
      #@tip_seq #NA
      #@anc_seq #NA
      #@seq_type #NA
      #@tipseq_file #NA
    

  # creating the jplace2 object for the output
  #### convert to phylo class of treeio package
    
  info <- c(json_data$metadata, version=json_data$version)
  data_out <- new("jplace2",
             treetext   = json_data$tree,
             phylo      = Ref.Phylo,
             placements = placements,
             info       = info,
             file       = normalizePath(jplace_file)
  )

  return(data_out)
}

#############################################################################
###    putting placements into a tree
#############################################################################

get.placements.phylo<-function(jplace2_data="", verbose=FALSE){
  # function that takes a jplace2 object
  # returns a phylo object of the reference tree that includes placements based on the "best" position
  # this cannot be done with bind.tree of ape because nodes are renumbered. 
  # Here, nodes are also renumbered, but their original osition is remembered internally
  
  # 1) checks before getting started and getting data
  if(is(jplace2_data,"jplace2")){    # check c.1: input must be of class jplace2
    Ref.Phylo <- jplace2_data@phylo
    placements <- jplace2_data@placements
    if(!is.null(placements)){    # check c.2: there must be placement data to work with
      placements<-get.placements.best(placements) #get the placements with the higest likelihood
    }else{
      stop('input does not contain placement data:
         could not execute get.placements.phylo()')
    }
  } else{stop('input must be of class jplace2:
         could not execute get.placements.phylo()')
  }
  
  print("started processing, this may take a while...")
  
  # 2) placing the placements in the reference tree
  # 2.1) build a matrices to keep track of the original edge numbers
  #   taking the following assumptions:
  #   - numbers 1:ntips are tips
  #   - number ntips+1 is the root (and the tree must be rooted)
  #   - numbers (ntips+1):nnodes are the internal nodes
  #   - negative values indicate new placements
  #   - edges have the same number as the child node they support
  n.tips <- length(Ref.Phylo$tip.label) # number of tips in the reference tree
  n.nodes <- nrow(Ref.Phylo$edge) # number of nodes in the reference tree
  tree_edges <- Ref.Phylo$edge # the phylo edge matrix that will be updated
  orig_tip <- data.frame(orig_edge = (1:n.tips), 
                         prev_edge = (1:n.tips), 
                         new_edge = (1:n.tips),
                         label = c(Ref.Phylo$tip.label),
                         stringsAsFactors=FALSE) # matrix keeping track of all the tips
  current_root <- previous_root <- orig_root <- n.tips + 1L # the root node
  orig_node <- data.frame(orig_edge = ((orig_root+1):(n.nodes+1)), 
                          prev_edge = ((orig_root+1):(n.nodes+1)), 
                          new_edge = ((orig_root+1):(n.nodes+1))) # matrix keeping track of all the internal nodes
  if(!is.null(Ref.Phylo$edge.length)){
    BL <- TRUE
    if(verbose){print("including branch lengths")}
    branch_lengths <- data.frame(orig_edge = c(1:(orig_root-1), (orig_root+1):(n.nodes+1)), 
                                 prev_edge = c(1:(orig_root-1), (orig_root+1):(n.nodes+1)),
                                 new_edge = c(1:(orig_root-1), (orig_root+1):(n.nodes+1)), 
                                 length = Ref.Phylo$edge.length) 
  } else{
    BL <- FALSE
    if(verbose){print("the reference tree has no edge lengths.")}
  }
  
  # remove double placements, keep just the ones with higest probability
  placements <- get.placements.best(placements)
  target_locations<-unique(placements$edge_num)
  if(verbose){print(paste("There are", nrow(placements), "placements for", 
                          length(target_locations), "edges"))}

  # function to update edges
  update.edge <- function(e){
    # updates an edge 
    # input:  e = an edge number to be updated
    # inherits following variebles from above:
    #         previous_root
    #         current_root
    #         p_loc_prev
    #         orig_tip
    #         orig_node
    # note: for the edge that will be split, of the 3 possible outcomes this function will give the name of the internal edge
    
    if(e == p_loc_prev){ 
      # special case: if it's the edge to be split
      if(e < previous_root){ # when it's a tip: return the (first) internal edge
        e_updated <- max(orig_node[orig_node$prev_edge==e,]$new_edge)
      } else if(e > previous_root){ # when it's an internal edge: return the edge furthest to the root (largest number)
        e_updated <- orig_node[orig_node$prev_edge==e & orig_node$orig_edge > 0,]$new_edge
      } else if(e == previous_root){
        e_updated <- max(orig_node$prev_edge) + 2
        }
    } else if(e < previous_root & length(orig_tip[orig_tip$prev_edge==e,]$new_edge) > 1){ 
      # special case: conflict with a new tip to an internal edge
      # return the edge that was in the tree first
      e_updated <- orig_tip[orig_tip$prev_edge==e & orig_tip$orig_edge != -p_loc_orig,]$new_edge
    } else{ 
      # any other case: the edge was already part of the tree in the previous itteration
      if(e < previous_root){ # case 1: changing a tip
        e_updated <- orig_tip[orig_tip$prev_edge==e,]$new_edge
      } else if(e == previous_root){ # case 2: changing the root
        e_updated <- current_root
      } else if(e > previous_root){ # case 3: changing an internal node
        e_updated <- orig_node[orig_node$prev_edge==e & orig_node$orig_edge != -p_loc_prev,]$new_edge
      } else{e_updated=NULL}
    }
    return(e_updated)
  }
  
  #---------------------------------------------
  # 2.2) running through the placements (p) and inserting them in the tree
  # approach taken: first look how many placements an edge will cary
  # if it is more than one, a subtree of the placements needs to be created first
  # the order of the subtree is based on the distal lengths, 
  # with the largest distal length splitting of first
  #---------------------------------------------
  for(p_loc_orig in target_locations){ 
    if(verbose){print(paste("processing placements for edge", p_loc_orig))}
    # p_edge is the target edge for thisparticular (set of) placement(s)
    p_placements <- placements[placements$edge_num==p_loc_orig,]
    
    # 2.2.1 assess the type of the placement 
    # there are three different cases were a placement can be addded: tip, node or root,
    # each one requires an individual treatment
    if(p_loc_orig < orig_root){  
      ### the target edge is a tip
      case_is_tip <- TRUE
      case_is_node <- case_is_root <- FALSE
    } else if(p_loc_orig > orig_root){
      ### the target edge is an internal node
      case_is_node <- TRUE
      case_is_tip <- case_is_root <- FALSE
    } else{
      ### the target edge is the root
      case_is_root <- TRUE
      case_is_node <- case_is_tip <- FALSE
    }
    
    # 2.2.2) check if there are multiple placements for one edge
    n_placements <- nrow(p_placements)
    
    # 2.2.3) any addition will shift the current root down by the number of placements
    previous_root <- current_root
    current_root <- current_root + n_placements
    
    # 2.2.4) to add any number of new tips/nodes, the previous settings must be saved
    orig_tip$prev_edge <- orig_tip$new_edge
    orig_node$prev_edge <- orig_node$new_edge
    if(BL){branch_lengths$prev_edge <- branch_lengths$new_edge}
   
    
    # order based on distal length, largest length last (because it will split of first)
    p_placements <- p_placements[with(p_placements, order(p_placements$distal_length, decreasing = FALSE)),]
    
    #######################################
    if(case_is_tip){   
    #######################################
      ##### I  adding a new tip to the list
      #---------------------------------------
      p_loc_prev <- orig_tip[orig_tip$orig_edge==p_loc_orig,]$prev_edge
      for(p in 1:nrow(p_placements)){
        orig_tip[nrow(orig_tip) + 1,] = list(-p_loc_orig, p_loc_prev, p_loc_prev + p, as.character(p_placements[p,]$name))
      }
      # re-number the tips, reference always before plcement
      orig_tip <- orig_tip[with(orig_tip, order(orig_tip$prev_edge, orig_tip$new_edge, decreasing = FALSE)),]
      orig_tip$new_edge <- c(1:nrow(orig_tip))
      
      ##### II split the edge to which the placement was added
      #---------------------------------------
      # this will create a one or more new internal node(s)/edge(s)
      max_edge <- max(abs(orig_node$orig_edge))
      for(p in 1:nrow(p_placements)){
        new_edge <- max_edge + p
        orig_node[nrow(orig_node) + 1,] = list(-new_edge, p_loc_prev, -new_edge)
      }
      p_parent <- get.ancestors(p_loc_prev, tree_edges, type="parent")
      # re-number the internal edges
      for(r in 1:nrow(orig_node)){
        if(orig_node[r,]$prev_edge < previous_root){ #new edges
          e <- (abs(orig_node[r,]$orig_edge) - max_edge - 1)
          orig_node[r,]$new_edge <- p_parent + n_placements + e
        } else if(orig_node[r,]$prev_edge >= p_parent){ #edges behind the placement
          orig_node[r,]$new_edge <- orig_node[r,]$prev_edge + (2*n_placements)
        } else if(orig_node[r,]$prev_edge < p_parent){ #edges before the placement
          orig_node[r,]$new_edge <- orig_node[r,]$prev_edge  + (1*n_placements)
        }
      }
      
      ##### III update the edgeList
      #---------------------------------------
      for(row in 1:nrow(tree_edges)) {
        for(col in 1:ncol(tree_edges)) {
          tree_edges[row, col] <- update.edge(tree_edges[row, col])
        }
      }
      
      ##### IV adding the new edge(s) 
      #---------------------------------------
      #(edge previousParent->p_loc_prev was split)
      # edge previousParent->newParent (largest) is already in the tree
      newParents <- orig_node[orig_node$prev_edge==p_loc_prev,]$new_edge
      newParents <- newParents[order(newParents, decreasing = TRUE)]
      newTips <- orig_tip[orig_tip$prev_edge==p_loc_prev & orig_tip$orig_edge < 0,]$new_edge
      newTips <- newTips[order(newTips, decreasing = TRUE)]
      for(b in 1:n_placements){
        newParent <- newParents[b]
        if(b < n_placements){
          newTip1 <- newParents[b+1] 
        } else{
          newTip1 <- p_loc_prev
        }
        newTip2 <- newTips[b]
        tree_edges <- rbind(tree_edges, c(newParent, newTip1)) # parent to reference tip
        tree_edges <- rbind(tree_edges, c(newParent, newTip2)) # parent to placement tip
      }
      tree_edges <- tree_edges[order(tree_edges[,2]),]

      if(BL){
        ##### V update the branch length matrix
        #---------------------------------------
        for(b in 1:nrow(branch_lengths)) {
          branch_lengths[b,]$new_edge <- update.edge(branch_lengths[b,]$prev_edge)
        }
        ##### VI correcting the length of the edge that was split
        #---------------------------------------
        p_distal <- max(p_placements$distal_length)
        prev_length <- branch_lengths[branch_lengths$prev_edge==p_loc_prev,]$length
        new_length <- prev_length - p_distal
        branch_lengths[branch_lengths$prev_edge==p_loc_prev,]$orig_edge <- -branch_lengths[branch_lengths$prev_edge==p_loc_prev,]$orig_edge
        if(new_length >= 0){
            branch_lengths[branch_lengths$prev_edge==p_loc_prev,]$length <- new_length
          } else{
            branch_lengths[branch_lengths$prev_edge==p_loc_prev,]$length = 0
            warning("branch length < 0 set to length 0")
          }
        
        ##### VII adding the new tips
        #---------------------------------------
        p_placements2 <- p_placements[with(p_placements, order(p_placements$distal_length, decreasing = TRUE)),]
        remainder_length <- prev_length
        for(b in 1:n_placements){
          p_pendant <- p_placements2[b,]$pendant_length
          # the tip edge
          branch_lengths[nrow(branch_lengths)+1,] = list(-p_loc_orig, p_loc_prev, newTips[b], p_pendant)
          p_distal <-  p_placements2[b,]$distal_length
          # the node edge
          if(b>1){
            node_distal <- remainder_length - p_distal
            if(node_distal != 0){
              branch_lengths[nrow(branch_lengths)+1,] = list(-p_loc_orig, p_loc_prev, newParents[b], node_distal)
            } else{
              ### special case: same branch length, so polytomy with previous tip
              parent_tip1 <- get.ancestors(newTips[b], edgeList=tree_edges, type="parent")
              parent_tip2 <- get.ancestors(newTips[b-1], edgeList=tree_edges, type="parent")
              if((parent_tip2 - parent_tip1) == 1){ #should be impossible to be anything else than 1
                ### adapt the orig_node matrix
                orig_node <- orig_node[-which(orig_node$new_edge == parent_tip2),]
                orig_node$new_edge<-unlist(sapply(orig_node$new_edge, function(t){if(t >= parent_tip2){t = t-1}else{t}}))
                ### adapt the edge_tree matrix
                grandparent_tip2 <- get.ancestors(parent_tip2, edgeList=tree_edges, type="parent")
                tree_edges[which(tree_edges[,1] == parent_tip2 & tree_edges[,2] == parent_tip1),1]<-grandparent_tip2
                tree_edges<-tree_edges[-which(tree_edges[,1] == grandparent_tip2 & tree_edges[,2] == parent_tip2),]
                tree_edges[which(tree_edges[,1] == parent_tip2 & tree_edges[,2] == newTips[b-1]),1]<-parent_tip1
                for(rw in 1:nrow(tree_edges)){
                  for(cl in 1:ncol(tree_edges)){
                    if(tree_edges[rw,cl] > parent_tip1){
                      tree_edges[rw,cl] <- tree_edges[rw,cl]-1
                    }
                  }
                }
                ### adapt the branch matrix
                for(rw in 1:nrow(branch_lengths)){
                  if(branch_lengths[rw,]$new_edge > parent_tip1){
                    branch_lengths[rw,]$new_edge <- branch_lengths[rw,]$new_edge-1
                  }
                }
                
              } else{
                print(paste("w1", parent_tip1, parent_tip2, newTips[b]))
                warning("in introducing a polytomy, something went wrong")
              }
            }
          }
          remainder_length <- p_distal
        }
        
        ##### VIII adding the original reference tip back
        #---------------------------------------
        newRef <- orig_tip[orig_tip$prev_edge==p_loc_prev & orig_tip$orig_edge > 0,]$new_edge
        branch_lengths[nrow(branch_lengths)+1,] = list(p_loc_orig, p_loc_prev, newRef, remainder_length)
        branch_lengths <- branch_lengths[order(branch_lengths$new_edge),]
      } else if(!BL & n_placements > 1){
        ##### V-bis collapsing polytomies when there are no branch lenths
        #---------------------------------------
        for(b in 2:n_placements){
          parent_tip1 <- get.ancestors(newTips[b], edgeList=tree_edges, type="parent")
          parent_tip2 <- get.ancestors(newTips[b-1], edgeList=tree_edges, type="parent")
          if((parent_tip2 - parent_tip1) == 1){ #should be impossible to be anything else than 1
            ### adapt the orig_node matrix
            orig_node <- orig_node[-which(orig_node$new_edge == parent_tip2),]
            orig_node$new_edge<-unlist(sapply(orig_node$new_edge, function(t){if(t >= parent_tip2){t = t-1}else{t}}))
            ### adapt the edge_tree matrix
            grandparent_tip2 <- get.ancestors(parent_tip2, edgeList=tree_edges, type="parent")
            tree_edges[which(tree_edges[,1] == parent_tip2 & tree_edges[,2] == parent_tip1),1]<-grandparent_tip2
            tree_edges<-tree_edges[-which(tree_edges[,1] == grandparent_tip2 & tree_edges[,2] == parent_tip2),]
            tree_edges[which(tree_edges[,1] == parent_tip2 & tree_edges[,2] == newTips[b-1]),1]<-parent_tip1
            for(rw in 1:nrow(tree_edges)){
              for(cl in 1:ncol(tree_edges)){
                if(tree_edges[rw,cl] > parent_tip1){
                  tree_edges[rw,cl] <- tree_edges[rw,cl]-1
                }
              }
            }
          } else{
            print(paste("w1", parent_tip1, parent_tip2, newTips[b]))
            warning("in introducing a polytomy, something went wrong")
          }
        }
      }
    #######################################
    } else if(case_is_node){ 
    #######################################
      ##### I  adding a new tip to the list
      #---------------------------------------
      p_loc_prev <- orig_node[orig_node$orig_edge==p_loc_orig,]$prev_edge
      # first look up the most logical place within the list of tips
      # here defined as one higher than the highest descendant-tip from the branch where the tip will be placed
      dsc <- get.descendants(p_loc_prev,tree_edges, type="all")
      dsc <- max(dsc[dsc < previous_root]) #place after which to insert the new node
      dsc <- update.edge(dsc)
      for(p in 1:nrow(p_placements)){
        orig_tip[nrow(orig_tip) + 1,] = list(-p_loc_orig, dsc, 0, as.character(p_placements[p,]$name))
      }
      # re-number the tips, reference always before plcement
      orig_tip <- orig_tip[with(orig_tip, order(orig_tip$prev_edge, -orig_tip$new_edge, decreasing = FALSE)),]
      orig_tip$new_edge <- c(1:nrow(orig_tip))

      ##### II split the edge to which the placement was added
      #---------------------------------------
      # this will create a one or more new internal node(s)/edge(s)
      max_edge <- max(abs(orig_node$orig_edge))
      for(p in 1:nrow(p_placements)){
        new_edge <- max_edge + p
        orig_node[nrow(orig_node) + 1,] = list(-new_edge, p_loc_prev, -new_edge)
      }
      # re-number the internal edges
      for(r in 1:nrow(orig_node)){
        if(orig_node[r,]$prev_edge == p_loc_prev){ 
          if(orig_node[r,]$orig_edge > 0){ #target edge
            orig_node[r,]$new_edge <- p_loc_prev + n_placements
          } else{ #new edges
            e <- (abs(orig_node[r,]$orig_edge) - max_edge)
            orig_node[r,]$new_edge <- p_loc_prev + n_placements + e
          }
        } else if(orig_node[r,]$prev_edge > p_loc_prev){
          orig_node[r,]$new_edge <- orig_node[r,]$prev_edge + 2*(n_placements)
        } else if(orig_node[r,]$prev_edge < p_loc_prev){
          orig_node[r,]$new_edge <- orig_node[r,]$prev_edge + 1*(n_placements)
        }
      }

      ##### III update the edgeList
      #---------------------------------------
      for(row in 1:nrow(tree_edges)) {
        for(col in 1:ncol(tree_edges)) {
          if(tree_edges[row, col] == p_loc_prev){
            if(col==1){ # special case: branch to split: new distal part (away from root)
              tree_edges[row, col] <- p_loc_prev + (1*n_placements)
            } else{ # special case: branch to split: new proximal part (close to root)
              tree_edges[row, col] <- p_loc_prev + (2*n_placements)
            }
          } else{
            tree_edges[row, col] <- update.edge(tree_edges[row, col])
          }
        }
      }
      
      ##### IV adding the new edge(s) 
      #---------------------------------------
      #(edge previousParent->p_loc_prev was split)
      # edge previousParent->newParent (largest) is already in the tree
      newParents <- orig_node[orig_node$prev_edge==p_loc_prev & orig_node$orig_edge < 0,]$new_edge
      newParents <- newParents[order(newParents, decreasing = TRUE)]
      newTips <- orig_tip[orig_tip$orig_edge==-p_loc_orig,]$new_edge
      newTips <- newTips[order(newTips, decreasing = TRUE)]
      for(b in 1:n_placements){
        newParent <- newParents[b]
        if(b < n_placements){
          newTip1 <- newParents[b+1] 
        } else{
          newTip1 <- p_loc_prev + n_placements
        }
        newTip2 <- newTips[b]
        tree_edges <- rbind(tree_edges, c(newParent, newTip1)) # parent to reference tip
        tree_edges <- rbind(tree_edges, c(newParent, newTip2)) # parent to placement tip
      }
      tree_edges <- tree_edges[order(tree_edges[,2]),]
      
      if(BL){
        ##### V update the branch length matrix
        #---------------------------------------
        for(b in 1:nrow(branch_lengths)) {
          branch_lengths[b,]$new_edge <- update.edge(branch_lengths[b,]$prev_edge)
        }
        ##### VI correcting the length of the edge that was split
        #---------------------------------------
        p_distal <- max(p_placements$distal_length)
        prev_length <- branch_lengths[branch_lengths$prev_edge==p_loc_prev,]$length
        new_length <- prev_length - p_distal
        if(new_length >= 0){
          branch_lengths[branch_lengths$prev_edge==p_loc_prev,]$length <- new_length
        } else{
          branch_lengths[branch_lengths$prev_edge==p_loc_prev,]$length = 0
          warning("branch length < 0 set to length 0")
        }
        
        ##### VII adding the new tips
        #---------------------------------------
        p_placements2 <- p_placements[with(p_placements, order(p_placements$distal_length, decreasing = TRUE)),]
        remainder_length <- prev_length
        for(b in 1:n_placements){
          p_pendant <- p_placements2[b,]$pendant_length
          # the tip edge
          branch_lengths[nrow(branch_lengths)+1,] = list(-p_loc_orig, p_loc_prev, newTips[b], p_pendant)
          p_distal <-  p_placements2[b,]$distal_length
          # the node edge
          if(b>1){
            node_distal <- remainder_length - p_distal
            if(node_distal != 0){
              branch_lengths[nrow(branch_lengths)+1,] = list(-p_loc_orig, p_loc_prev, newParents[b], node_distal)
            } else{
              ### special case: same branch length, so polytomy with previous tip
              parent_tip1 <- get.ancestors(newTips[b], edgeList=tree_edges, type="parent")
              parent_tip2 <- get.ancestors(newTips[b-1], edgeList=tree_edges, type="parent")
              if((parent_tip2 - parent_tip1) == 1){ #should be impossible to be anything else than 1
                ### adapt the orig_node matrix
                orig_node <- orig_node[-which(orig_node$new_edge == parent_tip2),]
                orig_node$new_edge<-unlist(sapply(orig_node$new_edge, function(t){if(t >= parent_tip2){t = t-1}else{t}}))
                ### adapt the edge_tree matrix
                grandparent_tip2 <- get.ancestors(parent_tip2, edgeList=tree_edges, type="parent")
                tree_edges[which(tree_edges[,1] == parent_tip2 & tree_edges[,2] == parent_tip1),1]<-grandparent_tip2
                tree_edges<-tree_edges[-which(tree_edges[,1] == grandparent_tip2 & tree_edges[,2] == parent_tip2),]
                tree_edges[which(tree_edges[,1] == parent_tip2 & tree_edges[,2] == newTips[b-1]),1]<-parent_tip1
                for(rw in 1:nrow(tree_edges)){
                  for(cl in 1:ncol(tree_edges)){
                    if(tree_edges[rw,cl] > parent_tip1){
                      tree_edges[rw,cl] <- tree_edges[rw,cl]-1
                    }
                  }
                }
                ### adapt the branch matrix
                for(rw in 1:nrow(branch_lengths)){
                  if(branch_lengths[rw,]$new_edge > parent_tip1){
                    branch_lengths[rw,]$new_edge <- branch_lengths[rw,]$new_edge-1
                  }
                }
                
              } else{
                print(paste("w1", parent_tip1, parent_tip2, newTips[b]))
                warning("in introducing a polytomy, something went wrong")
              }
            }
            
          }
          remainder_length <- p_distal
        }
        ##### VIII adding the remaining edge leading back
        #---------------------------------------
        newRef <- orig_node[!orig_node$new_edge %in% branch_lengths$new_edge,]$new_edge
        branch_lengths[nrow(branch_lengths)+1,] = list(-p_loc_orig, p_loc_prev, newRef, remainder_length)
        branch_lengths <- branch_lengths[order(branch_lengths$new_edge),]
      } else if(!BL & n_placements > 1){
        ##### V-bis collapsing polytomies when there are no branch lenths
        #---------------------------------------
        for(b in 2:n_placements){
          parent_tip1 <- get.ancestors(newTips[b], edgeList=tree_edges, type="parent")
          parent_tip2 <- get.ancestors(newTips[b-1], edgeList=tree_edges, type="parent")
          if((parent_tip2 - parent_tip1) == 1){ #should be impossible to be anything else than 1
            ### adapt the orig_node matrix
            orig_node <- orig_node[-which(orig_node$new_edge == parent_tip2),]
            orig_node$new_edge<-unlist(sapply(orig_node$new_edge, function(t){if(t >= parent_tip2){t = t-1}else{t}}))
            ### adapt the edge_tree matrix
            grandparent_tip2 <- get.ancestors(parent_tip2, edgeList=tree_edges, type="parent")
            tree_edges[which(tree_edges[,1] == parent_tip2 & tree_edges[,2] == parent_tip1),1]<-grandparent_tip2
            tree_edges<-tree_edges[-which(tree_edges[,1] == grandparent_tip2 & tree_edges[,2] == parent_tip2),]
            tree_edges[which(tree_edges[,1] == parent_tip2 & tree_edges[,2] == newTips[b-1]),1]<-parent_tip1
            for(rw in 1:nrow(tree_edges)){
              for(cl in 1:ncol(tree_edges)){
                if(tree_edges[rw,cl] > parent_tip1){
                  tree_edges[rw,cl] <- tree_edges[rw,cl]-1
                }
              }
            }
          } else{
            print(paste("w1", parent_tip1, parent_tip2, newTips[b]))
            warning("in introducing a polytomy, something went wrong")
          }
        }
      }
    #######################################
    }  else if(case_is_root){
    #######################################
      ##### I  adding a new tip to the list
      #---------------------------------------
      p_loc_prev <- previous_root
      for(p in 1:nrow(p_placements)){
        orig_tip[nrow(orig_tip) + 1,] = list(-p_loc_orig, p_loc_prev, p_loc_prev + (p-1), as.character(p_placements[p,]$name))
      }
      
      ##### II split the edge to which the placement was added
      #---------------------------------------
      # this will create a one or more new internal node(s)/edge(s)
      max_edge <- max(abs(orig_node$orig_edge))
      max_edge_prev <- max(abs(orig_node$prev_edge))
      for(p in 1:nrow(p_placements)){
        orig_node[nrow(orig_node) + 1,] = list(-(max_edge + p), p_loc_prev, max_edge_prev + p)
      }
      # re-number the internal edges
      for(r in 1:nrow(orig_node)){
        orig_node[r,]$new_edge <- orig_node[r,]$new_edge + 1*(n_placements)
      }
      
      ##### III update the edgeList
      #---------------------------------------
      for(row in 1:nrow(tree_edges)) {
        for(col in 1:ncol(tree_edges)) {
          if(tree_edges[row, col] == p_loc_prev){
            tree_edges[row, col] <- min(orig_node[orig_node$prev_edge == p_loc_prev,]$new_edge)
          } else{
            tree_edges[row, col] <- update.edge(tree_edges[row, col])
          }
        }
      }
      
      ##### IV adding the new edge(s) 
      #---------------------------------------
      #(edge previousParent->p_loc_prev was split)
      # edge previousParent->newParent (largest) is already in the tree
      newParents <- orig_node[orig_node$prev_edge==p_loc_prev & orig_node$orig_edge < 0,]$new_edge
      newParents <- newParents[order(newParents, decreasing = TRUE)]
      newTips <- orig_tip[orig_tip$orig_edge==-orig_root,]$new_edge
      newTips <- newTips[order(newTips, decreasing = FALSE)]
      for(b in 1:n_placements){
        newTip1 <- newTips[b]
        if(b < n_placements){
          newParent <- newParents[b]
          newTip2 <- newParents[b+1]
        } else{
          newParent <- current_root
          newTip2 <- max(newParents)
        }
        tree_edges <- rbind(tree_edges, c(newParent, newTip1)) # parent to reference tip
        tree_edges <- rbind(tree_edges, c(newParent, newTip2)) # parent to placement tip
      }
      tree_edges <- tree_edges[order(tree_edges[,2]),]
      
      if(BL){
        ##### V update the branch length matrix
        #---------------------------------------
        for(b in 1:nrow(branch_lengths)) {
          branch_lengths[b,]$new_edge <- update.edge(branch_lengths[b,]$prev_edge)
        }
        
        ##### VI adding the new tips
        #---------------------------------------
        p_placements2 <- p_placements[with(p_placements, order(p_placements$distal_length, decreasing = TRUE)),]
        remainder_length <- max(p_placements2$distal_length)
        for(b in 1:n_placements){
          p_pendant <- p_placements2[b,]$pendant_length
          # the tip edge
          branch_lengths[nrow(branch_lengths)+1,] = list(-p_loc_orig, p_loc_prev, newTips[b], p_pendant)
          p_distal <-  p_placements2[b,]$distal_length
          # the node edge
          if(b>1){
            node_distal <- remainder_length - p_distal
            if(node_distal != 0){
              branch_lengths[nrow(branch_lengths)+1,] = list(-p_loc_orig, p_loc_prev, newParents[b], node_distal)
            } else{
              ### special case: same branch length, so polytomy with previous tip
              parent_tip1 <- get.ancestors(newTips[b], edgeList=tree_edges, type="parent")
              parent_tip2 <- get.ancestors(newTips[b-1], edgeList=tree_edges, type="parent")
              ### adapt the orig_node matrix
              orig_node <- orig_node[-which(orig_node$new_edge == parent_tip2),]
              orig_node$new_edge<-unlist(sapply(orig_node$new_edge, function(t){if(t >= parent_tip2){t = t-1}else{t}}))
              ### adapt the edge_tree matrix
              if(parent_tip1==current_root){
                tree_edges<-tree_edges[-which(tree_edges[,1] == current_root & tree_edges[,2] == parent_tip2),]
                tree_edges[which(tree_edges[,1] == parent_tip2),1]<- current_root
              } else{
                grandparent_tip2 <- get.ancestors(parent_tip2, edgeList=tree_edges, type="parent")
                tree_edges[which(tree_edges[,1] == parent_tip2 & tree_edges[,2] == parent_tip1),1]<-grandparent_tip2
                tree_edges<-tree_edges[-which(tree_edges[,1] == grandparent_tip2 & tree_edges[,2] == parent_tip2),]
                tree_edges[which(tree_edges[,1] == parent_tip2 & tree_edges[,2] == newTips[b-1]),1]<-parent_tip1
                for(rw in 1:nrow(tree_edges)){
                  for(cl in 1:ncol(tree_edges)){
                    if(tree_edges[rw,cl] > parent_tip1){
                      tree_edges[rw,cl] <- tree_edges[rw,cl]-1
                    }
                  }
                }
                ### adapt the branch matrix
                for(rw in 1:nrow(branch_lengths)){
                  if(branch_lengths[rw,]$new_edge > parent_tip1){
                    branch_lengths[rw,]$new_edge <- branch_lengths[rw,]$new_edge-1
                  }
                }
              }
            }
            
          }
          remainder_length <- p_distal
        }
        ##### VIII adding the remaining edge leading back
        #---------------------------------------
        newRef <- orig_node[!orig_node$new_edge %in% branch_lengths$new_edge,]$new_edge
        branch_lengths[nrow(branch_lengths)+1,] = list(-p_loc_orig, p_loc_prev, newRef, remainder_length)
        branch_lengths <- branch_lengths[order(branch_lengths$new_edge),]
      } else if(!BL & n_placements > 1){
        ##### V-bis collapsing polytomies when there are no branch lenths
        #---------------------------------------
        for(b in 2:n_placements){
          parent_tip1 <- get.ancestors(newTips[b], edgeList=tree_edges, type="parent")
          parent_tip2 <- get.ancestors(newTips[b-1], edgeList=tree_edges, type="parent")
          ### adapt the orig_node matrix
          orig_node <- orig_node[-which(orig_node$new_edge == parent_tip2),]
          orig_node$new_edge<-unlist(sapply(orig_node$new_edge, function(t){if(t >= parent_tip2){t = t-1}else{t}}))
          ### adapt the edge_tree matrix
          grandparent_tip2 <- get.ancestors(parent_tip2, edgeList=tree_edges, type="parent")
          tree_edges[which(tree_edges[,1] == parent_tip2 & tree_edges[,2] == parent_tip1),1]<-grandparent_tip2
          tree_edges<-tree_edges[-which(tree_edges[,1] == grandparent_tip2 & tree_edges[,2] == parent_tip2),]
          tree_edges[which(tree_edges[,1] == parent_tip2 & tree_edges[,2] == newTips[b-1]),1]<-parent_tip1
          for(rw in 1:nrow(tree_edges)){
            for(cl in 1:ncol(tree_edges)){
              if(tree_edges[rw,cl] > parent_tip1){
                tree_edges[rw,cl] <- tree_edges[rw,cl]-1
              }
            }
          }
        }
      }
    } 
  }
  
  #---------------------------------------------
  # 3) putting everyting back into a phylo object
  tree_edges <- cbind(sapply(tree_edges[,1],function(x){x<-as.integer(x)}),
                sapply(tree_edges[,2],function(x){x<-as.integer(x)}))
  tree_edges <- tree_edges[order(tree_edges[,2]),]
  tip.label <- orig_tip$label
  Nnode <- as.integer(length(table(tree_edges[,1])))
  if(BL){
    edge.length <- as.vector(branch_lengths$length)
    Ref.Phylo.updated <- list(edge = tree_edges, edge.length = edge.length, Nnode = Nnode, tip.label = tip.label)
  } else{
    Ref.Phylo.updated <- list(edge = tree_edges, Nnode = Nnode, tip.label = tip.label)
  }
  class(Ref.Phylo.updated) <- "phylo"
  attr(Ref.Phylo.updated, "order") <- "cladewise"

  return(Ref.Phylo.updated)
}


test<-get.placements.phylo(test_jplace, verbose=TRUE)

checkValidPhylo(test)
table(test$edge[,1])

edge=data.frame(test$edge, edge_num=test$edge[,2])
colnames(edge)=c("parent", "node", "edge_num")
t <- ggtree(test, ladderize=FALSE) + geom_tiplab() + theme_tree2() + xlim(0, 100) + geom_text2(aes(subset = !isTip, label=label)) 
t %<+% edge + geom_label(aes(x=branch, label=edge_num))





    
