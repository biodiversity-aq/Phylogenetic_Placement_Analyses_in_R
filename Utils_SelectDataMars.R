# Functions to work with large mars datasets

require(tibble)
require(data.table)
require(ggplot2)
require(sf)
require(rnaturalearth)
require(rnaturalearthdata)


create.dataset <- function(SeqSetPath = NA){
  if(is.na(SeqSetPath)){
    stop()
  }
  
  DataOutput = data.frame('datasetID'=NA,'datasetName'=NA,
                          'parentEventID'=NA,'eventID'=NA,'eventName'=NA,
                          'decimalLatitude'=NA,'decimalLongitude'=NA,
                          'eventDate'=NA,'locality'=NA,
                          'targetGroup'=NA,'marker'=NA,'markerSubfragment'=NA,
                          'seqMethod'=NA,'investigation_type'=NA,'submitted_to_insdc'=NA, 'seqdata_URL'=NA,
                          'studyNumber'=NA, 'geneticAccessionNumber'=NA,'BioprojectNumber'=NA,
                          'primerSequenceForward'=NA,'primerSequenceReverse'=NA,'primerNameForward'=NA,'primerNameReverse'=NA,'mid'=NA,
                          'lib_reads_seqd'=NA,
                          'env_biome'=NA, 'env_feature'=NA, 'env_material'=NA,
                          'depth'=NA,'elev'=NA,'conduc'=NA,'salinity'=NA,'temp'=NA)
  
  
  
  totalCount = 0
  datasetCount = 0
  totalEventCount = 0
  
  fileNameList <- dir(SeqSetPath, pattern = "SeqSet_")
  for(fileName in fileNameList){
    print(fileName)
    sampleCount = 0
    eventList = c(); eventCount = 0
    additionalDataList <- c('depth', 'elev', 'conduc', 'salinity', 'temp')
    
    datasetCount = datasetCount + 1
    
    SeqSet = fileName
    MIxS = gsub("SeqSet_", "MiMARKS_", fileName)
    SeqSetData = read.csv(paste(SeqSetPath, SeqSet, sep='/'), header = FALSE, sep=',')
    names(SeqSetData) <- names(read.csv(paste(SeqSetPath, SeqSet, sep='/'), header = TRUE, sep=','))
    
    MIxSData = tryCatch({read.csv(paste(SeqSetPath, MIxS, sep='/'), header = TRUE, sep=',')},
                        error = function (e){read.csv(paste(SeqSetPath, gsub("SeqSet_", "MIxS_", fileName), sep='/')
                                                      , header = TRUE, sep=',')})
    
    #remove empty trailing columns
    MIxSData <- MIxSData[,colSums(is.na(MIxSData))<nrow(MIxSData)] 
    SeqSetData <- SeqSetData[,colSums(is.na(SeqSetData))<nrow(SeqSetData)] 
    
    
    #add a 2nd and 3th column to SeqSet, to be able to combine them with the MiMARKS
    SeqSetData <- add_column(SeqSetData, section = "SeqSet", .after = 0)
    SeqSetData <- add_column(SeqSetData, units = NA, .after = 2)
    names(SeqSetData)[2]<-"Structured.Comment.Name"
    
    #combine SeqSet and MIxS data
    combData <- tryCatch({rbind(SeqSetData, MIxSData)}, 
                         error = function(e){print('error in rbind');
                           print(paste('    colname not in MIxSData:', setdiff(names(SeqSetData), names(MIxSData))));
                           print(paste('    colname not in SeqSetData:', setdiff(names(MIxSData), names(SeqSetData))))})
    
    for(ns in 4:ncol(combData)){
      totalCount = totalCount + 1
      sampleCount = sampleCount + 1
      
      DataOutput[totalCount,'datasetID'] <- sprintf("D%04d", datasetCount)
      DataOutput[totalCount,'datasetName'] <- paste(sprintf("D%04d", datasetCount), gsub("SeqSet_", "", fileName),sep='_')
      
      if("DwC_event" %in% combData$Structured.Comment.Name){
        Par_eventName <- as.character(combData[combData$Structured.Comment.Name == "DwC_event",ns])
        if (Par_eventName %in% eventList){
          Par_eventID <- sprintf("PE%04d", match(Par_eventName,eventList))
        } else{
          eventCount <- eventCount+1
          Par_eventID = sprintf("PE%04d", eventCount)
          eventList[eventCount] <- Par_eventName
        }
      } else{
        eventCount <- eventCount+1
        Par_eventName <- as.character(combData[combData$Structured.Comment.Name == "unique_sequence_set_id",ns])
        Par_eventID <- paste(sprintf("PE%04d", eventCount))
      }
      
      DataOutput[totalCount,'parentEventID'] <- paste(Par_eventID, Par_eventName, sep='_')
      
      DataOutput[totalCount,'eventID'] <- paste(Par_eventID, sprintf("S%04d", sampleCount), sep='_')
      DataOutput[totalCount,'eventName'] <- as.character(combData[combData$Structured.Comment.Name == "unique_sequence_set_id",ns])
      
      lat_lon<-tryCatch({as.character(combData[combData$Structured.Comment.Name == "lat_lon",ns])},
                        error=function(e){"NA NA"})
      DataOutput[totalCount,'decimalLatitude'] <- as.numeric(strsplit(lat_lon, ' ')[[1]][1])
      DataOutput[totalCount,'decimalLongitude'] <- as.numeric(strsplit(lat_lon, ' ')[[1]][2])
      DataOutput[totalCount,'eventDate'] <- as.character(combData[combData$Structured.Comment.Name == "collection_date",ns][[1]])
      
      if("DwC_locality" %in% combData$Structured.Comment.Name){
        DataOutput[totalCount,'locality'] <- as.character(combData[combData$Structured.Comment.Name == "DwC_locality",ns])
      } else{
        DataOutput[totalCount,'locality'] <- as.character(combData[combData$Structured.Comment.Name == "geo_loc_name",ns])
      }
      
      DataOutput[totalCount,'targetGroup'] <- as.character(combData[combData$Structured.Comment.Name == "target_taxa",ns])
      DataOutput[totalCount,'marker'] <- as.character(combData[combData$Structured.Comment.Name == "target_gene",ns][[1]])
      DataOutput[totalCount,'markerSubfragment'] <- as.character(combData[combData$Structured.Comment.Name == "region_targeted",ns])
      DataOutput[totalCount,'seqMethod'] <- as.character(combData[combData$Structured.Comment.Name == "sequencing_technology",ns])
      DataOutput[totalCount,'investigation_type'] <- as.character(combData[combData$Structured.Comment.Name == "investigation_type",ns])
      DataOutput[totalCount,'submitted_to_insdc'] <- as.character(combData[combData$Structured.Comment.Name == "submitted_to_insdc",ns])
      DataOutput[totalCount,'seqdata_URL'] <- as.character(combData[combData$Structured.Comment.Name == "url_data_repository",ns])
      
      if(as.character(combData[combData$Structured.Comment.Name == "url_data_repository",ns]) != ""){
        SRR_nr <- as.character(combData[combData$Structured.Comment.Name == "url_data_repository",ns])
        SRR_nr <- strsplit(SRR_nr , '/')[[1]][length(strsplit(SRR_nr , '/')[[1]])]
        DataOutput[totalCount,'geneticAccessionNumber'] <- SRR_nr
      }
      
      if(as.character(combData[combData$Structured.Comment.Name == "sra_accession_number",ns]) != ""){
        DataOutput[totalCount,'studyNumber'] <- as.character(combData[combData$Structured.Comment.Name == "sra_accession_number",ns])
      } else if(as.character(combData[combData$Structured.Comment.Name == "genbank_accession_numbers",ns]) != ""){
        DataOutput[totalCount,'studyNumber'] <- as.character(combData[combData$Structured.Comment.Name == "genbank_accession_numbers",ns])
      }
      DataOutput[totalCount,'BioprojectNumber'] <- as.character(combData[combData$Structured.Comment.Name == "sra_project_number",ns])
      DataOutput[totalCount,'primerSequenceForward'] <- as.character(combData[combData$Structured.Comment.Name == "forward_primer_sequence",ns])
      DataOutput[totalCount,'primerSequenceReverse'] <- as.character(combData[combData$Structured.Comment.Name == "reverse_primer_sequence",ns])
      DataOutput[totalCount,'primerNameForward'] <- as.character(combData[combData$Structured.Comment.Name == "forward_primer",ns])
      DataOutput[totalCount,'primerNameReverse'] <- as.character(combData[combData$Structured.Comment.Name == "reverse_primer",ns])
      DataOutput[totalCount,'lib_reads_seqd'] <- as.character(combData[combData$Structured.Comment.Name == "number_of_spots_in_sra",ns])
      
      DataOutput[totalCount,'env_biome'] <- as.character(combData[combData$Structured.Comment.Name == "env_biome",ns])
      DataOutput[totalCount,'env_feature'] <- as.character(combData[combData$Structured.Comment.Name == "env_feature",ns])
      DataOutput[totalCount,'env_material'] <- as.character(combData[combData$Structured.Comment.Name == "env_material",ns])
      
      
      for(additionalDataName in additionalDataList){
        if(additionalDataName %in% combData$Structured.Comment.Name){
          DataOutput[totalCount,additionalDataName] <- as.character(combData[combData$Structured.Comment.Name == additionalDataName,ns])
        }
      }
      
    }
    totalEventCount = totalEventCount + eventCount
  }
  
  #print stats of process
  writeLines(paste('\n\n', datasetCount, 'datasets processed', '\n', 
                   totalEventCount, 'physical samples', '\n',
                   totalCount, 'sequence sets','\n'))
  
  #plot points on map
  world <- ne_countries(scale = "medium", returnclass = "sf")
  DataOnWorld<- ggplot(data = world) +
    geom_sf(color = "antiquewhite2", fill = "antiquewhite") + 
    xlab("Longitude") + ylab("Latitude") +
    theme(panel.background = element_rect(fill = "aliceblue")) +
    ggtitle("sampled points (with data) represented on mARS.biodiversity.aq") +
    geom_point(data = DataOutput, aes(x=decimalLongitude, y=decimalLatitude), color = "red", size = 0.8) 
  
  plot(DataOnWorld)
  
  return(DataOutput)
}

combine.data.frame <- function(df1, df2, fill=NA){
  if(ncol(df1) > 0 & nrow(df1) > 0){
    df_UpRight <- data.frame(matrix(nrow = nrow(df1), ncol = ncol(df2), data=fill))
    colnames(df_UpRight) <- colnames(df2)
    rownames(df_UpRight) <- rownames(df1)
    df1_b <- cbind(df1, df_UpRight)
    
    df_DownLeft <- data.frame(matrix(nrow = nrow(df2), ncol = ncol(df1), data=fill))
    rownames(df_DownLeft) <- rownames(df2)
    colnames(df_DownLeft) <- colnames(df1)
    df2_b <- cbind(df_DownLeft, df2)
    
    df_out <- data.frame(rbind(df1_b, df2_b))
  } else(
    df_out <- data.frame(df2)
  )

  return(df_out)
}

scan.sequences <- function(fasta = NULL, SeqNames = NULL){
  require(seqinr)
  ### function that runs through a fasta file or DNAStringset, and gets out the desired sequences
  # @param:fasta = a fasta file or DNAStringSet
  # @param:SeqNames = a vector with the names to be extracted
  # @output = a DNAStringSet of the selected sequences
  if(is.null(fasta) | is.null(SeqNames)){
    stop("Please provide a valid input for both fasta and SeqNames arguments")
  }
  if(!class(fasta)[1] == "DNAStringSet"){
    if(class(fasta)[1] == "characer" && file.exists(fasta)){
      fasta <-  readDNAStringSet(fasta)
    } else{
      stop("fasta argument must be either DNAStringSet or a valid file path")
    }
  }
  fastaOut <- fasta[SeqNames]
  return(fastaOut)
}

get.taxon <- function(fasta = NULL, taxtab = NULL, taxon = NULL, taxrank=c("phylum")){
  require(seqinr)
  ### function extracts all sequences from a given taxon from a fasta file
  # @param:fasta = a fasta file or DNAStringSet
  # @param:taxtab = a table with sequence names and their taxonomic identity
  # @param:taxon = the taxon name to be extracted
  # @param:rank = the rank of the taxon name to be extracted (either: kingdom, phylum, class, order, family, genus, species)
  # @output = a DNAStringSet of the sequences of the selected taxon
  if(is.null(fasta) | is.null(taxtab) | is.null(taxon)){
    stop("Please provide a valid input for fasta, taxtab and taxon arguments")
  }
  tax_sub <- taxtab[apply(taxtab[,taxrank], 1, function(x) any(x==taxon)), ]
  #taxnames <- as.character(tax_sub$X)
  taxnames <- rownames(tax_sub)
  #taxnames <- rownames(taxtab[eval(parse(text=paste("taxtab$", taxrank, sep=""))) == taxon,])
  if(length(taxnames)==0){
    warning("no ASVs selected, NULL returned")
    fastaOut <- NULL
  } else{
    fastaOut <- scan.sequences(fasta = fasta, SeqNames = taxnames)
  }
  return(fastaOut)
}
