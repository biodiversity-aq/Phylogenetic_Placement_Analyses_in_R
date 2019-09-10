# This script contains several functions to deal with alignments
# author Maxime Sweetlove 2019
# licence CC 4.0

#--------------------------------------------------------------
# functions to align short primer sequences to references, and get the position
#--------------------------------------------------------------
primerPosition<-function(primer, reference, start=TRUE){
  ### function that takes a primer sequence and a reference sequence
  # @param primer: a (short) primer DNA sequence, IUPAC notation
  # @param reference: a (long) reference DNA sequence, IUPAC notation
  # @param start: [boolean], TRUE: return the starting position of the alignment; FALSE end position
  ### returns the starting (start=TRUE) or ending (start=FALSE) position of the alignment
  mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -1, baseOnly = FALSE)
  
  locAlign <- tryCatch({pairwiseAlignment(pattern=toupper(primer), subject=reference, type= "local", 
                                          substitutionMatrix = mat, gapOpening = 8, 
                                          gapExtension = 4)}, 
                       error = function(e){print('error')})
  
  
  if(start==TRUE){
    return(locAlign@subject@range@start)
  } else{
    return(locAlign@subject@range@start+nchar(primer))
  }
}

split_DNAstringList_unknownSep<-function(DNAstringList_unknownSep){
  ### function that takes at string composed of DNA sequences with unkown separator
  # @param DNAstringList_unknownSep: a string with one or more DNA sequences (IUPAC notation), separated by an unknown separator
  ### returns a vector with the individual DNA sequences
  
  #trim trailing white spaces
  DNAstringList_unknownSep<-StrTrim(DNAstringList_unknownSep, pattern = " \t\n", method = "both")
  
  #assume space as default sep
  if(grepl(",", DNAstringList_unknownSep, fixed=TRUE)){
    DNAstringList_unknownSep<-gsub(" ", "", DNAstringList_unknownSep, fixed = TRUE)
    DNAstringList_unknownSep<-strsplit(DNAstringList_unknownSep, ',')[[1]]
  } else if(grepl(";", DNAstringList_unknownSep, fixed=TRUE)){
    DNAstringList_unknownSep<-gsub(" ", "", DNAstringList_unknownSep, fixed = TRUE)
    DNAstringList_unknownSep<-strsplit(DNAstringList_unknownSep, ';')[[1]]
  } else if(grepl("|", DNAstringList_unknownSep, fixed=TRUE)){
    DNAstringList_unknownSep<-gsub(" ", "", DNAstringList_unknownSep, fixed = TRUE)
    DNAstringList_unknownSep<-strsplit(DNAstringList_unknownSep, '|')[[1]]
  } else{
    DNAstringList_unknownSep<-strsplit(DNAstringList_unknownSep, ' ')[[1]]
  }
  
  return(DNAstringList_unknownSep)
}


#--------------------------------------------------------------
# functions xxx
#--------------------------------------------------------------