require(dplyr)
require(tidyr)
require(reshape)
require(reshape2)



#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

uniqueMappingFile <- args[1]
readCutoff <- args[2]
outputPath <- args[3]


# uniqueMappingFile <- "/home/schakraborty/Documents/HomeOffice/LadderSeq/FinalPipelines/FinalKallisto/uniquelyMappingReads.txt"
# percentageCutoff <- 0
# readCutoff <- 0
# # transcriptMartFile <- "/gcm-lfs1/shounak/LadderSeq/Quantification/configData/transLengthFile.txt"
# outputPath <- "/home/schakraborty/Documents/HomeOffice/LadderSeq/FinalPipelines/FinalKallisto/"


mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "hsapiens_gene_ensembl",
                         host = 'www.ensembl.org')


betaEstimator <- function(uniqueMappingFile, readCutoff, outputPath){

  ### Counting the number of reads that are mapped to a transcript
  bamdata <- data.table::fread(uniqueMappingFile, header = FALSE, strip.white = TRUE)
  bamdata <- as.data.frame(bamdata)
  names(bamdata) <- c("ReadName","TranscriptName","NN")


  #### Extracting the read band names and making a new column out of them ####
  bamdata["ReadBand"] <- as.numeric(substr(bamdata$ReadName,5, 5))
  names(bamdata) <- c("ReadName","TranscriptName","NN","ReadBand")
  bamdata <- dplyr::select(bamdata,1,2,4)

  ### Grouping the bam file data based on the transcript name and the read band
  groupedBamData <- bamdata %>%
    group_by(TranscriptName,ReadBand) %>%
    summarise(ReadCount = length(TranscriptName))







  ### Summing the read counts for each transcript ###
  normReadCount <- groupedBamData %>%
    group_by(TranscriptName) %>%
    summarise(ReadCountSum = sum(ReadCount))

  groupedBamData <- merge(groupedBamData,normReadCount,by="TranscriptName")

  groupedBamData <- groupedBamData %>% subset(grepl("ENST",TranscriptName))

  groupedBamData["PercentageReadCount"] <- groupedBamData$ReadCount/groupedBamData$ReadCountSum




  ### read transcript file
  # tLength <- read.table(transcriptMartFile, header = T)
  
  transcriptNames <- lapply(groupedBamData$TranscriptName, function(x) strsplit(as.character(x),'\\.')[[1]][1])
  transcriptNames <- as.data.frame(unlist(transcriptNames))
  names(transcriptNames) <- c("transcriptNames")
  groupedBamData$transcriptNames <- transcriptNames$transcriptNames
  tLength <- biomaRt::getBM(attributes = c("transcript_length","ensembl_transcript_id","ensembl_gene_id","transcript_biotype"),
                            filters='ensembl_transcript_id', values = groupedBamData$transcriptNames ,
                            mart = mart)
  names(tLength) <- c("TranscriptLength","transcriptNames","GeneName","transcriptBiotype")
  
  ## adding 200 to transcript length for poly a tails
  tLength["TranscriptLength"] <- tLength["TranscriptLength"] + 200
  
  # mergedTranscripts <- merge(groupedBamData,tLength, by= "transcriptNames")
  # mergedTranscripts <- mergedTranscripts[,c(2,7,8,9)]
  # names(mergedTranscripts) <- c("TranscriptName","TranscriptLength","GeneName","transcriptBiotype")

  groupedBamData <- merge(groupedBamData,tLength,by = "transcriptNames")

  ### Filtering out transcripts which are not protein coding
  groupedBamData <- groupedBamData %>%
    subset(grepl("protein_coding",groupedBamData$transcriptBiotype))


  ### Filtering out transcripts of very high length
  groupedBamData <- groupedBamData %>% subset(groupedBamData$TranscriptLength<8000)

  ### Assigning actual band where the transcript should have been ###
  groupedBamData["ActualBand"] <- groupedBamData$TranscriptLength %>% cut(c(0,1000,1500,2000,3000,4000,6000,10000),labels=FALSE)

  ### assigning correctness ###
  groupedBamData["Correctness"] <- (groupedBamData$ReadBand == groupedBamData$ActualBand)

  # ## filtering transcripts with low read counts
  filterdGroupedBamData <- groupedBamData %>%  subset(ReadCountSum>readCutoff)




#   ### Getting the low resolution band probabilities from the read counts in each band ###
#   metaAnalysis <- filterdGroupedBamData %>%
#     group_by(ActualBand,ReadBand) %>%
#     summarise(ReadCount = sum(ReadCount),NumTranscript = n())


#   bandProb <- metaAnalysis %>%
#     group_by(ActualBand) %>%
#     mutate(sum(ReadCount))

#   bandProb[3] <- round(bandProb[3]/bandProb[5],digits = 100)

#   bandProb <- bandProb[,c(1:3)]
#   betaOutputFile <- paste(outputPath,"estimatedBeta_low_text.txt",sep = "")
#   bandProbabilities <- dcast(data = bandProb,formula = ActualBand~ReadBand,value.var = "ReadCount", fill = 0)
#   bandProbabilities <- bandProbabilities[-1]
#   bandProbabilities <- format(bandProbabilities, scientific = FALSE)
#   write.table(bandProbabilities,file=betaOutputFile, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
#   ### finished with low resolution band probabilites ###






  ### Getting the high resolution band probabilities from the read counts in each band ###
  highResolutionBeta <- filterdGroupedBamData
  highResolutionBeta["ActualBand"] <-  highResolutionBeta$TranscriptLength %>%
    cut(c(0,910,1083,1268,1440,1590,1757,1920,2070,2273,2462,2701,2932,3230,3554,3893,4525,5324,6378,10000),labels=FALSE)

  metaAnalysis <- highResolutionBeta %>%
    group_by(ActualBand,ReadBand) %>%
    summarise(ReadCount = sum(ReadCount),NumTranscript = n())
  ## Getting the band probabilities from the read counts in each band
  bandProb <- metaAnalysis %>%
    group_by(ActualBand) %>%
    mutate(sum(ReadCount))

  bandProb[3] <- round(bandProb[3]/bandProb[5],digits = 10)
  bandProb <- as.data.frame(bandProb)
  #print(bandProb)

  bandProb <- bandProb[c(1:3)]
  betaOutputFile <- paste(outputPath,"migration_probs.txt",sep = "")
  bandProbabilities <- as.data.frame(cast(data = bandProb,formula = ActualBand~ReadBand,value.var = "ReadCount", na.rm = TRUE, fill = 0))
  #print(bandProbabilities)
  bandProbabilities <- bandProbabilities[-1]
  bandProbabilities <- format(bandProbabilities, scientific = FALSE)
  write.table(bandProbabilities,file=betaOutputFile, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
  ### finished with high resolution band probabilites ###


}


betaEstimator(uniqueMappingFile, readCutoff, outputPath)
