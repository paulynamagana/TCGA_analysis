### DOWNLOAD METHYLATION

library(TCGAbiolinks)
# Obs: The data in the legacy database has been aligned to hg19


query.general <- GDCquery(
  project = "TCGA-LUSC", 
  legacy = TRUE,
  data.category = "DNA methylation")

general = getResults(query.general)

summary(factor(general$data_type))


################## more specific

query.raw <- GDCquery(
  project = "TCGA-LUSC",
  data.category = "Raw microarray data", 
  data.type = "Raw intensities", 
  legacy = TRUE,
  platform = "Illumina Human Methylation 450", 
  file.type = "idat",
  experimental.strategy = "Methylation array")


View(getResults(query.raw))

results <- getResults(query.raw)

#how many data_types we have
summary(factor(results$data_type))

summary(factor(results$sample_type))


# download results query.raw
GDCdownload(query = query.raw)



library(tidyverse)
library(dplyr)


sheet <- query.raw[[1]][[1]]


colnames(sheet)[colnames(sheet) == 'sample_type'] <- 'Sample_Group'
targets <- sheet %>%
  mutate(Sample_Group = str_replace(Sample_Group, " ", "_"))

#create Basename column it should not include _Red.idat
library(stringr)

dataDirectory <- "./GDCdata/TCGA-LUSC/legacy/Raw_microarray_data/Raw_intensities"



##### create basename
basename_folder <- sub("_...\\.idat", "", targets$file_name)
targets$Basename <- paste(dataDirectory, basename_folder, sep="/")

##extract strings
sentrixid <- sub("_R.....", "", basename_folder)
sentrix_position<- sub( ".*_", "", basename_folder)

#######
targets$Slide <- sentrixid
targets$Array <- sentrix_position
targets$Sample_Label <- targets$Sample_Group


duplicated(targets$Basename)
targets <-targets[!duplicated(targets$Basename), ]


##review this code########################
##GDC does not have a SampleSheet.csv file, so have to create one
SampleSheet <- targets[,c("id", "file_name", "data_category", "platform", "file_id", "experimental_strategy", "project", "Sample_Label", "Array", "Slide", "Basename")]

####save samplesheet
write.csv(SampleSheet, "./GDCdata/TCGA-LUSC/legacy/Raw_microarray_data/Raw_intensities/SampleSheet.csv", row.names = FALSE)







