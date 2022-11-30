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



##GDC does not have a SampleSheet.csv file, so have to create one
SampleSheet <- results[,c("id", "data_format", "file_name", "data_category", "platform", "file_id", "experimental_strategy", "project", "sample_type")]



summary(factor(SampleSheet$platform))
summary(factor(SampleSheet$sample_type))

write.csv(SampleSheet, "./GDCdata/TCGA-LUSC/legacy/Raw_microarray_data/Raw_intensities/SampleSheet.csv", row.names = TRUE)







