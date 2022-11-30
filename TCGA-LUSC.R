library(data.table)


generate_count_mat <- function(path, pattern) {
  files = list.files(path, pattern, full.names = TRUE, recursive=TRUE, include.dirs=TRUE)
  mat = as.data.frame(do.call(cbind, lapply(files, function(x) fread(x, stringsAsFactors = FALSE))))
  mat <- mat[-c(1:4),]
  #gene_type <- as.character(mat[,3])
  rownames(mat) = mat[,1]
  mat = as.data.frame(mat[, seq(4, ncol(mat), 9)])
  #mat$gene_type <- gene_type
  #mat <- subset(mat, mat$gene_type == "protein_coding")
  #mat <- mat[,-c(ncol(mat))]
  return(mat)
}




# stage raw counts
dir <- "./GDCdata/TCGA-LUSC/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification"
mrna_counts <- generate_count_mat(dir, "//.rna_seq.augmented_star_gene_counts.tsv$")





library(tidyverse)

file_names <- list.files(dir, "//.rna_seq.augmented_star_gene_counts.tsv$", full.names = FALSE, recursive = TRUE, include.dirs = FALSE)
file_names <- sub(".*/", "", file_names)
samplesheet <- read.table("gdc_sample_sheet.2022-04-27.tsv", header=T, sep="\t")
samplesheet <- samplesheet[match(file_names, samplesheet$File.Name),]
colnames(mrna_counts) <- samplesheet$Sample.ID
meta <- subset(samplesheet, select=c(Sample.ID, Sample.Type))
rownames(meta) <- NULL
meta <- column_to_rownames(mrna_meta, var="Sample.ID")




























