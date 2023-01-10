#______________loading packages__________________#
library("TCGAbiolinks")
library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library("IlluminaHumanMethylation450kmanifest")
library("minfi")
library("limma")
library("missMethyl") # Can take a short time...
library("DMRcate")
library("Gviz")
library("ggplot2")
library("RColorBrewer")
library("edgeR")
library("stringr")
library("readr")
library("tidyverse")



# loading saved session: Once you saved your data, you can load it into your environment: 
data.met = readRDS(file = "data.met.RDS")
# met matrix
met <- as.data.frame(SummarizedExperiment::assay(data.met))
# clinical data
clinical <- data.frame(data.met@colData)








































