### Analyse data

library("limma")
library("minfi")
library("RColorBrewer")
library("missMethyl") # Can take a short time...
library("minfiData")
library("Gviz")
library("sesame") #for analysis
library("DMRcate")
library("DMRcatedata")
library("stringr")
library("mCSEA")
library("readr")
library("tidyverse")
library("edgeR")
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(stringr)

# look for all files
dataDirectory <- "./GDCdata/TCGA-LUSC/legacy/Raw_microarray_data/Raw_intensities"
list.files(dataDirectory, recursive = TRUE, full.names = TRUE)


dir.create("./meth_plots")
save_dir <- "./meth_plots/"


####################### SAMPLE SHEET##################################
##GDC does not have a SampleSheet.csv file, so created one when downloaded the files
# look for the file
list.files(dataDirectory, pattern = "*.csv")

#load SampleSheet.csv file
targets <- read_csv("GDCdata/TCGA-LUSC/legacy/Raw_microarray_data/Raw_intensities/SampleSheet.csv", 
                        col_types = cols(...1 = col_skip()))

targets <- targets %>%
  mutate(sample_type = str_replace(sample_type, " ", "_"))
colnames(targets)[colnames(targets) == 'sample_type'] <- 'Sample_Group'

#print head
head(targets)


#create Basename column it should not include _Red.idat
library(stringr)

targets$Basename <- paste(dataDirectory, targets$file_name, sep="/")
basename_folder <- sub("_...\\.idat", "", targets$Basename)
sentrix <- sub("_...\\.idat", "", targets$file_name)
sentrixid <- sub("_R.....", "", sentrix)
sentrix_position<- sub( ".*_", "", sentrix)



targets$Basename <- basename_folder
targets$Slide <- sentrixid
targets$Array <- sentrix_position
targets$Sample_Label <- targets$Sample_Group

head(targets$Basename)

#file
#6264509100/6264509100_R01C01_Grn.idat
#Basename
## 1  C:/Users/Chrit/Documents/R/win-library/3.4/methylationArrayAnalysis/extdata/6264509100/6264509100_R01C01

met <- read.metharray.exp(targets=targets, recursive =TRUE)

duplicated(targets$Basename)


############################################################
##### downloading from TCGA ######################
#_________ increasing the memory limit____________#
memory.limit(size = 28000)

# DNA methylation aligned to hg19
query_met <- GDCquery(project= "TCGA-LUSC", 
                      data.category = "DNA methylation", 
                      platform = "Illumina Human Methylation 450", 
                      legacy = TRUE)
GDCdownload(query_met)



#putting files togathers
data.met <- GDCprepare(query.raw)

#if you get error
##Error in stopAndCache(title):
#| File idatSignature needs to be cached to be used in sesame.
#| Please make sure you have updated ExperimentHub and try
#| > sesameDataCache()
#| to retrieve and cache needed sesame data.

## solution
##Each sesame datum from ExperimentHub is accessible through the sesameDataGet interface.
#It should be noted that all data must be pre-cached to local disk before they can be used.
#This design is to prevent conflict in annotation data caching and remove internet dependency.
#Caching needs only be done once per sesame/sesameData installation. One can cache data using

#sesameDataCache()


###save file ###
saveRDS(object = data.met,
        file = "data.met.RDS",
        compress = FALSE)



# loading saved session: Once you saved your data, you can load it into your environment: 
data.met = readRDS(file = "data.met.RDS")
# met matrix
met <- as.data.frame(SummarizedExperiment::assay(data.met))
# clinical data
clinical <- data.frame(data.met@colData)



#___________inspecting methylation data_______________#

# get the 450k annotation data
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

## remove probes with NA
probe.na <- rowSums(is.na(met))
table(probe.na == 0)
#FALSE   TRUE 
# 170732 314845 

# chose those has no NA values in rows
probe <- probe.na[probe.na == 0]
met <- met[row.names(met) %in% names(probe), ]


## remove probes that match chromosome  X and Y 
keep <- !(row.names(met) %in% ann450k$Name[ann450k$chr %in% c("chrX","chrY")])
table(keep)
met <- met[keep, ]
rm(keep) # remove no further needed probes.



## remove SNPs overlapped probe
table (is.na(ann450k$Probe_rs))
# probes without snp
no.snp.probe <- ann450k$Name[is.na(ann450k$Probe_rs)]

snp.probe <- ann450k[!is.na(ann450k$Probe_rs), ]
#snps with maf <= 0.05
snp5.probe <- snp.probe$Name[snp.probe$Probe_maf <= 0.05]


# filter met
met <- met[row.names(met) %in% c(no.snp.probe, snp5.probe), ]

#remove no-further needed dataset
rm(no.snp.probe, probe, probe.na, snp.probe, snp5.probe)

## Removing probes that have been demonstrated to map to multiple places in the genome.
# list adapted from https://www.tandfonline.com/doi/full/10.4161/epi.23470


crs.reac <- read.csv("cross_reactive_probe.chen2013.csv")
crs.reac <- crs.reac$TargetID[-1]

# filter met
met <- met[ -which(row.names(met) %in% crs.reac), ]
bval <- met

## converting beta values to m_values
mval <- t(apply(met, 1, function(x) log2(x/(1-x))))


#______________saving/loading_____________________#
# save data sets
#saveRDS(mval, file = "mval.RDS", compress = FALSE)
#saveRDS (bval, file = "bval.RDS", compress = FALSE)
#mval <- readRDS("mval.RDS")
#bval <- readRDS("bval.RDS")

head(bval[,1:5])

densityPlot(mval, sampGroups=clinical$sample_type, main="Beta values", 
            legend=FALSE, xlab="Beta values")
legend("top", legend = levels(factor(targets$Sample_Group)), 
       text.col=brewer.pal(8,"Dark2"))












##############--------- Differential methylation analysis -------#########
table(clinical$sample_type, clinical$tissue_or_organ_of_origin)

##                     Lower lobe, lung Lung, NOS Main bronchus Middle lobe, lung Overlapping lesion of lung Upper lobe, lung
#Primary Tumor                    128        19             5                11                          8              199
#Solid Tissue Normal               17         3             0                 0                          1               21


# filtering and grouping
clinical <- clinical[, c("barcode", "sample_type", "tissue_or_organ_of_origin")]
clinical <- na.omit(clinical)
barcode <- clinical$barcode

# removing samples from meth matrixes
bval <- bval[, colnames(bval) %in% barcode]
mval <- mval[, colnames(mval) %in% barcode]

# Making sure about samples in clinical and matrixes and their order
table(colnames(mval) %in% row.names(clinical))
table(colnames(bval) %in% row.names(clinical))
#
all(row.names(clinical) == colnames(bval))
all(row.names(clinical) == colnames(mval))


#Making grouping variable
clinical$sample_type <- as.factor(clinical$sample_type)
#levels(clinical$sample_type)
clinical$sample_type <- relevel(clinical$sample_type, ref = "Solid Tissue Normal")


###########_____________ DMC analysis________________###########
design <- model.matrix(~ sample_type, data = clinical)
# fit the linear model 
fit <- lmFit(mval, design)
fit2 <- eBayes(fit)



#extracting significantly methylated probes
deff.meth = topTable(fit2, coef=ncol(design), sort.by="p",number = nrow(mval), adjust.method = "BY")


# Visualization
# plot the top 10 most significantly differentially methylated CpGs 
png(paste0(save_dir, "top 10 most significantly differentially methylated CpGs.png"), width=1000, height=750)
par(mfrow=c(2,5))
sapply(rownames(deff.meth)[1:10], function(cpg){
  plotCpg(bval, cpg=cpg, pheno=clinical$sample_type, ylab = "Beta values")
})
dev.off()

# if Error in .Call.graphics(C_palette2, .Call(C_palette2, NULL)) : 
#invalid graphics state
#In addition: Warning messages:

#then run dev.off()  to clear the plots#


#making a volcano plot
#making dataset
dat <- data.frame(foldchange = fit[["coefficients"]][,2], logPvalue =  -log10(fit2[["p.value"]][,2]))
dat$threshold <- as.factor(abs(dat$foldchange) < 0.4)

#Visualization
png(paste0(save_dir, "volcano_plot_meth.png"), width=1000, height=750)
cols <- c("TRUE" = "grey", "FALSE" = "blue")
ggplot(data=dat, aes(x=foldchange, y = logPvalue, color=threshold)) +
  geom_point(alpha=.6, size=1.2) +
  scale_colour_manual(values = cols) +
  geom_vline(xintercept = 0.4, colour="#990000", linetype="dashed") + 
  geom_vline(xintercept = - 0.4, colour="#990000", linetype="dashed") +
  theme(legend.position="none") +
  xlab("Fold Change") +
  ylab("-log10 p value") +
  theme_bw() +
  theme(legend.position = "none")
dev.off()




#######------------------------- Differentially methylated regions (DMR) analysis
# setting some annotation
myAnnotation <- cpg.annotate(object = mval,
                             datatype = "array", 
                             what = "M", 
                             analysis.type = "differential", 
                             design = design, 
                             contrasts = FALSE, 
                             coef = "sample_typePrimary Tumor", 
                             arraytype = "450K",
                             fdr = 0.001)
## Your contrast returned 184682 individually significant probes. We recommend the default setting of pcutoff in dmrcate().


str(myAnnotation)


















detP <- detectionP(data.met)
head(detP)
















