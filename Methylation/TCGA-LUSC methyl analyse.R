#### process can be chunked to four main parts
############1########  Processing raw data
#Quality check
#Alignment and post-alignment processing
#Methylation calling
#Filtering bases

##############2 ####### Exploratory analysis
#Clustering
#PCA

#########3############Finding interesting regions
#Differential methylation
#Methylation segmentation

##############4###################Annotating interesting regions
#Nearest genes
#Annotation with other genomic features
#Integration with other quantitative genomics data





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
library("readr")
library("tidyverse")
library("edgeR")
library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library("IlluminaHumanMethylation450kmanifest")
library("stringr")

###################### NOT CURRENTLY USING
# look for all files
dataDirectory <- "./GDCdata/TCGA-LUSC/legacy/Raw_microarray_data/Raw_intensities"
list.files(dataDirectory, recursive = TRUE, full.names = TRUE)

### folder and dir for plots
dir.create("./meth_plots")
save_dir <- "./meth_plots/"


####################### SAMPLE SHEET##################################
##GDC does not have a SampleSheet.csv file, so created one when downloaded the files
# look for the file
list.files(dataDirectory, pattern = "*.csv")

#load SampleSheet.csv file
SampleSheet <- read_csv("GDCdata/TCGA-LUSC/legacy/Raw_microarray_data/Raw_intensities/SampleSheet.csv")


#file
#6264509100/6264509100_R01C01_Grn.idat
#Basename
## 1  C:/Users/Chrit/Documents/R/win-library/3.4/methylationArrayAnalysis/extdata/6264509100/6264509100_R01C01
met <- read.metharray.exp(targets=SampleSheet, verbose=TRUE, recursive =TRUE)
met

dim(met)# 622399 412

############################################################
##### downloading from TCGA ######################
#_________ increasing the memory limit____________#
#memory.limit(size = 28000)

# DNA methylation aligned to hg19
#query_met <- GDCquery(project= "TCGA-LUSC", 
#                      data.category = "DNA methylation", 
#                      platform = "Illumina Human Methylation 450", 
#                      legacy = TRUE)
#GDCdownload(query_met)



#putting files togathers
#data.met <- GDCprepare(query.raw)
#Processing  IDATs with Sesame - http://bioconductor.org/packages/sesame/
#Running opensesame - applying quality masking and nondetection masking (threshold P-value 0.05)
#Please cite: doi: 10.1093/nar/gky691 and 10.1093/nar/gkt090

#message("Raw data has ", nrow(data.met), " probes")
#message("Raw data has ", ncol(data.met), " samples")

#GDCprepare uses openSesame from sesame package with default arguments.
#GDCprepare transforms the downloaded data into a summarizedExperiment object
#or a data frame. If SummarizedExperiment is set to TRUE, TCGAbiolinks
#will add to the object sub-type information, which was defined by The
#Cancer Genome Atlas (TCGA) Research Network reports (the full list of papers
#can be seen in TCGAquery_subtype section in TCGAbiolinks vignette),
#and clinical information



#if you get error
##Error in stopAndCache(title):
#| File idatSignature needs to be cached to be used in sesame.
#| Please make sure you have updated ExperimentHub and try
#| > sesameDataCache()
#| to retrieve and cache needed sesame data.

## solution
##Each sesame data from ExperimentHub is accessible through the sesameDataGet interface.
#It should be noted that all data must be pre-cached to local disk before they can be used.
#This design is to prevent conflict in annotation data caching and remove internet dependency.
#Caching needs only be done once per sesame/sesameData installation. One can cache data using

#sesameDataCache()


###save file ###
#saveRDS(object = data.met,
 #       file = "data.met.RDS",
  #      compress = FALSE)

#get the 450 annotation data
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
head(ann450k)

# give the samples descriptive names
##change names
targets$ID <- paste(targets$Sample_Group,targets$id,sep=".")
sampleNames(met) <- targets$ID

#check samplenames
met


#######################################QUALITY CONTROL###############
###
##
##


# calculate the detection p-values
detP <- detectionP(met)
head(detP)



# examine mean detection p-values across all samples to identify any failed samples
png(paste0(save_dir, "Mean_detection_pvalues.png"), width=1000, height=750)
pal <- brewer.pal(8,"Dark2")
par(mfrow=c(1,2))
barplot(colMeans(detP), col=pal[factor(targets$Sample_Group)], las=2, 
        cex.names=0.8, ylab="Mean detection p-values")
abline(h=0.05,col="red")
legend("topleft", legend=levels(factor(targets$Sample_Group)), fill=pal,
       bg="white")

barplot(colMeans(detP), col=pal[factor(targets$Sample_Group)], las=2, 
        cex.names=0.8, ylim=c(0,0.002), ylab="Mean detection p-values")
abline(h=0.05,col="red")
legend("topleft", legend=levels(factor(targets$Sample_Group)), fill=pal, 
       bg="white")
dev.off()


#### The minfi qcReport function generates many other useful quality control plots.
#The minfi vignette describes the various plots and how they should be interpreted in
#detail. Generally, samples that look poor based on mean detection p-value will
#also look poor using other metrics and it is usually advisable to exclude them from further analysis.
qcReport(met, sampNames=targets$ID, sampGroups=targets$Sample_Group, 
         pdf="qcReport.pdf")

#We will use a detection p-value cutoff of >0.05, which removes the birth sample from further analysis.
#remove poor quality samples from rgSet
keep <- colMeans(detP) < 0.05
rgSet <- met[,keep]
dim(rgSet) # 622399    412


# remove poor quality samples from targets data
targets <- targets[keep,]
targets[,1:5]

# remove poor quality samples from detection p-value table
detP <- detP[,keep]
dim(detP) # 485512    412




#####################################NORMALISATION ############################

#preprocessFunnorm (Fortin et al. 2014) function is most appropriate for datasets
#with global methylation differences such as cancer/normal or vastly different tissue types, 
#If there exist global biological methylation differences between your samples, as for
#instance a dataset with cancer and normal samples, or a dataset with different tissues/cell
#types, use the preprocessFunnorm function as it is aimed for such datasets.

mSetSq  <- preprocessFunnorm(rgSet)




































####################################################################START FROM HERE################

# loading saved session: Once you saved your data, you can load it into your environment: 
data.met = readRDS(file = "data.met.RDS")
# met matrix
met <- as.data.frame(SummarizedExperiment::assay(data.met))
# clinical data
clinical <- data.frame(data.met@colData)

### MDS ###
pal <- brewer.pal(8,"Dark2")
png(paste0(save_dir, "MDS plot after GDCprepare.png"), width=1000, height=750)
plotMDS(met, top=1000, gene.selection="common", 
        col=pal[factor(clinical$sample_type)])
legend("topright", legend=levels(factor(clinical$sample_type)), text.col=pal,
       bg="white", cex=0.7)
dev.off()



#density
png(paste0(save_dir, "densities_after_GDCprepare.png"), width=1200, height=850)
plotDensities(met, legend=FALSE, main= "Densities sesame data")
dev.off() 


#############plot mean b values
# remove probes with NA (similar to na.omit)
tcga_na_filtered <- data.met[rowSums(is.na(assay(data.met))) == 0,]


df <- data.frame(
  "Sample.mean" = colMeans(assay(tcga_na_filtered), na.rm = TRUE),
  "groups" = tcga_na_filtered$sample_type
)


library(ggpubr)
png(paste0(save_dir, "mean_methylation_groups.png"), width=1200, height=850)
ggpubr::ggboxplot(
  data = df,
  y = "Sample.mean",
  x = "groups",
  color = "groups",
  add = "jitter",
  ylab = expression(paste("Mean DNA methylation (", beta, "-values)")),
  xlab = ""
) + stat_compare_means() 
dev.off()

##############______inspecting methylation data______#############################

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

dim(met) # 314845    412

## remove probes that match chromosome  X and Y 
keep <- !(row.names(met) %in% ann450k$Name[ann450k$chr %in% c("chrX","chrY")])
table(keep)
met <- met[keep, ]

dim(met) ##  308731    412
rm(keep) # remove no further needed probes.


## remove SNPs overlapped probe
table(is.na(ann450k$Probe_rs))
# probes without snp
no.snp.probe <- ann450k$Name[is.na(ann450k$Probe_rs)]

snp.probe <- ann450k[!is.na(ann450k$Probe_rs), ]
#snps with maf <= 0.05
snp5.probe <- snp.probe$Name[snp.probe$Probe_maf <= 0.05]


# filter met
met <- met[row.names(met) %in% c(no.snp.probe, snp5.probe), ]

dim(met) #283815    412

#remove no-further needed dataset
rm(no.snp.probe, probe, probe.na, snp.probe, snp5.probe)

## Removing probes that have been demonstrated to map to multiple places in the genome.
# list adapted from https://www.tandfonline.com/doi/full/10.4161/epi.23470
crs.reac <- read.csv("cross_reactive_probe.chen2013.csv")
crs.reac <- crs.reac$TargetID[-1]

# filter met
met <- met[ -which(row.names(met) %in% crs.reac), ]
dim(met) # 282641    412

#density
png(paste0(save_dir, "densities_after_filtering.png"), width=1200, height=850)
plotDensities(met, legend=FALSE, main= "Densities after filtering data")
dev.off() 

##################save files#################################################
#save as bval
bval <- met
## converting beta values to m_values
mval <- t(apply(met, 1, function(x) log2(x/(1-x))))


#######______________saving/loading_____________________############
# save data sets
saveRDS(mval, file = "mval.RDS", compress = FALSE)
saveRDS (bval, file = "bval.RDS", compress = FALSE)
#load datasets
#mval <- readRDS("mval.RDS")
#bval <- readRDS("bval.RDS")

#densities Bval
png(paste0(save_dir, "densities_bval.png"), width=1200, height=850)
plotDensities(bval, legend=FALSE, main= "Densities Bval data")
dev.off()

#densities mval
png(paste0(save_dir, "densities_mval.png"), width=1200, height=850)
plotDensities(mval, legend=FALSE, main= "Densities Mval data")
dev.off()



##############--------- Differential methylation analysis -------#########
table(clinical$sample_type, clinical$tissue_or_organ_of_origin)

##                     Lower lobe, lung Lung, NOS Main bronchus Middle lobe, lung Overlapping lesion of lung Upper lobe, lung
#Primary Tumor                    128        19             5                11                          8              199
#Solid Tissue Normal               17         3             0                 0                          1               21

ddddd <- qcReport(bval)





# filtering and grouping
clinical <- clinical[, c("barcode", "sample_type", "site_of_resection_or_biopsy")]
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



















