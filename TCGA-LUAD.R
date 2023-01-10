

# Load packages
library("TCGAbiolinks")
library("limma")
library("edgeR")
library("caret")
library("SummarizedExperiment")
library("gplots")
library("survival")
library("survminer")
library("RColorBrewer")
library("genefilter")
library("DESeq2")

######################
##############download ############ #####TCGA ##########

#extract all projects
GDCprojects = getGDCprojects()
head(GDCprojects[c("project_id", "name")] )

dplyr::filter(GDCprojects, grepl('Lung', name)) #projects for lung 
rm(GDCprojects) #get rid of data

TCGAbiolinks:::getProjectSummary("TCGA-LUAD") #get summary for LUAD


##query the data
query_TCGA = GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq") # some cases there might be more than 1 data type = gene expression, splice junstion quantification


#To visualize the query results in a more readable way, we can use the command getResults.
lusc_res = getResults(query_TCGA)
colnames(lusc_res)

#CHECK DATA
head(lusc_res)
summary(factor(lusc_res$sample_type))

###query all the data
query_TCGA_raw = GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  access= "open",
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  data.type = "Gene Expression Quantification",
  sample.type = c("Primary Tumor", "Solid Tissue Normal")) #double check


##download the data
GDCdownload(query_TCGA_raw, method = "api")
## samples have been downloaded


luad.count  <- GDCprepare(query = query_TCGA_raw,
                          save = TRUE,
                          save.filename = "TCGA_LUAD.rda", #save file as
                          summarizedExperiment = TRUE)
#Starting to add information to samples
#=> Add clinical information to samples
#=> Adding TCGA molecular information from marker papers
#=> Information will have prefix 'paper_' 









#get rid of unused data
#free up space
rm(query_TCGA, GDCprojects, lusc_res, query_TCGA_raw) #change if needed


#see dimensions
dim(luad.count)
#] 60660   598

##access associated clinical data
colnames(colData(luad.count))

#vital status
table(luad.count@colData$vital_status)
#Alive  Dead 
#379   219 


#tumour stage
table(luad.count@colData$ajcc_pathologic_stage)

##tumour classification
table(luad.count@colData$definition)
#Primary solid Solid Tissue Normal 
#539                  59 


#tissue of origin
table(luad.count@colData$tissue_or_organ_of_origin)

#gender
table(luad.count@colData$gender)

#race
table(luad.count@colData$race) # mostly white


#####gene expression eset
dim(assay(luad.count))     # gene expression matrices.
# 60660   598
head(assay(luad.count)[,1:10]) # expression of first 6 genes and first 10 samples
head(rowData(luad.count))     # ensembl id and gene id of the first 6 genes.


# Save the data as a file, if you need it later, you can just load this file
# instead of having to run the whole pipeline again

saveRDS(object = luad.count,
        file = "tcga_data.RDS",
        compress = FALSE)

#The data can be loaded with the following command
#tcga_data = readRDS(file = "tcga_data.RDS")

#extract clinical data
clinical_data = colData(luad.count)

#define groups
group = factor(clinical_data$definition)
levels(group)

#define Solid Tissue Normal as being the base or reference level.
group = relevel(group, ref="Solid Tissue Normal")

#create design matrix
design = model.matrix(~group)
head(design)


##remove genes, which have low amount of counts. 
dge = DGEList( # creating a DGEList object
  counts=assay(luad.count),
  samples=colData(luad.count),
  genes=as.data.frame(rowData(luad.count)))

##filtering counts
#at least 20 samples with a cpm of 0.5 or higher
#isexpr <- rowSums(cpm(dge)>0.5) >= 20
#table(isexpr)

#dge <- dge[isexpr,]
#dim(dge) # 26076   598




# filtering
keep = filterByExpr(dge, group=design) # defining which genes to keep
dge = dge[keep,,keep.lib.sizes=FALSE] # filtering the dge object
rm(keep) #  use rm() to remove objects from memory if you don't need them anymore
dim(dge) # 26757   598 without doing the cpm


####visualizations
dir.create("./plots") #create folder
save_dir <- "./plots/"
# density plot of raw read counts (log10)
dge$samples$lib.size


##library sizes
png(paste0(save_dir, "lib_sizes.png"), width=1000, height=900)
barplot(dge$samples$lib.size, names=colnames(dge), las=2)
# Add a title to the plot
title("Barplot of library sizes")
dev.off()


#analyse distributions
# Get log2 counts per million
logcounts <- cpm(dge,log=TRUE)
# Check distributions of samples using boxplots
png(paste0(save_dir, "boxplotslogcpms.png"), width=1000, height=900)
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (unnormalised)")
dev.off()




######################
#########################NORMALISE ################
#calculate normalization factors between libraries
dge = calcNormFactors(dge, method="TMM")
dge$samples$norm.factors


#normalisation by the method of trimmed mean of M-values (TMM)
#is performed using the calcNormFactors function in edgeR.
#The normalisation factors calculated here are used as a scaling factor for the library sizes. 
colnames(design)
colnames(design)<- gsub(" ", "_", colnames(design))
colnames(design) #verify underscores

# normalise the read counts with 'voom' function
#incorporate sample-level weights together with the abundance dependent weights estimated by voom
v = voom(dge, design, plot=TRUE)
v


#fit model
vfit = lmFit(v, design)
efit = eBayes(vfit, trend= TRUE)

#sumary tests
#Significance is defined using an adjusted p-value cutoff that is set at 5% by default. 
summary(decideTests(efit))





png(paste0(save_dir, "plotSA.png"), width=1000, height=900)
plotSA(efit, main="Final model: Mean-variance trend")
dev.off()


#topTable we can check the top10 genes classified as being differentially expressed. ’
# Show top genes
full_results <- topTable(efit, number=Inf)

#save data table
library(readr) 
write_csv(full_results, paste0(save_dir, path="_filtered_DE_results.csv"))




topGenes = topTable(efit, coef=2, sort.by="p", number=100)
write.csv(topGenes, paste0(save_dir, " top100genes.csv"))



#######################################################SURVIVAL #########
library("gplots")
library("survival")


clinical_data <- luad.count@colData
# 598 87



clin_df = clinical_data[clinical_data$definition == "Primary solid Tumor",
                   c("patient",
                     "vital_status",
                     "days_to_death",
                     "days_to_last_follow_up",
                     "gender",
                     "ajcc_pathologic_stage")]

names(clin_df)[names(clin_df) == "ajcc_pathologic_stage"] <- "tumor_stage"

# create a new boolean variable that has TRUE for dead patients
# and FALSE for live patients
clin_df$deceased = clin_df$vital_status == "Dead"

# create an "overall survival" variable that is equal to days_to_death
# for dead patients, and to days_to_last_follow_up for patients who
# are still alive
clin_df$overall_survival = ifelse(clin_df$deceased,
                                  clin_df$days_to_death,
                                  clin_df$days_to_last_follow_up)

# show first 10 samples
head(clin_df)


##The survival package provides us with an object,
##Surv, to form a dependent variable out of the overall_survival
##and deceased information:
Surv(clin_df$overall_survival, clin_df$deceased)
Surv(clin_df$overall_survival, clin_df$deceased) ~ clin_df$gender

# fit a survival model
fit = survfit(Surv(overall_survival, deceased) ~ gender, data=clin_df)
print(fit)



png(paste0(save_dir, "survival_gender.png"), width=700, height=700)
ggsurvplot(fit, data=clin_df)
dev.off()

##add pval
png(paste0(save_dir, "survival_gender_pval.png"), width=700, height=700)
ggsurvplot(fit,
           data=clin_df,
           xlab = "Time in days",   # customize X axis label.
           pval = TRUE, # add pvalue
           conf.int = TRUE) # Add confidence interval
dev.off()




###########################################
########survival genes
# let's extract the table of differential expression we got earlier
# print the first row, to see the gene name, the logFC value and the p-value
print(full_results[1, ])

#extract expression from v
d_mat = as.matrix(t(v$E))




# get the ensembl gene id of the first row
gene_id = "ENSG00000197208.6"
#SLC22A1_id = "ENSG00000175003.15"
#SLC22A4_id = "ENSG00000197208.6"
#SLC22A5_id = "ENSG00000197375.13"


# also get the common gene name of the first row
gene_name = "SLC22A4"


# visualize the gene expression distribution on the diseased samples
# versus the healthy samples
expr_diseased = d_mat[rownames(clin_df), gene_id]
expr_healthy = d_mat[setdiff(rownames(d_mat), rownames(clin_df)), gene_id]


##plot boxplot expression
png(paste0(save_dir, gene_name, "expression_gene.png"), width=700, height=700)
boxplot(expr_diseased, expr_healthy,
        names=c("Diseased", "Healthy"), main="Distribution of gene expression")
dev.off()



# get the expression values for the selected gene
clin_df$gene_value = d_mat[rownames(clin_df), gene_id]

# find the median value of the gene and print it
median_value = median(clin_df$gene_value)
print(median_value)



# divide patients in two groups, up and down regulated.
# if the patient expression is greater or equal to them median we put it
# among the "up-regulated", otherwise among the "down-regulated"
clin_df$gene = ifelse(clin_df$gene_value >= median_value, "UP", "DOWN")

# we can fit a survival model, like we did in the previous section
fit = survfit(Surv(overall_survival, deceased) ~ gene, data=clin_df)

# we can extract the survival p-value and print it
pval = surv_pvalue(fit, data=clin_df)$pval
print(pval)

# and finally, we produce a Kaplan-Meier plot
png(paste0(save_dir, gene_name, "survival_gene.png"), width=800, height=700)
ggsurvplot(fit,
           data=clin_df,
           pval = TRUE, # add pvalue
           risk.table=TRUE,
           xlab = "Time in days",   # customize X axis label.
           title=paste(gene_name),
           border = 'full',
           borderWidth = 1.0,
           borderColour = 'black')
dev.off()










#################new volcano ENHANCED ###########################
library(EnhancedVolcano)
#The default cut-off for log2FC is >|2|; the default cut-off for P value is 10e-6.


keyvals.colour <- ifelse(
  full_results$logFC < -1.5, 'royalblue',
  ifelse(full_results$logFC > 1.5, 'red',"grey"))
names(keyvals.colour)[keyvals.colour == 'red'] <- 'Up-regulated'
names(keyvals.colour)[keyvals.colour == 'grey'] <- 'Not-Significant'
names(keyvals.colour)[keyvals.colour == 'royalblue'] <- 'Down-regulated'



png(paste0(save_dir, "enhanced_volcano.png"), width=700, height=900)
EnhancedVolcano(full_results,
                lab = "",
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.001,
                FCcutoff = 1.5,
                title = "TCGA: cohort: LUNG LUAD (n=558)",
                subtitle = "Differential expression: Primary solid Solid (n=539) Tissue Normal (n=59) ",
                cutoffLineType = 'twodash',
                cutoffLineWidth = 0.8,
                pointSize = 4.0,
                labSize = 6.0,
                colAlpha = 1,
                colCustom = keyvals.colour,
                legendLabels=c('Not sig.','Log (base 2) FC','p-value',
                               'p-value & Log (base 2) FC'),
                col = c('grey', 'grey', 'grey', 'red3'),
                legendPosition = 'bottom',
                legendLabSize = 16,
                legendIconSize = 5.0,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                border = 'full',
                borderWidth = 1.0,
                borderColour = 'black')
dev.off()
###############################




###MORE VISUALISATIONS
##############plotMD
png(paste0(save_dir, "plot_MD.png"), width=700, height=900)
plotMD(v)
dev.off()


##plot dispersion
y <- estimateDisp(dge, design, robust=TRUE)
png(paste0(save_dir, "plot_BC_dispersion.png"), width=700, height=900)
plotBCV(y)
dev.off()









#####enrichment for pathways
## The pull function from dplyr is used to extract a particular column
library(org.Hs.eg.db); library(dplyr)
pathway_genes <- AnnotationDbi::select(org.Hs.eg.db,
                                       keys = "GO:0030198",
                                       keytype = "GO",
                                       columns="ENSEMBL") %>% pull(ENSEMBL)


go_table <- mutate(full_results,
                   inPathway = gene_id %in% pathway_genes,
                   isDE = adj.P.Val < 0.05 & abs(logFC) >1)
go_table

#inspect the genes in pathways
table(go_table$inPathway, go_table$isDE)

#chisquare
chisq.test(table(go_table$inPathway, go_table$isDE))





##analysis with cluterprofile
library(clusterProfiler)
library(dplyr)


##extract data for slc22a plot
gene_id_Ensemble <- full_results$gene_id
length(gene_id_Ensemble)

eset <- v$E
head(eset)
head(colnames(eset))
head(rownames(eset))

###check dimensiones
dim(eset)
# 26757 genes
# 598 samples

as_data <- as.data.frame(eset, SKIP=0)
dim(as_data) # check for no changes

library(tidyr); library(ggplot2); library(stringr); library(dplyr) #packages

pdata <- clinical_data

##fix table with genesnames and probes
as_data$probe <- gene_id_Ensemble
as_data$Symbol <- v$genes$gene_name
rownames(as_data) <- NULL


data_long <- gather(as_data, IDATfile, log_fold,
                    "TCGA-78-7156-01A-11R-2039-07":"TCGA-73-4675-01A-01R-1206-07", factor_key=FALSE)
head(data_long)


#######

#change names to match the functions
names(pdata)[names(pdata) == "barcode"] <- "IDATfile"
pdata <- as.data.frame(pdata) #convert
targetinfo <- dplyr::select(pdata, "IDATfile", "definition") #extract columns

#merge dataong and targetinfo for groups
data_long <- merge(data_long, targetinfo, by  = "IDATfile")
names(data_long)[names(data_long) == "definition"] <- "Group" # change name


#save data
write.csv(data_long, "data_long.csv")



###filter for genes of interest
SLC22A1_table <- data_long %>%
  filter(Symbol == "SLC22A1")

SLC22A4_table <- data_long %>%
  filter(Symbol == "SLC22A4")

SLC22A5_table <- data_long %>%
  filter(Symbol == "SLC22A5")


#function for plotting per gene
plot_gene <- function(data, title){
  ggplot(data, aes(x= Group, log_fold)) +
    geom_boxplot(outlier.shape = NA, color= "black",fill= c("gray60", "gray33")) +
    geom_jitter(width=0.08, height = 0.5, size= 2.0, color= "black") +
    ggtitle(data$genes, subtitle = "TCGA LUAD") +# We'll make this a jitter plot
    ylab("Expression log2 scale") + 
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size=24), 
          axis.title = element_text(size = 20), # font size of axis
          axis.text.x = element_text(size=16), #font size of x ticks
          axis.text.y = element_text(size=14),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=1))+ # font size of y ticks
    scale_y_continuous(breaks = round(seq(min(data$log_fold), max(data$log_fold), by = 0.5),1))
  ggsave(title, width = 15,height = 25, units="cm")
}




## save plots ##
plot_gene(SLC22A1_table, paste0(save_dir,"SLC22A1_expression.png"))
plot_gene(SLC22A4_table, paste0(save_dir,"SLC22A4_expression.png"))
plot_gene(SLC22A5_table, paste0(save_dir,"SLC22A5_expression.png"))






######################plot all SLC22 genes in boxplot
#table
SLC22 <- data_long %>%
  filter(stringr::str_detect(Symbol, "SLC22A"))



plot_all_genes <- function(data, title){
  ggplot(data, aes(x= Symbol ,log_fold, fill=Group))  +
    geom_boxplot(outlier.shape = NA, color= "black", position="dodge") +
    ggtitle("TCGA-LUAD ", subtitle="cohort: Primary solid Solid (n=539) Tissue Normal (n=59)") +
    ylab("Log2(counts+1)") + 
    xlab("") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size=20), 
          axis.title = element_text(size = 18), # font size of axis
          axis.text.x = element_text(size=16, angle = 90), #font size of x ticks
          axis.text.y = element_text(size=12))+ # font size of y ticks
    scale_y_continuous(breaks = round(seq(min(data$log_fold), max(data$log_fold), by = 1.5),1))+
    scale_fill_manual(values=c("gray60", "gray33"))
  ggsave(title, width = 35,height = 20, units="cm")
}


#plot and save
plot_all_genes(SLC22, paste0(save_dir,"SLC22_expression.png"))









############################ Go enrichment
gene_id<-sub("\\..*","", gene_id_Ensemble)
full_results$ENSEMBL <- gene_id


sigGenes <- full_results %>% 
  filter(adj.P.Val < 0.05, !is.na(ENSEMBL)) %>% pull(ENSEMBL)

enrich_go <- enrichGO(
  gene= sigGenes,
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL",
  ont = "BP",
  universe = gene_id,
  qvalueCutoff = 0.05,
  readable=TRUE
)

#turn into datafram
enrich_go %>% data.frame

##top 20 categories
png(paste0(save_dir, "ennrich_GO_top20.png"), width=700, height=900)
dotplot(enrich_go,
        showCategory=20)
dev.off()



enrich_go <- enrichplot::pairwise_termsim(enrich_go)
png(paste0(save_dir, "ennrich_GO_top20_pairwise.png"), width=700, height=900)
emapplot(enrich_go)
dev.off()


png(paste0(save_dir, "ennrich_GO_top20_pairwiseupset.png"), width=700, height=900)
enrichplot::upsetplot(enrich_go)
dev.off()



##################adding GENENAME

anno <- AnnotationDbi::select(org.Hs.eg.db, 
                              columns=c("ENSEMBL","GENENAME"),
                              keys=gene_id,
                              keytype = "ENSEMBL")

names(anno)[names(anno) == "ENSEMBL"] <- "gene_id"

#get rid of the extra for ensembl ids
full_results$gene_id <- sub("\\..*","", full_results$gene_id)

##export results
export_results <- merge(full_results, anno, by="gene_id")
write.csv(export_results, "exportedresults.csv")
