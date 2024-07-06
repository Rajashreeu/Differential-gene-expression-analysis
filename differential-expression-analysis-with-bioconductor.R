# installing Biconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

n#install packages from Bioconductor
BiocManager::install(c('Biobase','limma','geneplotter','enrichplot'))

BiocManager::install('EnhancedVolcano')

BiocManager::install('clusterProfiler')

#install pheatmap
install.packages('pheatmap')

# loading CRAN and Bioconductor packages
library(Biobase)
library(limma)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(geneplotter)
library(pheatmap)
library(enrichplot)
library(tidyr)
library(EnhancedVolcano)
library(clusterProfiler)
#Step 7 -
#First, we load the normalized expression assay, the phenotype data and the feature annotation data for this dataset.
exprsData <- read.delim("~/R PROGRAMMING/data_R/GSE27272Norm_exprs.txt")
phenoData <- read.delim("~/R PROGRAMMING/data_R/GSE27272Norm_phenoData.txt")
featureData <- read.delim("~/R PROGRAMMING/data_R/GSE27272Norm_featureData.txt")

#Step 8 
View(head(exprsData))
View(head(phenoData))
View(head(featureData))

#Step 9 -
#Variable for smoking status
phenoData['characteristics_ch1.1']

#Step 10 -
#Creating the an ExpressionSet object with all attributes

GSE27272_Eset<-ExpressionSet(as.matrix(exprsData))
View(head(GSE27272_Eset))
pData(GSE27272_Eset)<-phenoData
View(head(pData(GSE27272_Eset)))
featureData(GSE27272_Eset) <- as(featureData,"AnnotatedDataFrame")
View(head(featureData(GSE27272_Eset)))
#Exploratory Graph
#Step 11 -
#Before applying hypothesis testing on the data we should examine exploratory graphs like PCA and heatmaps to assesss our data.
fig <- function(width, heigth){
  options(repr.plot.width = width, repr.plot.height = heigth)
}

#Step 12 -
#Create a PCA plot examine the variation in the data by the phenotype variable of interest.
# Creates PCA Plots
fig(12,8)
GSE27272_exprs <- Biobase::exprs(GSE27272_Eset)
View(GSE27272_exprs)
PCA <- prcomp(t(GSE27272_exprs), scale = FALSE)
percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])
dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],Phenotype = Biobase::pData(GSE27272_Eset)$sex)

ggplot(dataGG, aes(PC1, PC2)) +
  geom_point(aes(colour = Phenotype)) +
  ggtitle("PCA plot of the GSE27272") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5,size=25, face='bold'),
        axis.text.x = element_text(size=12, face='bold'),
        axis.text.y = element_text(size=12, face='bold'),
        axis.title.x = element_text(size=18, face='bold'),
        axis.title.y = element_text(size=18, face='bold'),
        legend.title = element_text(size=18, face='bold'),
        legend.text = element_text(size=18) ) +
  scale_color_manual(values = c("hotpink", "deepskyblue"))

#Step 13 -
#Plotting a heatmap to examine the sample to sample distances and to see how well the samples cluster to sex.
annotation_for_heatmap <- data.frame(Phenotype = Biobase::pData(GSE27272_Eset)$sex)
row.names(annotation_for_heatmap) <- row.names(pData(GSE27272_Eset))

dists <- as.matrix(dist(t(GSE27272_exprs), method = "manhattan"))
rownames(dists) <- row.names(pData(GSE27272_Eset))
hmcol <- rev(colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))(255))
colnames(dists) <- NULL
diag(dists) <- NA
ann_colors <- list(
  Phenotype = c(female = "hotpink", male = "deepskyblue"))
pheatmap(dists, col = (hmcol),
         annotation_row = annotation_for_heatmap,
         annotation_colors = ann_colors,
         legend = TRUE,
         treeheight_row = 0,
         legend_breaks = c(min(dists, na.rm = TRUE),
                           max(dists, na.rm = TRUE)),
         
         legend_labels = (c("small distance", "large distance")),
         main = "Clustering heatmap for the GSE27272 samples")

#Step 14 -
#Filtering Data

# Filters the ExpressionSet (which includes the feature data and the expression data)
# to the genes that are not present in the Y chromosome
GSE27272_noY <-GSE27272_Eset[GSE27272_Eset@featureData@data$CHR!="Y",]
View(GSE27272_noY)

#Step 15 -
#Hypothesis Testing
design <- model.matrix(~0+phenoData$sex)
colnames(design) <- c("female","male")
GSE27272_samples <-
  as.character(phenoData$geo_accession)
rownames(design) <- GSE27272_samples
design <- model.matrix(~0+phenoData$sex)
colnames(design) <- c("female","male")
GSE27272_samples <-
  as.character(phenoData$geo_accession)
rownames(design) <- GSE27272_samples


#Step 16 -
contrast_matrix <- makeContrasts(female-male, levels= design)
#contrast_matrix <- makeContrasts(non_smoker-smoker, levels=design)

GSE27272_fit <- eBayes(contrasts.fit(lmFit(GSE27272_noY,
                                           
                                           design = design ),
                                     contrast_matrix))

#Step 17 -

table_GSE27272 <- topTable(GSE27272_fit, number = Inf,confint = TRUE)
hist(table_GSE27272$P.Value, col = brewer.pal(3, name = "Set2")[1],
     main = "Female vs Male - GSE27272", xlab = "p-values")
View(table_GSE27272)
#Step 18 -

GSE27272_Results <- data.frame(Ensembl_IDs= table_GSE27272$Ensembl_IDs,
                               
                               Entrez_IDs=table_GSE27272$Entrez_IDs,
                               Symbol=table_GSE27272$Symbols,Log2FC= table_GSE27272$logFC,
                               pvalue=table_GSE27272$P.Value,
                               adj.pvalue=table_GSE27272$adj.P.Val,
                               CI.L=table_GSE27272$CI.L,
                               CI.R = table_GSE27272$CI.R,
                               t=table_GSE27272$t,stringsAsFactors = FALSE)

head(GSE27272_Results)


#Step 19 -
#Volcano Plots

volcano_names <- ifelse(abs(GSE27272_fit$coefficients)>=1,
                        as.character(GSE27272_fit$genes$Symbols), NA)

volcanoplot(GSE27272_fit, coef = 1L, style = "p-value", highlight = 100,
            names = volcano_names,
            xlab = "Log2 Fold Change", ylab = NULL, pch=16, cex=0.35)

#Step 20 -
EnhancedVolcano(GSE27272_Results,
                lab = as.character(table_GSE27272$Symbols),
                x = 'Log2FC',
                title=" GSE27272 Volcano Plot Female vs Male",
                y = 'pvalue')

#Step 21 -
sigGenes <- GSE27272_Results[ GSE27272_Results$adj.pvalue < 0.05 &
                                
                                !is.na(GSE27272_Results$adj.pvalue) &
                                
                                abs(GSE27272_Results$Log2FC) > 1, ]

sigGenes

#Step 22 -
#We use the function enrichKEGG that will perform an enrichment analysis using the Kyoto Encyclopedia of Genes and Genomes (KEGG) database.
sigGenes <- GSE27272_Results$Entrez_IDs[ GSE27272_Results$adj.pvalue < 0.05 &
                                           
                                           !is.na(GSE27272_Results$adj.pvalue) &
                                           abs(GSE27272_Results$Log2FC) > 1 ]

sigGenes <- na.exclude(sigGenes)
kk <- enrichKEGG(gene = sigGenes, organism = 'hsa')
head(kk, n=10)







