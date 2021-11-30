#GUID: 2503736R Name: Alejandra Rodriguez Sosa
#ASSIGNMENT RNA-seq PART 3

setwd("E:/RNA-NGS/ASSIGNMENT/Thomas/Assessment")

library(Seurat)
library(ggplot2)
library(magrittr)
library(dplyr)

#Data represents counts (log transformed)
Expression_file <- read.table("Expression.txt", sep = "\t", header = TRUE)

#Generate two different plots that represent the expression of TREM2 between the two conditions, showing
#if the gene has different expression levels between the two conditions

#Creating my own theme to be used in graphics later on
theme_Alej <- function() {
  theme (
    plot.title = element_text(size = 18, face = "bold", colour = "purple"),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16)
  )
}
Plot_TREM = subset(Expression_file, Gene == "TREM2")

mycolours = c("red","light blue")

#The boxplot was first used to check which kind of graphic would it lead to but it was too clumped to assess the magnitude of the difference
#in the expression

png("E:/RNA-NGS/ASSIGNMENT/Thomas/Assessment/TREM2_box.png", height = 500, width = 600)
TREM2_boxplot = ggplot(Plot_TREM, aes(x= Condition, y = Expression, fill = Condition)) + 
  geom_boxplot() + 
  stat_summary(fun=mean, geom="point", shape=17, size=2, color="magenta") +
  scale_fill_manual(values=mycolours) +
  labs(title = "Expression of TREM2 Remission vs ActiveRA")+
  theme_Alej()
print(TREM2_boxplot)
dev.off()

#The violin plot seemed better than the boxplot but it was not clear enough to assess the difference in expression

#Violin plot for assessing the difference in expression of the gene TREM2
png("E:/RNA-NGS/ASSIGNMENT/Thomas/Assessment/TREM2_violin.png", height = 500, width = 600)
TREM2_violinplot = ggplot(Plot_TREM, aes(x= Condition, y = Expression, fill = Condition)) + 
  geom_violin() + 
  scale_fill_manual(values=mycolours) + 
  stat_summary(fun=mean, geom="point", shape=17, size=2, color="green")+ 
  geom_jitter(shape = 19, aes(color = Condition), alpha=0.2)+
  labs(title = "Expression of TREM2 Remission vs ActiveRA")+
  theme_Alej()
print(TREM2_violinplot)
dev.off()

#Frequency polygon plot (similar to histogram) for assessing the difference in expression of the gene TREM2
png("E:/RNA-NGS/ASSIGNMENT/Thomas/Assessment/TREM2_freq.png", height = 500, width = 600)
TREM2_freq = ggplot(Plot_TREM, aes(x= Expression, colour = Condition)) + 
  geom_freqpoly(binwidth=0.1, size = 0.65) + 
  scale_fill_manual(values=mycolours) +
  labs(title = "Expression of TREM2 Remission vs ActiveRA", subtitle = "Frequency polygon plot")+
  theme_Alej()
print(TREM2_freq)
dev.off()

#density plot for assessing the difference in expression of the gene TREM2
png("E:/RNA-NGS/ASSIGNMENT/Thomas/Assessment/TREM2_density.png", height = 500, width = 600)
TREM2_density = ggplot(Plot_TREM, aes(x=Expression, fill = Condition)) + 
  geom_density(alpha = 0.5)+
  labs(title = "Expression of TREM2 Remission vs ActiveRA", subtitle = "Density plot")+
  theme_Alej()
print(TREM2_density)
dev.off()

#t-test and obtained result
activeRA <- subset(Plot_TREM, Condition == "ActiveRA")
Remission <- subset(Plot_TREM, Condition == "Remission")

ttest = t.test(activeRA$Expression, Remission$Expression)

#Welch Two Sample t-test
#
#data:  activeRA$Expression and Remission$Expression
#t = -26.756, df = 7195, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.2395823 -0.2068721
#sample estimates:
#  mean of x  mean of y 
#0.07437191 0.29759911 

#SEURAT INTEGRATION ANALYSIS STARTS HERE#
#PART A) TASK 3
library(Seurat)
library(SeuratData)
library(cowplot)
library(ggplot2)

##ACTIVERA##
#load the data
acRAN.data <- Read10X(data.dir = "E:/RNA-NGS/ASSIGNMENT/Thomas/Assessment/acRAN_Seurat/GRCh38" )
# Initialize the Seurat object with the raw (non-normalized data).
acRAN <- CreateSeuratObject( counts = acRAN.data, project = "acRAN", min.cells = 3)
acRAN$sample <- "acRAN1"
acRAN$group <- "acRAN"
#QC and filtering cells by mitochondrial contamination (to filter low quality / dying cells)
acRAN[["percent.mt"]] <- PercentageFeatureSet(acRAN, pattern = "^MT-")
#VlnPlot(acRAN, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
acRAN <- subset(acRAN, subset = nFeature_RNA > 500 & nFeature_RNA < 2200 & percent.mt < 18)

# store mitochondrial percentage in object meta data
acRAN <- PercentageFeatureSet(acRAN, pattern = "^MT-", col.name = "percent.mt")

# run sctransform
acRAN <- SCTransform(acRAN, vars.to.regress = "percent.mt", verbose = FALSE)
#
#
##REMISSION##
#
#
#load the data
REM.data <- Read10X(data.dir = "E:/RNA-NGS/ASSIGNMENT/Thomas/Assessment/REM_Seurat/GRCh38")
# Initialize the Seurat object with the raw (non-normalized data).
REM <- CreateSeuratObject( counts = REM.data, project = "REM", min.cells = 3)
REM$sample <- "REM1"
REM$group <- "REM"
#QC and filtering cells by mitochondrial contamination (to filter low quality / dying cells)
REM[["percent.mt"]] <- PercentageFeatureSet(REM, pattern = "^MT-")
#VlnPlot(REM, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
REM <- subset(REM, subset = nFeature_RNA > 500 & nFeature_RNA < 2200 & percent.mt < 18)

# store mitochondrial percentage in object meta data
REM <- PercentageFeatureSet(REM, pattern = "^MT-", col.name = "percent.mt")

# run sctransform
REM <- SCTransform(REM, vars.to.regress = "percent.mt", verbose = FALSE)
#
#
# Combine the two objects
#
#
dim=15

my.anchors <- FindIntegrationAnchors(object.list = list(acRAN,REM), dims = 1:dim)
Combined <- IntegrateData(anchorset = my.anchors, dims = 1:dim)
DefaultAssay(Combined) <- "integrated"
Combined <- ScaleData(Combined, verbose = FALSE)
Combined <- RunPCA(Combined, npcs = dim, verbose = FALSE)
Combined <- RunUMAP(Combined, reduction = "pca", dims = 1:dim)
Combined <- FindNeighbors(Combined, reduction = "pca", dims = 1:dim)
Combined <- FindClusters(Combined, resolution = 0.2)

png("E:/RNA-NGS/ASSIGNMENT/Thomas/Assessment/UMAP1.png", height = 500, width = 600)
p1 <- DimPlot(Combined, reduction = "umap", group.by = "sample")
print(p1)
dev.off()

png("E:/RNA-NGS/ASSIGNMENT/Thomas/Assessment/UMAP2.png", height = 500, width = 600)
p2 <- DimPlot(Combined, reduction = "umap", label = TRUE)
print(p2)
dev.off()

png("E:/RNA-NGS/ASSIGNMENT/Thomas/Assessment/CombUMAP.png", height = 500, width = 900)
CombUMAP <- CombinePlots(plots = list(p1, p2))
print(CombUMAP)
dev.off()

## find conserved markers for the clusters now
Combined.filt.markers <- FindAllMarkers(Combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Combined.filt.markers  %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

top2 <- Combined.filt.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
top2_list <- as.data.frame(top2[c(7)])
write.table(top2_list, file = "top2_list.txt", sep = "",row.names = FALSE, col.names = FALSE, quote = FALSE)

top10 <- Combined.filt.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
top10_list <- as.data.frame(top10[c(7)])
write.table(top10_list, file = "top10_list.txt", sep = "",row.names = FALSE, col.names = FALSE, quote = FALSE)

#top10_heat <- DoHeatmap(Combined.filt, features = top10$gene) + NoLegend()
png("E:/RNA-NGS/ASSIGNMENT/Thomas/Assessment/top10.png", height = 500, width = 900)
print(top10_heat)
dev.off()
#View(top10)

#identify celltypes
png("E:/RNA-NGS/ASSIGNMENT/Thomas/Assessment/annotatedUMAP.png", height = 500, width = 900)
Combined.anno <- RenameIdents(Combined,'0'="macrophage_C1QA",'1'="neutrophil",'2'="macrophage_LYZ",
                              '3'="T-reg",'4'="macrophage_ACKR1",'5'="eosinophil",
                              '6'="T-reg",'7'="T-reg",'8'="Lymphocyte B",
                              '9'="naive B-cell",'10'="basophil",'11'="macrophage_TAGLN", '12'="T-reg")

annotated_umap <- DimPlot(Combined.anno,label =TRUE,split.by = "sample")
print(annotated_umap)
dev.off


#PART B) TASK 3
# Filtering away other clusters
Combined.filt <- subset(Combined.anno, idents = c("macrophage_C1QA", "macrophage_LYZ","macrophage_ACKR1","macrophage_TAGLN"), invert = FALSE)
DimPlot(Combined.filt, reduction = "umap", split.by = "sample")

# re-clustering, without the other cell types
DefaultAssay(object = Combined.filt) <- "integrated"
# Run the standard workflow for visualization and clustering
Combined.filt <- ScaleData(object = Combined.filt, verbose = FALSE)
Combined.filt <- RunPCA(object = Combined.filt, npcs = dim, verbose = FALSE)
# t-SNE and Clustering
Combined.filt <- RunUMAP(object = Combined.filt, reduction = "pca", dims = 1:dim)
Combined.filt <- FindNeighbors(object = Combined.filt, reduction = "pca", dims = 1:dim)
Combined.filt <- FindClusters(Combined.filt, resolution = 0.15)

DimPlot(Combined.filt, reduction = "umap", split.by = "sample")
#DimPlot(object = Combined.filt, reduction = "umap")


# filtering CLUSTERS 0,3,4 away other clusters
Combined.filt <- subset(Combined.filt, idents = c("1","2"), invert = FALSE)
DimPlot(Combined.filt, reduction = "umap", split.by = "sample")

# re-clustering, without the other cell types
DefaultAssay(object = Combined.filt) <- "integrated"
# Run the standard workflow for visualization and clustering
Combined.filt <- ScaleData(object = Combined.filt, verbose = FALSE)
Combined.filt <- RunPCA(object = Combined.filt, npcs = dim, verbose = FALSE)
# t-SNE and Clustering
Combined.filt <- RunUMAP(object = Combined.filt, reduction = "pca", dims = 1:dim)
Combined.filt <- FindNeighbors(object = Combined.filt, reduction = "pca", dims = 1:dim)
Combined.filt <- FindClusters(Combined.filt, resolution = 0.15)

png("E:/RNA-NGS/ASSIGNMENT/Thomas/Assessment/macrophage_plot.png", height = 500, width = 900)
Macrophage_plot = DimPlot(Combined.filt, reduction = "umap", split.by = "sample")
#DimPlot(object = Combined.filt, reduction = "umap")
print(Macrophage_plot)
dev.off()

#### GENERATING GENE LIST FROM MACROPHAGES FOR GENE ENRICHMENT STUDY

saveRDS(Combined.filt,"2503736R_RA_Assign1.rds")

DefaultAssay(Combined.filt) <- "RNA"
theme_set(theme_cowplot())

Combined.filt$celltype.group <- paste(Idents(Combined.filt), Combined.filt$sample, sep = "_")
Combined.filt$celltype <- Idents(Combined.filt)

n="0"
check_sanity = FALSE
## diff expression\
Idents(Combined.filt) <- "celltype.group"
response_2 <- FindMarkers(Combined.filt, ident.1 = paste(n,"_acRAN1", sep=""), ident.2 = paste(n,"_REM1", sep=""), verbose = FALSE, min.pct = 0.25)

#This table stores the number of genes that are up-regulated in Remission with corrected p-values < 0.01
remission_2_sig <- subset(response_2,response_2$avg_logFC > 0 & response_2$p_val_adj < 0.01)
remission_list = row.names(remission_2_sig)
write.table(remission_list, file = "remission_list.txt", sep = "",row.names = FALSE, col.names = FALSE, quote = FALSE)

