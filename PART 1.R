#GUID: 2503736R Name: Alejandra Rodriguez Sosa
#ASSIGNMENT RNA-seq PART 1

#Input the data
#set experiment design
# Make dispersion plots for both objects and compare
# Make rlog-based PCA plots for both objects and compare
# Perform differential expression on gene-level for all contrasts
# Compare MA (MvA) plots for standard NULL hypothesis (LFC=0) and NULL hypothesis of LFC<2

setwd("E:/RNA-NGS/ASSIGNMENT/part 1/REPORT")
library(ggplot2)
library(DESeq2)
library(magrittr)

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

exp_design=data.frame(id=c("s1","s2","s3","s4","s5","s6","s7","s8","s9","s10","s11","s12"), group = rep(c("A","B","C"),each=4))
write.csv(exp_design, "Design_Exp.csv", row.names=F)

countTable <- read.csv("gene_count_matrix.csv",row.names=1)
countTranscripts <- read.csv("transcript_count_matrix.csv", row.names = 1)
colTable <- read.csv("Design_Exp.csv",row.names=1)

#Ordering the count gene and transcript tables by name
col_order <- c("s1","s2","s3","s4","s5","s6","s7","s8","s9","s10","s11","s12")
countTableO <- countTable[, col_order]
countTranscriptsO <- countTranscripts[, col_order]

ddsG <- DESeq2::DESeqDataSetFromMatrix(countData=countTableO, colData=colTable, design=~group)
ddsT <- DESeq2::DESeqDataSetFromMatrix(countData=countTranscriptsO, colData=colTable, design=~group)
#DESeq2::counts()

#remove all - 0 columns to reduce memory size of dds and increase speed of transformation
ddsG_filtered = ddsG[which(rowSums(DESeq2::counts(ddsG))>0),]
ddsT_filtered = ddsT[which(rowSums(DESeq2::counts(ddsT))>0),]
#
##GENE##
#
#Differential expression analysis over the count-matrix input data with the function DESeq()
ddsG_DA <- DESeq2::DESeq(ddsG_filtered)

#dispersion plot to see the effect of dispersion shrinkage and export to png
png("E:/RNA-NGS/ASSIGNMENT/part 1/dispersionG.png", height = 500, width = 600)
dispersionG <-estimateDispersions(ddsG_DA, fitType = c("parametric"))
print(plotDispEsts(dispersionG))
dev.off()

#regularised-log transformation produces transformed log2 data to remove dependence of the variance 
#within-group variability (if blind=FALSE)
#assay is used to extract the matrix of normalised values (and head to inspect)
rld <- rlog(ddsG_DA, blind = FALSE)

head(assay(vsd), 3)

#               s1       s2        s3        s4        s5        s6        s7        s8        s9       s10       s11       s12
#Cers6   11.127149 10.90619 11.017347 10.919659 11.394174 11.179202 11.460337 11.573910 11.081421 11.087375 11.130647 10.900590
#Gm14023  8.221945  8.10786  8.069994  8.083956  8.171544  7.889029  8.210672  8.122517  8.271664  8.196575  8.135350  8.233409
#Bpifb1   7.889029  8.10786  9.107108  8.265720  8.543654  8.701101  8.074986  8.184178  8.361031  8.184553  8.379898  7.889029

png("E:/RNA-NGS/ASSIGNMENT/part 1/PCA_G.png", height = 500, width = 600)
PCA_G = plotPCA(rld, intgroup="group")
print(PCA_G)
dev.off()

#results extracts the table with the log2fold, p and p.adj for each group comparison
ddsG_resAB <- DESeq2::results(ddsG_DA, contrast = c("group", "A","B"))
ddsG_resAC <- DESeq2::results(ddsG_DA, contrast = c("group", "A","C"))
ddsG_resBC <- DESeq2::results(ddsG_DA, contrast = c("group", "B","C"))

#MAplot#
#
# A vs B
#
resAvsB2 <- results(ddsG_DA,contrast=c("group","A","B"),lfcThreshold=2, altHypothesis = "greaterAbs")
resAvsB_sort2 = resAvsB2[order(resAvsB2$padj),]
resAvsB_sig2 = subset(resAvsB_sort2, resAvsB_sort2$padj<0.001)
summary(resAvsB_sig2)

resAvsB0 <- results(ddsG_DA,contrast=c("group","A","B"),lfcThreshold=0)
resAvsB_sort0 = resAvsB0[order(resAvsB0$padj),]
resAvsB_sig0 = subset(resAvsB_sort0, resAvsB_sort0$padj<0.001)
summary(resAvsB_sig0)

#The DESeq2 package results in an outputted MA plot that is not flexible to customise, so I've implemented it using ggplot2.
#DESeq2::plotMA(resAvsB_sig0)
#DESeq2::plotMA(resAvsB_sig2)

#In order to create a ggplot we need a data frame and not a dds object. We obtain two different data frames for our MA plot.
AvB_0 = as.data.frame(resAvsB_sig0)
AvB_2 = as.data.frame(resAvsB_sig2)

#The data come from two different data frames so we specify with the ggplot aes() the axis, then we add the data and colours in different geoms.
png("E:/RNA-NGS/ASSIGNMENT/part 1/geneMAplot_AvB.png", height = 500, width = 600)
MAplot_AvsB = ggplot(NULL,aes(x = log10(baseMean) , y= log2FoldChange)) + 
  geom_point(data = AvB_0, color = "turquoise") +
  geom_point(data = AvB_2, color = "purple") +  
  geom_hline(yintercept = 0, color = "red", size=0.5) + 
  geom_hline(yintercept=2,size = 1, color="dark green") +
  geom_hline(yintercept=-2, size = 1, color="dark green") +
  labs(title = "MA plot: Significant genes, A vs B", subtitle = "LFC=0(blue)and LFC<2(purple)")+
  theme_Alej()
print(MAplot_AvsB)
dev.off()
#
# A vs C
#
resAvsC2 <- results(ddsG_DA,contrast=c("group","A","C"),lfcThreshold=2, altHypothesis = "greaterAbs")
resAvsC_sort2 = resAvsC2[order(resAvsC2$padj),]
resAvsC_sig2 = subset(resAvsC_sort2, resAvsC_sort2$padj<0.001)
summary(resAvsC_sig2)

resAvsC0 <- results(ddsG_DA,contrast=c("group","A","C"),lfcThreshold=0)
resAvsC_sort0 = resAvsC0[order(resAvsC0$padj),]
resAvsC_sig0 = subset(resAvsC_sort0, resAvsC_sort0$padj<0.001)
summary(resAvsC_sig0)

AvC_0 = as.data.frame(resAvsC_sig0)
AvC_2 = as.data.frame(resAvsC_sig2)

png("E:/RNA-NGS/ASSIGNMENT/part 1/geneMAplot_AvC.png", height = 500, width = 600)
MAplot_AvsC = ggplot(NULL,aes(x = log10(baseMean) , y= log2FoldChange)) + 
  geom_point(data = AvC_0, color = "turquoise") +
  geom_point(data = AvC_2, color = "purple") +  
  geom_hline(yintercept = 0, color = "red", size=0.5) + 
  geom_hline(yintercept=2,size = 1, color="dark green") +
  geom_hline(yintercept=-2, size = 1, color="dark green") +
  labs(title = "MA plot: Significant genes, A vs C", subtitle = "LFC=0(blue)and LFC<2(purple)")+
  theme_Alej()
print(MAplot_AvsC)
dev.off()
#
# B vs C
#
resBvsC2 <- results(ddsG_DA,contrast=c("group","B","C"),lfcThreshold=2, altHypothesis = "greaterAbs")
resBvsC_sort2 = resBvsC2[order(resBvsC2$padj),]
resBvsC_sig2 = subset(resBvsC_sort2, resBvsC_sort2$padj<0.001)
summary(resBvsC_sig2)

resBvsC0 <- results(ddsG_DA,contrast=c("group","B","C"),lfcThreshold=0)
resBvsC_sort0 = resBvsC0[order(resBvsC0$padj),]
resBvsC_sig0 = subset(resBvsC_sort0, resBvsC_sort0$padj<0.001)
summary(resBvsC_sig0)

BvC_0 = as.data.frame(resBvsC_sig0)
BvC_2 = as.data.frame(resBvsC_sig2)

png("E:/RNA-NGS/ASSIGNMENT/part 1/geneMAplot_BvC.png", height = 500, width = 600)
MAplot_BvsC = ggplot(NULL,aes(x = log10(baseMean) , y= log2FoldChange)) + 
  geom_point(data = BvC_0, color = "turquoise") +
  geom_point(data = BvC_2, color = "purple") +  
  geom_hline(yintercept = 0, color = "red", size=0.5) + 
  geom_hline(yintercept=2,size = 1, color="dark green") +
  geom_hline(yintercept=-2, size = 1, color="dark green") +
  labs(title = "MA plot: Significant genes, B vs C", subtitle = "LFC=0(blue)and LFC<2(purple)")+
  theme_Alej()
print(MAplot_BvsC)
dev.off()


##TRANSCRIPT##
#running the differential expression analysis over the count-matrix input data with the function DESeq()
ddsT_DA <- DESeq2::DESeq(ddsT_filtered)

#dispersion plot to see the effect of dispersion shrinkage
png("E:/RNA-NGS/ASSIGNMENT/part 1/dispersionT.png", height = 500, width = 600)
dispersionT <- estimateDispersions(ddsT_DA, fitType = c("parametric"))
print(plotDispEsts(dispersionT))
dev.off()

#regularised-log transformation, transformed log2 data. Within-group variability (if blind=FALSE)
rldT <- rlog(ddsT_DA, blind = FALSE)

head(assay(vsdT), 3)

#               s1       s2       s3       s4       s5       s6       s7       s8       s9      s10      s11      s12
#NM_026615 9.678745 9.764990 9.666272 9.571129 9.875623 9.827339 9.746903 9.801209 9.703627 9.734380 9.786080 9.644320
#NM_026613 9.025492 8.999296 9.071709 9.014954 9.152245 9.352904 9.198882 9.062500 9.045766 9.085110 9.080844 9.060745
#NM_016866 9.153159 9.280506 9.294834 9.316602 9.271032 9.351581 9.445668 9.263372 9.247042 9.294964 9.293570 9.285126

png("E:/RNA-NGS/ASSIGNMENT/part 1/PCA_T.png", height = 500, width = 600)
PCA_T = plotPCA(rldT, intgroup="group")
print(PCA_T)
dev.off()

#results extracts the table with the log2fold, p and p.adj for each group comparison
ddsT_resAB <- DESeq2::results(ddsT_DA, contrast = c("group", "A","B"))
ddsT_resAC <- DESeq2::results(ddsT_DA, contrast = c("group", "A","C"))
ddsT_resBC <- DESeq2::results(ddsT_DA, contrast = c("group", "B","C"))

#MAplots#
#
#Group A vs B
#
TresAvsB2 <- results(ddsT_DA,contrast=c("group","A","B"),lfcThreshold=2, altHypothesis = "greaterAbs")
TresAvsB_sort2 = TresAvsB2[order(TresAvsB2$padj),]
TresAvsB_sig2 = subset(TresAvsB_sort2, TresAvsB_sort2$pvalue<0.001)

TresAvsB0 <- results(ddsT_DA,contrast=c("group","A","B"),lfcThreshold=0)
TresAvsB_sort0 = TresAvsB0[order(TresAvsB0$padj),]
TresAvsB_sig0 = subset(TresAvsB_sort0, TresAvsB_sort0$pvalue<0.001)

#The DESeq2 package results in an outputted MA plot that is not flexible to customise, so I've implemented it using ggplot2.
#DESeq2::plotMA(TresAvsB_sig0)
#DESeq2::plotMA(TresAvsB_sig2)

#In order to create a ggplot we need a data frame and not a dds object. We obtain two different data frames for our MA plot.
TAvB_0 = as.data.frame(TresAvsB_sig0)
TAvB_2 = as.data.frame(TresAvsB_sig2)

#The data come from two different data frames so we specify with the ggplot aes() the axis, then we add the data and colours in different geoms.
png("E:/RNA-NGS/ASSIGNMENT/part 1/transcriptMAplot_AvB.png", height = 500, width = 600)
TMAplot_AvsB = ggplot(NULL,aes(x = log10(baseMean) , y= log2FoldChange)) + 
  geom_point(data = TAvB_0, color = "turquoise") +
  geom_point(data = TAvB_2, color = "purple") +  
  geom_hline(yintercept = 0, color = "red", size=0.5) + 
  geom_hline(yintercept=2,size = 1, color="dark green") +
  geom_hline(yintercept=-2, size = 1, color="dark green") +
  labs(title = "MA plot: Significant transcripts, A vs B", subtitle = "LFC=0(blue)and LFC<2(purple)")+
  theme_Alej()
print(TMAplot_AvsB)
dev.off()
#
#Group A vs C
#
TresAvsC2 <- results(ddsT_DA,contrast=c("group","A","C"),lfcThreshold=2, altHypothesis = "greaterAbs")
TresAvsC_sort2 = TresAvsC2[order(TresAvsC2$padj),]
TresAvsC_sig2 = subset(TresAvsC_sort2, TresAvsC_sort2$pvalue<0.001)

TresAvsC0 <- results(ddsT_DA,contrast=c("group","A","C"),lfcThreshold=0)
TresAvsC_sort0 = TresAvsC0[order(resAvsC0$padj),]
TresAvsC_sig0 = subset(TresAvsC_sort0, TresAvsC_sort0$pvalue<0.001)

TAvC_0 = as.data.frame(TresAvsC_sig0)
TAvC_2 = as.data.frame(TresAvsC_sig2)

png("E:/RNA-NGS/ASSIGNMENT/part 1/transcriptMAplot_AvC.png", height = 500, width = 600)
TMAplot_AvsC = ggplot(NULL,aes(x = log10(baseMean) , y= log2FoldChange)) + 
  geom_point(data = TAvC_0, color = "turquoise") +
  geom_point(data = TAvC_2, color = "purple") +  
  geom_hline(yintercept = 0, color = "red", size=0.5) + 
  geom_hline(yintercept=2,size = 1, color="dark green") +
  geom_hline(yintercept=-2, size = 1, color="dark green") +
  labs(title = "MA plot: Significant transcripts, A vs C", subtitle = "LFC=0(blue)and LFC<2(purple)")+
  theme_Alej()
print(TMAplot_AvsC)
dev.off()
#
#Group B vs C
#
TresBvsC2 <- results(ddsT_DA,contrast=c("group","B","C"),lfcThreshold=2, altHypothesis = "greaterAbs")
TresBvsC_sort2 = TresBvsC2[order(TresBvsC2$padj),]
TresBvsC_sig2 = subset(TresBvsC_sort2, TresBvsC_sort2$pvalue<0.001)

TresBvsC0 <- results(ddsT_DA,contrast=c("group","B","C"),lfcThreshold=0)
TresBvsC_sort0 = TresBvsC0[order(TresBvsC0$padj),]
TresBvsC_sig0 = subset(TresBvsC_sort0, TresBvsC_sort0$pvalue<0.001)

TBvC_0 = as.data.frame(TresBvsC_sig0)
TBvC_2 = as.data.frame(TresBvsC_sig2)

png("E:/RNA-NGS/ASSIGNMENT/part 1/transcriptMAplot_BvC.png", height = 500, width = 600)
TMAplot_BvsC = ggplot(NULL,aes(x = log10(baseMean) , y= log2FoldChange)) + 
  geom_point(data = TBvC_0, color = "turquoise") +
  geom_point(data = TBvC_2, color = "purple") +  
  geom_hline(yintercept = 0, color = "red", size=0.5) + 
  geom_hline(yintercept=2,size = 1, color="dark green") +
  geom_hline(yintercept=-2, size = 1, color="dark green") +
  labs(title = "MA plot: Significant transcripts, B vs C", subtitle = "LFC=0(blue)and LFC<2(purple)")+
  theme_Alej()
print(TMAplot_BvsC)
dev.off()

#Resources:
#https://hbctraining.github.io/DGE_workshop/lessons/04_DGE_DESeq2_analysis.html

