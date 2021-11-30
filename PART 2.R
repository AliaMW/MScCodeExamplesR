#GUID: 2503736R Name: Alejandra Rodriguez Sosa
#ASSIGNMENT RNA-seq PART 2

setwd("A:/Dropbox (Personal)/MSc Bioinformatics/Semester 2/RNA-NGS/ASSIGNMENT/data")
library(ggplot2)
library(pheatmap)
library(dplyr)
library(ggfortify)#to use autoplot
library(factoextra)
library(amap)
library(reshape) #for melt
library(ggpubr)

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
#png("E:/RNA-NGS/ASSIGNMENT/part 1/PLOTNAME.png", height = 500, width = 600)
#print(THEPLOT)
#dev.off()

#Three groups with 14 samples each: Flu (FLU), Health control (HC), Common cold (CC)
#Calculated using GRCh38 release 91

CC_HC <- read.table("DE_CC_vs_HC.csv", header = TRUE, sep = '\t')  
FLU_CC <- read.table("DE_FLU_vs_CC.csv", header = TRUE, sep = '\t')
FLU_HC <- read.table("DE_FLU_vs_HC.csv", header = TRUE, sep = '\t')
Expression_matrix <- read.table("EM.csv", header = TRUE, sep = '\t')
sample_sheet <- read.table("sample_sheet.csv", header = TRUE, sep = '\t')

#Downloading gene annotations for GRCh38 release 91
#As Ensembls biomart is currently running at reduced capacity and the correct annotations cant be downloaded please either:
#1) Goto biomart not in archive and use the latest version of the annotation (99 I think). 
#2) As annotations (gene symbol) is mainly aesthetic its not strictly needed to be able to answer the assignment questions, technically you can skip them.

biomart <- read.table("biomart.tsv", header=TRUE, sep='\t')

#add mean expression column; the object where function is applied, the dimension 1 for rows 2 for cols, name of function
Expression_matrix$mean_values = apply(Expression_matrix[2:43],1,mean)
write.table(Expression_matrix, file="Expression_matrix.tsv", quote=FALSE, sep="\t",col.names = NA)

#add a -log10p column to each differential expression table

logCC_HC <- -log(CC_HC["p"],10)
CC_HC = data.frame(CC_HC, logCC_HC)
names(CC_HC)[5] <- "log10p"

logFLU_CC <- -log(FLU_CC["p"],10)
FLU_CC = data.frame(FLU_CC,logFLU_CC)
names(FLU_CC)[5] <- "log10p"

logFLU_HC <- -log(FLU_HC["p"],10)
FLU_HC = data.frame(FLU_HC,logFLU_HC)
names(FLU_HC)[5] <- "log10p"

#add a column to each d.e. table that flags significance: p.adj<0.05, log2fold>1
CC_HC$sig = as.factor(CC_HC$p.adj < 0.05 & abs(CC_HC$log2fold) > 1.0)
FLU_CC$sig = as.factor(FLU_CC$p.adj < 0.05 & abs(FLU_CC$log2fold) > 1.0)
FLU_HC$sig = as.factor(FLU_HC$p.adj < 0.05 & abs(FLU_HC$log2fold) > 1.0)

write.table(CC_HC,file="CC_HC.csv", quote=FALSE, sep="\t",col.names=NA)
write.table(FLU_CC,file="FLU_CC.csv", quote=FALSE, sep="\t",col.names=NA)
write.table(FLU_HC,file="FLU_HC.csv", quote=FALSE, sep="\t",col.names=NA)

#make a master file combining diff.exp.tables with the expression matrix table AND save it to disk
merge1 <- merge(CC_HC, FLU_CC, by="ID")
merge2 <- merge(merge1, FLU_HC, by="ID")
master_file = merge(Expression_matrix, merge2, by="ID")
#colnames(master_file) #just to check the whole list of names in our master file

names(master_file)[45:59] <- c(
  "log2fold_CC_HC","p_CC_HC", "p.adj_CC_HC","log10_CC_HC","sig_CC_HC",
  "log2fold_FLU_CC","p_FLU_CC", "p.adj_FLU_CC","log10_FLU_CC","sig_FLU_CC","log2fold_FLU_HC",
  "p_FLU_HC", "p.adj_FLU_HC","log10_FLU_HC","sig_FLU_HC"
)
#rename only the columns with repeated names (coming from the differential expression matrices)
write.table(master_file,file="master_file.tsv",quote=FALSE,sep="\t",col.names = NA)

#US:UP and SIGNIFICATIVE / DS:DOWN and SIGNIFICATIVE / NC:NO CHANGE
CCHC_US_FLUCC_US = subset(master_file, log2fold_CC_HC > 0 & sig_CC_HC == "TRUE" & log2fold_FLU_CC > 0 & sig_FLU_CC == "TRUE")
CCHC_US_FLUCC_DS = subset(master_file, log2fold_CC_HC > 0 & sig_CC_HC == "TRUE" & log2fold_FLU_CC < 0 & sig_FLU_CC == "TRUE")
CCHC_US_FLUCC_NC = subset(master_file, log2fold_CC_HC > 0 & sig_CC_HC == "TRUE" & sig_FLU_CC == "FALSE")

CCHC_DS_FLUCC_DS = subset(master_file, log2fold_CC_HC < 0 & sig_CC_HC == "TRUE" & log2fold_FLU_CC < 0 & sig_FLU_CC == "TRUE")
CCHC_DS_FLUCC_US = subset(master_file, log2fold_CC_HC < 0 & sig_CC_HC == "TRUE" & log2fold_FLU_CC > 0 & sig_FLU_CC == "TRUE")
CCHC_DS_FLUCC_NC = subset(master_file, log2fold_CC_HC < 0 & sig_CC_HC == "TRUE" & sig_FLU_CC == "FALSE")

CCHC_NC_FLUCC_NC = subset(master_file, sig_CC_HC == "FALSE" & sig_FLU_CC == "FALSE")
CCHC_NC_FLUCC_DS = subset(master_file, sig_CC_HC == "FALSE" & log2fold_FLU_CC < 0 & sig_FLU_CC == "TRUE")
CCHC_NC_FLUCC_US = subset(master_file, sig_CC_HC == "FALSE" & log2fold_FLU_CC > 0 & sig_FLU_CC == "TRUE")
##
CCHC_US_FLUHC_US = subset(master_file, log2fold_CC_HC > 0 & sig_CC_HC == "TRUE" & log2fold_FLU_HC > 0 & sig_FLU_HC == "TRUE")
CCHC_US_FLUHC_DS = subset(master_file, log2fold_CC_HC > 0 & sig_CC_HC == "TRUE" & log2fold_FLU_HC < 0 & sig_FLU_HC == "TRUE")
CCHC_US_FLUHC_NC = subset(master_file, log2fold_CC_HC > 0 & sig_CC_HC == "TRUE" & sig_FLU_HC == "FALSE")

CCHC_DS_FLUHC_DS = subset(master_file, log2fold_CC_HC < 0 & sig_CC_HC == "TRUE" & log2fold_FLU_HC < 0 & sig_FLU_HC == "TRUE")
CCHC_DS_FLUHC_US = subset(master_file, log2fold_CC_HC < 0 & sig_CC_HC == "TRUE" & log2fold_FLU_HC > 0 & sig_FLU_HC == "TRUE")
CCHC_DS_FLUHC_NC = subset(master_file, log2fold_CC_HC < 0 & sig_CC_HC == "TRUE" & sig_FLU_HC == "FALSE")

CCHC_NC_FLUHC_NC = subset(master_file, sig_CC_HC == "FALSE" & sig_FLU_HC == "FALSE")
CCHC_NC_FLUHC_DS = subset(master_file, sig_CC_HC == "FALSE" & log2fold_FLU_HC < 0 & sig_FLU_HC == "TRUE")
CCHC_NC_FLUHC_US = subset(master_file, sig_CC_HC == "FALSE" & log2fold_FLU_HC > 0 & sig_FLU_HC == "TRUE")
##
FLUCC_US_FLUHC_US = subset(master_file, log2fold_FLU_CC > 0 & sig_FLU_CC == "TRUE" & log2fold_FLU_HC > 0 & sig_FLU_HC == "TRUE")
FLUCC_US_FLUHC_DS = subset(master_file, log2fold_FLU_CC > 0 & sig_FLU_CC == "TRUE" & log2fold_FLU_HC < 0 & sig_FLU_HC == "TRUE")
FLUCC_US_FLUHC_NC = subset(master_file, log2fold_FLU_CC > 0 & sig_FLU_CC == "TRUE" & sig_FLU_HC == "FALSE")

FLUCC_DS_FLUHC_US = subset(master_file, log2fold_FLU_CC < 0 & sig_FLU_CC == "TRUE" & log2fold_FLU_HC > 0 & sig_FLU_HC == "TRUE")
FLUCC_DS_FLUHC_DS = subset(master_file, log2fold_FLU_CC < 0 & sig_FLU_CC == "TRUE" & log2fold_FLU_HC < 0 & sig_FLU_HC == "TRUE")
FLUCC_DS_FLUHC_NC = subset(master_file, log2fold_FLU_CC < 0 & sig_FLU_CC == "TRUE" & sig_FLU_HC == "FALSE")

FLUCC_NC_FLUHC_NC = subset(master_file, sig_FLU_CC == "FALSE" & sig_FLU_HC == "FALSE")
FLUCC_NC_FLUHC_DS = subset(master_file, sig_FLU_CC == "FALSE" & log2fold_FLU_HC < 0 & sig_FLU_HC == "TRUE")
FLUCC_NC_FLUHC_US = subset(master_file, sig_FLU_CC == "FALSE" & log2fold_FLU_HC > 0 & sig_FLU_HC == "TRUE")

##PLOTTING STARTS HERE##
#To assess if expression relates to differential expression, sequencing deep enough and which are most differential genes: MA PLOT

#MA plots
png("A:/Dropbox (Personal)/MSc Bioinformatics/Semester 2/RNA-NGS/ASSIGNMENT/data", height = 500, width = 600)
MAplot_CC_HC = ggplot(CC_HC, aes(x = log10(master_file$mean_values) , y= log2fold, color = sig))+ geom_point() + theme_Alej()
print(MAplot_CC_HC)
dev.off()

png("A:/Dropbox (Personal)/MSc Bioinformatics/Semester 2/RNA-NGS/ASSIGNMENT/data", height = 500, width = 600)
MAplot_FLU_HC = ggplot(FLU_HC, aes(x = log10(master_file$mean_values) , y= log2fold, color = sig))+ geom_point()
print(MAplot_FLU_HC)
dev.off()

png("A:/Dropbox (Personal)/MSc Bioinformatics/Semester 2/RNA-NGS/ASSIGNMENT/data", height = 500, width = 600)
MAplot_FLU_CC = ggplot(FLU_CC, aes(x = log10(master_file$mean_values) , y= log2fold, color = sig))+ geom_point()
print(MAplot_FLU_CC)
dev.off()


mycolours = c("black","magenta","grey")

#Volcano plots - statistically significant gene expression changes - INCLUDE IN REPORT

png("A:/Dropbox (Personal)/MSc Bioinformatics/Semester 2/RNA-NGS/ASSIGNMENT/data", height = 500, width = 600)
VOLC_CC_HC = ggplot(CC_HC, aes(x = log2fold, y = log10p, colour=sig, group = sig)) + 
geom_point(size = 1) + 
labs(title = "Volcano plot Common Cold vs Health Control", x = "log2fold change", y = "-log10 p-value") +
scale_color_manual(values=mycolours,na.translate = FALSE)
print(VOLC_CC_HC)
dev.off()

png("A:/Dropbox (Personal)/MSc Bioinformatics/Semester 2/RNA-NGS/ASSIGNMENT/data", height = 500, width = 600)
VOLC_FLU_HC = ggplot(FLU_HC, aes(x = log2fold, y = log10p, colour=sig)) + 
geom_point(size = 1) + 
labs(title = "Volcano plot Flu vs Health Control", x = "log2fold change", y = "-log10 p-value") +
  scale_color_manual(values=mycolours,na.translate = FALSE)
print(VOLC_FLU_HC)
dev.off()

png("A:/Dropbox (Personal)/MSc Bioinformatics/Semester 2/RNA-NGS/ASSIGNMENT/data", height = 500, width = 600)
VOLC_FLU_CC = ggplot(FLU_CC, aes(x = log2fold, y = log10p, colour=sig)) + 
geom_point(size = 1) + 
labs(title = "Volcano plot Flu vs Common Cold", x = "log2fold change", y = "-log10 p-value") +
  scale_color_manual(values=mycolours,na.translate = FALSE)
print(VOLC_FLU_CC)
dev.off()

#Expression boxplot for most differential expressed genes 
#first sort the data in master_file by increasing p-val
#then extract columns and transpose data
#make expression boxplot: 2 groups being compared on the X axis, expression of the gene on y axis; sort by ascending p

#TOP3 MOST DIFFERENTIALLY EXPRESSED GENES COLD_HEALTH: X3210 X10906 X9182

DE_CCHC_SORT1 <- master_file[order(master_file$p_CC_HC, decreasing = FALSE, na.last = TRUE),]
DE_CCHC_SORT1 = DE_CCHC_SORT1[c(2:29)]
DE_CCHC_SORT_TRANS1 = data.frame(t(DE_CCHC_SORT1))
png("A:/Dropbox (Personal)/MSc Bioinformatics/Semester 2/RNA-NGS/ASSIGNMENT/data/BOX_CCHC1.png", height = 500, width = 600)
BOX_CCHC1 = ggplot(DE_CCHC_SORT_TRANS1, aes(x=sample_sheet$GROUP[c(1:28)], y=X3210, fill = sample_sheet$GROUP[c(1:28)])) + 
  geom_boxplot() + theme_Alej() + labs(title = "CCvsHC gene ENSG00000108950", x = "Group", y = "Gene expression X3210")+ coord_cartesian(ylim = c(0, 1000))
print(BOX_CCHC1)
dev.off()

DE_CCHC_SORT2 <- master_file[order(master_file$p_CC_HC, decreasing = FALSE, na.last = TRUE),]
DE_CCHC_SORT2 = DE_CCHC_SORT2[c(2:29)]
DE_CCHC_SORT_TRANS2 = data.frame(t(DE_CCHC_SORT2))
png("A:/Dropbox (Personal)/MSc Bioinformatics/Semester 2/RNA-NGS/ASSIGNMENT/data/BOX_CCHC2.png", height = 500, width = 600)
BOX_CCHC2 = ggplot(DE_CCHC_SORT_TRANS2, aes(x=sample_sheet$GROUP[c(1:28)], y=X10906, fill = sample_sheet$GROUP[c(1:28)])) + 
  geom_boxplot() + theme_Alej() + labs(title = "CCvsHC gene ENSG00000170439",x = "Group", y = "Gene expression X10906")+ coord_cartesian(ylim = c(0, 300))
print(BOX_CCHC2)
dev.off()

DE_CCHC_SORT3 <- master_file[order(master_file$p_CC_HC, decreasing = FALSE, na.last = TRUE),]
DE_CCHC_SORT3 = DE_CCHC_SORT3[c(2:29)]
DE_CCHC_SORT_TRANS3 = data.frame(t(DE_CCHC_SORT3))
png("A:/Dropbox (Personal)/MSc Bioinformatics/Semester 2/RNA-NGS/ASSIGNMENT/data/BOX_CCHC3.png", height = 500, width = 600)
BOX_CCHC3 = ggplot(DE_CCHC_SORT_TRANS3, aes(x=sample_sheet$GROUP[c(1:28)], y=X9182, fill = sample_sheet$GROUP[c(1:28)])) + 
  geom_boxplot() + theme_Alej() + labs(title = "CCvsHC gene ENSG00000161944",x = "Group", y = "Gene expression X9182")+ coord_cartesian(ylim = c(0, 3000))
print(BOX_CCHC3)
dev.off()

#TOP3 MOST DIFFERENTIALLY EXPRESSED GENES FLU_HEALTH: X14932 X10906 X6172 	ENSG00000123119","ENSG00000204936", "ENSG00000135424"

DE_FLUHC_SORT1 <- master_file[order(master_file$p_FLU_HC, decreasing = FALSE, na.last = TRUE),]
DE_FLUHC_SORT1 =DE_FLUHC_SORT1[c(2:15,30:43)]
DE_FLUHC_SORT_TRANS1 = data.frame(t(DE_FLUHC_SORT1))
png("A:/Dropbox (Personal)/MSc Bioinformatics/Semester 2/RNA-NGS/ASSIGNMENT/data/BOX_FLUHC1.png", height = 500, width = 600)
BOX_FLUHC1 = ggplot(DE_FLUHC_SORT_TRANS1, aes(x=sample_sheet$GROUP[c(1:14,29:42)], y=X14932, fill = sample_sheet$GROUP[c(1:14,29:42)])) + 
  geom_boxplot() + theme_Alej() + labs(title = "FLUvsHC gene ENSG00000123119",x = "Group", y = "Gene expression ENSG00000123119")+coord_cartesian(ylim = c(0, 5000))
print(BOX_FLUHC1)
dev.off()

DE_FLUHC_SORT2 <- master_file[order(master_file$p_FLU_HC, decreasing = FALSE, na.last = TRUE),]
DE_FLUHC_SORT2 =DE_FLUHC_SORT2[c(2:15,30:43)]
DE_FLUHC_SORT_TRANS2 = data.frame(t(DE_FLUHC_SORT2))
png("A:/Dropbox (Personal)/MSc Bioinformatics/Semester 2/RNA-NGS/ASSIGNMENT/data/BOX_FLUHC2.png", height = 500, width = 600)
BOX_FLUHC2 = ggplot(DE_FLUHC_SORT_TRANS2, aes(x=sample_sheet$GROUP[c(1:14,29:42)], y=X10906, fill = sample_sheet$GROUP[c(1:14,29:42)])) + 
  geom_boxplot() + theme_Alej() + labs(title = "FLUvsHC gene ENSG00000204936", x = "Group", y = "Gene expression ENSG00000204936")+coord_cartesian(ylim = c(0, 600))
print(BOX_FLUHC2)
dev.off()

DE_FLUHC_SORT3 <- master_file[order(master_file$p_FLU_HC, decreasing = FALSE, na.last = TRUE),]
DE_FLUHC_SORT3 =DE_FLUHC_SORT3[c(2:15,30:43)]
DE_FLUHC_SORT_TRANS3 = data.frame(t(DE_FLUHC_SORT3))
png("A:/Dropbox (Personal)/MSc Bioinformatics/Semester 2/RNA-NGS/ASSIGNMENT/data/BOX_FLUHC3.png", height = 500, width = 600)
BOX_FLUHC3 = ggplot(DE_FLUHC_SORT_TRANS3, aes(x=sample_sheet$GROUP[c(1:14,29:42)], y=X6172, fill = sample_sheet$GROUP[c(1:14,29:42)])) + 
  geom_boxplot() + theme_Alej() + labs(title = "FLUvsHC gene ENSG00000135424",x = "Group", y = "Gene expression ENSG00000135424")+coord_cartesian(ylim = c(0, 400))
print(BOX_FLUHC3)
dev.off()

#TOP3 MOST DIFFERENTIALLY EXPRESSED GENES FLU_COLD: X4708 X14932 X6172

DE_FLUCC_SORT1 <- master_file[order(master_file$p_FLU_CC, decreasing = FALSE, na.last = TRUE),]
DE_FLUCC_SORT1 = DE_FLUCC_SORT1[c(2:43)]
DE_FLUCC_SORT_TRANS1 = data.frame(t(DE_FLUCC_SORT1)).
png("A:/Dropbox (Personal)/MSc Bioinformatics/Semester 2/RNA-NGS/ASSIGNMENT/data/BOX_FLUCC1.png", height = 500, width = 600)
BOX_FLUCC1 = ggplot(DE_FLUCC_SORT_TRANS1, aes(x=sample_sheet$GROUP, y=X4708, fill = sample_sheet$GROUP)) + 
  geom_boxplot() + theme_Alej() + labs(title = "FLUvsCC gene X4708",x = "Group", y = "Gene expression X4708")+coord_cartesian(ylim = c(0, 75))
print(BOX_FLUCC1)
dev.off()

DE_FLUCC_SORT2 <- master_file[order(master_file$p_FLU_CC, decreasing = FALSE, na.last = TRUE),]
DE_FLUCC_SORT2 = DE_FLUCC_SORT2[c(2:43)]
DE_FLUCC_SORT_TRANS2 = data.frame(t(DE_FLUCC_SORT2))
png("A:/Dropbox (Personal)/MSc Bioinformatics/Semester 2/RNA-NGS/ASSIGNMENT/data/BOX_FLUCC2.png", height = 500, width = 600)
BOX_FLUCC2 = ggplot(DE_FLUCC_SORT_TRANS2, aes(x=sample_sheet$GROUP, y=X14932, fill = sample_sheet$GROUP)) + 
  geom_boxplot() + theme_Alej() + labs(title = "FLUvsCC gene X14932",x = "Group", y = "Gene expression X14932")+coord_cartesian(ylim = c(0, 5000))
print(BOX_FLUCC2)
dev.off()

DE_FLUCC_SORT3 <- master_file[order(master_file$p_FLU_CC, decreasing = FALSE, na.last = TRUE),]
DE_FLUCC_SORT3 = DE_FLUCC_SORT3[c(2:43)]
DE_FLUCC_SORT_TRANS3 = data.frame(t(DE_FLUCC_SORT3))
png("A:/Dropbox (Personal)/MSc Bioinformatics/Semester 2/RNA-NGS/ASSIGNMENT/data/BOX_FLUCC3.png", height = 500, width = 600)
BOX_FLUCC3 = ggplot(DE_FLUCC_SORT_TRANS3, aes(x=sample_sheet$GROUP, y=X6172, fill = sample_sheet$GROUP)) + 
  geom_boxplot() + theme_Alej() + labs(title = "FLUvsCC gene X6172",x = "Group", y = "Gene expression X6172")+coord_cartesian(ylim = c(0, 400))
print(BOX_FLUCC3)
dev.off()

#VIOLIN PLUS JITTER PLOTS FOR TOP 3 MOST DIFF EXPRESSED GENES FOR EACH GROUP
#TOP 3 VIOLIN PLOTS FOR FLU_COLD X4708 X14932 X6172
#plots have not been saved to png because I chose to look at boxplots

VIOL_FLUCC1 <- ggplot(DE_FLUCC_SORT_TRANS1, aes(x=sample_sheet$GROUP, y=X4708, fill = sample_sheet$GROUP)) +  
  geom_violin(trim=FALSE) + 
  stat_summary(fun=mean, geom="point", shape=17, size=2, color="magenta")+ 
  geom_jitter(shape = 1) +
  labs(x="group",y="expression ENSG00000123119")+coord_cartesian(ylim = c(0, 1500))

VIOL_FLUCC2 <- ggplot(DE_FLUCC_SORT_TRANS2, aes(x=sample_sheet$GROUP, y=X14932, fill = sample_sheet$GROUP)) +  
  geom_violin(trim=FALSE) + 
  stat_summary(fun=mean, geom="point", shape=17, size=2, color="magenta")+ 
  geom_jitter(shape = 1) +
  labs(x="group",y="expression ENSG00000204936")

VIOL_FLUCC3 <- ggplot(DE_FLUCC_SORT_TRANS3, aes(x=sample_sheet$GROUP, y=X6172, fill = sample_sheet$GROUP)) +  
  geom_violin(trim=FALSE) + 
  stat_summary(fun=mean, geom="point", shape=17, size=2, color="magenta")+ 
  geom_jitter(shape = 1) +
  labs(x="group",y="expression ENSG00000135424")+coord_cartesian(ylim = c(0, 1500))

FLU_CC_sorted <- master_file[order(master_file$p.adj_FLU_CC, decreasing = FALSE, na.last = TRUE), ]
FLU_CC_sorted = FLU_CC_sorted[1:3,c(1:43),]
FLU_CC_sorted_scaled = as.matrix(sapply(FLU_CC_sorted[,c(15:42)], as.numeric))
FLU_CC_sorted_scaled = scale(FLU_CC_sorted_scaled, center = TRUE, scale = TRUE)
FLU_CC_sorted_scaled = data.frame(FLU_CC_sorted_scaled)
FLU_CC_sorted_scaled$ID = c("	ENSG00000123119","ENSG00000204936", "ENSG00000135424")
FLU_CC_sorted_melted = melt(FLU_CC_sorted_scaled)
FLU_CC_sorted_melted$variable = sample_sheet$GROUP [match(unlist(FLU_CC_sorted_melted$variable), sample_sheet$SAMPLE)]

All_FLUCC_violin = ggplot(FLU_CC_sorted_melted, aes(x=ID, y=value, fill = variable)) + 
  geom_violin(alpha = 0.4, trim = F)+ 
  theme_classic()+ 
  labs(x="group",y="expression of ") + 
  coord_cartesian(ylim = c(-1.5, 1.5))+
  geom_jitter(shape=1)

#TOP 3 VIOLIN PLOTS FOR FLU_HEALTH X14932 X10906 X6172

VIOL_FLUHC1 <- ggplot(DE_FLUHC_SORT_TRANS1, aes(x=sample_sheet$GROUP[c(1:14,29:42)], y=X14932, fill = sample_sheet$GROUP)) +  
  geom_violin(trim=FALSE) + 
  stat_summary(fun=mean, geom="point", shape=17, size=2, color="magenta")+ 
  geom_jitter(shape = 1) +
  labs(x="group",y="expression ENSG00000204936")

VIOL_FLUHC2 <- ggplot(DE_FLUHC_SORT_TRANS2, aes(x=sample_sheet$GROUP[c(1:14,29:42)], y=X10906, fill = sample_sheet$GROUP)) +  
  geom_violin(trim=FALSE) + 
  stat_summary(fun=mean, geom="point", shape=17, size=2, color="magenta")+ 
  geom_jitter(shape = 1) +
  labs(x="group",y="expression ENSG00000170439")+coord_cartesian(ylim = c(0, 1000))

VIOL_FLUHC3 <- ggplot(DE_FLUHC_SORT_TRANS3, aes(x=sample_sheet$GROUP[c(1:14,29:42)], y=X6172, fill = sample_sheet$GROUP)) +  
  geom_violin(trim=FALSE) + 
  stat_summary(fun=mean, geom="point", shape=17, size=2, color="magenta")+ 
  geom_jitter(shape = 1) +
  labs(x="group",y="expression ENSG00000135424")+coord_cartesian(ylim = c(0, 1500))

FLU_HC_sorted <- master_file[order(master_file$p.adj_FLU_HC, decreasing = FALSE, na.last = TRUE), ]
FLU_HC_sorted = FLU_HC_sorted[1:3,c(1:43),]
FLU_HC_sorted_scaled = as.matrix(sapply(FLU_HC_sorted[,c(1:14,29:42)], as.numeric))
FLU_HC_sorted_scaled = scale(FLU_HC_sorted_scaled, center = TRUE, scale = TRUE)
FLU_HC_sorted_scaled = data.frame(FLU_HC_sorted_scaled)
FLU_HC_sorted_scaled$ID = c("	ENSG00000204936","ENSG00000170439", "ENSG00000135424")
FLU_HC_sorted_melted = melt(FLU_HC_sorted_scaled)
FLU_HC_sorted_melted$variable = sample_sheet$GROUP [match(unlist(FLU_HC_sorted_melted$variable), sample_sheet$SAMPLE)]

All_FLUHC_violin = ggplot(FLU_HC_sorted_melted, aes(x=ID, y=value, fill = variable)) + 
  geom_violin(alpha = 0.4, trim = F)+ 
  theme_classic()+ 
  labs(x="group",y="expression of ") + 
  coord_cartesian(ylim = c(-1.5, 1.5))+
  geom_jitter(shape=1)

#TOP 3 VIOLIN PLOTS FOR COLD_HEALTH  X3210 X10906 X9182

VIOL_CCHC1 <- ggplot(DE_CCHC_SORT_TRANS1, aes(x=sample_sheet$GROUP[c(1:28)], y=X3210, fill = sample_sheet$GROUP)) +  
  geom_violin(trim=FALSE) + 
  stat_summary(fun=mean, geom="point", shape=17, size=2, color="magenta")+ 
  geom_jitter(shape = 1) +
  labs(x="group",y="expression ENSG00000108950")

VIOL_CCHC2 <- ggplot(DE_CCHC_SORT_TRANS2, aes(x=sample_sheet$GROUP[c(1:28)], y=X10906, fill = sample_sheet$GROUP)) +  
  geom_violin(trim=FALSE) + 
  stat_summary(fun=mean, geom="point", shape=17, size=2, color="magenta")+ 
  geom_jitter(shape = 1) +
  labs(x="group",y="expression ENSG00000170439")

VIOL_CCHC3 <- ggplot(DE_CCHC_SORT_TRANS3, aes(x=sample_sheet$GROUP[c(1:28)], y=X9182, fill = sample_sheet$GROUP)) +  
  geom_violin(trim=FALSE) + 
  stat_summary(fun=mean, geom="point", shape=17, size=2, color="magenta")+ 
  geom_jitter(shape = 1) +
  labs(x="group",y="expression ENSG00000161944")

CC_HC_sorted <- master_file[order(master_file$p.adj_CC_HC, decreasing = FALSE, na.last = TRUE), ]
CC_HC_sorted = CC_HC_sorted[1:3,c(1:43),]
CC_HC_sorted_scaled = as.matrix(sapply(CC_HC_sorted[,c(1:28)], as.numeric))
CC_HC_sorted_scaled = scale(CC_HC_sorted_scaled, center = TRUE, scale = TRUE)
CC_HC_sorted_scaled = data.frame(CC_HC_sorted_scaled)
CC_HC_sorted_scaled$ID = c("	ENSG00000108950","ENSG00000170439", "ENSG00000161944")
CC_HC_sorted_melted = melt(CC_HC_sorted_scaled)
CC_HC_sorted_melted$variable = sample_sheet$GROUP [match(unlist(CC_HC_sorted_melted$variable), sample_sheet$SAMPLE)]

All_CCHC_violin = ggplot(CC_HC_sorted_melted, aes(x=ID, y=value, fill = variable)) + 
  geom_box(alpha = 0.4)+ 
  theme_classic()+ 
  labs(x="group",y="expression of ") + 
  coord_cartesian(ylim = c(-1.5, 1.5))+
  geom_jitter(shape=1)

#PCA plot for PC1 vs PC2 and PC3 vs PC4
numeric_data <- as.matrix(sapply(Expression_matrix[2:43],as.numeric))
#because prcomp uses samples as row and genes as columns and we have the opposite: transpose matrix
pca = prcomp(t(numeric_data))
#to extract the component data, matrix called x, we apply the following 
pca_coordinates=data.frame(pca$x)
#if we don't transpose: graph that shows how genes are related to each other! Instead of how samples are or not related to each other
#prcomp returns x: principal components for drawing a graph. X samples, X PCs. First PCs show the most variation across data (gene expression across samples)

#autoplot generates the PCA plots without customisation - proper PCA plots with geoms are done for z-scores
PC1_PC2 = autoplot(prcomp(t(numeric_data)), data = sample_sheet, colour = "GROUP", label = TRUE, label.size = 2)
PC3_PC4 = autoplot(prcomp(t(numeric_data)), data = sample_sheet,x= 3, y=4, colour = "GROUP", label = TRUE, label.size = 2)

#with Z-scores
zNumeric <- as.matrix(sapply(Expression_matrix[2:43], as.numeric))
zTrans = t(zNumeric)
zScale = scale(zNumeric, center = TRUE, scale = TRUE)
zPrcomp = prcomp(na.omit(t(zScale)))
zCoor = data.frame(zPrcomp$x)

sdev <- zPrcomp$sdev #standard deviation of each principal component
variance <- sdev^2 #variance explained by each principal component
proportion_var <- variance / sum(variance) #compute the proportion of variance explained by each principal component
proportion <- round(proportion_var*100,2) #round numbers (aesthetically better looking)
percentage <- paste( colnames(zCoor), "(", paste( as.character(proportion), "%", ")", sep="") )

png("A:/Dropbox (Personal)/MSc Bioinformatics/Semester 2/RNA-NGS/ASSIGNMENT/data/PC1_PC2_Z.png", height = 500, width = 600)
PC1_PC2_Z=ggplot(zCoor,aes(x=PC1,y=PC2, color = sample_sheet$GROUP, label = row.names(zCoor))) +
  geom_point() + 
  geom_text(size = 2)+
  xlab(percentage[1]) + ylab(percentage[2]) + labs(title = "PC1 vs PC2")+
  theme_Alej()
print(PC1_PC2_Z)
dev.off()

PC3_PC4_Z=ggplot(zCoor,aes(x=PC3,y=PC4, color = sample_sheet$GROUP, label = row.names(zCoor))) +
  geom_point()+ 
  geom_text(size = 2)+
  xlab(percentage[3]) + ylab(percentage[4]) + labs(title = "PC3 vs PC4")+
  theme_Alej()

#Eigenvalues plot for PCA (shows the percentage of variances explained by each principal component)
png("A:/Dropbox (Personal)/MSc Bioinformatics/Semester 2/RNA-NGS/ASSIGNMENT/data/fviz_eig.png", height = 500, width = 600)
fviz_eig(zPrcomp, barfill = "light pink", barcolor = "light pink") + theme_Alej()
print(fviz_eig)
dev.off()

#CORRELATION HEATMAPS - ZSCORES

colours = c("magenta","white","dark blue")

Na_omit <- na.omit(zScale)
cors = cor(Na_omit, method="spearman")
#correlation function correlates columns across rows; samples by col and genes by row
melted=melt(cors)
#parsing of the heatmap with geom_tile - to generate the three column matrix with x and y samples and the correlation calc: melt

png("A:/Dropbox (Personal)/MSc Bioinformatics/Semester 2/RNA-NGS/ASSIGNMENT/data/heatmap_z.png", height = 500, width = 600)
heatmap_z = ggplot(melted, aes(x = X1, y = X2, fill = value)) + 
  geom_tile()+scale_fill_gradientn(limits=c(0.80, 1),colours =colorRampPalette(colours)(150))+
  theme(axis.text.x=element_text(angle = 90,hjust=1,vjust=0.5)) +
  labs(title = "Correlation heatmap")
print(heatmap_z)
dev.off()

#Heatmaps generated by pheatmap function show a lesser level of customization and have not been included in report
#pheatmap_Z = pheatmap(cors) 

#Expression heatmap for each of the three differential comparisons

coloursExp = c("dark green", "white", "magenta")

#COLD VS HEALTH CONTROL
SIG_CCHC = master_file[which(master_file$sig_CC_HC == "TRUE"),]
SIG_CCHC_TABLE = SIG_CCHC[c(2:29)]
row.names(SIG_CCHC_TABLE) = SIG_CCHC$ID
SIG_CCHC_ZTR = na.omit(t(SIG_CCHC_TABLE))
SIG_CCHC_ZSCALE = scale(SIG_CCHC_ZTR, center = TRUE, scale = TRUE)
SIG_CCHC_ZSCORES = t(SIG_CCHC_ZSCALE)

CCHC_dist = Dist(SIG_CCHC_ZSCORES, method = "spearman")
CCHC_cluster = hclust(CCHC_dist, method = "average")
CCHC_dd = as.dendrogram(CCHC_cluster)
CCHC_ddRe = reorder(CCHC_dd,0, FUN = "average")
CCHC_order = order.dendrogram(CCHC_ddRe)
CCHC_orderlab = data.frame(SIG_CCHC_ZSCORES[CCHC_order,])
CCHC_orderlab$ID = rownames(CCHC_orderlab)

SIG_CCHC_melt = melt(CCHC_orderlab, na.rm = T)

png("A:/Dropbox (Personal)/MSc Bioinformatics/Semester 2/RNA-NGS/ASSIGNMENT/data/Exp_heatmap_CCHC.png", height = 500, width = 600)
Exp_heatmap_CCHC = ggplot(SIG_CCHC_melt, aes(x=variable, y=factor(ID, level=CCHC_orderlab$ID), fill=value)) + 
  geom_tile() + scale_fill_gradientn(colours =colorRampPalette(coloursExp)(100))+ 
  theme(axis.text.x=element_text(angle = 90,hjust=1,vjust=0.5),axis.text.y = element_blank(), axis.ticks.y = element_blank())+
  labs(title = "Expression heatmap CC vs HC")
print(Exp_heatmap_CCHC)
dev.off()

# FLU VS HEALTH CONTROL
SIG_FLUHC = master_file[which(master_file$sig_FLU_HC == "TRUE"),]
SIG_FLUHC_TABLE = SIG_FLUHC[c(2:15,30:43)]
row.names(SIG_FLUHC_TABLE) = SIG_FLUHC$ID
SIG_FLUHC_ZTR = na.omit(t(SIG_FLUHC_TABLE))
SIG_FLUHC_ZSCALE = scale(SIG_FLUHC_ZTR, center = TRUE, scale = TRUE)
SIG_FLUHC_ZSCORES = t(SIG_FLUHC_ZSCALE)

FLUHC_dist = Dist(SIG_FLUHC_ZSCORES, method = "spearman")
FLUHC_cluster = hclust(FLUHC_dist, method = "average")
FLUHC_dd = as.dendrogram(FLUHC_cluster)
FLUHC_ddRe = reorder(FLUHC_dd,0, FUN = "average")
FLUHC_order = order.dendrogram(FLUHC_ddRe)
FLUHC_orderlab = data.frame(SIG_FLUHC_ZSCORES[FLUHC_order,])
FLUHC_orderlab$ID = rownames(FLUHC_orderlab)

SIG_FLUHC_melt = melt(FLUHC_orderlab, na.rm = T)

png("A:/Dropbox (Personal)/MSc Bioinformatics/Semester 2/RNA-NGS/ASSIGNMENT/data/Exp_heatmap_FLUHC.png", height = 500, width = 600)
Exp_heatmap_FLUHC = ggplot(SIG_FLUHC_melt, aes(x=variable, y=factor(ID, level=FLUHC_orderlab$ID), fill=value)) + geom_tile() + 
  scale_fill_gradientn(colours =colorRampPalette(coloursExp)(100))+ 
  theme(axis.text.x=element_text(angle = 90,hjust=1,vjust=0.5), axis.text.y = element_blank(), axis.ticks.y = element_blank())+
  labs(title = "Expression heatmap FLU vs HC")
print(Exp_heatmap_FLUHC)
dev.off()

# FLU VS COLD 
SIG_FLUCC = master_file[which(master_file$sig_FLU_CC == "TRUE"),]
SIG_FLUCC_TABLE = SIG_FLUCC[c(16:43)]
row.names(SIG_FLUCC_TABLE) = SIG_FLUCC$ID
SIG_FLUCC_ZTR = na.omit(t(SIG_FLUCC_TABLE))
SIG_FLUCC_ZSCALE = scale(SIG_FLUCC_ZTR, center = TRUE, scale = TRUE)
SIG_FLUCC_ZSCORES = t(SIG_FLUCC_ZSCALE)

FLUCC_dist = Dist(SIG_FLUCC_ZSCORES, method = "spearman")
FLUCC_cluster = hclust(FLUCC_dist, method = "average")
FLUCC_dd = as.dendrogram(FLUCC_cluster)
FLUCC_ddRe = reorder(FLUCC_dd,0, FUN = "average")
FLUCC_order = order.dendrogram(FLUCC_ddRe)
FLUCC_orderlab = data.frame(SIG_FLUCC_ZSCORES[FLUCC_order,])
FLUCC_orderlab$ID = rownames(FLUCC_orderlab)

SIG_FLUCC_melt = melt(FLUCC_orderlab, na.rm = T)

png("A:/Dropbox (Personal)/MSc Bioinformatics/Semester 2/RNA-NGS/ASSIGNMENT/data/Exp_heatmap_FLUHC.png", height = 500, width = 600)
Exp_heatmap_FLUCC = ggplot(SIG_FLUCC_melt, aes(x=variable, y=factor(ID, level=FLUCC_orderlab$ID), fill=value)) + 
  geom_tile() + scale_fill_gradientn(colours =colorRampPalette(coloursExp)(100))+ 
  theme(axis.text.x=element_text(angle = 90,hjust=1,vjust=0.5), axis.text.y = element_blank(), axis.ticks.y = element_blank())+
  labs(title = "Expression heatmap FLU vs CC")
print(Exp_heatmap_FLUCC)
dev.off()

#RNA-seq [AvsB]n Plots

#Expression Heatmap combining all three differential comparisons and all samples 
#Gene (y-axis) clustering
colours2 = c("dark blue", "white", "magenta")
SIG_table1 = master_file[which(master_file$sig_FLU_CC == "TRUE" | master_file$sig_FLU_HC == "TRUE" | master_file$sig_CC_HC == "TRUE" ),]

SIG_All_TABLE = SIG_table1[c(2:43)]
row.names(SIG_All_TABLE) = SIG_table1$ID
SIG_All_ZTR = na.omit(t(SIG_All_TABLE))
SIG_All_ZSCALE = scale(SIG_All_ZTR, center = TRUE, scale = TRUE)
SIG_All_ZSCORES = t(SIG_All_ZSCALE)

#pheatmap(SIG_All_ZSCORES, cluster_rows = T, cluster_cols = F)

All_dist = Dist(SIG_All_ZSCORES, method = "spearman")
All_cluster = hclust(All_dist, method = "average")
All_dd = as.dendrogram(All_cluster)
All_ddRe = reorder(All_dd,0, FUN = "average")
All_order = order.dendrogram(All_ddRe)
All_orderlab = data.frame(SIG_All_ZSCORES[All_order,])
All_orderlab$ID = rownames(All_orderlab)

SIG_All_melt = melt(All_orderlab, na.rm = T)

png("A:/Dropbox (Personal)/MSc Bioinformatics/Semester 2/RNA-NGS/ASSIGNMENT/data/Exp_heatmap_All.png", height = 500, width = 600)
Exp_heatmap_All = ggplot(SIG_All_melt, aes(x=variable, y=factor(ID, level=All_orderlab$ID), fill=value)) + 
  geom_tile() + 
  theme(axis.text.x=element_text(angle = 90,hjust=1,vjust=0.5), axis.text.y = element_blank(), axis.ticks.y = element_blank())+
  scale_fill_gradientn(colours =colorRampPalette(colours2)(100))+
  labs(title = "Expression heatmap of the three differential comparisons")
print(Exp_heatmap_All)
dev.off()

#Fold vs Fold Scatter-plot

#CCvsHC & FLUvsCC
png("E:/RNA-NGS/ASSIGNMENT/data/CCHC_FLUCC_scatter.png", height = 500, width = 600)
master_file$sig_All <- as.factor(master_file$p.adj_CC_HC < 0.05 & abs(master_file$log2fold_CC_HC) > 1.0 & master_file$p.adj_FLU_CC < 0.05 & abs(master_file$log2fold_FLU_CC) > 1.0) 
master_file$sig_None <- as.factor((master_file$p.adj_CC_HC > 0.05 | abs(master_file$log2fold_CC_HC) < 1.0) & (master_file$p.adj_FLU_CC > 0.05 | abs(master_file$log2fold_FLU_CC) < 1.0)) 
master_file$sig_CCvsHC <- as.factor((master_file$p.adj_CC_HC < 0.05 & abs(master_file$log2fold_CC_HC) > 1.0) & (master_file$p.adj_FLU_CC > 0.05 | abs(master_file$log2fold_FLU_CC) < 1.0)) 
master_file$sig_FLUvsCC  <- as.factor((master_file$p.adj_CC_HC > 0.05 | abs(master_file$log2fold_CC_HC) < 1.0) & (master_file$p.adj_FLU_CC < 0.05 & abs(master_file$log2fold_FLU_CC) > 1.0)) 

master_file=na.omit(master_file)
Significance_colour = case_when(master_file$sig_FLUvsCC=="TRUE" ~ "sig_only_FLUvsCC",
                      master_file$sig_CCvsHC=="TRUE" ~ "sig_only_CCvsHC",
                      master_file$sig_All=="TRUE" ~ "sig_both",
                      master_file$sig_None=="TRUE" ~ "not_sig")


master_file$sig_type
cor(master_file$log2fold_CC_HC,y= master_file$log2fold_FLU_CC, method = "spearman")
CCHC_FLUCC_scatter = ggplot(master_file, aes(x = log2fold_CC_HC, y = log2fold_FLU_CC, color=Significance_colour)) + 
  geom_point() + 
  stat_cor(method = "spearman")+ 
  labs(x="CC vs HC log2fold", y="FLU vs CC log2fold", title = "CCHC vs FLUCC scatter plot")+
  theme_Alej()
print(CCHC_FLUCC_scatter)
dev.off()

#CCvsHC & FLUvsHC
png("E:/RNA-NGS/ASSIGNMENT/data/CCHC_FLUHC_scatter.png", height = 500, width = 600)
master_file$sig_All2 <- as.factor(master_file$p.adj_CC_HC < 0.05 & abs(master_file$log2fold_CC_HC) > 1.0 & master_file$p.adj_FLU_HC < 0.05 & abs(master_file$log2fold_FLU_HC) > 1.0) 
master_file$sig_None2 <- as.factor((master_file$p.adj_CC_HC > 0.05 | abs(master_file$log2fold_CC_HC) < 1.0) & (master_file$p.adj_FLU_HC > 0.05 | abs(master_file$log2fold_FLU_HC) < 1.0)) 
master_file$sig_CCvsHC2 <- as.factor((master_file$p.adj_CC_HC < 0.05 & abs(master_file$log2fold_CC_HC) > 1.0) & (master_file$p.adj_FLU_HC > 0.05 | abs(master_file$log2fold_FLU_HC) < 1.0)) 
master_file$sig_FLUvsHC2  <- as.factor((master_file$p.adj_CC_HC > 0.05 | abs(master_file$log2fold_CC_HC) < 1.0) & (master_file$p.adj_FLU_HC < 0.05 & abs(master_file$log2fold_FLU_HC) > 1.0)) 

master_file=na.omit(master_file)
Significance_colour2 = case_when(master_file$sig_FLUvsHC=="TRUE" ~ "sig_only_FLUvsHC",
                                master_file$sig_CCvsHC=="TRUE" ~ "sig_only_CCvsHC",
                                master_file$sig_All=="TRUE" ~ "sig_both",
                                master_file$sig_None=="TRUE" ~ "not_sig")


master_file$sig_type
cor(master_file$log2fold_CC_HC,y= master_file$log2fold_FLU_HC, method = "spearman")
CCHC_FLUHC_scatter = ggplot(master_file, aes(x = log2fold_CC_HC, y = log2fold_FLU_HC)) + 
  geom_point() + 
  stat_cor(method = "spearman") + 
  labs(title = "CCHC vs FLUHC scatter plot")+
  theme_Alej()
print(CCHC_FLUHC_scatter)
dev.off()

#FLUvsHC & FLUvsCC
FLUHC_FLUCC_scatter = ggplot(master_file, aes(x = log2fold_FLU_HC, y = log2fold_FLU_CC)) + 
  geom_point() + 
  stat_cor(method = "spearman", label.x = 3, label.y = 30) + 
  labs(title = "FLUCC vs FLUHC scatter plot")

#Expression matrices for each differential expression profile: expression heatmap for each of the gene lists

# CCHC_US_FLUCC_DS, CCHC_DS_FLUCC_DS, CCHC_US_FLUHC_DS, CCHC_DS_FLUHC_US, FLUCC_US_FLUHC_DS, FLUCC_DS_FLUHC_US, have 0 variables

#FLU VS COLD BOTH UPREGULATED CCHC_US_FLUCC_US heatmap
SIG_CCHC_US_FLUCC_US_TABLE = CCHC_US_FLUCC_US[c(2:43)]
row.names(SIG_CCHC_US_FLUCC_US_TABLE) = CCHC_US_FLUCC_US$ID
NC_ZTR_CCHC_US_FLUCC_US = na.omit(t(SIG_CCHC_US_FLUCC_US_TABLE))
NC_ZSCALE_CCHC_US_FLUCC_US = scale(SIG_ZTR_CCHC_US_FLUCC_US, center = TRUE, scale = TRUE)
NC_ZSCORES_CCHC_US_FLUCC_US = t(SIG_ZSCALE_CCHC_US_FLUCC_US)

CCHC_UP_FLUCC_UP_HEATMAP = pheatmap(NC_ZSCORES_CCHC_US_FLUCC_US, cluster_rows = T, cluster_cols = F)

# FLU VS COLD UP REGULATED AND NO CHANGE CCHC_US_FLUCC_NC
SIG_CCHC_US_FLUCC_NC_TABLE = CCHC_US_FLUCC_NC[c(2:43)]
row.names(SIG_CCHC_US_FLUCC_NC_TABLE) = CCHC_US_FLUCC_NC$ID
NC_ZTR_CCHC_US_FLUCC_NC = na.omit(t(SIG_CCHC_US_FLUCC_NC_TABLE))
NC_ZSCALE_CCHC_US_FLUCC_NC = scale(NC_ZTR_CCHC_US_FLUCC_NC, center = TRUE, scale = TRUE)
NC_ZSCORES_CCHC_US_FLUCC_NC = t(NC_ZSCALE_CCHC_US_FLUCC_NC)

CCHC_UP_FLUCC_NC_HEATMAP = pheatmap(NC_ZSCORES_CCHC_US_FLUCC_NC, cluster_rows = T, cluster_cols = F)

# FLU UP COLD DOWN REGULATED CCHC_DS_FLUCC_US

SIG_CCHC_DS_FLUCC_US_TABLE = CCHC_DS_FLUCC_US[c(2:43)]
row.names(SIG_CCHC_DS_FLUCC_US_TABLE) = CCHC_DS_FLUCC_US$ID
NC_ZTR_CCHC_DS_FLUCC_US = na.omit(t(SIG_CCHC_DS_FLUCC_US_TABLE))
NC_ZSCALE_CCHC_DS_FLUCC_US = scale(NC_ZTR_CCHC_DS_FLUCC_US, center = TRUE, scale = TRUE)
NC_ZSCORES_CCHC_DS_FLUCC_US = t(NC_ZSCALE_CCHC_DS_FLUCC_US)

CCHC_DN_FLUCC_UP_HEATMAP = pheatmap(NC_ZSCORES_CCHC_DS_FLUCC_US, cluster_rows = T, cluster_cols = F)

# FLU NO CHANGE COLD DOWN REGULTATED CCHC_DS_FLUCC_NC

SIG_CCHC_DS_FLUCC_NC_TABLE = CCHC_DS_FLUCC_NC[c(2:43)]
row.names(SIG_CCHC_DS_FLUCC_NC_TABLE) = CCHC_DS_FLUCC_NC$ID
NC_ZTR_CCHC_DS_FLUCC_NC = na.omit(t(SIG_CCHC_DS_FLUCC_NC_TABLE))
NC_ZSCALE_CCHC_DS_FLUCC_NC = scale(NC_ZTR_CCHC_DS_FLUCC_NC, center = TRUE, scale = TRUE)
NC_ZSCORES_CCHC_DS_FLUCC_NC = t(NC_ZSCALE_CCHC_DS_FLUCC_NC)

CCHC_DN_FLUCC_NC_HEATMAP = pheatmap(NC_ZSCORES_CCHC_DS_FLUCC_NC, cluster_rows = T, cluster_cols = F)

# COLD VS FLU NO CHANGE CCHC_NC_FLUCC_NC

SIG_CCHC_NC_FLUCC_NC_TABLE = CCHC_NC_FLUCC_NC[c(2:43)]
row.names(SIG_CCHC_NC_FLUCC_NC_TABLE) = CCHC_NC_FLUCC_NC$ID
NC_ZTR_CCHC_NC_FLUCC_NC = na.omit(t(SIG_CCHC_NC_FLUCC_NC_TABLE))
NC_ZSCALE_CCHC_NC_FLUCC_NC = scale(NC_ZTR_CCHC_NC_FLUCC_NC, center = TRUE, scale = TRUE)
NC_ZSCORES_CCHC_NC_FLUCC_NC = t(NC_ZSCALE_CCHC_NC_FLUCC_NC)

CCHC_NC_FLUCC_NC_HEATMAP = pheatmap(NC_ZSCORES_CCHC_NC_FLUCC_NC, cluster_rows = T, cluster_cols = F)

# COLD NO CHANGE FLU DOWN CCHC_NC_FLUCC_DS

SIG_CCHC_NC_FLUCC_DS_TABLE = CCHC_NC_FLUCC_DS[c(2:43)]
row.names(SIG_CCHC_NC_FLUCC_DS_TABLE) = CCHC_NC_FLUCC_DS$ID
NC_ZTR_CCHC_NC_FLUCC_DS = na.omit(t(SIG_CCHC_NC_FLUCC_DS_TABLE))
NC_ZSCALE_CCHC_NC_FLUCC_DS = scale(NC_ZTR_CCHC_NC_FLUCC_DS, center = TRUE, scale = TRUE)
NC_ZSCORES_CCHC_NC_FLUCC_DS = t(NC_ZSCALE_CCHC_NC_FLUCC_DS)

CCHC_NC_FLUCC_DN_HEATMAP = pheatmap(NC_ZSCORES_CCHC_NC_FLUCC_DS, cluster_rows = T, cluster_cols = F)

# COLD NO CHANGE FLU UP CCHC_NC_FLUCC_US

SIG_CCHC_NC_FLUCC_US_TABLE = CCHC_NC_FLUCC_US[c(2:43)]
row.names(SIG_CCHC_NC_FLUCC_US_TABLE) = CCHC_NC_FLUCC_US$ID
NC_ZTR_CCHC_NC_FLUCC_US = na.omit(t(SIG_CCHC_NC_FLUCC_US_TABLE))
NC_ZSCALE_CCHC_NC_FLUCC_US = scale(NC_ZTR_CCHC_NC_FLUCC_US, center = TRUE, scale = TRUE)
NC_ZSCORES_CCHC_NC_FLUCC_US = t(NC_ZSCALE_CCHC_NC_FLUCC_US)

CCHC_NC_FLUCC_UP_HEATMAP = pheatmap(NC_ZSCORES_CCHC_NC_FLUCC_US, cluster_rows = T, cluster_cols = F)

# COLD UP FLU HEALTH UP CCHC_US_FLUHC_US

SIG_CCHC_US_FLUHC_US_TABLE = CCHC_US_FLUHC_US[c(2:43)]
row.names(SIG_CCHC_US_FLUHC_US_TABLE) = CCHC_US_FLUHC_US$ID
NC_ZTR_CCHC_US_FLUHC_US = na.omit(t(SIG_CCHC_US_FLUHC_US_TABLE))
NC_ZSCALE_CCHC_US_FLUHC_US = scale(NC_ZTR_CCHC_US_FLUHC_US, center = TRUE, scale = TRUE)
NC_ZSCORES_CCHC_US_FLUHC_US = t(NC_ZSCALE_CCHC_US_FLUHC_US)

CCHC_UP_FLUHC_UP_HEATMAP = pheatmap(NC_ZSCORES_CCHC_US_FLUHC_US, cluster_rows = T, cluster_cols = F)

# COLD UP FLU HEALTH NO CHANGE CCHC_US_FLUHC_NC

SIG_CCHC_US_FLUHC_NC_TABLE = CCHC_US_FLUHC_NC[c(2:43)]
row.names(SIG_CCHC_US_FLUHC_NC_TABLE) = CCHC_US_FLUHC_NC$ID
NC_ZTR_CCHC_US_FLUHC_NC = na.omit(t(SIG_CCHC_US_FLUHC_NC_TABLE))
NC_ZSCALE_CCHC_US_FLUHC_NC = scale(NC_ZTR_CCHC_US_FLUHC_NC, center = TRUE, scale = TRUE)
NC_ZSCORES_CCHC_US_FLUHC_NC = t(NC_ZSCALE_CCHC_US_FLUHC_NC)

CCHC_UP_FLUHC_NC_HEATMAP = pheatmap(NC_ZSCORES_CCHC_US_FLUHC_NC, cluster_rows = T, cluster_cols = F)

# COLD DOWN FLU HEALTH DOWN CCHC_DS_FLUHC_DS

SIG_CCHC_DS_FLUHC_DS_TABLE = CCHC_DS_FLUHC_DS[c(2:43)]
row.names(SIG_CCHC_DS_FLUHC_DS_TABLE) = CCHC_DS_FLUHC_DS$ID
NC_ZTR_CCHC_DS_FLUHC_DS = na.omit(t(SIG_CCHC_DS_FLUHC_DS_TABLE))
NC_ZSCALE_CCHC_DS_FLUHC_DS = scale(NC_ZTR_CCHC_DS_FLUHC_DS, center = TRUE, scale = TRUE)
NC_ZSCORES_CCHC_DS_FLUHC_DS = t(NC_ZSCALE_CCHC_DS_FLUHC_DS)

CCHC_DN_FLUHC_DN_HEATMAP = pheatmap(NC_ZSCORES_CCHC_DS_FLUHC_DS, cluster_rows = T, cluster_cols = F)

# COLD DOWN FLU HEALTH NO CHANGE CCHC_DS_FLUHC_NC

SIG_CCHC_DS_FLUHC_NC_TABLE = CCHC_DS_FLUHC_NC[c(2:43)]
row.names(SIG_CCHC_DS_FLUHC_NC_TABLE) = CCHC_DS_FLUHC_NC$ID
NC_ZTR_CCHC_DS_FLUHC_NC = na.omit(t(SIG_CCHC_DS_FLUHC_NC_TABLE))
NC_ZSCALE_CCHC_DS_FLUHC_NC = scale(NC_ZTR_CCHC_DS_FLUHC_NC, center = TRUE, scale = TRUE)
NC_ZSCORES_CCHC_DS_FLUHC_NC = t(NC_ZSCALE_CCHC_DS_FLUHC_NC)

CCHC_DN_FLUHC_NC_HEATMAP = pheatmap(NC_ZSCORES_CCHC_DS_FLUHC_NC, cluster_rows = T, cluster_cols = F)

# COLD NO CHANGE FLU NO CHANGE CCHC_NC_FLUHC_NC

SIG_CCHC_NC_FLUHC_NC_TABLE = CCHC_NC_FLUHC_NC[c(2:43)]
row.names(SIG_CCHC_NC_FLUHC_NC_TABLE) = CCHC_NC_FLUHC_NC$ID
NC_ZTR_CCHC_NC_FLUHC_NC = na.omit(t(SIG_CCHC_NC_FLUHC_NC_TABLE))
NC_ZSCALE_CCHC_NC_FLUHC_NC = scale(NC_ZTR_CCHC_NC_FLUHC_NC, center = TRUE, scale = TRUE)
NC_ZSCORES_CCHC_NC_FLUHC_NC = t(NC_ZSCALE_CCHC_NC_FLUHC_NC)

CCHC_NC_FLUHC_NC_HEATMAP = pheatmap(NC_ZSCORES_CCHC_NC_FLUHC_NC, cluster_rows = T, cluster_cols = F)

# COLD NO CHANGE FLU DOWN CCHC_NC_FLUHC_DS

SIG_CCHC_NC_FLUHC_DS_TABLE = CCHC_NC_FLUHC_DS[c(2:43)]
row.names(SIG_CCHC_NC_FLUHC_DS_TABLE) = CCHC_NC_FLUHC_DS$ID
NC_ZTR_CCHC_NC_FLUHC_DS = na.omit(t(SIG_CCHC_NC_FLUHC_DS_TABLE))
NC_ZSCALE_CCHC_NC_FLUHC_DS = scale(NC_ZTR_CCHC_NC_FLUHC_DS, center = TRUE, scale = TRUE)
NC_ZSCORES_CCHC_NC_FLUHC_DS = t(NC_ZSCALE_CCHC_NC_FLUHC_DS)

CCHC_NC_FLUHC_DN_HEATMAP = pheatmap(NC_ZSCORES_CCHC_NC_FLUHC_DS, cluster_rows = T, cluster_cols = F)

# COLD NO CHANGE FLU UP CCHC_NC_FLUHC_US

SIG_CCHC_NC_FLUHC_US_TABLE = CCHC_NC_FLUHC_US[c(2:43)]
row.names(SIG_CCHC_NC_FLUHC_US_TABLE) = CCHC_NC_FLUHC_US$ID
NC_ZTR_CCHC_NC_FLUHC_US = na.omit(t(SIG_CCHC_NC_FLUHC_US_TABLE))
NC_ZSCALE_CCHC_NC_FLUHC_US = scale(NC_ZTR_CCHC_NC_FLUHC_US, center = TRUE, scale = TRUE)
NC_ZSCORES_CCHC_NC_FLUHC_US = t(NC_ZSCALE_CCHC_NC_FLUHC_US)

CCHC_NC_FLUHC_UP_HEATMAP = pheatmap(NC_ZSCORES_CCHC_NC_FLUHC_US, cluster_rows = T, cluster_cols = F)

# FLU-COLD UP FLU HEALTH UP FLUCC_US_FLUHC_US

SIG_FLUCC_US_FLUHC_US_TABLE = FLUCC_US_FLUHC_US[c(2:43)]
row.names(SIG_FLUCC_US_FLUHC_US_TABLE) = FLUCC_US_FLUHC_US$ID
NC_ZTR_FLUCC_US_FLUHC_US = na.omit(t(SIG_FLUCC_US_FLUHC_US_TABLE))
NC_ZSCALE_FLUCC_US_FLUHC_US = scale(NC_ZTR_FLUCC_US_FLUHC_US, center = TRUE, scale = TRUE)
NC_ZSCORES_FLUCC_US_FLUHC_US = t(NC_ZSCALE_FLUCC_US_FLUHC_US)

FLUCC_UP_FLUHC_UP_HEATMAP = pheatmap(NC_ZSCORES_FLUCC_US_FLUHC_US, cluster_rows = T, cluster_cols = F)

# FLU-COLD UP FLU HEALTH NO CHANGE FLUCC_US_FLUHC_NC

SIG_FLUCC_US_FLUHC_NC_TABLE = FLUCC_US_FLUHC_NC[c(2:43)]
row.names(SIG_FLUCC_US_FLUHC_NC_TABLE) = FLUCC_US_FLUHC_NC$ID
NC_ZTR_FLUCC_US_FLUHC_NC = na.omit(t(SIG_FLUCC_US_FLUHC_NC_TABLE))
NC_ZSCALE_FLUCC_US_FLUHC_NC = scale(NC_ZTR_FLUCC_US_FLUHC_NC, center = TRUE, scale = TRUE)
NC_ZSCORES_FLUCC_US_FLUHC_NC = t(NC_ZSCALE_FLUCC_US_FLUHC_NC)

FLUCC_UP_FLUHC_NC_HEATMAP = pheatmap(NC_ZSCORES_FLUCC_US_FLUHC_NC, cluster_rows = T, cluster_cols = F)

# FLU-COLD DOWN FLU HEALTH NO CHANGE FLUCC_DS_FLUHC_NC

SIG_FLUCC_DS_FLUHC_NC_TABLE = FLUCC_DS_FLUHC_NC[c(2:43)]
row.names(SIG_FLUCC_DS_FLUHC_NC_TABLE) = FLUCC_DS_FLUHC_NC$ID
NC_ZTR_FLUCC_DS_FLUHC_NC = na.omit(t(SIG_FLUCC_DS_FLUHC_NC_TABLE))
NC_ZSCALE_FLUCC_DS_FLUHC_NC = scale(NC_ZTR_FLUCC_DS_FLUHC_NC, center = TRUE, scale = TRUE)
NC_ZSCORES_FLUCC_DS_FLUHC_NC = t(NC_ZSCALE_FLUCC_DS_FLUHC_NC)

FLUCC_DN_FLUHC_NC_HEATMAP = pheatmap(NC_ZSCORES_FLUCC_DS_FLUHC_NC, cluster_rows = T, cluster_cols = F)

# FLU-COLD NO CHANGE FLU HEALTH NO CHANGE FLUCC_NC_FLUHC_NC

SIG_FLUCC_NC_FLUHC_NC_TABLE = FLUCC_NC_FLUHC_NC[c(2:43)]
row.names(SIG_FLUCC_NC_FLUHC_NC_TABLE) = FLUCC_NC_FLUHC_NC$ID
NC_ZTR_FLUCC_NC_FLUHC_NC = na.omit(t(SIG_FLUCC_NC_FLUHC_NC_TABLE))
NC_ZSCALE_FLUCC_NC_FLUHC_NC = scale(NC_ZTR_FLUCC_NC_FLUHC_NC, center = TRUE, scale = TRUE)
NC_ZSCORES_FLUCC_NC_FLUHC_NC = t(NC_ZSCALE_FLUCC_NC_FLUHC_NC)

FLUCC_NC_FLUHC_NC_HEATMAP = pheatmap(NC_ZSCORES_FLUCC_NC_FLUHC_NC, cluster_rows = T, cluster_cols = F)

# FLU-COLD NO CHANGE FLU HEALTH DOWN  FLUCC_NC_FLUHC_DS

SIG_FLUCC_NC_FLUHC_DS_TABLE = FLUCC_NC_FLUHC_DS[c(2:43)]
row.names(SIG_FLUCC_NC_FLUHC_DS_TABLE) = FLUCC_NC_FLUHC_DS$ID
NC_ZTR_FLUCC_NC_FLUHC_DS = na.omit(t(SIG_FLUCC_NC_FLUHC_DS_TABLE))
NC_ZSCALE_FLUCC_NC_FLUHC_DS = scale(NC_ZTR_FLUCC_NC_FLUHC_DS, center = TRUE, scale = TRUE)
NC_ZSCORES_FLUCC_NC_FLUHC_DS = t(NC_ZSCALE_FLUCC_NC_FLUHC_DS)

FLUCC_NC_FLUHC_DN_HEATMAP = pheatmap(NC_ZSCORES_FLUCC_NC_FLUHC_DS, cluster_rows = T, cluster_cols = F)

# FLU-COLD NO CHANGE FLU HEALTH UP FLUCC_NC_FLUHC_US

SIG_FLUCC_NC_FLUHC_US_TABLE = FLUCC_NC_FLUHC_US[c(2:43)]
row.names(SIG_FLUCC_NC_FLUHC_US_TABLE) = FLUCC_NC_FLUHC_US$ID
NC_ZTR_FLUCC_NC_FLUHC_US = na.omit(t(SIG_FLUCC_NC_FLUHC_US_TABLE))
NC_ZSCALE_FLUCC_NC_FLUHC_US = scale(NC_ZTR_FLUCC_NC_FLUHC_US, center = TRUE, scale = TRUE)
NC_ZSCORES_FLUCC_NC_FLUHC_US = t(NC_ZSCALE_FLUCC_NC_FLUHC_US)

FLUCC_NC_FLUHC_UP_HEATMAP = pheatmap(NC_ZSCORES_FLUCC_NC_FLUHC_US, cluster_rows = T, cluster_cols = F)

#METAGENE BOXPLOT: average expression accross several genes per sample. Average expression z-score per column (task 3). X:sample groups, y:exp
#meta-gene expression level colMeans() on the per gene expression z-scores for each profile

# CCHC_US_FLUCC_DS, CCHC_DS_FLUCC_DS,
# CCHC_US_FLUHC_DS, CCHC_DS_FLUHC_US, 
# FLUCC_US_FLUHC_DS, FLUCC_DS_FLUHC_US, have 0 variables

# METAGENES OF SIGNATURES: CCHC_UP_FLUCC_NC, CCHC_NC_FLUCC_UP, CCHC_UP_FLUHC_UP, CCHC_NC_FLUHC_UP, FLUCC_UP_FLUHC_UP, FLUCC_NC_FLUHC_UP

# Plots of metagenes are not included in the report as previous plots answer sufficiently the questions and thus they are not saved to disk.

cbPalette=c("dark red","dark green", "dark blue")

Metagene_CCHC_US_FLUCC_NC = data.frame(colMeans(NC_ZSCORES_CCHC_US_FLUCC_NC))
BOX_CCHC_US_FLUCC_NC = ggplot(Metagene_CCHC_US_FLUCC_NC, aes(x=sample_sheet$GROUP, y = colMeans.NC_ZSCORES_CCHC_US_FLUCC_NC., fill = sample_sheet$GROUP))+
  geom_boxplot() + 
  geom_jitter(shape = 19, aes(color = sample_sheet$GROUP)) +
  scale_color_manual(values = cbPalette)+ # Jitter color palette
  theme_Alej() + 
  labs(title = "Metagene CCvsHC up / FLUvsCC no change",x = "Group", y = "Gene exp")

Metagene_CCHC_NC_FLUCC_US = data.frame(colMeans(NC_ZSCORES_CCHC_NC_FLUCC_US))
BOX_CCHC_NC_FLUCC_US = ggplot(Metagene_CCHC_NC_FLUCC_US, aes(x=sample_sheet$GROUP, y = colMeans.NC_ZSCORES_CCHC_NC_FLUCC_US.,fill = sample_sheet$GROUP))+
  geom_boxplot() + 
  geom_jitter(shape = 19, aes(color = sample_sheet$GROUP)) +
  scale_color_manual(values = cbPalette)+ # Jitter color palette
  theme_Alej() + 
  labs(x = "Group", y = "Gene exp",title = "Metagene CCvsHC no change / FLUvsCC up")

Metagene_CCHC_US_FLUHC_US = data.frame(colMeans(NC_ZSCORES_CCHC_US_FLUHC_US))
BOX_CCHC_US_FLUHC_US = ggplot(Metagene_CCHC_US_FLUHC_US, aes(x=sample_sheet$GROUP, y = colMeans.NC_ZSCORES_CCHC_US_FLUHC_US.,fill = sample_sheet$GROUP))+
  geom_boxplot() + 
  geom_jitter(shape = 19, aes(color = sample_sheet$GROUP)) +
  scale_color_manual(values = cbPalette)+ # Jitter color palette
  theme_Alej() + 
  labs(x = "Group", y = "Gene exp",title = "Metagene CCvsHC up / FLUvsHC up")

Metagene_CCHC_NC_FLUHC_US = data.frame(colMeans(NC_ZSCORES_CCHC_NC_FLUHC_US))
BOX_CCHC_NC_FLUHC_US = ggplot(Metagene_CCHC_NC_FLUHC_US, aes(x=sample_sheet$GROUP, y = colMeans.NC_ZSCORES_CCHC_NC_FLUHC_US.,fill = sample_sheet$GROUP))+
  geom_boxplot() + 
  geom_jitter(shape = 19, aes(color = sample_sheet$GROUP)) +
  scale_color_manual(values = cbPalette)+ # Jitter color palette
  theme_Alej() + 
  labs(x = "Group", y = "Gene exp",title = "Metagene CCvsHC no change / FLUvsHC up")

Metagene_FLUCC_US_FLUHC_US = data.frame(colMeans(NC_ZSCORES_FLUCC_US_FLUHC_US))
BOX_FLUCC_US_FLUHC_US = ggplot(Metagene_FLUCC_US_FLUHC_US, aes(x=sample_sheet$GROUP, y = colMeans.NC_ZSCORES_FLUCC_US_FLUHC_US.,fill = sample_sheet$GROUP))+
  geom_boxplot() + 
  geom_jitter(shape = 19, aes(color = sample_sheet$GROUP)) +
  scale_color_manual(values = cbPalette)+ # Jitter color palette
  theme_Alej() + 
  labs(x = "Group", y = "Gene exp",title = "Metagene FLUvsCC up / FLUvsHC up")

Metagene_FLUCC_NC_FLUHC_US = data.frame(colMeans(NC_ZSCORES_FLUCC_NC_FLUHC_US))
BOX_FLUCC_NC_FLUHC_US = ggplot(Metagene_FLUCC_NC_FLUHC_US, aes(x=sample_sheet$GROUP, y = colMeans.NC_ZSCORES_FLUCC_NC_FLUHC_US.,fill = sample_sheet$GROUP))+
  geom_boxplot() + 
  geom_jitter(shape = 19, aes(color = sample_sheet$GROUP)) +
  scale_color_manual(values = cbPalette)+ # Jitter color palette
  theme_Alej() + 
  labs(x = "Group", y = "Gene exp",title = "Metagene FLUvsCC no change / FLUvsHC up")

#GENE ONTOLOGY - ENRICHMENT with webgestalt - genome as reference gene list

#Gene lists for each comparison 
#Each gene list is then exported to webgestalt to perform the Over Representation Analysis with the "genome" reference gene list option
#to obtain a plot depicting the biological activity of each comparison.

geneList_CCHC=master_file[which(master_file$sig_CC_HC=="TRUE"),]
geneList_CCHC_table=geneList_CCHC[c(1:29)]
geneList_CCHC_ID = geneList_CCHC_table[1]
write.table(geneList_CCHC_ID, file = "geneList_CCHC_ID.txt", sep = "",row.names = FALSE, col.names = FALSE, quote = FALSE)

geneList_FLUHC=master_file[which(master_file$sig_FLU_HC=="TRUE"),]
geneList_FLUHC_table=geneList_FLUHC[c(1:15,29:43)]
geneList_FLUHC_ID = geneList_FLUHC_table[1]
write.table(geneList_FLUHC_ID, file = "geneList_FLUHC_ID.txt", sep = "",row.names = FALSE, col.names = FALSE, quote = FALSE)

geneList_FLUCC=master_file[which(master_file$sig_FLU_CC=="TRUE"),]
geneList_FLUCC_table=geneList_FLUCC[c(1,16:43)]
geneList_FLUCC_ID = geneList_FLUCC_table[1]
write.table(geneList_FLUCC_ID, file = "geneList_FLUCC_ID.txt", sep = "",row.names = FALSE, col.names = FALSE, quote = FALSE)
