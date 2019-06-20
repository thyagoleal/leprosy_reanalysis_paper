##  Re-analysis of two-channel microarray GSE40950
##  2017-02-26 - Thyago Leal Calvo
## Reviewed July, 18th 2017.

root <- getwd()
options(digits = 3, sciphen = 999)
set.seed(128192)

## Load packages ---------------

library(GEOquery)
library(marray)
library(limma)

## Loading RAW files ---------------

gpr_files <- list.files(path = "raw_files/",pattern = "gpr", full.names = TRUE)
gpr_files

# 1 - Galfile -----

galinfo <- read.Galfile("raw_files/GPL16068_hoae.gal", 
                        path = getwd(),
                        info.id = c("ID", "Name"),
                        layout.id = c(Block = "Block", Row = "Row", Column = "Column"),
                        labels = "ID", 
                        sep = "\t", 
                        skip = 55)

# 2 - target info (I made it using info published in GEO) -------

targets <- readTargets("raw_files/samplesInfo.txt")

# Intensity files  ----------

columns <- c("F647 Median", "F555 Median", "B647 Mean", "B555 Mean")
names(columns) <- c("R", "G", "Rb", "Gb")

RG <- read.maimages(gpr_files,
                    source = "genepix.custom", 
                    columns = columns, 
                    green.only = FALSE,
                    wt.fun = NULL,
                    verbose = TRUE)

# Checking file quality --------

png("figs/Diag_plot.png", antialias = "default")
par(mfrow = c(1,2))
imageplot(log2(RG$Rb[,1]), RG$printer, low="white", high="red")
imageplot(log2(RG$Gb[,1]), RG$printer, low="white", high="green")
dev.off()

png("figs/RG_channel_densities.png", antialias = "default")
par(mfrow = c(1,1))
plotDensities(RG)
dev.off()

## Background correction --------

RG_bg <- backgroundCorrect(RG, 
                           method = "normexp", offset = 0)

summary(log2(RG_bg$G))

## Assessing quality

# Before and after BG correction:
# before

png("figs/BG_correction.png", width = 1200, height = 800)
par(mfrow=c(1,2))
plotMD(RG, main = "RG Without BG Correction")

# after
plotMD(RG_bg, main =  "RG With BG Correction")
dev.off()
par(mfrow=c(1,1))

RG_bg$genes <- galinfo$gnames@maInfo
names(RG_bg$printer)
RG_bg$printer <- getLayout(RG$genes)


## Checking distribution of Red and Green Values 
library(geneplotter)

pdf(file = "figs/Distribution_density_boxplot.pdf", paper = "a4r")
par(mfrow=c(2,1))

plotformula_Green = log2(RG_bg$G)~col(RG_bg$G)
boxplot(plotformula_Green, outline=FALSE,
          col="limegreen", xlab="Arrays",
          ylab=expression(log[2]~G), main="Boxplot of the Green Intensities from All Arrays")

plotformula_Red = log2(RG_bg$R)~col(RG_bg$R)
boxplot(plotformula_Red, outline=FALSE,
        col="red2", xlab="Arrays",
        ylab=expression(log[2]~R), main="Boxplot of the Red Intensities from All Arrays")

dev.off()

pdf(file = "figs/Distribution_density.pdf", paper = "a4r")
par(mfrow=c(2,1))

multidensity(plotformula_Green,
             main="Green Densities", xlab = expression(log[2]~Green))

multidensity(plotformula_Red,
             main="Red Densities", xlab = expression(log[2]~Red))

dev.off()

pdf(file = "figs/RG_densities.pdf", paper = "a4r")
plotDensities(RG_bg, main = "Red and Green Densities combined")
dev.off()

## Normalization ------- 

# 1 Robust Spline

MA <- normalizeWithinArrays(RG_bg, method = "robustspline")

# 2 Between Array Norm quantile method 

MA2 <- normalizeBetweenArrays(MA, method = "quantile")

pdf(file = "figs/Normalizatiom_MAplots.pdf", width = 1200, height = 600, paper = "a4r")
par(mfcol=c(1,3), oma=c(1,1,1,1))

limma::plotMA(RG_bg, main = "Before normalization")
limma::plotMA(MA, main = "After Robust spline Normalization")
limma::plotMA(MA2, main = "After Quantile Normalization")
dev.off()

## Removing uninformative probes and duplicates 

dim(MA2)

MA2 = MA2[!grepl("suid[.]*", MA2$genes$ID),]
MA2 = MA2[!grepl("EST[_]*", MA2$genes$ID),]
MA2 = MA2[!grepl("^MJ", MA2$genes$ID),]
MA2 = MA2[!grepl("AmbionSpike*", MA2$genes$ID),]
MA2 = MA2[!grepl("Random_*", MA2$genes$ID),]
MA2 = MA2[!grepl("empty", MA2$genes$ID),]
MA2 = MA2[!grepl("*[_]REVCOMP", MA2$genes$ID),]
MA2 = MA2[!grepl("BCR([\\])*", MA2$genes$ID),]
dim(MA2)

## Setting up study design ---------------------

design <- modelMatrix(targets, ref = "non-infected")

### Fitting linear model 

fit <- lmFit(MA2, design)

### Empirical Bayes analysis 
# The moderated t-statistics use sample standard deviations which have been squeezed
# towards a pooled standard deviation value

fit <- eBayes(fit)

## Summary table of some key statistics for the top genes. 

tt <- topTable(fit, adjust = "fdr", number = Inf, confint = TRUE)

tt_2 <-
        subset(
                tt,
                select = c(
                        "ID",
                        "Name",
                        "adj.P.Val",
                        "P.Value",
                        "logFC",
                        "CI.L",
                        "CI.R",
                        "t",
                        "B",
                        "AveExpr"
                )
        )
head(tt_2)

## Annotating with Entrez ID -------------------

library(Homo.sapiens)
library(AnnotationDbi)

entrez <- select(Homo.sapiens, keys = tt_2$ID, columns = "ENTREZID", keytype = "SYMBOL")
full_deg <- merge.data.frame(tt_2, entrez, by.x = "ID", by.y = "SYMBOL", all = FALSE)
full_deg <- full_deg[!duplicated(full_deg$ENTREZID),]
full_deg <- full_deg[!is.na(full_deg$ENTREZID),]

## Saving tables ------

library(WriteXLS)

WriteXLS(full_deg, ExcelFileName = "Results/Results_GSE40950.xls")

###### Extracting channels 

exprs_40950 = MA2$M
colnames(exprs_40950) <- targets$Cy5

entrez <- select(org.Hs.eg.db, keys = MA2$genes$ID, columns = c("ENTREZID"), keytype = "SYMBOL")
rownames(exprs_40950) <- rownames(MA2$genes)

write.table(exprs_40950, file = "Results/exp_40950.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

### Volcano plot =================

library(RColorBrewer)
library(gplots)
library(ggplot2)
library(limma)

## Grab all genes with their p-values and Log2FC to a new data frame:
data_for_plot <- topTable(fit, adjust.method = "fdr", number = Inf)

Significant = data_for_plot$P.Value < 0.05 & abs(data_for_plot$logFC) > 1

plot <- ggplot(data = data_for_plot, 
               aes(x = logFC,
                   y = -log10(P.Value), 
                   colour = Significant)) +
        geom_point(alpha = 0.8, size = .8) +
        xlab("log2 Fold Change") + ylab("-log10 P-value") +
        theme_bw() 

gA <- plot + ggtitle("Differentially Expressed Genes from GSE40950 - Infected vs Control") + theme(legend.position = "none") + scale_color_manual(values=c("lightgrey", "black"))+ geom_hline(yintercept = -log10(0.05), colour = "dodgerblue", linetype = 3) + geom_vline(xintercept = c(-1,1), colour = c("green", "red"), linetype = 3)

pdf("figs_PDF/VolcanoPlot.pdf", paper = "a4r", pagecentre = TRUE)
gA
dev.off()

## Heatmap --------------------

library(pheatmap)
library(RColorBrewer)

plot = tt_2$ID[tt_2$P.Value < 0.05]

label = data.frame(Condition = colnames(exprs_40950))
row.names(label) <- make.names(label$Condition, unique = T)
ann_colors = list(Condition = c(infected = "firebrick1", non.infected = "turquoise2"))

pdf("figs_PDF/heatmap_deg.pdf", paper = "a4", pagecentre = TRUE, onefile = TRUE)
pheatmap(exprs_40950[rownames(exprs_40950) %in% plot,],
         color = brewer.pal(12, "RdYlBu"),
         scale = "row",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "euclidean",
         clustering_method = "average",
         legend = TRUE,
         show_rownames = FALSE,
         main = "Heatmap of GSE40950 DEGs (703) with P-value < 0.05",
         border_color = NA, 
         annotation_col = label,
         annotation_legend = TRUE,
         annotation_colors = ann_colors, show_colnames = FALSE)
dev.off()

save.image(compress = TRUE)

### Venn Diagram  ---------------------
 
 results <- decideTests(fit, adjust.method = "fdr", p.value = 0.05)
 
 setEPS()
 postscript("VennDiagram_GSE22784.eps", paper = "a4", pagecentre = TRUE)
 
 vennDiagram(results, include = c("up", "down"), 
             counts.col = c("red", "blue"), 
             circle.col = c("green3", "yellow", "darkblue", "red"),
             names = "Infected vs Control, p < 0.05", cex = c(1,1,1))
 
dev.off()

#### Extracting separate channels in order to get matrix with raw samples per channel/dye  ------------

green <- c("control", "infected", "infected", "control", "infected", "control")
red <- c("infected", "control", "control", "infected", "control", "infected")

df_green <- RG_bg$G
df_red <- RG_bg$R

colnames(df_green) <- green
colnames(df_red) <- red

expr_matrix_dual_channel <- cbind(df_green, df_red)
dim(expr_matrix_dual_channel)

head(RG_bg$genes$ID)
length(RG_bg$genes$ID)

genes <- as.character(RG_bg$genes$ID)
expr_matrix_dual_channel <- cbind(genes, expr_matrix_dual_channel)

cutoff <- tt_2$ID[tt_2$P.Value <= 0.05]

expr_matrix_dual_channel <- expr_matrix_dual_channel[expr_matrix_dual_channel[,1] %in% cutoff,]

input <- as.data.frame(t(expr_matrix_dual_channel))

input <- input[-1,]

colnames(input) <- genes[genes %in% cutoff]

row.names(input) <- NULL

samples <- colnames(expr_matrix_dual_channel[, 2:13])
input <- cbind(samples, input)

dup_ges <- !duplicated(colnames(input))
table(dup_ges)

input <- input[,dup_ges]

input[,-1] <- sapply(input[,-1], function(x){as.numeric(as.character(x))})

str(input)

input[,-1] <- log2(input[,-1])


### mSVM-RFE Analysis --------------------------

# Set SVM-RFE directory. Then work through these commands.

set.seed(12345)
library(e1071)
source('/home/thyago/R/x86_64-pc-linux-gnu-library/3.4/SVM-RFE-master/msvmRFE.R')

# Basic usage
svmRFE(input, k = 10, halve.above = 10)

# Set up cross validation
nfold = 10
nrows = nrow(input)
folds = rep(1:nfold, len=nrows)[sample(nrows)]
folds
folds = lapply(1:nfold, function(x) which(folds == x))
folds

# Perform feature ranking on all training sets
results = lapply(folds, svmRFE.wrap, input, k = 10, halve.above = 10)
length(results)
results

# Obtain top features across ALL folds
top.features = WriteFeatures(results, 
                             input, 
                             save = F)

write.table(top.features, file = "Results/top.features_40950.txt", quote = FALSE, sep = "\t", row.names = FALSE)

head(top.features)

# Estimate generalization error using a varying number of top features 
featsweep = lapply(1:5, FeatSweep.wrap, results, input)
featsweep

# Make plot
no.info = min(prop.table(table(input[,1])))
errors = sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))

PlotErrors(errors, no.info = no.info)

## Test the features using PCA  --------------------
require(FactoMineR)
require(factoextra)

tmp = t(exprs_40950[rownames(exprs_40950) %in% top.features$FeatureName[1:10],])
class <- row.names(tmp)
row.names(tmp) <- NULL


pca1 <- PCA(tmp,  scale.unit = TRUE)

pdf("figs_PDF/PCA_top10.pdf")
fviz_pca_biplot(pca1, geom = "point", repel = TRUE,
                point.shape = 19,
                habillage = as.factor(class),
                pointsize = 3,
                title = "PCA Top 10 Genes from Feature Selection",
                col.var = "#A3A3A3") + theme_bw()
dev.off()

## Heatmap with the top 10 Rank from Feature Selection --------------------

pdf("figs_PDF/Heatmap_top50.pdf", onefile = TRUE, paper = "a4")
pheatmap(exprs_40950[rownames(exprs_40950) %in% top.features$FeatureName[1:50],],
         color = brewer.pal(11, "RdYlBu"),
         scale = "row",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "euclidean",
         clustering_method = "average",
         legend = TRUE,
         show_rownames = TRUE,
         main = "Heatmap of Top 50 DEG from Feature Selection Step",
         border_color = NA, 
         annotation_col = label,
         annotation_legend = TRUE,
         fontsize = 8,
         annotation_colors = ann_colors,
         show_colnames = FALSE)
dev.off()

# Gene Ontology Enrichment Analysis  ------------------

### Loading required packages ==============

library(clusterProfiler)
library(org.Hs.eg.db)
library(pathview)

#### Getting Universe genes ==================

universe <- full_deg$ENTREZID

# grab only gene symbols

## Getting changed genes ============

# Since we've got no genes with the desired adj.P-value threshold, 0.05, I'm gonna select based on p-value only.

table(full_deg$P.Value < 0.05)
deg <- full_deg[full_deg$P.Value <= 0.05, ]
dim(deg)

deg <- deg[, c(9,3,6)] # selecting only GENE SYMBOL, log2FC and p-value columns.

unique <- unique(deg$ENTREZID)
deg <- deg[deg$ENTREZID %in% unique, ]
dim(deg)
table(is.na(deg$ENTREZID))

## Converting gene SYMBOL to ENTREZ ID using annotate and organisms.

# DEG

genes <- deg$ID

genes = bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
universe = bitr(universe, fromType = "SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

deg <- deg[deg$ID %in% genes$SYMBOL, ]
deg <- deg[!duplicated(deg$ID),]

deg <- cbind(deg, genes$ENTREZID)

colnames(deg)[1] <- "Entrez"

# Creating geneList object

## feature 1: numeric vector
geneList = deg[, 2]

## feature 2: named vector
names(geneList) = as.character(deg[, 1])

## feature 3: decreasing order
geneList = sort(geneList, decreasing = TRUE)
head(geneList)

## Cellular Component ------------------

go.cc <-  enrichGO(gene = names(geneList),
                   universe      = universe,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "CC",
                   pAdjustMethod = "fdr",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

go.cc@result[,1:7]

go.cc <- simplify(go.cc, cutoff = 0.7, by = "p.adjust", select_fun = min)
pdf("figs_PDF/go.cc.pdf", paper = "a4r", pagecentre = TRUE)
barplot(go.cc, 
        colorBy = "p.adjust",
        title = "Gene Ontology over-representation test for Cellular Component",
        showCategory = length(go.cc@result$ID), font.size = 8)
dev.off()

cnetplot(go.cc, categorySize="pvalue", foldChange=geneList, fixed = FALSE)

write.csv(go.cc@result[,1:7], file = "Results/GO_cc.csv", row.names = FALSE)

## Biological Process  ------------------

go.bp <-  enrichGO(gene = names(geneList),
                   universe      = universe,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "BP",
                   pAdjustMethod = "fdr",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

go.bp

go.bp <- simplify(go.bp, cutoff = 0.5, by = "p.adjust", select_fun = min)

View(go.bp@result[,1:5])

pdf("figs_PDF/go.bp.pdf", paper = "a4r")
barplot(go.bp, 
        colorBy = "p.adjust",
        title = "Gene Ontology over-representation test for Cellular Component",
        showCategory = length(go.bp@result$ID), font.size = 8)

dev.off()
write.csv(go.bp@result[,1:7], file = "Results/GO_bp.csv", row.names = FALSE)


## Molecular function  ------------------

go.mf <-  enrichGO(gene = names(geneList),
                   universe      = universe,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "MF",
                   pAdjustMethod = "fdr",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

go.mf@result

go.mf <- simplify(go.mf, cutoff = 0.3, by = "p.adjust", select_fun = min)

View(go.mf@result[,1:7])

pdf("figs_PDF/go.mf.pdf", paper = "a4r")

barplot(go.mf, 
        colorBy = "p.adjust",
        title = "Gene Ontology over-representation test for Cellular Component",
        showCategory = length(go.mf@result$ID), font.size = 8)

dev.off()
write.csv(go.mf@result[,1:7], file = "Results/GO_mf.csv", row.names = FALSE)

## Reactome Enrichmment Analysis ---------------

reac <- enrichPathway(names(geneList), pAdjustMethod = "fdr", universe = universe$ENTREZID, readable = TRUE)

View(reac@result[,1:7])
dotplot(reac, showCategory = 10, color)
enrichMap(reac, layout=igraph::layout.kamada.kawai, vertex.label.cex = 1, fixed = FALSE)

cnetplot(reac, categorySize="pvalue", foldChange=geneList, fixed = FALSE)


### Exporting DEG table -----------

DEG <- topTable(fit, number = Inf, genelist = fit$genes, adjust.method = "fdr")

library(Homo.sapiens)
library(AnnotationDbi)

DEG <- DEG[DEG$P.Value < 0.005,]

entrez <- select(Homo.sapiens, keys = DEG$ID, columns = "ENTREZID", keytype = "SYMBOL")

DEG <- merge.data.frame(entrez, DEG, by.x = "SYMBOL", by.y = "ID", all = FALSE)

DEG <- DEG[!is.na(DEG$ENTREZID),]
DEG <- DEG[!duplicated(DEG$ENTREZID),]

DEG <- DEG[, c("ENTREZID", "logFC", "adj.P.Val")]

write.table(DEG, file = "Results/DEG.GSE40950.InfvsCon.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


########################
####      END       ####
########################

rm(list=ls())
q(save = "no")




