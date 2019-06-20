##########################################################
## GSE443 - Re-analysis of affymetrix one-channel array
## Thyago Leal Calvo - thyagoleal@yahoo.com 
## Last update: June 9th 2018
##########################################################

## Settings -------------------------------

set.seed(1618)
options(digits = 3, download.file.method.GEOquery = NULL, scipen = 999)
root <- getwd()

## Loading packages -------------------------------

require(GEOquery)
require(limma)
require(org.Hs.eg.db)
require(AnnotationDbi)
require(WriteXLS)
require(hgu95av2.db)

## Loading data from GEO -------------------------------

gset <- getGEO(GEO = "GSE443",
               destdir = "raw/",
               GSEMatrix = TRUE, getGPL = TRUE, AnnotGPL = TRUE)

gset <- gset$GSE443_series_matrix.txt.gz

dim(gset)
boxplot(log2(exprs(gset)))

summary(exprs(gset))

exp <- exprs(gset)

## Annotating features -------------------------

annt <- select(hgu95av2.db, keys = row.names(exp), columns = c("ENTREZID", "SYMBOL", "GENENAME", "REFSEQ"))

annt <- annt[!duplicated(annt$PROBEID),]
row.names(annt) <- annt$PROBEID

annt <- annt[row.names(exp),]

# Check

identical(row.names(annt), row.names(exp))

featureData(gset) <- new("AnnotatedDataFrame", annt)

head(fData(gset))

## Fixing PhenoData ------------------------------

pd <- data.frame(Sample = sampleNames(gset), Group = factor(sub(".*\\((.*)\\).*", "\\1", pData(gset)[,1])), row.names = sampleNames(gset))

phenoData(gset) <- new("AnnotatedDataFrame", pd)

gset

plotMDS(gset, labels = pData(gset)$Group, col = as.integer(pData(gset)$Group))

# Removing genes with negative raw expression
keep <- apply(exp, 1,  function(row) all(row > 0))
table(keep)

gset <- gset[keep, ]

rm(keep)

exp <- log2(exprs(gset))

gset <- ExpressionSet(assayData = exp, phenoData = gset@phenoData, featureData = gset@featureData)

plotMDS(gset, labels = pData(gset)$Group, col = as.integer(pData(gset)$Group))

saveRDS(gset, "objects/eSet.Rds", compress = "bzip2")
rm(exp, exp2, annt)

## Analysis ----------------------------

groups <- pData(gset)$Group
groups

design <- model.matrix(~groups)

fit <- lmFit(gset, design)

## Removing duplicate ENTREZ ---------------------------

fit <- fit[order(fit$Amean, decreasing = TRUE), ]
fit <- fit[!duplicated(fit$genes$ENTREZID), ]
fit <- fit[!is.na(fit$genes$ENTREZID),]
dim(fit)

fit2 <- eBayes(fit)

tt <- topTable(fit2, coef = 2, number = Inf, adjust.method = "fdr", confint = TRUE)

## Saving results ------------

WriteXLS(tt, "results/DEG_full_results.xls", row.names = FALSE)

exp_matrix <- exprs(gset)

exp_matrix <- exp_matrix[row.names(fit2),]
exp_matrix <- as.data.frame(exp_matrix)
colnames(exp_matrix) <- make.names(pData(gset)$Group, unique = TRUE)
exp_matrix$entrez <- fit2$genes$ENTREZID
exp_matrix$symbol <- fit2$genes$SYMBOL

write.table(exp_matrix, "results/exp_matrix_gse443.txt", quote = FALSE, sep = "\t", row.names = TRUE)

### Volcano plot =================

library(RColorBrewer)
library(gplots)
library(ggplot2)
library(limma)

Significant = tt$P.Value <= 0.05 & abs(tt$logFC) > 1

plot <- ggplot(data = tt, 
                     aes(x = logFC,
                         y = -log10(P.Value), 
                         colour = Significant)) +
        geom_point(alpha = 1, size = .8) +
        xlab("log2 Fold Change") + ylab("-log10 P-value") + theme_bw() +
        theme(legend.position = "none") + scale_color_manual(values=c("lightgrey", "black"))

plot <- plot + ggtitle("Differentially Expressed Genes from GSE443 - Leprosy Lesion from LL vs BT patients") + geom_hline(yintercept = -log10(0.05), colour = "dodgerblue", linetype = 3)+geom_vline(xintercept = c(1,-1), colour = c("red", "limegreen"), linetype = 3) 

pdf("figs/VolcanoPlot.pdf")
plot
dev.off()

# Heatmap ----------------------------

library(pheatmap)
library(RColorBrewer)

plot = tt$PROBEID[tt$P.Value< 0.01]

label = data.frame(Class = gset@phenoData@data$Group)
row.names(label) <- gset@phenoData@data$Sample

pdf("figs/heatmap_deg.pdf", onefile = FALSE)
pheatmap(exprs(gset[featureNames(gset) %in% plot,]),
         color = brewer.pal(11, "RdYlBu"),
         scale = "row",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         legend = TRUE,
         show_rownames = FALSE,
         main = paste("Heatmap of GSE443 DEGs (", length(unique(plot)),  ") with P-value <= 0.01", sep = ""),
         border_color = NA, 
         annotation_col = label,
         annotation_legend = TRUE,
         show_colnames = FALSE)
dev.off()

###  Gene Ontology Enrichment Analysis  -----------------------

universe <- tt$Gene.ID

deg <- tt[tt$P.Value < 0.05,]
deg <- deg[,c("Gene.ID", "logFC", "P.Value")]
deg <- deg[!duplicated(deg$Gene.ID),]

geneList <- deg$logFC
names(geneList) <- as.character(deg$Gene.ID)

geneList <- sort(geneList, decreasing = TRUE) ## geneList with Log2FC + Entrez for DEG with P-val < 0.05
head(geneList)

library(clusterProfiler)
library(org.Hs.eg.db)
library(pathview)

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

pdf("figs_pdf/CellComponent.pdf", width = 10, pagecentre = TRUE, onefile = TRUE)

barplot(go.cc, 
        colorBy = "p.adjust",
        title = "Gene Ontology over-representation test for Cellular Component",
        showCategory = length(go.cc@result$ID), font.size = 10)
dev.off()

cnetplot(go.cc, categorySize="pvalue", foldChange=geneList, fixed = FALSE)

write.csv(go.cc@result[,1:7], file = "results/go.cc.csv", row.names = FALSE)

## Biological Process  ------------------

go.bp <-  enrichGO(gene = names(geneList),
                   universe      = universe,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "BP",
                   pAdjustMethod = "fdr",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

go.bp@result[,1:7]

pdf("figs_pdf/BioProcess.pdf", width = 10, pagecentre = TRUE, onefile = TRUE)

barplot(go.bp, 
        colorBy = "p.adjust",
        title = "Gene Ontology over-representation test for Biological Process",
        showCategory = length(go.bp@result$ID), font.size = 10)
dev.off()

cnetplot(go.bp, categorySize="pvalue", foldChange=geneList, fixed = FALSE)

write.csv(go.bp@result[,1:7], file = "results/go.bp.csv", row.names = FALSE)

## Molecular Function  ------------------

go.mf <-  enrichGO(gene = names(geneList),
                   universe      = universe,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "Mf",
                   pAdjustMethod = "fdr",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

go.mf@result[,1:7]

pdf("figs_pdf/MolecularFunction.pdf", width = 10, pagecentre = TRUE, onefile = TRUE)

barplot(go.mf, 
        colorBy = "p.adjust",
        title = "Gene Ontology over-representation test for Molecular Function",
        showCategory = length(go.mf@result$ID), font.size = 10)
dev.off()

cnetplot(go.mf, categorySize="pvalue", foldChange=geneList, fixed = FALSE)

write.csv(go.mf@result[,1:7], file = "results/go.mf.csv", row.names = FALSE)

## Saving Table DEG --------------

DEG <- topTable(fit2, number = Inf, genelist = fit$genes, adjust.method = "fdr", p.value = 0.05, lfc = 0 )

DEG$logFC <- with(DEG, round(logFC, 3))

### END  ########
rm(list=ls())
q(save = "no")
        