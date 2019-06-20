##############################################################################
##  Re-analysis of single-channel microarray of Anna Beatriz - ABI - 2012    #
##  2017-02-27
##  Reviwed May 2018.
##############################################################################

set.seed(210920)
options(digits = 3)

### Loading packages ----------------------------------------------

library(GEOquery)
library(limma)
library(WriteXLS)
library(pheatmap)
library(RColorBrewer)

#### Reading the expression matrix already normalized  ---------------------------

gset <- getGEO(GEO = "GSE35423", GSEMatrix = TRUE, destdir = "raw/")
gset <- gset$GSE35423_series_matrix.txt.gz

boxplot(exprs(gset))

## Creating the phenotype object from data published  ------------------------

conditions <- pData(gset)["source_name_ch1"]
array <- row.names(pData(gset))
conditions <- data.frame(infection = as.factor(c("infected", "infected", "control", "control", "infected", "infected", "control", "control")), time = as.factor(c(rep(24, 4), rep(48, 4))))
p.data <- data.frame(Array = array, Cond = conditions)
p.data
row.names(p.data) <- p.data$Array
colnames(exprs(gset)) == row.names(p.data)

boxplot(exprs(gset))

saveRDS(gset, file = "obj/eSet.rds", compress = "bzip2")

## Getting feature info ------------------------

features <- gset@featureData@data

## Seting up comparisons ------------------------

factors <- paste(p.data$Cond.infection, p.data$Cond.time, sep=".")
f <- factor(factors)

design <- model.matrix(~0+f)
colnames(design) <- c("control24", "control48", "infected24", "infected48")

fit <- lmFit(gset, design = design)

## Removing duplicated ENTREZID and useless probes --------------------

# Removing obsolete and pseudogenes

table(fit$genes$Status)

fit <- fit[fit$genes$Status == "current",]
fit <- fit[order(fit$Amean, decreasing = TRUE), ]
fit <- fit[!duplicated(fit$genes$GENE), ]
fit <- fit[!is.na(fit$genes$GENE), ]
dim(fit)

cont.matrix <-
        makeContrasts(
                Inf24 = "infected24-control24",
                Inf48 = "infected48-control48",
                Inf48vsInf24 = "infected48-infected24",
                levels = design
        )

fit2 <- contrasts.fit(fit, contrasts = cont.matrix)

fit2 <- eBayes(fit2)

inf24vsCon24 <-topTable(fit2, coef = "Inf24", number = Inf, adjust.method = "fdr", confint = TRUE) 
inf48vsCon48 <-topTable(fit2, coef = "Inf48", number = Inf, adjust.method = "fdr", confint = TRUE) 
inf48vsInf24 <-topTable(fit2, coef = "Inf48vsInf24", number = Inf, adjust.method = "fdr", confint = TRUE) 

### Saving Table with results ---------- 

list <- list("inf24vsCon24", "inf48vsCon48", "inf48vsInf24")

for(i in list){
        filename <- paste("results/",i,".xls", sep = "")
        WriteXLS(x = i, ExcelFileName = filename)
}
rm(i, filename)

## Venn Diagram ---------------------

results <- decideTests(fit2, adjust.method = "fdr", p.value = 0.05, lfc = 1)


pdf("figs_pdf/VennDiagram.pdf", paper = "a4", pagecentre = TRUE)

vennDiagram(results, include = c("up", "down"), 
            counts.col = c("firebrick1", "dodgerblue2"), 
            circle.col = c("limegreen", "yellow", "dodgerblue2", "firebrick1"))

dev.off()

### Volcano plot =================

library(RColorBrewer)
library(gplots)
library(ggplot2)
library(limma)

## Infection 24 h vs Control 24 h -------------------------------

pdf("figs_pdf/VolcanoPlots.pdf", onefile = TRUE, paper = "a4r", pagecentre = TRUE)

coef1 <- topTable(fit2,
                  coef = 1,
                  adjust.method = "fdr",
                  number = Inf)

sig = coef1$P.Value < 0.05 & abs(coef1$logFC) > 1
table(sig)

names(coef1)

plot_coef1 <- ggplot(data = coef1, 
                     aes(x = logFC,
                         y = -log10(P.Value), 
                         colour = sig)) +
        geom_point(alpha = 1, size = .5) +
        xlab("log2 Fold Change") + ylab("-log10 P-value") + theme_bw() +
        theme(legend.position = "none") + scale_color_manual(values=c("lightgrey", "black"))

plot_coef1 <- plot_coef1 + ggtitle("Differentially Expressed Genes from GSE35423- Infection 24 h vs Control 24 h") + geom_hline(yintercept = -log10(0.05), colour = "dodgerblue", linetype = 3) + geom_vline(xintercept = c(-1,1), colour = c("green", "red"), linetype = 3)

plot_coef1

## Infection 48 h vs Control 48 h -------------------------------

coef2 <- topTable(fit2,
                  coef = 2,
                  adjust.method = "fdr",
                  number = Inf)

sig = coef2$P.Value < 0.05 & abs(coef2$logFC) > 1

names(coef2)

plot_coef2 <- ggplot(data = coef2, 
                     aes(x = logFC,
                         y = -log10(P.Value), 
                         colour = sig)) +
        geom_point(alpha = 1, size = .5) +
        xlab("log2 Fold Change") + ylab("-log10 P-value") + theme_bw() +
        theme(legend.position = "none") + scale_color_manual(values=c("lightgrey", "black"))

plot_coef2 <- plot_coef2 + ggtitle("Differentially Expressed Genes from GSE35423- Infection 48 h vs Control 48 h") + geom_hline(yintercept = -log10(0.05), colour = "dodgerblue", linetype = 3) + geom_vline(xintercept = c(-1,1), colour = c("green", "red"), linetype = 3)

plot_coef2

## Infection 48 h vs Control 48 h -------------------------------

coef3 <- topTable(fit2,
                  coef = 3,
                  adjust.method = "fdr",
                  number = Inf)

sig = coef3$P.Value < 0.05 & abs(coef3$logFC) > 1

names(coef3)

plot_coef3 <- ggplot(data = coef3, 
                     aes(x = logFC,
                         y = -log10(P.Value), 
                         colour = sig)) +
        geom_point(alpha = 1, size = .5) +
        xlab("log2 Fold Change") + ylab("-log10 P-value") + theme_bw() +
        theme(legend.position = "none") + scale_color_manual(values=c("lightgrey", "black"))

plot_coef3 <- plot_coef3 + ggtitle("Differentially Expressed Genes from GSE35423- Infection 48 h vs Infection 24 h") + geom_hline(yintercept = -log10(0.05), colour = "dodgerblue", linetype = 3) + geom_vline(xintercept = c(-1,1), colour = c("green", "red"), linetype = 3)

plot_coef3
dev.off()

## Heatmap ---------------------------

plot = inf48vsCon48$adj.P.Val < 1e-3 & abs(inf48vsCon48$logFC)> 1

label = data.frame(Treatment = p.data$Cond.infection, Time = p.data$Cond.time)
row.names(label) <- make.names(p.data$Array, unique = T)

exp <- exprs(gset)
exp <- exp[row.names(exp) %in% row.names(fit2$genes), ]

pdf("figs_pdf/heatmap_deg.pdf", paper = "a4", onefile = FALSE)
pheatmap(exp[plot,],
         color = brewer.pal(9, "RdYlBu"),
         scale = "row",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "euclidean",
         clustering_method = "average",
         legend = TRUE,
         show_rownames = FALSE,
         main = paste("Heatmap of GSE35423 DEGs (", unname(table(plot)[2]),  ") with Adj.P-value < 0.001 & |logFC| > 2", sep = ""),
         border_color = NA, 
         annotation_col = label,
         annotation_legend = TRUE,
         show_colnames = FALSE)
dev.off()

## Preparing the expression matrix  ------------------------------------

input = exprs(gset)
colnames(input) <- paste(p.data$Cond.infection, p.data$Cond.time, sep = ".")
input <- as.data.frame(input)

input<- cbind(input, gset@featureData@data[, c("GENE", "Gene Symbol")])
colnames(input)[9:10] <- c("Entrez", "Symbol")

input <- input[!is.na(input$Entrez),]
input <- input[input$Entrez != "",]
input <- input[order(rowMeans(input[,1:8], na.rm = TRUE), decreasing = TRUE), ]
input <- input[!duplicated(input$Entrez), ]
dim(input)

row.names(input) <- as.character(input$Entrez)

## Exporting ExprMatrix ---------------

write.table(input,file = "results/expr_matrix_35423.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# Gene Ontology Enrichment Analysis  ------------------

### Loading required packages ==============

library(clusterProfiler)
library(org.Hs.eg.db)
library(pathview)

#### Getting Universe genes ==================

universe <- fit$genes$GENE
# grab only gene entrez

## Getting changed genes for infection 24 h -----------

tt <- topTable(fit2, coef = 1, number = Inf, genelist = fit$genes, adjust.method = "fdr")

table(tt$P.Value < 0.005)
deg_inf24 <- tt[tt$adj.P.Val <= 0.05, ]
dim(deg_inf24)

deg_inf24 <- deg_inf24[, c(3,18,22)] # selecting only Entrez, log2FC and p-value columns.
head(deg_inf24)

deg_inf24 <- deg_inf24[!duplicated(deg_inf24$GENE),]
dim(deg_inf24)
table(is.na(deg_inf24$GENE))

deg_inf24 <- deg_inf24[!is.na(deg_inf24$GENE),]
dim(deg_inf24)

# Creating geneList object

## feature 1: numeric vector
head(deg_inf24)

geneList_inf24 = deg_inf24[, 2]

## feature 2: named vector
names(geneList_inf24) = as.character(deg_inf24[, 1])

## feature 3: decreasing order
geneList_inf24 = sort(geneList_inf24, decreasing = TRUE)
head(geneList_inf24)

## Getting changed genes for infection 48 h -----------

tt <- topTable(fit2, coef = 2, number = Inf, genelist = fit$genes, adjust.method = "fdr")

table(tt$P.Value < 0.005)
deg_inf48 <- tt[tt$adj.P.Val <= 0.05, ]
dim(deg_inf48)

deg_inf48 <- deg_inf48[, c(3,18,22)] # selecting only Entrez, log2FC and p-value columns.
head(deg_inf48)

deg_inf48 <- deg_inf48[!duplicated(deg_inf48$GENE),]
dim(deg_inf48)
table(is.na(deg_inf48$GENE))

deg_inf48 <- deg_inf48[!is.na(deg_inf48$GENE),]
dim(deg_inf48)

# Creating geneList object

## feature 1: numeric vector
head(deg_inf48)

geneList_inf48 = deg_inf48[, 2]

## feature 2: named vector
names(geneList_inf48) = as.character(deg_inf48[, 1])

## feature 3: decreasing order
geneList_inf48 = sort(geneList_inf48, decreasing = TRUE)
head(geneList_inf48)

##### GO Enrichment Analysis Infection 24 h --------------------------

## Cellular Component 24h ------------------

go.cc.24 <-  enrichGO(gene = names(geneList_inf24),
                   universe      = universe,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "CC",
                   pAdjustMethod = "fdr",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

go.cc.24@result[,1:7]


pdf("Figs_pdf/CellComponent_Inf24h.pdf", width = 10, pagecentre = TRUE, onefile = TRUE)

barplot(go.cc.24, 
        colorBy = "p.adjust",
        title = "Gene Ontology over-representation test for Cellular Component 24 h Infection",
        showCategory = length(go.cc.24@result$ID), font.size = 8)
dev.off()

cnetplot(go.cc.24, categorySize="pvalue", foldChange=geneList_inf24, fixed = FALSE)

write.csv(go.cc.24@result[,1:7], file = "results/go.cc.24.csv", row.names = FALSE)

## Biological Process 24h ------------------

go.bp.24 <-  enrichGO(gene = names(geneList_inf24),
                      universe      = universe,
                      OrgDb         = org.Hs.eg.db,
                      ont           = "BP",
                      pAdjustMethod = "fdr",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)

go.bp.24@result[,1:7]

go.bp.24 <- simplify(go.bp.24, cutoff = 0.5, by = "p.adjust", select_fun = min)

pdf("Figs_pdf/BioProcess_Inf24h.pdf", width = 10, pagecentre = TRUE, onefile = TRUE)

barplot(go.bp.24, 
        colorBy = "p.adjust",
        title = "Gene Ontology over-representation test for Biological Process 24 h Infection",
        showCategory = length(go.bp.24@result$ID), font.size = 8)
dev.off()

cnetplot(go.bp.24, categorySize="pvalue", foldChange=geneList_inf24, fixed = FALSE)

write.csv(go.bp.24@result[,1:7], file = "results/go.bp.24.csv", row.names = FALSE)


## Molecular Function24h ------------------

go.mf.24 <-  enrichGO(gene = names(geneList_inf24),
                      universe      = universe,
                      OrgDb         = org.Hs.eg.db,
                      ont           = "MF",
                      pAdjustMethod = "fdr",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)

go.mf.24@result[,1:7]

pdf("Figs_pdf/MolFun_Inf24h.pdf", width = 10, pagecentre = TRUE, onefile = TRUE)

barplot(go.mf.24, 
        colorBy = "p.adjust",
        title = "Gene Ontology over-representation test for Molecular Function 24 h Infection",
        showCategory = length(go.mf.24@result$ID), font.size = 8)
dev.off()

cnetplot(go.mf.24, categorySize="pvalue", foldChange=geneList_inf24, fixed = FALSE)

write.csv(go.mf.24@result[,1:7], file = "results/go.mf.24.csv", row.names = FALSE)

##### GO Enrichment Analysis Infection 48 h --------------------------

## Cellular Component 48h ------------------

go.cc.48 <- enrichGO(gene = names(geneList_inf48),
                      universe      = universe,
                      OrgDb         = org.Hs.eg.db,
                      ont           = "CC",
                      pAdjustMethod = "fdr",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)

go.cc.48@result[,1:7]

pdf("Figs_pdf/CellComponent_Inf48h.pdf", width = 10, pagecentre = TRUE, onefile = TRUE)

barplot(go.cc.48, 
        colorBy = "p.adjust",
        title = "Gene Ontology over-representation test for Cellular Component 48 h Infection",
        showCategory = length(go.cc.48@result$ID), font.size = 8)
dev.off()

cnetplot(go.cc.48, categorySize="pvalue", foldChange=geneList_inf48, fixed = FALSE)

write.csv(go.cc.48@result[,1:7], file = "results/go.cc.48.csv", row.names = FALSE)

## Biological Process 48 h ------------------

go.bp.48 <-  enrichGO(gene = names(geneList_inf48),
                      universe      = universe,
                      OrgDb         = org.Hs.eg.db,
                      ont           = "BP",
                      pAdjustMethod = "fdr",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)

go.bp.48@result[,1:7]

go.bp.48 <- simplify(go.bp.48, cutoff = 0.5, by = "p.adjust", select_fun = min)

pdf("Figs_pdf/BioProcess_Inf48h.pdf", width = 10, pagecentre = TRUE, onefile = TRUE)

barplot(go.bp.48, 
        colorBy = "p.adjust",
        title = "Gene Ontology over-representation test for Biological Process 48 h Infection",
        showCategory = length(go.bp.48@result$ID), font.size = 8)
dev.off()

cnetplot(go.bp.48, categorySize="pvalue", foldChange=geneList_inf48, fixed = FALSE)

write.csv(go.bp.48@result[,1:7], file = "results/go.bp.48.csv", row.names = FALSE)

## Molecular Function 48h ------------------

go.mf.48 <-  enrichGO(gene = names(geneList_inf48),
                      universe      = universe,
                      OrgDb         = org.Hs.eg.db,
                      ont           = "MF",
                      pAdjustMethod = "fdr",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)

go.mf.48@result[,1:7]

pdf("Figs_pdf/MolFun_Inf48h.pdf", width = 10, pagecentre = TRUE, onefile = TRUE)

barplot(go.mf.48, 
        colorBy = "p.adjust",
        title = "Gene Ontology over-representation test for Molecular Function 48 h Infection",
        showCategory = length(go.mf.48@result$ID), font.size = 8)
dev.off()

cnetplot(go.mf.48, categorySize="pvalue", foldChange=geneList_inf48, fixed = FALSE)

write.csv(go.mf.48@result[,1:7], file = "results/go.mf.48.csv", row.names = FALSE)

## END --------------
rm(list=ls())
q(save = "no")
