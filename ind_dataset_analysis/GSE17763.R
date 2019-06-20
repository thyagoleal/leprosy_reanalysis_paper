## GSE17763 - Divergence of macrophage phagocytic and antimicrobial programs in leprosy
## Re-analysis of affymetrix one-channel array - 2017-03-01  Thyago Leal 
## Reviewed: July, 21st 2017.

## Default code ------

root <- getwd()
set.seed(139)
options(digits = 3)

## Loading libraries -------

library("GEOquery")
library("affy")
library("affyiof")
library("simpleaffy")
library("affyQCReport")
library("arrayQualityMetrics")
library("limma")
library("WriteXLS")
library("hgu133plus2.db")

### Importing annotation from GEO -------

gset <- getGEO("GSE17763", destdir = "raw_files/")

## Importing celfiles ---

celfiles = list.files("raw_files/GSE17763_RAW/", pattern = ".CEL", full.names = TRUE)

raw = read.affybatch(celfiles)
sampleNames(raw)

sampleNames(raw) = sub("\\.CEL$", "",sampleNames(raw))

## Importing the file describind the samples (created manually)

samples <- as.data.frame(read.table("raw_files/sampleinformation", sep = "\t", header = T))
row.names(samples) <- samples$samples
row.names(samples) == sampleNames(raw)

pData(raw) <- samples
pData(raw)

### Quality assessment -----------

png("figs/BoxplotRaw.png")
boxplot(raw, main="Raw Boxplot",
                    outline = TRUE, col="lightblue")
dev.off()

quality <- qc(raw)

pdf("figs_pdf/QCStatistics.pdf")
plot(quality)
dev.off()

## Preprocessing and normalization ---------------------

norm.data <- rma(raw)

png("figs/BoxplotRMA.png")
boxplot(norm.data)
dev.off()

dim(norm.data)

## Annotating -----------

probes <- row.names(norm.data)

annt <-  select(hgu133plus2.db, probes, columns = c("ENTREZID", "SYMBOL", "GENENAME"))
dim(norm.data)

annt <- annt[!duplicated(annt$PROBEID),]
row.names(annt) <- annt$PROBEID

phenodata = new("AnnotatedDataFrame",
                data = samples)

featuredata <- new("AnnotatedDataFrame",
                   data = annt)
norm.data <-
        ExpressionSet(
                assayData = exprs(norm.data),
                featureData = featuredata,
                phenoData = phenodata,
                annotation = "hgu95av2"
        )

saveRDS(norm.data, file = "objs/eSet.Rds", compress = "bzip2")

## Performing t-test with Limma ------------------

f <- factor(samples$Classification, levels = c("BT", "LL", "RR"))
design <- model.matrix( ~ 0 + f)
colnames(design) <- c("BT", "LL", "RR")
fit <- lmFit(norm.data, design)

contrast.matrix <- makeContrasts(LL - BT, LL - RR, BT - RR,
                                 levels = design)

## Removing duplicated by leaving the one with maximum expressions across samples -------

o <- order(fit$Amean, decreasing = TRUE)
fit <- fit[o, ]
d <- duplicated(fit$genes$ENTREZID)
fit <- fit[!d, ]

fit <- contrasts.fit(fit, contrast.matrix)

fit2 <- eBayes(fit)

results_LLvsBT <-
        topTable(fit2,
                 coef = 1,
                 adjust.method = "fdr",
                 number = Inf, 
                 confint = TRUE)
results_LLvsRR  <-
        topTable(fit2,
                 coef = 2,
                 adjust.method = "fdr",
                 number = Inf,
                 confint = TRUE)
results_BTvsRR <-
        topTable(fit2,
                 coef = 3,
                 adjust.method = "fdr",
                 number = Inf, 
                 confint = TRUE)

print("Does probe match?")
ifelse(all(row.names(results_LLvsBT) == results_LLvsBT$PROBEID),
       "Yes, probes match!",
       "Something is wrong!")

## Saving results -----------------------------

list <- list("results_LLvsBT", "results_LLvsRR", "results_BTvsRR")

for(i in list){
        filename <- paste("results/",i,".xls", sep = "")
        WriteXLS(x = i, ExcelFileName = filename)
}
rm(i, filename)

## Expression Matrix -------------------------

exp_matrix <- exprs(norm.data)
head(exp_matrix)
colnames(exp_matrix) <- samples$Classification
exp_matrix <- as.data.frame(exp_matrix)
annotation <-  select(hgu133plus2.db, keys = row.names(exp_matrix), columns = c("ENTREZID", "SYMBOL"))
annotation <- annotation[!duplicated(annotation$PROBEID), ]

exp_matrix <- merge.data.frame(exp_matrix, annotation, by.x = 'row.names', by.y = "PROBEID")

o <- order(rowMeans(exp_matrix[,2:25]), decreasing = T)
exp_matrix <- exp_matrix[o, ] 
exp_matrix <- exp_matrix[!duplicated(exp_matrix$ENTREZID), ]
exp_matrix <- exp_matrix[!is.na(exp_matrix$ENTREZID),]

write.table(exp_matrix, file = "results/exp_matrix_gse17763.txt", quote = F, sep = "\t", row.names = FALSE)

### Volcano plot =================

library(RColorBrewer)
library(gplots)
library(ggplot2)
library(limma)

## LL vs BT -------

sig_1 = results_LLvsBT$P.Value < 0.05 & abs(results_LLvsBT$logFC) > 1

names(results_LLvsBT)

plot_coef1 <- ggplot(data = results_LLvsBT, 
                     aes(x = logFC,
                         y = -log10(P.Value), 
                         colour = sig_1)) +
        geom_point(alpha = 1, size = .8) +
        xlab("log2 Fold Change") + ylab("-log10 P-value") + theme_bw() +
        theme(legend.position = "none") + scale_color_manual(values=c("lightgrey", "black")) + geom_hline(yintercept = -log10(0.05), colour = "dodgerblue", linetype = 3) + geom_vline(xintercept = c(-1,1), colour = c("green", "red"), linetype = 3)

plot_coef1 <- plot_coef1 + ggtitle("Differentially Expressed Genes from GSE17763 - LL vs BT")

pdf(file = "figs/VolcanoPlots.pdf", onefile = TRUE)

plot_coef1


## LL vs RR -------

sig_2 = results_LLvsRR$P.Value < 0.05 & abs(results_LLvsRR$logFC) > 1

names(results_LLvsRR)

plot_coef2 <- ggplot(data = results_LLvsRR, 
                     aes(x = logFC,
                         y = -log10(P.Value), 
                         colour = sig_2)) +
        geom_point(alpha = 1, size = .8) +
        xlab("log2 Fold Change") + ylab("-log10 P-value") + theme_bw() +
        theme(legend.position = "none") + scale_color_manual(values=c("lightgrey", "black")) + geom_hline(yintercept = -log10(0.05), colour = "dodgerblue", linetype = 3) + geom_vline(xintercept = c(-1,1), colour = c("green", "red"), linetype = 3)

plot_coef2 <- plot_coef2 + ggtitle("Differentially Expressed Genes from GSE17763 - LL vs RR")

plot_coef2

## BT vs RR  -------

sig_3 = results_BTvsRR$P.Value < 0.05 & abs(results_BTvsRR$logFC) > 1

names(results_BTvsRR)

plot_coef3 <- ggplot(data = results_BTvsRR, 
                     aes(x = logFC,
                         y = -log10(P.Value), 
                         colour = sig_3)) +
        geom_point(alpha = 1, size = .8) +
        xlab("log2 Fold Change") + ylab("-log10 P-value") + theme_bw() +
        theme(legend.position = "none") + scale_color_manual(values=c("lightgrey", "black"))+ geom_hline(yintercept = -log10(0.05), colour = "dodgerblue", linetype = 3) + geom_vline(xintercept = c(-1,1), colour = c("green", "red"), linetype = 3)

plot_coef3 <- plot_coef3 + ggtitle("Differentially Expressed Genes from GSE17763 - BT vs RR")

plot_coef3
dev.off()

# Heatmap ----------------------------

library(pheatmap)
library(RColorBrewer)

plot <- topTable(fit2, number = Inf, coef = 1, adjust.method = "fdr")
plot <- plot$ENTREZID[plot$adj.P.Val < 0.05]

label = data.frame(Class = samples$Classification)
row.names(label) <- make.names(samples$Classification, unique = TRUE)
ann_colors = list(Class = c(BT = "limegreen",
                            LL = "firebrick1",
                            RR = "dodgerblue"))

pdf("figs/heatmap_deg.pdf", paper = "a4", pagecentre = TRUE, onefile = FALSE)
pheatmap(exp_matrix[rownames(exp_matrix) %in% plot,2:25],
         color = brewer.pal(9,"RdYlBu"),
         scale = "row",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "euclidean",
         clustering_method = "average",
         legend = TRUE,
         show_rownames = FALSE,
         main = paste("Heatmap of GSE17763 DEGs (", length(plot), ") with Adj.P-value < 0.05", sep = ""),
         border_color = NA, 
         annotation_col = label,
         annotation_legend = TRUE,
         annotation_colors = ann_colors, show_colnames = FALSE)
dev.off()

## Venn Diagram  ---------------------

require(limma)

res <- decideTests(fit2, adjust.method = "fdr", p.value = 0.05, lfc = 0)

pdf("figs/VennDiagram_GSE17763.pdf", paper = "a4", pagecentre = TRUE)

vennDiagram(res, include = c("up", "down"), 
            counts.col = c("firebrick1", "dodgerblue3"), 
            circle.col = c("firebrick1", "yellow", "dodgerblue3"),
            names = c("LL vs BT", "LL vs RR", "BT vs RR"))

dev.off()

###  Gene Ontology Enrichment Analysis  -----------------------

## Separate Contrasts

## Tissue LL vs BT ----------------------------

universe <- results_LLvsBT$ENTREZID

coef1 <- results_LLvsBT[results_LLvsBT$adj.P.Val < 0.05,]
deg1 <- coef1$logFC
names(deg1) <- as.character(coef1$ENTREZID)
deg1 <- deg1[!duplicated(names(deg1))]
deg1 <- sort(deg1, decreasing = TRUE) ## geneList with Log2FC + Entrez for DEG with P-val < 0.05

library(clusterProfiler)
library(org.Hs.eg.db)
library(pathview)

## LL vs BT Cell component ---

deg1.cc <-  enrichGO(gene = names(deg1),
                     universe      = universe,
                     OrgDb         = org.Hs.eg.db,
                     ont           = "CC",
                     pAdjustMethod = "fdr",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)
head(deg1.cc)

deg1.cc <- simplify(deg1.cc, cutoff = 0.5, by = "p.adjust", select_fun = min)

pdf("figs_pdf/CellComponent_LLvsBT.pdf", width = 10, pagecentre = TRUE, onefile = TRUE)

barplot(deg1.cc, 
        colorBy = "p.adjust",
        title = "Gene Ontology over-representation test for Cellular Component (LL vs BT)",
        showCategory = length(deg1.cc@result$ID),
        font.size = 12)

dev.off()
cnetplot(deg1.cc, showCategory = 10, foldChange = deg1, fixed = FALSE)

write.csv(deg1.cc@result[,1:7], file = "results/GO_cc_LLvsBT.csv", row.names = FALSE)

## LL vs BT Biological Process  ---

deg1.bp <- enrichGO(gene = names(deg1),
                     universe      = universe,
                     OrgDb         = org.Hs.eg.db,
                     ont           = "BP",
                     pAdjustMethod = "fdr",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)
head(deg1.bp)

deg1.bp <- simplify(deg1.bp, cutoff = 0.5, by = "p.adjust", select_fun = min)

pdf("figs_pdf/BioProcess_LLvsBT.pdf", width = 10, pagecentre = TRUE, onefile = TRUE)

barplot(deg1.bp, 
        colorBy = "p.adjust",
        title = "Gene Ontology over-representation test for Biological Process (LL vs BT)",
        showCategory = length(deg1.bp@result$ID),
        font.size = 12)

dev.off()

cnetplot(deg1.bp, showCategory = 4, foldChange = deg1, fixed = FALSE)

write.csv(deg1.bp@result[,1:7], file = "results/GO_BP_LLvsBT.csv", row.names = FALSE)


## LL vs BT Molecular Function  ---

deg1.mf <- enrichGO(gene = names(deg1),
                    universe      = universe,
                    OrgDb         = org.Hs.eg.db,
                    ont           = "MF",
                    pAdjustMethod = "fdr",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)
head(deg1.mf)

deg1.mf <- simplify(deg1.mf, cutoff = 0.5, by = "p.adjust", select_fun = min)

pdf("figs_pdf/MolFunction_LLvsBT.pdf", width = 10, pagecentre = TRUE, onefile = TRUE)

barplot(deg1.mf, 
        colorBy = "p.adjust",
        title = "Gene Ontology over-representation test for Molecular Function (LL vs BT)",
        showCategory = length(deg1.mf@result$ID),
        font.size = 12)

dev.off()

cnetplot(deg1.mf, showCategory = 10, foldChange = deg1, fixed = FALSE)

write.csv(deg1.mf@result[,1:7], file = "results/GO_MF_LLvsBT.csv", row.names = FALSE)




## Tissue LL vs RR ----------------------------

universe <- results_LLvsRR$ENTREZID

coef2 <- results_LLvsRR[results_LLvsRR$adj.P.Val < 0.05,]
deg2 <- coef2$logFC
names(deg2) <- as.character(coef2$ENTREZID)
deg2 <- deg2[!duplicated(names(deg2))]
deg2 <- sort(deg2, decreasing = TRUE) ## geneList with Log2FC + Entrez for DEG with P-val < 0.05


## LL vs RR Cell component ---

deg2.cc <-  enrichGO(gene = names(deg2),
                     universe      = universe,
                     OrgDb         = org.Hs.eg.db,
                     ont           = "CC",
                     pAdjustMethod = "fdr",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)
head(deg2.cc)

deg2.cc <- simplify(deg2.cc, cutoff = 0.5, by = "p.adjust", select_fun = min)

pdf("figs_pdf/CellComponent_LLvsRR.pdf", width = 10, pagecentre = TRUE, onefile = TRUE)

barplot(deg2.cc, 
        colorBy = "p.adjust",
        title = "Gene Ontology over-representation test for Cellular Component (LL vs RR)",
        showCategory = length(deg2.cc@result$ID),
        font.size = 12)

dev.off()
cnetplot(deg2.cc, showCategory = 10, foldChange = deg2, fixed = FALSE)

write.csv(deg2.cc@result[,1:7], file = "results/GO_cc_LLvsRR.csv", row.names = FALSE)

## LL vs RR Biological Process  ---

deg2.bp <- enrichGO(gene = names(deg2),
                    universe      = universe,
                    OrgDb         = org.Hs.eg.db,
                    ont           = "BP",
                    pAdjustMethod = "fdr",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)
head(deg2.bp)

deg2.bp <- simplify(deg2.bp, cutoff = 0.5, by = "p.adjust", select_fun = min)

pdf("figs_pdf/BioProcess_LLvsRR.pdf", width = 10, pagecentre = TRUE, onefile = TRUE)

barplot(deg2.bp, 
        colorBy = "p.adjust",
        title = "Gene Ontology over-representation test for Biological Process (LL vs RR)",
        showCategory = length(deg2.bp@result$ID),
        font.size = 12)

dev.off()

cnetplot(deg2.bp, showCategory = 10, foldChange = deg2, fixed = FALSE)

write.csv(deg2.bp@result[,1:7], file = "results/GO_BP_LLvsRR.csv", row.names = FALSE)


## LL vs RR Molecular Function  ---

deg2.mf <- enrichGO(gene = names(deg2),
                    universe      = universe,
                    OrgDb         = org.Hs.eg.db,
                    ont           = "MF",
                    pAdjustMethod = "fdr",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)
head(deg2.mf)

deg2.mf <- simplify(deg2.mf, cutoff = 0.5, by = "p.adjust", select_fun = min)

pdf("figs_pdf/MolFunction_LLvsRR.pdf", width = 10, pagecentre = TRUE, onefile = TRUE)

barplot(deg2.mf, 
        colorBy = "p.adjust",
        title = "Gene Ontology over-representation test for Molecular Function (LL vs RR)",
        showCategory = length(deg2.mf@result$ID),
        font.size = 12)

dev.off()

cnetplot(deg2.mf, showCategory = 10, foldChange = deg2, fixed = FALSE)

write.csv(deg2.mf@result[,1:7], file = "results/GO_MF_LLvsRR.csv", row.names = FALSE)

## Tissue BT vs RR ----------------------------

universe <- results_LLvsRR$ENTREZID

coef3 <- results_LLvsRR[results_LLvsRR$adj.P.Val < 0.05,]
deg3 <- coef3$logFC
names(deg3) <- as.character(coef3$ENTREZID)
deg3 <- deg3[!duplicated(names(deg3))]
deg3 <- sort(deg3, decreasing = TRUE) ## geneList with Log2FC + Entrez for DEG with P-val < 0.05


## BT vs RR Cell component ---

deg3.cc <-  enrichGO(gene = names(deg3),
                     universe      = universe,
                     OrgDb         = org.Hs.eg.db,
                     ont           = "CC",
                     pAdjustMethod = "fdr",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)
head(deg3.cc)

deg3.cc <- simplify(deg3.cc, cutoff = 0.5, by = "p.adjust", select_fun = min)

pdf("figs_pdf/CellComponent_LLvsRR.pdf", width = 10, pagecentre = TRUE, onefile = TRUE)

barplot(deg3.cc, 
        colorBy = "p.adjust",
        title = "Gene Ontology over-representation test for Cellular Component (BT vs RR)",
        showCategory = length(deg3.cc@result$ID),
        font.size = 12)

dev.off()
cnetplot(deg3.cc, showCategory = 10, foldChange = deg3, fixed = FALSE)

write.csv(deg3.cc@result[,1:7], file = "results/GO_cc_LLvsRR.csv", row.names = FALSE)

## BT vs RR Biological Process  ---

deg3.bp <- enrichGO(gene = names(deg3),
                    universe      = universe,
                    OrgDb         = org.Hs.eg.db,
                    ont           = "BP",
                    pAdjustMethod = "fdr",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)
head(deg3.bp)

deg3.bp <- simplify(deg3.bp, cutoff = 0.5, by = "p.adjust", select_fun = min)

pdf("figs_pdf/BioProcess_LLvsRR.pdf", width = 10, pagecentre = TRUE, onefile = TRUE)

barplot(deg3.bp, 
        colorBy = "p.adjust",
        title = "Gene Ontology over-representation test for Biological Process (BT vs RR)",
        showCategory = length(deg3.bp@result$ID),
        font.size = 12)

dev.off()

cnetplot(deg3.bp, showCategory = 10, foldChange = deg3, fixed = FALSE)

write.csv(deg3.bp@result[,1:7], file = "results/GO_BP_LLvsRR.csv", row.names = FALSE)


## BT vs RR Molecular Function  ---

deg3.mf <- enrichGO(gene = names(deg3),
                    universe      = universe,
                    OrgDb         = org.Hs.eg.db,
                    ont           = "MF",
                    pAdjustMethod = "fdr",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)
head(deg3.mf)

deg3.mf <- simplify(deg3.mf, cutoff = 0.5, by = "p.adjust", select_fun = min)

pdf("figs_pdf/MolFunction_LLvsRR.pdf", width = 10, pagecentre = TRUE, onefile = TRUE)

barplot(deg3.mf, 
        colorBy = "p.adjust",
        title = "Gene Ontology over-representation test for Molecular Function (BT vs RR)",
        showCategory = length(deg3.mf@result$ID),
        font.size = 12)

dev.off()

cnetplot(deg3.mf, showCategory = 10, foldChange = deg3, fixed = FALSE)

write.csv(deg3.mf@result[,1:7], file = "results/GO_MF_LLvsRR.csv", row.names = FALSE)

## Exporting All DEGs table ------

cont <-colnames(contrast.matrix)

for(i in cont){
        deg <- topTable(fit2, coef = i, genelist = fit$genes, number = Inf, adjust.method = "fdr", p.value = 0.05)
        deg <- deg[, c("ENTREZID", "logFC", "adj.P.Val")]
        deg$logFC <- with(deg, round(logFC, 3))
        write.table(x = deg, 
                    file = paste("results/", "DEG.", gsub(" ", "", i), ".txt", sep = ""),
                    quote = FALSE, 
                    sep = "\t", 
                    row.names = FALSE,
                    col.names = TRUE
                    )
}

rm(i,  deg)

### END #### 
rm(list=ls())
q(save = "no")