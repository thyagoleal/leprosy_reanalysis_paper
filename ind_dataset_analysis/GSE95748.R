### Reanalysis GSE95748
## Thyago Leal Calvo - thyagoleal@yahoo.com
## September 5th, 2017.
####################################################

## Settings ---------------------------

root <- getwd()
set.seed(129)
options(digits = 3, scipen = 999)

## Loading packages 

require("limma")
require("Biobase")
require("mouse4302.db")
require("affy")
library("simpleaffy")
library("affyQCReport")
library("arrayQualityMetrics")

## Loading raw .CEL files -------------------------

celfiles <- list.files("raw/GSE95748_RAW/", pattern = ".CEL$", full.names = TRUE)

raw = read.affybatch(celfiles)
sampleNames(raw)
sampleNames(raw) = sub("\\.CEL$", "",sampleNames(raw))

### Quality assessment -----------

png("figs/BoxplotRaw.png")
boxplot(raw, main="Raw Boxplot",
        outline = TRUE, col="lightblue")
dev.off()

quality <- qc(raw)

pdf("figs/QCStatistics.pdf")
plot(quality)
dev.off()

## Preprocessing and normalization ---------------------

norm.data <- rma(raw)

png("figs/BoxplotRMA.png")
boxplot(norm.data)
dev.off()

dim(norm.data)

## Fixing Phenotype - Samples Information -----------------

samples <- data.frame(samples = sampleNames(norm.data), Group = c("Control", "Control", "Infected_14d", "Infected_14d", "Infected_28", "Infected_28", "pSLCs", "pSLCs"))

samples

## Fixing Annotation --------------------

annot <- select(x = mouse4302.db, keys = row.names(norm.data), columns = c("ENTREZID", "SYMBOL", "GENENAME"))

annot <- annot[!duplicated(annot$PROBEID), ]
row.names(annot) <- annot$PROBEID

annot <- annot[row.names(norm.data),]

# Check
all(row.names(norm.data) == row.names(annot))

annot <- new("AnnotatedDataFrame", data =  annot)

norm.data@featureData <- annot

saveRDS(norm.data, "objs/eSet.rds", compress = "bzip2")

## Performing t-test with Limma ------------------

f <- factor(samples$Group)
design <- model.matrix(~0+f)
colnames(design) <- c("Control", "Infected14", "Infected28", "pSLCs")

fit <- lmFit(norm.data, design)

## Removing duplicated ENTREZID 

fit <- fit[order(fit$Amean, decreasing = TRUE), ]
fit <- fit[!duplicated(fit$genes$ENTREZID), ]

contrast.matrix <-
        makeContrasts(
                Infected14 - Control,
                Infected28 - Control,
                (Infected14 + Infected28) / 2 - Control,
                levels = design
        )

fit <- contrasts.fit(fit, contrast.matrix)

fit2 <- eBayes(fit)

tt <- topTable(fit2, coef ="Infected14 - Control", adjust.method = "fdr", number = Inf, confint = TRUE)

tt2 <- topTable(fit2, coef ="Infected28 - Control", adjust.method = "fdr", number = Inf, confint = TRUE)

tt3 <- topTable(fit2, coef ="(Infected14 + Infected28)/2 - Control", adjust.method = "fdr", number = Inf, confint = TRUE)

## Saving results ------------------------

require(WriteXLS)

WriteXLS(tt, ExcelFileName = "results/Results_14days.xls")
WriteXLS(tt2, ExcelFileName = "results/Results_28days.xls")
WriteXLS(tt3, ExcelFileName = "results/Results_MeanTime.xls")

# Exporting expression matrix -----------------

exp_matrix <- exprs(norm.data)
head(exp_matrix)
colnames(exp_matrix) <- samples$Group
exp_matrix <- exp_matrix[row.names(fit2),]

head(exp_matrix)

tmp <- merge.data.frame(exp_matrix, annot@data, by.x = "row.names", by.y = "PROBEID", all.x = TRUE)

head(tmp)

exp_matrix <- tmp[, 2:10]
head(exp_matrix)
rm(tmp)

write.table(exp_matrix, file = "results/exp_matrix_GSE95748.txt", sep = "\t", col.names = TRUE)

## Converting ids to human homologs --------------------------------

library(biomaRt)

mus <- as.character(unique(annot$ENTREZID))

human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

homolog_human <- getLDS(
        attributes = "entrezgene",
        filters = "entrezgene",
        values = mus,
        mart = mouse,
        attributesL = c("entrezgene", "hgnc_symbol"),
        martL = human
)

tt <- topTable(fit2, coef ="Infected14 - Control", adjust.method = "fdr", number = Inf, confint = TRUE)

results_infected14 <- merge.data.frame(tt, homolog_human, by.x = "ENTREZID", by.y = "NCBI.gene.ID")

results_infected14 <- results_infected14[, c(2:15)]

names(results_infected14)[13] <- "ENTREZID"

o <- order(results_infected14$adj.P.Val, decreasing = FALSE)

WriteXLS(results_infected14, ExcelFileName = "results/Results_GSE95748_inf14d_human_ids.xls")

tt2 <- topTable(fit2, coef ="Infected28 - Control", adjust.method = "fdr", number = Inf, confint = TRUE)

results_infected28 <- merge.data.frame(tt2, homolog_human, by.x = "ENTREZID", by.y = "NCBI.gene.ID")

results_infected28 <- results_infected28[, c(2:14)]
names(results_infected28)[12] <- "ENTREZID"
results_infected28 <- results_infected28[, c(12:13, 1:10)]
WriteXLS(results_infected28, ExcelFileName = "results/Results_GSE95748_inf28d_human_ids.xls")

## Venn Diagram ---------------------

require(limma)

results <- decideTests(fit2, adjust.method = "fdr", p.value = 0.05)

pdf("figs/VennDiagram.pdf", paper = "a4", pagecentre = TRUE)

vennDiagram(results, include = c("up", "down"), 
            counts.col = c("firebrick1", "dodgerblue2"), 
            circle.col = c("limegreen", "yellow", "dodgerblue2", "firebrick1"), cex = 1)

dev.off()

### Volcano plot =================

library(RColorBrewer)
library(gplots)
library(ggplot2)
library(limma)

## Infection 14 days  vs Control -------------------------------

coef1 <- topTable(fit2,
                  coef = 1,
                  adjust.method = "fdr",
                  number = Inf)

Significant = coef1$P.Value < 0.05 & abs(coef1$logFC) > 1

names(coef1)

plot_coef1 <- ggplot(data = coef1, 
                     aes(x = logFC,
                         y = -log10(P.Value), 
                         colour = Significant)) +
        geom_point(alpha = 1, size = .5) +
        xlab("log2 Fold Change") + ylab("-log10 P-value") + theme_bw() +
        theme(legend.position = "none") + scale_color_manual(values=c("lightgrey", "black"))

plot_coef1 <- plot_coef1 + ggtitle("Differentially Expressed Genes from GSE95748- Infection 14 days vs Control") + geom_hline(yintercept = -log10(0.05), colour = "dodgerblue", linetype = 3) + geom_vline(xintercept = c(-1,1), colour = c("green", "red"), linetype = 3)

pdf("figs/VolcanoPlots.pdf", onefile = TRUE, paper = "a4r", pagecentre = TRUE)
plot_coef1


## Infection 28 days  vs Control -------------------------------

coef2 <- topTable(fit2,
                  coef = 2,
                  adjust.method = "fdr",
                  number = Inf)

Significant = coef2$P.Value < 0.05 & abs(coef2$logFC) > 1

names(coef2)

plot_coef2 <- ggplot(data = coef2, 
                     aes(x = logFC,
                         y = -log10(P.Value), 
                         colour = Significant)) +
        geom_point(alpha = 1, size = .5) +
        xlab("log2 Fold Change") + ylab("-log10 P-value") + theme_bw() +
        theme(legend.position = "none") + scale_color_manual(values=c("lightgrey", "black"))

plot_coef2 <- plot_coef2 + ggtitle("Differentially Expressed Genes from GSE95748- Infection 28 days vs Control") + geom_hline(yintercept = -log10(0.05), colour = "dodgerblue", linetype = 3) + geom_vline(xintercept = c(-1,1), colour = c("green", "red"), linetype = 3)

plot_coef2
dev.off()

## Heatmap ----------------------------

library(pheatmap)
library(RColorBrewer)

plot = tt$ENTREZID[tt3$adj.P.Val < 0.005]

label = data.frame(Group = samples$Group)
row.names(label) <- make.names(label$Group, unique = TRUE)

pdf("figs/heatmap_deg.pdf", paper = "a4", pagecentre = TRUE, onefile = FALSE)
pheatmap(exp_matrix[exp_matrix$ENTREZID %in% plot, 1:8],
         color = brewer.pal(8, "RdYlBu"),
         scale = "row",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method = "average",
         legend = TRUE,
         show_rownames = FALSE,
         main = paste("Heatmap of GSE95748 DEGs (", length(plot),  ") with Adj.P-value < 0.005", sep = ""),
         border_color = NA, 
         annotation_col = label,
         annotation_legend = TRUE,
         show_colnames = FALSE)
dev.off()

# Gene Ontology Enrichment Analysis  ------------------

### Loading required packages ==============

library(clusterProfiler)
library(org.Mm.eg.db)
library(pathview)

#### Getting Universe genes ==================

universe <- results_infected14$ENTREZID
# grab only gene entrez

## Getting changed genes for infection 14 days -----------

table(results_infected14$P.Value < 0.005)
deg_inf14 <- results_infected14[results_infected14$adj.P.Val <= 0.005, ]
dim(deg_inf14)

deg_inf14 <- deg_inf14[, c(9,2,6)] # selecting only Entrez, log2FC and p-value columns.
head(deg_inf14)

deg_inf14 <- deg_inf14[!duplicated(deg_inf14$ENTREZID),]
dim(deg_inf14)
table(is.na(deg_inf14$ENTREZID))

deg_inf14 <- deg_inf14[!is.na(deg_inf14$ENTREZID),]
dim(deg_inf14)

# Creating geneList object

## feature 1: numeric vector
head(deg_inf14)

geneList_inf14 = deg_inf14[, 2]

## feature 2: named vector
names(geneList_inf14) = as.character(deg_inf14[, 1])

## feature 3: decreasing order
geneList_inf14 = sort(geneList_inf14, decreasing = TRUE)
head(geneList_inf14)

##### GO Enrichment Analysis Infection 14 d --------------------------

## Cellular Component 14 d ------------------

go.cc.14 <-  enrichGO(gene = names(geneList_inf14),
                      universe      = universe,
                      OrgDb         = org.Mm.eg.db,
                      ont           = "CC",
                      pAdjustMethod = "fdr",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)

go.cc.14@result[,1:7]

go.cc.14 <- simplify(go.cc.24, cutoff = 0.7, by = "p.adjust", select_fun = min)

pdf("figs_pdf/CellComponent_inf14d.pdf", width = 10, pagecentre = TRUE, onefile = TRUE)

barplot(go.cc.14, 
        colorBy = "p.adjust",
        title = "Gene Ontology over-representation test for Cellular Component 14 d Infection",
        showCategory = length(go.cc.24@result$ID), font.size = 8)
dev.off()

cnetplot(go.cc.14, categorySize="pvalue", foldChange=geneList_inf14, fixed = FALSE)

write.csv(go.cc.14@result[,1:7], file = "results/go.cc.14.csv", row.names = FALSE)

## Biological Process 14h ------------------

go.bp.14 <-  enrichGO(gene = names(geneList_inf14),
                      universe      = universe,
                      OrgDb         = org.Mm.eg.db,
                      ont           = "BP",
                      pAdjustMethod = "fdr",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)

go.bp.14@result[,1:7]

go.bp.14 <- simplify(go.bp.14, cutoff = 0.5, by = "p.adjust", select_fun = min)

pdf("figs_pdf/BioProcess_inf14d.pdf", width = 10, pagecentre = TRUE, onefile = TRUE)


barplot(
        go.bp.14,
        colorBy = "p.adjust",
        title = "Gene Ontology over-representation test for Biological Process 14 d Infection",
        showCategory = length(go.bp.14@result$ID),
        font.size = 8
)
dev.off()

cnetplot(
        go.bp.14,
        categorySize = "pvalue",
        foldChange = geneList_inf14,
        fixed = FALSE
)

write.csv(go.bp.14@result[, 1:7], file = "results/go.bp.14.csv", row.names = FALSE)

## Molecular Function 14h ------------------

go.mf.14 <-  enrichGO(gene = names(geneList_inf14),
                      universe      = universe,
                      OrgDb         = org.Mm.eg.db,
                      ont           = "MF",
                      pAdjustMethod = "fdr",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)

go.mf.14@result[,1:7]

go.mf.14 <-
        simplify(go.mf.14,
                 cutoff = 0.7,
                 by = "p.adjust",
                 select_fun = min)

pdf(
        "figs_pdf/MolFunction_inf14d.pdf",
        width = 10,
        pagecentre = TRUE,
        onefile = TRUE
)

barplot(
        go.mf.14,
        colorBy = "p.adjust",
        title = "Gene Ontology over-representation test for Molecular Function 14 d Infection",
        showCategory = length(go.mf.14@result$ID),
        font.size = 8
)
dev.off()

cnetplot(
        go.mf.14,
        categorySize = "pvalue",
        foldChange = geneList_inf14,
        fixed = FALSE
)

write.csv(go.mf.14@result[, 1:7], file = "results/go.mf.14.csv", row.names = FALSE)




























