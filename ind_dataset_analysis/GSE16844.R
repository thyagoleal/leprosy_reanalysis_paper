## GSE16844 - neutrophil recruitment 
## Re-analysis of affymetrix one-channel array - 2017-03-02  Thyago Leal 
## Reviewed: July, 21st 207

## Default code -----------

options(digits = 3)
root <- getwd()
set.seed(2107)
dir()

## Loading libraries ---------------------

library("affy")
library("hgu133plus2.db")
library("limma")
library("simpleaffy")
library("affyQCReport")
library("arrayQualityMetrics")
library("hgu133plus2.db")
library("AnnotationDbi")
library("WriteXLS")
library("pheatmap")
library("RColorBrewer")

## Importing celfiles ---------------------

celfiles = list.files("raw_files/", pattern = ".CEL", full.names = TRUE)

raw = read.affybatch(celfiles)
sampleNames(raw)
sampleNames(raw) = sub("\\.CEL$", "",sampleNames(raw))

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

png("figs/BoxplotRAW.png")
boxplot(raw, col = "lightblue", outline = F, ylab = "Log2 Intensidade", main = "Distribuição das intensidades antes da normalização", xlab = "Arranjos")
dev.off()

norm.data <- rma(raw)

png("figs/BoxplotRMA.png")
boxplot(exprs(norm.data), col = "lightblue", outline = F, ylab = "Log2 Intensidade", main = "Distribuição das intensidades após normalização", xlab = "Arranjos")
dev.off()

dim(norm.data)

## Fixing Phenotype - Samples Information -----------------

samples <- data.frame(samples = sampleNames(norm.data), Classification = c(rep("ENL",2), rep("LL",6), rep("ENL", 4), "LL"))

row.names(samples) <- samples$samples

## Annotating and creating new eset -----------

## Annotation 

probes <- row.names(norm.data)

annt <- select(hgu133plus2.db, probes, columns = c("ENTREZID", "SYMBOL", "GENENAME"))

annt <- annt[match(probes, annt$PROBEID),]

row.names(annt) <- annt$PROBEID

featuredata <- new("AnnotatedDataFrame",
                                 data=annt)
phenoData <- new("AnnotatedDataFrame",
                 data=samples)

eset <- ExpressionSet(assayData=exprs(norm.data),
                            phenoData=phenoData, 
                      featureData = featuredata,
                            annotation="hgu133plus2.db")


ifelse(all(row.names(eset) == eset@featureData@data$PROBEID),
       "All Set",
       "Something IS WRONG, mate!")

saveRDS(eset, file = "objs/eSet.Rds", compress = "bzip2")

## Performing t-test with Limma ------------------

f <- factor(samples$Classification, levels = c("ENL", "LL"))
design <- model.matrix( ~ 0 + f)
colnames(design) <- c("ENL", "LL")
fit <- lmFit(eset, design)

## Filtering out duplicated probes ---------------

o <- order(fit$Amean, decreasing = TRUE)
fit <- fit[o,]
d <- duplicated(fit$genes$ENTREZID)

fit <- fit[!d,]

contrast.matrix <- makeContrasts(LL-ENL,
                                 levels=design)

fit <- contrasts.fit(fit, contrast.matrix)

fit2 <- eBayes(fit)

results_LLvsENL <- topTable(fit2, adjust.method = "fdr", number = Inf, confint = TRUE)
head(results_LLvsENL)

## counting the Adjusted P-values obtained (just for fun) ---------

pvals <- c(0.05, 0.005, 0.0005, 0.0001, 0.00001)

for(p in pvals){
        print(paste(
                "Number of DEG with Adj.P-Val less than or equal to ", p,": ",
                length(which(results_LLvsENL$adj.P.Val <= p)),
                sep = ""))
}

rm(p, pvals)

## Saving results -----------------------------

WriteXLS(results_LLvsENL, ExcelFileName = "results/results_GSE16844_LLvsENL.xls", verbose = TRUE)

exp_matrix <- as.data.frame(exprs(norm.data))
exp_matrix$probes <- row.names(exp_matrix)
rownames(exp_matrix)

tmp <- select(
                hgu133plus2.db,
                keys = row.names(exp_matrix),
                columns = c("ENTREZID", "SYMBOL")
        )

head(exp_matrix)

tmp <- tmp[match(row.names(exp_matrix), tmp$PROBEID),]
exp_matrix <- merge.data.frame(tmp, exp_matrix, by.x = "PROBEID", by.y = 0, all.x = TRUE, sort = F)
all(exp_matrix$probes == exp_matrix$PROBEID)
rm(tmp)
row.names(exp_matrix) <- exp_matrix$probes
o <- order(rowMeans(exp_matrix[,4:16]), decreasing = TRUE)
exp_matrix <- exp_matrix[o, ]
exp_matrix <- exp_matrix[!duplicated(exp_matrix$ENTREZID), ]
exp_matrix <- exp_matrix[!is.na(exp_matrix$ENTREZID), ]
colnames(exp_matrix)[4:16] <- make.names(samples$Classification, unique = T)

exp_matrix <- exp_matrix[, 1:16]

write.table(
        exp_matrix,
        file = "results/exp_matrix_gse16844.txt",
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
)


### Volcano plot ------------------

library(RColorBrewer)
library(gplots)
library(ggplot2)

significant <-
        results_LLvsENL$P.Value < 0.05 & abs(results_LLvsENL$logFC) > 1
table(significant)


names(results_LLvsENL)

plot1 <- ggplot(data = results_LLvsENL,
                aes(
                        x = logFC,
                        y = -log10(P.Value),
                        colour = significant
                )) +
        geom_point(alpha = 1, size = .8) +
        xlab("log2 Fold Change") + ylab("-log10 P-value") + theme_bw() +
        theme(legend.position = "none") + scale_color_manual(values = c("lightgrey", "black"))+ geom_hline(yintercept = -log10(0.05), colour = "dodgerblue", linetype = 3) + geom_vline(xintercept = c(-1,1), colour = c("green", "red"), linetype = 3) 

plot1 <-
        plot1 + ggtitle("Differentially Expressed Genes (1314) GSE16844 - LL vs ENL")

pdf(
        "figs/VolcanoPlot_LLvsENL.pdf",
        onefile = TRUE,
        paper = "a4r",
        pagecentre = TRUE
)
plot1
dev.off()
rm(plot1, significant)

# Heatmap ----------------------------

plot <- results_LLvsENL$ENTREZID[results_LLvsENL$adj.P.Val < 0.05 & abs(results_LLvsENL$logFC) > 1]

label = data.frame(Class = samples$Classification)
row.names(label) <- make.names(samples$Classification, unique = T)
ann_colors = list(Class = c(ENL = "limegreen",
                            LL = "firebrick1"))

pdf(
        "figs/heatmap_deg.pdf",
        paper = "a4",
        onefile = FALSE
)
pheatmap(
        exp_matrix[exp_matrix$ENTREZID %in% plot, 4:16],
        color = brewer.pal(11, "RdYlBu"),
        scale = "row",
        clustering_distance_rows = "correlation",
        clustering_distance_cols = "euclidean",
        clustering_method = "average",
        legend = TRUE,
        show_rownames = FALSE,
        main = paste("Heatmap of GSE17763 DEGs (", length(unique(plot)),") with Adj.P-value < 0.05 & |log2FC| >1", sep = ""),
        border_color = NA,
        annotation_col = label,
        annotation_legend = TRUE,
        annotation_colors = ann_colors,
        show_colnames = FALSE
)
dev.off()

## PCA with all the DEG and all samples -----------------------------------------

require(FactoMineR)
require(factoextra)

pca.all <- PCA(t(exp_matrix[, 4:16]),  scale.unit = TRUE)

pdf("figs/PCA_all_DEG_samples.pdf", onefile = FALSE)
fviz_pca_ind(
        pca.all,
        geom = "point",
        habillage = as.factor(samples$Classification),
        pointsize = 3,
        title = "PCA from All DEG (5028) and all samples"
) + theme_bw()
dev.off()

###  Gene Ontology Enrichment Analysis  -----------------------

# LL vs ENL 

universe <- results_LLvsENL$ENTREZID

deg <- results_LLvsENL[results_LLvsENL$adj.P.Val < 0.05, ]
deg <- deg[, c("logFC", "ENTREZID")]
deg1 <- deg$logFC
names(deg1) <- as.character(deg$ENTREZID)
deg1 <- deg1[!duplicated(names(deg1))]
deg1 <- sort(deg1, decreasing = TRUE) ## geneList with Log2FC + Entrez for DEG with P-val < 0.05

library(clusterProfiler)
library(org.Hs.eg.db)
library(pathview)

## LL vs ENL Cell component ---

deg1.cc <-  enrichGO(gene = names(deg1),
                     universe      = universe,
                     OrgDb         = org.Hs.eg.db,
                     ont           = "CC",
                     pAdjustMethod = "fdr",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)
head(deg1.cc)

pdf("figs_pdf/CellComponent_LLvsENL.pdf", width = 10, pagecentre = TRUE, onefile = TRUE)

barplot(
        deg1.cc,
        colorBy = "p.adjust",
        title = "Ontologias Gênicas (GO) super-representadas para Componentes Celulares (LL vs ENL)",
        showCategory = length(deg1.cc@result$ID),
        font.size = 12
)

dev.off()
cnetplot(deg1.cc, showCategory = 10, foldChange = deg1, fixed = FALSE)

write.csv(deg1.cc@result[, 1:7], file = "results/GO_cc_LLvsENL.csv", row.names = FALSE)

## LL vs ENL Biological Process  ---

deg1.bp <- enrichGO(gene = names(deg1),
                    universe      = universe,
                    OrgDb         = org.Hs.eg.db,
                    ont           = "BP",
                    pAdjustMethod = "fdr",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)
head(deg1.bp)

pdf("figs_pdf/BioProcess_LLvsENL.pdf", width = 10, pagecentre = TRUE, onefile = TRUE)

barplot(deg1.bp, 
        colorBy = "p.adjust",
        title = "Ontologias Gênicas (GO) super-representadas para Processos Biológicos (LL vs ENL",
        showCategory = length(deg1.bp@result$ID),
        font.size = 12)

dev.off()

cnetplot(deg1.bp, showCategory = 4, foldChange = deg1, fixed = FALSE)

write.csv(deg1.bp@result[,1:7], file = "results/GO_BP_LLvsENL.csv", row.names = FALSE)


## LL vs ENL Molecular Function  ---

deg1.mf <- enrichGO(
        gene = names(deg1),
        universe      = universe,
        OrgDb         = org.Hs.eg.db,
        ont           = "MF",
        pAdjustMethod = "fdr",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.05,
        readable      = TRUE
)
head(deg1.mf)

pdf("figs_pdf/MolFunction_LLvsENL.pdf", width = 10, pagecentre = TRUE, onefile = TRUE)

barplot(deg1.mf, 
        colorBy = "p.adjust",
        title = "Ontologias Gênicas (GO) super-representadas para Funções Moleculares (LL vs ENL)",
        showCategory = length(deg1.mf@result$ID),
        font.size = 12)

dev.off()

write.csv(deg1.mf@result[,1:7], file = "results/GO_MF_LLvsENL.csv", row.names = FALSE)

### END ### 
rm(list=ls())
q(save = "no")
