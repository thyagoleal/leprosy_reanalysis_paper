## Analysis of the GSE74481 data set
## Thyago Leal Calvo. May 25th, 2018
## thyagoleal@yahoo.com

## Settings ---------------

options(digits = 3)
set.seed(2396)
root <- getwd()

library(ggplot2)
library(GEOquery) 
library(limma)
library(convert)
library(geneplotter)
library(AnnotationDbi)
library(readxl)
library(HsAgilentDesign026652.db)
library(WriteXLS)
library(FactoMineR)
library(factoextra)

## Reading RAW Files   ------------------------------------------------

dir("raw_files/")

files <- list.files(path = "raw_files/", pattern = "^GSM.*.txt$", full.names = TRUE)
files

## Reading RAW Files   ------------------------------------------------

RG <- read.maimages(files, source = "agilent")
head(RG$genes)

# Checking file quality --------

pdf("figs/Diag_plot.pdf")
par(mfrow = c(1, 2))
imageplot(log2(RG$Rb[, 1]), RG$printer, low = "white", high = "red")
imageplot(log2(RG$Gb[, 1]), RG$printer, low = "white", high = "green")
dev.off()

pdf("figs/RG_channel_densities.pdf")
par(mfrow = c(1, 1))
plotDensities(RG)
dev.off()

postscript("figs/boxplot_before_norm.eps", width = 1000, height = 300)
boxplot(log2(RG$G), outline = FALSE, xlab = "Arranjos", ylab = "Intensidade (Log2)", main = "Gráfico de caixas com a distribuição da intensidade de cada arranjo antes da normalização", xaxt = "n", col = "lightblue")
axis(1, at = 1:ncol(RG), labels = 1:ncol(RG))
dev.off()

## Background correction --------

RG_bg <- backgroundCorrect(RG, method = "normexp")

## Assessing quality

# Before and after BG correction:
# before

pdf("figs/BG_correction.pdf", width = 1200, height = 800)
par(mfrow=c(1,2))
plotMD(RG, main = "RG Without BG Correction")

# after
plotMD(RG_bg, main =  "RG With BG Correction")
dev.off()

## Checking distribution of Red and Green Values 

pdf(file = "figs_pdf/Distribution_density_boxplot.pdf", paper = "a4r")
par(mfrow = c(2, 1))

plotformula_Green = log2(RG_bg$G) ~ col(RG_bg$G)
boxplot(
        plotformula_Green,
        outline = FALSE,
        col = "limegreen",
        xlab = "Arrays",
        ylab = expression(log[2] ~ G),
        main = "Boxplot of the Green Intensities from All Arrays"
)

plotformula_Red = log2(RG_bg$R) ~ col(RG_bg$R)
boxplot(
        plotformula_Red,
        outline = FALSE,
        col = "red2",
        xlab = "Arrays",
        ylab = expression(log[2] ~ R),
        main = "Boxplot of the Red Intensities from All Arrays"
)

dev.off()

pdf(file = "figs_pdf/Distribution_density.pdf", paper = "a4r")
par(mfrow = c(2, 1))

multidensity(plotformula_Green,
             main = "Green Densities",
             xlab = expression(log[2] ~ Green))

multidensity(plotformula_Red,
             main = "Red Densities", xlab = expression(log[2] ~ Red))

dev.off()

pdf(file = "figs_pdf/RG_densities.pdf", paper = "a4r")
par(mfrow = c(1, 1))
plotDensities(RG_bg, main = "Red and Green Densities combined")
dev.off()

## Normalization ---------------------

# 1 Robust Spline

MA <- normalizeWithinArrays(RG_bg, method = "robustspline")

# 2 Between Array Norm quantile method 

MA2 <- normalizeBetweenArrays(MA, method = "quantile")

pdf(
        file = "figs_pdf/Normalizatiom_MAplots.pdf",
        width = 1200,
        height = 600,
        paper = "a4r"
)
par(mfcol = c(1, 3), oma = c(1, 1, 1, 1))

limma::plotMA(RG_bg, main = "Before normalization")
limma::plotMA(MA, main = "After Robust spline Normalization")
limma::plotMA(MA2, main = "After Quantile Normalization")
dev.off()

postscript("figs/boxplot_after_norm.eps",
           width = 1000,
           height = 300)
boxplot(
        MA2$A,
        outline = FALSE,
        xlab = "Arranjos",
        ylab = "Intensidade (Log2)",
        main = "Gráfico de caixas com a distribuição da intensidade de cada arranjo após etapa de normalização",
        xaxt = "n",
        col = "lightblue"
)
axis(1, at = 1:ncol(RG), labels = 1:ncol(RG))
dev.off()

## Assembly of ExpressionSet ----------------------------

mdata <- as(MA2, "ExpressionSet")
class(mdata)

dim(exprs(mdata))
par(mfrow=c(1,1))
boxplot(exprs(mdata), outline = FALSE)

fData(mdata) <- MA2$genes

## Fixing sample Info and Gene Info (Annotations ) ------------------

## Phenotype data

exprs <- exprs(mdata)
colnames(exprs) <- strsplit2(colnames(exprs), split = "[_]")[,3]

sampinfo <- read.table("raw_files/sampInfo.txt", header = TRUE)
all(colnames(exprs) %in% sampinfo$Sample.ID)
row.names(sampinfo) <- sampinfo$Sample.ID

sampinfo <- sampinfo[colnames(exprs),]
sampinfo$Group <- gsub("[*]", "", sampinfo$Clinical.Form)
sampinfo$Group <- gsub("C", "Healthy", sampinfo$Group)
sampinfo$Group <- factor(sampinfo$Group)
sampinfo$MDT <- ifelse(grepl("[*]", sampinfo$Clinical.Form) == TRUE, 'Yes', 'No')

## Feature data

genes <- fData(mdata)
annotation <- select(HsAgilentDesign026652.db, keys = genes$ProbeName, columns = c("ENTREZID", "SYMBOL"))

annotation <- cbind(genes$ControlType, annotation)

pheno = new("AnnotatedDataFrame", data = sampinfo)
feature = new("AnnotatedDataFrame", data = annotation )
eset <- ExpressionSet(assayData = exprs, 
                      phenoData = pheno,
                      annotation = "HsAgilentDesign026652.db", 
                      featureData = feature)

rm(list = c("annotation", "feature", "genes", "exprs", "sampinfo"))

## Filtering out Control probes and probes not annotated with ENTREZID

dim(eset)

table(fData(eset)$`genes$ControlType`)
eset <- eset[!fData(eset)$`genes$ControlType` != 0, ]
eset <- eset[!is.na(fData(eset)$ENTREZID), ]

dim(eset)

saveRDS(eset, file = "objects/eset.rds", compress = "bzip2")

## Quality Control ----------------------------

plotMDS(eset, labels = pData(eset)$Group, pch = 19)

pca1 <- prcomp(t(exprs(eset)), center = TRUE, scale. = TRUE)

qplot(
        x = pca1$x[, 1],
        pca1$x[, 2],
        geom = "point",
        col = eset@phenoData@data$Group
) + theme_bw() + geom_point(size = 3) 

## Differential Expression Analysis with Limma -----------------

design <- model.matrix(~ 0 + Group + Gender + Ethnicity, data = pData(eset))
design
colnames(design)[1:8] <- levels(pData(eset)$Group)

contrasts <-
        makeContrasts(
                DiseaseVsCon = (TT + BT + BB + BL + LL + R1 + R2) / 7 - Healthy,
                AllFormsvsCon = (TT + BT + BB + BL + LL) / 5 - Healthy,
                IntermedvsCon = (BT + BB + BL) / 3 - Healthy,
                PoleTvsCon = (TT + BT) / 2 - Healthy,
                PoleLvsCon = (LL + BL) / 2 - Healthy,
                R1vsCon = R1 - Healthy,
                R2vsCon = R2 - Healthy,
                TTvsHealthy = TT - Healthy,
                BTvsHealthy = BT - Healthy,
                BLvsHealthy = BL - Healthy,
                LLvsHealthy = LL - Healthy,
                BBvsHealthy = BB - Healthy, 
                LLvsBT = LL - BT,
                LLvsR1 = LL - R1,
                LLvsR2 = LL - R2, 
                PoleLvsPoleT = (LL + BL) / 2 - (TT + BT) / 2,
                NonTubvsTub = (LL + BL + BB) / 3 - (TT + BT) / 2,
                PolevsPole = LL - TT,
                ClinvsRea1 = (TT + BT + BB + BL) / 4 - R1,
                ClinvsRea2 = (BL + LL) / 2 - R2,
                Reactions = R1-R2,
                TT_BT_BBvsR1 = (TT + BT + BB )/3 - R1,
                
                levels = design)

fit <- lmFit(eset, design)

### Removing duplicate probes --------------

o <- order(fit$Amean, decreasing = TRUE)
fit <- fit[o, ]
d <- duplicated(fit$genes$ENTREZID)

fit <- fit[!d, ]

fit2 <- contrasts.fit(fit, contrasts)
fit3 <- eBayes(fit2)

## Results --------------

tt = topTable(fit3, coef = "PoleLvsPoleT", number = Inf, adjust.method = "fdr", sort.by = "logFC", confint = TRUE)

cont <- colnames(contrasts)
results <- list()

for(i in cont) {
        results[[i]] <-topTable(fit3,
                        coef = i,
                        number = Inf,
                        genelist = fit$genes,
                        adjust.method = "fdr", confint = TRUE)
}
rm(i)

## PCA with only DEG according to omnibus FDR p-value selected from all contrasts (F, p-value, corrected by BH) #####

ind <- row.names(fit3)[p.adjust(fit3$F.p.value, method = "bonferroni") <= 1e-3]
        
pca2 <- prcomp(t(exprs(eset)[ind, ]), center = TRUE, scale. = TRUE)

qplot(
        x = pca2$x[, 1],
        pca2$x[, 2],
        geom = "point",
        xlab = "PC1",
        ylab = 'PC2',
        col = eset@phenoData@data$Group
) + geom_point(size = 3) + labs(colour = "Grupos") + scale_color_brewer(type = "div", palette = "Set1") + theme_bw() 
rm(pca2)

pca2 <- FactoMineR::PCA(X = t(exprs(eset)[ind,]), scale.unit = TRUE)

factoextra::fviz_pca_ind(
        pca2,
        axes = c(1, 2),
        geom = "point",
        pointshape = 19,
        pointsize = 3,
        repel = TRUE,
        palette = "Set1",
        habillage = eset@phenoData@data$Group,
        mean.point = FALSE,
        ggtheme = theme_bw(),
        title = ""
)
ggsave(
        filename = "figs/PCA_topVar.pdf",
        device = "pdf",
        units = 'cm',
        width = 16,
        height = 16
)

row.names(pca2$rotation) <- fit3$genes[ind,]$SYMBOL

fviz_contrib(pca2, choice = "var", axes = 1, fill = "lightgrey", color = "black", top = 20)

ggsave(
        filename = "figs/PCA_ExpVarAx1Genes.pdf",
        device = "pdf",
        units = 'cm',
        width = 15,
        height = 10
)

fviz_contrib(pca2, choice = "var", axes = 2, fill = "lightgrey", color = "black", top = 20)

ggsave(
        filename = "figs/PCA_ExpVarAx2Genes.pdf",
        device = "pdf",
        units = 'cm',
        width = 15,
        height = 10
)


## Exporing DEG tables  --------------------------------

for (i in names(results)) {
        WriteXLS(
                x = results[i],
                ExcelFileName = paste("results/", "results.", i, ".xls", sep = "")
        )
}

rm(i)

## Exporting Expression Matrix ------------------

exp_matrix <- exprs(eset)

exp_matrix <- cbind(exp_matrix, eset@featureData@data)

o <- order(rowMeans(exp_matrix[, 1:76]), decreasing = TRUE)
exp_matrix <- exp_matrix[o, ]
exp_matrix <- exp_matrix[!duplicated(exp_matrix$ENTREZID), ]
exp_matrix <- exp_matrix[!is.na(exp_matrix$ENTREZID), ]

rownames(exp_matrix) <- as.character(exp_matrix$ENTREZID)
exp_matrix <- exp_matrix[, -c(77, 80)]
exp_matrix <- as.data.frame(exp_matrix)
colnames(exp_matrix)[1:76] <- make.names(pData(eset)$Group, unique = TRUE)
write.table(exp_matrix, file = "results/exp_74481.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

### Volcano plot =================

library(RColorBrewer)
library(gplots)
library(ggplot2)
library(limma)

## Volcano Plots with loop

pdf(
        file = "figs_pdf/VolcanoPlots.pdf",
        onefile = TRUE,
        paper = "a4r",
        pagecentre = TRUE
)

for(i in names(results)) {
        df1 <- as.data.frame(results[[i]])
        sig <- with(df1, P.Value < 0.05 & abs(logFC) > 1)
                plot.title <-
                paste(
                        "DEG from GSE74481",
                        i,
                        "P.Val < 0.05 & LogFC > 1",
                        sep = " "
                )
        plot <-     ggplot(data = df1,
                           aes(
                                   x = logFC,
                                   y = -log10(P.Value),
                                   colour = sig
                           )) +
                geom_point(alpha = 1, size = .8) +
                xlab("log2 Fold Change") + ylab("-log10 P-value") + theme_bw() +
                theme(legend.position = "none") + scale_color_manual(values =
                                                                             c("lightgrey", "black")) + ggtitle(plot.title) + geom_hline(yintercept = -log10(0.05), colour = "dodgerblue", linetype = 3) + geom_vline(xintercept = c(-1,1), colour = c("green", "red"), linetype = 3)
        
        print(plot)
      
}
  dev.off()
  
rm(df1, sig, plot.title, plot)

###  Heatmap ----------------------------

# I'm gonna plot just the DEG with Adj.P.val < 0.01 and logFC > 2

library(pheatmap)
library(RColorBrewer)

plot = tt[with(tt, adj.P.Val < 0.001 & logFC > 1),]$ENTREZID

label = data.frame(Group = eset@phenoData@data$Group, MDT = eset@phenoData@data$MDT)
row.names(label) <- make.names(eset@phenoData@data$Group, unique = T)
ann_colors = list(Group = c(BB = "#E41A1C",
                            BL = "#377EB8",
                            BT = "#F781BF",
                            LL = "#984EA3",
                            R1 = "#FF7F00",
                            R2 = "#FFFF33",
                            TT = "#A65628",
                            Healthy = "#4DAF4A"), 
                  MDT = c(No = "lightgrey", Yes = "darkgrey"))

pdf("figs/heatmap_deg.pdf", paper = "a4", pagecentre = TRUE, onefile = FALSE)

pheatmap(exp_matrix[plot, 1:76],
         color = brewer.pal(11, "RdYlBu"),
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         legend = TRUE,
         show_rownames = FALSE,
         main = paste("Heatmap of GSE74481 DEGs (", length(unique(plot)), ") with Adj.P-value < 0.001 & logFC > 1", sep = ""),
         border_color = NA, 
         annotation_col = label,
         annotation_legend = TRUE,
         annotation_colors = ann_colors,
         show_colnames = FALSE)
dev.off()
# 
# 
# ###  Gene Ontology Enrichment Analysis  -----------------------
# 
# # Since there's a lot of different contrasts, I'm gonna perform the enrichment analysis only on some of them, which I consider most informative...
# 
# library(clusterProfiler)
# library(org.Hs.eg.db)
# library(pathview)
# 
# ##
# 
# ## All disease forms against healthy  -------
# 
# universe <- fit3$genes$ENTREZID
# 
# deg1 <-
#         topTable(
#                 fit3,
#                 coef = 'DiseaseVsCon',
#                 adjust.method = "fdr",
#                 p.value = 0.05,
#                 number = Inf
#         )
# 
# genelist1 <- deg1$logFC
# names(genelist1) <- deg1$ENTREZID
# genelist1 <- genelist1[!duplicated(names(genelist1))]
# genelist1 <- sort(genelist1, decreasing = TRUE)
# 
# ## Cellular Component All (-reactions) vs Healthy
# 
# gsea1.cc <-  enrichGO(
#         gene = names(genelist1),
#         universe      = universe,
#         OrgDb         = org.Hs.eg.db,
#         ont           = "CC",
#         pAdjustMethod = "fdr",
#         pvalueCutoff  = 0.05,
#         qvalueCutoff  = 0.05,
#         readable      = TRUE
# )
# head(gsea1.cc)
# write.csv(gsea1.cc@result[, 1:7], file = "results/GO_cc_ALLvsHealthy.csv", row.names = FALSE)
# 
# gsea1.cc <-
#         simplify(gsea1.cc,
#                  cutoff = 0.7,
#                  by = "p.adjust",
#                  select_fun = min)
# 
# pdf(
#         "figs_pdf/CellComponent_ALLvsHealthy.pdf",
#         width = 10,
#         pagecentre = TRUE,
#         onefile = TRUE
# )
# 
# barplot(
#         gsea1.cc,
#         colorBy = "p.adjust",
#         title = "Ontologias gênicas (GO) super-representadas para Componentes Celulares (Doença vs Controles",
#         showCategory = length(gsea1.cc@result$ID),
#         font.size = 12
# )
# 
# dev.off()
# cnetplot(
#         gsea1.cc,
#         showCategory = 3,
#         foldChange = genelist1,
#         fixed = FALSE
# )
# 
# ## Biological Process All vs Healthy
# 
# gsea1.BP <-  enrichGO(
#         gene = names(genelist1),
#         universe      = universe,
#         OrgDb         = org.Hs.eg.db,
#         ont           = "BP",
#         pAdjustMethod = "fdr",
#         pvalueCutoff  = 0.05,
#         qvalueCutoff  = 0.05,
#         readable      = TRUE
# )
# head(gsea1.BP)
# 
# gsea1.BP <-
#         simplify(gsea1.BP,
#                  cutoff = 0.5,
#                  by = "p.adjust",
#                  select_fun = min)
# 
# pdf(
#         "figs_pdf/BioProcess_ALLvsHealthy.pdf",
#         width = 10,
#         pagecentre = TRUE,
#         onefile = TRUE
# )
# 
# barplot(
#         gsea1.BP,
#         colorBy = "p.adjust",
#         title = "Ontologias gênicas (GO) super-representadas para Processos Biológicos (Doença vs Controles",
#         showCategory = length(gsea1.BP@result$ID),
#         font.size = 12
# )
# 
# dev.off()
# cnetplot(
#         gsea1.BP,
#         showCategory = 6,
#         foldChange = genelist1,
#         fixed = FALSE
# )
# 
# write.csv(gsea1.BP@result[, 1:7], file = "results/GO_BP_ALLvsHealthy.csv", row.names = FALSE)
# 
# ## Molecular function All (-reactions) vs Healthy
# 
# gsea1.MF <-  enrichGO(
#         gene = names(genelist1),
#         universe      = universe,
#         OrgDb         = org.Hs.eg.db,
#         ont           = "MF",
#         pAdjustMethod = "fdr",
#         pvalueCutoff  = 0.05,
#         qvalueCutoff  = 0.05,
#         readable      = TRUE
# )
# head(gsea1.MF)
# 
# gsea1.MF <-
#         simplify(gsea1.MF,
#                  cutoff = 0.6,
#                  by = "p.adjust",
#                  select_fun = min)
# 
# pdf(
#         "figs_pdf/MolFunction_ALLvsHealthy.pdf",
#         width = 10,
#         pagecentre = TRUE,
#         onefile = TRUE
# )
# 
# barplot(
#         gsea1.MF,
#         colorBy = "p.adjust",
#         title = "Ontologias gênicas (GO) super-representadas para Funções Moleculares (Doença vs Controles",
#         showCategory = length(gsea1.MF@result$ID),
#         font.size = 12
# )
# 
# dev.off()
# cnetplot(
#         gsea1.MF,
#         showCategory = 6,
#         foldChange = genelist1,
#         fixed = FALSE
# )
# 
# write.csv(gsea1.MF@result[, 1:7], file = "results/GO_MF_ALLvsHealthy.csv", row.names = FALSE)
# 
# ## TT vs Healthy (except reactional states) -------
# 
# deg2 <-
#         topTable(
#                 fit3,
#                 coef = 'TTvsHealthy',
#                 adjust.method = "fdr",
#                 p.value = 0.05,
#                 number = Inf
#         )
# 
# genelist2 <- deg2$logFC
# names(genelist2) <- deg2$ENTREZID
# genelist2 <- genelist2[!duplicated(names(genelist2))]
# genelist2 <- sort(genelist2, decreasing = TRUE)
# 
# ## Cellular Component TT  vs Healthy
# 
# gsea2.cc <-  enrichGO(
#         gene = names(genelist2),
#         universe      = universe,
#         OrgDb         = org.Hs.eg.db,
#         ont           = "CC",
#         pAdjustMethod = "fdr",
#         pvalueCutoff  = 0.05,
#         qvalueCutoff  = 0.05,
#         readable      = TRUE
# )
# head(gsea2.cc)
# 
# gsea2.cc <-
#         simplify(gsea2.cc,
#                  cutoff = 0.7,
#                  by = "p.adjust",
#                  select_fun = min)
# 
# 
# pdf(
#         "figs_pdf/CellComponent_TTvsHealthy.pdf",
#         width = 10,
#         pagecentre = TRUE,
#         onefile = TRUE
# )
# 
# barplot(
#         gsea2.cc,
#         colorBy = "p.adjust",
#         title = "Gene Ontology over-representation test for Cellular Component (TTvsHealthy)",
#         showCategory = length(gsea2.cc@result$ID),
#         font.size = 12
# )
# 
# dev.off()
# cnetplot(
#         gsea2.cc,
#         showCategory = 3,
#         foldChange = genelist2,
#         fixed = FALSE
# )
# 
# write.csv(gsea2.cc@result[, 1:7], file = "results/GO_cc_TTvsHealthy.csv", row.names = FALSE)
# 
# ## Biological Process TT  vs Healthy
# 
# gsea2.BP <-  enrichGO(
#         gene = names(genelist2),
#         universe      = universe,
#         OrgDb         = org.Hs.eg.db,
#         ont           = "BP",
#         pAdjustMethod = "fdr",
#         pvalueCutoff  = 0.05,
#         qvalueCutoff  = 0.05,
#         readable      = TRUE
# )
# head(gsea2.BP)
# 
# gsea2.BP <-
#         simplify(gsea2.BP,
#                  cutoff = 0.5,
#                  by = "p.adjust",
#                  select_fun = min)
# 
# pdf(
#         "figs_pdf/BioProcess_TTvsHealthy.pdf",
#         width = 10,
#         pagecentre = TRUE,
#         onefile = TRUE
# )
# 
# barplot(
#         gsea2.BP,
#         colorBy = "p.adjust",
#         title = "Gene Ontology over-representation test for Biological Process (TTvsHealthy)",
#         showCategory = length(gsea2.BP@result$ID),
#         font.size = 12
# )
# 
# dev.off()
# cnetplot(
#         gsea2.BP,
#         showCategory = 6,
#         foldChange = genelist2,
#         fixed = FALSE
# )
# 
# write.csv(gsea2.BP@result[, 1:7], file = "results/GO_BP_TTvsHealthy.csv", row.names = FALSE)
# 
# ## Molecular function TT vs Healthy
# 
# gsea2.MF <-  enrichGO(
#         gene = names(genelist2),
#         universe      = universe,
#         OrgDb         = org.Hs.eg.db,
#         ont           = "MF",
#         pAdjustMethod = "fdr",
#         pvalueCutoff  = 0.05,
#         qvalueCutoff  = 0.05,
#         readable      = TRUE
# )
# head(gsea2.MF)
# 
# write.csv(gsea2.MF@result[, 1:7], file = "results/GO_MF_ALL-rvsHealthy.csv", row.names = FALSE)
# 
# gsea2.MF <-
#         simplify(gsea2.MF,
#                  cutoff = 0.7,
#                  by = "p.adjust",
#                  select_fun = min)
# 
# pdf(
#         "figs_pdf/MolFunction_TTvsHealthy.pdf",
#         width = 10,
#         pagecentre = TRUE,
#         onefile = TRUE
# )
# 
# barplot(
#         gsea2.MF,
#         colorBy = "p.adjust",
#         title = "Gene Ontology over-representation test for Molecular function (TTvsHealthy)",
#         showCategory = length(gsea2.MF@result$ID),
#         font.size = 12
# )
# 
# dev.off()
# cnetplot(
#         gsea2.MF,
#         showCategory = 6,
#         foldChange = genelist2,
#         fixed = FALSE
# )
# 
# ## Borderline vs Healthy (except reactional states) -------
# 
# deg3 <-
#         topTable(
#                 fit3,
#                 coef = 'IntermedvsCon',
#                 adjust.method = "fdr",
#                 p.value = 0.05,
#                 number = Inf
#         )
# 
# genelist3 <- deg3$logFC
# names(genelist3) <- as.character(deg3$ENTREZID)
# genelist3 <- genelist3[!duplicated(names(genelist3))]
# genelist3 <- sort(genelist3, decreasing = TRUE)
# 
# ## Cellular Component Borderline vs Healthy
# 
# gsea3.cc <-  enrichGO(
#         gene = names(genelist3),
#         universe      = universe,
#         OrgDb         = org.Hs.eg.db,
#         ont           = "CC",
#         pAdjustMethod = "fdr",
#         pvalueCutoff  = 0.05,
#         qvalueCutoff  = 0.05,
#         readable      = TRUE
# )
# head(gsea3.cc)
# 
# gsea3.cc <-
#         simplify(gsea3.cc,
#                  cutoff = 0.8,
#                  by = "p.adjust",
#                  select_fun = min)
# 
# 
# pdf(
#         "figs_pdf/CellComponent_BorderlinevsHealthy.pdf",
#         width = 10,
#         pagecentre = TRUE,
#         onefile = TRUE
# )
# 
# barplot(
#         gsea3.cc,
#         colorBy = "p.adjust",
#         title = "Gene Ontology over-representation test for Cellular Component (BorderlinevsHealthy)",
#         showCategory = length(gsea3.cc@result$ID),
#         font.size = 12
# )
# 
# dev.off()
# cnetplot(
#         gsea3.cc,
#         showCategory = 8,
#         foldChange = genelist3,
#         fixed = FALSE
# )
# 
# write.csv(gsea3.cc@result[, 1:7], file = "results/GO_cc_BorderlinevsHealthy.csv", row.names = FALSE)
# 
# ## Biological Process Borderline vs Healthy
# 
# gsea3.BP <-  enrichGO(
#         gene = names(genelist3),
#         universe      = universe,
#         OrgDb         = org.Hs.eg.db,
#         ont           = "BP",
#         pAdjustMethod = "fdr",
#         pvalueCutoff  = 0.05,
#         qvalueCutoff  = 0.05,
#         readable      = TRUE
# )
# head(gsea3.BP)
# 
# gsea3.BP <-
#         simplify(gsea3.BP,
#                  cutoff = 0.4,
#                  by = "p.adjust",
#                  select_fun = min)
# 
# pdf(
#         "figs_pdf/BioProcess_BorderlinevsHealthy.pdf",
#         width = 10,
#         pagecentre = TRUE,
#         onefile = TRUE
# )
# 
# barplot(
#         gsea3.BP,
#         colorBy = "p.adjust",
#         title = "Gene Ontology over-representation test for Biological Process (BorderlinevsHealthy)",
#         showCategory = length(gsea3.BP@result$ID),
#         font.size = 12
# )
# 
# dev.off()
# cnetplot(
#         gsea3.BP,
#         showCategory = 6,
#         foldChange = genelist3,
#         fixed = FALSE
# )
# 
# write.csv(gsea3.BP@result[, 1:7], file = "results/GO_BP_BorderlinevsHealthy.csv", row.names = FALSE)
# 
# ## Molecular function Borderline vs Healthy
# 
# gsea3.MF <-  enrichGO(
#         gene = names(genelist3),
#         universe      = universe,
#         OrgDb         = org.Hs.eg.db,
#         ont           = "MF",
#         pAdjustMethod = "fdr",
#         pvalueCutoff  = 0.05,
#         qvalueCutoff  = 0.05,
#         readable      = TRUE
# )
# head(gsea3.MF)
# 
# pdf(
#         "figs_pdf/MolFunction_BorderlinevsHealthy.pdf",
#         width = 10,
#         pagecentre = TRUE,
#         onefile = TRUE
# )
# 
# barplot(
#         gsea3.MF,
#         colorBy = "p.adjust",
#         title = "Gene Ontology over-representation test for Molecular function (BorderlinevsHealthy)",
#         showCategory = length(gsea3.MF@result$ID),
#         font.size = 12
# )
# 
# dev.off()
# cnetplot(
#         gsea3.MF,
#         showCategory = 15,
#         foldChange = genelist3,
#         fixed = FALSE
# )
# 
# write.csv(gsea3.MF@result[, 1:7], file = "results/GO_MF_ALL-rvsHealthy.csv", row.names = FALSE)

#### END ###
rm(list=ls())
q(save = "no")
