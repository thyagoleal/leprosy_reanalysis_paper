### Reanalysis of GSE100853 
## Thyago Leal - thyagoleal@yahoo.com
## September June 1st 2018
########################################################################

## Settings ------------------------

options(digits = 3, scipen = 0, download.file.method = "libcurl")
set.seed(1209)
root <- getwd()

## Libraries ------------------------

require(limma)
require(illuminaHumanv4.db)
require(GEOquery)
require(Biobase)

## I decided to use the normalized data from GEO instead of the "raw", because I can make mistakes in pairing the samples with the headings of the raw data. Although, I'm gonna use the P-values from the RAW data to subset my expression set from GEO.

raw <- read.ilmn("raw/GSE100853_non-normalized.txt",
                 probeid = "ID_REF", expr = "PS", other.columns = "Detection Pval", sep = "\t")

raw.p <- neqc(raw, detection.p = raw$other$`Detection Pval`)

expressed <- row.names(raw.p$other$`Detection Pval`)[rowSums(raw.p$other$`Detection Pval` < 0.05) >= 3]

## Grabbing data from GEO ------------------------------------

## Since the experimental data are only available through GEO, I'll get data from there.

geo.data <- getGEO("GSE100853", destdir = "raw/", parseCharacteristics = FALSE)
geo.data <- geo.data$GSE100853_series_matrix.txt.gz

pData(geo.data)

## creating phenotype data -------------------------

targets <-
        data.frame(
                samples = colnames(geo.data@assayData$exprs),
                treatment = ifelse(
                        grepl("non-stimulated", geo.data@phenoData@data$source_name_ch1),
                        yes = "Control",
                        no = "Treated"),
                Sex = as.factor(strsplit2(geo.data@phenoData@data$characteristics_ch1.1, " ")[,2]),
                Age = as.integer(strsplit2(geo.data@phenoData@data$characteristics_ch1, " ")[,2]),
                Stimulation.Time = strsplit2(geo.data@phenoData@data$characteristics_ch1.2, " ")[,4],
                Dose = strsplit2(geo.data@phenoData@data$characteristics_ch1.3, " ")[,2])

patients = unname(geo.data@phenoData@data$title)
patients = strsplit2(patients, split = " ")[,1]
patients = strsplit2(patients, split = "_")[,1]
patients = strsplit2(patients, split = "PS")[,2]

targets <- cbind(patients, targets)
targets

row.names(targets) <- targets$samples

targets$Stimulation.Time <- as.integer(with(targets, gsub(pattern = "h", x = Stimulation.Time, replacement = "")))

targets$Dose <- factor(with(targets, gsub(pattern = "ug/mL", x = Dose, replacement = "")))

str(targets)

## Removing probes and log2 ----------------------------------------------

# The detection values contain P values for testing whether each probe is more intense than the negative control probes. Small values are evidence that the probe corresponds to a truly expressed gene. So, I'm gonna use this to filter out probes that are not expressed in at least 10 arrays, according to a detection p-values of 5%.

eSet <- geo.data[expressed,]
plotMDS(eSet, labels = targets$treatment)
pData(eSet) <- targets

exprs(eSet) <- log2(exprs(eSet))

## Fixing annotation - Getting ENTREZID --------------

entrez <- select(illuminaHumanv4.db, keys = fData(eSet)$GB_ACC, columns = c("ENTREZID", "SYMBOL", "GENENAME"), keytype =  "REFSEQ")

entrez <- entrez[!is.na(entrez$ENTREZID), ]

eSet <- eSet[fData(eSet)$GB_ACC %in% entrez$REFSEQ,]
all(fData(eSet)$GB_ACC == entrez$REFSEQ)
eSet@featureData <- new("AnnotatedDataFrame", data = entrez)

saveRDS(eSet, file = "obj/eSet.Rds", compress = "bzip2")

## Differential Expression Analysis ----------------------------------------

patients <- factor(targets$patients)
sex <- factor(targets$Sex)
age <- as.integer(targets$Age)
time <- as.factor(targets$Stimulation.Time)
dose <- factor(targets$Dose)

treatment <- factor(targets$treatment, levels = c("Control", "Treated"))

design <- model.matrix(~treatment + dose + sex+time)

cor.dup <- duplicateCorrelation(eSet, design, block = patients)
cor.dup$consensus

fit <- lmFit(eSet, design, block = patients, correlation = cor.dup$consensus)

## Filtering out duplicated probes ---------------

o <- order(fit$Amean, decreasing = TRUE)
fit <- fit[o,]
d <- duplicated(fit$genes$ENTREZID)
fit <- fit[!d,]

fit2 <- eBayes(fit)

top <- topTable(fit2, coef = "treatmentTreated", adjust = "fdr", genelist = fit$genes, number = Inf, confint = TRUE)

## counting the Adjusted P-values obtained (just for fun) ---------

pvals <- c(0.05, 0.005, 0.0005, 0.0001, 0.00001)

for(p in pvals){
        print(paste(
                "Number of DEG with Adj.P-Val less than or equal to ", p,": ",
                length(which(top$adj.P.Val <= p  & abs(top$logFC)>1)),
                sep = ""))
}

rm(p, pvals)

## Saving results -----------------------------

library(WriteXLS)

WriteXLS(top, ExcelFileName = "results/results_GSE100583.xls", verbose = TRUE)

## Expression matrix 

exp_matrix <- as.data.frame(exprs(eSet))
exp_matrix <- exp_matrix[row.names(fit), ]

exp_matrix <- cbind(fit$genes[,2:3], exp_matrix)

exp_matrix <- exp_matrix[!is.na(exp_matrix$ENTREZID), ]

write.table(exp_matrix, file = "results/exp.matrix.GSE100853.txt", sep = "\t", col.names = TRUE, row.names = FALSE)

### Volcano plot ------------------

library(RColorBrewer)
library(gplots)
library(ggplot2)

significant <-
        top$P.Value < 0.05 & abs(top$logFC) > 1
table(significant)


names(top)

plot1 <- ggplot(data = top,
                aes(
                        x = logFC,
                        y = -log10(P.Value),
                        colour = significant
                )) +
        geom_point(alpha = 1, size = .8) +
        xlab("log2 Fold Change") + ylab("-log10 P-value") + theme_bw() +
        theme(legend.position = "none") + scale_color_manual(values = c("lightgrey", "black"))+ geom_hline(yintercept = -log10(0.05), colour = "dodgerblue", linetype = 3) + geom_vline(xintercept = c(-1,1), colour = c("green", "red"), linetype = 3) 

plot1 <-
        plot1 + ggtitle("Differentially Expressed Genes (196) from GSE100853 - Treated vs Control")

pdf(
        "figs/VolcanoPlot.pdf",
        onefile = TRUE,
        paper = "a4r",
        pagecentre = TRUE
)
plot1
dev.off()

# Heatmap ----------------------------

library(pheatmap)
library(RColorBrewer)

plot <- na.omit(top$ENTREZID[top$adj.P.Val < 0.05 & abs(top$logFC)>1 ])

label = data.frame(Treatment = targets$treatment, Dose = targets$Dose, Stimulus.Time = targets$Stimulation.Time)
row.names(label) <- targets$samples
ann_colors = list(Treatment = c(Control = "dodgerblue",
                            Treated = "firebrick1"),
                  Dose = c("0" = "gray90",
                           "5" = "gray50",
                           "20" = "gray10"),
                  Stimulus.Time = c("26" = "lightyellow",
                                    "32" = "darkorange"))

pdf(
        "figs/heatmap_deg.pdf",
        onefile = FALSE
)
pheatmap(
        exp_matrix[exp_matrix$ENTREZID %in% plot, 3:104],
        color = brewer.pal(11, "RdYlBu"),
        scale = "row",
        clustering_distance_rows = "correlation",
        clustering_distance_cols = "euclidean",
        clustering_method = "average",
        legend = TRUE,
        show_rownames = FALSE,
        main = paste("Heatmap of GSE100853 DEGs (", length(plot), ") with Adj.P-value < 0.05", sep = ""),
        border_color = NA,
        annotation_col = label,
        annotation_legend = TRUE,
        annotation_colors = ann_colors,
        show_colnames = FALSE
)
dev.off()

###  Gene Ontology Enrichment Analysis  -----------------------

require(clusterProfiler)  
library(org.Hs.eg.db)
library(pathview)

universe <- results$ENTREZID

deg <- results[results$adj.P.Val < 0.05, ]
deg <- deg[, c("logFC", "ENTREZID")]
deg1 <- deg$logFC
names(deg1) <- as.character(deg$ENTREZID)
deg1 <- deg1[!duplicated(names(deg1))]
deg1 <- sort(deg1, decreasing = TRUE) ## geneList with Log2FC + Entrez for DEG with P-val < 0.05

### Cell component ------
        
deg1.cc <-  enrichGO(gene = names(deg1),
                             universe      = universe,
                             OrgDb         = org.Hs.eg.db,
                             ont           = "CC",
                             pAdjustMethod = "fdr",
                             pvalueCutoff  = 0.05,
                             qvalueCutoff  = 0.05,
                             readable      = TRUE)
head(deg1.cc)

pdf("figs_pdf/CellComponent.pdf", width = 10, pagecentre = TRUE, onefile = TRUE)

barplot(
        deg1.cc,
        colorBy = "p.adjust",
        title = "Ontologias Gênicas (GO) super-representadas para Componentes Celulares (Tratado vs Controle)",
        showCategory = length(deg1.cc@result$ID),
        font.size = 12
)

dev.off()
cnetplot(deg1.cc, showCategory = 10, foldChange = deg1, fixed = FALSE)

write.csv(deg1.cc@result[, 1:7], file = "results/GO_cc_LLvsENL.csv", row.names = FALSE)


## Biological Process  -------
        
deg1.bp <- enrichGO(gene = names(deg1),
                            universe      = universe,
                            OrgDb         = org.Hs.eg.db,
                            ont           = "BP",
                            pAdjustMethod = "fdr",
                            pvalueCutoff  = 0.05,
                            qvalueCutoff  = 0.05,
                            readable      = TRUE)
head(deg1.bp)

deg1.bp@result$Description[36] <- "adaptive immune response based on somatic recombination" 

pdf("figs_pdf/BioProcess.pdf", width = 10, pagecentre = TRUE, onefile = TRUE)

barplot(deg1.bp, 
        colorBy = "p.adjust",
        title = "Ontologias Gênicas (GO) super-representadas para Processos Biológicos (Tratado vs Controle)",
        showCategory = length(deg1.bp@result$ID),
        font.size = 11)

dev.off()

cnetplot(deg1.bp, showCategory = 4, foldChange = deg1, fixed = FALSE)

write.csv(deg1.bp@result[,1:7], file = "results/GO_BP_LLvsENL.csv", row.names = FALSE)

## Molecular Function  -----------
        
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
        title = "Ontologias Gênicas (GO) super-representadas para Funções Moleculares (Tratado vs Controle)",
        showCategory = length(deg1.mf@result$ID),
        font.size = 12)

dev.off()

write.csv(deg1.mf@result[,1:7], file = "results/GO_MF_LLvsENL.csv", row.names = FALSE)

### END ### 
rm(list=ls())
q(save = "no")



