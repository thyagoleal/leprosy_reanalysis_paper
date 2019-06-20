## Analysis Microarray GSE24280 - July, 18th 2017.
## Last update: June 1st 2018.
## Thyago Leal Calvo - thyagoleal@yahoo.com
############################################################################

## Settings -----------------------------

root <- getwd()
set.seed(1290)
options(digits = 3)

### Loading Packages ----------

require(oligo)
require(limma)
require(Biobase)
require(geneplotter)
library(WriteXLS)
library(GEOquery)

### Converting pairs to xys ---------

pairfiles = list.files("raw_files/pairs/", pattern = ".pair", full.names = TRUE)

## CONVERSION TOOL - Benilton Carvalho - June/2013

pair2xys <- function(pairFiles, outdir=getwd(), verbose=TRUE){
        if (verbose) message('Output directory: ', outdir)
        for (pairFile in pairFiles){
                if (verbose) message('Processing ', basename(pairFile))
                header <- readLines(pairFile, n=1)
                pair <- read.delim(pairFile, header=TRUE, sep='\t',
                                   stringsAsFactors=FALSE, comment.char='#')
                maxX <- max(pair$X)
                maxY <- max(pair$Y)
                xys <- expand.grid(X=1:maxX, Y=1:maxY)
                xys <- merge(xys, pair[, c('X', 'Y', 'PM')], all.x=TRUE)
                xys <- pair[, c('X', 'Y', 'PM')]
                names(xys) <- c('X', 'Y', 'SIGNAL')
                xys$COUNT <- ifelse(is.na(xys$SIGNAL), NA_integer_, 1L)
                xys <- xys[with(xys, order(Y, X)),]
                rownames(xys) <- NULL
                xysFile <- file.path(outdir, gsub('\\.pair$', '\\.xys', basename(pairFile)))
                if (verbose) message('Writing ', basename(xysFile))
                writeLines(header, con=xysFile)
                suppressWarnings(write.table(xys, file=xysFile, sep='\t',
                                             row.names=FALSE, quote=FALSE, append=TRUE))
        }
}
### END CONVERSION TOOL

## Converting  --------------------------------------

pair2xys(pairFiles = pairfiles, outdir = "raw_files/", verbose = TRUE)
rm(pairfiles)

## Reading files with oligo --------------------------------------

xysfiles <- list.files(path = "raw_files/", pattern = ".xys",
                       full.names = TRUE)

### Reading XysFiles  ------------------------------------------

install.packages("pd.gpl10191.090828.hg18.opt.exp/", repos = NULL, type = "source")

expressionset <- read.xysfiles(xysfiles, pkgname = "pd.gpl10191.090828.hg18.opt.expr")
boxplot(exprs(expressionset))

## RMA BG correction and normalization =========

eset <- rma(expressionset)
boxplot(eset, outline=FALSE)

## Getting GEO for annotation -----------------------

gse <- getGEO("GSE24280", AnnotGPL = TRUE, destdir = "raw_files/")
gse <- gse$GSE24280_series_matrix.txt.gz

pdata <- gse@phenoData@data[,c("source_name_ch1", "characteristics_ch1")]

head(pdata)
condition = base::paste(pdata$source_name_ch1, pdata$characteristics_ch1, sep = ".")
condition <- gsub(" ", "", condition)
condition <- gsub("*patient[.]tissue*", "", condition)
condition <- gsub("[(:)]", ".", condition)
condition <- gsub("Skinbiopsy", "t", condition)
condition <- gsub("*Blood", "b", condition)
condition[25] <- "healthy.t"
pdata = data.frame(samples = as.factor(rownames(pdata)), Condition = condition)
row.names(pdata) <- pdata$samples
pdata

pdata$samples == colnames(exprs(eset))

dimnames(eset)[[2]] <- strsplit2(colnames(exprs(eset)), split = "_")[,1]

dim(eset)
rm(gse)

## Creating Annotation Package for this chip --------------------------
## Run only if the package is not already installed 

if(!"hg18.opt.db" %in% row.names(installed.packages())) {
        require(AnnotationForge)
        require(AnnotationDbi)
        
        makeDBPackage(
                "HUMANCHIP_DB",
                affy = FALSE,
                prefix = "hg18.opt",
                fileName = "raw_files/annotation.txt",
                baseMapType = "eg",
                outputDir = getwd(),
                version = "1.0.0",
                manufacturer = "Nimblegene"
        )
        
        install.packages("hg18.opt.db/", repos = NULL, type = "source")
} else{
        cat("Package is already installed, loading it")
        require(hg18.opt.db)
}

annotation(eset) <- "hg18.opt.db"

## Creating ExpressionSet with Annotated ---------------

require(hg18.opt.db)

exprs <- exprs(eset)

featdata <-
        select(
                hg18.opt.db,
                keys = as.character(row.names(exprs)),
                columns = c("ENTREZID", "SYMBOL", "GENENAME")
        )


colnames(exprs) <- pdata$samples
pdata <- new("AnnotatedDataFrame", data = pdata)
featuredata <- new("AnnotatedDataFrame", data = featdata)
row.names(featuredata) <- featuredata@data$PROBEID
eset <- ExpressionSet(
        assayData = exprs,
        phenoData = pdata,
        featureData = featuredata,
        annotation = "hg18.opt.db"
)

dim(eset)

saveRDS(eset, file = "objs/eSet.Rds", compress = "bzip2")

# Analysis --------------------------

design <- model.matrix( ~ 0 + pdata$Condition)
colnames(design) <- levels(pdata$Condition)

fit <- lmFit(eset, design)
contrast.matrix <- makeContrasts(MB.b - BT.b,
                                 MB.t - healthy.t,
                                 MB.t - BT.t,
                                 BT.t - healthy.t,
                                 levels = design)

# Removing Duplicate probes --------

o <- order(fit$Amean, decreasing = TRUE)
fit <- fit[o, ]
d <- duplicated(fit$genes$ENTREZID)
fit <- fit[!d,]

fit <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit)

tt <- topTable(fit2,
               adjust.method = "fdr",
               number = Inf, confint = TRUE)

BloodMBvsBloodBT <- topTable(fit2,
                             coef = 1,
                             adjust.method = "fdr",
                             number = Inf, 
                             confint = TRUE)

TissueMBvsHealthy <- topTable(fit2,
                              coef = 2,
                              adjust.method = "fdr",
                              number = Inf, 
                              confint = TRUE)


TissueMBvsTissueBT <- topTable(fit2,
                               coef = 3,
                               adjust.method = "fdr",
                               number = Inf, 
                               confint = TRUE)

TissueBTvsHealthy <- topTable(fit2,
                              coef = 4,
                              adjust.method = "fdr",
                              number = Inf,
                              confint = TRUE)

## Exporing DEG table  --------------------------------

list <- list(BloodMBvsBloodBT, TissueMBvsHealthy, TissueMBvsTissueBT, TissueBTvsHealthy)
names(list) <- c("BloodMBvsBloodBT", "TissueMBvsHealthy", "TissueMBvsTissueBT", "TissueBTvsHealthy")

## Saving all the tables for each contrasts
for(i in names(list)){
        filename <- paste("results/",i,".xls", sep = "")
        WriteXLS(x = i, ExcelFileName = filename)
}
rm(i, filename)

## Exporting Expression Matrix ------------------

expr_matrix <- exprs(eset)
colnames(expr_matrix) <- make.names(pdata$Condition, unique = TRUE)
expr_matrix <- merge.data.frame(expr_matrix, featdata, by.x = "row.names", by.y = "PROBEID")
row.names(expr_matrix) <- expr_matrix$Row.names
expr_matrix <- expr_matrix[row.names(fit2),]

expr_matrix <- expr_matrix[!is.na(expr_matrix$ENTREZID),]
expr_matrix <- expr_matrix[!duplicated(expr_matrix$ENTREZ),-1]

write.table(expr_matrix, file = "results/exp_24280.txt", sep = "\t")

### Volcano plot =================

library(RColorBrewer)
library(gplots)
library(ggplot2)
library(limma)

## Blood MB vs BT -------------------------------

pdf("figs/VolcanoPlots.pdf", onefile= TRUE)

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
        geom_point(alpha = 1, size = .8) +
        xlab("log2 Fold Change") + ylab("-log10 P-value") + theme_bw() +
        theme(legend.position = "none") + scale_color_manual(values=c("lightgrey", "black"))

plot_coef1 <- plot_coef1 + ggtitle("Differentially Expressed Genes from GSE24280 - Blood MB vs BT") + geom_hline(yintercept = -log10(0.05), colour = "dodgerblue", linetype = 3) + geom_vline(xintercept = c(-1,1), colour = c("green", "red"), linetype = 3)


plot_coef1

## Tissue Mb vs Control  -------------------------------

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
        geom_point(alpha = 1, size = .8) +
        xlab("log2 Fold Change") + ylab("-log10 P-value") + theme_bw() +
        theme(legend.position = "none") + scale_color_manual(values=c("lightgrey", "black"))

plot_coef2 <- plot_coef2 + ggtitle("Differentially Expressed Genes from GSE24280 - Tissue MB vs Control") + geom_hline(yintercept = -log10(0.05), colour = "dodgerblue", linetype = 3) + geom_vline(xintercept = c(-1,1), colour = c("green", "red"), linetype = 3)

plot_coef2

## Tissue Mb vs BT  -------------------------------

coef3 <- topTable(fit2,
                  coef = 3,
                  adjust.method = "fdr",
                  number = Inf)

Significant = coef3$P.Value < 0.05 & abs(coef3$logFC) > 1

names(coef3)

plot_coef3 <- ggplot(data = coef3, 
                     aes(x = logFC,
                         y = -log10(P.Value), 
                         colour = Significant)) +
        geom_point(alpha = 1, size = .8) +
        xlab("log2 Fold Change") +
        ylab("-log10 P-value") +
        theme_bw() + 
        theme(legend.position = "none") + scale_color_manual(values=c("lightgrey", "black"))

plot_coef3 <- plot_coef3 + ggtitle("Differentially Expressed Genes from GSE24280 - Tissue MB vs BT") + geom_hline(yintercept = -log10(0.05), colour = "dodgerblue", linetype = 3) + geom_vline(xintercept = c(-1,1), colour = c("green", "red"), linetype = 3)

plot_coef3

## Tissue BT vs Control  -------------------------------

coef4 <- topTable(fit2,
                  coef = 4,
                  adjust.method = "fdr",
                  number = Inf)

Significant = coef4$P.Value < 0.05 & abs(coef4$logFC) > 1

names(coef4)

plot_coef4 <- ggplot(data = coef4, 
                     aes(x = logFC,
                         y = -log10(P.Value), 
                         colour = Significant)) +
        geom_point(alpha = 1, size = .8) +
        xlab("log2 Fold Change") +
        ylab("-log10 P-value") +
        theme_bw() + 
        theme(legend.position = "none") + scale_color_manual(values=c("lightgrey", "black"))

plot_coef4 <- plot_coef4 + ggtitle("Differentially Expressed Genes from GSE24280 - Tissue BT vs Control") + geom_hline(yintercept = -log10(0.05), colour = "dodgerblue", linetype = 3) + geom_vline(xintercept = c(-1,1), colour = c("green", "red"), linetype = 3)

plot_coef4
dev.off()

# Heatmap ----------------------------

library(pheatmap)
library(RColorBrewer)

DEG <- topTable(fit2, number = Inf, coef = 3, adjust.method = "BH")

plot = DEG$PROBEID[DEG$adj.P.Val < 0.05]

label = data.frame(Class = pdata$Condition)
row.names(label) <- make.names(pdata$Condition, unique = T)
ann_colors = list(Class = c(BT.b = "darkgrey",
                            MB.b = "firebrick1",
                            BT.t = "yellow",
                            MB.t = "purple",
                            healthy.t = "limegreen"))

pdf("figs/heatmap_deg.pdf", onefile = FALSE)
pheatmap(na.omit(expr_matrix[plot, 1:25]),
         color = brewer.pal(10, "RdYlBu"),
         scale = "row",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "euclidean",
         clustering_method = "average",
         legend = TRUE,
         show_rownames = FALSE,
         main = paste("Heatmap of GSE24280 DEGs (", length(plot), ") with Adj.P-value < 0.05", sep = ""),
         border_color = NA, 
         annotation_col = label,
         annotation_legend = TRUE,
         annotation_colors = ann_colors, show_colnames = FALSE)
dev.off()

## Venn Diagram ---------------------

require(limma)

results <- decideTests(fit2, adjust.method = "fdr", p.value = 0.05)


pdf("figs/VennDiagram.pdf", paper = "a4", pagecentre = TRUE)

vennDiagram(results, include = c("up", "down"), 
            counts.col = c("firebrick1", "dodgerblue2"), 
            circle.col = c("limegreen", "yellow", "dodgerblue2", "firebrick1"))

dev.off()

## PCA with all the DEG and all samples -----------------------------------------

p <- DEG$PROBEID[DEG$adj.P.Val < 0.05]

require(FactoMineR)
require(factoextra)

pca.all <- PCA(t(expr_matrix[rownames(expr_matrix) %in% p,1:25]),  scale.unit = TRUE)

pdf("Figs_pdf/PCA_all_DEG_samples.pdf")
fviz_pca_ind(pca.all, geom = "point", 
             habillage = pdata@data$Condition,
             pointsize = 3, title = "PCA from All DEG (187) and all samples") + theme_bw()
dev.off()

## END
rm(list=ls())
q(save = "no")
