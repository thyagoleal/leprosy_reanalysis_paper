## sdef - synthesize lists of significant features in related experiments - LLvsControl
## Last update: June 14th, 2019.
## Thyago Leal Calvo, thyagoleal@yahoo.com
## ****************************************

## Settings ------------

options(digits = 3, scipen = 999)
root <- getwd()
dir()
set.seed(1918)
require(sdef)
require(readxl)
require(dplyr)
require(ggplot2)
require(AnnotationDbi)
require(org.Hs.eg.db)
require(WriteXLS)

## LL vs Control ------------
## datasets used: 24280, 74481

GSE24280_LLvscontrol <- read_xls("../GSE24280/results/TissueMBvsHealthy.xls")
GSE24280_LLvscontrol <- as.data.frame(GSE24280_LLvscontrol)

GSE74481_LLvscontrol<- read_xls("../GSE74481/results/results.LLvsHealthy.xls")
GSE74481_LLvscontrol <- as.data.frame(GSE74481_LLvscontrol)

tmp.LLvsCon <- list(GSE24280_LLvscontrol, GSE74481_LLvscontrol)
names(tmp.LLvsCon) <- c("GSE24280_LLvscontrol", "GSE74481_LLvscontrol")

mylist <- list()

for(i in names(tmp.LLvsCon)) {
        df <- tmp.LLvsCon[[i]]
        df <- df[!duplicated(df$ENTREZID), ]
        mylist[[i]] <- df
}

rm(i, df)

merged <-
        Reduce(function(x, y)
                merge.data.frame(x, y, all = FALSE, by = "ENTREZID"),
               mylist)
merged <- merged[, c(1, grep("^P.Value.*", colnames(merged)))]
row.names(merged) <- as.character(merged$ENTREZID)
colnames(merged)[2:3] <- names(mylist)
merged <- merged[, 2:3]

# Remove ribosomal genes 
rib <- read.table("ribosomal_proteins.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
head(rib, 5)
merged <- merged[!row.names(merged) %in% as.character(rib$GeneID), ]

if(!file.exists("LLvsControl/objects/LLvsCon_Pvals.Rds")) {
        saveRDS(merged, file = "LLvsControl/objects/LLvsCon_Pvals.Rds", 
                compress = "bzip2")
        message("Rds object saved!\n")
}

setwd("LLvsControl/")

Th_LLvsCon <- ratio(data = merged, pvalue = TRUE, interval = 0.01)
Rh_LLvsCon <- baymod(iter = 500, output.ratio = Th_LLvsCon)
MC_LLvsCon <- Tmc(iter = 200, output.ratio = Th_LLvsCon)
feat.names <- row.names(merged)
feat.lists.Bayesian <-
        extractFeatures.R(output.ratio = Th_LLvsCon,
                          output.bay = Rh_LLvsCon,
                          feat.names = feat.names,
                          h = Rh_LLvsCon)

createTable(output.ratio = Th_LLvsCon, output.bay = Rh_LLvsCon, h = 0.02)

annt_LLvsCon <-
        select(
                org.Hs.eg.db,
                keys = as.character(feat.lists.Bayesian$max$Names),
                columns = c("SYMBOL", "GENENAME"),
                keytype = "ENTREZID"
        )

## Grabbing Adjusted P-values ----------------------------------

merged <-
        Reduce(function(x, y)
                merge.data.frame(x, y, all = FALSE, by = "ENTREZID"),
               mylist)

merged <- merged[, c(1, grep("^adj.P.Val|^logFC", x = colnames(merged)))]
colnames(merged)[-1] <- paste(rep(c("logFC", "adj.P.Val"), 2), rep(names(mylist), each = 2), sep = "_")
row.names(merged) <- as.character(merged$ENTREZID)
results_LLvsCon_sdef <- merge.data.frame(merged[annt_LLvsCon$ENTREZID, ], annt_LLvsCon, by.x = "row.names", by.y = "ENTREZID") 
names(results_LLvsCon_sdef)

## Summary Statistics across comparisons --------------------------

# Median log2FC
results_LLvsCon_sdef$Median_Log2FC <-
        apply(results_LLvsCon_sdef[, grep("^logFC", colnames(results_LLvsCon_sdef))], MARGIN = 1, FUN = median)

# Absolute Median log2FC
results_LLvsCon_sdef$`|medianLogFC|` <-
        abs(results_LLvsCon_sdef$Median_Log2FC)

# Number of DEG with FDR < 0.1 (10%)
results_LLvsCon_sdef$Significant <-
        apply(results_LLvsCon_sdef[, grep("^adj.P.Val", colnames(results_LLvsCon_sdef))], 1, function(x)
                sum(abs(x < 0.1)))

# Unicode Arrows (symbols) indicating up- or down-regulation
results_LLvsCon_sdef$Direction <-
        as.character(apply(
                apply(
                        results_LLvsCon_sdef[, grep("^logFC", colnames(results_LLvsCon_sdef))],
                        MARGIN = 2,
                        FUN = function(x) {
                                ifelse(x < 0, "\U2193", "\U2191")
                        }
                ),
                1,
                paste,
                collapse = " "
        ))

# Ordering columns by index (I know, I should be using names...)
colnames(results_LLvsCon_sdef)
results_LLvsCon_sdef <- results_LLvsCon_sdef[, c(2, 7, 8, 3:6, 9:12)]

## Exporting *.XLS file --------------------------------

WriteXLS(x = results_LLvsCon_sdef, ExcelFileName = "Results_sdef_LLvsCon.xls")

## Organizing files ---------------------------------

file.rename(list.files(pattern = "*.csv", full.names = TRUE),
            to = gsub("./", "output/", x = list.files(
                    pattern = "*.csv", full.names = TRUE
            )))

file.rename(list.files(pattern = "*.ps", full.names = TRUE),
            to = gsub("./", "figs/", x = list.files(
                    pattern = "*.ps", full.names = TRUE
            )))

file.rename("dataratio .Rdata", "objects/dataratio.Rdata")

if(getwd() != root) setwd(root)

## End program -----------------
rm(list = ls())
q(save = 'no')
