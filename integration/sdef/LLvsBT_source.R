## sdef - synthesize lists of significant features in related experiments - LLvsBT
## Last update: 14th June, 2019.
## Thyago Leal Calvo, thyagoleal@yahoo.com
# ****************************************

## Settings ------------

root <- getwd()
set.seed(1918)
library(sdef)
library(readxl)
library(dplyr)
library(ggplot2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(WriteXLS)

##### 1st Context  - Experiments involving LL vs BT comparisons ############################

# datasets needed: GSE24280, 17763, 74481, 443.

GSE24280_LLvsBT <- as.data.frame(read_xls("../GSE24280/results/TissueMBvsTissueBT.xls"))

GSE17763_LLvsBT <- as.data.frame(read_xls("../GSE17763/results/results_LLvsBT.xls"))

GSE74481_LLvsBT <- as.data.frame(read_xls("../GSE74481/results/results.LLvsBT.xls"))

GSE443_LLvsBT <- as.data.frame(read_xls("../GSE443/results/DEG_full_results.xls"))

tmp.LLvsBT <- list(GSE24280_LLvsBT, GSE17763_LLvsBT, GSE74481_LLvsBT, GSE443_LLvsBT)
names(tmp.LLvsBT) <- c("GSE24280_LLvsBT", "GSE17763_LLvsBT", "GSE74481_LLvsBT", "GSE443_LLvsBT")

mylist <- list()

for(i in names(tmp.LLvsBT)) {
        df <- tmp.LLvsBT[[i]]
        df <- df[!duplicated(df$ENTREZID), ]
        df <- df[!is.na(df$ENTREZID), ]
        mylist[[i]] <- df
}

rm(i, df)

merged <-
        Reduce(function(x, y)
                merge.data.frame(x, y, all = TRUE, by = "ENTREZID"),
               mylist)

merged <- merged[, c(1, grep("P.Value.*", colnames(merged)))]
row.names(merged) <- as.character(merged$ENTREZID)
colnames(merged)[2:5] <- names(mylist)
merged <- merged[, 2:5]

## Keeping only genes present in 3 out of 4 studies. 
keep <- (apply(merged, 1, function(x)sum(abs(is.na(x)))) <= 1)
table(keep)
merged <- merged[keep, ]

# Replacing missing p-values with 0.5
merged[is.na(merged)] <- 0.5

# Remove ribosomal genes 
rib <- read.table("ribosomal_proteins.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
head(rib, 5)
merged <- merged[!row.names(merged) %in% as.character(rib$GeneID), ]

setwd("LLvsBT/")

Th_LLvsBT <- ratio(data = merged,
                   pvalue = TRUE,
                   interval = 0.1)

Rh_LLvsBT <- baymod(iter = 500, output.ratio = Th_LLvsBT)

MC_LLvsBT <- Tmc(iter = 200, output.ratio = Th_LLvsBT)

feat.names <- row.names(merged)

feat.lists.Bayesian <-
        extractFeatures.R(
                output.ratio = Th_LLvsBT,
                output.bay = Rh_LLvsBT,
                feat.names = feat.names,
                h = 0.1
        )

createTable(output.ratio = Th_LLvsBT, output.bay = Rh_LLvsBT)

annt_LLvsBT <-
        select(
                org.Hs.eg.db,
                keys = as.character(feat.lists.Bayesian$User$`h= 0.1`$Names),
                columns = c("SYMBOL", "GENENAME"),
                keytype = "ENTREZID"
        )

results_LLvsBT_sdef <-
        merge.data.frame(feat.lists.Bayesian$User$`h= 0.1`,
                         annt_LLvsBT,
                         by.x = "Names",
                         by.y = "ENTREZID") 

## Grabbing Adjusted P-values ----------------------------------

merged <-
        Reduce(function(x, y)
                merge.data.frame(x, y, all = TRUE, by = "ENTREZID"),
               mylist)

merged <- merged[, c(1, grep("^adj.P.Val.|^logFC.*", colnames(merged)))]
colnames(merged)[-1] <- paste(rep(c("logFC", "adj.P.Val."), times = 4), rep(names(mylist), each = 2), sep = "_")
row.names(merged) <- as.character(merged$ENTREZID)
results_LLvsBT_sdef <- merge.data.frame(merged[annt_LLvsBT$ENTREZID, ], annt_LLvsBT, by.x = "row.names", by.y = "ENTREZID") 
names(results_LLvsBT_sdef)

## Summary Statistics across comparisons --------------------------

# Median log2FC
results_LLvsBT_sdef$Median_Log2FC <-
        apply(results_LLvsBT_sdef[, grep("^logFC", colnames(results_LLvsBT_sdef))], MARGIN = 1, FUN = median)

# Absolute Median log2FC
results_LLvsBT_sdef$`|medianLogFC|` <-
        abs(results_LLvsBT_sdef$Median_Log2FC)

# Number of DEG with FDR < 0.1 (10%)
results_LLvsBT_sdef$Significant <-
        apply(results_LLvsBT_sdef[, grep("^adj.P.Val", colnames(results_LLvsBT_sdef))], 1, function(x)
                sum(abs(x < 0.1)))

# Unicode Arrows (symbols) indicating up- or down-regulation
results_LLvsBT_sdef$Direction <-
        as.character(apply(
                apply(
                        results_LLvsBT_sdef[, grep("^logFC", colnames(results_LLvsBT_sdef))],
                        MARGIN = 2,
                        FUN = function(x) {
                                ifelse(x < 0, "\U2193", "\U2191")
                        }
                ),
                1,
                paste,
                collapse = " "
        ))

# Ordering columns
colnames(results_LLvsBT_sdef)
results_LLvsBT_sdef <- results_LLvsBT_sdef[, c("ENTREZID", "SYMBOL", "GENENAME", "logFC_GSE24280_LLvsBT", "adj.P.Val._GSE24280_LLvsBT", "logFC_GSE17763_LLvsBT", "adj.P.Val._GSE17763_LLvsBT", "logFC_GSE74481_LLvsBT", "adj.P.Val._GSE74481_LLvsBT", "logFC_GSE443_LLvsBT", "adj.P.Val._GSE443_LLvsBT", "Median_Log2FC", "|medianLogFC|", "Significant", "Direction" )]

## Exporting *.XLS file --------------------------------

WriteXLS(x = results_LLvsBT_sdef, ExcelFileName = "Results_sdef_LLvsBT.xls")

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

## End program -----------------
rm(list=ls())
q(save="no")
