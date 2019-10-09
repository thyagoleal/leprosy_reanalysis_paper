## sdef - synthesize lists of significant features in related experiments - LLvsR2
## Last updated: June 14th, 2019.
## Thyago Leal Calvo, thyagoleal@yahoo.com
## ******************************************

## Settings ------------

options(digits = 3, scipen = 999)
dir()
root <- getwd()
set.seed(1918)
require(sdef)
require(readxl)
require(dplyr)
require(ggplot2)
require(AnnotationDbi)
require(org.Hs.eg.db)
require(WriteXLS)

## LL vs R2 ------------
## datasets used: 16844, 74481

GSE16844_LLvsR2 <- as.data.frame(read_xls("../GSE16844/results/results_GSE16844_LLvsENL.xls"))

GSE74481_LLvsR2 <-  as.data.frame(read_xls("../GSE74481/results/results.LLvsR2.xls"))

tmp.LLvsR2 <- list(GSE16844_LLvsR2, GSE74481_LLvsR2)
names(tmp.LLvsR2) <- c("GSE16844_LLvsR2", "GSE74481_LLvsR2")

mylist <- list()

for(i in names(tmp.LLvsR2)) {
        df <- tmp.LLvsR2[[i]]
        df <- df[!duplicated(df$ENTREZID), ]
        df <- df[!is.na(df$ENTREZID),]
        mylist[[i]] <- df
}

rm(i, df)

merged <-
        Reduce(function(x, y)
                merge.data.frame(x, y, all = FALSE, by = "ENTREZID"),
               mylist)

merged <- merged[, c(1, grep("P.Value.", colnames(merged)))]
row.names(merged) <- as.character(merged$ENTREZID)
colnames(merged)[2:3] <- names(mylist)
merged <- merged[, 2:3]

# Remove ribosomal genes 
rib <- read.table("ribosomal_proteins.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
head(rib, 5)
merged <- merged[!row.names(merged) %in% as.character(rib$GeneID), ]

setwd("LLvsR2/")

Th_LLvsR2 <- ratio(data = merged, pvalue = TRUE)
Rh_LLvsR2 <- baymod(iter = 500, output.ratio = Th_LLvsR2)
MC_LLvsR2 <- Tmc(iter = 200, output.ratio = Th_LLvsR2)
feat.names <- row.names(merged)
feat.lists.Bayesian <-
        extractFeatures.R(output.ratio = Th_LLvsR2,
                          output.bay = Rh_LLvsR2,
                          feat.names = feat.names)

createTable(output.ratio = Th_LLvsR2, output.bay = Rh_LLvsR2)

annt_LLvsR2 <-
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

merged <- merged[, c(1, grep("^adj.P.Val.|^logFC.*", colnames(merged)))]
colnames(merged)[-1] <- paste(rep(c("logFC", "adj.P.Val."), times = 2), rep(names(mylist), each = 2), sep = "_")
row.names(merged) <- as.character(merged$ENTREZID)
results_LLvsR2_sdef <- merge.data.frame(merged[annt_LLvsR2$ENTREZID, ], annt_LLvsR2, by.x = "row.names", by.y = "ENTREZID") 

names(results_LLvsR2_sdef)

## Summary Statistics across comparisons --------------------------

# Median log2FC
results_LLvsR2_sdef$Median_Log2FC <-
        apply(results_LLvsR2_sdef[, grep("^logFC", colnames(results_LLvsR2_sdef))], MARGIN = 1, FUN = median)

# Absolute Median log2FC
results_LLvsR2_sdef$`|medianLogFC|` <-
        abs(results_LLvsR2_sdef$Median_Log2FC)

# Number of DEG with FDR < 0.1 (10%)
results_LLvsR2_sdef$Significant <-
        apply(results_LLvsR2_sdef[, grep("^adj.P.Val", colnames(results_LLvsR2_sdef))], 1, function(x)
                sum(abs(x < 0.1)))

# Unicode Arrows (symbols) indicating up- or down-regulation
results_LLvsR2_sdef$Direction <-
        as.character(apply(
                apply(
                        results_LLvsR2_sdef[, grep("^logFC", colnames(results_LLvsR2_sdef))],
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
colnames(results_LLvsR2_sdef)
results_LLvsR2_sdef <- results_LLvsR2_sdef[, c(2, 7, 8, 3:6, 9:12)]

## Exporting *.XLS file --------------------------------

WriteXLS(results_LLvsR2_sdef, ExcelFileName = "Results_sdef_LLvsR2.xls")

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

## End program-----------
setwd(root)
rm(list=ls())
q(save = "no")
