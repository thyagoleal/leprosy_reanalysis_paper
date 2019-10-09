## sdef - synthesize lists of significant features in related experiments - LL vs In vitro
## Last update: June 14th, 2019.
## Thyago Leal Calvo, thyagoleal@yahoo.com
## ****************************************

## Settings ------------

dir()
options(digits = 3, scipen = 999)
root <- getwd()
set.seed(1918)
require(sdef)
require(readxl)
require(dplyr)
require(ggplot2)
require(AnnotationDbi)
require(org.Hs.eg.db)
require(WriteXLS)

### In vitro Experiments -------------------------------------

GSE35423_Inf24vsCon24 <-
        as.data.frame(read_xlsx("../GSE35423/results/inf48vsInf24.xlsx"))
names(GSE35423_Inf24vsCon24)[2] <- "ENTREZID"

GSE35423_Inf48vsCon48 <-
        as.data.frame(read_xls("../GSE35423/results/inf48vsCon48.xls"))
names(GSE35423_Inf48vsCon48)[2] <- "ENTREZID"

GSE100853_StimulatedvsControl <-
        as.data.frame(read_xls("../GSE100853/results/results_GSE100583.xls"))
names(GSE100853_StimulatedvsControl)

GSE95748_inf14d <-
        as.data.frame(read_xls(
                "../GSE95748/results/Results_GSE95748_inf14d_human_ids.xls"
        ))
names(GSE95748_inf14d)

GSE95748_in28d <-
        as.data.frame(read_xls(
                "../GSE95748/results/Results_GSE95748_inf28d_human_ids.xls"
        ))
names(GSE95748_in28d)

tmp.invitro <-
        list(
                GSE35423_Inf24vsCon24,
                GSE35423_Inf48vsCon48,
                GSE100853_StimulatedvsControl,
                GSE95748_inf14d,
                GSE95748_in28d
        )
names(tmp.invitro) <-
        c(
                "GSE35423_Inf24vsCon24",
                "GSE35423_Inf48vsCon48",
                "GSE100853_StimulatedvsControl",
                "GSE95748_inf14d",
                "GSE95748_in28d"
        )

mylist <- list()

for (i in names(tmp.invitro)) {
        df <- tmp.invitro[[i]]
        df <- df[!duplicated(df$ENTREZID), ]
        df <- df[!is.na(df$ENTREZID), ]
        mylist[[i]] <- df
}

rm(i, df)

merged <-
        Reduce(function(x, y)
                merge.data.frame(x, y, all = FALSE, by = "ENTREZID"),
               mylist)

merged <- merged[, c(1, grep("P.Value.*", colnames(merged)))]
row.names(merged) <- as.character(merged$ENTREZID)
colnames(merged)[2:6] <- names(mylist)
merged <- merged[, 2:6]

# Remove ribosomal genes 
rib <- read.table("ribosomal_proteins.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
head(rib, 5)
merged <- merged[!row.names(merged) %in% as.character(rib$GeneID), ]

## save pvals

if (!file.exists("invitro/objects/Invitro_Pvals.Rds")) {
        saveRDS(merged, file = "invitro/objects/Invitro_Pvals.Rds", compress = "bzip2")
        cat("Rds file saved!\n")
}

setwd("invitro/")

Th_models <- ratio(data = merged,
                   pvalue = TRUE,
                   interval = 0.05)

Rh_models <- baymod(iter = 500, output.ratio = Th_models)

MC_models <- Tmc(iter = 200, output.ratio = Th_models)
feat.names <- row.names(merged)
feat.lists.Bayesian <-
        extractFeatures.R(
                output.ratio = Th_models,
                output.bay = Rh_models,
                feat.names = feat.names,
                h = NULL
        )

createTable(output.ratio = Th_models, output.bay = Rh_models)

annt_models <-
        select(
                org.Hs.eg.db,
                keys = as.character(feat.lists.Bayesian$max$Names),
                columns = c("SYMBOL", "GENENAME"),
                keytype = "ENTREZID"
        )

## Grabbing Adjusted P-values ---------------------

merged <-
        Reduce(function(x, y)
                merge.data.frame(x, y, all = FALSE, by = "ENTREZID"),
               mylist)

merged <-
        merged[, c(1, grep("^adj.P.Val|^logFC.*", colnames(merged)))]
colnames(merged)[-1] <-
        paste(rep(c("logFC", "adj.P.Val."), times = 5), rep(names(mylist), each = 2), sep = "_")
row.names(merged) <- as.character(merged$ENTREZID)

results_models_sdef <-
        merge.data.frame(merged[annt_models$ENTREZID, ], annt_models, by.x = "row.names", by.y = "ENTREZID")

names(results_models_sdef)

## Summary Statistics across comparisons --------------------------

# Median log2FC
results_models_sdef$Median_Log2FC <-
        apply(results_models_sdef[, grep("^logFC", colnames(results_models_sdef))], MARGIN = 1, FUN = median)

# Absolute Median log2FC
results_models_sdef$`|medianLogFC|` <-
        abs(results_models_sdef$Median_Log2FC)

# Number of DEG with FDR < 0.1 (10%)
results_models_sdef$Significant <-
        apply(results_models_sdef[, grep("^adj.P.Val", colnames(results_models_sdef))], 1, function(x)
                sum(abs(x < 0.1)))

# Unicode Arrows (symbols) indicating up- or down-regulation
results_models_sdef$Direction <-
        as.character(apply(
                apply(
                        results_models_sdef[, grep("^logFC", colnames(results_models_sdef))],
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
colnames(results_models_sdef)
results_models_sdef <-
        results_models_sdef[, c(2, 13, 14, 3:12, 15:18)]

## Exporting *.XLS file --------------------------------

WriteXLS::WriteXLS(results_models_sdef, ExcelFileName = "Final.DEG.invitro.xls")

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

## End program -----------
setwd(root)
rm(list = ls())
quit(save = "no")
