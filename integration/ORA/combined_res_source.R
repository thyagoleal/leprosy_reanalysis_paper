## GO Enrichment Analysis from meta-analysis/integration results
## Thyago Leal Calvo - thyagoleal@yahoo.com
## Last updated: June 17th, 2019
## ********************************

## Settings --------------
set.seed(12909)
options(digits = 3)
dir()
root <- getwd()

library(clusterProfiler)
library(readxl)
library(WriteXLS)
library(org.Hs.eg.db)
library(ggplot2)
source("src/Functions.R")

## Reading in data ---------------------

## LL vs BT *************

## REM, only genes with mean ES combined >= 1 & FDR < 0.1 are kept

llvsbt_REM <- as.data.frame(read_xls("../meta-analysis/LLvsBT/results/LLvsBT_ES.xls"))
llvsbt_REM_sub <- base::subset(llvsbt_REM, abs(muhat)>= 1 & FDR < 0.1)

## LL vs R2  *************

## REM, only genes with mean ES combined >= 1 & FDR < 0.1 are kept

llvsr2_REM <- as.data.frame(read_xls("../meta-analysis/LLvsR2/results/llvsr2_ES.xls"))
llvsr2_REM_sub <- base::subset(llvsr2_REM, abs(muhat)>= 1 & FDR < 0.01)

## LL vs Control *************

# sdef, only genes with concordant expression direction

llvscon_sdef <- as.data.frame(read_xls("../sdef/LLvsControl/Results_sdef_LLvsCon.xls"))
llvscon_sdef$Direction2 <- as.integer(factor(llvscon_sdef$Direction))
llvscon_sdef_sub <- llvscon_sdef[llvscon_sdef$Direction2 %in% c(1, 4),]

llvscon_universe <- readRDS("../sdef/LLvsControl/objects/LLvsCon_Pvals.Rds")
llvscon_universe <- as.character(row.names(llvscon_universe))

## In vitro studies *************

# sdef, only genes with concordant expression direction

invitro_sdef <- as.data.frame(read_xls("../sdef/invitro/Final.DEG.invitro.xls"))
invitro_sdef <- invitro_sdef[invitro_sdef$`|medianLogFC|` >= 0.5, ]

invitro_univ <- readRDS("../sdef/invitro/objects/Invitro_Pvals.Rds")
invitro_univ <- as.character(row.names(invitro_univ))

## Enrichment analysis ---------------------------

## LLvsBT ## 

llvsbt.go <-
        enrichGO(
               as.character(llvsbt_REM_sub$Row.names), # only DE REM genes
                OrgDb = "org.Hs.eg.db",
                ont = "BP",
                universe = as.character(llvsbt_REM$Row.names), # all REM genes
                pvalueCutoff = 0.1,
                pAdjustMethod = "fdr",
                qvalueCutoff = 0.1,
                readable = TRUE,
                pool = FALSE)

# z-score 
llvsbt.go@result$Zscore <- gozscore2(llvsbt.go, llvsbt_REM)

llvsbt.go@result$ZscoreSign <- ifelse(llvsbt.go@result$Zscore >= 1, "Mostly up-reg.", ifelse(llvsbt.go@result$Zscore <= -1, "Mostly down-reg.", "Up/Down reg."))

#llvsbt.go@result <- llvsbt.go@result[order(llvsbt.go@result$qvalue, decreasing = TRUE), ]

#llvsbt.go <- clusterProfiler::simplify(llvsbt.go, cutoff = 0.7, by = "p.adjust")

dotplot(llvsbt.go, showCategory = 20, size = NULL, x = "Count") + 
        ggplot2::facet_wrap(~llvsbt.go@result$ZscoreSign[1:20]) +
        ggplot2::scale_color_viridis_c(option = "C")

dev.copy2eps(
        device = "postscript",
        file = "results/combined/llvsbt_REM.eps",
        family = "ArialMT",
        colormodel = "cmyk",
        width = 8.5,
        height = 4.5
)

# subset desired ontologies for heatplot
llvsbt.go2 <- llvsbt.go

llvsbt.go2@result <-
        llvsbt.go2@result[llvsbt.go2@result$ID %in% c(
                "GO:0045087",
                "GO:0050900",
                "GO:0032640",
                "GO:0002685",
                "GO:0006935",
                "GO:0030593",
                "GO:0034340", 
                "GO:0002683", 
                "GO:0055094", 
                "GO:0009913",
                "GO:0060337",
                "GO:0030216",
                "GO:0043588",
                "GO:0031424",
                "GO:0070268"
        ),]

enrichplot::heatplot(llvsbt.go2, showCategory = 15, foldChange = structure(llvsbt_REM$muhat, names = llvsbt_REM$Row.names)) + ggplot2::scale_fill_continuous(type = "viridis", option = "C", values = c(0, 0.45, 1))

dev.copy2pdf(
        device = "pdf", 
        file = "results/combined/llvsbt_heatplot.pdf", 
        family = "ArialMT", 
        colormodel = "cmyk", 
        width = 9, 
        height = 2
)

cnetplot(llvsbt.go2, foldChange = structure(llvsbt_REM$muhat, names = llvsbt_REM$Row.names), showCategory = 14)

dev.copy2pdf(
        device = "pdf",
        file = "results/combined/llvsbt_cnetplot.pdf",
        family = "ArialMT",
        colormodel = "cmyk",
        width = 9,
        height = 9
)

WriteXLS(llvsbt.go@result, ExcelFileName = "results/combined/LLvsBT_enrichment.xls")

rm(list = ls(pattern = "llvsbt.*"))

## LLvsR2 --------------

llvsr2.go <-
        enrichGO(
                as.character(llvsr2_REM_sub$Row.names),
                universe = as.character(llvsr2_REM$Row.names),
                OrgDb = "org.Hs.eg.db",
                ont = "BP",
                pvalueCutoff = 0.1,
                pAdjustMethod = "fdr",
                qvalueCutoff = 0.1,
                readable = TRUE,
                pool = FALSE
        )

# z-score 
llvsr2.go@result$Zscore <- gozscore2(llvsr2.go, llvsr2_REM)

llvsr2.go@result$ZscoreSign <- ifelse(llvsr2.go@result$Zscore >= 1, "Mostly up-reg.", ifelse(llvsr2.go@result$Zscore <= -1, "Mostly down-reg.", "Up/Down reg."))

dotplot(llvsr2.go, showCategory = 9, x = "Count") + ggplot2::facet_wrap(~llvsr2.go@result$ZscoreSign[1:9])+ ggplot2::scale_color_viridis_c(option = "C")

dev.copy2eps(
        device = "postscript",
        file = "results/combined/llvsr2_REM.eps",
        family = "ArialMT",
        colormodel = "cmyk",
        width = 7,
        height = 3.2
)

# subset desired ontologies for heatplot
llvsr2.go2 <- llvsr2.go

llvsr2.go2@result <-
        llvsr2.go2@result[llvsr2.go2@result$ID %in% c(
                "GO:0043312",
                "GO:0042119",
                "GO:0036230",
                "GO:0048771",
                "GO:0055072",
                "GO:0042742",
                "GO:0090382",
                "GO:0051452",
                "GO:0019730",
                "GO:0030198",
                "GO:0006959",
                "GO:0015991"
        ),]

enrichplot::heatplot(llvsr2.go2, showCategory = 6, foldChange = structure(llvsr2_REM$muhat, names = llvsr2_REM$Row.names)) + ggplot2::scale_fill_continuous(type = "viridis", option = "C", values = c(0, 0.45, 1))

dev.copy2pdf(
        device = "pdf",
        file = "results/combined/llvsr2_heatplot.pdf",
        family = "ArialMT",
        colormodel = "cmyk",
        width = 6,
        height = 2
)

cnetplot(llvsr2.go2, foldChange = structure(llvsr2_REM$muhat, names = llvsr2_REM$Row.names), showCategory = 8)
dev.copy2pdf(
        device = "pdf",
        file = "results/combined/llvsr2_cnetplot.pdf",
        family = "ArialMT",
        colormodel = "cmyk",
        width = 8,
        height = 7
)

WriteXLS(llvsr2.go@result, ExcelFileName = "results/combined/LLvsR2_comb.xls")

rm(list = ls(pattern = "llvsr2.*"))

## LLvsCon -----------------

llvscon.go <-
        enrichGO(
                as.character(llvscon_sdef_sub$ENTREZID),
                universe = llvscon_universe,
                OrgDb = "org.Hs.eg.db",
                ont = "BP",
                pvalueCutoff = 0.1,
                pAdjustMethod = "fdr",
                qvalueCutoff = 0.1,
                readable = TRUE,
                pool = FALSE
        )

llvscon.go_simp <- simplify(llvscon.go)

## Abbreviate long strings

llvscon.go_simp@result["GO:0002460", "Description"] <- "Adaptative immune response based on somatic recomb."
llvscon.go_simp@result["GO:0002822", "Description"] <- "Regulation of adaptive immune response based on somatic recomb."

dotplot(llvscon.go_simp, x = "Count", showCategory = 25) + ggplot2::scale_color_viridis_c(option = "C")

dev.copy2eps(
        device = "postscript",
        file = "results/combined/llvscon_sdef.eps",
        family = "ArialMT",
        colormodel = "cmyk",
        width = 9,
        height = 5
)

## getting mean log2fold change between two datasets to use in heatplot

exp1 <- read_xls("../GSE74481/results/results.LLvsHealthy.xls")
exp1 <- exp1[, c("ENTREZID", "logFC")]

exp2 <- read_xls("../GSE24280/results/TissueMBvsHealthy.xls")
exp2 <- exp2[, c("ENTREZID", "logFC")]

llvscon_exp <- merge(exp1, exp2, by = "ENTREZID")
llvscon_exp <- llvscon_exp[llvscon_exp$ENTREZID %in% llvscon_sdef_sub$ENTREZID, ]

llvscon_exp <- structure(rowMeans(llvscon_exp[, c("logFC.x", "logFC.y")]), names = as.character(llvscon_exp$ENTREZID))
rm(exp1, exp2)

## choose ontologies for the heatplot
llvscon.go2 <- llvscon.go_simp
llvscon.go2@result <-
        llvscon.go2@result[llvscon.go2@result$ID %in% c(
                "GO:0071346",
                "GO:0008544",
                "GO:0043312",
                "GO:0042110",
                "GO:0006007",
                "GO:0002209", 
                "GO:0006897", 
                "GO:0032365",
                "GO:0010872",
                "GO:0001837",
                "GO:0048863",
                "GO:0030216"
        ),]

heatplot(llvscon.go2, showCategory = 12, foldChange = llvscon_exp) +
        ggplot2::scale_fill_continuous(type = "viridis", option = "C", values = c(0, 0.45, 1))

dev.copy2pdf(
        device = "pdf",
        file = "results/combined/llvscon_heatplot.pdf",
        family = "ArialMT",
        colormodel = "cmyk",
        width = 9,
        height = 3
)

cnetplot(llvscon.go2, foldChange = llvscon_exp, showCategory = 10)

dev.copy2pdf(
        device = "pdf",
        file = "results/combined/llvscon_cnetplot.pdf",
        family = "ArialMT",
        colormodel = "cmyk",
        width = 8,
        height = 7
)

WriteXLS(llvscon.go_simp@result, ExcelFileName = "results/combined/LLvsCon_comb.xls")

rm(list = ls(pattern = "llvscon.*"))

## In vitro studies ---------------

invitro.go <-
        enrichGO(
               as.character(invitro_sdef$ENTREZID),
                universe = invitro_univ,
                OrgDb = "org.Hs.eg.db",
                ont = "BP",
                pvalueCutoff = 0.1,
                pAdjustMethod = "fdr",
                qvalueCutoff = 0.1,
                readable = TRUE,
                pool = FALSE
        )

# renaming long string
invitro.go@result["GO:0070059", "Description"] <- "intrinsic apoptotic sig. path. in response to endo. reticulum stress"

dotplot(invitro.go, showCategory = 13, x = "Count") + ggplot2::scale_color_viridis_c(option = "C")

dev.copy2eps(
        device = "postscript",
        file = "results/combined/invitro_sdef.eps",
        family = "ArialMT",
        colormodel = "cmyk",
        width = 8,
        height = 4
)

# read mean log2fc

invitro_eff <- readRDS("../sdef/invitro/objects/Invitro_log2FC.Rds")

invitro_foldchange <- apply(invitro_eff, MARGIN = 1, FUN = median)
head(invitro_foldchange)

heatplot(invitro.go, showCategory = 13, foldChange = invitro_foldchange) + ggplot2::scale_fill_continuous(type = "viridis", option = "C", values = c(0, 0.45, 1))

dev.copy2pdf(
        device = "pdf",
        file = "results/combined/invitro_heatplot.pdf",
        family = "ArialMT",
        colormodel = "cmyk",
        width = 9,
        height = 2.5
)

cnetplot(invitro.go, showCategory = 13, foldChange = invitro_foldchange)

dev.copy2pdf(
        device = "pdf",
        file = "results/combined/invitro_cnetplot.pdf",
        family = "ArialMT",
        colormodel = "cmyk",
        width = 8,
        height = 7.5
)

rm(list = ls(pattern = "^invitro."))

WriteXLS(invitro.go@result, ExcelFileName = "results/combined/invitro.xls")

## End
setwd(root)
rm(list = ls())
q(save = "no")
