## Meta-analysis comparing tools for gene expression data from published datasets 
## Author: Thyago Leal Calvo - thyagoleal@yahoo.com                        
## Last updated: June 14th 2019
## ********************************

## Settings ----------------------------

set.seed(13130)
options(digits = 3)
options(download.file.method = "wget")

library("MetaDE")
library("readxl")
library("org.Hs.eg.db")
library("AnnotationDbi")
source("src/source.R")
library("WriteXLS")
library("VennDiagram")

## Read Expression Matrix from individual analysis step #########
# This should be already normalized and in log2 space

# Since the study GSE24280 only contains one sample in one group it is impossible to 
# calculate the effect size adequately. Therefore, for this context only the p-value-based
# methods are going to be employed. 

## Getting P-values calculated in limma ----------------------------

p_gse24280 <- as.data.frame(read_xls("../../GSE24280/results/TissueMBvsTissueBT.xls"))
p_gse24280 <- p_gse24280[,c("ENTREZID", "P.Value", "adj.P.Val")]
p_gse24280 <- p_gse24280[!is.na(p_gse24280$ENTREZID),]
row.names(p_gse24280) <- p_gse24280$ENTREZID

p_gse74481 <- as.data.frame(read_xls("../../GSE74481/results/results.LLvsHealthy.xls"))
p_gse74481 <- p_gse74481[,c("ENTREZID", "P.Value", "adj.P.Val")]
p_gse74481 <- p_gse74481[!is.na(p_gse74481$ENTREZID),]
row.names(p_gse74481) <- p_gse74481$ENTREZID

## Genes common to all studies
i <- Reduce(intersect, list(row.names(p_gse24280), row.names(p_gse74481)))

# remove ribosomal genes
rib <- read.table("../ribosomal_proteins.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
head(rib, 2)

i <- setdiff(i, rib$GeneID)

tmp <- cbind(p_gse24280[i,]$P.Value, p_gse74481[i,]$P.Value)
row.names(tmp) <- i
colnames(tmp) <- c("GSE24280","GSE74481")

rm(i)

## List holding p-values for all studies
llvscon_pvals <- list(p = tmp)
rm(tmp)

## Meta-analyses using different algorithms ----------------------

llvscon.mul <-
        MetaDE.pvalue(
                llvscon_pvals,
                meta.method = c("roP","maxP", "SR"),
                rth = 2, 
                parametric = FALSE)

class(llvscon.mul) <- "MetaDE.pvalue"

## Merging multiple meta-analyses results for visualization -----------------------

pdf("fig/DEGs_Meta.pdf", width = 6, height = 6)
draw.DEnumber(llvscon.mul, maxcut = 0.1, FDR = TRUE)
dev.off()

## Grabbing FDR according to all methods and saving ------------------------

## Grabbing study specific FDR 

i <- Reduce(intersect, list(row.names(p_gse24280), row.names(p_gse74481)))
# remove ribosomal genes

i <- setdiff(i, rib$GeneID)

tmp <- cbind(p_gse24280[i,]$adj.P.Val, p_gse74481[i,]$adj.P.Val)
row.names(tmp) <- i
colnames(tmp) <- c("GSE24280", "GSE74481")

res <- cbind(tmp[i, ], llvscon.mul$meta.analysis$FDR[i, ])
rm(tmp)

res <- as.data.frame(res)
res$Symbol <-
        select(
                org.Hs.eg.db,
                keys = row.names(res),
                columns = "SYMBOL",
                "ENTREZID"
        )$SYMBOL 

# Saving results from other methods

WriteXLS(res, ExcelFileName = "results/llvscon_allMethods.xls", row.names = TRUE)

# # NOT RUN
# ## Venn Diagram by meta-analysis methods ------------------------------
# 
# # Sources
# 
# venn.diagram(
#         x = list(
#                 SDEF = sdef,
#                 roP = rop,
#                 maxP = maxp,
#                 SR = sr
#         ),
#         imagetype = "tiff",
#         resolution = 500,
#         filename = "fig/Venn_llvscon.tiff",
#         col = "black",
#         width = 3000,
#         height = 3000,
#         fill = c("firebrick1", "dodgerblue", "yellow", "limegreen"),
#         alpha = rep(.3, 4),
#         ext.line.lwd = 2,
#         scaled = FALSE,
#         ext.text = TRUE,
#         lwd = rep(2, 4),
#         force.unique = TRUE,
#         units = 'px',
#         cex = 0.9,
#         cat.fontface = "bold",
#         margin = 0.1,
#         cat.cex = 1,
#         main.cex = 1.1,
#         na = "remove",
#         fontfamily = "Arial",
#         cat.fontfamily = "Arial",
#         main.fontfamily = "Arial",
#         at.default.pos = "outer"
# )
#  
# file.remove(list.files(path = "fig/", pattern = "*.log$", full.names = TRUE))
 
# ## saving the intersection between all methods --------------

sdef <- read_xls("../../sdef/LLvsControl/Results_sdef_LLvsCon.xls")
sdef <- as.character(sdef$ENTREZID)

rop <- as.character(row.names(res[res$roP < 0.1,]))

maxp <- as.character(row.names(res[res$maxP < 0.1,]))

sr <- as.character(row.names(res[res$SR < 0.1,]))

intersection <- data.frame(
        ENTREZID = Reduce(intersect,
                        list(
                                sdef,
                                rop,
                                maxp
                        )))

intersection$SYMBOL <- select(org.Hs.eg.db, keys = as.character(intersection$ENTREZID), columns = "SYMBOL", "ENTREZID")$SYMBOL
 
WriteXLS(intersection, ExcelFileName = "results/Intersect_fdr0.1.xls")

### END 
rm(list = ls())
q(save = "no", status = 0L)
