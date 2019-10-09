## Meta-analysis (Invitro exp) comparing tools for gene expression data from published datasets 
## Author: Thyago Leal Calvo - thyagoleal@yahoo.com                        
## Last updated: June 14th, 2019                                                                
## ******************************

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

## Since this category has a lot of different studies with distincts stimuli types I'm only using the P-value-based tools. 

## Getting P-values calculated in limma ----------------------------

## GSE35423 48h ------------
p_gse35423_48h <- as.data.frame(read_xls("../../GSE35423/results/inf48vsCon48.xls"))
p_gse35423_48h <- p_gse35423_48h[,c("GENE", "P.Value", "adj.P.Val")]
p_gse35423_48h <- p_gse35423_48h[!is.na(p_gse35423_48h$GENE),]
p_gse35423_48h <- p_gse35423_48h[order(p_gse35423_48h$adj.P.Val, decreasing = FALSE), ]
p_gse35423_48h <- p_gse35423_48h[!duplicated(p_gse35423_48h$GENE),]
row.names(p_gse35423_48h) <- p_gse35423_48h$GENE

## GSE35423 24h ------------
p_gse35423_24h <- as.data.frame(read_xlsx("../../GSE35423/results/inf48vsInf24.xlsx"))
p_gse35423_24h <- p_gse35423_24h[,c("GENE", "P.Value", "adj.P.Val")]
p_gse35423_24h <- p_gse35423_24h[!is.na(p_gse35423_24h$GENE),]
p_gse35423_24h <- p_gse35423_24h[order(p_gse35423_24h$adj.P.Val, decreasing = FALSE), ]
p_gse35423_24h <- p_gse35423_24h[!duplicated(p_gse35423_24h$GENE),]
row.names(p_gse35423_24h) <- p_gse35423_24h$GENE

## GSE100853 Stimulated vs Con -------------------

p_gse100853 <- as.data.frame(read_xls("../../GSE100853/results/results_GSE100583.xls"))
p_gse100853 <- p_gse100853[,c("ENTREZID", "P.Value", "adj.P.Val")]
p_gse100853 <- p_gse100853[order(p_gse100853$adj.P.Val, decreasing = FALSE), ]
p_gse100853 <- p_gse100853[!duplicated(p_gse100853$ENTREZID), ]
row.names(p_gse100853) <- p_gse100853$ENTREZID

## GSE95748_inf14d --------------------

p_gse95748_inf14d <- as.data.frame(read_xls("../../GSE95748/results/Results_GSE95748_inf14d_human_ids.xls"))
p_gse95748_inf14d <- p_gse95748_inf14d[,c("ENTREZID", "P.Value", "adj.P.Val")]
p_gse95748_inf14d <- p_gse95748_inf14d[order(p_gse95748_inf14d$adj.P.Val,decreasing = FALSE), ]
p_gse95748_inf14d <- p_gse95748_inf14d[!duplicated(p_gse95748_inf14d$ENTREZID),]
p_gse95748_inf14d <- p_gse95748_inf14d[!is.na(p_gse95748_inf14d$ENTREZID), ]
row.names(p_gse95748_inf14d) <- p_gse95748_inf14d$ENTREZID

## GSE95748_inf28d --------------------

p_gse95748_inf28d <- as.data.frame(read_xls("../../GSE95748/results/Results_GSE95748_inf28d_human_ids.xls"))
p_gse95748_inf28d <- p_gse95748_inf28d[,c("ENTREZID", "P.Value", "adj.P.Val")]
p_gse95748_inf28d <- p_gse95748_inf28d[order(p_gse95748_inf28d$adj.P.Val,decreasing = FALSE), ]
p_gse95748_inf28d <- p_gse95748_inf28d[!duplicated(p_gse95748_inf28d$ENTREZID),]
p_gse95748_inf28d <- p_gse95748_inf28d[!is.na(p_gse95748_inf28d$ENTREZID), ]
row.names(p_gse95748_inf28d) <- p_gse95748_inf28d$ENTREZID

## Genes common to all studies
i <- Reduce(intersect, list(row.names(p_gse100853), row.names(p_gse35423_24h), row.names(p_gse35423_48h), row.names(p_gse95748_inf14d), row.names(p_gse95748_inf28d)))

# remove ribosomal genes
rib <- read.table("../ribosomal_proteins.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
head(rib, 2)

i <- setdiff(i, rib$GeneID)

tmp <- cbind(p_gse100853[i,]$P.Value, p_gse35423_24h[i,]$P.Value, p_gse35423_48h[i, ]$P.Value, p_gse95748_inf14d[i,]$P.Value,  p_gse95748_inf28d[i,]$P.Value)

row.names(tmp) <- i
colnames(tmp) <- c("GSE100853","GSE35423_24h", "GSE35423_48h", "GSE95748_14d", "GSE95748_28d")

rm(i)

## List holding p-values for all studies
invitro_pvals <- list(p = tmp)
rm(tmp)

## Meta-analyses using different algorithms ----------------------

invitro.mul <-
        MetaDE.pvalue(
                invitro_pvals,
                meta.method = c("roP","maxP", "SR"),
                rth = 4,
                parametric = FALSE)

class(invitro.mul) <- "MetaDE.pvalue"

## Merging multiple meta-analyses results for visualization -----------------------

pdf("fig/DEGs_Meta.pdf", width = 6, height = 6)
draw.DEnumber(result = invitro.mul, maxcut = 0.1, FDR = TRUE)
dev.off()

## Grabbing FDR according to all methods and saving ------------------------

## Grabbing study specific FDR 

i <- Reduce(intersect, list(row.names(p_gse100853), row.names(p_gse35423_24h), row.names(p_gse35423_48h), row.names(p_gse95748_inf14d), row.names(p_gse95748_inf28d)))

# remove ribosomal genes
i <- setdiff(i, rib$GeneID)

tmp <- cbind(p_gse100853[i,]$adj.P.Val, p_gse35423_24h[i,]$adj.P.Val, p_gse35423_48h[i, ]$adj.P.Val, p_gse95748_inf14d[i,]$adj.P.Val,  p_gse95748_inf28d[i,]$adj.P.Val)
row.names(tmp) <- i
colnames(tmp) <- c("GSE100853","GSE35423_24h", "GSE35423_48h", "GSE95748_14d", "GSE95748_28d")

res <- cbind(tmp[i, ], invitro.mul$meta.analysis$FDR[i, ])
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

WriteXLS(res, ExcelFileName = "results/invitro_allMethods.xls", row.names = TRUE)
 
# # NOT RUN
# ## Venn Diagram by meta-analysis methods ------------------------------
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
#         filename = "fig/Venn_invitro.tiff",
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
# file.remove(list.files(
#         path = "fig/",
#         pattern = "*.log$",
#         full.names = TRUE
# ))
 
## saving the intersection between all methods --------------

# Sources

sdef <- read_xls("../../sdef/invitro/results/Final.DEG.invitro.xls")
sdef <- as.character(sdef$ENTREZID)

rop <- row.names(res[res$roP < 0.1,])

maxp <- row.names(res[res$maxP < 0.1,])

sr <- row.names(res[res$SR < 0.1,])
 
intersection <- data.frame(
         ENTREZID = Reduce(intersect,
                        list(
                                sdef,
                                rop,
                                maxp, sr
                        )))
 
intersection$SYMBOL <- select(org.Hs.eg.db, keys = as.character(intersection$ENTREZID), columns = "SYMBOL", "ENTREZID")$SYMBOL
 
WriteXLS(intersection, ExcelFileName = "results/Intersect_fdr0.1_ALL.xls")
 
## saving the intersection between all methods except SR --------------
 
intersection <- data.frame(
        ENTREZID = Reduce(intersect,
                          list(
                                  sdef,
                                  rop,
                                  maxp
                          )))

intersection$SYMBOL <- select(org.Hs.eg.db, keys = as.character(intersection$ENTREZID), columns = "SYMBOL", "ENTREZID")$SYMBOL

WriteXLS(intersection, ExcelFileName = "results/Intersect_fdr0.1_exceptSR.xls")

### END 
rm(list = ls())
q(save = "no", status = 0L)
