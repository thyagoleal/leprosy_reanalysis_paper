## Meta-analysis comparing tools for gene expression data from published datasets 
## Author: Thyago Leal Calvo - thyagoleal@yahoo.com                        
## Last updated: June 14th 2019                                                                

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

### GSE24280 ------------------------------

gse24280 <- read.table("../../GSE24280/results/exp_24280.txt", header = TRUE, sep = "\t")
head(gse24280, 2)

## Formatting

table(duplicated(gse24280$ENTREZ))

row.names(gse24280) <- as.character(gse24280$ENTREZ)
gse24280 <- gse24280[,c(13:24)]
head(gse24280, 2)
label_24280 <- ifelse(grepl("^MB.", colnames(gse24280)) == TRUE, '1', '0')
label_24280

## GSE17763  ------------------------------

gse17763 <- read.table("../../GSE17763/results/exp_matrix_gse17763.txt", header = TRUE, sep = "\t")
head(gse17763, 2)

## Formatting 

table(duplicated(gse17763$ENTREZ))

row.names(gse17763) <- as.character(gse17763$ENTREZ)
gse17763 <- gse17763[, c(2:18)]
head(gse17763, 2)
label_17763 <- ifelse(grepl("LL", colnames(gse17763)) == TRUE, '1', '0') 
label_17763

## GSE74481 ------------------------------

gse74481 <- read.table("../../GSE74481/results/exp_74481.txt", header = TRUE, sep = "\t")

table(duplicated(gse74481$ENTREZ))

row.names(gse74481) <- as.character(gse74481$ENTREZ)
gse74481 <- gse74481[, grep("^LL|^BT", colnames(gse74481))]
head(gse74481, 2)
label_74481 <- ifelse(grepl("LL", colnames(gse74481)) == TRUE, '1', '0')
label_74481

## transforming ratios to positive expression values by adding a constant
gse74481 <- gse74481 + abs(floor(min(gse74481)))

## GSE443 ------------------------------

gse443 <- read.table("../../GSE443/results/exp_matrix_gse443.txt", header = TRUE, sep = "\t")

table(duplicated(gse443$entrez))
names(gse443)
row.names(gse443) <- as.character(gse443$entrez)
gse443 <- gse443[, grep("^LL|^BT", colnames(gse443))]
head(gse443, 2)
label_443 <- ifelse(grepl("LL", colnames(gse443)) == TRUE, '1', '0')
label_443

## Constructing the list -----------------------------

i <- Reduce(intersect, x = list(row.names(gse24280), rownames(gse17763), row.names(gse74481), row.names(gse443)))

## exclude genes from ribosomal proteins
rib <- read.table("../ribosomal_proteins.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
head(rib, 3)

i <- setdiff(i, rib$GeneID)

llvsbt <-
        list(
                data.matrix(gse24280[i, ]),
                data.matrix(gse17763[i, ]),
                data.matrix(gse74481[i, ]),
                data.matrix(gse443[i, ])
        )
rm(i)

names(llvsbt) <- c("GSE2480", "GSE17763", "GSE74481", "GSE443")

K <- length(llvsbt)

label <- list(GSE24280 = label_24280, GSE17763 = label_17763, GSE74481 = label_74481, GSE443 = label_443)

clin.data <- lapply(label, function(x) {data.frame(x)} )
for (k in 1:length(clin.data)){
        colnames(clin.data[[k]]) <- "label"
}

# meta analysis REM
meta.res <- MetaDE(
        data = llvsbt,
        clin.data = clin.data,
        data.type = "continuous",
        resp.type = "twoclass",
        response = 'label',
        ind.method = rep('limma', 4),
        meta.method = "REM",
        select.group = c('0', '1'),
        ref.level = c('0'),
        paired = rep(FALSE, length(llvsbt)),
        REM.type = "HO",
        tail = 'abs'
)

## saving main holder object 

saveRDS(meta.res, file = "obj/llvsbt_meta.rds", compress = "bzip2")

# Constructing a matrix to hold the results

llvsbt.rem2 <- matrix(data = c(
        meta.res$meta.analysis$mu.hat,
        meta.res$meta.analysis$mu.var,
        meta.res$meta.analysis$tau2,
        meta.res$meta.analysis$FDR),
        nrow = length(meta.res$meta.analysis$mu.hat),
        ncol = 4,
        dimnames = list(
                row.names(meta.res$meta.analysis$FDR), c("muhat", "muvar", "tau2", "FDR")), byrow = FALSE)

llvsbt.rem2 <- cbind(meta.res$ind.ES[row.names(llvsbt.rem2), ], llvsbt.rem2)

llvsbt.rem2 <- as.data.frame(llvsbt.rem2)
llvsbt.rem2$Symbol <- select(org.Hs.eg.db, keys = row.names(llvsbt.rem2),columns = c("SYMBOL"), keytype = "ENTREZID")$SYMBOL

colnames(llvsbt.rem2)[1:4] <- names(llvsbt)

llvsbt.rem2 <- merge(llvsbt.rem2, meta.res$ind.Var, by = "row.names")
colnames(llvsbt.rem2)[11:14] <- paste("var", names(llvsbt), sep = "_")

WriteXLS(llvsbt.rem2, ExcelFileName = "results/LLvsBT_ES.xls")

rm(list = setdiff(ls(), c("meta.res", 'draw.DEnumber_all')))

## Getting P-values calculated in limma ----------------------------

p_gse24280 <- as.data.frame(read_xls("../../GSE24280/results/TissueMBvsTissueBT.xls"))
p_gse24280 <- p_gse24280[,c("ENTREZID", "P.Value", "adj.P.Val")]
p_gse24280 <- p_gse24280[!is.na(p_gse24280$ENTREZID),]
row.names(p_gse24280) <- p_gse24280$ENTREZID

p_gse17763 <- as.data.frame(read_xls("../../GSE17763/results/results_LLvsBT.xls"))
p_gse17763 <- p_gse17763[,c("ENTREZID", "P.Value", "adj.P.Val")]
p_gse17763 <- p_gse17763[!is.na(p_gse17763$ENTREZID),]
row.names(p_gse17763) <- p_gse17763$ENTREZID

p_gse443 <- as.data.frame(read_xls("../../GSE443/results/DEG_full_results.xls"))
p_gse443 <- p_gse443[,c("ENTREZID", "P.Value", "adj.P.Val")]
p_gse443 <- p_gse443[!duplicated(p_gse443$ENTREZID),]
row.names(p_gse443) <- p_gse443$ENTREZID

p_gse74481 <- as.data.frame(read_xls("../../GSE74481/results/results.LLvsBT.xls"))
p_gse74481 <- p_gse74481[,c("ENTREZID", "P.Value", "adj.P.Val")]
p_gse74481 <- p_gse74481[!is.na(p_gse74481$ENTREZID),]
row.names(p_gse74481) <- p_gse74481$ENTREZID

## Genes common to all studies
i <- Reduce(intersect, list(row.names(p_gse24280), row.names(p_gse17763), row.names(p_gse74481), row.names(p_gse443)))

## exclude genes from ribosomal proteins
rib <- read.table("../ribosomal_proteins.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
i <- setdiff(i, rib$GeneID)

tmp <- cbind(p_gse24280[i,]$P.Value, p_gse17763[i,]$P.Value, p_gse74481[i,]$P.Value, p_gse443[i,]$P.Value)
row.names(tmp) <- i
colnames(tmp) <- c("GSE24280","GSE17763", "GSE74481", "GSE443")

rm(i)

## List holding p-values for all studies
llvsbt_pvals <- list(p = tmp)
rm(tmp)

## Meta-analyses using different algorithms ----------------------

llvsbt.mul <-
        MetaDE.pvalue(
                llvsbt_pvals,
                meta.method = c("roP","maxP", "SR"),
                rth = 3, parametric = FALSE)

## Merging multiple meta-analyses results for visualization ----------------

res <- cbind(llvsbt.mul$ind.p, llvsbt.mul$meta.analysis$pval)
tmp <- cbind(meta.res$meta.analysis$pval, meta.res$meta.analysis$FDR)
res <- cbind(res, tmp[row.names(res),][,1])
colnames(res)[8] <- "REM" 
rm(tmp)

pdf("fig/DEGs_Meta.pdf", width = 6, height = 6)
draw.DEnumber_all(result = res, maxcut = 0.1, FDR = TRUE, main = "Number of DEG by individual and meta-analyses methods \n(LL vs. BT)")
dev.off()
rm(res)

## Grabbing FDR according to all methods and saving ------------------------

## Grabbing study specific FDR 

i <- Reduce(intersect, list(row.names(p_gse24280), row.names(p_gse17763), row.names(p_gse74481), row.names(p_gse443)))

## exclude genes from ribosomal proteins
i <- setdiff(i, rib$GeneID)

tmp <- cbind(p_gse24280[i,]$adj.P.Val, p_gse17763[i,]$adj.P.Val,  p_gse74481[i,]$adj.P.Val, p_gse443[i,]$adj.P.Val)

row.names(tmp) <- i

colnames(tmp) <- c("GSE24280","GSE17763", "GSE74481", "GSE443")

res <- cbind(tmp[i, ], llvsbt.mul$meta.analysis$FDR[i, ])

res <- cbind(res[i, ], meta.res$meta.analysis$FDR[i,])
colnames(res)[8] <- "REM" 
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

WriteXLS(res, ExcelFileName = "results/LLvsBT_allMethods.xls", row.names = TRUE)

## NOT RUN! 

## Venn Diagram by meta-analysis methods ------------------------------
# 
# # Sources
# 
# sdef <- read_xls("../../sdef/LLvsBT/Results_sdef_LLvsBT.xls", sheet = 1)
# sdef <- as.character(sdef$ENTREZID)
# 
# rem <- as.character(row.names(res[res$REM < 0.1,]))
# 
# rop <- row.names(res[res$roP < 0.1,])
# maxp <- row.names(res[res$maxP < 0.1,])
# sr <- row.names(res[res$SR < 0.1,])
# 
# venn.diagram(
#          x = list(
#                  SDEF = sdef,
#                  REM = rem,
#                  roP = rop,
#                  maxP = maxp,
#                  SR = sr
#          ),
#          imagetype = "tiff",
#          resolution = 500,
#          filename = "fig/Venn_LLvsBT.tiff",
#          col = "black",
#          width = 3000,
#          height = 3000,
#          fill = c("firebrick1", "dodgerblue", "yellow", "limegreen", "grey"),
#          alpha = rep(.3, 5),
#          ext.line.lwd = 2,
#          scaled = FALSE,
#          ext.text = TRUE,
#          lwd = rep(2, 5),
#          force.unique = TRUE,
#          units = 'px',
#          cex = 0.9,
#          cat.fontface = "bold",
#          margin = 0.1,
#          cat.cex = 1,
#          main.cex = 1.1,
#          na = "remove",
#          fontfamily = "Arial",
#          cat.fontfamily = "Arial",
#          main.fontfamily = "Arial",
#          at.default.pos = "outer"
# )
# 
# file.remove(list.files(path = "fig/", pattern = "*.log$", full.names = TRUE))

# ## saving the intersection between all methods --------------

intersection <- data.frame(
         ENTREZID = Reduce(intersect,
                         list(
                                 sdef,
                                 rem,
                                 rop,
                                 maxp,
                                 sr
                         )))

intersection$SYMBOL <- select(org.Hs.eg.db, keys = as.character(intersection$ENTREZID), columns = "SYMBOL", "ENTREZID")$SYMBOL

WriteXLS(intersection, ExcelFileName = "results/Intersect_fdr0.1.xls")

## saving the intersection between all methods except SR --------------

intersection <- data.frame(
        ENTREZID = Reduce(intersect,
                          list(
                                  sdef,
                                  rem,
                                  rop,
                                  maxp
                          )))

intersection$SYMBOL <- select(org.Hs.eg.db, keys = as.character(intersection$ENTREZID), columns = "SYMBOL", "ENTREZID")$SYMBOL

WriteXLS(intersection, ExcelFileName = "results/Intersect_fdr0.1_ExceptSR.xls")

### END
rm(list = ls())
q(save = "no", status = 0L)
