## Meta-analysis (LLvsR2) using different methods on gene expression data from published datasets
## Author: Thyago Leal Calvo - thyagoleal@yahoo.com                        
## Last updated: June 14th 2019
## **********************************

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

### GSE16844 ------------------------------

gse16844 <- read.table("../../GSE16844/results/exp_matrix_gse16844.txt", header = TRUE, sep = "\t")

## Formatting

table(duplicated(gse16844$ENTREZ))
gse16844 <- gse16844[!is.na(gse16844$ENTREZID),]

row.names(gse16844) <- as.character(gse16844$ENTREZ)
gse16844 <- gse16844[,c(4:ncol(gse16844))]
head(gse16844, 2)
label_16844 <- ifelse(grepl("LL", colnames(gse16844)) == TRUE, '1', '0')
label_16844

## GSE74481 ------------------------------

gse74481 <- read.table("../../GSE74481/results/exp_74481.txt", header = TRUE, sep = "\t")

table(duplicated(gse74481$ENTREZ))

row.names(gse74481) <- as.character(gse74481$ENTREZ)
gse74481 <- gse74481[, grep("^LL|^R2", colnames(gse74481))]
head(gse74481, 2)
label_74481 <- ifelse(grepl("LL", colnames(gse74481)) == TRUE, '1', '0')
label_74481

## transforming ratios to positive expression values by adding a constant (c = 8)
gse74481 <- gse74481 + abs(floor(min(gse74481)))

## Constructing the list -----------------------------

i <- Reduce(intersect, x = list(row.names(gse16844), row.names(gse74481)))

## removing ribosomal genes
rib <- read.table("../ribosomal_proteins.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
head(rib, 2)

i <- setdiff(i, rib$GeneID)

llvsr2 <-
        list(data.matrix(gse16844[i,]),
             data.matrix(gse74481[i, ]))

rm(i)

names(llvsr2) <- c("GSE16844", "GSE74481")

K <- length(llvsr2)

label <- list(GSE16844 = label_16844, GSE74481 = label_74481)

clin.data <- lapply(label, function(x) {data.frame(x)} )
for (k in 1:length(clin.data)){
        colnames(clin.data[[k]]) <- "label"
}

clin.data

# meta analysis REM
meta.res <- MetaDE(
        data = llvsr2,
        clin.data = clin.data,
        data.type = "continuous",
        resp.type = "twoclass",
        response = 'label',
        ind.method = rep('limma', 2),
        meta.method = "REM",
        select.group = c('0', '1'),
        ref.level = c('0'),
        paired = rep(FALSE, length(llvsr2)),
        REM.type = "HO",
        tail = 'abs'
)

## saving main holder object 

saveRDS(meta.res, file = "obj/llvsr2_meta.rds", compress = "bzip2")

# Constructing a matrix to hold the results

llvsr2.rem2 <- matrix(data = c(
        meta.res$meta.analysis$mu.hat,
        meta.res$meta.analysis$mu.var,
        meta.res$meta.analysis$tau2,
        meta.res$meta.analysis$FDR),
        nrow = length(meta.res$meta.analysis$mu.hat),
        ncol = 4,
        dimnames = list(
                row.names(meta.res$meta.analysis$FDR), c("muhat", "muvar", "tau2", "FDR")), byrow = FALSE)

llvsr2.rem2 <- cbind(meta.res$ind.ES[row.names(llvsr2.rem2), ], llvsr2.rem2)

llvsr2.rem2 <- as.data.frame(llvsr2.rem2)
llvsr2.rem2$Symbol <- select(org.Hs.eg.db, keys = row.names(llvsr2.rem2),columns = c("SYMBOL"), keytype = "ENTREZID")$SYMBOL

colnames(llvsr2.rem2)[1:2] <- names(llvsr2)

llvsr2.rem2 <- merge(llvsr2.rem2, meta.res$ind.Var, by = "row.names")
colnames(llvsr2.rem2)[9:10] <- paste("var", names(llvsr2), sep = "_")

WriteXLS(llvsr2.rem2, ExcelFileName = "results/llvsr2_ES.xls")

rm(list = setdiff(ls(), c("meta.res", 'draw.DEnumber_all')))

## Getting P-values calculated in limma ----------------------------

p_gse16844 <- as.data.frame(read_xls("../../GSE16844/results/results_GSE16844_LLvsENL.xls"))
p_gse16844 <- p_gse16844[,c("ENTREZID", "P.Value", "adj.P.Val")]
p_gse16844 <-p_gse16844[!is.na(p_gse16844$ENTREZID),]
row.names(p_gse16844) <- p_gse16844$ENTREZID

p_gse74481 <- as.data.frame(read_xls("../../GSE74481/results/results.LLvsR2.xls"))
p_gse74481 <- p_gse74481[,c("ENTREZID", "P.Value", "adj.P.Val")]
p_gse74481 <- p_gse74481[!is.na(p_gse74481$ENTREZID),]
row.names(p_gse74481) <- p_gse74481$ENTREZID

## Genes common to all studies
i <- Reduce(intersect, list(row.names(p_gse16844), row.names(p_gse74481)))

# removing ribosomal genes
rib <- read.table("../ribosomal_proteins.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
i <- setdiff(i, rib$GeneID)

tmp <- cbind(p_gse16844[i,]$P.Value, p_gse74481[i,]$P.Value)
row.names(tmp) <- i
colnames(tmp) <- c("GSE16844","GSE74481")

rm(i)

## List holding p-values for all studies
llvsr2_pvals <- list(p = tmp)
rm(tmp)

## Meta-analyses using different algorithms ----------------------

llvsr2.mul <-
        MetaDE.pvalue(
                llvsr2_pvals,
                meta.method = c("roP","maxP", "SR"),
                rth = 2, 
                parametric = FALSE
        )

## Merging multiple meta-analyses results for visualization -----------------------

res <- cbind(llvsr2.mul$ind.p, llvsr2.mul$meta.analysis$pval)
tmp <- cbind(meta.res$meta.analysis$pval, meta.res$meta.analysis$FDR)
tmp <- tmp[row.names(res),]
res <- cbind(res, tmp[,1])
rm(tmp)

colnames(res)[6] <- "REM" 

pdf("fig/DEGs_Meta.pdf", width = 6, height = 6)
draw.DEnumber_all(result = res, maxcut = 0.1, FDR = TRUE, main = "Number of DEG by individual and meta-analyses methods \n(LL vs. BT)")
dev.off()
rm(res)

## Grabbing FDR according to all methods and saving ------------------------

## Grabbing study specific FDR 

i <- Reduce(intersect, list(row.names(p_gse16844), row.names(p_gse74481)))

# removing ribosomal genes
rib <- read.table("../ribosomal_proteins.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
i <- setdiff(i, rib$GeneID)

tmp <- cbind(p_gse16844[i,]$adj.P.Val, p_gse74481[i,]$adj.P.Val)
row.names(tmp) <- i
colnames(tmp) <- c("GSE16844", "GSE74481")

res <- cbind(tmp[i, ], llvsr2.mul$meta.analysis$FDR[i, ])

res <- cbind(res[i, ], meta.res$meta.analysis$FDR[i,])
colnames(res)[6] <- "REM" 
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

WriteXLS(res, ExcelFileName = "results/llvsr2_allMethods.xls", row.names = TRUE)

## NOT RUN
# 
# ## Venn Diagram by meta-analysis methods ------------------------------
# 
# # Sources
# 
 
# venn.diagram(
#         x = list(
#                 SDEF = sdef,
#                 REM = rem,
#                 roP = rop,
#                 maxP = maxp,
#                 SR = sr
#         ),
#         imagetype = "tiff",
#         resolution = 500,
#         filename = "fig/Venn_llvsr2.tiff",
#         col = "black",
#         width = 3000,
#         height = 3000,
#         fill = c("firebrick1", "dodgerblue", "yellow", "limegreen", "grey"),
#         alpha = rep(.3, 5),
#         ext.line.lwd = 2,
#         scaled = FALSE,
#         ext.text = TRUE,
#         lwd = rep(2, 5),
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

sdef <- read_xls("../../sdef/LLvsR2/Results_sdef_LLvsR2.xls")
sdef <- as.character(sdef$ENTREZID)
 
rem <- row.names(res[res$REM < 0.1, ])
rop <- row.names(res[res$roP < 0.1,])
maxp <- row.names(res[res$maxP < 0.1,])
sr <- row.names(res[res$SR < 0.1,])

intersection <- data.frame(
        ENTREZID = Reduce(intersect,
                        list(
                                sdef,
                                rem,
                                rop,
                                maxp, sr
                        )))

intersection$SYMBOL <- select(org.Hs.eg.db, keys = as.character(intersection$ENTREZID), columns = "SYMBOL", "ENTREZID")$SYMBOL

WriteXLS(intersection, ExcelFileName = "results/Intersect_fdr0.1.All.xls")
 
# ## saving the intersection between all methods except SR --------------
 
intersection <- data.frame(
        ENTREZID = Reduce(intersect,
                          list(
                                  sdef,
                                  rem,
                                  rop,
                                  maxp
                          )))

intersection$SYMBOL <- select(org.Hs.eg.db, keys = as.character(intersection$ENTREZID), columns = "SYMBOL", "ENTREZID")$SYMBOL
 
WriteXLS(intersection, ExcelFileName = "results/Intersect_fdr0.1_Except_SR.xls")
 
# ### END 

rm(list = ls())
q(save = "no", status = 0L)
