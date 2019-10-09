## Code for drawing forest plots for random effects meta-analysis of gene expression
## Thyago Leal Calvo - <thyagoleal@yahoo.com>
## Last update: June 14th, 2019
## *****************************

## Settings ------------------------

set.seed(120)
library(readxl)
library(forestplot)
library(WriteXLS)
source("src/source.R")

## LL vs. BT forest plots --------------------------------

# load data

llvsbt <- as.data.frame(read_xls("../meta-analysis/LLvsBT/results/LLvsBT_ES.xls"))
genes <- as.character(llvsbt$Entrezid[llvsbt$FDR < 0.1 & abs(llvsbt$muhat)>= 1])
genes <- as.character(sort(as.integer(genes), decreasing = FALSE))

# format data

llvsbt <- llvsbt[,c("Entrezid", "GSE2480", "GSE17763", "GSE74481", "GSE443", "muhat", "muvar", "Symbol", "var_GSE2480", "var_GSE17763", "var_GSE74481", "var_GSE443")]

names(llvsbt)[1] <- "ENTREZID"
row.names(llvsbt) <- as.character(llvsbt$ENTREZID)

llvsbt$ci_overall <- ci(llvsbt$muvar)
llvsbt <- cbind(llvsbt, ci(llvsbt[, c("var_GSE2480", "var_GSE17763", "var_GSE74481", "var_GSE443")]))

colnames(llvsbt)[14:17] <- paste("ci", colnames(llvsbt)[2:5], sep = "_")

# draw plot

pdf(
        file = "figs/fp_llvsbt.pdf",
        width = 5,
        height = 4,
        family = "ArialMT",
        title = "Foresplot",
        onefile = TRUE
)

for(i in genes){
                forestplot(
                labeltext = c(colnames(llvsbt)[2:5], "Summary"),
                mean = unlist(llvsbt[i, c(2:6)]),
                lower = unlist(llvsbt[i, c(2:6)] - llvsbt[i, c(14:17, 13)]),
                upper = unlist(llvsbt[i, c(2:6)] + llvsbt[i, c(14:17, 13)]),
                is.summary = c(rep(FALSE, 4), TRUE),
                xlab = expression(paste("Standardized log"[2], "FC", sep = "")),
                lwd.zero = 2,
                col = fpColors(
                        box = "grey30",
                        line = "grey50",
                        summary = "grey10",
                        zero = c("firebrick2")),
                title = llvsbt[i, ]$Symbol,
                #boxsize = 0.25, 
                txt_gp = fpTxtGp(
                        label = gpar(cex = 1.3), 
                        xlab = gpar(cex = 1),
                        title = gpar(cex = 1.5)))
}

dev.off()
rm(i, genes, llvsbt)

## Save LL vs. BT table with CI --------------

llvsbt <- as.data.frame(read_xls("../meta-analysis/LLvsBT/results/LLvsBT_ES.xls"))

names(llvsbt)[1] <- "ENTREZID"

llvsbt$ci_overall <- ci(llvsbt$muvar)

llvsbt$ci95 <- paste("[", round(llvsbt$muhat-llvsbt$ci_overall, digits = 3),", ", round(llvsbt$muhat+llvsbt$ci_overall, digits = 3),"]", sep = "")

llvsbt <- llvsbt[, c("ENTREZID", "Symbol", names(llvsbt)[2:6], "ci95", "FDR", "tau2")]

WriteXLS(llvsbt, ExcelFileName = "../meta-analysis/LLvsBT/results/LLvsBT_ES_95CI.xls", BoldHeaderRow = TRUE)

rm(list = c('llvsbt', 'i', 'genes'))

## LL vs. R2 forest plots --------------------------------

# load data

llvsr2 <- as.data.frame(read_xls("../meta-analysis/LLvsR2/results/llvsr2_ES.xls"))
genes <- as.character(llvsr2$Row.names[llvsr2$FDR < 0.1 & abs(llvsr2$muhat)>= 1])
genes <- as.character(sort(as.integer(genes), decreasing = FALSE))

# format data

llvsr2 <- llvsr2[,c("Row.names", "GSE16844", "GSE74481", "muhat", "muvar", "Symbol", "var_GSE16844", "var_GSE74481")]

names(llvsr2)[1] <- "ENTREZID"
row.names(llvsr2) <- as.character(llvsr2$ENTREZID)

llvsr2$ci_overall <- ci(llvsr2$muvar)
llvsr2 <- cbind(llvsr2, ci(llvsr2[, c("var_GSE16844", "var_GSE74481")]))

colnames(llvsr2)[10:11] <- paste("ci", colnames(llvsr2)[2:3], sep = "_")

# draw plot
pdf(
        file = "figs/fp_llvsr2.pdf",
        width = 4,
        height = 3,
        family = "ArialMT",
        title = "Foresplot",
        onefile = TRUE
)

for(i in genes){
        forestplot(
                labeltext = c(colnames(llvsr2)[2:3], "Summary"),
                mean = unlist(llvsr2[i, c(2:4)]),
                lower = unlist(llvsr2[i, c(2:4)] - llvsr2[i, c(10, 11, 9)]),
                upper = unlist(llvsr2[i, c(2:4)] + llvsr2[i, c(10, 11, 9)]),
                is.summary = c(rep(FALSE, 2), TRUE),
                xlab = expression(paste("Standardized log"[2], "FC", sep = "")),
                lwd.zero = 2,
                col = fpColors(
                        box = "grey30",
                        line = "grey50",
                        summary = "grey10",
                        zero = c("firebrick2")),
                title = llvsr2[i, ]$Symbol,
                #boxsize = 0.25, 
                txt_gp = fpTxtGp(
                        label = gpar(cex = 1.3), 
                        xlab = gpar(cex = 1),
                        title = gpar(cex = 1.5)))
}

dev.off()

## Save LL vs. R2 table with 95ci  -------------- 

llvsr2 <- as.data.frame(read_xls("../meta-analysis/LLvsR2/results/llvsr2_ES.xls"))

names(llvsr2)[1] <- "ENTREZID"

llvsr2$ci_overall <- ci(llvsr2$muvar)

llvsr2$ci95 <- paste("[", round(llvsr2$muhat-llvsr2$ci_overall, digits = 3),", ", round(llvsr2$muhat+llvsr2$ci_overall, digits = 3),"]", sep = "")

llvsr2 <- llvsr2[, c("ENTREZID", "Symbol", names(llvsr2)[2:3], "ci95", "FDR", "tau2")]

WriteXLS(llvsr2, ExcelFileName = "../meta-analysis/LLvsR2/results/llvsr2_ES_95CI.xls", BoldHeaderRow = TRUE)

## end
rm(list = ls())
q(save = "no")
