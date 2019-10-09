## Upset plot to show intersection of genes differentially expressed selected by different integration methods 
## Author: Thyago Leal Calvo - thyagoleal@yahoo.com
## Last update: June 14th, 2019
## ********************************

## Settings ---------------------

set.seed(1018)
library(UpSetR)
library(readxl)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(WriteXLS)

## Reading results from meta-analytical methods and SDEF ----------------

## LL vs BT -----------

# SDEF
sdef <- as.character(read_xls("../sdef/LLvsBT/Results_sdef_LLvsBT.xls")$ENTREZID)

# Meta-analytical

meta <- read_xls("../meta-analysis/LLvsBT/results/LLvsBT_allMethods.xls", sheet = "res")
colnames(meta)[1] <- "ENTREZID"
if ("REM" %in% colnames(meta)) {
  # this uses FDR
        rem <- as.character(meta[meta$REM <= 0.1, ]$ENTREZID)
        
} else{
        cat("REM column not found!")
}

## subsetting by FDR

rop <- as.character(meta[meta$roP <= 0.1,]$ENTREZID)
maxp <- as.character(meta[meta$maxP <= 0.1, ]$ENTREZID)
sr <- as.character(meta[meta$SR <= 0.1, ]$ENTREZID)

if (!is.null(rem)) {
        list <-
                list(
                        "REM" = rem,
                        "RoP" = rop,
                        "MaxP" = maxp,
                        "SR" = sr,
                        "SDEF" = sdef
                )
        
} else {
        cat("REM column not present, running with remaining categories")
        list <-
                list(
                        "RoP" = rop,
                        "MaxP" = maxp,
                        "SR" = sr,
                        "SDEF" = sdef
                )
}

df <- fromList(list)
upset(df, order.by = "freq")
dev.copy2pdf(out.type = "pdf",
             file = paste("figs/LLvsBT_upset.pdf"))

## saving intersection lists

# all methods

intersection1 <- data.frame(
                ENTREZID = Reduce(intersect,
                                list(
                                        sdef,
                                        rem,
                                        rop,
                                        maxp,
                                        sr
                                )))

intersection1$SYMBOL <- select(org.Hs.eg.db, keys = as.character(intersection1$ENTREZID), columns = "SYMBOL", "ENTREZID")$SYMBOL

WriteXLS(intersection1, ExcelFileName = "../meta-analysis/LLvsBT/results/Intersect_fdr0.1.xls")

# all except sr

intersection2 <- data.frame(
        ENTREZID = Reduce(intersect,
                          list(
                                  sdef,
                                  rem,
                                  rop,
                                  maxp
                                  
                          )))

intersection2$SYMBOL <- select(org.Hs.eg.db, keys = as.character(intersection2$ENTREZID), columns = "SYMBOL", "ENTREZID")$SYMBOL

WriteXLS(intersection2, ExcelFileName = "../meta-analysis/LLvsBT/results/Intersect_fdr0.1_ExceptSR.xls")

rm(list = ls())

## LL vs R2 ------------

# SDEF
sdef <- as.character(read_xls("../sdef/LLvsR2/Results_sdef_LLvsR2.xls")$ENTREZID)

# Meta-analytical

meta <- read_xls("../meta-analysis/LLvsR2/results/llvsr2_allMethods.xls", sheet = "res")
colnames(meta)[1] <- "ENTREZID"
if ("REM" %in% colnames(meta)) {
        rem <- as.character(meta[meta$REM <= 0.1, ]$ENTREZID)
        
} else{
        cat("REM column not found!")
}
rop <- as.character(meta[meta$roP <= 0.1,]$ENTREZID)
maxp <- as.character(meta[meta$maxP <= 0.1, ]$ENTREZID)
sr <- as.character(meta[meta$SR <= 0.1, ]$ENTREZID)

if (!is.null(rem)) {
        list <-
                list(
                        "REM" = rem,
                        "RoP" = rop,
                        "MaxP" = maxp,
                        "SR" = sr,
                        "SDEF" = sdef
                )
        
} else {
        cat("REM column not present, running with remaining categories")
        list <-
                list(
                        "RoP" = rop,
                        "MaxP" = maxp,
                        "SR" = sr,
                        "SDEF" = sdef
                )
}

df <- fromList(list)
upset(df, order.by = "freq")
dev.copy2pdf(out.type = "pdf",
             file = paste("figs/LLvsR2_upset.pdf"))

## saving intersection lists

# all methods

intersection1 <- data.frame(
        ENTREZID = Reduce(intersect,
                          list(
                                  sdef,
                                  rem,
                                  rop,
                                  maxp,
                                  sr
                          )))

intersection1$SYMBOL <- select(org.Hs.eg.db, keys = as.character(intersection1$ENTREZID), columns = "SYMBOL", "ENTREZID")$SYMBOL

WriteXLS(intersection1, ExcelFileName = "../meta-analysis/LLvsR2/results/Intersect_fdr0.1.xls")


# all methods except sr

intersection2 <- data.frame(
        ENTREZID = Reduce(intersect,
                          list(
                                  sdef,
                                  rem,
                                  rop,
                                  maxp
                          )))

intersection2$SYMBOL <- select(org.Hs.eg.db, keys = as.character(intersection2$ENTREZID), columns = "SYMBOL", "ENTREZID")$SYMBOL

WriteXLS(intersection2, ExcelFileName = "../meta-analysis/LLvsR2/results/Intersect_fdr0.1_Except_SR.xls")

rm(list = ls())

## LL vs Control ------------

# SDEF
sdef <- as.character(read_xls("../sdef/LLvsControl/Results_sdef_LLvsCon.xls")$ENTREZID)

# Meta-analytical

meta <- read_xls("../meta-analysis/LLvsCon/results/llvscon_allMethods.xls", sheet = "res")
colnames(meta)[1] <- "ENTREZID"

rop <- as.character(meta[meta$roP <= 0.1,]$ENTREZID)
maxp <- as.character(meta[meta$maxP <= 0.1, ]$ENTREZID)
sr <- as.character(meta[meta$SR <= 0.1, ]$ENTREZID)

list <- list(
                        "RoP" = rop,
                        "MaxP" = maxp,
                        "SR" = sr,
                        "SDEF" = sdef
                )

df <- fromList(list)
upset(df, order.by = "freq")
dev.copy2pdf(out.type = "pdf",
             file = paste("figs/LLvsControl_upset.pdf"))

# all methods

intersection1 <- data.frame(
        ENTREZID = Reduce(intersect,
                          list(
                                  sdef,
                                  rop,
                                  maxp,
                                  sr
                          )))

intersection1$SYMBOL <- select(org.Hs.eg.db, keys = as.character(intersection1$ENTREZID), columns = "SYMBOL", "ENTREZID")$SYMBOL

WriteXLS(intersection1, ExcelFileName = "../meta-analysis/LLvsCon/results/Intersect_fdr0.1.xls")

# all methods except sr

intersection2 <- data.frame(
        ENTREZID = Reduce(intersect,
                          list(
                                  sdef,
                                  rop,
                                  maxp
                          )))

intersection2$SYMBOL <- select(org.Hs.eg.db, keys = as.character(intersection2$ENTREZID), columns = "SYMBOL", "ENTREZID")$SYMBOL

WriteXLS(intersection2, ExcelFileName = "../meta-analysis/LLvsCon/results/Intersect_fdr0.1_Except_SR.xls")

rm(list = ls())

## In vitro: Inf vs. Control ------------

# SDEF
sdef <- as.character(read_xls("../sdef/invitro/results/Final.DEG.invitro.xls")$ENTREZID)

# Meta-analytical

meta <- read_xls("../meta-analysis/invitro/results/invitro_allMethods.xls", sheet = "res")
colnames(meta)[1] <- "ENTREZID"

rop <- as.character(meta[meta$roP <= 0.1, ]$ENTREZID)
maxp <- as.character(meta[meta$maxP <= 0.1,]$ENTREZID)
sr <- as.character(meta[meta$SR <= 0.1,]$ENTREZID)

list <-
        list(
                "RoP" = rop,
                "MaxP" = maxp,
                "SR" = sr,
                "SDEF" = sdef
        )

df <- fromList(list)
upset(df, order.by = "freq")
dev.copy2pdf(out.type = "pdf",
             file = paste("figs/invitro_upset.pdf"))

# all methods

intersection1 <- data.frame(
        ENTREZID = Reduce(intersect,
                          list(
                                  sdef,
                                  rop,
                                  maxp,
                                  sr
                          )))

intersection1$SYMBOL <- select(org.Hs.eg.db, keys = as.character(intersection1$ENTREZID), columns = "SYMBOL", "ENTREZID")$SYMBOL

WriteXLS(intersection1, ExcelFileName = "../meta-analysis/invitro/results/Intersect_fdr0.1.xls")

# all methods except sr

intersection2 <- data.frame(
        ENTREZID = Reduce(intersect,
                          list(
                                  sdef,
                                  rop,
                                  maxp
                          )))

intersection2$SYMBOL <- select(org.Hs.eg.db, keys = as.character(intersection2$ENTREZID), columns = "SYMBOL", "ENTREZID")$SYMBOL

WriteXLS(intersection2, ExcelFileName = "../meta-analysis/invitro/results/Intersect_fdr0.1_Except_SR.xls")

# exit
rm(list = ls())
q(status = 0, save = "no")
