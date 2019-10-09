## Thyago Leal Calvo - thyagoleal@yahoo.com
## May 2018.

## Custom version of MetaDE::draw.DEnumber() function for some situations
## v.0.1

# Z-score function based on t distribution

ci <- function(var,
               df,
               es = NULL,
               conf.level = 0.95) {
        
       stopifnot(!is.null(df))
                if (is.null(es)) {
                        t <- qt(1 - (1 - conf.level) / 2, df)
                        ci <- t * sqrt(var)
                        return(ci)
                        } else {
                                t <- qt(1 - (1 - conf.level) / 2, df)
                                ci.l <- es - t * sqrt(var)
                                ci.r <- es + t * sqrt(var)
                                return(data.frame(es, ci.l, ci.r))
                        }
                }


draw.DEnumber2 <-
        function (result,
                  maxcut,
                  mlty = NULL,
                  mcol = NULL,
                  mlwd = NULL,
                  mpch = NULL,
                  FDR = TRUE)
        {
                if (class(result) == "MetaDE.pvalue") {
                        pm <- cbind(result$ind.p, result$meta.analysis$pval)
                }
                else if (class(result) == "MetaDE.ES") {
                        pm <- cbind(result$pval)
                        colnames(pm) <- attr(result, "meta.method")
                }
                else if (class(result) == "MetaDE.minMCC") {
                        pm <- cbind(result$meta.analysis$pval)
                        colnames(pm) <-
                                attr(result$meta.analysis, "meta.method")
                }
                else {
                        pm <- result
                }
                method <- colnames(pm)
                if (FDR)
                        pm <-
                        cbind(apply(pm, 2, function(x)
                                p.adjust(x, method = "BH")))
                maxp <- max(pm, na.rm = T)
                if (maxcut > maxp) {
                        cat(
                                "Your maximum cut point exceeds the maximum of observed p-value/FDR\n",
                                "we will use",
                                maxp,
                                "as the maximum cut point\n"
                        )
                        maxcut <- maxp
                }
                ns <- ncol(pm)
                ymax <-
                        max(apply(pm, 2, function(x)
                                sum(x <= maxcut, na.rm = T)))
                if (is.null(mlty))
                        mlty = 1:ns
                if (is.null(mcol))
                        mcol = 1:ns
                if (is.null(mlwd))
                        mlwd = rep(2, ns)
                if (is.null(mpch))
                        mpch = 1:ns
                xlab0 <- ifelse(FDR, "FDR cut-off", "p-value cut-off")
                get.c <- function(cut, pm) {
                        s <- apply(pm, 2, function(x, y)
                                sum(x <= y, na.rm = T),
                                y = cut)
                        return(sum(dist(cbind(
                                cut, s
                        ))))
                }
                mycut <- as.matrix(seq(0, maxcut, length = 20))
                dis <- apply(mycut, 1, get.c, pm = pm)
                minx.pos <- mycut[which.max(dis)]
                plot(
                        c(0, maxcut),
                        c(1, ymax),
                        type = "n",
                        xlab = xlab0,
                        ylab = "Significant genes"
                )
                for (i in 1:ns) {
                        y.pos <- sum(pm[, i] <= minx.pos, na.rm = T)
                        if (y.pos == 0) {
                                x.pos <- minx.pos
                        }
                        else {
                                x.pos <- sort(pm[, i])[y.pos]
                        }
                        points(
                                x.pos,
                                y.pos,
                                pch = mpch[i],
                                col = mcol[i],
                                lwd = 3
                        )
                        lines(
                                sort(pm[, i]),
                                rank(sort(pm[, i]), ties.method = "max"),
                                lty = mlty[i],
                                col = mcol[i],
                                lwd = mlwd[i]
                        )
                }
                legend(
                        "topleft",
                        method,
                        lty = mlty,
                        lwd = mlwd,
                        col = mcol,
                        bty = "n",
                        pch = mpch
                )
        }

### draw.DEnumber_all ###################
# This mod function will draw the DE plot for any tests considering they're formated adequately

draw.DEnumber_all <- function(result,
                               maxcut,
                               mlty = NULL,
                               mcol = NULL,
                               mlwd = NULL,
                               mpch = NULL,
                               FDR = TRUE, 
                               main = NULL)
{
        
       pm <- result
        method <- colnames(pm)
        if (FDR)
                pm <-
                cbind(apply(pm, 2, function(x)
                        p.adjust(x, method = "BH")))
        maxp <- max(pm, na.rm = T)
        if (maxcut > maxp) {
                cat(
                        "Your maximum cut point exceeds the maximum of observed p-value/FDR\n",
                        "we will use",
                        maxp,
                        "as the maximum cut point\n"
                )
                maxcut <- maxp
        }
        ns <- ncol(pm)
        ymax <-
                max(apply(pm, 2, function(x)
                        sum(x <= maxcut, na.rm = T)))
        if (is.null(mlty))
                mlty = 1:ns
        if (is.null(mcol))
                mcol = 1:ns
        if (is.null(mlwd))
                mlwd = rep(2, ns)
        if (is.null(mpch))
                mpch = 1:ns
        xlab0 <- ifelse(FDR, "FDR cut-off", "p-value cut-off")
        get.c <- function(cut, pm) {
                s <- apply(pm, 2, function(x, y)
                        sum(x <= y, na.rm = T),
                        y = cut)
                return(sum(dist(cbind(
                        cut, s
                ))))
        }
        mycut <- as.matrix(seq(0, maxcut, length = 20))
        dis <- apply(mycut, 1, get.c, pm = pm)
        minx.pos <- mycut[which.max(dis)]
        plot(
                c(0, maxcut),
                c(1, ymax),
                type = "n",
                xlab = xlab0,
                ylab = "Significant genes", main = main)
        for (i in 1:ns) {
                y.pos <- sum(pm[, i] <= minx.pos, na.rm = T)
                if (y.pos == 0) {
                        x.pos <- minx.pos
                }
                else {
                        x.pos <- sort(pm[, i])[y.pos]
                }
                points(
                        x.pos,
                        y.pos,
                        pch = mpch[i],
                        col = mcol[i],
                        lwd = 3
                )
                lines(
                        sort(pm[, i]),
                        rank(sort(pm[, i]), ties.method = "max"),
                        lty = mlty[i],
                        col = mcol[i],
                        lwd = mlwd[i]
                )
        }
        legend(
                "topleft",
                method,
                lty = mlty,
                lwd = mlwd,
                col = mcol,
                bty = "n",
                pch = mpch, y.intersp = 0.8
        )
}

### Function to draw the heatmap

heatmap <- function(list, genes) {
        
        ## Part 1 - Formatting and extracting creating labels 
        tmp <- list()
        for (j in seq_along(1:length(list))) {
                tmp[[j]] <- unlist(list[[j]][[1]])
        }
        rm(j)
        mat <- do.call(cbind, tmp)
        b <- vector(mode = "character")
        for (j in seq_along(1:length(list))) {
                b <- append(b, unlist(list[[j]][[2]]))
        }
        colnames(mat) <- make.names(b, unique = TRUE)
        rm(j)
        lab <- vector(mode = "character")
        for (j in seq_along(1:length(list))) {
                lab <- append(studies, values = rep(paste(j, sep = ""), length(list[[j]][[2]])))
        }
        rm(j)
        
        ## Part 2 - Standardizing matrix
        
        mat <- data.matrix(mat)
        leg <- data.frame(Group = b, Study = lab, row.names = make.names(b, unique = TRUE))
        pheatmap::pheatmap(mat[genes, ], 
                           scale = "row", 
                           border_color = NA, 
                           cluster_rows = FALSE,
                           cluster_cols = FALSE,
                           show_colnames = FALSE, 
                           clustering_distance_rows = "correlation",
                           legend = FALSE,
                           annotation_col = NULL
        )
        rm(b, mat, leg)
}



