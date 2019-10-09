gozscore <- function(x, y) {
        stopifnot(all(!missing(x), !missing(y)))
        res <- x@result
        exp <-
                structure(y$Median_Log2FC,
                          names = y$SYMBOL,
                          class = "numeric")
        genes <- strsplit(res$geneID, split = "\\/")
        stopifnot(all(unlist(genes) %in% names(exp)))
        zscore <- list()
        for (i in seq_along(genes)) {
                zscore[i] <-
                        (sum(exp[genes[[i]]] > 0) - sum(exp[genes[[i]]] < 0)) / sqrt(length(exp[genes[[i]]]))
        }
        return(unlist(zscore))
}

gozscore2 <- function(x, y) {
        stopifnot(all(!missing(x), !missing(y)))
        res <- x@result
        exp <-
                structure(y$muhat,
                          names = y$Symbol,
                          class = "numeric")
        genes <- strsplit(res$geneID, split = "\\/")
        stopifnot(all(unlist(genes) %in% names(exp)))
        zscore <- list()
        for (i in seq_along(genes)) {
                zscore[i] <-
                        (sum(exp[genes[[i]]] > 0) - sum(exp[genes[[i]]] < 0)) / sqrt(length(exp[genes[[i]]]))
        }
        return(unlist(zscore))
}
