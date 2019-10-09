## Thyago Leal Calvo - thyagoleal@yahoo.com
## February 2019.

# ****************************************************
# Confidence interval function based on Z distribution
# ****************************************************

ci <- function(var, conf.level = 0.95) {
        z <- qnorm(p = 1 - (1 - conf.level) / 2)
        ci <- z * sqrt(var)
        return(ci)
}

# ****************************************************
# Confidence interval function based on students't t dist
# ****************************************************

cist <- function(var, conf.level = 0.95, df = NULL) {
        stopifnot(!is.null(df))
        z <- qt(p = 1 - (1 - conf.level) / 2, df = df)
        ci <- z * sqrt(var)
        return(ci)
}
