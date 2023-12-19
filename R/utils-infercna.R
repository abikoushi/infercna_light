rowcenter = function(m, by = 'mean') {
    m = as.matrix(m)
    if (by == 'mean')  by = rowMeans(m, na.rm = T)
    else if (by == 'median') by = matrixStats::rowMedians(m, na.rm = T)
    else stopifnot(is.numeric(by) & length(by) == nrow(m))
    t(scale(t(m), center = by, scale = F))
}

colcenter = function(m, by = 'mean') {
    m = as.matrix(m)
    if (by == 'mean')  by = colMeans(m, na.rm = T)
    else if (by == 'median') by = matrixStats::colMedians(m, na.rm = T)
    else stopifnot(is.numeric(by) & length(by) == ncol(m))
    scale(m, center = by, scale = F)
}

unlogtpm = function(m, bulk = F) {
    # wrapper around scalop::tpm since scalop::tpm is confusing..
    # in that it does not generate tpm from counts, but rather removes log, scaling and pseudocount
    # if (has_dim(m)) 
    m = as.matrix(m)
    if (bulk) x = 1
    else x = 10
    #(2^m) * x - 1 
    (2^(m) - 1) * x
}

logtpm = function(m, bulk = F) {
    #if (has_dim(m))
    m = as.matrix(m)
    if (bulk) x = 1
    else x = 10
    log2((m/x) + 1)
}

clip <- function(m, range = c(-3, 3)) {
    m = as.matrix(m)
    m[m < range[[1]]] <- range[[1]]
    m[m > range[[2]]] <- range[[2]]
    m
}

#' @title Rolling Means
#' @description Apply a rolling window mean to a matrix or vector.
#' @param m a numeric vector or matrix. If the latter, each column will be processed separately.
#' @param k width of rolling window. Default: 100
#' @param endrule character string indicating how the values at the beginning and the end of the data should be treated. One of "mean", "trim", "keep", "constant". See caTools::runmean for more details. Default: 'mean'
#' @param align specifies whether result should be centered (default), left-aligned or right-aligned. See caTools::runmean for more details. Default: 'center'
#' @param verbose print progress messages. Default: TRUE
#' @return a numeric vector or matrix of the same size as <m>. Only in case of endrule=trim, the output vectors will be shorter and output matrices will have fewer rows. 
#' @seealso 
#'  \code{\link[caTools]{runmean}}
#' @rdname runMean
#' @export 
#' @importFrom caTools runmean
runMean <- function(m,
                    k = 100,
                    endrule = 'mean',
                    align = 'center',
                    verbose = T) {

    if (!is.null(dim(m))) {
        m = as.matrix(m)
    }

    if (is.null(dim(m))) return(FALSE)
    if (nrow(m) == 0) return(FALSE)

    if (nrow(m) < k) {
        k = nrow(m)
        if (verbose) message('Adjusting window to the max. number of genes in chromosome (', k, ')')
    }

    mout = caTools::runmean(m,
                            k = k,
                            endrule = endrule,
                            align = align)

    if (!is.null(dim(m))) {
        colnames(mout) = colnames(m)
        rownames(mout) = rownames(m)
        mout = as.data.frame(mout)
    } else {
        names(mout) = names(m)
    }

    mout
}

