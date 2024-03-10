check_positive <- function(x, len) {
    length(x)==len && is.numeric(x) && all(x>=0)
}

kernel_type <- function(kern) {
    kernels <- c("optimal", "uniform", "triangular", "epanechnikov")
    if (is.function(kern)) {
        "function"
    } else if (is.character(kern) && kern %in% kernels) {
        kern
    } else {
        stop(paste(c("'kern' must be a function or one of:", kernels),
                   collapse=" "))
    }
}

process_options <- function(M, se.method, method, d, kern) {
    if (kernel_type(kern)=="optimal" && method!="SRD")
        stop("Optimal kernel requires sharp RD design.")
    if (!missing(M)) {
        m_len <- if (method=="FRD") 2 else 1
        if (!check_positive(M, m_len))
            stop(paste0("M must be a non-negative numeric vector of length",
                        m_len, "."))
    }
    if (!(se.method %in% c("nn", "EHW", "supplied.var"))) {
        stop("Unsupported se.method")
    }
    if (se.method=="nn" && !is.null(d$clusterid))
        stop(paste0("'se.method=\"nn\"' not allowed with clustered standard",
                    "  errors.\nUse 'se.method=\"EHW\"'"))

    if (ncol(d$covs)>0 && method=="IP")
        stop("Covariates not allowed whem method is 'IP'.")
}
