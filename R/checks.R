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
    m_len <- if (method=="FRD") 2 else 1
    if (!missing(M)) {
        if (!check_positive(M, m_len))
            stop(paste0("M must be a non-negative numeric vector of length",
                        m_len, "."))
    }

    if (method=="IP") {
        if (!is.null(d$covs))
            stop("Covariates not allowed whem method is 'IP'.")
    } else if (min(sum(d$p), sum(d$m))==0) {
        stop("No observations on one side of the cutoff")
    }

    if (!(se.method %in% c("nn", "EHW", "supplied.var"))) {
        stop("Unsupported se.method")
    }
    if (se.method=="nn" && !is.null(d$clusterid))
        stop(paste0("'se.method=\"nn\"' not allowed with clustered standard",
                    "  errors.\nUse 'se.method=\"EHW\"'"))
}
