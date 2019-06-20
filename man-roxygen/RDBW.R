#' @param h bandwidth, a scalar parameter. For fuzzy or sharp RD, it can be a
#'     named vector of length two with names \code{"p"} and \code{"m"}, in which
#'     case the bandwidth \code{h["m"]} is used for observations below the
#'     cutoff, and the bandwidth \code{h["p"]} is used for observations above
#'     the cutoff. If not supplied, optimal bandwidth is computed according to
#'     criterion given by \code{opt.criterion}.
