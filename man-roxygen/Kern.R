#' @param kern specifies kernel function used in the local regression. It can
#'     either be a string equal to \code{"triangular"} (\eqn{k(u)=(1-|u|)_{+}}),
#'     \code{"epanechnikov"} (\eqn{k(u)=(3/4)(1-u^2)_{+}}), or \code{"uniform"}
#'     (\eqn{k(u)= (|u|<1)/2}), or else a kernel function.
#' @param order Order of local regression 1 for linear, 2 for quadratic.
