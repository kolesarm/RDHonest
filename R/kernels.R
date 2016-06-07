#' Equivalent kernel for local linear regression.
#'
#' Calculates equivalent kernel for local polynomial regression.
#' @param kernel kernel type. Can be a function supported on \eqn{[-1, 1]}
#' (boundary kernel) or \eqn{[0,1]} (interior kernel), or else one of "uniform"
#' (\eqn{I(|u|<1)}), "triangular" (\eqn{(1-|u|)I(|u|<1)}), or "epanechnikov"
#' (k(u)=(3/4)(1-u^2)I(|u|<1)).
#' @param boundary Logical scalar, specifying whether we are at a boundary.
#' @param order Order of local polynomial: \code{0} means local constant,
#'     \code{1} local linear, \code{2} local quadratic etc.
#' @return Equivalent kernel function.
#' @examples
#' EqKern(kernel = "uniform", order = 2)
#' @export
EqKern <- function(kernel = "uniform", boundary = TRUE, order = 0) {
    ## support
    su <- function(u) (u <= 1) * (u >= -1 + boundary)
    if(is.function(kernel)) {
        EqKernN(kernel, boundary = boundary, order = order)
    } else if (order > 2) {
        K <- EqKern(kernel = kernel, boundary=boundary, order = 0)
        EqKernN(K, boundary = boundary, order = order)
    } else if (order == 0 && boundary == TRUE) {
        switch(kernel,
               uniform = function(u) su(u),
               triangular = function(u) 2 * (1 - u) * su(u),
               epanechnikov = function(u) (3 / 2) * (1 - u^2) * su(u))
    } else if (order == 1 && boundary == TRUE) {
        switch(kernel,
               uniform = function(u) (4 - 6*u) * su(u),
               triangular = function(u) 6*(1 - 2*u) * (1 - u) * su(u),
               epanechnikov = function(u) 6/19 * (16-30*u) * (1-u^2) * su(u))
    } else if (order == 2 && boundary == TRUE) {
        switch(kernel,
               uniform = function(u) (9 - 36*u + 30*u^2) * su(u),
               triangular = function(u) 12 * (1-5*u+5*u^2) * (1-u) * su(u),
               epanechnikov = function(u) 1/8 * (85 - 400*u + 385*u^2) *
                   (1-u^2) * su(u))
    } else if (order == 0 && boundary == FALSE) {
        switch(kernel,
               uniform = function(u) su(u) / 2,
               triangular = function(u) (1 - abs(u)) * su(u),
               epanechnikov = function(u) (3 / 4) * (1 - u^2) * su(u))
    } else if (order == 1 && boundary == FALSE) {
        switch(kernel,
               uniform = function(u) su(u) / 2,
               triangular = function(u) (1 - abs(u)) * su(u),
               epanechnikov = function(u) 3/4 * (1 - u^2) * su(u))
    } else if (order == 2 && boundary == FALSE) {
        switch(kernel,
               uniform = function(u) (9 - 15 * u^2) / 8 * su(u),
               triangular = function(u) 6/7 * (2-5*u^2) * (1-abs(u)) * su(u),
               epanechnikov = function(u) 15/32 * (3-7*u^2) * (1-u^2) * su(u))
    }
}

#' Moments of a kernel.
#'
#' Computes moments of a kernel over \eqn{X=[0, 1]} (boundary case), or
#' \eqn{X=[-1, 1]} (interior case),
#'
#' @param K kernel function.
#' @param moment order \eqn{j} of moment to compute.
#' @param type Type of moment. "raw" computes \eqn{\int_X u^j k(u)}, "absolute"
#'     computes \eqn{\int_X |u^j k(u)|}, and "raw2" computes
#'     \eqn{\int_X u^j k(u)^2}.
#' @inheritParams EqKern
#' @return Integral value (a scalar).
#' @examples
#' KernMoment(function(u) u<1, moment = 3, boundary = FALSE)
#' KernMoment(EqKern(kernel = "triangular", order = 2),
#'            moment = 3)
#' @export
KernMoment <- function(K, moment = 0, boundary = TRUE, type = "raw") {
    fkt <- switch(type,
                  raw=function(u) u^moment*K(u),
                  absolute=function(u) abs(u^moment*K(u)),
                  raw2=function(u) u^moment*K(u)^2)

    stats::integrate(fkt, lower=-1+boundary, upper=1,
                     rel.tol = .Machine$double.eps^0.75)$value
}

#' M matrix
#' @keywords internal
KernM <- function(K, order = 2, boundary = boundary) {
    M <- outer(0:order, 0:order, "+")

    matrix(sapply(M, function(j)
        KernMoment(K, moment = j, boundary = boundary)),
           nrow = (order + 1))
}

#' Compute Equivalent kernel numerically
#' @inheritParams EqKern
#' @param K original kernel
#' @keywords internal
EqKernN <- function(K, boundary = TRUE, order = 0) {
    s <- drop(solve(KernM(K, order = order, boundary = boundary))[1, ])

    function(u)
        Reduce("+", lapply(seq_along(s), function(j) s[j] * u^(j - 1))) * K(u)
}