#' @section Note:
#' \code{subset} is evaluated in the same way as variables in \code{formula},
#' that is first in \code{data} and then in the environment of \code{formula}.

#' @param formula object of class \code{"formula"} (or one that can be coerced
#'     to that class) of the form \code{outcome ~ running_variable}
#' @param data optional data frame, list or environment (or object coercible by
#'     \code{as.data.frame} to a data frame) containing the outcome and running
#'     variables in the model. If not found in \code{data}, the variables are
#'     taken from \code{environment(formula)}, typically the environment from
#'     which the function is called.
#' @param subset optional vector specifying a subset of observations to be used
#'     in the fitting process.
#' @param cutoff specifies the RD cutoff in the running variable.
#' @param weights Optional vector of weights to weight the observations
#'     (useful for aggregated data). Disregarded if optimal kernel is used.
#' @param na.action function which indicates what should happen when the data
#'     contain \code{NA}s. The default is set by the \code{na.action} setting of
#'     \code{options} (usually \code{na.omit}).
