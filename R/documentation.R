#' Head Start data from Ludwig and Miller (2007)
#'
#' Subset of Ludwig-Miller data. Counties with missing poverty rate, or with
#' both outcomes missing (\code{hs} and \code{mortality}) were removed. In the
#' original dataset, Yellowstone County, MT (\code{oldcode = 27056}) was entered
#' twice, here the duplicate is removed. Yellowstone National Park, MT
#' (\code{oldcode = 27057}) is also removed due to it being an outlier for both
#' outcomes. Counties with \code{oldcode} equal to (3014, 32032, 47010, 47040,
#' 47074, 47074, 47078, 47079, 47096) matched more than one FIPS entry, so the
#' county labels may not be correct. Mortality data is missing for Alaska.
#' @format A data frame with 3,127 rows and 9 variables:
#'
#' \describe{
#' \item{statefp}{State FIPS code}
#' \item{countyfp}{County FIPS code}
#' \item{oldcode}{ID in Ludwig-Miller dataset}
#' \item{povrate60}{Poverty rate in 1960 relative to 300th poorest county (which
#'      had poverty rate 59.1984)}
#'
#' \item{mortHS}{Average Mortality rate per 100,000 for children aged 5-9 over
#' 1973--83 due to causes addressed as part of Head Start's health services.}
#'
#' \item{mortInj}{Average Mortality rate per 100,000 for children aged 5-9 over
#' 1973--83 due to injury.}
#'
#' \item{highSchool}{High school completion rate in 1990 census, ages 18-24}
#' \item{statepc}{State postal code}
#' \item{county}{County name}
#' }
#' @source Douglas Miller's website,
#' \url{http://faculty.econ.ucdavis.edu/faculty/dlmiller/statafiles/}
#' @references{
#'
#' \cite{Ludwig , J., and D. L. Miller (2007): "Does Head Start improve
#'       children's life chances? Evidence from a regression discontinuity
#'       design," Quarterly Journal of Economics, 122 (1), 159-208.}
#'
#' }
"headst"


#' Austrian unemployment duration data from Lalive (2008)
#'
#' Subset of Lalive data for individuals in the regions affected by the REBP
#' program
#'
#' @format A data frame with 29,371 rows and 4 variables:
#' \describe{
#'   \item{age}{Age in years, at monthly accuracy}
#'   \item{period}{Indicator for whether REBP is in place}
#'   \item{female}{Indicator for female}
#'   \item{duration}{unemployment duration in weeks}
#'  }
#' @references{
#'
#' \cite{Lalive, R. (2008): "How Do Extended Benefits Affect Unemployment
#'       Duration? A Regression Discontinuity Approach." Journal of Econometrics,
#'       142 (2): 785-806.}
#'
#' }
#' @source Rafael Lalive's website,
#'     \url{https://sites.google.com/site/rafaellalive/}
"rebp"


#' Lee (2008) US House elections dataset
#'
#' @format A data frame with 6,558 rows and 2 variables:
#' \describe{
#'   \item{voteshare}{Vote share in next election}
#'   \item{margin}{Democratic margin of victory}
#'  }
#' @source Mostly Harmless Econometrics data archive,
#' \url{https://economics.mit.edu/faculty/angrist/data1/mhe}
#' @references{
#'
#' \cite{Lee , D. S. (2008): "Randomized experiments from non-random selection in
#'       U.S. House elections," Journal of Econometrics, 142 (2), 675-697.}
#'
#' }
"lee08"

#' Oreopoulos (2006) UK general household survey dataset
#'
#' @format A data frame with 73,954 rows and 2 variables:
#' \describe{
#'   \item{earnings}{Annual earnings in 1998 (UK pounds)}
#'   \item{yearat14}{Year individual turned 14}
#'  }
#' @references{
#'
#' \cite{Oreopoulos, P. (2006): "Estimating Average and Local Average Treatment
#'       Effects When Compulsory Education Schooling Laws Really Matter",
#'       American Economic Review, 96 (1), 152-175}
#'
#' }
#' @source American Economic Review data archive,
#' \url{https://doi.org/10.1257/000282806776157641}
"cghs"

#' Battistin et al. (2009) retirement consumption puzzle dataset
#'
#' @format A data frame with 30,006 rows and 6 variables:
#' \describe{
#'   \item{survey_year}{Survey year}
#'   \item{elig_year}{Years to/from eligibility (males)}
#'   \item{retired}{Retirement status (males)}
#'   \item{food}{Total household food expenditure}
#'   \item{c}{Total household consumption}
#'   \item{cn}{Total household expenditure on non-durable goods}
#' }
#'
#' @references{
#' \cite{Battistin, Erich, Agar Brugiavini, Enrico Rettore, and
#'     Guglielmo Weber. 2009. “The Retirement Consumption Puzzle: Evidence from
#'     a Regression Discontinuity Approach.” American Economic Review 99 (5):
#'     2209–26.}
#' }
#' @source American Economic Review data archive,
#' \url{https://doi.org/10.1257/aer.99.5.2209}
"rcp"

#' Constants for common kernels.
#'
#' First four moments of uniform, triangular, and Epanechnikov equivalent
#' kernels. Up to numerical integration precision, these moments are matched by
#' \code{KernMoment()}. See vignette \code{lpkernels}
#'
#' @format A data frame with 18 rows and 19 variables:
#' \describe{
#'   \item{kernel}{Kernel type.}
#'   \item{order}{Order of local polynomial.}
#'   \item{boundary}{Boundary regression?}
#'   \item{mu0, mu1, mu2, mu3, mu4}{\eqn{\int_X u^j k(u) d u}, raw moments}
#'   \item{nu0, nu1, nu2, nu3, nu4}{\eqn{\int_X u^j k^2(u) d u}, raw moments of
#'         kernel squared}
#'   \item{pi0, pi1, pi2, pi3, pi4}{\eqn{\int_X |u^j k(u)| d u}, absolute moments}
#'   \item{pMSE}{constant for pointwise MSE optimal bandwidth,
#'        \eqn{((p+1)!^2\nu_0 / (2(p+1)\mu_{p+1}^2))^{1/(2p+3)}}, see page 67 in
#'        Fan and Gijbels (1996)}
#'  }
#' @source Computed analytically using symbolic math software
#' @references{
#'
#' \cite{Fan , J., and I. Gijbels (1996): Local Polynomial Modelling and Its
#' Applications, Monographs on Statistics and Applied Probability. Chapman &
#' Hall/CRC, New York, NY.}
#'
#' }
"kernC"
