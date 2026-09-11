#' Calculate Global Moran's I
#'
#' Computes Global Moran's I for a numeric variable using a supplied spatial
#' weights object. Both an analytical Moran's I test and a Monte Carlo
#' permutation test are performed.
#'
#' The function is intended for use on the DataSHIELD server side. Before
#' calculating Moran's I, the number of finite observations is checked against
#' the server-side `nfilter.tab` disclosure threshold. If the number of valid
#' spatial units is below this threshold, the analysis is not performed.
#'
#' @param X A numeric vector containing one observed value for each spatial
#'   unit. The order of values in `X` must correspond to the order of spatial
#'   units represented in `LISTW`.
#'
#' @param LISTW A spatial weights object of class `listw`, typically created
#'   using `spdep::nb2listw()`. The number of spatial units represented by
#'   `LISTW` must equal `length(X)`.
#'
#' @param nsim Integer specifying the number of Monte Carlo permutations used
#'   by `spdep::moran.mc()`. Default is 199.
#'
#' @details
#' Global Moran's I measures spatial autocorrelation in `X` according to the
#' neighbour relationships and weights defined by `LISTW`.
#'
#' Positive values indicate that spatial units with similar values tend to be
#' located near one another. Negative values indicate that neighbouring units
#' tend to have dissimilar values. Values close to the expectation under
#' spatial randomness indicate little evidence of global spatial
#' autocorrelation.
#'
#' The function performs:
#'
#' \itemize{
#'   \item an analytical Moran's I test using `spdep::moran.test()`;
#'   \item a Monte Carlo permutation test using `spdep::moran.mc()`.
#' }
#'
#' Spatial units with no neighbours are permitted by setting
#' `zero.policy = TRUE`.
#'
#' Only aggregate test results are returned. Individual spatial-unit values,
#' neighbour identities, spatial lags, and cell-specific contributions to
#' Moran's I are not returned.
#'
#' @return
#' A list.
#'
#' If the disclosure threshold is not satisfied, the list contains:
#'
#' \describe{
#'   \item{Table_valid}{Logical indicating whether the analysis was permitted.}
#'   \item{Nvalid}{Number of finite observations.}
#'   \item{Nmissing}{Number of missing or non-finite observations.}
#'   \item{Ntotal}{Total number of spatial units.}
#'   \item{Warning}{Disclosure-control warning message.}
#' }
#'
#' If the analysis is permitted, the list contains:
#'
#' \describe{
#'   \item{Table_valid}{TRUE.}
#'   \item{Nvalid}{Number of finite observations.}
#'   \item{Nmissing}{Number of missing or non-finite observations.}
#'   \item{Ntotal}{Total number of spatial units.}
#'   \item{moran_I}{Observed Global Moran's I statistic.}
#'   \item{expectation}{Expected Moran's I under the null hypothesis.}
#'   \item{variance}{Variance of Moran's I under the analytical test.}
#'   \item{z}{Standardised Moran's I test statistic.}
#'   \item{p_analytic}{Analytical p-value from `spdep::moran.test()`.}
#'   \item{p_mc}{Monte Carlo permutation p-value from `spdep::moran.mc()`.}
#' }
#'
#' @seealso
#' `spdep::moran.test()`, `spdep::moran.mc()`
#' @importFrom spdep moran.test moran.mc
#' @keywords internal
#' @export
MoransI <- function (X, LISTW, nsim) {
  #############################################################
  # MODULE 1: CAPTURE THE nfilter SETTINGS
  #thr <- dsBase::listDisclosureSettingsDS()
  nfilter.tab <- as.numeric(thr$nfilter.tab)
  #nfilter.glm <- as.numeric(thr$nfilter.glm)
  #nfilter.subset <- as.numeric(thr$nfilter.subset)
  #nfilter.string <- as.numeric(thr$nfilter.string)
  #############################################################
  
  results <- list()
  
  for (v in 1:length(X)) {
    res.table <- data.frame()
    for (w in names(listw)) {
      x <- X[[v]]
      x.length <- sum(is.finite(x))
      x.na <- sum(!is.finite(x))
      x.total <- length(x)

      if (x.length < nfilter.tab) {
        table.valid <- FALSE
        cell.count.warning <- paste0("Number of cells is less than ", nfilter.tab, ".")
        result <- list(table.valid,
                       X.length,
                       X.na,
                       X.total,
                       cell.count.warning)
        names(result) <- list("Table_valid",
                              "Nvalid",
                              "Nmissing",
                              "Ntotal",
                              "Warning")
        return(result)
      } else{
        mt <- spdep::moran.test(x,
                                listw[[w]],
                                zero.policy = TRUE,
                                na.action = na.omit)
        mc <- spdep::moran.mc(
          x,
          listw[[w]],
          nsim = nsim,
          zero.policy = TRUE,
          na.action = na.omit
        )
    
        res <- data.table(
          variable <- v,
          weights <- w,
          moran_I = unname(mt$estimate[1]),
          z = unname(mt$statistic),
          p_mc = mc$p.value
        )
        res.table <- rbind(res.table, res)
      }
    }
    results[[v]] <- res.table
  }
  return(results)
}
