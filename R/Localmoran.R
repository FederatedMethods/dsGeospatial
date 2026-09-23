#' Calculate Local Moran's I
#'
#' Computes Local Moran's I for a numeric variable using a supplied spatial
#' weights object. 
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
#'   \item{moran_I}{Observed Global Moran's I statistic.}
#'   \item{z}{Standardised Moran's I test statistic.}
#' }
#'
#' @seealso
#' `spdep::moranlocal()`
#' @importFrom spdep moran.test moran.mc
#' @keywords internal
#' @export
Localmoran <- function (X, LISTW) {
  #############################################################
  # MODULE 1: CAPTURE THE nfilter SETTINGS
  thr <- dsBase::listDisclosureSettingsDS()
  nfilter.tab <- as.numeric(thr$nfilter.tab)
  #nfilter.glm <- as.numeric(thr$nfilter.glm)
  #nfilter.subset <- as.numeric(thr$nfilter.subset)
  #nfilter.string <- as.numeric(thr$nfilter.string)
  #############################################################
  
  
  x.length <- sum(is.finite(X))
  x.na <- sum(!is.finite(X))
  x.total <- length(X)
  
  if (x.length < nfilter.tab) {
    table.valid <- FALSE
    cell.count.warning <- paste0("Number of cells is less than ", nfilter.tab, ".")
    result <- list(table.valid,
                   x.length,
                   x.na,
                   x.total,
                   cell.count.warning)
    names(result) <- list("Table_valid",
                          "Nvalid",
                          "Nmissing",
                          "Ntotal",
                          "Warning")
    return(result)
  } 
  else{
    mt <- spdep::localmoran(X,
                            LISTW,
                            zero.policy = TRUE,
                            na.action = na.omit)
    
    
    res <- list(
      moran_I = as.data.frame(mt),
      quadr = as.data.frame(attr(mt, "quadr"))
    )
    
    return(res)
    
  }
  
}
