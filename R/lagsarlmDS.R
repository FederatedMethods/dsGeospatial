#' @title Fit a Spatial Lag Regression Model
#'
#' @description
#' Fits a spatial autoregressive lag (SAR) model on the server side using
#' \code{spatialreg::lagsarlm()}. The function incorporates spatial dependence
#' in the outcome through a spatial weights object while applying DataSHIELD
#' disclosure controls before model fitting.
#'
#' @details
#' The function constructs the model matrix specified by \code{formula} and
#' applies the DataSHIELD \code{nfilter.tab} and \code{nfilter.glm} disclosure
#' settings before fitting the spatial lag model.
#'
#' The \code{nfilter.tab} setting is used to ensure that a sufficient number
#' of valid observations are available. The \code{nfilter.glm} setting prevents
#' fitting an oversaturated model in which the number of model parameters is
#' too large relative to the number of observations.
#'
#' If the disclosure checks are satisfied, a spatial lag model of the general
#' form
#'
#' \deqn{y = \rho Wy + X\beta + \epsilon}
#'
#' is fitted using \code{spatialreg::lagsarlm()}. Here, \eqn{W} represents the
#' spatial weights structure and \eqn{\rho} is the spatial autoregressive
#' parameter describing dependence between the outcome at a location and the
#' outcomes at neighbouring locations.
#'
#' The model is fitted using \code{method = "Matrix"} to allow sparse matrix
#' methods to be used for the spatial weights structure. Locations without
#' neighbours are permitted through \code{zero.policy = TRUE}.
#'
#' Only selected model outputs are returned to the client.
#'
#' @param formula An object of class \code{formula} specifying the dependent
#' variable and explanatory variables to be included in the spatial lag model.
#'
#' @param LISTW A spatial weights object of class \code{listw}, defining the
#' neighbourhood structure and spatial weights between observations.
#'
#' @param data A \code{data.frame} containing the variables specified in
#' \code{formula}. The ordering of observations must correspond to the ordering
#' used in \code{LISTW}.
#'
#' @return
#' If the disclosure checks are satisfied, a list containing:
#' \itemize{
#'   \item{\code{rho}}{The estimated spatial autoregressive parameter.}
#'   \item{\code{coef}}{The estimated regression coefficients.}
#' }
#'
#' If a disclosure check fails, a list is returned containing:
#' \itemize{
#'   \item{\code{Table_valid}}{Logical value indicating that the disclosure
#'   check was not satisfied.}
#'   \item{\code{Nvalid}}{Number of valid values used in the disclosure check.}
#'   \item{\code{Nmissing}}{Number of non-finite values.}
#'   \item{\code{Ntotal}}{Total number of values considered.}
#'   \item{\code{Warning}}{A message describing the disclosure-control failure.}
#' }
#'
#' @seealso
#' \code{\link[spatialreg]{lagsarlm}},
#' \code{\link[dsBase]{listDisclosureSettingsDS}}
#'
#' @export



lagsarlmDS <- function (formula,  LISTW, data) {
  #############################################################
  # MODULE 1: CAPTURE THE nfilter SETTINGS
  thr <- dsBase::listDisclosureSettingsDS()
  nfilter.tab <- as.numeric(thr$nfilter.tab)
  nfilter.glm <- as.numeric(thr$nfilter.glm)
  #nfilter.subset <- as.numeric(thr$nfilter.subset)
  #nfilter.string <- as.numeric(thr$nfilter.string)
  #############################################################
  
  X <- stats::model.response(
    stats::model.frame(formula, data = data)
  )
  
  x.length <- sum(is.finite(X))
  x.na <- sum(!is.finite(X))
  x.total <- length(X)
  X.matrix <- stats::model.matrix(formula, data = data)
  N <- nrow(X.matrix)
  nparams <- ncol(X.matrix)
  
  if (x.length < nfilter.tab | nparams > nfilter.glm * N) {
    table.valid <- FALSE
    if (x.length < nfilter.tab){
    cell.count.warning <- paste0("Number of cells is less than ", nfilter.tab, ".")
    } else{
      cell.count.warning <- paste0("ERROR: Model is oversaturated (too many model parameters relative to sample size)",
                             "leading to a possible risk of disclosure - please simplify model.")
    }
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
    model <- spatialreg::lagsarlm(
      formula,
      data = data,
      listw = LISTW,
      zero.policy = TRUE,
      method = "Matrix"
    )
    
    res <- list(
      rho = unname(model$rho),
      coef = model$coefficients)
    
    
    return(res)
    
  }
  
}
