#' @title Fit a Spatial Lag of X (SLX) Regression Model
#'
#' @description
#' Fits a Spatial Lag of X (SLX) regression model using
#' \code{spatialreg::lmSLX()}, applying DataSHIELD disclosure checks
#' before fitting the model.
#'
#' @details
#' The SLX model incorporates spatially lagged explanatory variables
#' to examine the association between an outcome and the characteristics
#' of both local and neighbouring locations.
#'
#' The model is defined as:
#' \deqn{y = X\beta + WX\theta + \epsilon}
#'
#' where \eqn{y} is the dependent variable, \eqn{X} is the matrix of
#' explanatory variables, \eqn{W} is the spatial weights matrix,
#' \eqn{\beta} represents the coefficients of local explanatory
#' variables, and \eqn{\theta} represents the coefficients of
#' spatially lagged explanatory variables.
#'
#' Before fitting the model, the function retrieves the DataSHIELD
#' disclosure thresholds \code{nfilter.tab} and \code{nfilter.glm}.
#' It checks:
#' \itemize{
#'   \item Whether the number of complete observations is at least
#'   \code{nfilter.tab}.
#'   \item Whether the number of model parameters satisfies the
#'   \code{nfilter.glm} threshold.
#' }
#'
#' The model is fitted with \code{zero.policy = TRUE}, allowing
#' observations with no neighbours in the spatial weights object.
#'
#' @param formula An R formula specifying the dependent variable
#'   and explanatory variables.
#' @param LISTW A spatial weights object of class \code{listw},
#'   defining the spatial relationships between observations.
#' @param data A server-side data frame containing the variables
#'   specified in \code{formula}.
#'
#' @return
#' If the disclosure checks are satisfied, a list containing:
#' \describe{
#'   \item{\code{coef}}{Estimated regression coefficients, including
#'   coefficients for spatially lagged explanatory variables.}
#' }
#'
#' If the disclosure checks fail, a list containing:
#' \describe{
#'   \item{\code{Table_valid}}{Logical value \code{FALSE}.}
#'   \item{\code{Nvalid}}{Number of complete observations.}
#'   \item{\code{Nmissing}}{Number of incomplete observations.}
#'   \item{\code{Ntotal}}{Total number of observations.}
#'   \item{\code{Warning}}{Message explaining why the disclosure
#'   checks failed.}
#' }
#'
#' @seealso \code{\link[spatialreg]{lmSLX}}
#'
#' @export
## next function 
# lmSLX 
# SLX: Is obesity prevalence associated with fast-food 
# accessibility in neighbouring areas?
# SAR: Is obesity prevalence associated with obesity prevalence in neighbouring 
# areas, after accounting for local predictors?
lmSLXDS <- function (formula,  LISTW, data) {
  #############################################################
  # MODULE 1: CAPTURE THE nfilter SETTINGS
  thr <- dsBase::listDisclosureSettingsDS()
  nfilter.tab <- as.numeric(thr$nfilter.tab)
  nfilter.glm <- as.numeric(thr$nfilter.glm)
  #nfilter.subset <- as.numeric(thr$nfilter.subset)
  #nfilter.string <- as.numeric(thr$nfilter.string)
  #############################################################
  model_data <- stats::model.frame(
    formula,
    data = data,
    na.action = stats::na.pass
  )
  
  complete <- stats::complete.cases(model_data)
  
  x.length <- sum(complete)
  x.na <- sum(!complete)
  x.total <- nrow(model_data)
  
  X <- stats::model.response(
    stats::model.frame(formula, data = data)
  )
  
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
    model <- spatialreg::lmSLX(
      formula,
      data = data,
      listw = LISTW,
      zero.policy = TRUE
    )
    
    res <- list(
      coef = model$coefficients)
    return(res)
    
  }
  
}
