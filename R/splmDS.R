
#' @title Spatial Linear Regression with Disclosure Controls
#'
#' @description
#' Fits a spatial linear regression model on server-side spatial data
#' using \code{spmodel::splm()}, while applying DataSHIELD disclosure
#' checks before model estimation.
#'
#' The model accounts for spatial dependence between observations by
#' explicitly modelling the spatial covariance structure of the
#' regression errors.
#'
#' @details
#' The spatial linear model is expressed as:
#'
#' \deqn{y_i = \beta_0 + \sum_{k=1}^{p}\beta_k x_{ik} + \epsilon_i}
#'
#' where \eqn{y_i} is the dependent variable, \eqn{x_{ik}} represents
#' the explanatory variables, and \eqn{\epsilon_i} is a spatially
#' correlated error term.
#'
#' Unlike ordinary least squares regression, the model incorporates
#' spatial correlation into the estimation of regression coefficients
#' and their uncertainty.
#'
#' The function uses the default arguments of \code{spmodel::splm()},
#' including its default spatial covariance structure and estimation
#' settings. For large datasets, \code{spmodel::splm()} may
#' automatically use a computationally efficient local approximation.
#'
#' If the input spatial data contain polygon geometries, the
#' underlying modelling function may convert these to centroids
#' to represent observation locations.
#'
#' Before fitting the model, the function applies the following
#' DataSHIELD disclosure checks:
#'
#' \itemize{
#'   \item The number of complete observations must be at least
#'   \code{nfilter.tab}.
#'   \item The number of regression parameters must not exceed
#'   \code{nfilter.glm} multiplied by the sample size.
#' }
#'
#' If either check fails, model estimation is skipped and a
#' disclosure warning is returned.
#'
#' These checks provide baseline disclosure protection but do not
#' independently guarantee that all spatial model outputs are
#' disclosure-safe.
#'
#' @param formula An R formula specifying the dependent variable
#'   and one or more independent variables for the spatial
#'   linear regression model.
#'
#' @param data A server-side spatial data object, typically an
#'   \code{sf} object, containing the variables specified in
#'   \code{formula} and spatial geometry or coordinates required
#'   for estimating spatial dependence.
#'
#' @param dist_cutoff Currently unused. Included in the function
#'   signature but not passed to \code{spmodel::splm()}.
#'   This argument does not affect model estimation.
#'
#' @return
#' A named list containing the following elements when the model
#' passes the disclosure checks:
#'
#' \describe{
#'   \item{\code{fixed}}{
#'     Estimated fixed-effect regression coefficients.
#'   }
#'   \item{\code{spcov}}{
#'     Estimated spatial covariance parameters describing
#'     the spatial dependence structure.
#'   }
#'   \item{\code{randcov}}{
#'     Estimated random-effect covariance parameters, where
#'     applicable. May be empty or \code{NULL} when no random
#'     effects are specified.
#'   }
#' }
#'
#' If the disclosure checks fail, a named list is returned
#' containing:
#'
#' \describe{
#'   \item{\code{Table_valid}}{
#'     Logical value \code{FALSE}, indicating that the model
#'     failed the disclosure checks.
#'   }
#'   \item{\code{Nvalid}}{
#'     Number of complete observations.
#'   }
#'   \item{\code{Nmissing}}{
#'     Number of observations containing missing values.
#'   }
#'   \item{\code{Ntotal}}{
#'     Total number of observations.
#'   }
#'   \item{\code{Warning}}{
#'     Message explaining why the disclosure checks failed.
#'   }
#' }
#'
#' @seealso
#' \code{\link[spmodel]{splm}},
#' \code{\link[dsBase]{listDisclosureSettingsDS}}
#'
#' @export
#' 

splmDS <- function (formula, data) {
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
    model <- spmodel::splm(
      formula, 
      data = data)
    
    res <- list(
      fixed = model$coefficients$fixed,
      spcov = model$coefficients$spcov,
      randcov = model$coefficients$randcov
      )
    
    return(res)
    
  }
  
}
