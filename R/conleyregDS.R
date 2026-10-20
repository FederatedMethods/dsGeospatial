#' @title Linear Regression with Conley Spatial Standard Errors
#'
#' @description
#' Performs linear regression with spatial heteroskedasticity and
#' autocorrelation consistent (Conley) standard errors on server-side
#' spatial data using \code{conleyreg::conleyreg()}.
#'
#' The function accounts for spatial correlation between observations
#' when estimating coefficient standard errors, while retaining the
#' ordinary least squares (OLS) coefficient estimates.
#'
#' DataSHIELD disclosure checks are applied before fitting the model.
#'
#' @details
#' The regression model is expressed as:
#'
#' \deqn{y_i = \beta_0 + \sum_{k=1}^{p}\beta_k x_{ik} + \epsilon_i}
#'
#' where \eqn{y_i} is the dependent variable, \eqn{x_{ik}} represents
#' the independent variables, and \eqn{\epsilon_i} is the error term.
#'
#' Unlike ordinary linear regression, Conley regression adjusts the
#' estimated standard errors to account for spatial dependence between
#' observations within a specified distance cutoff.
#'
#' The distance cutoff determines the maximum distance over which
#' spatial correlation is considered when estimating the standard errors.
#' It does not restrict the observations used to estimate the
#' regression coefficients.
#'
#' The function applies the following DataSHIELD disclosure checks:
#' \itemize{
#'   \item The number of complete observations must be at least
#'   \code{nfilter.tab}.
#'   \item The number of regression parameters must not exceed
#'   \code{nfilter.glm} multiplied by the sample size.
#' }
#'
#' If the checks are satisfied, the model is fitted using
#' \code{conleyreg::conleyreg()} and the coefficient table is returned.
#'
#' These checks provide baseline disclosure protection but do not
#' independently guarantee that all model outputs are disclosure-safe.
#'
#' @param formula An R formula specifying the dependent variable
#'   and one or more independent variables for the regression model.
#'
#' @param data A server-side spatial data object containing the
#'   variables specified in \code{formula} and spatial coordinates
#'   required for calculating spatial distances.
#'
#' @param dist_cutoff A positive numeric value specifying the maximum
#'   distance over which spatial correlation is accounted for.
#'   The units depend on the spatial distance configuration used
#'   by \code{conleyreg::conleyreg()}.
#'
#' @return
#' A named list containing the following element when the model
#' passes the disclosure checks:
#'
#' \describe{
#'   \item{\code{coef}}{
#'     Matrix containing the regression coefficient estimates,
#'     Conley-adjusted standard errors, t-statistics, and p-values.
#'   }
#' }
#'
#' If the disclosure checks fail, a named list is returned containing:
#'
#' \describe{
#'   \item{\code{Table_valid}}{
#'     Logical value \code{FALSE}, indicating that the model
#'     did not satisfy the disclosure requirements.
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
#'     Message explaining why the disclosure check failed.
#'   }
#' }
#'
#' @seealso
#' \code{\link[conleyreg]{conleyreg}},
#' \code{\link[dsBase]{listDisclosureSettingsDS}}
#'
#' @export
conleyregDS <- function (formula, data, dist_cutoff) {
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
    model <- conleyreg::conleyreg(
      formula, 
      data = data, 
      dist_cutoff = dist_cutoff)
    
    res <- list(
      coef = as.matrix(model)
    )
    return(res)
    
  }
  
}
