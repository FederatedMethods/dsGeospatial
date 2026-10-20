#' @title Geographically Weighted Regression (GWR)
#'
#' @description
#' Performs Geographically Weighted Regression (GWR) on server-side spatial
#' data using \code{GWmodel::gwr.basic()}. GWR estimates spatially varying
#' regression coefficients by fitting locally weighted linear regressions
#' at each observation location.
#'
#' The function applies DataSHIELD disclosure checks based on
#' \code{nfilter.tab} and \code{nfilter.glm} before fitting the model.
#' Instead of returning individual local regression coefficients, it returns
#' global regression results, GWR diagnostics, and aggregated summaries
#' of the local coefficients.
#'
#' @details
#' The GWR model is expressed as:
#'
#' \deqn{y_i = \beta_0(u_i,v_i) +
#' \sum_{k=1}^{p}\beta_k(u_i,v_i)x_{ik} + \epsilon_i}
#'
#' where \eqn{(u_i,v_i)} represents the spatial coordinates of observation
#' \eqn{i}, and \eqn{\beta_k(u_i,v_i)} represents the location-specific
#' regression coefficient for predictor \eqn{k}.
#'
#' The function uses the following fixed model settings:
#' \itemize{
#'   \item Bandwidth: 1000 nearest neighbours.
#'   \item Kernel: bisquare.
#'   \item Adaptive bandwidth: TRUE.
#' }
#'
#' If \code{dMat} is NULL, spatial distances are calculated internally
#' by \code{GWmodel::gwr.basic()}.
#'
#' If \code{dMat} is provided, the function checks its dimensions against
#' the number of observations and requires a \code{cutoff_m} attribute.
#' The supplied sparse distance matrix is converted into a dense matrix,
#' with zero entries replaced by the cutoff distance and diagonal entries
#' restored to zero.
#'
#' Disclosure checks include:
#' \itemize{
#'   \item The number of complete observations must be at least
#'   \code{nfilter.tab}.
#'   \item The number of regression parameters must not exceed
#'   \code{nfilter.glm} multiplied by the sample size.
#' }
#'
#' These are global disclosure checks and do not independently establish
#' the disclosure safety of local GWR estimates or their summaries.
#'
#' To reduce disclosure risk, the function returns only the first quartile,
#' median, mean, and third quartile of each local regression coefficient.
#' Minimum and maximum coefficient values are excluded.
#'
#' @param formula An R formula specifying the dependent variable and
#'   independent variables for the GWR model.
#'
#' @param LISTW A spatial weights object. Currently included for
#'   compatibility with other spatial regression functions, but not used
#'   in the GWR calculation.
#'
#' @param data A server-side spatial data object containing the variables
#'   specified in \code{formula} and spatial coordinates required for GWR.
#'
#' @param dMat An optional spatial distance matrix. If \code{NULL},
#'   distances are calculated internally. If provided, the matrix must
#'   have dimensions matching the observations in \code{data} and contain
#'   a valid \code{cutoff_m} attribute.
#'
#' @return
#' A named list containing the following elements when the model passes
#' the disclosure checks:
#'
#' \describe{
#'   \item{\code{global_coefficients}}{
#'     Estimated coefficients from the global linear regression model.
#'   }
#'   \item{\code{global_summary}}{
#'     Global regression coefficient table containing estimates,
#'     standard errors, t-statistics, and p-values.
#'   }
#'   \item{\code{gwr_diagnostics}}{
#'     GWR diagnostic statistics returned by \code{GWmodel}, including
#'     model fit and information criteria.
#'   }
#'   \item{\code{gwr_coefficient_summary}}{
#'     Matrix containing the first quartile (Q1), median, mean, and
#'     third quartile (Q3) of the local GWR coefficient estimates.
#'   }
#' }
#'
#' If the disclosure checks fail, a named list is returned containing:
#' \describe{
#'   \item{\code{Table_valid}}{Logical value \code{FALSE}.}
#'   \item{\code{Nvalid}}{Number of complete observations.}
#'   \item{\code{Nmissing}}{Number of incomplete observations.}
#'   \item{\code{Ntotal}}{Total number of observations.}
#'   \item{\code{Warning}}{Explanation of the failed disclosure check.}
#' }
#'
#' @seealso
#' \code{\link[GWmodel]{gwr.basic}},
#' \code{\link[dsBase]{listDisclosureSettingsDS}}
#'
#' @export
gwrbasicDS <- function (formula,  data, dMat) {
  #############################################################
  # MODULE 1: CAPTURE THE nfilter SETTINGS
  thr <- dsBase::listDisclosureSettingsDS()
  nfilter.tab <- as.numeric(thr$nfilter.tab)
  nfilter.glm <- as.numeric(thr$nfilter.glm)
  #nfilter.subset <- as.numeric(thr$nfilter.subset)
  #nfilter.string <- as.numeric(thr$nfilter.string)
  #############################################################
  args <- list(formula = formula,
               data = data,
               bw = 1000,
               kernel = "bisquare",
               adaptive = TRUE)
  
  
  if(!is.null(dMat)){
    if (nrow(dMat) != nrow(data) ||
        ncol(dMat) != nrow(data)) {
      stop("Dimensions of dMat do not match data.")
    }
    
    cutoff <- attr(dMat, "cutoff_m")
    
    if (is.null(cutoff) || !is.finite(cutoff)) {
      stop("dMat must contain a valid cutoff_m attribute.")
    }
    
    M <- as.matrix(dMat)
    M[M == 0] <- cutoff_m
    diag(M) <- 0
    args$dMat <- M
  }
  
  
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
    
    model <- do.call(GWmodel::gwr.basic, args)
    
    coef_names <- names(coef(model$lm))
    coef_names[coef_names == "(Intercept)"] <- "Intercept"
    gwr_coef <- as.data.frame(model$SDF)[, coef_names, drop = FALSE]
    gwr_summary <- sapply(gwr_coef, function(x) {
      c(
        Q1 = unname(quantile(x, 0.25, na.rm = TRUE)),
        Median = median(x, na.rm = TRUE),
        Mean = mean(x, na.rm = TRUE),
        Q3 = unname(quantile(x, 0.75, na.rm = TRUE))
      )
    })
    
    res <- list(
      global_coefficients = coef(model$lm),
      global_summary = summary(model$lm)$coefficients,
      gwr_diagnostics = model$GW.diagnostic,
      gwr_coefficient_summary = gwr_summary
    )
    
    return(res)
    
  }
  
}
