
#' @title Bayesian Gaussian Regression Using INLA in DataSHIELD
#'
#' @description
#' Fits a Bayesian Gaussian regression model using the
#' \code{INLA::inla()} function on data held within a DataSHIELD
#' server. The function applies server-side disclosure checks
#' before model fitting and returns posterior summaries of the
#' model parameters and model fit statistics.
#'
#' @details
#' Integrated Nested Laplace Approximation (INLA) is a computational
#' approach for approximate Bayesian inference in latent Gaussian
#' models. It provides posterior distributions of model parameters
#' without requiring conventional Markov Chain Monte Carlo (MCMC)
#' sampling.
#'
#' The function first evaluates the specified model formula and
#' determines the number of complete and incomplete observations.
#' It then constructs the model matrix to determine the number of
#' regression parameters.
#'
#' Two DataSHIELD disclosure checks are applied:
#' \itemize{
#'   \item The number of complete observations must be greater than
#'   or equal to the server-side \code{nfilter.tab} threshold.
#'   \item The number of regression parameters must not exceed
#'   \code{nfilter.glm} multiplied by the number of observations
#'   in the model matrix.
#' }
#'
#' If either disclosure check fails, the model is not fitted and
#' a warning is returned.
#'
#' If the checks are satisfied, the function fits a Bayesian
#' Gaussian regression model using \code{INLA::inla()} with
#' \code{family = "gaussian"}.
#'
#' The Deviance Information Criterion (DIC) and Widely Applicable
#' Information Criterion (WAIC) are calculated to support model
#' comparison.
#'
#' The function returns posterior summaries for fixed effects
#' and hyperparameters, together with DIC and WAIC values.
#'
#' Spatial or other random effects may be specified through
#' INLA model formulas, but the current disclosure checks only
#' assess the conventional model matrix and do not account for
#' the number of latent random effects.
#'
#' @param formula An R formula specifying the Bayesian regression
#'   model to be fitted. The left-hand side represents the response
#'   variable, and the right-hand side specifies the explanatory
#'   variables. The response is assumed to follow a Gaussian
#'   distribution.
#'
#' @param data A server-side data frame or compatible data object
#'   containing the variables referenced in \code{formula}.
#'
#' @return
#' A list containing the results of the Bayesian regression model.
#'
#' If the disclosure checks are satisfied, the list contains:
#' \describe{
#'   \item{\code{fixed}}{
#'     Data frame containing posterior summaries of the fixed-effect
#'     coefficients, including posterior means, standard deviations,
#'     quantiles, and other statistics provided by INLA.
#'   }
#'   \item{\code{hyperparameters}}{
#'     Data frame containing posterior summaries of the model
#'     hyperparameters, such as the Gaussian observation precision.
#'   }
#'   \item{\code{dic}}{
#'     Numeric value of the Deviance Information Criterion (DIC).
#'   }
#'   \item{\code{waic}}{
#'     Numeric value of the Widely Applicable Information
#'     Criterion (WAIC).
#'   }
#' }
#'
#' If a disclosure check fails, the function returns:
#' \describe{
#'   \item{\code{Table_valid}}{
#'     Logical value \code{FALSE}, indicating that the disclosure
#'     requirements were not satisfied.
#'   }
#'   \item{\code{Nvalid}}{
#'     Number of complete observations for the specified model.
#'   }
#'   \item{\code{Nmissing}}{
#'     Number of observations containing missing model variables.
#'   }
#'   \item{\code{Ntotal}}{
#'     Total number of observations evaluated.
#'   }
#'   \item{\code{Warning}}{
#'     Character string describing the disclosure-check failure.
#'   }
#' }
#'
#' @note
#' This function is intended for execution on a DataSHIELD server.
#' Individual-level observations are not returned to the client.
#'
#' The function currently supports a Gaussian response distribution
#' only. DIC and WAIC are intended for comparing models fitted to
#' comparable datasets.
#'
#' The disclosure checks are based on the minimum observation count
#' and the ratio of model parameters to observations. These checks
#' alone do not guarantee that the returned Bayesian summaries are
#' disclosure-safe, particularly for models containing spatial
#' random effects or other latent variables.
#'
#' @seealso
#' \code{\link[INLA]{inla}},
#' \code{\link[dsBase]{listDisclosureSettingsDS}}
#'
#' @examples
#' \dontrun{
#' # Example executed within an authorised DataSHIELD server
#'
#' result <- inlaDS(
#'   formula = obesity_qof_2024_25 ~
#'     dist_food_fastfood + pct_deprived_any,
#'   data = cells_object
#' )
#'
#' # Posterior summaries of fixed effects
#' result$fixed
#'
#' # Posterior summaries of hyperparameters
#' result$hyperparameters
#'
#' # Model comparison criteria
#' result$dic
#' result$waic
#' }
#'
#' @export
#' 
inlaDS <- function (formula, data) {
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
    
    model <- INLA::inla(
      formula,
      family = "gaussian",
      data = data,
      control.compute = list(
        dic = TRUE,
        waic = TRUE
        )
      )
    
    res <- list(
      fixed = model$summary.fixed,
      hyperparameters = model$summary.hyperpar,
      dic = model$dic$dic,
      waic = model$waic$waic
    )
    return(res)
    
  }
  
}



