#' @title Calculate Spatial Eigenvectors in DataSHIELD
#'
#' @description
#' Computes spatial eigenvectors from geographical coordinates using
#' \code{spmoran::weigen()}. The function extracts coordinates from a
#' server-side spatial dataset, calculates spatial eigenvectors, and
#' returns summary statistics describing the resulting eigenvalue
#' distribution.
#'
#' Disclosure checks are applied before performing the analysis to
#' ensure that the number of valid spatial observations meets the
#' minimum threshold specified by the DataSHIELD server.
#'
#' @details
#' Spatial eigenvectors represent spatial patterns across geographical
#' observations and can be used to account for spatial autocorrelation
#' in regression models.
#'
#' The function first checks the number of observations with non-empty
#' geometries and non-missing geometry-validity results against the
#' server-side \code{nfilter.tab} disclosure threshold.
#'
#' If the disclosure check is satisfied, the function extracts the
#' centroid coordinates of the spatial observations and calculates
#' spatial eigenvectors using \code{spmoran::weigen()}.
#'
#' The spatial neighbourhood structure can be defined using
#' k-nearest neighbours or Delaunay triangulation, depending on the
#' specified \code{type}.
#'
#' The function returns the number of retained eigenvectors, summary
#' statistics of their eigenvalues, and the eigenvalue spectrum.
#' The observation-level eigenvector matrix is not returned.
#'
#' @param X A server-side \code{sf} object containing spatial geometries.
#'   Each row represents a spatial observation.
#'
#' @param type A character string specifying the method used to
#'   construct the spatial neighbourhood structure. Supported options
#'   include \code{"knn"} for k-nearest neighbours and \code{"tri"}
#'   for Delaunay triangulation.
#'
#' @param k An integer specifying the number of nearest neighbours
#'   when \code{type = "knn"}.
#'
#' @param threshold A numeric value specifying the relative eigenvalue
#'   threshold used to select spatial eigenvectors.
#'
#' @param enum An integer specifying the maximum number of spatial
#'   eigenvectors to extract.
#'
#' @return
#' A list containing the results of the spatial eigenvector analysis.
#'
#' If the disclosure check is satisfied, the list contains:
#' \describe{
#'   \item{\code{n_eigenvectors}}{
#'     Number of retained spatial eigenvectors.
#'   }
#'   \item{\code{eigenvalue_summary}}{
#'     Named numeric vector containing the first quartile
#'     (\code{Q1}), median (\code{Median}), mean (\code{Mean}),
#'     and third quartile (\code{Q3}) of the retained eigenvalues.
#'   }
#'   \item{\code{spectrum}}{
#'     Numeric vector containing the retained eigenvalues.
#'   }
#' }
#'
#' If the disclosure check fails, the function returns:
#' \describe{
#'   \item{\code{Table_valid}}{
#'     Logical value \code{FALSE}, indicating that the analysis
#'     did not satisfy the disclosure threshold.
#'   }
#'   \item{\code{Nvalid}}{
#'     Number of observations with non-empty geometries and
#'     non-missing geometry-validity results.
#'   }
#'   \item{\code{Nmissing}}{
#'     Number of observations failing the geometry checks.
#'   }
#'   \item{\code{Ntotal}}{
#'     Total number of spatial observations.
#'   }
#'   \item{\code{Warning}}{
#'     Character string describing the disclosure-check failure.
#'   }
#' }
#'
#' @note
#' This function is intended for execution on a DataSHIELD server.
#' The spatial data remain on the server, and only the specified
#' analytical summaries are returned to the client.
#'
#' Passing the minimum observation-count check does not, by itself,
#' guarantee that all returned statistics are disclosure-safe.
#' In particular, releasing the full eigenvalue spectrum or exact
#' counts may require additional disclosure assessment.
#'
#' @seealso
#' \code{\link[spmoran]{weigen}},
#' \code{\link[sf]{st_centroid}},
#' \code{\link[dsBase]{listDisclosureSettingsDS}}
#'
#' @export
weigenDS <- function (X, weights, k, threshold, enum) {
  #############################################################
  # MODULE 1: CAPTURE THE nfilter SETTINGS
  thr <- dsBase::listDisclosureSettingsDS()
  nfilter.tab <- as.numeric(thr$nfilter.tab)
  #nfilter.glm <- as.numeric(thr$nfilter.glm)
  #nfilter.subset <- as.numeric(thr$nfilter.subset)
  #nfilter.string <- as.numeric(thr$nfilter.string)
  #############################################################
  
  
  
  x.total <- nrow(X)
  
  valid <- !sf::st_is_empty(X) &
    !is.na(sf::st_is_valid(X))
  
  x.length <- sum(valid)
  x.na <- sum(!valid)
  
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
    coords <- sf::st_coordinates(sf::st_centroid(X))
    
    weig <- spmoran::weigen(x = coords,
                            type = weights,
                            k = k, 
                            threshold = threshold, 
                            enum = enum)
    res <- list(
      n_eigenvectors = ncol(weig$sf),
      eigenvalue_summary = c(
        Q1 = unname(quantile(weig$ev, 0.25)),
        Median = median(weig$ev),
        Mean = mean(weig$ev),
        Q3 = unname(quantile(weig$ev, 0.75))
      ),
      spectrum = weig$ev
    )
    
    return(res)
    
  }
  
}