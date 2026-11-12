#' Prepare Spatial Data for Case and Negative Mapping
#'
#' @description
#' Prepares spatial data for mapping by linking study observations to
#' precomputed display coordinates and filtering cases according to the
#' number of neighbouring cases within a network distance threshold.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Selects precomputed display coordinates for the specified area.
#'   \item Links study observations to distance matrix row identifiers
#'   using the supplied mapping table, if required.
#'   \item Merges study observations with their display coordinates.
#'   \item Counts neighbouring cases using the supplied network
#'   distance matrix.
#'   \item Retains cases with at least \code{k} other cases within
#'   the network distance threshold.
#'   \item Separates retained cases and negative observations and
#'   returns them together with the basemap.
#' }
#'
#' The distance threshold is determined by the construction of \code{D}.
#' For example, when \code{D} represents a 5 km network distance matrix,
#' the neighbouring case counts are calculated within that network.
#'
#' The function assumes that non-zero entries in \code{D} represent
#' connections between locations and that \code{matrix_row_id}
#' corresponds to the row and column indices of \code{D}.
#'
#' @param D A square network distance matrix, potentially sparse,
#'   containing distances between spatial locations. Non-zero entries
#'   indicate locations included in the neighbourhood calculation.
#'
#' @param study A \code{data.frame} or \code{data.table} containing
#'   study observations, including \code{cell_id} or
#'   \code{matrix_row_id}, and a binary outcome variable.
#'
#' @param mapping A \code{data.table} containing \code{cell_id} and
#'   \code{matrix_row_id}, used to associate study observations with
#'   rows of the distance matrix.
#'
#' @param pos A \code{data.table} containing precomputed display
#'   coordinates, with columns \code{matrix_row_id}, \code{area},
#'   \code{x}, and \code{y}.
#'
#' @param basemap A list containing the basemap information required
#'   for plotting, such as the raster image and its spatial bounds.
#'
#' @param area Character string specifying the area to display.
#'   Defaults to \code{"east"}.
#'
#' @param case_var Character string specifying the binary outcome
#'   column in \code{study}, where 1 indicates a case and 0
#'   indicates a negative observation. Defaults to \code{"case"}.
#'
#' @param k Integer specifying the minimum number of other cases
#'   required within the network neighbourhood for a case to be
#'   retained. Defaults to 5.
#'
#' @return
#' A list containing:
#' \describe{
#'   \item{\code{cases}}{A data frame containing cases that satisfy
#'   the neighbouring-case threshold, including their display
#'   coordinates.}
#'   \item{\code{negatives}}{A data frame containing all negative
#'   observations with available display coordinates in the
#'   selected area.}
#'   \item{\code{basemap}}{The supplied basemap object.}
#' }
#'
#' @note
#' This function returns observation-level data, including display
#' coordinates. The neighbouring-case threshold is not, by itself,
#' a sufficient disclosure-control mechanism for DataSHIELD.
#' The returned data must remain within an appropriately controlled
#' environment unless a separate disclosure assessment authorises
#' its release.
#'
#' @keywords internal
#'
#' @export
#' 
plotSpatialDS <- function(
    D,
    study,
    mapping,
    pos,
    basemap,
    area,
    case_var,
    k = 5
) {
  
  pos <- pos[pos$area == area, ]
  
  # Ensure study has matrix_row_id
  if (!"matrix_row_id" %in% names(study)) {
    study <- merge(
      study,
      mapping[, c("cell_id", "matrix_row_id"), drop = FALSE],
      by = "cell_id"
    )
  }
  
  # Add map coordinates
  study <- merge(
    study,
    pos[, c("matrix_row_id", "x", "y"), drop = FALSE],
    by = "matrix_row_id"
  )
  
  # 4. Suppress cases without sufficient nearby cases
  case_rows <- study$matrix_row_id[study[[case_var]] == 1]
  
  cases_per_row <- tabulate(
    case_rows,
    nbins = ncol(D)
  )
  
  cases_within_5km <- as.vector(
    (D != 0) %*% cases_per_row
  )
  
  other_cases <- cases_within_5km[case_rows] +
    cases_per_row[case_rows] - 1
  
  keep_case <- other_cases >= k
  
  cases <- study[study[[case_var]] == 1, ]
  cases <- cases[keep_case, ]
  
  negatives <- study[study[[case_var]] == 0, ]
  
  result <- list(cases = cases,
                 negatives = negatives,
                 basemap = basemap)
}




