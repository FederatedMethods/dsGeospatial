## Client-side function
#'
#' @title Plot a geographic heatmap
#' @description Creates a geographic heatmap showing the group-wise mean of a
#' numeric variable across Lower Layer Super Output Areas (LSOAs).
#'
#' @details This function calls the server-side function \code{MeanSdGpDS} to
#' calculate the mean of \code{x} within groups defined by \code{y}. The grouping
#' variable is expected to contain LSOA codes, such as \code{lsoa11cd}. These
#' results are then joined to LSOA boundary data obtained using \code{boundr} and
#' plotted as a choropleth map using \code{ggplot2::geom_sf}.
#'
#' Depending on \code{type}, the heatmap can be produced in different ways:
#' \cr
#' (1) \code{"combine"}: pooled results are calculated across all studies and
#' plotted as a single combined map. This requires the LSOA levels to be identical
#' across all studies. \cr
#' (2) \code{"split"}: study-specific maps are produced separately and displayed
#' using facets. \cr
#' (3) \code{"both"}: both combined and study-specific outputs are prepared where
#' applicable.
#'
#' For \code{type = "combine"}, \code{names_region} must be of length one. For
#' \code{type = "split"}, \code{names_region} should correspond to the regions
#' required for each study.
#'
#' Disclosure control is applied on the server side. If any table contains a cell
#' count between 1 and the minimum filter threshold, the function stops and returns
#' a warning message asking the user to regroup the data.
#'
#' Server function called: \code{MeanSdGpDS}
#'
#' @param x A character string specifying the name of the numeric variable for
#' which the group-wise mean is calculated.
#' @param y A character string specifying the name of the grouping variable.
#' This should usually contain geographic area identifiers such as LSOA 2011 codes.
#' @param type A character string specifying the type of output to produce.
#' This can be \code{"combine"}, \code{"split"} or \code{"both"}.
#' Default is \code{"combine"}.
#' @param do.checks Logical. If \code{TRUE}, checks are performed to confirm that
#' the input objects are defined in all studies and are of equivalent class across
#' studies. Default is \code{FALSE} to save time.
#' @param names_region A character vector specifying the local authority district
#' name or names used to obtain LSOA boundary data from \code{boundr}. Default is
#' \code{"Liverpool"}.
#' @param shapefile Currently unused. Reserved for future support for a user-supplied
#' spatial object or shapefile.
#' @param datasources A list of \code{\link[DSI]{DSConnection-class}} objects
#' obtained after login. If \code{NULL}, the function uses the currently available
#' DataSHIELD connections returned by \code{datashield.connections_find()}.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object containing a geographic heatmap.
#' The map shows the group-wise mean of \code{x} for each LSOA, either pooled
#' across studies or split by study depending on \code{type}.
#'
#' @importFrom grDevices colorRampPalette
#' @export
#'
#' @examples
#' \dontrun{
#'   library(DSI)
#'   library(DSOpal)
#'   library(dsBaseClient)
#'
#'   builder <- DSI::newDSLoginBuilder()
#'   builder$append(server = "study1",
#'                  url = "https://opal-demo.obiba.org",
#'                  user = "administrator",
#'                  password = "password",
#'                  table = "UPRN.synthetic_asthma_uprn_level",
#'                  driver = "OpalDriver")
#'
#'   logindata <- builder$build()
#'
#'   connections <- DSI::datashield.login(logins = logindata,
#'                                        assign = TRUE,
#'                                        symbol = "D")
#'
#'   ds.geoheatmapPlot(
#'     x = "D$has_asthma",
#'     y = "D$lsoa11cd",
#'     type = "combine",
#'     names_region = "Liverpool",
#'     datasources = connections
#'   )
#'
#'   datashield.logout(connections)
#' }
#' @importFrom grDevices colorRampPalette
#' @export



ds.geoheatmapPlot <- function(x=NULL, y=NULL, type='combine', do.checks=FALSE,
                              names_region = "Liverpool", shapefile = NULL,
                              datasources=NULL){

  # look for DS connections
  # look for DS connections
  if(is.null(datasources)){
    datasources <- datashield.connections_find()
  }

  # ensure datasources is a list of DSConnection-class
  if(!(is.list(datasources) && all(unlist(lapply(datasources, function(d) {methods::is(d,"DSConnection")}))))){
    stop("The 'datasources' were expected to be a list of DSConnection-class objects", call.=FALSE)
  }

  if(is.null(x)){
    stop("Please provide the name of the input vector!", call.=FALSE)
  }

  if(is.null(y)){
    stop("Please provide the name of the input vector!", call.=FALSE)
  }

  if(do.checks){

    # check if the input objects are defined in all the studies
    isDefined(datasources, x)
    isDefined(datasources, y)

    # call the internal function that checks the input object is of the same class in all studies.
    typ1 <- checkClass(datasources, x)
    typ2 <- checkClass(datasources, y)
  }

  xvarname <- gsub("_", " ", sub(".*\\$", "", x))


  # call the server side function that calculates mean and standard deviation
  # by group in each study

  calltext <- paste0("MeanSdGpDS(", x, ",", y, ")")
  output <- DSI::datashield.aggregate(datasources, as.symbol(calltext))


  # Work around : introducing output of server side directly
  # output$UPRN <- UPRN
  #

  numsources <- length(output)

  all.tables.valid<-1

  for(i in 1:numsources){
    if(unlist(output[[i]]$Table_valid==FALSE)) all.tables.valid <- all.tables.valid-1
  }

  if(all.tables.valid<1){
    warning.message<-"At least 1 cell count is 1-nfilter, please regroup"
    return(warning.message)
  }

  plot.matrix <- list()

  ##### Load shapefile
  # if type is combine names_region can be only of length 1

  if(type == 'split' && length(names_region) != numsources){
    names_region <- rep(names_region, times = numsources)
  }

  shape_list <- lapply(names_region, function(region) {
    boundr::bounds(
      "lsoa",
      within_level = "lad",
      within_names = region,
      lookup_year = 2011,
      opts = boundr::boundr_options(resolution = "BFC")
    ) |>
      dplyr::select(lsoa11cd, geometry)
  })


  if(all.tables.valid==1){
    # COMBINE OR BOTH
    if(type=="combine"){

      if(length(names_region)!=1){
        stop("expected length of names_region is 1")
      }

      # if type is combine check whether lsoa levels are same in the two studies
       lsoa_names <- lapply(output, function(x) names(x$Mean_gp))
       names.identical <- all(vapply(lsoa_names, identical, logical(1),y = lsoa_names[[1]]))


       if(names.identical) {
         lsoanames <- lsoa_names[[1]]
       } else {
         stop("lsoa levels across studies not sae, can't proceed with type = 'combine'")
       }

      numsources <- length(output)
      mean.matrix <- NULL
      n.matrix <- NULL

      for(j in 1:numsources){
        mean.matrix <- rbind(mean.matrix,as.numeric(unlist(output[[j]][2])))
        n.matrix <- rbind(n.matrix,as.numeric(unlist(output[[j]][4])))
      }

      nsum.vector <- rep(1,numsources)
      # Calculate weighted means across studies in each group
      mean.gp <- (diag(t(mean.matrix)%*%n.matrix))/(t(n.matrix)%*%nsum.vector)

      mean.data <- data.frame(lsoa11cd = lsoanames,
      mean.study = as.numeric(mean.gp),
      study = 'combine')

      shape_sf <- shape_list[[1]] |>
        dplyr::left_join(mean.data, by = "lsoa11cd") |>
        sf::st_as_sf()

      plot.matrix <- shape_sf

    }

  ##### Type Split

  if(type=="split"){

    numsources <- length(output)
    mean.matrix <- NULL

    for(j in 1:numsources){

      mean.data <- data.frame(
        lsoa11cd = names(output[[j]][2]$Mean_gp),
        mean.study = as.numeric(output[[j]][2]$Mean_gp),
        study = names(datasources[j])
      )

      shape_sf <- shape_list[[j]] |>
        dplyr::left_join(mean.data, by = "lsoa11cd") |>
        sf::st_as_sf()

      plot.matrix[[j]] <- shape_sf

    }
    plot.matrix <- do.call(rbind, plot.matrix)
  }

  }


  # Plotting

plotresult <-  ggplot2::ggplot(plot.matrix) +
      ggplot2::geom_sf(ggplot2::aes(fill = mean.study)) +
      ggplot2::scale_fill_gradientn(
        colours = grDevices::colorRampPalette(c(
          "#440154",
          "#414487",
          "#2A788E",
          "#22A884",
          "#7AD151",
          "#FDE725"
        ))(100),
        na.value = "grey90"
      ) +
      ggplot2::labs(fill = xvarname) +
      ggplot2::facet_wrap(~study) +
      ggplot2::theme_minimal()

return(plotresult)
}

