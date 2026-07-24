#' @title Calculates the correlation of R objects in the server-side
#' @description This function calculates the correlation of two variables or the correlation
#' matrix for the variables of an input data frame.
#' @details In addition to computing correlations; this function produces a table outlining the
#' number of complete cases and a table outlining the number of missing values to allow the
#' user to decide the 'relevance' of the correlation based on the number of complete
#' cases included in the correlation calculations.
#'
#' If the argument \code{y} is not NULL, the dimensions of the object have to be
#' compatible with the argument \code{x}.
#'
#' The function calculates the pairwise correlations based on casewise complete cases which means that
#' it omits all the rows in the input data frame that include at least one cell with a missing value,
#' before the calculation of correlations.
#'
#' If \code{type} is set to \code{'split'} (default), the correlation of two variables or the
#' variance-correlation matrix of an input data frame and the number of complete cases and missing
#' values are returned for every single study. If type is set to \code{'combine'}, the pooled
#' correlation, the total number of complete cases and the total number of missing values aggregated
#' from all the involved studies, are returned.
#'
#' Server function called: \code{corDS}
#'
#' @param x a character string providing the name of the input vector, data frame or matrix.
#' @param y a character string providing the name of the input vector, data frame or matrix.
#' Default NULL.
#' @param type a character string that represents the type of analysis to carry out.
#' This must be set to \code{'split'} or \code{'combine'}.  Default \code{'split'}. For more information see details.
#' @param datasources a list of \code{\link[DSI]{DSConnection-class}} objects obtained after login.
#' If the \code{datasources} argument is not specified
#' the default set of connections will be used: see \code{\link[DSI]{datashield.connections_default}}.
#' @return \code{ds.cor} returns a list containing the number of missing values in each variable,
#' the number of missing variables casewise, the correlation matrix,
#' the number of used complete cases. The function applies two disclosure controls. The first disclosure
#' control checks that the number of variables is not bigger than a percentage of the individual-level records (the allowed
#' percentage is pre-specified by the 'nfilter.glm'). The second disclosure control checks that none of them is dichotomous
#' with a level having fewer counts than the pre-specified 'nfilter.tab' threshold.
#' @author DataSHIELD Development Team
#' @examples
#' \dontrun{
#'
#' ## Version 6, for version 5 see the Wiki
#'   # Connecting to the Opal servers
#'
#'   require('DSI')
#'   require('DSOpal')
#'   require('dsBaseClient')
#'
#'   builder <- DSI::newDSLoginBuilder()
#'   builder$append(server = "study1",
#'                  url = "http://192.168.56.100:8080/",
#'                  user = "administrator", password = "datashield_test&",
#'                  table = "CNSIM.CNSIM1", driver = "OpalDriver")
#'   builder$append(server = "study2",
#'                  url = "http://192.168.56.100:8080/",
#'                  user = "administrator", password = "datashield_test&",
#'                  table = "CNSIM.CNSIM2", driver = "OpalDriver")
#'   builder$append(server = "study3",
#'                  url = "http://192.168.56.100:8080/",
#'                  user = "administrator", password = "datashield_test&",
#'                  table = "CNSIM.CNSIM3", driver = "OpalDriver")
#'   logindata <- builder$build()
#'
#'   # Log onto the remote Opal training servers
#'   connections <- DSI::datashield.login(logins = logindata, assign = TRUE, symbol = "D")
#'
#'   # Example 1: Get the correlation matrix of two continuous variables
#'   ds.cor(x="D$LAB_TSC", y="D$LAB_TRIG", type="combine", datasources = connections)
#'
#'   # Example 2: Get the correlation matrix of the variables in a dataframe
#'   ds.dataFrame(x=c("D$LAB_TSC", "D$LAB_TRIG", "D$LAB_HDL", "D$PM_BMI_CONTINUOUS"),
#'                newobj="D.new", check.names=FALSE, datasources=connections)
#'   ds.cor("D.new", type="combine", datasources = connections)
#'
#'   # clear the Datashield R sessions and logout
#'   datashield.logout(connections)
#'
#' }
#' @export
#'
ds.groupcor <- function(x=NULL, y= NULL, index = NULL, type="split", datasources=NULL,
                        do.checks = TRUE, names_region = "Liverpool"){

  # look for DS connections
  if(is.null(datasources)){
    datasources <- datashield.connections_find()
  }

  # ensure datasources is a list of DSConnection-class
  if(!(is.list(datasources) && all(unlist(lapply(datasources, function(d) {methods::is(d,"DSConnection")}))))){
    stop("The 'datasources' were expected to be a list of DSConnection-class objects", call.=FALSE)
  }

  if(is.null(x)){
    stop("x=NULL. Please provide the name of a matrix or dataframe or the names of two numeric vectors!", call.=FALSE)
  }


  if(is.null(y)){
    stop("Please provide the name of the input vector!", call.=FALSE)
  }

  if(is.null(index)){
    stop("Please provide the name of the input vector!", call.=FALSE)
  }

  # check the type of the input objects
  if(do.checks){

    # check if the input objects are defined in all the studies
    isDefined(datasources, x)
    isDefined(datasources, y)

    # call the internal function that checks the input object is of the same class in all studies.
    typ1 <- checkClass(datasources, x)
    typ2 <- checkClass(datasources, y)
  }



  calltext <- paste0("groupcorDS(", x, ",", y, ",",index,")")

  output <- DSI::datashield.aggregate(datasources, as.symbol(calltext))
  nfilter.tab <- output[[1]][["Nfilter.tab"]]
  output <- lapply(output, function(x) {
    x[-length(x)]
  })

  # name of the studies to be used in the output
  stdnames <- names(datasources)



  ## Loading map
  if(type == 'split' && length(names_region) != length(stdnames)){
    names_region <- rep(names_region, times = length(stdnames))
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



  # call the server side function
  if (type=="split"){
    names(shape_list) <- stdnames
    covariance <- list()
    sqrt.diag <- list()
    correlation <- list()
    splitresults <- list()
    results <- list()
    for(i in 1:length(output)){
      for(j in 1: length(output[[i]])){ # -1 to exclude [[nfiltertab]]
        covariance[[j]] <- matrix(0, ncol=dim(output[[i]][[j]][[1]])[2], nrow=dim(output[[i]][[j]][[1]])[1])
        correlation[[j]] <- matrix(0, ncol=dim(output[[i]][[j]][[1]])[2], nrow=dim(output[[i]][[j]][[1]])[1])
        colnames(correlation[[j]]) <- colnames(output[[i]][[j]][[1]])
        rownames(correlation[[j]]) <- colnames(output[[i]][[j]][[1]])
        for(m in 1:dim(output[[i]][[j]][[1]])[1]){
          for(n in 1:dim(output[[i]][[j]][[1]])[2]){
            covariance[[j]][m,n] <- (1/(output[[i]][[j]][[3]][m,n]-1))*(output[[i]][[j]][[1]][m,n])-(1/(output[[i]][[j]][[3]][m,n]*(output[[i]][[j]][[3]][m,n]-1)))*output[[i]][[j]][[2]][m,n]*output[[i]][[j]][[2]][n,m]
          }
        }
        correlation[[j]] <- stats::cov2cor(covariance[[j]])
        splitresults[[j]] <- list(output[[i]][[j]][[4]][[1]], output[[i]][[j]][[4]][[2]], correlation[[j]], output[[i]][[j]][[3]])
        n1 <- "Number of missing values in each variable"
        n2 <- "Number of missing values casewise"
        n3 <- "Correlation Matrix"
        n4 <- "Number of complete cases used"
        names(splitresults[[j]]) <- c(n1, n2, n3, n4)
      }
      names(splitresults) <- names(output[[i]])
      results[[i]] <- splitresults
    }
    names(results) <- stdnames

    results_data <- do.call(
      rbind,
      lapply(names(results), function(server_name) {

        server_results <- results[[server_name]]

        correlation_data <- data.frame(
          lsoa11cd = names(server_results),

          complete.counts = vapply(
            server_results,
            function(x) {
              x[["Number of complete cases used"]][1, 2]
            },
            numeric(1)
          ),

          correlation.coef = vapply(
            server_results,
            function(x) {
              x[["Correlation Matrix"]][1, 2]
            },
            numeric(1)
          ),

          server = server_name,
          row.names = NULL
        )

        dplyr::left_join(
          shape_list[[server_name]],
          correlation_data,
          by = "lsoa11cd"
        )
      })
    )

    results_data <- sf::st_as_sf(results_data)
    rownames(results_data) <- NULL

  }
  else if (type=="combine"){
      common_lsoas <- sort(Reduce(intersect, lapply(output, names)))

      if(is.null(common_lsoas) || length(common_lsoas) < output[[1]][["Nfilter.tab"]]){
        stop("Common lsoa in servers is less than the safe limit")
      }
      # lsoa in all servers should match or warning that common lsoa is only analysed
      # in that case check if the matrix size is not disclosive
      combined.results <- list()
      for (lsoa in common_lsoas) {
        reference <- output[[1]][[lsoa]]

         combined.sums.of.products <- matrix(
            0,
            nrow = nrow(reference[[1]]),
            ncol = ncol(reference[[1]]),
            dimnames = dimnames(reference[[1]])
          )

          combined.sums <- matrix(
            0,
            nrow = nrow(reference[[2]]),
            ncol = ncol(reference[[2]]),
            dimnames = dimnames(reference[[2]])
          )

          combined.complete.cases <- matrix(
            0,
            nrow = nrow(reference[[3]]),
            ncol = ncol(reference[[3]]),
            dimnames = dimnames(reference[[3]])
          )

          combined.missing.cases.vector <- matrix(
            0,
            nrow = nrow(reference[[4]][[1]]),
            ncol = ncol(reference[[4]][[1]]),
            dimnames = dimnames(reference[[4]][[1]])
          )

          combined.missing.cases.matrix <- matrix(
            0,
            nrow = nrow(reference[[4]][[2]]),
            ncol = ncol(reference[[4]][[2]]),
            dimnames = dimnames(reference[[4]][[2]])
          )

          combined.sums.of.squares <- matrix(
            0,
            nrow = nrow(reference[[5]]),
            ncol = ncol(reference[[5]]),
            dimnames = dimnames(reference[[5]])
          )

          # Combine sufficient statistics across relevant servers
          for (i in stdnames) {

            lsoa_output <- output[[i]][[lsoa]]

            combined.sums.of.products <-
              combined.sums.of.products + lsoa_output[[1]]

            combined.sums <-
              combined.sums + lsoa_output[[2]]

            combined.complete.cases <-
              combined.complete.cases + lsoa_output[[3]]

            combined.missing.cases.vector <-
              combined.missing.cases.vector + lsoa_output[[4]][[1]]

            combined.missing.cases.matrix <-
              combined.missing.cases.matrix + lsoa_output[[4]][[2]]

            combined.sums.of.squares <-
              combined.sums.of.squares + lsoa_output[[5]]
          }

          # Calculate covariance matrix for this LSOA
          combined.covariance <- matrix(
            NA_real_,
            nrow = nrow(combined.sums.of.products),
            ncol = ncol(combined.sums.of.products),
            dimnames = dimnames(combined.sums.of.products)
          )

          for (m in seq_len(nrow(combined.sums.of.products))) {
            for (n in seq_len(ncol(combined.sums.of.products))) {

              n_complete <- combined.complete.cases[m, n]

              if (!is.na(n_complete) && n_complete > 1) {

                combined.covariance[m, n] <-
                  combined.sums.of.products[m, n] / (n_complete - 1) -
                  (
                    combined.sums[m, n] *
                      combined.sums[n, m]
                  ) / (
                    n_complete * (n_complete - 1)
                  )
              }
            }
          }

          # Convert covariance to correlation
          combined.correlation <- tryCatch(
            stats::cov2cor(combined.covariance),
            error = function(e) {
              matrix(
                NA_real_,
                nrow = nrow(combined.covariance),
                ncol = ncol(combined.covariance),
                dimnames = dimnames(combined.covariance)
              )
            }
          )

          combined.results[[lsoa]] <- list(
            "Number of missing values in each variable" =
              combined.missing.cases.vector,

            "Number of missing values casewise" =
              combined.missing.cases.matrix,

            "Number of complete cases used" =
              combined.complete.cases,

            "Covariance Matrix" =
              combined.covariance,

            "Correlation Matrix" =
              combined.correlation
          )
        }

        results <- combined.results

        correlation_data <- data.frame(
          lsoa11cd = names(results),
          complete.counts = sapply(
            results,
            function(x) x$`Number of complete cases used`[1, 2]
          ),
          correlation.coef = sapply(
            results,
            function(x) x$`Correlation Matrix`[1, 2]
          ),
          server = if (type == "split") names(results) else rep("combine", length(results)),
          row.names = NULL
        )

        results_data <- dplyr::left_join(
            shape_list[[1]],
            correlation_data,
            by = "lsoa11cd"
          ) %>% dplyr::mutate(server = "combine")
      }

    else{
      stop('Function argument "type" has to be either "combine" or "split"', call.=FALSE)
    }

# Plotting
##### Load shapefile
 corplot <- ggplot2::ggplot(results_data) +
      ggplot2::geom_sf(ggplot2::aes(fill = correlation.coef)) +
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
      ggplot2::labs(fill = "Correlation Coefficient") +
      ggplot2::facet_wrap(~server) +
      ggplot2::theme_minimal()


return(corplot)
}



