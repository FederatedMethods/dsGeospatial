geospatialdev_1.2
================
2026-07-03

## Function summary

Plots heatmap for lsoa level aggregate of selected variables on map.
Uses modified ds.meanSdGp and MeanSdGpDS functions for lsoa level
aggrgation. Default map setting from boundr package, future version
accepts shapefile readable with st_read. Plotting done with library sf
and ggplot2.

Current testing done on synthetic health outcome data created by RVD.

## Upload data and load back with DSLite.

``` r
library(DSLite)
```

    ## Loading required package: DSI

    ## Loading required package: progress

    ## Loading required package: R6

    ## Loading required package: rly

``` r
library(devtools)
```

    ## Loading required package: usethis

``` r
library(dsBase)
```

    ## Registered S3 method overwritten by 'gamlss':
    ##   method   from
    ##   print.ri bit

``` r
library(dsBaseClient)
library(dsGeospatial)

dat <- qof_asthma_cohort_dist_gs_sf |> 
  dplyr::select("lsoa11cd", "has_asthma","uprn","distance_local_greenspace") |> 
  sf::st_drop_geometry()

devtools::load_all("/Users/tkariya/Documents/DataSHIELD/dsExample")
```

    ## ℹ Loading dsExample

``` r
#devtools::load_all("/Users/tkariya/Documents/DataSHIELD/dsExampleClient")


dslite.server1 <- newDSLiteServer(
  tables = list(
    uprn_asthma_green1 = dat
  )
)

dslite.server2 <- newDSLiteServer(
  tables = list(
    uprn_asthma_green2 = dat
  )
)


dslite.server1$config(defaultDSConfiguration(include=c("dsBase", "dsExample")))
dslite.server1$aggregateMethod("MeanSdGpDS", "MeanSdGpDS")

dslite.server2$config(defaultDSConfiguration(include=c("dsBase", "dsExample")))
dslite.server2$aggregateMethod("MeanSdGpDS", "MeanSdGpDS")


options(
  datashield.privacyControlLevel = "banana",
  nfilter.tab = 3,
  nfilter.subset = 3,
  nfilter.glm = 0.33,
  nfilter.string = 80,
  nfilter.stringShort = 20,
  nfilter.kNN = 3,
  nfilter.levels.density = 0.33,
  nfilter.levels.max = 40,
  nfilter.noise = 0.25,
  nfilter.privacy.old = 5
)


builder <- DSI::newDSLoginBuilder()

builder$append(
  server = "server_1",
  url = "dslite.server1",
  table = "uprn_asthma_green1",
  driver = "DSLiteDriver"
)

builder$append(
  server = "server_2",
  url = "dslite.server2",
  table = "uprn_asthma_green2",
  driver = "DSLiteDriver"
)


logindata <- builder$build()


conns <- DSI::datashield.login(logins = logindata, assign = FALSE)
```

    ## 
    ## Logging into the collaborating servers

``` r
datashield.assign.table(conns, "D", logindata)
```

Set nfilter.tab

``` r
nfilter.tab  <- 3
```

## Modified meanSDGPDS

modification: Convert object to numeric first and then find mean and
var.

``` r
## Server side function

MeanSdGpDS <- function (X, INDEX){
  
  FUN.mean <- function(x) {mean(as.numeric(x),na.rm=TRUE)} # changed lines 
  FUN.var <- function(x)  {stats::var(as.numeric(x),na.rm=TRUE)} # changed lines 
  
  analysis.matrix <- data.frame(X,INDEX) # changed lines 
  
  data.complete<-stats::complete.cases(analysis.matrix) 
  
  Ntotal<-dim(analysis.matrix)[1]
  Nmissing<-sum(!data.complete)
  Nvalid<-sum(data.complete)
  
  simplify<-TRUE
  
  analysis.matrix.no.miss<-analysis.matrix[data.complete,]
  nv<-dim(analysis.matrix)[2]
  
  X<-as.vector(analysis.matrix.no.miss[,1])
  INDEX<-analysis.matrix.no.miss[,2:nv]
  
  
  if (!is.list(INDEX)) 
    INDEX <- list(INDEX)
  nI <- length(INDEX)
  if (!nI) 
    stop("'INDEX' is of length zero")
  namelist <- vector("list", nI)
  names(namelist) <- names(INDEX)
  extent <- integer(nI)
  nx <- length(X)
  one <- 1L
  group <- rep.int(one, nx)
  ngroup <- one
  for (i in seq_along(INDEX)) {
    index <- as.factor(INDEX[[i]])
    if (length(index) != nx) 
      stop("arguments must have same length")
    namelist[[i]] <- levels(index)
    extent[i] <- nlevels(index)
    group <- group + ngroup * (as.integer(index) - one)
    ngroup <- ngroup * nlevels(index)
  }
  #    if (is.null(FUN.mean)) 
  #        return(group)
  
  #CALCULATE GROUP MEANS
  ans <- lapply(X = split(X, group), FUN = FUN.mean)
  index <- as.integer(names(ans))
  if (simplify && all(unlist(lapply(ans, length)) == 1L)) {
    ansmat <- array(dim = extent, dimnames = namelist)
    ans <- unlist(ans, recursive = FALSE)
  }
  else {
    ansmat <- array(vector("list", prod(extent)), dim = extent, 
                    dimnames = namelist)
  }
  if (length(index)) {
    names(ans) <- NULL
    ansmat[index] <- ans
  }
  ansmat.mean<-ansmat
  
  #CALCULATE GROUP SDs
  ans <- lapply(X = split(X, group), FUN = FUN.var)
  index <- as.integer(names(ans))
  if (simplify && all(unlist(lapply(ans, length)) == 1L)) {
    ansmat <- array(dim = extent, dimnames = namelist)
    ans <- unlist(ans, recursive = FALSE)
  }
  else {
    ansmat <- array(vector("list", prod(extent)), dim = extent, 
                    dimnames = namelist)
  }
  if (length(index)) {
    names(ans) <- NULL
    ansmat[index] <- ans
  }
  ansmat.sd<-sqrt(ansmat)
  
  
  #CALCULATE GROUP SIZES AND CHECK VALID
  
  ansmat.count<-table(group)
  
  # Set filter for cell sizes that are too small
  # the minimum number of observations that are allowed (the below function gets the value from opal)
  any.invalid.cell<-(sum(ansmat.count<nfilter.tab&ansmat.count>0)>=1)
  if(!any.invalid.cell)
  {
    table.valid<-TRUE
    cell.count.warning<-paste0("All tables valid") 
    result<-list(table.valid,ansmat.mean,ansmat.sd,ansmat.count,Nvalid,Nmissing,Ntotal,cell.count.warning)
    names(result)<-list("Table_valid","Mean_gp","StDev_gp", "N_gp","Nvalid","Nmissing","Ntotal","Message")
    return(result)
  }
  
  if(any.invalid.cell)
  {
    table.valid<-FALSE
    cell.count.warning<-paste0("At least one group has between 1 and ", nfilter.tab-1, " observations. Please change groups") 
    result<-list(table.valid,Nvalid,Nmissing,Ntotal,cell.count.warning)
    names(result)<-list("Table_valid","Nvalid","Nmissing","Ntotal","Warning")
    return(result)
  }
  
}
```

## geoheatmap function

Inputs : inputs of ds.meanSdGp and shapefile and names_region (when
using boundr) Return : mean plot of selected variable. mean can be
‘split’ or ‘combined’ across studies.

``` r
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



ds.geoheatmapPlot <- function(x = NULL, y = NULL, type = 'combine', 
                              do.checks = FALSE, stat = "mean",
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
      sd.matrix <- NULL
      Nvalid <- 0
      Nmissing <- 0
      Ntotal <- 0
      
      for(j in 1:numsources){
        mean.matrix <- rbind(mean.matrix,as.numeric(unlist(output[[j]][2])))
        n.matrix <- rbind(n.matrix,as.numeric(unlist(output[[j]][4])))
        sd.matrix <- rbind(sd.matrix,as.numeric(unlist(output[[j]][3])))
        Nvalid <- Nvalid+as.numeric(unlist(output[[j]][5]))
        Nmissing <- Nmissing+as.numeric(unlist(output[[j]][6]))
        Ntotal <- Ntotal+as.numeric(unlist(output[[j]][7]))
      }
    
      var.matrix <- sd.matrix^2
      nsum.vector <- rep(1,numsources)
      # Calculate weighted means across studies in each group
      mean.gp <- (diag(t(mean.matrix)%*%n.matrix))/(t(n.matrix)%*%nsum.vector)
      var.gp <- (diag(t(var.matrix)%*%n.matrix))/(t(n.matrix)%*%nsum.vector)
      SD.gp <- sqrt(var.gp)
      N.gp <- (t(n.matrix) %*% nsum.vector)
      SEM.gp <- SD.gp/sqrt(N.gp)
      
      mean.data <- data.frame(lsoa11cd = lsoanames,
                              mean.study = as.numeric(mean.gp),
                              study = 'combine')
      
      var.data <- data.frame(lsoa11cd = lsoanames,
                             mean.study = as.numeric(var.gp),
                             study = 'combine')
      
      SD.data <- data.frame(lsoa11cd = lsoanames,
                           mean.study = as.numeric(SD.gp),
                           study = 'combine')
      
      N.data <- data.frame(lsoa11cd = lsoanames,
                           mean.study = as.numeric(N.gp),
                           study = 'combine')
      
      SEM.data <- data.frame(lsoa11cd = lsoanames,
                             mean.study = as.numeric(SEM.gp),
                             study = 'combine')

      shape_sf <- shape_list[[1]] |>
        dplyr::left_join(get(paste0(stat,".data")), by = "lsoa11cd") |>
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
  
  ggplot2::ggplot(plot.matrix) +
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
  
}
```

## Outputs

``` r
ds.geoheatmapPlot(x="D$has_asthma", y = "D$lsoa11cd", datasources = conns)
```

![](geospatial_dev_1.2_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
ds.geoheatmapPlot(x="D$has_asthma", y = "D$lsoa11cd", datasources = conns, 
                  type='split')
```

![](geospatial_dev_1.2_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->

``` r
ds.geoheatmapPlot(x="D$has_asthma", y = "D$lsoa11cd", datasources = conns,
                  type='split', names_region = c("Liverpool", "Sefton"))
```

![](geospatial_dev_1.2_files/figure-gfm/unnamed-chunk-6-3.png)<!-- -->

``` r
ds.geoheatmapPlot(x="D$distance_local_greenspace", y = "D$lsoa11cd", 
                  datasources = conns)
```

![](geospatial_dev_1.2_files/figure-gfm/unnamed-chunk-6-4.png)<!-- -->

``` r
ds.geoheatmapPlot(x="D$has_asthma", y = "D$lsoa11cd", datasources = conns,
                  type = 'split', names_region = c("Liverpool", "Sefton"))
```

![](geospatial_dev_1.2_files/figure-gfm/unnamed-chunk-6-5.png)<!-- -->
