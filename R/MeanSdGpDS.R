#' Calculate Grouped Summary Statistics on the DataSHIELD Server
#'
#' @description
#' Computes summary statistics for a numerical variable grouped by one
#' or more categorical variables on the DataSHIELD server. The function
#' calculates group means, standard deviations, quantiles and sample
#' sizes while enforcing DataSHIELD disclosure-control thresholds.
#'
#' @details
#' This server-side function is designed to be called by the client-side
#' function \code{ds.summarystat()} through a DataSHIELD aggregate request.
#'
#' The function first removes observations with missing values in either
#' the numerical variable or any grouping variable. It then calculates
#' the following statistics for each group:
#'
#' \itemize{
#'   \item Arithmetic mean.
#'   \item Standard deviation (SD).
#'   \item Quantiles at 5%, 10%, 25%, 50%, 75%, 90% and 95%.
#'   \item Number of observations.
#' }
#'
#' When multiple grouping variables are supplied, statistics are
#' calculated for each combination of their factor levels.
#'
#' Before returning group-level statistics, the function checks
#' disclosure-control settings using
#' \code{dsBase::listDisclosureSettingsDS()}. The minimum permitted
#' cell size is determined by \code{nfilter.tab}.
#'
#' If any non-empty group contains fewer observations than
#' \code{nfilter.tab}, the function suppresses all group-level
#' statistics and returns a warning message instead.
#'
#' If all non-empty groups satisfy the disclosure-control threshold,
#' the calculated summary statistics are returned to the client.
#'
#' @param X A numerical vector containing the variable for which
#'   grouped summary statistics are calculated. Values are coerced
#'   to numeric during calculation.
#'
#' @param INDEX A grouping variable or a data frame containing
#'   multiple grouping variables. Each grouping variable is converted
#'   to a factor. When multiple variables are supplied, groups are
#'   defined by combinations of their factor levels.
#'
#' @return
#' A named list whose contents depend on whether the
#' disclosure-control check is satisfied.
#'
#' If all group counts are valid, the returned list contains:
#' \describe{
#'   \item{Table_valid}{Logical value \code{TRUE}, indicating that
#'     all non-empty groups satisfy the disclosure threshold.}
#'   \item{Mean_gp}{Array containing the mean for each group.}
#'   \item{StDev_gp}{Array containing the standard deviation
#'     for each group.}
#'   \item{Q_gp}{Array or list-array containing the seven quantiles
#'     for each group.}
#'   \item{N_gp}{Frequency table containing group counts.}
#'   \item{Nvalid}{Total number of complete observations.}
#'   \item{Nmissing}{Total number of observations with missing values
#'     in X or INDEX.}
#'   \item{Ntotal}{Total number of input observations.}
#'   \item{Message}{Character string indicating that all tables
#'     passed disclosure-control checks.}
#' }
#'
#' If at least one non-empty group fails the disclosure-control
#' check, the returned list contains:
#' \describe{
#'   \item{Table_valid}{Logical value \code{FALSE}.}
#'   \item{Nvalid}{Total number of complete observations.}
#'   \item{Nmissing}{Total number of incomplete observations.}
#'   \item{Ntotal}{Total number of input observations.}
#'   \item{Warning}{Character string describing the disclosure-control
#'     violation and requesting regrouping.}
#' }
#'
#'
#' @seealso
#' \code{\link{ds.summarystat}}
#'
#' @export
MeanSdGpDS <- function (X, INDEX){
  
  
  #############################################################
  # MODULE 1: CAPTURE THE nfilter SETTINGS
  thr <- dsBase::listDisclosureSettingsDS()
  nfilter.tab <- as.numeric(thr$nfilter.tab)
  #nfilter.glm <- as.numeric(thr$nfilter.glm)
  #nfilter.subset <- as.numeric(thr$nfilter.subset)
  #nfilter.string <- as.numeric(thr$nfilter.string)
  #############################################################

  FUN.mean <- function(x) {mean(as.numeric(x),na.rm=TRUE)} # changed lines
  FUN.var <- function(x)  {stats::var(as.numeric(x),na.rm=TRUE)} # changed lines
  FUN.qq <- function(x)  {stats::quantile(as.numeric(x), 
                                          c(0.05,0.1,0.25,0.5,0.75,0.9,0.95), 
                                          na.rm=TRUE)}

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
  
 
  #CALCULATE GROUP QUANTILES
  
  ans <- lapply(X = split(X, group), FUN = FUN.qq)
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
  ansmat.qq<-ansmat

  #CALCULATE GROUP SIZES AND CHECK VALID

  ansmat.count<-table(group)

  # Set filter for cell sizes that are too small
  # the minimum number of observations that are allowed (the below function gets the value from opal)
  any.invalid.cell<-(sum(ansmat.count<nfilter.tab&ansmat.count>0)>=1)

  if(!any.invalid.cell)
  {
    table.valid<-TRUE
    cell.count.warning<-paste0("All tables valid")
    result<-list(table.valid,ansmat.mean,ansmat.sd,ansmat.qq, ansmat.count,Nvalid,Nmissing,Ntotal,cell.count.warning)
    names(result)<-list("Table_valid","Mean_gp","StDev_gp", "Q_gp","N_gp","Nvalid","Nmissing","Ntotal","Message")
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

