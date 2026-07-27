
## Server side function
#' Calculate grouped summary statistics
#'
#' Server-side aggregate function that computes the mean, standard deviation,
#' and sample size of a numeric variable grouped by one or more categorical
#' variables. Missing observations are removed prior to analysis and DataSHIELD
#' disclosure controls are applied to ensure that group sizes satisfy the
#' configured minimum cell count threshold (`nfilter.tab`).
#'
#' @param X A numeric vector.
#' @param INDEX A factor, character, or list of grouping variables of the same
#'   length as `X`.
#'
#' @return A named list.
#' If all groups satisfy the disclosure threshold, the returned list contains:
#' \describe{
#'   \item{Table_valid}{Logical indicating whether all groups passed disclosure checks.}
#'   \item{Mean_gp}{Array of group means.}
#'   \item{StDev_gp}{Array of group standard deviations.}
#'   \item{N_gp}{Table of group sample sizes.}
#'   \item{Nvalid}{Number of complete observations used in the analysis.}
#'   \item{Nmissing}{Number of observations removed due to missing values.}
#'   \item{Ntotal}{Total number of observations supplied.}
#'   \item{Message}{Disclosure status message.}
#' }
#'
#' If one or more groups fail the disclosure threshold, only the validity flag,
#' observation counts and warning message are returned.
#'
#' @details
#' This function is intended for server-side execution within the DataSHIELD
#' framework and should not normally be called directly by end users.
#'
#' @keywords internal
#' @export
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

