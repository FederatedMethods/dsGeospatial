# correlation matrix

#' DSLite config
#'#' Calculate grouped correlation summary statistics
#'
#' Server-side aggregate function that computes the summary statistics required
#' to construct covariance and correlation matrices within each level of one or
#' more grouping variables. The function returns group-specific sums, sums of
#' products, sums of squares, complete-case counts and missing-value summaries,
#' allowing the client to calculate covariance and correlation coefficients
#' without disclosing individual-level data.
#'
#' Disclosure controls are applied using the DataSHIELD disclosure settings.
#' The function checks for oversaturated covariance matrices, binary variables
#' with insufficient cell counts, and grouped tables with counts below the
#' configured disclosure threshold (`nfilter.tab`).
#'
#' @param X A numeric vector.
#' @param Y A numeric vector.
#' @param INDEX A factor, character vector, or list of grouping variables of the
#'   same length as `X` and `Y`.
#'
#' @return
#' If all disclosure checks are satisfied, a named list containing one element
#' per group. Each group contains:
#' \describe{
#'   \item{sums.of.products}{Matrix of pairwise sums of products.}
#'   \item{sums}{Matrix of variable sums.}
#'   \item{complete.counts}{Matrix of complete-case counts.}
#'   \item{na.counts}{List containing column-wise and case-wise missing-value counts.}
#'   \item{sums.of.squares}{Matrix of sums of squares.}
#' }
#' In addition, the returned list contains:
#' \describe{
#'   \item{Nfilter.tab}{The disclosure threshold used for grouped tables.}
#' }
#'
#' If disclosure checks fail, a reduced list is returned containing:
#' \describe{
#'   \item{Table_valid}{Logical indicating whether disclosure checks passed.}
#'   \item{Nvalid}{Number of complete observations.}
#'   \item{Nmissing}{Number of incomplete observations.}
#'   \item{Ntotal}{Total number of observations.}
#'   \item{Warning}{Disclosure warning message.}
#' }
#'
#' @details
#' This function is intended for server-side execution within the DataSHIELD
#' framework and should not normally be called directly by end users.
#'
#' @keywords internal
#' @export
groupcorDS <- function (X, Y, INDEX){

  #############################################################
  #MODULE 1: CAPTURE THE nfilter SETTINGS
  thr <- dsBase::listDisclosureSettingsDS()
  nfilter.tab <- as.numeric(thr$nfilter.tab)
  nfilter.glm <- as.numeric(thr$nfilter.glm)
  #############################################################

  ## define corDS function

  FUN <- function(x){
    X <- x$X
    Y <- x$Y
    dataframe <- as.data.frame(cbind(X,Y))


    # names of the variables
    cls <- colnames(dataframe)

    # number of the input variables
    N.vars <- ncol(dataframe)

    ######################
    # DISCLOSURE CONTROLS
    ######################

    ##############################################################
    # FIRST TYPE OF DISCLOSURE TRAP - TEST FOR OVERSATURATION
    # TEST AGAINST nfilter.glm
    ##############################################################

    varcov.saturation.invalid <- 0

    if(N.vars > (nfilter.glm * nrow(dataframe))){

      varcov.saturation.invalid <- 1

      studysideMessage <- "ERROR: The ratio of the number of variables over the number of individual-level
                          records exceeds the allowed threshold, there is a possible risk of disclosure"
      stop(studysideMessage, call. = FALSE)

    }

    # CHECK X MATRIX VALIDITY
    # Check no dichotomous X vectors with between 1 and filter.threshold
    # observations at either level

    X.mat <- as.matrix(dataframe)

    dimX <- dim(X.mat)

    num.Xpar <- dimX[2]

    Xpar.invalid <- rep(0, num.Xpar)

    for(pj in 1:num.Xpar){
      unique.values.noNA <- unique((X.mat[,pj])[stats::complete.cases(X.mat[,pj])])
      if(length(unique.values.noNA)==2){
        tabvar <- table(X.mat[,pj])[table(X.mat[,pj])>=1] #tabvar COUNTS N IN ALL CATEGORIES WITH AT LEAST ONE OBSERVATION
        min.category <- min(tabvar)
        if(min.category < nfilter.tab){
          Xpar.invalid[pj] <- 1
        }
      }
    }

    # if any of the vectors in X matrix is invalid then the function returns an error

    if(is.element('1', Xpar.invalid)==TRUE & varcov.saturation.invalid==0){

      studysideMessage <- "ERROR: at least one variable is binary with one category less than the filter threshold for table cell size"
      stop(studysideMessage, call. = FALSE)

    }

    # if all vectors in X matrix are valid then the output matrices are calculated

    if(is.element('1', Xpar.invalid)==FALSE & varcov.saturation.invalid==0){

      # calculate the number of NAs in each variable separately
      column.NAs <- matrix(ncol=N.vars, nrow=1)
      colnames(column.NAs) <- cls
      for(i in 1:N.vars){
        column.NAs[1,i] <- length(dataframe[,i])-length(dataframe[stats::complete.cases(dataframe[,i]),i])
      }

      # remove any rows from the dataframe that include NAs
      casewise.dataframe <- dataframe[stats::complete.cases(dataframe),]

      # calculate the number of NAs casewise
      casewise.NAs.all <- dim(dataframe)[1]-dim(casewise.dataframe)[1]
      casewise.NAs <- matrix(casewise.NAs.all, ncol=N.vars, nrow=N.vars)
      rownames(casewise.NAs) <- cls
      colnames(casewise.NAs) <- cls

      # counts for NAs to be returned to the client:
      # This is a list with (a) a vector with the number of NAs in each variable (i.e. in each column)
      # separately and (b) the number of NAs casewise (i.e. the number of rows deleted from the input dataframe
      # which are the rows that at least one of their cells include a missing value)
      na.counts <- list(column.NAs, casewise.NAs)
      names(na.counts) <- list(paste0("Number of NAs in each column"), paste0("Number of NAs casewise"))

      # A matrix with elements the sum of products between each two variables
      sums.of.products <- matrix(ncol=N.vars, nrow=N.vars)
      rownames(sums.of.products) <- cls
      colnames(sums.of.products) <- cls
      for(m in 1:N.vars){
        for(p in 1:N.vars){
          sums.of.products[m,p] <- sum(as.numeric(as.character(casewise.dataframe[,m]))*as.numeric(as.character(casewise.dataframe[,p])))
        }
      }

      # A matrix with elements the sum of each variable
      sums <- matrix(ncol=N.vars, nrow=N.vars)
      rownames(sums) <- cls
      colnames(sums) <- cls
      for(m in 1:N.vars){
        for(p in 1:N.vars){
          sums[m,p] <- sum(as.numeric(as.character(casewise.dataframe[,m])))
        }
      }

      # A matrix with elements the sum of squares of each variable after removing missing values casewise
      sums.of.squares <- matrix(ncol=N.vars, nrow=N.vars)
      rownames(sums.of.squares) <- cls
      colnames(sums.of.squares) <- cls
      for(m in 1:N.vars){
        for(p in 1:N.vars){
          sums.of.squares[m,p] <- sum(as.numeric(as.character(casewise.dataframe[,m]))*as.numeric(as.character(casewise.dataframe[,m])))
        }
      }

      # A natrix with elements the number of complete cases casewise
      complete.counts <- matrix(dim(casewise.dataframe)[1], ncol=N.vars, nrow=N.vars)
      rownames(complete.counts) <- cls
      colnames(complete.counts) <- cls

    }

    return(list(sums.of.products=sums.of.products, sums=sums, complete.counts=complete.counts, na.counts=na.counts, sums.of.squares=sums.of.squares))
  }

  if(!identical(length(INDEX), length(X), length(Y))){
    stop("X, Y and INDEX should be of the same length")
  }


  analysis.matrix <- data.frame(X, Y, INDEX) # changed lines

  data.complete<-stats::complete.cases(analysis.matrix)

  Ntotal<-dim(analysis.matrix)[1]
  Nmissing<-sum(!data.complete)
  Nvalid<-sum(data.complete)

  simplify<-TRUE

  analysis.matrix.no.miss<-analysis.matrix[data.complete,]
  #nv<-dim(analysis.matrix)[2]

  # X<-as.vector(analysis.matrix.no.miss[,1])
  # Y<-as.vector(analysis.matrix.no.miss[,2])
  # INDEX <- analysis.matrix.no.miss[,3]

  dataframe <- as.data.frame(cbind(X,Y))

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



  #CALCULATE GROUP SIZES AND CHECK VALID

  ansmat.count<-table(group)

  # Set filter for cell sizes that are too small
  # the minimum number of observations that are allowed (the below function gets the value from opal)
  any.invalid.cell<-(sum(ansmat.count < nfilter.tab & ansmat.count > 0 ) >= 1)


  if(any.invalid.cell){
    table.valid<-FALSE
    cell.count.warning<-paste0("At least one group has between 1 and ", nfilter.tab-1, " observations. Please change groups")
    result<-list(table.valid,Nvalid,Nmissing,Ntotal,cell.count.warning)
    names(result)<-list("Table_valid","Nvalid","Nmissing","Ntotal","Warning")
    return(result)
  }


  if(!any.invalid.cell){
    table.valid<-TRUE
    cell.count.warning<-paste0("All tables valid")
    result <- lapply(split(dataframe, group), FUN) # check if to introduce lsoa names
    names(result) <- unlist(namelist)
    result[["Nfilter.tab"]] <- nfilter.tab
    return(result)
  }


}

