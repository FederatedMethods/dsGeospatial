
groupglmDS1 <- function(formula, family, weights, offset, data, index){

    
    #############################################################
    # MODULE 1: CAPTURE DISCLOSURE SETTINGS
    #############################################################
    
    thr <- dsBase::listDisclosureSettingsDS()
    
    nfilter.tab <- as.numeric(thr$nfilter.tab)
    nfilter.glm <- as.numeric(thr$nfilter.glm)
    
    #############################################################
    # MODULE 2: PREPARE DATA AND GROUPING VARIABLE
    #############################################################
    
    if (is.null(data)) {
      stop("'data' must be specified for grouped GLM fitting.")
    }
    
    dataTable <- eval(parse(text = data), envir = parent.frame())
    
    if (!is.data.frame(dataTable)) {
      stop("'data' must resolve to a data frame.")
    }
    
    # index is a column name within dataTable
    if (!is.character(index) || length(index) != 1L) {
      stop("'index' must be a single character string.")
    }
    
    if (!index %in% names(dataTable)) {
      stop("Grouping variable not found in data.")
    }
    
    group.index <- as.factor(dataTable[[index]])
    
    group.rows <- split(
      seq_len(nrow(dataTable)),
      group.index,
      drop = TRUE
    )
    
    group.names <- names(group.rows)
    
    #############################################################
    # MODULE 3: INITIALISE GROUP RESULTS
    #############################################################
    
    group.results <- vector("list", length(group.rows))
    names(group.results) <- group.names
    
    #############################################################
    # MODULE 4: VALIDATE EACH GROUP
    #############################################################
    
    for (g in seq_along(group.rows)) {
      
      rows <- group.rows[[g]]
      
      group.data <- dataTable[rows, , drop = FALSE]
      
      errorMessage <- "No errors"
      
      y.invalid <- 0
      w.invalid <- 0
      o.invalid <- 0
      glm.saturation.invalid <- 0
      
      ###########################################################
      # CHECK MINIMUM GROUP SIZE
      ###########################################################
      
      if (nrow(group.data) < nfilter.tab) {
        
        group.results[[g]] <- list(
          valid = FALSE,
          errorMessage = "Group does not satisfy disclosure requirements"
        )
        
        next
      }
      
      ###########################################################
      # FIT INITIAL GLM (ONE ITERATION)
      ###########################################################
      
      formula2use <- stats::as.formula(formula)
      
      mod.glm.ds <- tryCatch(
        stats::glm(
          formula = formula2use,
          family = family,
          data = group.data,
          x = TRUE,
          y = TRUE,
          control = stats::glm.control(maxit = 1)
        ),
        error = function(e) NULL
      )
      
      if (is.null(mod.glm.ds)) {
        
        group.results[[g]] <- list(
          valid = FALSE,
          errorMessage = "Unable to initialise GLM for this group"
        )
        
        next
      }
      
      X.mat <- as.matrix(mod.glm.ds$x)
      y.vect <- as.vector(mod.glm.ds$y)
      
      dimX <- dim(X.mat)
      
      num.p <- dimX[2]
      num.N <- dimX[1]
      
      ###########################################################
      # FIRST DISCLOSURE TRAP: MODEL SATURATION
      ###########################################################
      
      if (num.p > nfilter.glm * num.N) {
        glm.saturation.invalid <- 1
      }
      
      ###########################################################
      # SECOND DISCLOSURE TRAP: RESPONSE VARIABLE
      ###########################################################
      
      unique.values.noNA.y <- unique(
        y.vect[stats::complete.cases(y.vect)]
      )
      
      if (length(unique.values.noNA.y) == 2) {
        
        tabvar <- table(y.vect)
        tabvar <- tabvar[tabvar >= 1]
        
        if (min(tabvar) < nfilter.tab) {
          y.invalid <- 1
        }
      }
      
      ###########################################################
      # THIRD DISCLOSURE TRAP: DESIGN MATRIX
      ###########################################################
      
      num.Xpar <- ncol(X.mat)
      
      Xpar.invalid <- rep(0, num.Xpar)
      
      for (pj in seq_len(num.Xpar)) {
        
        xcol <- X.mat[, pj]
        
        unique.values.noNA <- unique(
          xcol[stats::complete.cases(xcol)]
        )
        
        if (length(unique.values.noNA) == 2) {
          
          tabvar <- table(xcol)
          tabvar <- tabvar[tabvar >= 1]
          
          if (min(tabvar) < nfilter.tab) {
            Xpar.invalid[pj] <- 1
          }
        }
      }
      
      ###########################################################
      # FOURTH DISCLOSURE TRAP: WEIGHTS
      ###########################################################
      
      if (!is.null(weights)) {
        
        w.vect <- group.data[[weights]]
        
        if (is.null(w.vect)) {
          w.invalid <- 1
        } else {
          
          unique.values.noNA.w <- unique(
            w.vect[stats::complete.cases(w.vect)]
          )
          
          if (length(unique.values.noNA.w) == 2) {
            
            tabvar <- table(w.vect)
            tabvar <- tabvar[tabvar >= 1]
            
            if (min(tabvar) < nfilter.tab) {
              w.invalid <- 1
            }
          }
        }
      }
      
      ###########################################################
      # FIFTH DISCLOSURE TRAP: OFFSET
      ###########################################################
      
      if (!is.null(offset)) {
        
        o.vect <- group.data[[offset]]
        
        if (is.null(o.vect)) {
          o.invalid <- 1
        } else {
          
          unique.values.noNA.o <- unique(
            o.vect[stats::complete.cases(o.vect)]
          )
          
          if (length(unique.values.noNA.o) == 2) {
            
            tabvar <- table(o.vect)
            tabvar <- tabvar[tabvar >= 1]
            
            if (min(tabvar) < nfilter.tab) {
              o.invalid <- 1
            }
          }
        }
      }
      
      ###########################################################
      # FINAL VALIDITY CHECK
      ###########################################################
      
      valid <- !(
        y.invalid > 0 ||
          w.invalid > 0 ||
          o.invalid > 0 ||
          any(Xpar.invalid > 0) ||
          glm.saturation.invalid > 0
      )
      
      if (!valid) {
        errorMessage <- "STUDY DATA OR APPLIED MODEL INVALID FOR THIS GROUP"
      }
      
      group.results[[g]] <- list(
        valid = valid,
        dimX = dimX,
        coef.names = names(mod.glm.ds$coefficients),
        y.invalid = y.invalid,
        Xpar.invalid = Xpar.invalid,
        w.invalid = w.invalid,
        o.invalid = o.invalid,
        glm.saturation.invalid = glm.saturation.invalid,
        errorMessage = errorMessage
      )
    }
    
    #############################################################
    # MODULE 5: RETURN GROUP VALIDATION RESULTS
    #############################################################
    
    return(group.results)
  
}