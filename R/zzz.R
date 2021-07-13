
#' @import digest
# Auxiliary for getSplash.R so we can use the original file and don't have to change anything there

.onLoad <- function(libname, pkgname) {
  RMassBank.env <<- new.env()
  RMassBank.env$ReadAnnotation <- FALSE
  RMassBank.env$testnumber <- 1
  ## new variables
  RMassBank.env$verbose.output <- FALSE
  RMassBank.env$logging.file <- NULL
  RMassBank.env$export.invalid <- FALSE
  RMassBank.env$export.molfiles <- TRUE
  RMassBank.env$strictMsMsSpectraSelection <- FALSE
  
  mb <- list()
  attach(RMassBank.env)
}

utils::globalVariables(c("cpdID", "isotopes","mzCalc"))


