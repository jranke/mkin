# $Id: doRUnit.R 96 2011-04-29 11:10:40Z jranke $
# Adapted from a version around 2.9 of the rcdk package by Rajarshi Guha
if(require("RUnit", quietly=TRUE)) {
 
  pkg <- "mkin"
  path <- system.file(package=pkg, "unitTests")

  cat("\nRunning unit tests\n")
  print(list(pkg=pkg, getwd=getwd(), pathToUnitTests=path))
 
  library(package=pkg, character.only=TRUE)
 
  ## Define tests
  testSuite <- defineTestSuite(name=paste(pkg, " Unit Tests"),
                                          dirs=path)
  ## Run
  tests <- runTestSuite(testSuite)
 
  ## Default report name
  pathReport <- file.path(path, "report")
 
  ## Report to stdout and text files
  cat("------------------- UNIT TEST SUMMARY ---------------------\n\n")
  printTextProtocol(tests, showDetails=FALSE)
  printTextProtocol(tests, showDetails=FALSE,
                    fileName=paste(pathReport, "Summary.txt", sep=""))
  printTextProtocol(tests, showDetails=TRUE,
                    fileName=paste(pathReport, ".txt", sep=""))
 
  ## Report to HTML file
  printHTMLProtocol(tests, fileName=paste(pathReport, ".html", sep=""))
 
  ## Return stop() to cause R CMD check stop in case of
  ##  - failures i.e. FALSE to unit tests or
  ##  - errors i.e. R errors
  tmp <- getErrors(tests)
  if(tmp$nFail > 0 | tmp$nErr > 0) {
    stop(paste("\n\nunit testing failed (#test failures: ", tmp$nFail,
               ", #R errors: ",  tmp$nErr, ")\n\n", sep=""))
  }
} else {
  warning("cannot run unit tests -- package RUnit is not available")
}
