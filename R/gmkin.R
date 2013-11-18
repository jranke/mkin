gmkin <- function() {
  if (require(gWidgetsWWW2)) load_app(system.file("GUI/gmkin.R", package = "mkin"))
  else {
    message(
       "\nYou need to install gWidgetsWWW2 in order to run the mkin GUI gmkin.\n",
       "This  package is not currently on CRAN but can be installed from github\n",
       "using the devtools package:\n",
       "library(devtools)\n",
       "install_github('gWidgetsWWW2', 'jverzani')")
  }
}
