# Work in progress - simple Gtk2 GUI for mkin
require(gWidgets)
options("guiToolkit"="RGtk2")

w <- gwindow("Simple R GUI for kinetic evaluations",
  width=800, height=500)

# Project definition expanding group
pr <- gexpandgroup("Project definition", container=w)

n.observed <- 1
max.n.observed <- 9
observed.names = c("parent", paste("M", 1:(max.n.observed - 1), sep=""))

prg <- ggroup(horizontal=FALSE, cont = w)
prl <- glayout(cont = prg)
prl[1,1] <- glabel("Number of observed variables", cont=prl)
prl[1,2] <- (n.observed.gw = gcombobox(
  1:max.n.observed, 
  handler = function(h, ...) {
    n.observed <- svalue(n.observed.gw)
    visible(observed.gw) <- c(
      rep(TRUE, svalue(n.observed)), 
      rep(FALSE, max.n.observed - n.observed))
  },
  cont=prl))

observed.gw <- gdf(
  items = data.frame(Index = 1:max.n.observed, 
    Name = observed.names, stringsAsFactors=FALSE),
  name = "Names of observed variables",
  cont=prg)

