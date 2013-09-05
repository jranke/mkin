# Simple gWidgetsWWW2 GUI for mkin
# Set the GUI title and create the parent frame {{{1
GUI_title <- "Simple Browser based GUI for kinetic evaluations using mkin"
w <- gwindow(GUI_title)
sb <- gstatusbar("Powered by gWidgetsWWW2 and Rook", cont = w)
g <- gframe(GUI_title, cont = w, use.scrollwindow = TRUE, horizontal = FALSE)

# Set default values for project data objects {{{1
project_file <- "mkin_project_1.RData"
# Observed variables {{{2
n.observed <- 2
observed.names = c("parent", paste("M", 1:(n.observed - 1), sep=""))
observed.df = data.frame(Index = 1:n.observed, 
                         Name = observed.names,
                         stringsAsFactors = FALSE)

# Project definition and file management {{{1
upload_file_handler <- function(h, ...)  # {{{2
{
  tmpfile <- normalizePath(svalue(h$obj), winslash = "/")
  try(load(tmpfile))
  project_file <<- pr.gf$filename
  svalue(wf.ge) <- project_file
  observed.gdf[,] <- observed.df
}
save_to_file_handler <- function(h, ...) # {{{2
{
   observed.df <- observed.gdf[,]
   save(observed.df, file = project_file)
   galert(paste("Saved project contents to", project_file), parent = w)
}

prg <- gexpandgroup("Project definition", cont = g)

pr.vg <- ggroup(cont = prg, horizontal = FALSE)
pr.hg <- ggroup(cont = pr.vg, horizontal = TRUE)
pr.gf <- gfile(text = "Select project file", cont = pr.hg,
               handler = upload_file_handler)
pr.vg2 <- ggroup(cont = pr.hg, horizontal = FALSE)
pr.hg2 <- ggroup(cont = pr.vg2, horizontal = TRUE)
glabel("Current project file name is", cont = pr.hg2)
change_project_file_handler = function(h, ...) {
  project_file <<- as.character(svalue(h$obj))
}
wf.ge <- gedit(project_file, cont = pr.hg2, 
               handler = change_project_file_handler)

gbutton("Save current project contents to this file", cont = pr.vg2,
        handler = save_to_file_handler)

observed.gdf <- gdf(observed.df, name = "Names of observed variables", 
                    width = 500, height = 250, cont = pr.vg)
observed.gdf$set_column_width(1, 40)

# vim: set foldmethod=marker ts=2 sw=2 expandtab:
