# Simple gWidgetsWWW2 GUI for mkin
# Set the GUI title and create the parent frame {{{1
require("mkin")
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
                         Chemical = "NA",
                         stringsAsFactors = FALSE)
# Studies {{{2
studies.df <- data.frame(Index = as.integer(1), 
                         Author = "FOCUS kinetics workgroup",
                         Year = "2006", Title = "FOCUS Kinetics",
                         Datasets = as.integer(3), 
                         stringsAsFactors = FALSE)

# Datasets {{{2
ds <- list()
# FOCUS 2006 A {{{3
ds[[1]] <- list(
  study_nr = 1,
  dataset_nr = 1,
  title = "FOCUS example dataset A",
  sampling_times = unique(FOCUS_2006_A$time),
  time_unit = "NA",
  observed = "1",
  unit = "% AR",
  replicates = 1,
  data = FOCUS_2006_A
)
ds[[1]]$data$override = "NA"
ds[[1]]$data$weight = 1

# Set the initial dataset number
ds.cur = 1
# Project data management {{{1
upload_file_handler <- function(h, ...)  # {{{2
{
  tmpfile <- normalizePath(svalue(h$obj), winslash = "/")
  try(load(tmpfile))
  project_file <<- pr.gf$filename
  svalue(wf.ge) <- project_file
  observed.gdf[,] <- observed.df
  studies.gdf[,] <- studies.df 
}
save_to_file_handler <- function(h, ...) # {{{2
{
   observed.df <- observed.gdf[,]
   studies.df <- studies.gdf[,]
   save(observed.df, studies.df, file = project_file)
   galert(paste("Saved project contents to", project_file), parent = w)
}

# Add widgets for project data management to an expandable group {{{1
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

# Expandable group for observed variables {{{1
ovg <- gexpandgroup("Observed variables", cont = g)
observed.gdf <- gdf(observed.df, name = "Names of observed variables", 
                    width = 500, height = 250, cont = ovg)
observed.gdf$set_column_width(1, 40)

# Expandable group for studies {{{1
stg <- gexpandgroup("Studies", cont = g)
studies.gdf <- gdf(studies.df, name = "Studies in the project", 
                    width = 500, height = 250, cont = stg)
studies.gdf$set_column_width(1, 40)
studies.gdf$set_column_width(2, 200)

# Expandable group for datasets {{{1
dsg <- gexpandgroup("Datasets", cont = g)
# The following will be generated from the dataset object list later
ds.df <- data.frame(Index = 1:3,
                   Study = as.integer(1),
                   Title = paste("FOCUS dataset", c("A", "B", "C")), 
                   icon = asIcon(c("editor", "editor", "editor")),
                   stringsAsFactors = FALSE)
ds.gtable <- gtable(ds.df, width = "auto", cont = dsg)
ds.gtable$set_column_width(1, 40)
ds.gtable$set_column_width(2, 40)
ds.gtable$set_column_width(3, 200)

# Dataset editor {{{2
update_dataset_handler <- function(h, ...) {
}
ds.editor <- gframe(paste("Dataset", ds.cur), cont = g)
ds.e.head <- ggroup(cont = ds.editor, horizontal = FALSE)
ds.e.1 <- ggroup(cont = ds.e.head, horizontal = TRUE)
glabel(paste0("Study ", ds[[ds.cur]]$study_nr, ": "), cont = ds.e.1) 
ds.title.ge <- gedit(ds[[ds.cur]]$title, cont = ds.e.1, 
                     handler = update_dataset_handler)

ds.e.2 <- glayout(cont = ds.e.head)
ds.e.2[1, 1] <- glabel("Sampling times: ", cont = ds.e.2) 
ds.e.2[1, 2] <- gedit(paste(ds[[ds.cur]]$sampling_times, collapse = ", "), 
                      cont = ds.e.2)
ds.e.2[1, 3] <- glabel("Unit: ", cont = ds.e.2)
ds.e.2[1, 4] <- gedit(ds[[ds.cur]]$time_unit, width = 8, cont = ds.e.2)
ds.e.2[2, 1] <- glabel("Observed variables: ", cont = ds.e.2) 
ds.e.2[2, 2] <- gedit(ds[[ds.cur]]$observed, cont = ds.e.2) 
ds.e.2[2, 3] <- glabel("Unit: ", cont = ds.e.2)
ds.e.2[2, 4] <- gedit(ds[[ds.cur]]$unit, width = 8, cont = ds.e.2)
ds.e.2[3, 1] <- glabel("Replicates: ", cont = ds.e.2) 
ds.e.2[3, 2] <- gedit(ds[[ds.cur]]$replicates, width = 2, cont = ds.e.2)
ds.e.2[3, 3:4] <- gbutton("Generate empty grid", cont = ds.e.2, 
                          handler = update_dataset_handler)
visible(ds.e.2) <- TRUE
# 1}}}
# vim: set foldmethod=marker ts=2 sw=2 expandtab:
