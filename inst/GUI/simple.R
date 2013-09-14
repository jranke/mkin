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
                         stringsAsFactors = FALSE)

# Datasets {{{2
ds <- list()
# FOCUS 2006 datasets {{{3
for (i in 1:6) {
  ds.letter = LETTERS[i]
  ds.name = paste0("FOCUS_2006_", ds.letter)
  ds[[i]] <- list(
    study_nr = 1,
    dataset_nr = i,
    title = paste("FOCUS example dataset", ds.letter),
    sampling_times = unique(get(ds.name)$time),
    time_unit = "NA",
    observed = length(levels(get(ds.name)$name)),
    unit = "% AR",
    replicates = 1,
    data = get(ds.name)
  )
  ds[[i]]$data$override = "NA"
  ds[[i]]$data$weight = 1
}
# Dataframe with datasets for selecting them with the gtable widget {{{2
ds.df <- data.frame()
update_ds.df <- function() {
  ds.n <- length(ds)
  ds.df <<- data.frame(Index = 1:ds.n, Study = character(ds.n), Title = character(ds.n), 
                       icon = asIcon(rep("editor", ds.n)), stringsAsFactors = FALSE)
  for (i in 1:ds.n)
  {
    ds.df[i, "Study"] <<- ds[[i]]$study_nr
    ds.df[i, "Title"] <<- ds[[i]]$title
  }
}
update_ds.df()

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
#  ds <- ds
}
save_to_file_handler <- function(h, ...) # {{{2
{
   observed.df <- data.frame(observed.gdf[,], stringsAsFactors = FALSE)
   studies.df <- data.frame(studies.gdf[,], stringsAsFactors = FALSE)
   save(observed.df, studies.df, ds, file = project_file)
   galert(paste("Saved project contents to", project_file), parent = w)
}

# Add widgets for project file management to an expandable group {{{1
prg <- gexpandgroup("Project file management", cont = g)

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
                   width = 500, height = 200, cont = stg)
studies.gdf$set_column_width(1, 40)
studies.gdf$set_column_width(2, 200)

# Expandable group for datasets {{{1
dsg <- gexpandgroup("Datasets", cont = g, horizontal = FALSE)
ds.handler <- function(h, ...) {
  ds.cur <<- svalue(h$obj, index = FALSE)
  delete(ds.editor, ds.e.head)
  delete(ds.editor, ds.e.data)
  show_ds_editor()
}
ds.gtable <- gtable(ds.df, handler = ds.handler, cont = dsg)
size(ds.gtable) <- list(columnWidths = c(40, 40, 200, 40))

# Dataset editor {{{2
update_dataset_handler <- function(h, ...) {
  galert("test", parent = w)
}

new_dataset_handler <- function(h, ...) {
  galert("test", parent = w)
}

ds.editor <- gframe("Dataset 1", horizontal = FALSE, cont = dsg)
ds.e.head <- ggroup(cont = ds.editor, horizontal = FALSE)
ds.e.data <- ggroup(cont = ds.editor, horizontal = FALSE)

show_ds_editor <- function() {
  svalue(ds.editor) <- paste("Dataset", ds.cur)
  ds.e.head <<- ggroup(cont = ds.editor, horizontal = FALSE)
  ds.e.1 <- ggroup(cont = ds.e.head, horizontal = TRUE)
  glabel("Title: ", cont = ds.e.1) 
  ds.title.ge <- gedit(ds[[ds.cur]]$title, cont = ds.e.1, 
    handler = function(h, ...) { 
      ds[[ds.cur]]$title <<- svalue(h$obj)
      update_ds.df()
      ds.gtable[,] <- ds.df
    })
  glabel(" from ", cont = ds.e.1) 
  ds.study.gc <- gcombobox(paste("Study", studies.gdf[,1]), cont = ds.e.1) 
  ds.e.save <- gbutton("Save changes", cont = ds.e.1, 
                            handler = update_dataset_handler)
  ds.e.new <- gbutton("New dataset", cont = ds.e.1, 
                            handler = new_dataset_handler)

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
  ds.e.data <<- ggroup(cont = ds.editor, horizontal = FALSE)
  ds.e.gdf <- gdf(ds[[ds.cur]]$data, name = "Kinetic data", 
                   width = 700, height = 700, cont = ds.e.data)
  ds.e.gdf$set_column_width(2, 50)
  ds.e.gdf$set_column_width(3, 50)
  ds.e.gdf$set_column_width(4, 50)
}
show_ds_editor()
# 1}}}
# vim: set foldmethod=marker ts=2 sw=2 expandtab:
