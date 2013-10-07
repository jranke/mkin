# Simple gWidgetsWWW2 GUI for mkin
# Set the GUI title and create the parent frame {{{1
require("mkin")
GUI_title <- "Simple Browser based GUI for kinetic evaluations using mkin"
w <- gwindow(GUI_title)
sb <- gstatusbar("Powered by gWidgetsWWW2 and Rook", cont = w)
g <- gframe(GUI_title, cont = w, use.scrollwindow = TRUE, horizontal = FALSE)

# Set default values for project data objects {{{1
project_file <- "mkin_project_1.RData"
# Studies {{{2
studies.df <- data.frame(Index = as.integer(1), 
                         Author = "FOCUS kinetics workgroup",
                         Year = "2006", 
                         Title = "FOCUS Kinetics",
                         stringsAsFactors = FALSE)

# Datasets {{{2
ds <- list()
# FOCUS 2006 datasets {{{3
for (i in 1:5) {
  ds.letter = LETTERS[i]
  ds.index <- as.character(i)
  ds.name = paste0("FOCUS_2006_", ds.letter)
  ds[[ds.index]] <- list(
    study_nr = 1,
    title = paste("FOCUS example dataset", ds.letter),
    sampling_times = unique(get(ds.name)$time),
    time_unit = "NA",
    observed = as.character(unique(get(ds.name)$name)),
    unit = "% AR",
    replicates = 1,
    data = get(ds.name)
  )
  ds[[ds.index]]$data$name <- as.character(ds[[ds.index]]$data$name)
  ds[[ds.index]]$data$override = "NA"
  ds[[ds.index]]$data$weight = 1
}
# Dataframe with datasets for selection with the gtable widget {{{2
update_ds.df <- function() { # {{{3
  ds.n <- length(ds)
  ds.df <<- data.frame(Index = 1:ds.n, 
                       Study = character(ds.n), 
                       Title = character(ds.n),
                       stringsAsFactors = FALSE)
  for (i in 1:ds.n)
  {
    ds.index <- names(ds)[[i]]
    ds.df[i, "Study"] <<- ds[[ds.index]]$study_nr
    ds.df[i, "Title"] <<- ds[[ds.index]]$title
  }
}
ds.df <- data.frame()
update_ds.df()

# Set the initial dataset number globally
ds.cur = "1"
# Project data management {{{1
upload_file_handler <- function(h, ...)  # {{{2
{
  tmpfile <- normalizePath(svalue(h$obj), winslash = "/")
  try(load(tmpfile))
  project_file <<- pr.gf$filename
  svalue(wf.ge) <- project_file
  studies.gdf[,] <- studies.df 
  ds.cur <<- "1"
  ds <<- ds
  update_ds.df()
  ds.gtable[,] <- ds.df
  update_ds_editor()
}
save_to_file_handler <- function(h, ...) # {{{2
{
   studies.df <- data.frame(studies.gdf[,], stringsAsFactors = FALSE)
   save(studies.df, ds, file = project_file)
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

# Expandable group for studies {{{1
stg <- gexpandgroup("Studies", cont = g)
update_study_selector <- function(h, ...) {
  delete(ds.e.1, ds.study.gc)
  ds.study.gc <<- gcombobox(paste("Study", studies.gdf[,1]), cont = ds.e.1) 
  svalue(ds.study.gc, index = TRUE) <- ds[[ds.cur]]$study_nr
}
studies.gdf <- gdf(studies.df, name = "Studies in the project",
                   width = 500, height = 200, cont = stg)
studies.gdf$set_column_width(1, 40)
studies.gdf$set_column_width(2, 200)
addHandlerChanged(studies.gdf, update_study_selector)

# Expandable group for datasets {{{1
dsg <- gexpandgroup("Dataset selector", cont = g, horizontal = FALSE)
 
# Dataset table with handler {{{2
ds.switcher <- function(h, ...) {
  ds.cur <<- as.character(svalue(h$obj))
  update_ds_editor()
}
ds.gtable <- gtable(ds.df, cont = dsg)
addHandlerClicked(ds.gtable, ds.switcher)
size(ds.gtable) <- list(columnWidths = c(40, 40, 200))
 
# Dataset editor {{{2
# Handler functions {{{3
copy_dataset_handler <- function(h, ...) {
  ds.old <- ds.cur
  ds.cur <<- as.character(1 + length(ds))
  svalue(ds.editor) <- paste("Dataset", ds.cur)
  ds[[ds.cur]] <<- ds[[ds.old]]
  update_ds.df()
  ds.gtable[,] <- ds.df
}
 
delete_dataset_handler <- function(h, ...) {
  ds[[ds.cur]] <<- NULL
  names(ds) <<- as.character(1:length(ds))
  ds.cur <<- names(ds)[[1]]
  update_ds.df()
  ds.gtable[,] <- ds.df
  update_ds_editor()
}
 
new_dataset_handler <- function(h, ...) {
  ds.cur <<- as.character(1 + length(ds))
  ds[[ds.cur]] <<- list(
                        study_nr = 1,
                        title = "",
                        sampling_times = c(0, 1),
                        time_unit = "NA",
                        observed = "parent",
                        unit = "NA",
                        replicates = 1,
                        data = data.frame(
                                          name = "parent",
                                          time = c(0, 1),
                                          value = c(100, NA),
                                          override = "NA",
                                          weight = 1,
                                          stringsAsFactors = FALSE
                                          )
                        )
  update_ds.df()
  ds.gtable[,] <- ds.df
  update_ds_editor()
}

empty_grid_handler <- function(h, ...) {
  obs <- strsplit(svalue(ds.e.obs), ", ")[[1]]
  sampling_times <- strsplit(svalue(ds.e.st), ", ")[[1]]
  replicates <- as.numeric(svalue(ds.e.rep))
  new.data = data.frame(
    name = rep(obs, each = replicates * length(sampling_times)),
    time = rep(sampling_times, each = replicates, times = length(obs)),
    value = "NA",
    override = "NA",
    weight = 1
  )
  ds.e.gdf[,] <- new.data
}

save_ds_changes_handler <- function(h, ...) {
  ds[[ds.cur]]$title <<- svalue(ds.title.ge)
  ds[[ds.cur]]$study_nr <<- as.numeric(gsub("Study ", "", svalue(ds.study.gc)))
  update_ds.df()
  ds.gtable[,] <- ds.df
  tmpd <- ds.e.gdf[,]
  ds[[ds.cur]]$data <<- tmpd
  ds[[ds.cur]]$sampling_times <<- sort(unique(tmpd$time))
  ds[[ds.cur]]$time_unit <<- svalue(ds.e.stu)
  ds[[ds.cur]]$observed <<- unique(tmpd$name)
  ds[[ds.cur]]$unit <<- svalue(ds.e.obu)
  ds[[ds.cur]]$replicates <<- max(aggregate(tmpd$time, 
                                            list(tmpd$time, tmpd$name), length)$x)
  update_ds_editor()
}
 

# Widget setup {{{3
ds.editor <- gframe("Dataset 1", horizontal = FALSE, cont = dsg)
# Line 1 {{{4
ds.e.1 <- ggroup(cont = ds.editor, horizontal = TRUE)
glabel("Title: ", cont = ds.e.1) 
ds.title.ge <- gedit(ds[[ds.cur]]$title, cont = ds.e.1) 
glabel(" from ", cont = ds.e.1) 
ds.study.gc <- gcombobox(paste("Study", studies.gdf[,1]), cont = ds.e.1) 

# Line 2 {{{4
ds.e.2 <- ggroup(cont = ds.editor, horizontal = TRUE)
gbutton("Copy dataset", cont = ds.e.2, handler = copy_dataset_handler)
gbutton("Delete dataset", cont = ds.e.2, handler = delete_dataset_handler)
gbutton("New dataset", cont = ds.e.2, handler = new_dataset_handler)

# Line 3 with forms {{{4
ds.e.forms <- ggroup(cont= ds.editor, horizontal = TRUE)

ds.e.3a <- gvbox(cont = ds.e.forms)
ds.e.3a.gfl <- gformlayout(cont = ds.e.3a)
ds.e.st  <- gedit(paste(ds[[ds.cur]]$sampling_times, collapse = ", "),
                  width = 50,
                  label = "Sampling times", 
                  cont = ds.e.3a.gfl)
ds.e.stu <- gedit(ds[[ds.cur]]$time_unit, 
                  width = 20,
                  label = "Unit", cont = ds.e.3a.gfl)
ds.e.rep <- gedit(ds[[ds.cur]]$replicates, 
                  width = 20,
                  label = "Replicates", cont = ds.e.3a.gfl)

ds.e.3b <- gvbox(cont = ds.e.forms)
ds.e.3b.gfl <- gformlayout(cont = ds.e.3b)
ds.e.obs <- gedit(paste(ds[[ds.cur]]$observed, collapse = ", "),
                  width = 50,
                  label = "Observed", cont = ds.e.3b.gfl)
ds.e.obu <- gedit(ds[[ds.cur]]$unit,
                  width = 20, label = "Unit", 
                  cont = ds.e.3b.gfl)
gbutton("Generate empty grid for kinetic data", cont = ds.e.3b, 
        handler = empty_grid_handler)

# Save button {{{4
gbutton("Save changes", cont = ds.editor, handler = save_ds_changes_handler)

# Kinetic Data {{{4
ds.e.gdf <- gdf(ds[[ds.cur]]$data, name = "Kinetic data", 
                width = 700, height = 700, cont = ds.editor)
ds.e.gdf$set_column_width(2, 70)

# Update the editor {{{3
update_ds_editor <- function() {
  svalue(ds.editor) <- paste("Dataset", ds.cur)
  svalue(ds.title.ge) <- ds[[ds.cur]]$title
  svalue(ds.study.gc, index = TRUE) <- ds[[ds.cur]]$study_nr

  svalue(ds.e.st) <- paste(ds[[ds.cur]]$sampling_times, collapse = ", ")
  svalue(ds.e.stu) <- ds[[ds.cur]]$time_unit
  svalue(ds.e.obs) <- paste(ds[[ds.cur]]$observed, collapse = ", ")
  svalue(ds.e.obu) <- ds[[ds.cur]]$unit
  svalue(ds.e.rep) <- ds[[ds.cur]]$replicates

  ds.e.gdf[,] <- ds[[ds.cur]]$data
}
# 3}}}
# 2}}}
# 1}}}
# vim: set foldmethod=marker ts=2 sw=2 expandtab:
