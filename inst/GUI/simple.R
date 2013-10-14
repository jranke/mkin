# $Id$ {{{1

# Simple gWidgetsWWW2 GUI for mkin

# Copyright (C) 2013 Johannes Ranke
# Contact: jranke@uni-bremen.de, johannesranke@eurofins.com

# This file is part of the R package mkin

# mkin is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.

# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.

# You should have received a copy of the GNU General Public License along with
# this program. If not, see <http://www.gnu.org/licenses/>
require(mkin); require(canvas) # {{{1
# Set the GUI title and create the parent frame {{{1
GUI_title <- "Simple Browser based GUI for kinetic evaluations using mkin"
w <- gwindow(GUI_title)
sb <- gstatusbar("Powered by gWidgetsWWW2 and Rook", cont = w)
g <- gframe(GUI_title, cont = w, use.scrollwindow = TRUE, horizontal = FALSE)
# Set default values for project data {{{1
# Initial project file name {{{2
project_file <- "mkin_FOCUS_2006.RData"
# Initial studies {{{2
studies.df <- data.frame(Index = as.integer(1), 
                         Author = "FOCUS kinetics workgroup",
                         Year = "2006", 
                         Title = "FOCUS Kinetics",
                         stringsAsFactors = FALSE)

# Initial datasets {{{2
ds <- list()
observed.all <- vector()
for (i in 1:5) {
  ds.letter = LETTERS[i]
  ds.index <- as.character(i)
  ds.name = paste0("FOCUS_2006_", ds.letter)
  ds[[ds.index]] <- list(
    study_nr = 1,
    title = paste("FOCUS example dataset", ds.letter),
    sampling_times = unique(get(ds.name)$time),
    time_unit = "",
    observed = as.character(unique(get(ds.name)$name)),
    unit = "% AR",
    replicates = 1,
    data = get(ds.name)
  )
  ds[[ds.index]]$data$name <- as.character(ds[[ds.index]]$data$name)
  ds[[ds.index]]$data$override = as.numeric(NA)
  ds[[ds.index]]$data$weight = 1
}
# Initial models {{{2
m <- list()
m[["1"]] <- mkinmod(parent = list(type = "SFO"))
m[["1"]]$name = "SFO"
m[["2"]] <- mkinmod(parent = list(type = "FOMC"))
m[["2"]]$name = "FOMC"
m[["3"]] <- mkinmod(parent = list(type = "DFOP"))
m[["3"]]$name = "DFOP"
m[["4"]] <- mkinmod(parent = list(type = "SFO", to = "m1"),
                          m1 = list(type = "SFO"),
                          use_of_ff = "max")
m[["4"]]$name = "SFO_SFO"
# Initial fits {{{2
f <- list()
override <- function(d) {
  data.frame(name = d$name, time = d$time, 
             value = ifelse(is.na(d$override), d$value, d$override),
             weight = d$weight)
}
f[["1"]] <- mkinfit(m[["1"]], override(ds[["1"]]$data), err = "weight")
f[["1"]]$dataset_title = ds[["1"]]$title
f[["1"]]$model_name = m[["1"]]$name
# Data frames with datasets, models and fits to be continuosly updated {{{1
# Dataframe with datasets for selection {{{2
update_ds.df <- function() {
  ds.n <- length(ds)
  ds.df <<- data.frame(Index = 1:ds.n, 
                       Title = character(ds.n),
                       Study = character(ds.n), 
                       stringsAsFactors = FALSE)
  for (i in 1:ds.n)
  {
    ds.index <- names(ds)[[i]]
    ds.df[i, "Title"] <<- ds[[ds.index]]$title
    ds.df[i, "Study"] <<- ds[[ds.index]]$study_nr
    observed = as.character(unique(ds[[ds.index]]$data$name))
    observed.all <<- union(observed, observed.all)
  }
}
ds.df <- data.frame()
update_ds.df()
ds.cur = "1"
# Dataframe with models for selection {{{2
update_m.df <- function() {
  m.n <- length(m)
  m.df <<- data.frame(Index = 1:m.n, 
                      Name = character(m.n),
                      stringsAsFactors = FALSE)
  for (i in 1:m.n) {
    m.index <- names(m)[[i]]
    m.df[i, "Name"] <<- m[[m.index]]$name
  }
}
m.df <- data.frame()
update_m.df()
m.cur = "1"
# Dataframe with fits for selection {{{2
update_f.df <- function() {
  f.n <- length(f)
  f.df <<- data.frame(Index = 1:f.n,
                      Dataset = character(f.n),
                      Model = character(f.n),
                      stringsAsFactors = FALSE)
  for (i in 1:f.n) {
    f.index <- names(f)[[i]]
    f.df[i, "Dataset"] <<- f[[f.index]]$dataset_title
    f.df[i, "Model"] <<- f[[f.index]]$model_name
  }
}
f.df <- data.frame()
update_f.df()
f.cur = "1"
# Expandable group for project data management {{{1
prg <- gexpandgroup("Project file management", cont = g)
# Project data management handler functions {{{2
upload_file_handler <- function(h, ...)
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
  m.cur <<- "1"
  m <<- m
  update_m.df()
  m.gtable[,] <- m.df
  update_m_editor()
}
save_to_file_handler <- function(h, ...)
{
   studies.df <- data.frame(studies.gdf[,], stringsAsFactors = FALSE)
   save(studies.df, ds, m, file = project_file)
   galert(paste("Saved project contents to", project_file), parent = w)
}
# Project data management GUI elements {{{2
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

# Datasets and models {{{1
dsm <- gframe("Datasets and models - double click to edit", cont = g,
              horizontal = TRUE)
 
# Dataset table with handler {{{2
ds.switcher <- function(h, ...) {
  ds.cur <<- as.character(svalue(h$obj))
  update_ds_editor()
  visible(dse) <- TRUE
  visible(me) <- FALSE
}
ds.gtable <- gtable(ds.df, multiple = TRUE, cont = dsm)
addHandlerDoubleClick(ds.gtable, ds.switcher)
size(ds.gtable) <- list(columnWidths = c(40, 200, 40))

# Model table with handler {{{2
m.switcher <- function(h, ...) {
  m.cur <<- as.character(svalue(h$obj))
  update_m_editor()
  visible(dse) <- FALSE
  visible(me) <- TRUE
}
m.gtable <- gtable(m.df, multiple = TRUE, cont = dsm)
addHandlerDoubleClick(m.gtable, m.switcher)
size(m.gtable) <- list(columnWidths = c(40, 200))

# Section for selecting datasets and model
dsmsel <- gvbox(cont = dsm)
ds_plot_handler <- function(h, ...) {
  ds.sel <- svalue(ds.gtable)
  n.ds.sel <- length(ds.sel)
  for (i in 1:n.ds.sel) {
    prows[[i]] <- ggroup(cont = pfv)
    d <- ds[[ds.sel[[i]]]]

    f <- tempfile()
    canvas(f, width = 500, height = 350)
      plot(0, type = "n",
           xlim = c(0, max(d$data$time, na.rm = TRUE)),
           ylim = c(0, max(d$data$value, na.rm = TRUE)),
           main = d$title)
      for (obs_var in d$observed) {
        points(subset(d$data, name == obs_var, c(time, value))) 
      }
    dev.off()
    plots[[i]] <- gcanvas(f, cont = prows[[i]], width = 500, 350)
  }

}
dsplot <- gbutton("Plot selected datasets", cont = dsmsel, 
                  handler = ds_plot_handler)
 
# Expandable group for the dataset editor {{{1
dse <- gexpandgroup("Dataset editor", cont = g, horizontal = FALSE)
visible(dse) <- FALSE

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
    value = NA,
    override = NA,
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
ds.editor <- gframe("Dataset 1", horizontal = FALSE, cont = dse)
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
enter_next_value_handler <- function(h, ...) galert("next value", parent = w)
addHandlerChanged(ds.e.gdf, enter_next_value_handler)

# Update the dataset editor {{{3
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

# Expandable group for the model editor {{{1
me <- gexpandgroup("Model editor", cont = g, horizontal = FALSE)
visible(me) <- FALSE

# Handler functions {{{3
copy_model_handler <- function(h, ...) {
  m.old <- m.cur
  m.cur <<- as.character(1 + length(m))
  svalue(m.editor) <- paste("Model", m.cur)
  m[[m.cur]] <<- m[[m.old]]
  update_m.df()
  m.gtable[,] <- m.df
}
 
delete_model_handler <- function(h, ...) {
  m[[m.cur]] <<- NULL
  names(m) <<- as.character(1:length(m))
  m.cur <<- "1"
  update_m.df()
  m.gtable[,] <- m.df
  update_m_editor()
}

add_observed_handler <- function(h, ...) {
  obs.i <- length(m.e.rows) + 1
  m.e.rows[[obs.i]] <<- ggroup(cont = m.editor, horizontal = TRUE)
  m.e.obs[[obs.i]] <<- gcombobox(observed.all, selected = obs.i, 
                                cont = m.e.rows[[obs.i]])
  m.e.type[[obs.i]] <<- gcombobox(c("SFO", "FOMC", "DFOP", "HS", "SFORB"),
                                 cont = m.e.rows[[obs.i]])
  svalue(m.e.type[[obs.i]]) <- "SFO"
  glabel("to", cont = m.e.rows[[obs.i]]) 
  m.e.to[[obs.i]] <<- gedit("", cont = m.e.rows[[obs.i]])
  m.e.sink[[obs.i]] <<- gcheckbox("Path to sink", 
                                  checked = TRUE, cont = m.e.rows[[obs.i]]) 
  gbutton("Remove compound", handler = remove_compound_handler, 
          action = obs.i, cont = m.e.rows[[obs.i]])
}

remove_compound_handler <- function(h, ...) {
  m[[m.cur]]$spec[[h$action]] <<- NULL
  update_m_editor()
}

save_m_changes_handler <- function(h, ...) {
  spec <- list()
  for (obs.i in 1:length(m.e.rows)) {
    spec[[obs.i]] <- list(type = svalue(m.e.type[[obs.i]]),
                          to = svalue(m.e.to[[obs.i]]),
                          sink = svalue(m.e.sink[[obs.i]]))
    if(spec[[obs.i]]$to == "") spec[[obs.i]]$to = NULL
    names(spec)[[obs.i]] <- svalue(m.e.obs[[obs.i]])
  }
  m[[m.cur]] <<- mkinmod(use_of_ff = svalue(m.ff.gc), 
                         speclist = spec)
  m[[m.cur]]$name <<- svalue(m.name.ge) 
  update_m.df()
  m.gtable[,] <- m.df
}
 
# Widget setup {{{3
m.editor <- gframe("Model 1", horizontal = FALSE, cont = me)
m.e.0 <- ggroup(cont = m.editor, horizontal = TRUE)
glabel("Model name: ", cont = m.e.0) 
m.name.ge <- gedit(m[[m.cur]]$name, cont = m.e.0) 
glabel("Use of formation fractions: ", cont = m.e.0) 
m.ff.gc <- gcombobox(c("min", "max"), cont = m.e.0)
svalue(m.ff.gc) <- m[[m.cur]]$use_of_ff

# Model handling buttons {{{4
m.e.b <- ggroup(cont = m.editor, horizontal = TRUE)
gbutton("Copy model", cont = m.e.b, handler = copy_model_handler)
gbutton("Delete model", cont = m.e.b, handler = delete_model_handler)
gbutton("Add transformation product", cont = m.e.b, 
        handler = add_observed_handler)
gbutton("Save changes", cont = m.e.b, handler = save_m_changes_handler)


m.observed <- names(m[[m.cur]]$spec)
m.e.rows <- m.e.obs <- m.e.type <- m.e.to <- m.e.sink <- list()
obs.to <- ""

show_m_spec <- function() {
  for (obs.i in 1:length(m.observed)) {
    m.e.rows[[obs.i]] <<- ggroup(cont = m.editor, horizontal = TRUE)
    m.e.obs[[obs.i]] <<- gcombobox(observed.all, selected = obs.i, 
                                  cont = m.e.rows[[obs.i]])
    m.e.type[[obs.i]] <<- gcombobox(c("SFO", "FOMC", "DFOP", "HS", "SFORB"),
                                   cont = m.e.rows[[obs.i]])
    svalue(m.e.type[[obs.i]]) <<- m[[m.cur]]$spec[[obs.i]]$type
    glabel("to", cont = m.e.rows[[obs.i]]) 
    obs.to <<- ifelse(is.null(m[[m.cur]]$spec[[obs.i]]$to), "",
                 m[[m.cur]]$spec[[obs.i]]$to)
    m.e.to[[obs.i]] <<- gedit(obs.to, cont = m.e.rows[[obs.i]])
    m.e.sink[[obs.i]] <<- gcheckbox("Path to sink", checked = m[[m.cur]]$spec[[obs.i]]$sink,
              cont = m.e.rows[[obs.i]]) 
    if (obs.i > 1) {
      gbutton("Remove compound", handler = remove_compound_handler, 
              action = obs.i, cont = m.e.rows[[obs.i]])
    }
  }
}
show_m_spec()

# Update the model editor {{{3
update_m_editor <- function() {
  svalue(m.editor) <- paste("Model", m.cur)
  svalue(m.name.ge) <- m[[m.cur]]$name
  svalue(m.ff.gc) <- m[[m.cur]]$use_of_ff
  for (oldrow.i in 1:length(m.e.rows)) {
    delete(m.editor, m.e.rows[[oldrow.i]])
  }
  m.observed <<- names(m[[m.cur]]$spec)
  m.e.rows <<- m.e.obs <<- m.e.type <<- m.e.to <<- m.e.sink <<- list()
  show_m_spec()
}

# 3}}}
# 2}}}
# Plots and fitting {{{1
pf <- gframe("Plots and fitting", cont = g)
pfv <- gvbox(cont = pf)
prows <- plots <- list()

# Show the fits {{{1
#mf <- gnotebook(cont = g)
# s <- s.gt <- list()
#s[[1]] <- summary(fits[[1]])
#for (i in 1:length(fits)) {
#  fits[[i]] <- gframe(fits[[1]]$name, cont = mf, label = i)
#  s.tmp <- capture.output(print(s[[i]]))
#  s.gt[[i]] <- gtext(s.tmp, width = 600, cont = fits[[i]],
#                     use.codemirror = TRUE)
#}
# 1}}}
# vim: set foldmethod=marker ts=2 sw=2 expandtab:
