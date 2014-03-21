# gWidgetsWWW2 GUI for mkin {{{1

# Copyright (C) 2013,2014 Johannes Ranke
# Portions of this file are copyright (C) 2013 Eurofins Regulatory AG, Switzerland
# Contact: jranke@uni-bremen.de

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

require(mkin) # {{{1
# Set the GUI title and create the basic widget layout {{{1
w      <- gwindow("gmkin - Browser based GUI for kinetic evaluations using mkin")
sb     <- gstatusbar(paste("Powered by gWidgetsWWW2, ExtJS, Rook, FME, deSolve and minpack.lm --- Working directory is", getwd()), cont = w)
pg     <- gpanedgroup(cont = w, default.size = 300)
center <- gnotebook(cont = pg)
left   <- gvbox(cont = pg)
# Helper functions {{{1
# Override function for making it possible to override original data in the GUI {{{2
override <- function(d) {
  data.frame(name = d$name, time = d$time, 
             value = ifelse(is.na(d$override), d$value, d$override),
             err = d$err)
}
# Set default values for project data {{{1
# Initial project file name {{{2
project_name <- "mkin_FOCUS_2006"
project_file <- paste0(project_name, ".RData")
# Initial studies {{{2
studies.df <- data.frame(Index = as.integer(1), 
                         Citation = "FOCUS (2006) Guidance on degradation kinetics",
                         stringsAsFactors = FALSE)
# Initial datasets {{{2
ds <- list()
observed.all <- vector()
for (i in 1:2) {
  ds.letter = LETTERS[i + 2]
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
  ds[[ds.index]]$data$err = 1
}
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
# Initial fit lists {{{2
# The fits and summaries are collected in lists of lists
f <- s <- list()
# Dataframe with fits for selection {{{2
update_f.df <- function() {
  f.df <<- data.frame(Fit = character(),
                      Dataset = character(),
                      Model = character(), 
                      stringsAsFactors = FALSE)
  f.count <- 0
  for (fit.index in names(f)) {
    f.count <- f.count + 1
    ftmp <- f[[fit.index]]
    f.df[f.count, ] <<- c(as.character(f.count), ftmp$ds.index, ftmp$mkinmod$name)
  }
}
f.df.empty <- f.df <- data.frame(Fit = "0", 
                               Dataset = "", 
                               Model = "",
                               stringsAsFactors = FALSE)
# Widgets and handlers for project data {{{1
prg <- gexpandgroup("Project file management", cont = left, horizontal = FALSE)
# Project data management handler functions {{{2
upload_file_handler <- function(h, ...)
{
  # General
  tmpfile <- normalizePath(svalue(h$obj), winslash = "/")
  project_file <<- pr.gf$filename
  project_name <<- try(load(tmpfile))
  if (inherits(project_name, "try-error")) {
    galert(paste("Failed to load", project_file, "from", getwd()), parent = w)
  } 

  svalue(sb) <- paste("Loaded project file", project_file, "from working directory", getwd())
  svalue(pr.ge) <- project_name
  workspace <- get(project_name)

  # Studies
  studies.gdf[,] <- studies.df <- workspace$studies.df

  # Datasets
  ds.cur <<- workspace$ds.cur
  ds <<- workspace$ds
  update_ds.df()
  ds.gtable[,] <- ds.df
  update_ds_editor()

  # Models
  m.cur <<- workspace$m.cur
  m <<- workspace$m
  update_m.df()
  m.gtable[,] <- m.df
  update_m_editor()

  # Fits
  f.cur <<- workspace$f.cur
  f <<- workspace$f
  s <<- workspace$s
  if (length(f) > 0) {
    update_f.df()
    ftmp <<- f[[f.cur]]
    stmp <<- s[[f.cur]]
    ds.i <<- ds.cur
    update_plotting_and_fitting()
  } else {
    f.df <<- f.df.empty
    update_ds_editor()
    svalue(center) <- 1
  }
  f.gtable[,] <- f.df
}
save_to_file_handler <- function(h, ...)
{
  studies.df <- data.frame(studies.gdf[,], stringsAsFactors = FALSE)
  workspace <- list(
                    studies.df = studies.df,

                    ds = ds,
                    ds.cur = ds.cur,

                    m = m,
                    m.cur = m.cur,

                    f = f,
                    f.cur = f.cur,

                    s = s)
  assign(project_name, workspace)
  save(list = project_name, file = project_file)
  svalue(sb) <- paste("Saved project contents to", project_file, "in working directory", getwd())
}
change_project_name_handler = function(h, ...) {
  project_name <<- as.character(svalue(h$obj))
  project_file <<- paste0(project_name, ".RData")
}
# Project data management GUI elements {{{2
pr.gf <- gfile(text = "Select project file", cont = prg,
               handler = upload_file_handler)
pr.ge <- gedit(project_name, cont = prg, label = "Project", width = 240,
               handler = change_project_name_handler)
# The save button is always visible {{{1
gbutton("Save current project contents", cont = left,
        handler = save_to_file_handler)

# GUI widgets and a function for Studies {{{1
stg <- gexpandgroup("Studies", cont = left)
visible(stg) <- FALSE
update_study_selector <- function(h, ...) {
  delete(ds.e.1, ds.study.gc)
  ds.study.gc <<- gcombobox(paste("Study", studies.gdf[,1]), cont = ds.e.1) 
  svalue(ds.study.gc, index = TRUE) <- ds[[ds.cur]]$study_nr
}
studies.gdf <- gdf(studies.df, name = "Edit studies in the project",
                   width = 290, height = 200, cont = stg)
studies.gdf$set_column_width(1, 40)
studies.gdf$set_column_width(2, 240)
addHandlerChanged(studies.gdf, update_study_selector)
# Datasets and models {{{1
dsm <- gframe("Datasets and models", cont = left, horizontal = FALSE)
# Dataset table with handler {{{2
ds.switcher <- function(h, ...) {
  ds.cur <<- as.character(svalue(h$obj))
  update_ds_editor()
  svalue(center) <- 1
}
ds.gtable <- gtable(ds.df, width = 290, cont = dsm)
addHandlerDoubleClick(ds.gtable, ds.switcher)
size(ds.gtable) <- list(columnWidths = c(40, 200, 40))
ds.gtable$value <- 1

# Model table with handler {{{2
m.switcher <- function(h, ...) {
  m.cur <<- as.character(svalue(h$obj))
  update_m_editor()
  svalue(center) <- 2
}
m.gtable <- gtable(m.df, width = 290, cont = dsm)
addHandlerDoubleClick(m.gtable, m.switcher)
size(m.gtable) <- list(columnWidths = c(40, 240))
m.gtable$value <- 1

# Button for setting up a fit for the selected dataset and model
gbutton("Configure fit for selected model and dataset", cont = dsm, 
        handler = function(h, ...) {
          ds.i <<- as.character(svalue(ds.gtable))
          m.i <<- as.character(svalue(m.gtable))
          ftmp <<- suppressWarnings(mkinfit(m[[m.i]],
                                            override(ds[[ds.i]]$data),
                                            err = "err", 
                                            control.modFit = list(maxiter = 0)))
          ftmp$ds.index <<- ds.i
          ftmp$ds <<- ds[[ds.i]]
          stmp <<- summary(ftmp)
          svalue(pf) <- paste0("Dataset ", ds.i, ", Model ", m[[m.i]]$name)
          show_plot("Initial", default = TRUE)
          svalue(f.gg.opts.st) <<- ftmp$solution_type
          svalue(f.gg.opts.weight) <<- ftmp$weight
          svalue(f.gg.opts.atol) <<- ftmp$atol
          svalue(f.gg.opts.rtol) <<- ftmp$rtol
          svalue(f.gg.opts.reweight.method) <<- ifelse(
                                         is.null(ftmp$reweight.method),
                                         "none", ftmp$reweight.method)
          svalue(f.gg.opts.reweight.tol) <<- ftmp$reweight.tol
          svalue(f.gg.opts.reweight.max.iter) <<- ftmp$reweight.max.iter
          f.gg.parms[,] <- get_Parameters(stmp, FALSE)
          svalue(f.gg.summary) <- capture.output(stmp)
          svalue(center) <- 3
        })

# Fits {{{1
f.gf <- gframe("Fits", cont = left, horizontal = FALSE)
# Fit table with handler {{{2
f.switcher <- function(h, ...) {
  if (svalue(h$obj) != "0") {
    f.cur <<- svalue(h$obj)
    ftmp <<- f[[f.cur]]
    stmp <<- s[[f.cur]]
    ds.i <<- ftmp$ds.index
    update_plotting_and_fitting()
  }
  svalue(center) <- 3
}
f.gtable <- gtable(f.df, width = 290, cont = f.gf)
addHandlerDoubleClick(f.gtable, f.switcher)
size(f.gtable) <- list(columnWidths = c(40, 60, 180))

# Dataset editor {{{1
ds.editor <- gframe("Dataset 1", horizontal = FALSE, cont = center, label = "Dataset editor")
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
                        time_unit = "",
                        observed = "parent",
                        unit = "",
                        replicates = 1,
                        data = data.frame(
                                          name = "parent",
                                          time = c(0, 1),
                                          value = c(100, NA),
                                          override = "NA",
                                          err = 1,
                                          stringsAsFactors = FALSE
                                          )
                        )
  update_ds.df()
  ds.gtable[,] <- ds.df
  update_ds_editor()
}

new_ds_from_csv_handler <- function(h, ...) {
  tmpfile <- normalizePath(svalue(h$obj), winslash = "/")
  tmpd <- try(read.table(tmpfile, sep = "\t", header = TRUE, stringsAsFactors = FALSE))
  tmpdw <- mkin_wide_to_long(tmpd)
  if (class(tmpd) != "try-error") {
    ds.cur <<- as.character(1 + length(ds))
    ds[[ds.cur]] <<- list(
                          study_nr = NA,
                          title = "New upload",
                          sampling_times = sort(unique(tmpd$t)),
                          time_unit = "",
                          observed = unique(tmpdw$name),
                          unit = "",
                          replicates = max(aggregate(tmpdw$time,
                                                       list(tmpdw$time,
                                                            tmpdw$name),
                                                     length)$x),
                          data = tmpdw)
    ds[[ds.cur]]$data$override <<- as.numeric(NA)
    ds[[ds.cur]]$data$err <<- 1
    update_ds.df()
    ds.gtable[,] <- ds.df
    update_ds_editor()
  } else {
    galert("Uploading failed", parent = "w")
  }
}

empty_grid_handler <- function(h, ...) {
  obs <- strsplit(svalue(ds.e.obs), ", ")[[1]]
  sampling_times <- strsplit(svalue(ds.e.st), ", ")[[1]]
  replicates <- as.numeric(svalue(ds.e.rep))
  new.data = data.frame(
    name = rep(obs, each = replicates * length(sampling_times)),
    time = as.numeric(rep(sampling_times, each = replicates, times = length(obs))),
    value = as.numeric(NA),
    override = as.numeric(NA),
    err = 1,
    stringsAsFactors = FALSE
  )
  ds.e.gdf[,] <- new.data
}

keep_ds_changes_handler <- function(h, ...) {
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
  observed.all <<- union(observed.all, ds[[ds.cur]]$observed)
}
 
# Widget setup {{{3
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

upload_dataset.gf <- gfile(text = "Upload text file", cont = ds.e.2, 
        handler = new_ds_from_csv_handler)

# Line 3 with forms {{{4
ds.e.forms <- ggroup(cont= ds.editor, horizontal = TRUE)

ds.e.3a <- gvbox(cont = ds.e.forms)
ds.e.3a.gfl <- gformlayout(cont = ds.e.3a)
ds.e.st  <- gedit(paste(ds[[ds.cur]]$sampling_times, collapse = ", "),
                  width = 40,
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
generate_grid.gb <- gbutton("Generate empty grid for kinetic data", cont = ds.e.3b, 
        handler = empty_grid_handler)
tooltip(generate_grid.gb) <- "Overwrites the kinetic data shown below"

# Keep button {{{4
gbutton("Keep changes", cont = ds.editor, handler = keep_ds_changes_handler)

# Kinetic Data {{{4
ds.e.gdf <- gdf(ds[[ds.cur]]$data, name = "Kinetic data", 
                width = 500, height = 700, cont = ds.editor)
ds.e.gdf$set_column_width(2, 70)

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
# Model editor {{{1
m.editor <- gframe("Model 1", horizontal = FALSE, cont = center, label = "Model editor")
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

keep_m_changes_handler <- function(h, ...) {
  spec <- list()
  for (obs.i in 1:length(m.e.rows)) {
    to_vector = strsplit(svalue(m.e.to[[obs.i]]), ", ")[[1]]
    if (length(to_vector) == 0) to_vector = ""
    spec[[obs.i]] <- list(type = svalue(m.e.type[[obs.i]]),
                          to = to_vector,
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
gbutton("Keep changes", cont = m.e.b, handler = keep_m_changes_handler)


m.observed <- names(m[[m.cur]]$spec)
m.e.rows <- m.e.obs <- m.e.type <- m.e.to <- m.e.sink <- list()
obs.to <- ""

# Show the model specification {{{4
show_m_spec <- function() {
  for (obs.i in 1:length(m[[m.cur]]$spec)) {
    obs.name <- names(m[[m.cur]]$spec)[[obs.i]]
    m.e.rows[[obs.i]] <<- ggroup(cont = m.editor, horizontal = TRUE)
    m.e.obs[[obs.i]] <<- gcombobox(observed.all, selected = 0, 
                                  cont = m.e.rows[[obs.i]])
    svalue(m.e.obs[[obs.i]]) <<- obs.name
    m.e.type[[obs.i]] <<- gcombobox(c("SFO", "FOMC", "DFOP", "HS", "SFORB"),
                                   cont = m.e.rows[[obs.i]])
    svalue(m.e.type[[obs.i]]) <<- m[[m.cur]]$spec[[obs.i]]$type
    glabel("to", cont = m.e.rows[[obs.i]]) 
    obs.to <<- ifelse(is.null(m[[m.cur]]$spec[[obs.i]]$to), "",
                 paste(m[[m.cur]]$spec[[obs.i]]$to, collapse = ", "))
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
# Plotting and fitting {{{1
show_plot <- function(type, default = FALSE) {
  Parameters <- f.gg.parms[,]
  Parameters.de <- subset(Parameters, Type == "deparm", type)
  stateparms <- subset(Parameters, Type == "state")[[type]]
  deparms <- as.numeric(Parameters.de[[type]])
  names(deparms) <- rownames(Parameters.de)
  if (type == "Initial" & default == FALSE) {
    ftmp <<- suppressWarnings(mkinfit(ftmp$mkinmod, 
                                      override(ds[[ds.i]]$data),
                                      state.ini = stateparms, 
                                      parms.ini = deparms,
                                      err = "err", 
                                      control.modFit = list(maxiter = 0)))
    ftmp$ds.index <<- ds.i
    ftmp$ds <<- ds[[ds.i]]
  } 

  svg(tf, width = 7, height = 5)
  plot(ftmp, main = ftmp$ds$title,
       xlab = ifelse(ftmp$ds$time_unit == "", "Time", 
                     paste("Time in", ftmp$ds$time_unit)),
       ylab = ifelse(ds[[ds.i]]$unit == "", "Observed", 
                     paste("Observed in", ftmp$ds$unit)))
  dev.off()
  svalue(plot.gs) <<- tf
}
get_Parameters <- function(stmp, optimised)
{
  pars <- rbind(stmp$start[1:2], stmp$fixed)

  pars$fixed <- c(rep(FALSE, length(stmp$start$value)),
                  rep(TRUE, length(stmp$fixed$value)))
  pars$name <- rownames(pars)
  Parameters <- data.frame(Name = pars$name,
                           Type = pars$type,
                           Initial = pars$value,
                           Fixed = pars$fixed,
                           Optimised = as.numeric(NA))
  Parameters <- rbind(subset(Parameters, Type == "state"),
                      subset(Parameters, Type == "deparm"))
  rownames(Parameters) <- Parameters$Name
  if (optimised) {
    Parameters[rownames(stmp$bpar), "Optimised"] <- stmp$bpar[, "Estimate"]
  }
  return(Parameters)
}
run_fit <- function() {
  Parameters <- f.gg.parms[,]
  Parameters.de <- subset(Parameters, Type == "deparm")
  deparms <- Parameters.de$Initial
  names(deparms) <- Parameters.de$Name
  defixed <- names(deparms[Parameters.de$Fixed])
  Parameters.ini <- subset(Parameters, Type == "state")
  iniparms <- Parameters.ini$Initial
  names(iniparms) <- sub("_0", "", Parameters.ini$Name)
  inifixed <- names(iniparms[Parameters.ini$Fixed])
  weight <- svalue(f.gg.opts.weight)
  if (weight == "manual") {
    err = "err"
  } else {
    err = NULL
  }
  reweight.method <- svalue(f.gg.opts.reweight.method)
  if (reweight.method == "none") reweight.method = NULL
  ftmp <<- mkinfit(ftmp$mkinmod, override(ds[[ds.i]]$data),
                   state.ini = iniparms,
                   fixed_initials = inifixed,
                   parms.ini = deparms, 
                   fixed_parms = defixed,
                   solution_type = svalue(f.gg.opts.st),
                   atol = as.numeric(svalue(f.gg.opts.atol)),
                   rtol = as.numeric(svalue(f.gg.opts.rtol)),
                   weight = weight,
                   err = err,
                   reweight.method = reweight.method,
                   reweight.tol = as.numeric(svalue(f.gg.opts.reweight.tol)),
                   reweight.max.iter = as.numeric(svalue(f.gg.opts.reweight.max.iter))
                   )
  ftmp$ds.index <<- ds.i
  ftmp$ds <<- ds[[ds.i]]
  stmp <<- summary(ftmp)
  show_plot("Optimised")
  svalue(f.gg.opts.st) <- ftmp$solution_type
  svalue(f.gg.opts.weight) <- ftmp$weight.ini
  f.gg.parms[,] <- get_Parameters(stmp, TRUE)
  svalue(f.gg.summary) <- capture.output(stmp)
}
ds.i <- m.i <- "1"
f.cur <- "0"

# GUI widgets {{{2
pf <- gframe("Dataset 1, Model SFO", horizontal = FALSE, 
             cont = center, label = "Plotting and fitting")

# Mid group with plot and options {{{3
f.gg.mid <- ggroup(cont = pf)
ftmp <- suppressWarnings(mkinfit(m[[m.cur]], override(ds[[ds.i]]$data), 
                                 err = "err", 
                                 control.modFit = list(maxiter = 0)))
ftmp$ds.index = ds.i
ftmp$ds = ds[[ds.i]]
stmp <- summary(ftmp)
Parameters <- get_Parameters(stmp, FALSE)
tf <- get_tempfile(ext=".svg")
svg(tf, width = 7, height = 5)
plot(ftmp)
dev.off()
plot.gs <- gsvg(tf, container = f.gg.mid, width = 420, height = 300)
f.gg.opts <- gformlayout(cont = f.gg.mid)
solution_types <- c("auto", "analytical", "eigen", "deSolve")
f.gg.opts.st <- gcombobox(solution_types, selected = 1, 
                          label = "solution_type", width = 200, 
                          cont = f.gg.opts)
f.gg.opts.atol <- gedit(ftmp$atol, label = "atol", width = 20, 
                         cont = f.gg.opts)
f.gg.opts.rtol <- gedit(ftmp$rtol, label = "rtol", width = 20, 
                         cont = f.gg.opts)
weights <- c("manual", "none", "std", "mean")
f.gg.opts.weight <- gcombobox(weights, selected = 1, label = "weight", 
                              width = 200, cont = f.gg.opts)
f.gg.opts.reweight.method <- gcombobox(c("none", "obs"), selected = 1,
                                       label = "reweight.method",
                                       width = 200,
                                       cont = f.gg.opts)
f.gg.opts.reweight.tol <- gedit(1e-8, label = "reweight.tol",
                                 width = 20, cont = f.gg.opts)
f.gg.opts.reweight.max.iter <- gedit(10, label = "reweight.max.iter",
                                 width = 20, cont = f.gg.opts)

# Dataframe with initial and optimised parameters {{{3
f.gg.parms <- gdf(Parameters, width = 420, height = 300, cont = pf, 
                   do_add_remove_buttons = FALSE)
f.gg.parms$set_column_width(1, 200)
f.gg.parms$set_column_width(2, 50)
f.gg.parms$set_column_width(3, 60)
f.gg.parms$set_column_width(4, 50)
f.gg.parms$set_column_width(5, 60)

# Row with buttons {{{3
f.gg.buttons <- ggroup(cont = pf)
show.initial.gb <- gbutton("Show initial", 
                           handler = function(h, ...) show_plot("Initial"),
                           cont = f.gg.buttons)
tooltip(show.initial.gb) <- "Show model with current inital settings for current dataset"
run.fit.gb <- gbutton("Run", 
                      handler = function(h, ...) run_fit(), cont =
                      f.gg.buttons)
tooltip(run.fit.gb) <- "Fit with current settings on the current dataset, with the original model"
keep.fit.gb <- gbutton("Keep", 
                       handler = function(h, ...) {
                            f.cur <<- as.character(length(f) + 1)
                            f[[f.cur]] <<- ftmp
                            s[[f.cur]] <<- stmp
                            update_f.df()
                            f.gtable[,] <<- f.df
                          }, cont = f.gg.buttons)
tooltip(keep.fit.gb) <- "Store the optimised model with all settings and the current dataset in the fit list"

show.plots.gb <- gbutton("Show plots", 
                      handler = function(h, ...) show_plot_window(), 
                      cont = f.gg.buttons)
tooltip(show.plots.gb) <- "Show a window with plots including residual plots"

# Summary {{{3
f.gg.summary <- gtext(capture.output(stmp), cont = pf,
                      use.codemirror = TRUE)


delete.fit.gb <- gbutton("Delete", handler = function(h, ...) {
          if (length(f) > 0) {
            f[[f.cur]] <<- NULL
            s[[f.cur]] <<- NULL
          }
          if (length(f) > 0) {
            names(f) <<- as.character(1:length(f))
            names(s) <<- as.character(1:length(f))
            update_f.df()
            f.cur <<- "1"
            ftmp <<- f[[f.cur]]
            stmp <<- s[[f.cur]]
            ds.i <<- ftmp$ds.index
            update_plotting_and_fitting()
          } else {
            f.df <<- f.df.empty
            f.cur <<- "0"
          }
          f.gtable[,] <<- f.df
        }, cont = f.gg.buttons)
tooltip(delete.fit.gb) <- "Delete the currently loaded fit from the fit list"

# Update the plotting and fitting area {{{3
update_plotting_and_fitting <- function() {
  svalue(pf) <- paste0("Fit ", f.cur, ": Dataset ", ftmp$ds.index, 
                       ", Model ", ftmp$mkinmod$name)
  show_plot("Optimised")  
  svalue(f.gg.opts.st) <- ftmp$solution_type
  svalue(f.gg.opts.weight) <- ftmp$weight.ini
  svalue(f.gg.opts.reweight.method) <- ifelse(is.null(ftmp$reweight.method),
                                                      "none", 
                                                      ftmp$reweight.method)
  svalue(f.gg.opts.reweight.tol) <- ftmp$reweight.tol
  svalue(f.gg.opts.reweight.max.iter) <- ftmp$reweight.max.iter
  f.gg.parms[,] <- get_Parameters(stmp, TRUE)
  svalue(f.gg.summary) <- capture.output(stmp)
}

# Show plot window with residual plots {{{3
show_plot_window <- function(h, ...) {
  n.obs = length(ftmp$mkinmod$spec)
  obs.vars = names(ftmp$mkinmod$spec)
  parent = obs.vars[1]
  if(n.obs == 1) {
    n.rows = 1
    ps = 7
  } else {
    n.rows = 1 + ceiling(n.obs / 2)
    ps = 10
  }
  imgwidth = 800
  imgheight = 360 * n.rows
  pw <- gwindow("Plot window", parent = w,
                width = imgwidth + 20, height = imgheight + 100)
  pwg <- ggroup(cont = pw, horizontal = FALSE)
  make_plots <- function() {
    par(mfrow = c(n.rows, 2))
    plot(ftmp, main = ftmp$ds$title,
         xlab = ifelse(ftmp$ds$time_unit == "", "Time", 
                       paste("Time in", ftmp$ds$time_unit)),
         ylab = ifelse(ds[[ds.i]]$unit == "", "Observed", 
                       paste("Observed in", ftmp$ds$unit)))
    if (n.obs > 1) {
      plot(ftmp, legend = FALSE,
           main = paste0("Zoomed in on metabolite", 
                         ifelse(n.obs > 2, "s", "")),
           xlab = ifelse(ftmp$ds$time_unit == "", "Time", 
                         paste("Time in", ftmp$ds$time_unit)),
           ylab = ifelse(ds[[ds.i]]$unit == "", "Observed", 
                         paste("Observed in", ftmp$ds$unit)),
           ylim = c(0, max(subset(ftmp$data, 
                        variable != parent)$observed)))
      for (met in obs.vars[-1]) {
        mkinresplot(ftmp, met, legend = FALSE,
                    main = paste("Residual plot for", met))
      }
    } else {
      mkinresplot(ftmp, parent, legend = FALSE,
                  main = paste("Residual plot for", parent),
           xlab = ifelse(ftmp$ds$time_unit == "", "Time", 
                         paste("Time in", ftmp$ds$time_unit)),
           ylab = ifelse(ds[[ds.i]]$unit == "", "Residuals", 
                         paste("Residuals in", ftmp$ds$unit)))
    }
  }

  tf2 <- get_tempfile(ext = ".png")
  png(tf2,  width = imgwidth / 50 , height = imgheight / 50,
    units = "cm", res = 300, pointsize = ps)
  make_plots()
  dev.off()

  ghtml(paste0("<img width='", imgwidth, "' height='", imgheight,
               "' src='", get_tempfile_url(tf2), "' />"), 
        cont = pwg)
}
# vim: set foldmethod=marker ts=2 sw=2 expandtab: {{{1
