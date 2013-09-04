# Simple gWidgetsWWW2 GUI for mkin
w <- gwindow("Simple R GUI for kinetic evaluations")
sb <- gstatusbar("Powered by gWidgetsWWW2 and Rook", cont = w)
g <- gframe("Simple R GUI for kinetic evaluations", cont = w,
            use.scrollwindow = TRUE, horizontal = FALSE)

# Project definition {{{1
prg <- gexpandgroup("Project definition", cont = g)
visible(prg) <- TRUE

make_observed <- function(observed.df) {
  if (!exists("observed.df")) {
    n.observed <- 2
    observed.names = c("parent", paste("M", 1:(n.observed - 1), sep=""))
    observed.df = data.frame(Index = 1:n.observed, 
                             Name = observed.names,
                             stringsAsFactors = FALSE)
  }

  observed.gdf <- gdf(observed.df, 
                     name = "Names of observed variables", 
                     width = 500, height = 250, cont = prg)
  observed.gdf$set_column_width(1, 40)
}

f <- gfile(text = "Set project file", cont = prg,
           handler = function(h, ...) 
           {
             tmpfile <- normalizePath(svalue(h$obj), winslash = "/")
             load(tmpfile)
             make_observed(observed.df)
           })

# Editable table of studies {{{1
stg <- gexpandgroup("Studies", cont = g)
visible(stg) <- FALSE
studies.df <- data.frame(Index = as.integer(1), 
                         Author = "FOCUS kinetics workgroup", 
                         Year = "2006", 
                         Title = "FOCUS Kinetics",
                         Datasets = as.integer(3),
                         stringsAsFactors = FALSE)
studies.gdf <- gdf(studies.df, 
                   name = "Editable list of studies in the project",
                   width = "auto", 
                   cont = stg)
studies.gdf$set_column_width(1, 40)
studies.gdf$set_column_width(2, 200)

# Table of datasets to select them for editing {{{1
dsg <- gexpandgroup("Datasets", cont = g)
visible(dsg) <- FALSE
ds.df <- data.frame(Index = 1:3, 
                    Study = as.integer(1), 
                    Title = paste("FOCUS dataset", c("A", "B", "C")),
                    icon = asIcon(c("editor", "editor", "editor")),
                    stringsAsFactors = FALSE)
ds.gtable <- gtable(ds.df, width = "auto", cont = dsg)
ds.gtable$set_column_width(1, 40)
ds.gtable$set_column_width(2, 40)
ds.gtable$set_column_width(3, 200)


# 1}}}
# vim: set foldmethod=marker ts=2 sw=2 expandtab:
