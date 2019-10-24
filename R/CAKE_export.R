#' Export a list of datasets format to a CAKE study file
#' 
#' In addition to the datasets, the pathways in the degradation model can be
#' specified as well.
#' 
#' @param ds A named list of datasets in long format as compatible with
#'   \code{\link{mkinfit}}.
#' @param map A character vector with CAKE compartment names (Parent, A1, ...),
#'   named with the names used in the list of datasets.
#' @param links An optional character vector of target compartments, named with
#'   the names of the source compartments. In order to make this easier, the
#'   names are used as in the datasets supplied.
#' @param filename Where to write the result. Should end in .csf in order to be
#'   compatible with CAKE.
#' @param path An optional path to the output file.
#' @param overwrite If TRUE, existing files are overwritten.
#' @param study The name of the study.
#' @param description An optional description.
#' @param time_unit The time unit for the residue data.
#' @param res_unit The unit used for the residues.
#' @param comment An optional comment.
#' @param date The date of file creation.
#' @param optimiser Can be OLS or IRLS.
#' @importFrom utils write.table
#' @return The function is called for its side effect.
#' @author Johannes Ranke
#' @export
CAKE_export <- function(ds, map = c(parent = "Parent"),
  links = NA,
  filename = "CAKE_export.csf", path = ".", overwrite = FALSE,
  study = "Codlemone aerobic soil degradation",
  description = "",
  time_unit = "days",
  res_unit = "% AR",
  comment = "Created using mkin::CAKE_export",
  date = Sys.Date(),
  optimiser = "IRLS")
{
  file <- file.path(path, filename)
  if (file.exists(file) & !overwrite) stop(file, " already exists, stopping")
  csf <- file(file, encoding = "latin1", open = "w+")
  on.exit(close(csf))

  add <- function(x) cat(paste0(x, "\r\n"), file = csf, append = TRUE)
  add0 <- function(x) cat(x, file = csf, append = TRUE)

  add("[FileInfo]")
  add("CAKE-Version: 3.3 (Release)")
  add(paste("Name:", study))
  add(paste("Description:", description))
  add(paste("MeasurementUnits:", res_unit))
  add(paste("TimeUnits:", time_unit))
  add(paste("Comments:", comment))
  add(paste("Date:", date))
  add(paste("Optimiser:", optimiser))
  add("")

  add("[Data]")

  for (i in seq_along(ds)) {
    add(paste("NewDataSet:", names(ds)[i]))
    d <- mkin_long_to_wide(ds[[i]])
    names(d) <- c("Time", map[names(d)[-1]])
    write.table(d, csf,
      sep = "\t", col.names = TRUE,
      row.names = FALSE,
      quote = FALSE, eol = "\r\n", na = "")
    add("")
  }

  if (!is.na(links)) {
    add("")
    add("[Model]")
    add(paste0("ParentCompartment: Parent\t", names(map)[1], "\t", names(map)[1]))
    for (name in names(map)[-1]) {
      add(paste0("Compartment: ", map[name], "\t", name, "\t", name))
    }
    for (li in names(links)) {
      add(paste0("Link: ", map[li], "\t", map[links[li]], "\t0.5\t0\t1\tFree\tExplicit"))
    }

  }

  add("")
  add("[ComponentNames]")
  for (name in names(map)) {
    add(paste0(map[name], ":", name))
  }

}
