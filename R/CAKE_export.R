# Copyright (C) 2019 Johannes Ranke
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
