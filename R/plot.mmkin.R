# Copyright (C) 2015-2016,2019 Johannes Ranke
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

plot.mmkin <- function(x, main = "auto", legends = 1, 
                       resplot = c("time", "errmod"),
                       show_errmin = TRUE,
                       errmin_var = "All data", errmin_digits = 3,
                       cex = 0.7, rel.height.middle = 0.9, ...) {
  n.m <- nrow(x)
  n.d <- ncol(x)

  resplot <- match.arg(resplot)

  # We can handle either a row (different models, same dataset)
  # or a column (same model, different datasets)
  if (n.m > 1 & n.d > 1) stop("Please select fits either for one model or for one dataset")
  if (n.m == 1 & n.d == 1) loop_over = "none"
  if (n.m > 1) loop_over <- "models"
  if (n.d > 1) loop_over <- "datasets"
  n.fits <- length(x)

  # Set the main plot titles from the names of the models or the datasets
  # Will be integer indexes if no other names are present in the mmkin object
  if (main == "auto") {
    main = switch(loop_over,
                  none = paste(rownames(x), colnames(x)),
                  models = colnames(x),
                  datasets = rownames(x))
  }

  oldpar <- par(no.readonly = TRUE)

  # Set relative plot heights, so the first and the last plot are the norm
  # and the middle plots (if n.fits >2) are smaller by rel.height.middle
  rel.heights <- if (n.fits > 2) c(1, rep(rel.height.middle, n.fits - 2), 1)
                 else rep(1, n.fits)
  layout(matrix(1:(2 * n.fits), n.fits, 2, byrow = TRUE), heights = rel.heights)

  par(cex = cex)

  for (i.fit in 1:n.fits) {

    # Margins for top row of plots when we have more than one row
    # Reduce bottom margin by 2.1 - hides x axis legend
    if (i.fit == 1 & n.fits > 1) {
      par(mar = c(3.0, 4.1, 4.1, 2.1))
    }

    # Margins for middle rows of plots, if any
    if (i.fit > 1 & i.fit < n.fits) {
      # Reduce top margin by 2 after the first plot as we have no main title,
      # reduced plot height, therefore we need rel.height.middle in the layout
      par(mar = c(3.0, 4.1, 2.1, 2.1))
    }

    # Margins for bottom row of plots when we have more than one row
    if (i.fit == n.fits & n.fits > 1) {
      # Restore bottom margin for last plot to show x axis legend
      par(mar = c(5.1, 4.1, 2.1, 2.1))
    }

    fit <- x[[i.fit]]
    plot(fit, legend = legends == i.fit, ...)

    title(main, outer = TRUE, line = -2)

    fit_name <- switch(loop_over,
                       models = rownames(x)[i.fit],
                       datasets = colnames(x)[i.fit],
                       none = "")

    if (show_errmin) {
      chi2 <- signif(100 * mkinerrmin(fit)[errmin_var, "err.min"], errmin_digits)

      # Use LateX if the current plotting device is tikz
      if (names(dev.cur()) == "tikz output") {
        chi2_text <- paste0(fit_name, " $\\chi^2$ error level = ", chi2, "\\%")
      } else {
        chi2_perc <- paste0(chi2, "%")
        chi2_text <- bquote(.(fit_name) ~ chi^2 ~ "error level" == .(chi2_perc))
      }
      mtext(chi2_text, cex = cex, line = 0.4)
    }

    if (resplot == "time") {
      mkinresplot(fit, legend = FALSE, ...)
    } else {
      mkinerrplot(fit, legend = FALSE, ...)
    }
    if (show_errmin) {
      mtext(paste(fit_name, "residuals"), cex = cex, line = 0.4)
    }
  }

  par(oldpar, no.readonly = TRUE)
}
