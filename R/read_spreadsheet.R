#' Read datasets and relevant meta information from a spreadsheet file
#'
#' This function imports one dataset from each sheet of a spreadsheet file.
#' These sheets are selected based on the contents of a sheet 'Datasets', with
#' a column called 'Dataset Number', containing numbers identifying the dataset
#' sheets to be read in. In the second column there must be a grouping
#' variable, which will often be named 'Soil'. Optionally, time normalization
#' factors can be given in columns named 'Temperature' and 'Moisture'.
#'
#' There must be a sheet 'Compounds', with columns 'Name' and 'Acronym'.
#' The first row read after the header read in from this sheet is assumed
#' to contain name and acronym of the parent compound.
#'
#' The dataset sheets should be named using the dataset numbers read in from
#' the 'Datasets' sheet, i.e. '1', '2', ... . In each dataset sheet, the name
#' of the observed variable (e.g. the acronym of the parent compound or
#' one of its transformation products) should be in the first column,
#' the time values should be in the second colum, and the observed value
#' in the third column.
#'
#' In case relevant covariate data are available, they should be given
#' in a sheet 'Covariates', containing one line for each value of the grouping
#' variable specified in 'Datasets'. These values should be in the first
#' column and the column must have the same name as the second column in
#' 'Datasets'. Covariates will be read in from columns four and higher.
#' Their names should preferably not contain special characters like spaces,
#' so they can be easily used for specifying covariate models.
#'
#' A similar data structure is defined as the R6 class [mkindsg], but
#' is probably more complicated to use.
#'
#' @param path Absolute or relative path to the spreadsheet file
#' @param valid_datasets Optional numeric index of the valid datasets, default is
#' to use all datasets
#' @param parent_only Should only the parent data be used?
#' @param normalize Should the time scale be normalized using temperature
#' and moisture normalisation factors in the sheet 'Datasets'?
#' @export
read_spreadsheet <- function(path, valid_datasets = "all",
  parent_only = FALSE, normalize = TRUE)
{
  if (!requireNamespace("readxl", quietly = TRUE))
    stop("Please install the readxl package to use this function")

  # Read the compound table
  compounds <- readxl::read_excel(path, sheet = "Compounds")
  parent <- compounds[1, ]$Acronym

  # Read in meta information
  ds_meta <- readxl::read_excel(path, sheet = "Datasets")
  ds_meta["Dataset Number"] <- as.character(ds_meta[["Dataset Number"]])

  # Select valid datasets
  if (valid_datasets[1] == "all") valid_datasets <- 1:nrow(ds_meta)
  ds_numbers_valid <- ds_meta[valid_datasets, ]$`Dataset Number`
  grouping_factor <- names(ds_meta[2]) # Often "Soil"

  # Read in valid datasets
  ds_raw <- lapply(ds_numbers_valid,
    function(dsn) readxl::read_excel(path, sheet = as.character(dsn)))

  # Make data frames compatible with mmkin
  ds_tmp <- lapply(ds_raw, function(x) {
    ds_ret <- x[1:3] |>
      rlang::set_names(c("name", "time", "value")) |>
      transform(value = as.numeric(value))
  })
  names(ds_tmp) <- ds_numbers_valid

  # Normalize with temperature and moisture correction factors
  if (normalize) {
    ds_norm <- lapply(ds_numbers_valid, function(ds_number) {
      f_corr <- as.numeric(ds_meta[ds_number, c("Temperature", "Moisture")])
      ds_corr <- ds_tmp[[ds_number]] |>
        transform(time = time * f_corr[1] * f_corr[2])
      return(ds_corr)
    })
  } else {
    ds_norm <- ds_tmp
  }
  names(ds_norm) <- ds_numbers_valid

  # Select parent data only if requested
  if (parent_only) {
    ds_norm <- lapply(ds_norm, function(x) subset(x, name == parent))
    compounds <- compounds[1, ]
  }

  # Create a single long table to combine datasets with the same group name
  ds_all <- vctrs::vec_rbind(!!!ds_norm, .names_to = "Dataset Number")
  ds_all_group <- merge(ds_all, ds_meta[c("Dataset Number", grouping_factor)])
  groups <- unique(ds_meta[valid_datasets, ][[grouping_factor]])

  ds <- lapply(groups, function(x) {
      ret <- ds_all_group[ds_all_group[[grouping_factor]] == x, ]
      ret[c("name", "time", "value")]
    }
  )
  names(ds) <- groups

  # Get covariates
  covariates_raw <- readxl::read_excel(path, sheet = "Covariates")
  covariates <- as.data.frame(covariates_raw[4:ncol(covariates_raw)])
  nocov <- setdiff(groups, covariates_raw[[1]])
  if (length(nocov) > 0) {
    message("Did not find covariate data for ", paste(nocov, collapse = ", "))
    message("Not returning covariate data")
    attr(ds, "covariates") <- NULL
  } else {
    rownames(covariates) <- covariates_raw[[1]]
    covariates <- covariates[which(colnames(covariates) != "Remarks")]
    # Attach covariate data if available
    attr(ds, "covariates") <- covariates[groups, , drop = FALSE]
  }

  # Attach the compound list to support automatic model building
  attr(ds, "compounds") <- as.data.frame(compounds)

  return(ds)
}
