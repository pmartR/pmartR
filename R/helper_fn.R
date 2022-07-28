# Functions to get omicsData attributes ----------------------------------------

#' Fetch the data_info attribute
#'
#' Retrieves the values in the data_info attribute from an omicsData object.
#'
#' @param omicsData An object of class pepData, proData, metabData, lipidData,
#'                  or nmrData.
#'
#' @return A list containing seven elements:
#'         \itemize{
#'
#'           \item data_scale -- A Character string indicating the scale of the
#'           data in \code{e_data}.
#'
#'           \item norm_info -- A list containing a single element indicating
#'           whether the data in \code{e_data} have been normalized.
#'
#'           \item num_edata -- The number of unique entries present in the
#'           \code{edata_cname} column in \code{e_data}.
#'
#'           \item num_miss_obs -- An integer. The number of missing
#'           observations in \code{e_data}.
#'
#'           \item prop_missing -- A number between 0 and 1. The proportion of
#'           missing observations in \code{e_data}.
#'
#'           \item num_samps -- An integer indicating the number of samples or
#'           columns (excluding the identifier column \code{edata_cname}) in
#'           \code{e_data}.
#'
#'           \item data_types -- A character string describing the type of data
#'           in \code{e_data}.
#'
#'         }
#'
#' @export
#'
get_data_info <- function (omicsData) {

  # Check class of omicsData.
  if (!inherits(omicsData, c("pepData",
                             "proData",
                             "metabData",
                             "lipidData",
                             "nmrData"))) {

    # Lay down an error in the console.
    stop (paste("omicsData must be of class 'pepData', 'proData',",
                "'metabData', 'lipidData', or 'nmrData'",
                sep = " "))

  }

  # Extract and return the data_info attribute.
  return (attr(omicsData, 'data_info'))

}

#' Fetch the normalization status of the data
#'
#' This function returns the norm_info element of the data_info attribute
#' indicating whether the data have been normalized.
#'
#' @param omicsObject an object of the class 'pepData', 'proData', 'metabData',
#'   'lipidData', 'nmrData', 'statRes', or 'trellData', usually created by
#'   \code{\link{as.pepData}}, \code{\link{as.proData}},
#'   \code{\link{as.metabData}}, \code{\link{as.lipidData}},
#'   \code{\link{as.nmrData}}, \code{\link{imd_anova}}, or
#'   \code{\link{format_data}} respectively.
#'
#' @return A logical value indicating whether the data have been normalized.
#'
#' @rdname get_data_norm
#'
#' @export
#'
get_data_norm <- function (omicsObject) {

  # check that omicsObject is of appropriate class #
  if(!inherits(omicsObject, c("pepData", "proData", "metabData", "lipidData",
                              "nmrData", "statRes", "trellData")))
    stop (paste("omicsData must be of class 'pepData', 'proData', 'metabData',",
                "'lipidData', 'nmrData', 'statRes', or 'trellData'",
                sep = " "))

  return (attr(omicsObject, "data_info")$norm_info$is_normalized)

}

#' Fetch the meta_info attribute
#'
#' Retrieves the values in the meta_info attribute from an omicsData object.
#'
#' @param omicsData An object of class pepData, proData, metabData, lipidData,
#'                  or nmrData.
#'
#' @return A list containing two elements:
#'         \itemize{
#'
#'           \item meta_data -- Logical. Indicates if the \code{e_meta} data
#'           frame was provided.
#'
#'           \item num_emeta -- The number of unique entries present in the
#'           \code{emeta_cname} column in \code{e_meta}.
#'
#'         }
#'
#' @export
#'
get_meta_info <- function (omicsData) {

  # Check class of omicsData.
  if (!inherits(omicsData, c("pepData",
                             "proData",
                             "metabData",
                             "lipidData",
                             "nmrData"))) {

    # Lay down an error in the console.
    stop (paste("omicsData must be of class 'pepData', 'proData',",
                "'metabData', 'lipidData', or 'nmrData'",
                sep = " "))

  }

  # Extract and return the meta_info attribute.
  return (attr(omicsData, 'meta_info'))

}

#' Fetch the filters attribute
#'
#' Retrieves the values in the filters attribute from an omicsData object.
#'
#' @param omicsData An object of class pepData, proData, metabData, lipidData,
#'                  or nmrData.
#'
#' @return A list containing filter class objects. Each element in this list
#'         corresponds to a filter applied to the data. The filters will be
#'         listed in the order they were applied. A filter object contains two
#'         elements:
#'         \itemize{
#'
#'           \item threshold -- The threshold used to filter \code{e_data}. This
#'           value depends on the type of filter applied.
#'
#'           \item filtered -- A vector containing the identifiers from the
#'           \code{edata_cname} column that will be filtered.
#'
#'           \item method -- A character string indicating the type of method
#'           used to filter. This only applies when imdanova_filter is used.
#'
#'         }
#'
#' @export
#'
get_filters <- function (omicsData) {

  # Check class of omicsData.
  if (!inherits(omicsData, c("pepData",
                             "proData",
                             "metabData",
                             "lipidData",
                             "nmrData"))) {

    # Lay down an error in the console.
    stop (paste("omicsData must be of class 'pepData', 'proData',",
                "'metabData', 'lipidData', or 'nmrData'",
                sep = " "))

  }

  # Send a message if the filters object is empty.
  if (length(attr(omicsData, 'filters')) == 0) {

    message("No filters have been applied.")

    # Return NULL if there are no filters for the pmart app people.
    return (NULL)

  } else {

    # Extract and return the filters attribute.
    return (attr(omicsData, 'filters'))

  }

}

# Extracts the types of filters that have been applied. This function will be
# used at the beginning of the applyFilt function to give a warning if the
# same type of filter has already been applied.
get_filter_type <- function (omicsData) {

  # Store the filters attribute from the previous filters (if any).
  pFilters <- attr(omicsData, 'filters')

  # Determine the length of the filters attribute. The length corresponds to the
  # number of filters applied previously.
  nFilters <- length(pFilters)

  # Check the length of filters.
  if (nFilters == 0) {

    # If no filters have been applied return NULL.
    return (NULL)

  } else {

    # Fabricate a vector to hold the filter types.
    pTypes <- vector(length = nFilters)

    # Loop through each filter.
    for (e in 1:nFilters) {

      # Extract the filter type.
      pTypes[[e]] <- pFilters[[e]]$type

    }

    # Return the character vector containing the filter types!!
    return (pTypes)

  }

}

#' Fetch the check.names attribute
#'
#' Retrieves the value in the check.names attribute from an omicsData object.
#'
#' @param omicsData An object of class pepData, proData, metabData, lipidData,
#'                  or nmrData.
#'
#' @return A logical value indicating if the syntax of the column names in a
#'         data frame should be checked. See \code{\link[base]{data.frame}} for
#'         more details.
#'
#' @export
#'
get_check_names <- function (omicsData) {

  # Check class of omicsData.
  if (!inherits(omicsData, c("pepData",
                             "proData",
                             "metabData",
                             "lipidData",
                             "nmrData"))) {

    # Lay down an error in the console.
    stop (paste("omicsData must be of class 'pepData', 'proData',",
                "'metabData', 'lipidData', or 'nmrData'",
                sep = " "))

  }

  # Extract and return the check.names attribute.
  return (attr(omicsData, 'check.names'))

}

#' Fetch the isobaric_info attribute
#'
#' Retrieves the values in the isobaric_info attribute from an omicsData object.
#'
#' @param omicsData An object of class pepData, proData, metabData, lipidData,
#'                  or nmrData.
#'
#' @return A list containing the following six elements:
#'         \itemize{
#'
#'           \item exp_cname --
#'
#'           \item channel_cname --
#'
#'           \item refpool_channel --
#'
#'           \item refpool_cname --
#'
#'           \item refpool_notation --
#'
#'           \item norm_info -- A list containing a single logical element that
#'           indicates whether the data have been normalized to a reference
#'           pool.
#'
#'         }
#'
#' @export
#'
get_isobaric_info <- function (omicsData) {

  # Check class of omicsData.
  if (!inherits(omicsData, c("isobaricpepData",
                             "pepData"))) {

    # Lay down an error in the console.
    stop (paste("omicsData must be of class 'isobaricpepData' and 'pepData'",
                sep = " "))

  }

  # Extract and return the isobaric_info attribute.
  return (attr(omicsData, 'isobaric_info'))

}

#' Fetch the isobaric normalization info
#'
#' This function returns the norm_info element of the isobaric_info attribute
#' which indicates if the data have been isobaric normalized.
#'
#' @param omicsData an object of the class 'pepData', 'isobaricpepData' or
#'   'proData', usually created by \code{\link{as.pepData}},
#'   \code{\link{as.isobaricpepData}}.
#'
#' @return A logical value indicating whether the data have been isobaric
#'   normalized.
#'
#' @rdname get_isobaric_norm
#'
#' @export
#'
get_isobaric_norm <- function (omicsData) {

  # check that omicsData is of appropriate class #
  if (!inherits(omicsData, c("pepData", "proData", "isobaricpepData")))
    stop("omicsData must be of class 'pepData', 'proData' or 'isobaricpepData'")

  return (attr(omicsData, "isobaric_info")$norm_info$is_normalized)

}

#' Fetch the nmr_info attribute
#'
#' Retrieves the values in the nmr_info attribute from an omicsData object.
#'
#' @param omicsData An object of class pepData, proData, metabData, lipidData,
#'                  or nmrData.
#'
#' @return A list containing the following three elements:
#'         \itemize{
#'
#'           \item metabolite_name --
#'
#'           \item sample_property_cname --
#'
#'           \item norm_info -- A list containing two logical elements that
#'           indicate i) whether the data have been normalized to a spiked in
#'           metabolite or to a property taking sample-specific values and ii)
#'           whether the data have been back transformed so the values are on a
#'           similar scale to the raw values before normalization.
#'
#'         }
#'
#' @export
#'
get_nmr_info <- function (omicsData) {

  # Check class of omicsData.
  if (!inherits(omicsData, c("pepData",
                             "proData",
                             "metabData",
                             "lipidData",
                             "nmrData"))) {

    # Lay down an error in the console.
    stop (paste("omicsData must be of class 'pepData', 'proData',",
                "'metabData', 'lipidData', or 'nmrData'",
                sep = " "))

  }

  # Extract and return the nmr_info attribute.
  return (attr(omicsData, 'nmr_info'))

}

#' Fetch the NMR normalization info
#'
#' This function returns the norm_info element of the nmr_info attribute which
#' indicates if the data have been NMR normalized.
#'
#' @param omicsData an object of the class 'nmrData'.
#'
#' @return A logical value indicating whether the data have been NMR normalized.
#'
#' @rdname get_nmr_norm
#'
#' @export
#'
get_nmr_norm <- function (omicsData) {

  # check that omicsData is of appropriate class #
  if (!inherits(omicsData, "nmrData"))
    stop ("omicsData must be of class 'nmrData'")

  return (attr(omicsData, "nmr_info")$norm_info$is_normalized)

}

#' Fetch the group_DF attribute
#'
#' Retrieves the values in the group_DF attribute from an omicsData object.
#'
#' @param omicsData An object of class pepData, proData, metabData, lipidData,
#'                  or nmrData.
#'
#' @return A data.frame with columns for sample ID and group. If two main
#'         effects are provided the original main effect levels for each sample
#'         are returned as the third and fourth columns of the data frame.
#'         Additionally, the covariates provided will be listed as attributes of
#'         this data frame.
#'
#' @export
#'
get_group_DF <- function (omicsData) {

  # Check class of omicsData.
  if (!inherits(omicsData, c("pepData",
                             "proData",
                             "metabData",
                             "lipidData",
                             "nmrData"))) {

    # Lay down an error in the console.
    stop (paste("omicsData must be of class 'pepData', 'proData',",
                "'metabData', 'lipidData', or 'nmrData'",
                sep = " "))

  }

  # Extract and return the group_DF attribute.
  return (attr(omicsData, 'group_DF'))

}

#' Fetch the original data scale
#'
#' Retrieves the character string indicating the scale the data was originally
#' on when read into R.
#'
#' @param omicsObject an object of class 'pepData', 'proData', 'metabData',
#'  'lipidData', or 'nmrData'.
#'
#' @return A character string.
#'
#' @export
#'
get_data_scale_orig <- function (omicsObject){

  # Check that the input is an appropriate class.
  if (!inherits(omicsObject, c("pepData", "proData", "metabData",
                               "lipidData", "nmrData"))) {

    # Halt! You are using an unholy input object. Come back when your data is
    # decent.
    stop ("omicsObject is not an accepted class.")

  }

  return (attr(omicsObject, "data_info")$data_scale_orig)

}

#' Fetch the current data scale
#'
#' This function returns current data scale which may be different from the
#' original data scale (if \code{edata_transform} was used).
#'
#' @param omicsObject an object of the class 'pepData', 'proData', 'metabData',
#'   'lipidData', 'nmrData', 'statRes', or 'trellData', usually created by
#'   \code{\link{as.pepData}}, \code{\link{as.proData}},
#'   \code{\link{as.metabData}}, \code{\link{as.lipidData}},
#'   \code{\link{as.nmrData}}, \code{\link{imd_anova}}, or
#'   \code{\link{format_data}} respectively.
#'
#' @return a character string describing data scale
#'
#' @rdname get_data_scale
#'
#' @export
#'
get_data_scale <- function (omicsObject) {

  # check that omicsObject is of appropriate class #
  if(!inherits(omicsObject, c("pepData", "proData", "metabData", "lipidData",
                              "nmrData", "statRes", "trellData")))
    stop (paste("omicsData must be of class 'pepData', 'proData', 'metabData',",
                "'lipidData', 'nmrData', 'statRes', or 'trellData'",
                sep = " "))

  return(attr(omicsObject, "data_info")$data_scale)

}

#' Fetch the e_data column name
#'
#' This function returns the name of the column in e_data that contains the
#' biomolecule IDs.
#'
#' @param omicsObject an object of the class 'pepData', 'proData', 'metabData',
#'   'lipidData', 'nmrData', 'statRes', or 'trellData', usually created by
#'   \code{\link{as.pepData}}, \code{\link{as.proData}},
#'   \code{\link{as.metabData}}, \code{\link{as.lipidData}},
#'   \code{\link{as.nmrData}}, \code{\link{imd_anova}}, or
#'   \code{\link{format_data}} respectively.
#'
#' @return a character string describing e_data cname
#'
#' @rdname get_edata_cname
#'
#' @export
#'
get_edata_cname <- function (omicsObject) {

  # check that omicsObject is of appropriate class #
  if(!inherits(omicsObject, c("pepData", "proData", "metabData", "lipidData",
                              "nmrData", "statRes", "trellData")))
    stop (paste("omicsData must be of class 'pepData', 'proData', 'metabData',",
                "'lipidData', 'nmrData', 'statRes', or 'trellData'",
                sep = " "))

  return (attr(omicsObject, "cnames")$edata_cname)

}

#' Fetch the f_data column name
#'
#' This function returns the name of the column in f_data that contains the
#' names of the samples.
#'
#' @param omicsObject an object of the class 'pepData', 'proData', 'metabData',
#'   'lipidData', 'nmrData', 'statRes', or 'trellData', usually created by
#'   \code{\link{as.pepData}}, \code{\link{as.proData}},
#'   \code{\link{as.metabData}}, \code{\link{as.lipidData}},
#'   \code{\link{as.nmrData}}, \code{\link{imd_anova}}, or
#'   \code{\link{format_data}} respectively.
#'
#' @return a character string describing f_data cname
#'
#' @rdname get_fdata_cname
#'
#' @export
#'
get_fdata_cname<- function (omicsObject) {

  # check that omicsObject is of appropriate class #
  if(!inherits(omicsObject, c("pepData", "proData", "metabData", "lipidData",
                              "nmrData", "statRes", "trellData")))
    stop (paste("omicsData must be of class 'pepData', 'proData', 'metabData',",
                "'lipidData', 'nmrData', 'statRes', or 'trellData'",
                sep = " "))

  return (attr(omicsObject, "cnames")$fdata_cname)

}

#' Fetch the e_meta column name
#'
#' This function returns the name of the column in e_meta that contains the
#' mapping variable IDs.
#'
#' @param omicsObject an object of the class 'pepData', 'proData', 'metabData',
#'   'lipidData', 'nmrData', 'statRes', or 'trellData', usually created by
#'   \code{\link{as.pepData}}, \code{\link{as.proData}},
#'   \code{\link{as.metabData}}, \code{\link{as.lipidData}},
#'   \code{\link{as.nmrData}}, \code{\link{imd_anova}}, or
#'   \code{\link{format_data}} respectively.
#'
#' @return a character string describing e_meta cname
#'
#' @rdname get_emeta_cname
#'
#' @export
#'
get_emeta_cname <- function (omicsObject) {

  # check that omicsObject is of appropriate class #
  if(!inherits(omicsObject, c("pepData", "proData", "metabData", "lipidData",
                              "nmrData", "statRes", "trellData")))
    stop (paste("omicsData must be of class 'pepData', 'proData', 'metabData',",
                "'lipidData', 'nmrData', 'statRes', or 'trellData'",
                sep = " "))

  return (attr(omicsObject, "cnames")$emeta_cname)

}

# Functions to set omicsData attributes ----------------------------------------

# Create a set function that will return the value for each attribute. For
# example, set_data_info will perform all of the calculations to fill in the
# data_info attribute. These functions will be called in the as.xxx functions
# to create an xxxData object but can be used individually to update any one of
# the attributes at a later time.

# Sets/updates the values in the data_info attribute
#
# @param e_data
#
# @param edata_cname
#
# @param data_scale_orig
#
# @param data_scale
#
# @param data_types
#
# @param norm_info
#
# @param is_normalized
#
# @return A list containing all the elements in the data_info attribute:
#         (list all attributes)
#
set_data_info <- function (e_data,
                           edata_cname,
                           data_scale_orig,
                           data_scale,
                           data_types,
                           norm_info,
                           is_normalized,
                           batch_info,
                           is_bc) {

  # Identify the column number that contains the IDs.
  id_col <- which(names(e_data) == edata_cname)

  # Count the number of missing values in e_data.
  num_miss_obs <- sum(is.na(e_data[, -id_col]))

  # Calculate the proportion of missing data in e_data.
  prop_missing <- num_miss_obs / prod(dim(e_data[, -id_col]))

  # Extract the number of things (e.g., peptides, lipids, metabolites, ...).
  num_edata <- nrow(e_data)

  # Procure the number of samples.
  num_samps <- ncol(e_data) - 1

  # Set data normalization information.
  norm_info$is_normalized <- is_normalized

  batch_info$is_bc <- is_bc

  # Return all of the information that belongs in the data_info attribute.
  return (list(data_scale_orig = data_scale_orig,
               data_scale = data_scale,
               norm_info = norm_info,
               num_edata = num_edata,
               num_miss_obs = num_miss_obs,
               prop_missing = prop_missing,
               num_samps = num_samps,
               data_types = data_types,
               batch_info = batch_info))

}

# Sets/updates the values in the meta_info attribute
#
# @param e_data
#
# @param emeta_cname
#
# @return A list containing all the elements in the meta_info attribute:
#         (list all attributes)
#
set_meta_info <- function (e_meta,
                           emeta_cname) {

  # Determine if meta data is present.
  meta_data <- ifelse(is.null(e_meta), FALSE, TRUE)

  # Test if emeta_cname is not null.
  if (!is.null(emeta_cname)) {

    # Enumerate the number of unique proteins that map to a peptide in e_data.
    # When using other data (e.g., lipid or metabolite) this counts the number
    # of unique mapping variables associated with a biomolecule.
    num_emeta <- length(unique(e_meta[, emeta_cname]))

  } else {

    # If emeta_cname is null set the number of proteins to null.
    num_emeta <- NULL

  }

  # Return the list of meta attributes.
  return (list(meta_data = meta_data,
               num_emeta = num_emeta))

}

# Sets/updates the values in the isobaric_info attribute.
#
# @param exp_cname
#
# @param channel_cname
#
# @param refpool_channel
#
# @param refpool_cname
#
# @param refpool_notation
#
# @param norm_info
#
# @param isobaric_norm
#
# @return A list containing all the elements in the isobaric_info attribute:
#         (list all attributes)
#
set_isobaric_info <- function (exp_cname,
                               channel_cname,
                               refpool_channel,
                               refpool_cname,
                               refpool_notation,
                               norm_info,
                               isobaric_norm) {

  # Set the elements of the norm_info list.
  norm_info$is_normalized <- isobaric_norm

  # Return the list of isobaric_info attributes.
  return (list(exp_cname = exp_cname,
               channel_cname = channel_cname,
               refpool_channel = refpool_channel,
               refpool_cname = refpool_cname,
               refpool_notation = refpool_notation,
               norm_info = norm_info))

}

# Sets/updates the values in the nmr_info attribute.
#
# @param metabolite_name
#
# @param sample_property_cname
#
# @param norm_info
#
# @param nmr_norm
#
# @param backtransform
#
# @return A list containing all the elements in the nmr_info attribute:
#         (list all attributes)
#
set_nmr_info <- function (metabolite_name,
                          sample_property_cname,
                          norm_info,
                          nmr_norm,
                          backtransform) {

  # Set the elements of the norm_info list.
  norm_info$is_normalized <- nmr_norm
  norm_info$backtransform <- backtransform

  # Return the list of objects belonging to this attribute.
  return(list(metabolite_name = metabolite_name,
              sample_property_cname = sample_property_cname,
              norm_info = norm_info))

}

# Sets/updates the filters attribute with a filter class object.
#
# This function will create a filter class object. The output will always have
# the same elements but not all of them will be used for every filter type. This
# object will be appended to the list in the filters attribute for an omicsData
# object.
#
# @param filter_type
#
# @param threshold
#
# @param filtered
#
# @param filter_method
#
# @return A list containing all the elements in the filters attribute:
#         (list all attributes)
#
set_filter <- function (type,
                        threshold,
                        filtered,
                        method) {

  # Create an object that will have the filter elements and class added to it
  # later (in the next 10 lines or so).
  filta <- list()

  # Add the class of the filter being applied to filta.
  filta$type <- type

  # Add threshold to filta.
  filta$threshold <- threshold

  # Append filtered to the filta list.
  filta$filtered <- filtered

  # Affix the filter method to filta. This is only used with imdanovaFilt.
  filta$method <- method

  # Create the filter class.
  class(filta) = "filter"

  return (filta)

}

#' Set check.names attribute of omicsData object
#'
#' This function sets the check.names attribute of an omicsData object
#'
#' @param omicsData an object of the class 'pepData', 'proData', 'metabData',
#'   'lipidData', or 'nmrData', usually created by \code{\link{as.pepData}},
#'   \code{\link{as.proData}}, \code{\link{as.metabData}},
#'   \code{\link{as.lipidData}}, or \code{\link{as.nmrData}}, respectively.
#' @param set_to logical indicating what to set check.names attribute to.
#'   Defaults to TRUE.
#'
#' @return omicsData object with updated check.names attribute
#'
#' @rdname set_check_names
#'
#' @export
#'
set_check_names <- function (omicsData, set_to = TRUE) {

  # check that omicsData is of appropriate class #
  if(!inherits(omicsData, c("pepData", "proData", "metabData",
                            "lipidData", "nmrData")))
    stop (paste("omicsData must be of class 'pepData', 'proData',",
                "'metabData', 'lipidData', or 'nmrData'",
                sep = " "))

  #check that set_to is logical
  if(!is.logical(set_to)) stop ("set_to must be of class 'logical' ")

  attr(omicsData, "check.names") <- set_to

  return (omicsData)

}

# Functions to get statRes attributes ------------------------------------------

#' Return comparisons of statRes object
#'
#' This function returns comparisons from statRes or trellData object
#'
#' @param compObj is an object with the comparison attribute; specifically
#'   objects of class 'statRes' and 'trellData' objects derived from 'statRes'
#'   objects in \code{\link{format_data}}
#' @return returns a data frame with comparisons and their indices
#' @examples
#' \dontrun{
#' library(pmartR)
#' library(pmartRdata)
#'
#' my_prodata = group_designation(omicsData = pro_object,
#'                                main_effects = c("Condition"))
#'
#' imdanova_Filt = imdanova_filter(omicsData = my_prodata)
#'
#' my_prodata = applyFilt(filter_object = imdanova_Filt,
#'                        omicsData = my_prodata,
#'                        min_nonmiss_anova=2)
#'
#' imd_anova_res = imd_anova(omicsData = my_prodata,
#'                           test_method = 'comb',
#'                           pval_adjust='bon')
#'
#' result = get_comparisons(imd_anova_res)
#' }
#'
#' @rdname get_comparisons
#'
#' @export
#'
get_comparisons<- function(compObj){

  #check that compObj object is of 'statRes' or 'trellData' class
  if(!inherits(compObj, c("statRes", "trellData")))
    stop("object must be of class 'statRes' or 'trellData'")

  #check that compObj object is of 'statRes' or 'trellData' class
  if(inherits(compObj, "trellData") && is.null(attr(compObj, "comparisons"))) {

    return (NULL)

  } else {

    #pull comparisons attribute
    comp = attr(compObj, "comparisons")

    result = data.frame("comparisons" = as.character(comp),
                        "index" = 1:length(comp),
                        stringsAsFactors = FALSE)

    return(result)

  }

}

#' Return data_class of statRes or trellData object
#'
#' This function returns data_class attribute from statRes or trellData object,
#' inherited from the omicsData used in \code{\link{imd_anova}} or
#' \code{\link{format_data}}
#'
#' @param dcObj an object of class 'statRes' or 'trellData'
#'
#' @return returns the data_class attribute from a 'statRes' or 'trellData'
#'   object
#'
#' @examples
#' \dontrun{
#' library(pmartRdata)
#'
#' my_prodata = group_designation(omicsData = pro_object,
#'                                main_effects = c("Condition"))
#'
#' imdanova_Filt = imdanova_filter(omicsData = my_prodata)
#'
#' my_prodata = applyFilt(filter_object = imdanova_Filt,
#'                        omicsData = my_prodata,
#'                        min_nonmiss_anova=2)
#'
#' imd_anova_res = imd_anova(omicsData = my_prodata,
#'                           test_method = 'comb',
#'                           pval_adjust='bon')
#'
#' result = get_data_class(imd_anova_res)
#' }
#'
#' @rdname get_data_class
#'
#' @export
#'
get_data_class<- function (dcObj) {

  #check that compObj object is of 'statRes' class
  if(!inherits(dcObj, c("statRes", "trellData")))
    stop("dcObj object must be of class 'statRes' or 'trellData'")

  return(attr(dcObj, "data_class"))

}

#' Get group table
#'
#' This function returns a table with number of samples per group
#'
#' @param omicsObject an object of the class 'pepData', 'proData', 'metabData',
#'   'lipidData', 'nmrData', 'statRes', or 'trellData', usually created by
#'   \code{\link{as.pepData}}, \code{\link{as.proData}},
#'   \code{\link{as.metabData}}, \code{\link{as.lipidData}},
#'   \code{\link{as.nmrData}}, \code{\link{imd_anova}}, or
#'   \code{\link{format_data}} respectively.
#'
#' @return a table containing number of samples per group
#'
#' @rdname get_group_table
#'
#' @export
#'
get_group_table <- function (omicsObject) {

  # check that omicsObject is of appropriate class #
  if(!inherits(omicsObject, c("pepData", "proData", "metabData", "lipidData",
                              "nmrData", "statRes", "trellData")))
    stop (paste("omicsData must be of class 'pepData', 'proData', 'metabData',",
                "'lipidData', 'nmrData', 'statRes', or 'trellData'",
                sep = " "))

  if(is.null(attr(omicsObject, "group_DF"))) {

    return (NULL)

  } else {

    # should the result be constructed from omicsObject$f_data or from the
    # group_DF attr? if so...
    group = attr(omicsObject, "group_DF")$Group

    return (table(group))

  }

}

#'Helper to find the names of columns of a data.frame that contain exactly all
#'the elements of an input column
#'
#'@param df A data.frame whose columns we want to match to some query column.
#'@param col Vector of values which will be compared to a column in df.
#'
#'@return vector of column names of df that contain exactly all the elements of
#'the input column
#'
#'@keywords internal
column_matches_exact <- function(df, col) {
  diffs = lapply(df, function(df_col) {
    length(setdiff(
      union(df_col, col),
      intersect(df_col, col)
    ))
  })

  matched_cnames = names(diffs)[which(diffs == 0)]

  return(matched_cnames)
}

#' Custom message function to pretty-print text with newlines so you can follow
#' character limit guidelines in source code.
#' @noRd
wrap_message <- function(..., prefix = " ", initial = ""){
  message(strwrap(..., prefix = prefix, initial = initial))
}
