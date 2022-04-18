#' Creates Attribute of omicsData object containing Group Membership Based on
#' Specified Main Effects
#'
#' The method assigns each sample to a group, for future analyses, based on the
#' variable(s) specified as main effects.
#'
#' @param omicsData an object of the class 'lipidData', 'metabData', 'pepData',
#'   'proData', 'isobaricpepData', or 'nmrData', usually created by
#'   \code{\link{as.lipidData}}, \code{\link{as.metabData}},
#'   \code{\link{as.pepData}}, \code{\link{as.proData}},
#'   \code{\link{as.isobaricpepData}}, or \code{\link{as.nmrData}},
#'   respectively.
#' @param main_effects a character vector with no more than two variable names
#'   that should be used as main effects to determine group membership of
#'   samples. The variable name must match a column name from \code{f_data}.
#' @param covariates a character vector of no more than two variable names that
#'   should be used as covariates in downstream analyses. Covariates are
#'   typically variables that a user wants to account for in the analysis but
#'   quantifying/examining the effect of the variable is not of interest.
#' @param cov_type An optional character vector (must be the same length as
#'   \code{covariates} if used) indicating the class or type of each covariate.
#'   For example, "numeric", "character", or "factor". Partial matching ("num"
#'   for "numeric") is NOT used and the entire class/type must be typed out. If
#'   the class of a covariate does not match the input to \code{cov_type} the
#'   covariate will be coerced to that type. For example, if the covariate is a
#'   numeric vector of 0s and 1s (indicating two categories) and the input to
#'   \code{cov_type} is a class other than numeric this vector will be coerced
#'   to a character vector. The default value is NULL. In this case the class of
#'   the covariates is neither checked nor altered.
#' @param pairs A character string indicating the column in \code{f_data} that
#'   contains the pairing information. This string must match the column name
#'   exactly.
#' @param batch_id an optional character vector of no more than one variable that
#'  should be used as batch information for downstream analyses. Batch ID is
#'  similar to covariates but unlike covariates it is specific to that of specific
#'  batch effects
#' @param time_course an optional character string specifying the variable name
#'   of a time course variable, if applicable to the experimental data.
#'   CURRENTLY NOT SUPPORTED
#'
#' @details Groups are formed based on the levels of the main effect variables.
#'   One or two main effect variables are allowed. In the case of two main
#'   effect variables, groups are formed based on unique combinations of the
#'   levels of the two main effect variables. Any samples with level NA for a
#'   main effect variable will be removed from the data and will not be included
#'   in the final group designation results. Groups with a single sample are
#'   allowed, as is a single group.
#'
#' @return An object of the same class as the input \code{omicsData} object -
#'   the provided object with the samples filtered out, if any NAs were produced
#'   in designating groups. An attribute 'group_DF', a data.frame with columns
#'   for sample id and group, is added to the object. If two main effects are
#'   provided the original main effect levels for each sample are returned as
#'   the third and fourth columns of the data.frame. Additionally, the
#'   covariates provided will be listed as attributes of this data.frame. Use of
#'   time_course is currently not supported, however the eventual functionality
#'   is such that if time_course is included, a column for 'TimeCourse' will be
#'   output as well.
#'
#' @examples
#' \dontrun{
#' library(pmartRdata)
#' data(lipid_object)
#' lipid_object2 <- group_designation(omicsData = lipid_object,
#'                                    main_effects = "Condition")
#' attr(lipid_object2, "group_DF")
#' }
#'
#' @author Lisa Bramer, Kelly Stratton
#'
#' @export
#'
group_designation <- function (omicsData,
                               main_effects = NULL,
                               covariates = NULL,
                               cov_type = NULL,
                               pairs = NULL,
                               batch_id = NULL,
                               time_course = NULL) {

  # Initial checks on input arguments ------------------------------------------

  # check that omicsData is of appropriate class #
  if (!inherits(omicsData, c("pepData", "proData", "metabData",
                             "isobaricpepData", "lipidData", "nmrData"))) {

    # Throw an error that the input for omicsData is not the appropriate class.
    stop(paste("omicsData must be of class 'pepData', 'proData', 'metabData',",
               "'isobaricpepData', 'lipidData', or 'nmrData'",
               sep = ' '))

  }

  # check that isobaric data has been normalized #
  if (inherits(omicsData, "isobaricpepData") &&
      (attr(omicsData, "isobaric_info")$norm_info$is_normalized != TRUE)) {

    # Lay down an isobaric error for the user to trip over.
    stop (paste("omicsData with class 'isobaricpepData' must be normalized",
                "using normalize_isobaric(omicsData, apply_norm = TRUE) prior",
                "to calling group_designation()",
                sep = " "))

  }
  # no need for an analogous check for nmrData objects because there are not
  # entire samples that are used for normalization of NMR data

  # Check main_effects ---------------

  # A main effect does not need to be supplied if a pairing variable is.
  if (is.null(main_effects) && is.null(pairs)) {

    # A main effect must be specified unless a pairs variable is specified.
    stop (
      paste("A main effect must be specified unless the data are paired.",
            "In the case of paired data 'pairs' must be specified if",
            "there are no main effects.",
            sep = " ")
    )

  }

  # Check if any main effects have been provided.
  if (!is.null(main_effects)) {

    # Make sure the user does not use the forbidden main effect name. We must
    # defend the secrets of pmartR!
    if ("no_main_effects" %in% main_effects) {

      # Stop! You shall not pass!
      stop ("The name 'no_main_effects' cannot be used as a main effect name.")

    }

    # Check that main_effects is a character vector #
    if (!is.character(main_effects)) {

      # Stop production with a character vector error.
      stop ("main_effects must be a character vector.")

    }

    # Check that main_effects is of an appropriate length. The point is to be
    # like Goldilocks.
    if (length(main_effects) < 1) {

      # Error out with too few main effects.
      stop("No main effects were provided")

    }
    if (length(main_effects) > 2) {

      # Error out with too many main effects.
      stop ("No more than two main effects can be provided")

    }

    # Check that main_effects given are in f_data #
    if (sum(main_effects %in% names(omicsData$f_data)) != length(main_effects)) {

      stop("One or more of the main_effects is not found in f_data of omicsData")

    }

  }

  # Check covariates ---------------

  # See if covariates are present.
  if (!is.null(covariates)) {

    # Check that covariates is a character vector #
    if (!is.character(covariates)) {

      # Stop production with a character vector error.
      stop ("covariates must be a character vector.")

    }

    # Check that covariates is an appropriate length. Think Goldilocks! #
    if (length(covariates) > 2) {

      # Error out with too many covariates.
      stop ("No more than two covariates can be provided.")

    }

    # Check that covariates given are in f_data #
    if (sum(covariates %in% names(omicsData$f_data)) != length(covariates)) {

      stop("One or more of the covariates is not found in f_data of omicsData")

    }

    # Check if the covariate types are specified.
    if (!is.null(cov_type)) {

      # Make sure the length of covariates and cov_type are the same.
      if (length(covariates) != length(cov_type)) {

        stop ("The length of covariates and cov_type must be the same.")

      }

      # Loop through all covariates and compare their class to the cov_type
      # input. If the classes do not match the covariate will be converted to a
      # character vector.
      for (e in 1:length(covariates)) {

        # Nab the class of the eth covariate.
        da_class <- class(omicsData$f_data[[covariates[[e]]]])

        # Check class and convert to character if the class does not match the
        # covariate type. We cannot simply check if da_class and cov_type are
        # the same and convert the covariate to a character vector if they are
        # not. This is an insufficient check because the type of vector could be
        # "double" but the cov_type is "numeric". In this case the covariate
        # would incorrectly be converted to a character vector.
        if (da_class %in% c("numeric", "integer", "double") &&
            !(cov_type[[e]] %in% c("numeric", "integer", "double"))) {

          # Convert the covariate to a character vector.
          omicsData$f_data[[covariates[[e]]]] <- as.character(
            omicsData$f_data[[covariates[[e]]]]
          )

        } else if (!(da_class %in% c("numeric", "integer", "double")) &&
                   da_class != cov_type[[e]]) {

          # Convert the covariate to a character vector.
          omicsData$f_data[[covariates[[e]]]] <- as.character(
            omicsData$f_data[[covariates[[e]]]]
          )

        }

      }

    }

  }

  # Check batch_id ---------------

  # See if batch_id is present
  if (!is.null(batch_id)) {

    # Check that batch_id is a character vector #
    if (!is.character(batch_id)) {

      # Stop production with a character vector error
      stop ("batch_id must be a character vector.")
    }

    # Check that batch_id is an appropriate length
    if (length(batch_id) > 1){

      # Error out with too many batch_ids
      stop ("No more than one batch_id can be provided.")
    }

    # Check that batch_id is given in f_data #
    if (!(batch_id %in% names(omicsData$f_data))) {

      stop ("batch_id is not found in f_data of omicsData")
    }


  }

  # Check pairs ---------------

  # Have a looksie at the pairs argument. If it is present put it through the
  # usual methods of information extraction and verification.
  if (!is.null(pairs)) {

    # Check that pairs is a character vector.
    if (!is.character(pairs)) {

      # Holy inappropriate input type, Batman.
      stop ("pairs must be a character vector.")

    }

    # Make sure there is only one character string specified.
    if (length(pairs) > 1) {

      # Holy too many pairs, Batman!
      stop ("Only one paired variable can be specified.")

    }

    # Make sure the paired variable exists in f_data. How can we subset by
    # something that doesn't exist?!
    if (!(pairs %in% names(omicsData$f_data))) {

      # Holy missing variable, Batman!
      stop ("The variable specified for pairs does not exist in f_data.")

    }

    # Count the number of pair IDs that do not have two entries in f_data.
    not_two <- which(unname(table(omicsData$f_data[, pairs])) < 2)

    # Ensure each pair has at least two observations in f_data.
    if (length(not_two) > 0) {

      pair_id <- names(table(omicsData$f_data[, pairs]))[not_two]

      stop (paste("The following pairs do not have at least two samples to ",
                  "form a pair: ",
                  knitr::combine_words(pair_id),
                  ".",
                  sep = ""))

    }

    # Check that the main effect(s) are the same for each pair.
    if (!is.null(main_effects)) {

      # If there are two main effects create a new variable where the main
      # effects are combined with a "_" between them.
      if (length(main_effects) == 2) {

        # Count the number of unique combined main effect values for each pair.
        # The format for combined main effects is me1_me2.
        reprobates <- omicsData$f_data %>%
          dplyr::rowwise() %>%
          dplyr::mutate(
            both_me = paste(!!rlang::sym(main_effects[[1]]),
                            !!rlang::sym(main_effects[[2]]),
                            sep = "_")
          ) %>%
          dplyr::ungroup() %>%
          dplyr::group_by(!!rlang::sym(pairs)) %>%
          dplyr::mutate(n_me = dplyr::n_distinct(both_me)) %>%
          dplyr::ungroup() %>%
          dplyr::filter(n_me > 1) %>%
          dplyr::pull(get_fdata_cname(omicsData))

        # Use the main effect variable directly (only one main effect exists).
      } else {

        # Count the number of unique main effect values for each pair.
        reprobates <- omicsData$f_data %>%
          dplyr::group_by(!!rlang::sym(pairs)) %>%
          dplyr::mutate(
            n_me = dplyr::n_distinct(!!rlang::sym(main_effects))
          ) %>%
          dplyr::ungroup() %>%
          dplyr::filter(n_me > 1) %>%
          dplyr::pull(get_fdata_cname(omicsData))

      }

      # Check if some main effects differ between pairs.
      if (length(reprobates) > 0) {

        # Let them have it for making my life miserable.
        stop (paste("The following samples have main effects that differ",
                    "between pairs:",
                    knitr::combine_words(reprobates),
                    sep = " "))

      }

    }

    # Check that the covariate(s) are the same for each pair.
    if (!is.null(covariates)) {

      # If there are two covariates create a new variable where the covariates
      # are combined with a "_" between them.
      if (length(covariates) == 2) {

        # Count the number of unique combined covariate values for each pair.
        # The format for combined covariates is cov1_cov2.
        reprobates <- omicsData$f_data %>%
          dplyr::rowwise() %>%
          dplyr::mutate(
            both_cov = paste(!!rlang::sym(covariates[[1]]),
                             !!rlang::sym(covariates[[2]]),
                             sep = "_")
          ) %>%
          dplyr::ungroup() %>%
          dplyr::group_by(!!rlang::sym(pairs)) %>%
          dplyr::mutate(n_cov = dplyr::n_distinct(both_cov)) %>%
          dplyr::ungroup() %>%
          dplyr::filter(n_cov > 1) %>%
          dplyr::pull(get_fdata_cname(omicsData))

        # Use the covariate variable directly (only one covariate exists).
      } else {

        # Count the number of unique covariate values for each pair.
        reprobates <- omicsData$f_data %>%
          dplyr::group_by(!!rlang::sym(pairs)) %>%
          dplyr::mutate(
            n_cov = dplyr::n_distinct(!!rlang::sym(covariates))
          ) %>%
          dplyr::ungroup() %>%
          dplyr::filter(n_cov > 1) %>%
          dplyr::pull(get_fdata_cname(omicsData))

      }

      # Check if some covariates differ between pairs.
      if (length(reprobates) > 0) {

        # Let them have it for making my life miserable.
        stop (paste("The following samples have covariates that differ",
                    "between pairs:",
                    knitr::combine_words(reprobates),
                    sep = " "))

      }

    }

  }

  # Check time_course ---------------

  # See if time_course is present.
  if (!is.null(time_course)) {

    # added Aug 2020, since use of this argument is currently not actively
    # supported #
    if (is.na(time_course)) {time_course = NULL}

    # Give an error and don't run the rest of the function because time_course
    # is not currently supported.
    stop (paste("Use of the time_course argument is currently not supported.",
                "In group_designation set time_course = NULL.",
                sep = " "))

    if (!is.character(time_course)) {
      # check that time_course is a character #
      stop (paste("time_course must be a character specifying the name of the",
                  "variable denoting time",
                  sep = " "))
    }

    # check that time_course is found in the data #
    if(!(time_course %in% names(omicsData$f_data))) {

      stop ("time_course is not found in f_data of omicsData")

    }

    # check that no more than 1 main_effect is specified when time_course is non
    # NULL #
    if(length(main_effects) > 1) {

      stop ("Only 1 main effect may be specified when time_course is provided")

    }

  }

  # Create the group data frame and attributes ---------------------------------

  # If no main effect was provided but the data are paired create a substitute
  # main effects variable with just one level.
  if (is.null(main_effects) && !is.null(pairs)) {

    # Add a main effect name. This will be used as a place holder because later
    # in this function and downstream functions expect a main effect variable to
    # be present.
    main_effects <- "no_main_effect"

    # Create a main effect with just one level. I tried to select something that
    # won't show up as an actual main effect level in someone's study. However,
    # this rarely goes well and someone at some point will come crying to us
    # because pmartR is giving them crazy output and/or errors.
    omicsData$f_data$no_main_effect <- "paired_diff"

  }

  # pull sample id column name #
  samp_id <- get_fdata_cname(omicsData)

  # Get number of main effects #
  n.maineffects = length(main_effects)

  # Case 1: 1 main effect variable #
  if (n.maineffects == 1) {

    # create output formatted with first column being sample id and second
    # column group id #
    output <- data.frame(Sample.ID = as.character(omicsData$f_data[, samp_id]),
                         Group = as.character(omicsData$f_data[, main_effects]),
                         stringsAsFactors = FALSE)
    names(output)[1] = samp_id

    # Case 2: 2 main effect variables #
  } else if (n.maineffects == 2) {

    # create a group variable and paste main effect levels together for samples
    # # samples with a value of NA for either main effect will have a Group
    # value of NA #
    Group = rep(NA, nrow(omicsData$f_data))

    # identify samples that will have a Group membership that is not missing #
    nonna.group <- (!is.na(omicsData$f_data[, main_effects[[1]]]) &
                      !is.na(omicsData$f_data[, main_effects[[2]]]))

    # Combine names of the first and second main effects to create a group
    # variable.
    Group[nonna.group] <- paste(
      omicsData$f_data[nonna.group, main_effects[[1]]],
      omicsData$f_data[nonna.group, main_effects[[2]]],
      sep = "_"
    )

    # create output formatted with first column being sample id and second
    # column group id # third and fourth columns are the original main effect
    # levels #
    output <- data.frame(
      Sample.ID = as.character(omicsData$f_data[, samp_id]),
      Group = Group,
      me1 = as.character(omicsData$f_data[, main_effects[[1]]]),
      me2 = as.character(omicsData$f_data[, main_effects[[2]]]),
      stringsAsFactors = FALSE
    )

    names(output) = c(samp_id, 'Group', main_effects[[1]], main_effects[[2]])
  }

  # if time_course is non NULL add a column to the output #
  ### commented out Aug 2020 ###
  # if(!is.null(time_course)){
  #   output$TimeCourse = temp_data[, which(colnames(temp_data) == time_course)]
  # }

  # check for group sizes <2; if there are any, set them to NA's # CHANGED BELOW

  # Check if any of the groups have NAs.
  if (any(is.na(output$Group))) {
    ### COMMENTED OUT AUG 2020 TO UPDATE FUNCTIONALITY TO ALLOW FOR GROUPS OF
    # SIZE 1 output$Group = as.character(output$Group) sm.groups =
    # which(table(output$Group)<2) sm.groups = names(sm.groups) n.sm.groups =
    # length(sm.groups) if(n.sm.groups>0){ for(i in 1:n.sm.groups){
    # output$Group[output$Group==sm.groups[i]] = NA } }

    # output a warning message listing the samples that are removed (7/8/2016
    # KS; updated Aug 2020 by KS)
    na.group <- as.character(output[is.na(output$Group), samp_id])
    n.samps = length(na.group)
    mystr <- paste("The following ",
                   n.samps,
                   " sample(s) has/have been removed from the dataset due to",
                   " missing group information: ",
                   sep = "")
    mystr2 <- paste(as.character(na.group), sep="' '", collapse=", ")
    warning(paste0(mystr, mystr2))
    # note: this doesn't actually remove them from the dataset, that is done
    # further below

  } else {

    # If no group samples have NA values set na.group to null. This will allow
    # all rows in the group data frame (output) to be selected.
    na.group <- NULL

  }

  # Set na.time to null because time course data is not supported.
  na.time <- NULL

  # check for time courses with <2 samples; if there are any, set them to NA's
  ### commented out Aug 2020 ###
  # if(!is.null(time_course)){
  #   sm.groups = which(table(output$TimeCourse)<2)
  #   sm.groups = names(sm.groups)
  #   n.sm.groups = length(sm.groups)
  #   if(n.sm.groups>0){
  #     for(i in 1:n.sm.groups){
  #       output$TimeCourse[output$TimeCourse==sm.groups[i]] = NA
  #     }
  #   }

  # output a warning message listing the groups and samples that are removed
  # (7/8/2016 KS) ## COMMENTED OUT AUG 2020 ### sm.samps =
  # as.character(output[is.na(output$TimeCourse), samp_id]) n.samps =
  # length(sm.samps) mystr <- paste("The following ", n.sm.groups, " time
  # courses have been removed from the dataset due to a time course size of less
  # than 2 samples: ", sep="") mystr2 <- paste(as.character(sm.groups), sep="'
  # '", collapse=", ") mystr3 <- paste(". This corresponds to the following ",
  # n.samps, " samples: ", sep="") mystr4 <- paste(as.character(sm.samps),
  # sep="' '", collapse=", ") warning(paste(mystr, mystr2, mystr3, mystr4,
  # sep=""))

  # find any samples with NA for time course #
  ### commented out Aug 2020 ###
  # if(!is.null(time_course)){
  #   na.time = as.character(output[which(is.na(output$TimeCourse)), samp_id])
  # }else{
  # na.time = NULL
  # }

  # Select all rows that have an NA for either group or time.
  all.na = union(na.group, na.time)

  # remove samples with group or time NA from output and data
  if (length(all.na) > 0) {

    # remove from output #
    output <- output[-which(output[, samp_id] %in% all.na),]

    # Find rows or columns that need to be removed from e_data and f_data.
    edat_ids = which(names(omicsData$e_data) %in% all.na)
    fdat_ids = which(omicsData$f_data[, samp_id] %in% all.na)

    # Remove columns from e_data that have groups with NAs.
    omicsData$e_data = omicsData$e_data[, -edat_ids]

    # Remove rows from f_data that have groups with NAs.
    omicsData$f_data = omicsData$f_data[-fdat_ids, ]

  }

  # Add attributes to the group data frame for the main effects, covariates,
  # time course, and non-singleton groups.
  attr(output, "main_effects") = main_effects

  # Set the covariates attribute according to the input (very complicated).
  if (is.null(covariates)) {

    holy_covariates_batman <- NULL

    # Change the covariates attribute to a data frame. Later this will be
    # updated to make sure the covaraite columns are in the correct format, but
    # for now they will just be added to the covariates data frame.
  } else if (length(covariates) == 1) {

    # Make the data frame with the sample ID as the first column and the one and
    # only covariate as the second column.
    holy_covariates_batman <- data.frame(
      sample_id = as.character(omicsData$f_data[, samp_id]),
      cov1 = omicsData$f_data[, covariates],
      stringsAsFactors = FALSE
    )

    # Rename the columns to match the names in f_data.
    names(holy_covariates_batman) <- c(samp_id, covariates)

    # Runs when there are just two covariates.
  } else {

    # Make the data frame with the sample ID as the first column and the second
    # and third columns as the covariates.
    holy_covariates_batman <- data.frame(
      sample_id = as.character(omicsData$f_data[, samp_id]),
      cov1 = omicsData$f_data[, covariates[[1]]],
      cov2 = omicsData$f_data[, covariates[[2]]],
      stringsAsFactors = FALSE
    )

    # Rename the columns to match the names in f_data.
    names(holy_covariates_batman) <- c(samp_id,
                                       covariates[[1]],
                                       covariates[[2]])

  }

  attr(output, "covariates") <- holy_covariates_batman

  # Include the name of the variable containing the paired information as an
  # attribute of group_DF. This will be used in the functions that compute the
  # difference between each pair.
  attr(output, "pairs") <- pairs

  # Set the batch_id attribute according to the input
  if (is.null(batch_id)) {

    holy_batch_robin <- NULL
  } else {

    # make the data frame with the Sample ID as the first column and the batch id
    # as the second column
    holy_batch_robin <- data.frame(
      sample_id = as.character(omicsData$f_data[, samp_id]),
      batch = omicsData$f_data[, batch_id],
      stringsAsFactors = FALSE
    )
    # rename columns to match the names in f_data
    names(holy_batch_robin) <- c(samp_id, batch_id)
  }

  attr(output, "batch_id") <- holy_batch_robin

  ### changed to NA Aug 2020 ###
  attr(output, "time_course") = NULL #time_course
  ## added attribute that lists groups with 2 or more samples Nov 2020 ##
  # Find groups with more than one sample in them.
  nonsingleton_groups <- names(which(table(output$Group) > 1))
  attr(output, "nonsingleton_groups") <- nonsingleton_groups

  # Add the group information to the group_DF attribute in the omicsData object.
  attr(omicsData, "group_DF") = output

  # Update the data_info attribute.
  attr(omicsData, 'data_info') <- set_data_info(
    e_data = omicsData$e_data,
    edata_cname = get_edata_cname(omicsData),
    data_scale_orig = get_data_scale_orig(omicsData),
    data_scale = get_data_scale(omicsData),
    data_types = get_data_info(omicsData)$data_types,
    norm_info = get_data_info(omicsData)$norm_info,
    is_normalized = get_data_info(omicsData)$norm_info$is_normalized,
    batch_info = get_data_info(omicsData)$batch_info,
    is_bc = get_data_info(omicsData)$batch_info$is_bc
  )

  # Update the meta_info attribute.
  attr(omicsData, 'meta_info') <- set_meta_info(
    e_meta = omicsData$e_meta,
    emeta_cname = get_emeta_cname(omicsData)
  )

  # Return the updated and polished omicsData object!!!
  return (omicsData)

}
