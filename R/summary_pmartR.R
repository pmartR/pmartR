#' Produce a basic summary of a pmartR omicsData S3 Object
#'
#' This function will provide basic summary statistics for omicsData objects
#' from the pmartR package.
#'
#' @param omicsData an object of the class 'lipidData', 'metabData', 'pepData',
#'   'proData', or 'nmrData' usually created by \code{\link{as.lipidData}},
#'   \code{\link{as.metabData}}, \code{\link{as.pepData}},
#'   \code{\link{as.proData}}, or \code{\link{as.nmrData}}, respectively.
#'
#' @return a summary table for the pmartR omicsData object. If assigned to a
#'   variable, the elements of the summary table are saved in a list format.
#'
#' @examples
#' \dontrun{
#' library(pmartRdata)
#' data("pep_object")
#' summary(pep_object)
#' }
#'
#' @author Lisa Bramer, Kelly Stratton, Thomas Johansen
#'
#' @export
#'
#' @rdname summary-omicsData
#' @name summary-omicsData
#'
summary.pepData <- function(omicsData) {

  # Army captain, "Nooooo! ... General, the enemy is about to overtake us!"
  # Army general, "CALL IN THE SUMMARIZER!!!"
  summarizer(omicsData)

}

#'@export
#'@rdname summary-omicsData
#'@name summary-omicsData
summary.proData <- function(omicsData) {

  # Army captain, "Nooooo! ... General, the enemy is about to overtake us!"
  # Army general, "CALL IN THE SUMMARIZER!!!"
  summarizer(omicsData)

}


#'@export
#'@rdname summary-omicsData
#'@name summary-omicsData
summary.lipidData <- function(omicsData) {

  # Army captain, "Nooooo! ... General, the enemy is about to overtake us!"
  # Army general, "CALL IN THE SUMMARIZER!!!"
  summarizer(omicsData)

}



#'@export
#'@rdname summary-omicsData
#'@name summary-omicsData
summary.metabData <- function(omicsData) {

  # Army captain, "Nooooo! ... General, the enemy is about to overtake us!"
  # Army general, "CALL IN THE SUMMARIZER!!!"
  summarizer(omicsData)

}


#'@export
#'@rdname summary-omicsData
#'@name summary-omicsData
summary.nmrData <- function(omicsData) {

  # Army captain, "Nooooo! ... General, the enemy is about to overtake us!"
  # Army general, "CALL IN THE SUMMARIZER!!!"
  summarizer(omicsData)

}

# A Neat little function that summarizes an omicsData object. This function adds
# covariate and pair counts to the original summary function.
#
# @author Evan A Martin
summarizer <- function (omicsData) {

  # The following values will always be present no matter what. They wont change
  # depending on the main effects, covariates, and whether the data are paired.
  # For the remainder of the script we will append values to this list.
  res <- list(
    class = class(omicsData),
    num_samps = attr(omicsData, "data_info")$num_samps,
    num_edata = attr(omicsData, "data_info")$num_edata,
    num_emeta = attr(omicsData, "meta_info")$num_emeta,
    num_miss_obs = attr(omicsData, "data_info")$num_miss_obs,
    prop_missing = round(attr(omicsData, "data_info")$prop_missing, 3)
  )

  # Create a row names vector that matches the elements of res. Because res will
  # always have the elements above the row names vector will also have a name
  # for each of those elements. This will be the base vector. Additional
  # elements will be appended to the row names vector depending on main effects,
  # covariates, and whether the data are paired.
  res_names <- c(
    "Class",
    paste("Unique ", get_edata_cname(omicsData), "s (e_data)", sep = ""),
    paste("Unique ", get_fdata_cname(omicsData), "s (f_data)", sep = ""),
    if (is.null(get_emeta_cname(omicsData)))
      "Rows (e_meta)" else
        paste("Unique ", get_emeta_cname(omicsData), "s (e_meta)", sep = ""),
    "Missing Observations",
    "Proportion Missing"
  )

  # Check if the data are grouped.
  if (!is.null(attr(omicsData, "group_DF"))) {

    # Add group counts ---------------

    group_vec <- attr(omicsData, "group_DF")$Group
    levels <- unique(group_vec)
    counts <- vector(mode = "list", length = length(levels))

    for (i in 1:length(levels)) {

      # Determine the number of samples in the current group.
      counts[[i]] <- length(which(group_vec == levels[[i]]))

    }

    # Name the counts with their corresponding group names.
    names(counts) <- levels

    # Add the group counts to res.
    res <- c(res, counts)

    # Update the names of res depening on the main effects of the groups. For
    # example, if the data are paired and there are no main effects we don't
    # want to print the number of samples in the paired_diff group. (We want to
    # keep this group a secret from the user. Shhhhhhh!)
    if ("paired_diff" %in% levels) {

      res_names <- c(res_names, "Total Samples")

    } else {

      # Add each group name to the res_names vector.
      res_names <- c(res_names, paste("Samples per group:", levels, sep = " "))

    }

    # Add covariate counts ---------------

    # Check for covariates. If they exist we will tell the user how many.
    if (!is.null(attr(get_group_DF(omicsData), "covariates"))) {

      # Nab the covariate data frame. This will be used later for sundry
      # purposes.
      covies <- attr(get_group_DF(omicsData), "covariates")

      # Determine the number of covariates. This is the number of columns in the
      # covariates attribute minus one (because of the sample name column).
      n_cov <- dim(covies)[[2]] - 1

      # Loop through each covariate, check its class (e.g., numeric, character,
      # factor), and count the number of samples in each level.
      for (e in 1:n_cov) {

        # Check if the covariate is categorical (i.e., the class is character or
        # factor).
        if (is.character(covies[, e + 1]) || is.factor(covies[, e + 1])) {

          # Count the number of levels in the current covariate and create a
          # list to hold the number of samples in each level.
          cov_level <- as.character(unique(covies[, e + 1]))
          level_count <- vector(mode = "list", length = length(cov_level))

          # Loop through each level and count the number of samples.
          for (v in 1:length(cov_level)) {

            # Determine the number of samples in the current level.
            level_count[[v]] <- length(which(covies[, e + 1] == cov_level[[v]]))

          }

          # Name the counts with their corresponding covariate/level names.
          names(level_count) <- paste(names(covies)[e + 1],
                                      cov_level,
                                      sep = "_")

          # Runs when the covariate is numeric.
        } else {

          # The covariate is continuous: We will just report the name of the
          # covaraite along with the total number of samples for this covariate.
          level_count <- list(nrow(covies))
          names(level_count) <- names(covies)[e + 1]

        }

        # Add the covariate counts to res.
        res <- c(res, level_count)

        # Add the name of the covariate levels to res_names.
        res_names <- c(res_names,
                       paste("Samples per covariate:",
                             names(level_count), sep = " "))

      }

    }

    # Add pair counts ------------

    # Check if the data are paired. If they are the number of pairs will be
    # added to the beautiful summary we are creating!
    if (!is.null(attr(get_group_DF(omicsData), "pairs"))) {

      # Nab the name of the pair variable. This will be used to subset f_data.
      pair_name <- attr(get_group_DF(omicsData), "pairs")

      # Add the number of pairs to res.
      res <- c(res,
               list(pairs = dplyr::n_distinct(omicsData$f_data[, pair_name])))

      # Update the names vector with the row name for the number of pairs.
      res_names <- c(res_names, "Pairs")

    }

  }

  # Convert all elements in res to a character.
  res <- lapply(res, function(x) ifelse(is.null(x), "NA", as.character(x)))

  # Turn res into a data frame, remove all column names, and include the
  # established row names.
  catmat <- data.frame(unlist(res, use.names = FALSE))
  colnames(catmat) <- NULL
  rownames(catmat) <- res_names

  return (catmat)

}
