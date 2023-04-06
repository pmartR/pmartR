#' Molecule Filter Summary
#'
#' Provide summary of a moleculeFilt S3 object
#'
#' @param object S3 object of class 'moleculeFilt' created by
#'   \code{\link{molecule_filter}}
#' @param min_num integer value specifying the minimum number of times each
#'   feature must be observed across all samples. Default value is NULL.
#' @return a summary table giving the number of biomolecules by number of
#'   observed values across all samples. If min_num is specified, the numbers of
#'   biomolecules to be filtered and to be retained based on the specified
#'   threshold are reported. If, upon creation of moleculeFilt object,
#'   \code{use_groups = TRUE} or \code{use_batches = TRUE} were specified, the
#'   numbers reported by the summary are based on groups and/or batches.
#'
#' @examples
#' library(pmartRdata)
#' myfilter <- molecule_filter(omicsData = pep_object)
#' summary(myfilter)
#' summary(myfilter, min_num = 2)
#'
#' @seealso \code{\link{molecule_filter}}
#'
#' @author Lisa Bramer, Kelly Stratton
#'
#' @export
#' @rdname summary.moleculeFilt
#' @name summary.moleculeFilt

summary.moleculeFilt <- function(object, min_num = NULL, ...) {
  filter_object <- object

  if (!is.null(min_num)) {
    # check that min_num is not a vector #
    if (length(min_num) > 1) stop("min_num must be of length 1")
    # check that min_num is numeric >= 0 #
    if (!is.numeric(min_num) | min_num < 0) stop("min_num must be an integer >= 0")
    # check that min_num is an integer #
    if (min_num %% 1 != 0) stop("min_num must be an integer >= 0")
    # check that min_num is less than the max number of observations #
    if (min_num > max(filter_object$Num_Observations)) stop("min_num cannot be greater than the number of samples")
  }

  # return the numeric version of plot, the threshold used, the number that would be tested and the number that would not be tested

  # how many peptides appear in the dataset once, twice, 3 times, etc.
  cut_data <- table(cut(filter_object$Num_Observations, breaks = -1:max(filter_object$Num_Observations)))
  cumcounts <- cumsum(cut_data)
  pep_observation_counts <- data.frame(num_observations = 0:(length(cumcounts) - 1), frequency_counts = cumcounts)

  if (!is.null(min_num)) {
    # get number molecules tested
    num_not_tested <- pep_observation_counts$frequency_counts[pep_observation_counts$num_observations == (min_num - 1)]

    # get number molecules not tested
    num_tested <- pep_observation_counts$frequency_counts[nrow(pep_observation_counts)] - num_not_tested
  } else {
    num_tested = NULL
    num_not_tested = NULL
    min_num = NULL
  }

  res <- list(pep_observation_counts = pep_observation_counts[-1, ], min_num = min_num, num_not_filtered = num_tested, num_filtered = num_not_tested)

  class(res) = c("moleculeFilterSummary", "list")

  attr(res, "use_batch") <- attr(filter_object, "use_batch")
  attr(res, "use_groups") <- attr(filter_object, "use_groups")

  return(res)
}



#' Molecule Filter Print Method
#'
#' Print method for moleculeFilt S3 object
#' @param x the moleculeFilt summary to print
#' @param ... further arguments passed to or from other methods
#'
#' @export
#' @name print.moleculeFilterSummary

print.moleculeFilterSummary <- function(x, ...) {
  object <- x

  # create output #
  cat("\nSummary of Molecule Filter \n----------------------------------\n")
  if (attr(object, "use_batch") & !attr(object, "use_groups")) {
    cat("\nMinimum Number is Minimum Number per Batch \n----------------------------------\n")
  }
  if (!attr(object, "use_batch") & attr(object, "use_groups")) {
    cat("\nMinimum Number is Minimum Number per Group \n----------------------------------\n")
  }
  if (attr(object, "use_batch") & attr(object, "use_groups")) {
    cat("\nMinimum Number is Minimum Number per Batch and Group \n----------------------------------\n")
  }

  if (!is.null(object$min_num)) {
    cat(sprintf("Minimum Number:%d\nFiltered:%d\nNot Filtered:%d\n\n", object$min_num, object$num_filtered, object$num_not_filtered))
  }
  print(object$pep_observation_counts, row.names = FALSE)
}

#' RNA Filter Summary
#'
#' Provide summary of a RNAFilt S3 object
#'
#' @param object S3 object of class 'RNAFilt' created by
#'   \code{\link{RNA_filter}}.
#' @param min_nonzero integer or float between 0 and 1. Cut-off for number of
#'   unique biomolecules with non-zero counts or as a proportion of total
#'   biomolecules. Defaults to NULL.
#' @param size_library integer cut-off for sample library size (i.e. number of
#'   reads). Defaults to NULL.
#'
#' @return a summary table giving the minimum, maximum, 1st and 3rd quartiles,
#'   mean and standard deviation for library size (the number of unique
#'   biomolecules with non-zero observations per sample), and the proportion of
#'   non-zero observations over the total number of biomolecules.
#'
#' @examples
#' library(pmartRdata)
#' myfilter <- RNA_filter(omicsData = rnaseq_object)
#' summary(myfilter)
#' summary(myfilter, min_nonzero = 2)
#'
#' @author Rachel Richardson
#'
#' @seealso \code{\link{RNA_filter}}
#'
#' @export
#' @rdname summary.RNAFilt
#' @name summary.RNAFilt

summary.RNAFilt <- function(object,
                            size_library = NULL,
                            min_nonzero = NULL,
                            ...) {
  filter_object <- object

  # Have a looksie at the class of the filter object.
  if (!inherits(filter_object, "RNAFilt")) {
    # Fezzik, tear his arms off.
    stop("filter_object must be of class RNAFilt")
  }

  ## Checks for size_library as single integer
  if (!is.null(size_library) &&
    (length(size_library) > 1 ||
      size_library %% 1 != 0 ||
      size_library > max(filter_object$LibrarySize)
    )
  ) stop(
    paste0(
      "size_library must be integer of length 1 less than max library size (",
      max(filter_object$LibrarySize),
      ")"
    )
  )

  ## Checks for min_nonzero as single numeric
  if (!is.null(min_nonzero)) {
    ## Length
    if (length(min_nonzero) > 1) stop("min_nonzero must be length 1")

    ## proportion or int
    if (!(min_nonzero %% 1 == 0 || (min_nonzero > 0 && min_nonzero < 1))) stop(
      "min_nonzero must be integer or numeric between 0 and 1."
    )

    ## Within appropriate bounds
    if (min_nonzero %% 1 == 0 && min_nonzero > max(filter_object$NonZero)) stop(
      paste0(
        "min_nonzero exceeds maximum number of non-zero biomolecules (",
        max(filter_object$NonZero),
        ")"
      )
    )

    if (min_nonzero %% 1 != 0 &&
      min_nonzero > max(filter_object$ProportionNonZero)) stop(
      paste0(
        "min_nonzero exceeds maximum proportion of non-zero biomolecules (",
        signif(max(filter_object$ProportionNonZero)),
        ")"
      )
    )
  }

  #---------------
  temp_obj <- filter_object
  if (!is.null(min_nonzero)) {
    column_use <- ifelse(min_nonzero %% 1 == 0,
      "NonZero", "ProportionNonZero"
    )
    temp_obj <- temp_obj[temp_obj[[column_use]] >= min_nonzero, ]
  }

  if (!is.null(size_library)) {
    temp_obj <- temp_obj[temp_obj[["LibrarySize"]] >= size_library, ]
  }

  ## Pre-filter
  df <- as.data.frame(apply(temp_obj[-1], 2, summary))
  df["SD", ] <- apply(temp_obj[-1], 2, sd)

  rmved <- !(filter_object[[colnames(filter_object)[[1]]]] %in% temp_obj[[colnames(filter_object)[[1]]]])

  res <- list(
    Summary = df,
    samples_filtered = filter_object[[colnames(filter_object)[[1]]]][rmved],
    num_filtered = sum(rmved),
    num_not_filtered = sum(!rmved),
    nonzero_thresh = min_nonzero,
    library_thresh = size_library
  )

  class(res) = c("RNAFiltSummary", class(res))

  return(res)
}


#' RNA Filter Print Method
#'
#' Print method for summary of RNAFilt
#' @param x the RNAFilt summary to print
#' @param ... further arguments passed to or from other methods
#'
#' @export
#' @name print.RNAFiltSummary
print.RNAFiltSummary <- function(x, ...) {
  object <- x

  # create output #
  cat("\nSummary of RNA Filter \n----------------------------------\n")

  if (!is.null(object$library_thresh)) {
    cat(sprintf(
      "Library size cut-off: %d\n\n",
      object$library_thresh
    ))
  }

  if (!is.null(object$nonzero_thresh)) {
    if (object$nonzero_thresh %% 1 == 0) {
      cat(sprintf(
        "N Non-zero biomolecule cut-off: %d\n\n",
        object$nonzero_thresh
      ))
    } else {
      cat(sprintf(
        "Proportion Non-zero biomolecule cut-off: %f\n\n",
        object$nonzero_thresh
      ))
    }
  }

  if (length(object$samples_filtered) > 0) {
    cat(sprintf(
      "Samples Filtered:%d\n Samples Not Filtered:%d\n\n",
      object$num_filtered,
      object$num_not_filtered
    ))
  }

  print(object$Summary)
}


#' Total Count Filter Summary
#'
#' Provide summary of a totalCountFilt S3 object
#'
#' @param object S3 object of class 'totalCountFilt' created by
#' \code{\link{total_count_filter}}.
#' @param min_count numeric value greater than 1 and less than the value
#' given by filter_object$Total_Count. Values below min_count are filtered out.
#' Default value is NULL.
#' @return a summary of the Total Count values, number of zero values, and
#' non-zero values. If a min_count is provided the biomolecules that would be
#' filtered at this threshold are reported.
#'
#' @examples
#' \dontrun{
#' library(pmartRdata)
#' myfilt <- total_count_filter(omicsData = rnaseq_object)
#' summary(myfilt, min_count = 15)
#' }
#'
#' @author Rachel Richardson
#'
#' @seealso \code{\link{total_count_filter}}
#'
#' @export
#' @rdname summary.totalCountFilt
#' @name summary.totalCountFilt
summary.totalCountFilt <- function(object, min_count = NULL, ...) {
  filter_object <- object

  # checks for cv_threshold if not null #
  if (!is.null(min_count)) {
    # check that cv_threshold is numeric
    if (!is.numeric(min_count)) stop("cv_threshold must be numeric of length 1")
    # chack that cv_threshold is of length 1
    if (length(min_count) > 1) stop("cv_threshold must be numeric of length 1")
    # check that cv_threshold is more than 1 and less than max CV value
    if (min_count <= 1 | min_count >= max(filter_object$Total_Counts, na.rm = TRUE)) stop("cv_threshold must be greater than 1 and less than the maximum CV value")
  }

  # get summary of Total_Counts #
  count_sum <- summary(filter_object$Total_Counts)

  # get biomolecules to filter if min_count is not NULL #
  if (!is.null(min_count)) {
    filt <- as.character(filter_object[filter_object$Total_Counts < min_count, 1])
    if (length(filt) == 0) filt <- NULL
  } else filt <- NULL

  res <- list(
    Summary_all = count_sum,
    filtered_biomolecules = length(filt)
  )

  class(res) = c("totalCountFiltSummary", "list")

  return(res)
}

#' Total Count Filter Print Method
#'
#' Print method for summary of Total Count filter
#' @param x the Total Count filter summary to print
#' @param ... further arguments passed to or from other methods
#'
#' @export
#' @name print.totalCountFiltSummary
#'
print.totalCountFiltSummary <- function(x, ...) {
  object <- x

  # create output #
  cat(c(
    "\nSummary of Total Count Filter\n----------------------\nCounts:\n",
    capture.output(object$Summary_all)
  ), sep = "\n")
  cat(c("\nNumber Filtered Biomolecules:", paste(object$filtered_biomolecules, collapse = ", "), "\n\n"))
}


#' Proteomics Filter Summary
#'
#' Provide summary of a proteomicsFilt S3 object
#'
#' @param object S3 object of class 'proteomicsFilt' created by
#'   \code{\link{proteomics_filter}}.
#' @param min_num_peps optional integer value between 1 and the maximum
#'   number of peptides that map to a protein in the data. The value specifies
#'   the minimum number of peptides that must map to a protein. Any protein with
#'   less than \code{min_num_peps} mapping to it will be returned as a protein
#'   that should be filtered. Default value is NULL.
#' @param degen_peps logical indicator of whether to filter out 'degenerate' or 'redundant'
#'   peptides (i.e. peptides mapping to multiple proteins) (TRUE) or not
#'   (FALSE). Default value is FALSE.
#' @return a summary table giving the number of Observed Proteins per Peptide
#'   and number of Observed Peptides per Protein. If min_num_peps is specified
#'   and/or degen_peps is TRUE, the number of biomolecules to be filtered with
#'   the specified threshold(s) are reported.
#'
#' @examples
#' library(pmartRdata)
#' myfilt <- proteomics_filter(omicsData = pep_object)
#' summary(myfilt, degen_peps = TRUE) # there are no degenerate peptides to filter out
#' summary(myfilt, min_num_peps = 2)
#'
#' @author Lisa Bramer
#'
#' @seealso \code{\link{proteomics_filter}}
#'
#' @rdname summary.proteomicsFilt
#' @name summary.proteomicsFilt
#' @export
summary.proteomicsFilt <- function(object, min_num_peps = NULL, degen_peps = FALSE,
                                   ...) {
  filter_object <- object

  # error checks for min_num_peps, if not NULL #
  if (!is.null(min_num_peps)) {
    # check that min_num_peps is not a vector #
    if (length(min_num_peps) > 1) stop("min_num_peps must be of length 1")
    # check that min_num_peps is numeric and >=1 #
    if (!is.numeric(min_num_peps) | min_num_peps < 1) stop("min_num_peps must be an integer greater than or equal to 1")
    # check that min_num_peps is an integer #
    if (min_num_peps %% 1 != 0) stop("min_num_peps must be an integer greater than or equal to 1")
    # check that min_num_peps is of length 1 #
    if (length(min_num_peps) != 1) stop("min_num_peps must be of length 1")
    # check that min_num_peps is less than the total number of peptides #
    if (min_num_peps > max(max(filter_object$counts_by_pep$n), max(filter_object$counts_by_pro$n))) stop("min_num_peps cannot be greater than the total number of peptides")
  }
  # check that degen_peps is logical #
  if (!is.logical(degen_peps)) stop("degen_peps must be either TRUE or FALSE")


  pep_sum <- summary(filter_object$counts_by_pep$n)
  pro_sum <- summary(filter_object$counts_by_pro$n)

  if (!is.null(min_num_peps)) {
    count_bypro <- filter_object$counts_by_pro
    count_bypep <- filter_object$counts_by_pep

    pro_id <- names(count_bypro)[names(count_bypro) != "n"]
    pep_id <- names(count_bypep)[names(count_bypep) != "n"]
    pro_filt <- as.character(data.frame(count_bypro[which(count_bypro$n < min_num_peps), ])[, pro_id])

    # determine which peptides no longer have a protein to map to  #
    ## find rows in peptide.info that correspond to proteins to be filtered ##
    protfilt.ids <- which(attr(filter_object, "e_meta")[, pro_id] %in% pro_filt)


    ## find the peptides that are in the filter list but are not in the unfiltered lists ##
    pep_filt <- as.character(setdiff(attr(filter_object, "e_meta")[protfilt.ids, pep_id], attr(filter_object, "e_meta")[-protfilt.ids, pep_id]))


    if (length(pep_filt) == 0) {
      pep_filt = NULL
    }
    if (length(pro_filt) == 0) {
      pro_filt = NULL
    }

    filter_object_new2 <- list(proteins_filt = pro_filt, peptides_filt = pep_filt)
  } else {
    filter_object_new2 <- list(peptides_filt = c(), proteins_filt = c())
  }

  if (degen_peps) {
    count_bypro <- filter_object$counts_by_pro
    count_bypep <- filter_object$counts_by_pep

    pro_id <- names(count_bypro)[names(count_bypro) != "n"]
    pep_id <- names(count_bypep)[names(count_bypep) != "n"]
    degen_peptides <- as.character(data.frame(count_bypep[which(count_bypep$n > 1), ])[, pep_id])

    ## identify any proteins that now will not have peptides mapping to them ##
    ## find rows in e_meta that correspond to peptides to be filtered ##
    pepfilt.ids <- which(attr(filter_object, "e_meta")[, pep_id] %in% degen_peptides)

    ## find the proteins that are in the filter list but are not in the unfiltered lists ##
    add_prots <- as.character(setdiff(attr(filter_object, "e_meta")[pepfilt.ids, pro_id], attr(filter_object, "e_meta")[-pepfilt.ids, pro_id]))

    if (length(add_prots) == 0) {
      add_prots = NULL
    }
    if (length(degen_peptides) == 0) {
      degen_peptides = NULL
    }
    filter_object_new1 <- list(peptides_filt = degen_peptides, proteins_filt = add_prots)
  } else {
    filter_object_new1 <- list(peptides_filt = c(), proteins_filt = c())
  }

  ## consolidate filter_object_new1 and filter_object_new2 ##
  filter_object_new <- list(proteins_filt = unique(c(filter_object_new1$proteins_filt, filter_object_new2$proteins_filt)), peptides_filt = unique(c(filter_object_new1$peptides_filt, filter_object_new2$peptides_filt)))

  num_filtered <- lapply(filter_object_new, length)
  # return the 5 number summary for both parts of the filter, give total number of peps and/or prots filtered
  res <- list(num_per_pep = pep_sum, num_per_pro = pro_sum, num_pep_filtered = num_filtered$peptides_filt, num_pro_filtered = num_filtered$proteins_filt)
  res$num_pro_notfiltered = nrow(filter_object$counts_by_pro) - res$num_pro_filtered
  res$num_pep_notfiltered = nrow(filter_object$counts_by_pep) - res$num_pep_filtered
  class(res) = c("proteomicsFilterSummary", "list")

  return(res)
}

#' Proteomics Filter Print Method
#'
#' Print method for summary of proteomics filter
#' @param x the proteomics filter summary to print
#' @param ... further arguments passed to or from other methods
#'
#' @export
#' @name print.proteomicsFilterSummary
#'
print.proteomicsFilterSummary <- function(x, ...) {
  object <- x

  # create output #
  catmat <- data.frame(
    c(object$num_per_pep, " ", object$num_pep_filtered, object$num_pep_notfiltered),
    c(object$num_per_pro, " ", object$num_pro_filtered, object$num_pro_notfiltered)
  )
  colnames(catmat) <- c("Obs Proteins Per Peptide", "Obs Peptides Per Protein")
  rownames(catmat) <- c(names(object$num_per_pep), " ", "Filtered", "Not Filtered")

  cat("\nSummary of Proteomics Filter \n\n")
  cat(capture.output(catmat), sep = "\n")
  cat("\n")
}

#' IMD-ANOVA Filter Summary
#'
#' Provide summary of a imdanovaFilt S3 object
#'
#' @param object S3 object of class 'imdanovaFilt' created by
#'   \code{\link{imdanova_filter}}.
#' @param min_nonmiss_gtest integer value specifying the minimum number of
#'   non-missing feature values allowed per group for \code{gtest_filter}.
#'   Defaults to NULL. Suggested value is 3.
#' @param min_nonmiss_anova integer value specifying the minimum number of
#'   non-missing feature values allowed per group for \code{anova_filter}.
#'   Defaults to NULL. Suggested value is 2.
#' @param comparisons data frame with columns for "Control" and "Test"
#'   containing the different comparisons of interest. Comparisons will be made
#'   between the Test and the corresponding Control (e.g. Control is the
#'   reference group).
#'
#' @return If min_nonmiss_gtest or min_nonmiss_anova is specified, the number of
#'   biomolecules to be filtered with the specified threshold are reported.
#'
#' @examples
#' library(pmartRdata)
#' mypep <- group_designation(omicsData = pep_object, main_effects = "Phenotype")
#' myfilt <- imdanova_filter(omicsData = mypep)
#' summary(myfilt, min_nonmiss_anova = 2, min_nonmiss_gtest = 3)
#'
#' @seealso \code{\link{imdanova_filter}}
#'
#' @author Lisa Bramer
#'
#' @rdname summary.imdanovaFilt
#' @name summary.imdanovaFilt
#' @export
summary.imdanovaFilt <- function(object, min_nonmiss_anova = NULL,
                                 min_nonmiss_gtest = NULL, comparisons = NULL, ...) {
  filter_object <- object

  ## initial checks ##

  # it is fine if both min_nonmiss_anova and min_nonmiss_gtest are NULL in the summary function #

  # check that if they aren't NULL, min_nonmiss_anova and min_nonmiss_gtest are numeric, >=2 and >=3, respectively,
  # and neither are bigger than the minimum group size (group_sizes in an attribute of the filter_object, see below) #
  if (!is.null(min_nonmiss_anova)) {
    # check that min_nonmiss_anova is not a vector #
    if (length(min_nonmiss_anova) > 1) stop("min_nonmiss_anova must be of length 1")
    # check that min_nonmiss_anova is numeric >= 2 #
    if (!is.numeric(min_nonmiss_anova) | min_nonmiss_anova < 2) stop("min_nonmiss_anova must be an integer >= 2")
    # check that min_nonmiss_anova is an integer #
    if (min_nonmiss_anova %% 1 != 0) stop("min_nonmiss_anova must be an integer >= 2")
    # check that min_nonmiss_anova is less than the minimum group size among non-singleton groups #
    nonsingleton_groups <- attributes(filter_object)$nonsingleton_groups
    if (min_nonmiss_anova > min(attributes(filter_object)$group_sizes$n_group[which(attributes(filter_object)$group_sizes$Group %in% nonsingleton_groups)])) stop("min_nonmiss_anova cannot be greater than the minimum group size")
  }
  if (!is.null(min_nonmiss_gtest)) {
    # check that min_nonmiss_gtest is not a vector #
    if (length(min_nonmiss_gtest) > 1) stop("min_nonmiss_gtest must be of length 1")
    # check that min_nonmiss_gtest is numeric >= 3 #
    if (!is.numeric(min_nonmiss_gtest) | min_nonmiss_gtest < 3) stop("min_nonmiss_gtest must be an integer >= 3")
    # check that min_nonmiss_gtest is an integer #
    if (min_nonmiss_gtest %% 1 != 0) stop("min_nonmiss_gtest must be an integer >= 3")
    # check that min_nonmiss_gtest is less than the minimum group size #
    nonsingleton_groups <- attributes(filter_object)$nonsingleton_groups
    if (min_nonmiss_gtest > min(attributes(filter_object)$group_sizes$n_group[which(attributes(filter_object)$group_sizes$Group %in% nonsingleton_groups)])) stop("min_nonmiss_gtest cannot be greater than the minimum group size")
  }

  ## end of initial checks ##

  # Check if data are paired. If they are we will filter biomolecules on paired
  # differences.
  if (!is.null(
    attr(attr(attr(filter_object, "omicsData"), "group_DF"), "pair_id")
  )) {
    # Create an omicsData object on the differences.
    diff_omicsData <- as.diffData(attr(filter_object, "omicsData"))

    # Create a new filter object with the differenced data. (For future readers:
    # I meant to use the word "differenced". I hope it made you chuckle.) We
    # need this object for the counts based on differences instead of pairs
    # because later we will filter and perform statistics on the differences.
    filter_object <- imdanova_filter(diff_omicsData)
  }

  my_names <- as.character(attr(filter_object, "group_sizes")$Group)

  # if min_nonmiss_anova is not provided, the vector of zeros will indicate nothing needs removing
  inds_rm_anova <- rep(0, nrow(filter_object))

  # if min_nonmiss_gtest is not provided, the vector of zeros will indicate nothing needs removing
  inds_rm_gtest <- rep(0, nrow(filter_object))

  # Determine what will be filtered: ANOVA -------------------------------------

  # if min_nonmiss_anova is provided
  if (!is.null(min_nonmiss_anova)) {
    # Check if there are no main effects. This can occur when data are paired.
    if ("paired_diff" %in% my_names) {
      inds_rm_anova <- filter_object$paired_diff < min_nonmiss_anova
    } else {
      # dplyr mumbo jumbo!
      # Determine which groups meet the min_nonmiss_anova criterion.
      inds_rm_anova <- filter_object %>%
        # Determine which groups have non-missing counts above the cutoff.
        dplyr::mutate(
          dplyr::across(dplyr::any_of(my_names), ~ . >= min_nonmiss_anova)
        ) %>%
        # Count the number of groups above the cutoff.
        dplyr::mutate(
          n_groups = rowSums(dplyr::across(dplyr::any_of(my_names)))
        )

      # Determine what will be filtered when comparisons is NULL.
      if (is.null(comparisons)) {
        # dplyr mumbo jumbo!
        # When comparisons are NULL we need at least two groups (across all
        # groups) that have counts above the min_nonmiss_anova threshold.
        inds_rm_anova <- inds_rm_anova %>%
          # Count the number of groups above the cutoff.
          dplyr::mutate(insufficient = n_groups < 2) %>%
          dplyr::pull(insufficient)
      } else {
        # Grab the groups in the Test and Control columns.
        testers <- unique(comparisons$Test)
        controllers <- unique(comparisons$Control)
        combiners <- unique(c(testers, controllers))

        # Sum across the unique groups in test, control, and combined (the
        # unique set of groups from both test and control). This will be used to
        # determine which rows need to be filtered depending on which scenario
        # we are in.
        inds_rm_anova <- inds_rm_anova %>%
          dplyr::mutate(
            n_test = rowSums(dplyr::across(dplyr::all_of(testers))),
            n_control = rowSums(dplyr::across(dplyr::all_of(controllers))),
            n_combine = rowSums(dplyr::across(dplyr::all_of(combiners)))
          )

        # Scenario 1: one group is compared to multiple other groups.
        if (length(controllers) == 1) {
          inds_rm_anova <- inds_rm_anova %>%
            # Count the number of groups above the cutoff.
            dplyr::mutate(insufficient = n_groups < 2) %>%
            # We can't filter the rows like the code does in applyFilt because
            # the length of inds_rm_anova and inds_rm_gest need to be the same.
            # Filtering either inds_rm_anova or inds_rm_gtest would break the
            # code at the end that counts the number of biomolecules that are
            # either filtered or not filtered.
            #
            # Change values in insufficient according to the values in n_test
            # and n_control. Rows in insufficient without any groups above the
            # cutoff in test or control must be TRUE (the corresponding
            # biomolecule is removed).
            #
            # We can't test one thing against nothing :)
            dplyr::mutate(
              insufficient = dplyr::case_when(
                n_test == 0 | n_control == 0 ~ TRUE,
                TRUE ~ insufficient
              )
            ) %>%
            dplyr::pull(insufficient)

          # Scenario 2: Some groups are compared to some other groups. In this
          # scenario a group can be in both test and control.
        } else {
          inds_rm_anova <- inds_rm_anova %>%
            # Count the number of groups above the cutoff.
            dplyr::mutate(insufficient = n_groups < 2) %>%
            # We can't filter the rows like the code does in applyFilt because
            # the length of inds_rm_anova and inds_rm_gest need to be the same.
            # Filtering either inds_rm_anova or inds_rm_gtest would break the
            # code at the end that counts the number of biomolecules that are
            # either filtered or not filtered.
            #
            # Change values in insufficient according to the values in n_test,
            # n_control and n_combine. Rows in insufficient without any groups
            # above the cutoff in control or test and rows where there is at
            # least one group in test or control but the count of combined
            # groups is less than two need to be TRUE (the corresponding
            # biomolecule is removed).
            #
            # We can't test one thing against nothing nor can we test one thing
            # against itself :)
            dplyr::mutate(
              insufficient = dplyr::case_when(
                (n_test == 0 | n_control == 0) |
                  (n_test > 0 & n_control > 0 & n_combine < 2) ~ TRUE,
                TRUE ~ insufficient
              )
            ) %>%
            dplyr::pull(insufficient)
        }
      }
    }
  }

  # Determine what will be filtered: G-Test ------------------------------------

  # if min_nonmiss_gtest is provided
  if (!is.null(min_nonmiss_gtest)) {
    # Check if there are no main effects. This can occur when data are paired.
    if ("paired_diff" %in% my_names) {
      # We cannot run a gtest when there are no main effects. Create a vector of
      # zeros so the gtest will not affect the counting below.
      inds_rm_gtest <- rep(0, dim(filter_object)[[1]])
    } else {
      # Determine what will be filtered when comparisons is NULL.
      if (is.null(comparisons)) {
        # dplyr mumbo jumbo!
        # When comparisons are NULL we need at least one group (across all
        # groups) that has a count above min_nonmiss_gtest.
        inds_rm_gtest <- filter_object %>%
          # Determine which groups have non-missing counts above the cutoff.
          dplyr::mutate(
            # Find groups above G-Test criterion.
            dplyr::across(dplyr::all_of(my_names), ~ . >= min_nonmiss_gtest),
            # Count number of groups above G-Test criterion.
            n_groups = rowSums(dplyr::across(dplyr::all_of(my_names))),
            # Determine which biomolecules have insufficient data.
            insufficient = n_groups == 0
          ) %>%
          dplyr::pull(insufficient)

        # Determine what will be filtered when custom comparisons are provided.
      } else {
        # Grab the groups in both the Test and Control columns.
        combiners <- unique(c(comparisons$Test, comparisons$Control))

        # dplyr mumbo jumbo!
        inds_rm_gtest <- filter_object %>%
          dplyr::mutate(
            # Find groups specified in comparisons above G-Test criterion.
            dplyr::across(dplyr::all_of(combiners), ~ . >= min_nonmiss_gtest),
            # Count number of groups above G-Test criterion.
            n_groups = rowSums(dplyr::across(dplyr::all_of(combiners))),
            # Determine which biomolecules have insufficient data.
            insufficient = n_groups == 0
          ) %>%
          dplyr::pull(insufficient)
      }
    }
  }

  # Count number of things that will be filtered -------------------------------

  # if-statement for the case where both min_nonmiss_anova and min_nonmiss_gtest
  # are non-NULL. Note: num_not_tested will be the count of (inds_rm_anova +
  # inds_rm_gtest) that sum to 2. That is the intersection of both filters.
  if (!is.null(min_nonmiss_anova) & !is.null(min_nonmiss_gtest)) {
    # We will count differently if the data are paired and there are no main
    # effects. In this case we only need the sum of inds_rm_anova and
    # inds_rm_gtest to be greater than zero.
    if ("paired_diff" %in% my_names) {
      # get number molecules not tested
      num_filtered <- sum((inds_rm_anova + inds_rm_gtest) > 0)
    } else {
      # get number molecules not tested
      num_filtered <- sum((inds_rm_anova + inds_rm_gtest) > 1)
    }

    # get number molecules tested
    num_not_filtered <- nrow(filter_object) - num_filtered
  } else if (!is.null(min_nonmiss_anova) | !is.null(min_nonmiss_gtest)) {
    # get number molecules not tested
    num_filtered <- sum((inds_rm_anova + inds_rm_gtest) > 0)

    # get number molecules tested
    num_not_filtered <- nrow(filter_object) - num_filtered
  } else {
    num_not_filtered = NULL
    num_filtered = NULL
  }

  res <- list(
    pep_observation_counts = nrow(filter_object),
    num_filtered = num_filtered,
    num_not_filtered = num_not_filtered
  )

  class(res) = c("imdanovaFilterSummary", "list")

  return(res)
}


#' IMD-ANOVA Filter Print Method
#'
#' Print method for summary of imdanova filter
#' @param x the imdanova filter summary to print
#' @param ... further arguments passed to or from other methods
#'
#' @export
#' @name print.imdanovaFilterSummary
#'
print.imdanovaFilterSummary <- function(x, ...) {
  object <- x

  # create output #
  if (is.null(object$num_filtered)) {
    object$num_filtered <- NA
  }

  if (is.null(object$num_not_filtered)) {
    object$num_not_filtered <- NA
  }

  catmat <- data.frame(c(object$pep_observation_counts, object$num_filtered, object$num_not_filtered))
  rownames(catmat) <- c("Total Observations: ", "Filtered: ", "Not Filtered: ")

  colnames(catmat) <- NULL

  cat("\nSummary of IMD-ANOVA Filter\n")
  cat(capture.output(catmat), sep = "\n")
}


#' RMD Filter Summary
#'
#' Provide summary of a rmdFilt S3 object
#'
#' @param object S3 object of class 'rmdFilt' created by
#'   \code{\link{rmd_filter}}.
#' @param pvalue_threshold A threshold for the Robust Mahalanobis Distance (RMD)
#'   p-value. All samples with p-values below the threshold will be filtered
#'   out. Default value is NULL. Suggested value is 0.0001
#' @return a summary of the p-values associated with running RMD-PAV across all
#'   samples. If a p-value threshold is provided the samples that would be
#'   filtered at this threshold are reported.
#'
#' @examples
#' library(pmartRdata)
#' mymetab <- group_designation(omicsData = metab_object, main_effects = "Phenotype")
#' mymetab <- edata_transform(omicsData = mymetab, data_scale = "log2")
#' myfilt <- rmd_filter(omicsData = mymetab)
#' summary(myfilt, pvalue_threshold = 0.001)
#'
#' @author Lisa Bramer, Kelly Stratton
#'
#' @seealso \code{\link{rmd_filter}}
#'
#' @export
#' @rdname summary.rmdFilt
#' @name summary.rmdFilt
summary.rmdFilt <- function(object, pvalue_threshold = NULL, ...) {
  filter_object <- object

  # check that pvalue_threshold is numeric [0,1] #
  if (!is.null(pvalue_threshold)) {
    if (!is.numeric(pvalue_threshold)) stop("pvalue_threshold must be numeric between 0 and 1")
    if (pvalue_threshold < 0 | pvalue_threshold > 1) stop("pvalue_threshold must be numeric between 0 and 1")
  }

  # get metrics used #
  samp_id <- names(attr(filter_object, "group_DF"))[1]
  metrics <- attributes(filter_object)$metrics

  # get samples filtered out, if pvalue_threshold is specified #
  filt <- NULL
  if (!is.null(pvalue_threshold)) {
    filt <- as.character(filter_object[filter_object$pvalue < pvalue_threshold, samp_id])
    if (length(filt) == 0) filt <- NULL
  }

  # Check if data are paired. If they are we will make sure no pairs are split
  # (e.g., one sample in a pair will be filtered but the other sample will not).
  if (!is.null(attr(attr(filter_object, "group_DF"), "pair_id"))) {
    # Snatch the sample name and pair name as they will be used in multiple
    # places. It will save some typing. ... However, all the typing I just saved
    # has probably been undone by writing this comment.
    sample_name <- attr(filter_object, "fdata_cname")
    pair_name <- attr(attr(filter_object, "group_DF"), "pair_id")

    # !#!#!#!#!#!#!#!#!#!
    # The following code checks if the samples in a pair will be split. For
    # example, one sample in a pair will be filtered and another sample in the
    # pair will not be filtered. If a pair is split remove ALL samples
    # associated with that pair.
    # !#!#!#!#!#!#!#!#!#!

    # Snag the associated pair IDs for the samples that will be filtered.
    filtered_pairs <- attr(filter_object, "fdata") %>%
      dplyr::filter(!!rlang::sym(sample_name) %in% filt) %>%
      dplyr::pull(!!rlang::sym(pair_name))

    # Go back to f_data and nab all the sample names corresponding to the pair
    # IDs associated with the original samples that were selected for removal.
    # These sample names will be used to filter the data. This vector will not
    # be the same as the original vector if any pairs were split. If there
    # were no pairs split the two vectors will be the same.
    filt <- attr(filter_object, "fdata") %>%
      dplyr::filter(!!rlang::sym(pair_name) %in% filtered_pairs) %>%
      dplyr::pull(!!rlang::sym(sample_name)) %>%
      as.character()
  }

  # for display #
  if (is.null(filt)) {
    filt <- "NULL"
  }

  res <- list(pvalue = summary(filter_object$pvalue), metrics = metrics, filtered_samples = filt)

  class(res) = c("rmdFilterSummary", "list")

  return(res)
}

#' RMD Filter Print Method
#'
#' Print method for summary of RMD filter
#' @param x the RMD filter summary to print
#' @param ... further arguments passed to or from other methods
#'
#' @export
#' @name print.rmdFilterSummary
#'
print.rmdFilterSummary <- function(x, ...) {
  object <- x

  # create output #
  cat(c("\nSummary of RMD Filter\n----------------------\nP-values:\n", capture.output(object$pvalue)), sep = "\n")
  cat(c("\nMetrics Used:", paste(object$metrics, collapse = ", "), "\n"))
  cat(c("\nFiltered Samples:", paste(object$filtered_samples, collapse = ", "), "\n\n"))
}


#' Coefficient of Variation (CV) Filter Summary
#'
#' Provide summary of a cvFilt S3 object
#'
#' @param filter_object S3 object of class 'cvFilt' created by
#'   \code{\link{cv_filter}}.
#' @param cv_threshold numeric value greater than 1 and less than the value
#'   given by filter_object$CV. CV values above cv_threshold are filtered out.
#'   Default value is NULL.
#' @return a summary of the CV values, number of NA values, and non-NA values.
#'   If a CV threshold is provided, the biomolecules that would be filtered
#'   based on this threshold are reported.
#'
#' @examples
#' library(pmartRdata)
#' mypep <- group_designation(omicsData = pep_object, main_effects = "Phenotype")
#' to_filter <- cv_filter(omicsData = mypep, use_groups = TRUE)
#' summary(to_filter, cv_threshold = 30)
#'
#' @author Lisa Bramer
#'
#' @seealso \code{\link{cv_filter}}
#'
#' @export
#' @rdname summary.cvFilt
#' @name summary.cvFilt
summary.cvFilt <- function(object, cv_threshold = NULL, ...) {
  filter_object <- object

  # checks for cv_threshold if not null #
  if (!is.null(cv_threshold)) {
    # check that cv_threshold is numeric
    if (!is.numeric(cv_threshold)) stop("cv_threshold must be numeric of length 1")
    # chack that cv_threshold is of length 1
    if (length(cv_threshold) > 1) stop("cv_threshold must be numeric of length 1")
    # check that cv_threshold is more than 1 and less than max CV value
    if (cv_threshold <= 1 | cv_threshold >= max(filter_object$CV, na.rm = TRUE)) stop("cv_threshold must be greater than 1 and less than the maximum CV value")
  }

  # get rid of NAs #
  new_object <- filter_object[!is.na(filter_object$CV), ]

  # get summary of CVs #
  CVs <- summary(new_object$CV)

  if (!is.null(attributes(filter_object)$tot_nas)) {
    # get total NAs, total non-NAs #
    tot_NAs <- attributes(filter_object)$tot_nas
    tot_non_NAs <- nrow(filter_object) - tot_NAs
  } else {
    # get total zeros, total non-zeros #
    tot_zeros <- attributes(filter_object)$tot_zeros
    tot_non_zeros <- nrow(filter_object) - tot_zeros
  }

  # get biomolecules to filter if cv_threshold is not NULL #
  filt <- NULL
  if (!is.null(cv_threshold)) {
    filt <- as.character(new_object[new_object$CV > cv_threshold, 1])
    if (length(filt) == 0) filt <- NULL
  }

  if (!is.null(attributes(filter_object)$tot_nas)) {
    res <- list(
      CVs = CVs,
      tot_NAs = tot_NAs,
      tot_non_NAs = tot_non_NAs,
      filtered_biomolecules = length(filt)
    )
  } else {
    res <- list(
      CVs = CVs,
      tot_zeros = tot_zeros,
      tot_non_zeros = tot_non_zeros,
      filtered_biomolecules = length(filt)
    )
  }

  class(res) = c("cvFilterSummary", "list")

  return(res)
}

#' CV Filter Print Method
#'
#' Print method for summary of CV filter
#' @param x the CV filter summary to print
#' @param ... further arguments passed to or from other methods
#'
#' @export
#' @name print.cvSummary
#'
print.cvFilterSummary <- function(x, ...) {
  object <- x

  # create output #
  cat(c("\nSummary of Coefficient of Variation (CV) Filter\n----------------------\nCVs:\n", capture.output(object$CVs)), sep = "\n")
  if (!is.null(object$tot_NAs)) {
    cat(c("\nTotal NAs:", object$tot_NAs))
    cat(c("\nTotal Non-NAs:", object$tot_non_NAs, "\n"))
  } else {
    cat(c("\nTotal zeros:", object$tot_zeros))
    cat(c("\nTotal Non-zeros:", object$tot_non_zeros, "\n"))
  }
  cat(c("\nNumber Filtered Biomolecules:", paste(object$filtered_biomolecules, collapse = ", "), "\n\n"))
}


#' Custom Filter Summary
#'
#' Provide summary of a customFilt S3 object
#'
#' @param object S3 object of class 'customFilt' created by
#'   \code{\link{custom_filter}}.
#'
#' @return a summary of the items in e_data, f_data, and e_meta that will be
#'   removed as a result of applying the custom filter.
#'
#' @examples
#' library(pmartRdata)
#' to_filter <- custom_filter(omicsData = metab_object, e_data_remove = "fumaric acid", f_data_remove = "Sample_1_Phenotype2_B")
#' summary(to_filter)
#'
#' to_filter2 <- custom_filter(omicsData = metab_object, f_data_keep = metab_object$f_data$SampleID[1:10])
#' summary(to_filter2)
#'
#' @author Lisa Bramer
#'
#' @seealso \code{\link{custom_filter}}
#'
#' @export
#' @rdname summary.customFilt
#' @name summary.customFilt
summary.customFilt <- function(object, ...) {
  filter_object <- object

  # get omicsData object #
  omicsData <- attr(filter_object, "omicsData")
  summary_orig <- summary(omicsData)

  # get names #
  edata_id <- get_edata_cname(omicsData)
  emeta_id <- get_emeta_cname(omicsData)
  samp_id <- get_fdata_cname(omicsData)

  # get counts #
  num_samples <- attr(filter_object, "num_samples")
  num_edata <- attr(filter_object, "num_edata")
  num_emeta <- attr(filter_object, "num_emeta") ## is null where not used

  # apply the filter #
  filtered_data <- applyFilt(filter_object, omicsData)
  summary_filt <- summary(filtered_data)

  mode_rmv <- !is.null(filter_object$e_data_remove) |
    !is.null(filter_object$f_data_remove) |
    !is.null(filter_object$e_meta_remove)

  mode_kp <- !is.null(filter_object$e_data_keep) |
    !is.null(filter_object$f_data_keep) |
    !is.null(filter_object$e_meta_keep)

  # if filter_object contains removes
  if (mode_rmv) {
    # samples #
    if (!is.null(filter_object$f_data_remove)) {
      samps_filt <- length(filter_object$f_data_remove)
      samps_left <- num_samples - samps_filt
    } else {
      samps_filt <- 0
      samps_left <- num_samples
    }

    # e_data #
    if (!is.null(filter_object$e_data_remove)) {
      edata_filt <- length(filter_object$e_data_remove)
      edata_left <- num_edata - edata_filt
    } else {
      edata_filt <- num_edata - nrow(filtered_data$e_data)
      edata_left <- nrow(filtered_data$e_data)
    }

    # e_meta #
    if (!is.null(num_emeta)) {
      if (!is.null(filter_object$e_meta_remove)) {
        emeta_filt <- length(filter_object$e_meta_remove)
        emeta_left <- num_emeta - emeta_filt
      } else {
        emeta_filt <- num_emeta - length(unique(filtered_data$e_meta[, emeta_id]))
        emeta_left <- length(unique(filtered_data$e_meta[, emeta_id]))
      }
    }

    disp_colnames <- c("Filtered", "Remaining", "Total")
  }

  # if filter_object contains keeps
  # tags of `_left` and `_filt` correspond to items kept and items not kept
  if (mode_kp) {
    # samples #
    if (!is.null(filter_object$f_data_keep)) {
      samps_left <- length(filter_object$f_data_keep)
      samps_filt <- num_samples - samps_left
    } else {
      samps_left <- num_samples
      samps_filt <- 0
    }

    # e_data #
    if (!is.null(filter_object$e_data_keep)) {
      edata_left <- length(filter_object$e_data_keep)
      edata_filt <- num_edata - edata_left
    } else {
      edata_left <- nrow(filtered_data$e_data)
      edata_filt <- num_edata - nrow(filtered_data$e_data)
    }

    # e_meta #
    if (!is.null(num_emeta)) {
      if (!is.null(filter_object$e_meta_keep)) {
        emeta_left <- length(filter_object$e_meta_keep)
        emeta_filt <- num_emeta - emeta_left
      } else {
        emeta_left <- length(unique(filtered_data$e_meta[, emeta_id]))
        emeta_filt <- num_emeta - length(unique(filtered_data$e_meta[, emeta_id]))
      }
    }

    disp_colnames <- c("Discarded", "Kept", "Total")
  }

  # Display #

  ## Display text ##
  edata_plural_text <- ifelse(length(grep("s$", edata_id)) > 0, " ", "s ")
  emeta_plural_text <- ifelse(length(grep("s$", emeta_id)) > 0, " ", "s ")
  samp_plural_text <- ifelse(length(grep("s$", samp_id)) > 0, " ", "s ")
  samp_id <- paste(samp_id, samp_plural_text, "(f_data)", sep = "")
  edata_id <- paste(edata_id, edata_plural_text, "(e_data)", sep = "")
  display_emeta_id <- paste(emeta_id, emeta_plural_text, "(e_meta)", sep = "")

  ## construct data frame ##
  if (is.null(emeta_id)) {
    disp <- data.frame(
      c(samps_filt, edata_filt),
      c(samps_left, edata_left),
      c(num_samples, num_edata)
    )
    rownames(disp) <- c(samp_id, edata_id)
  } else {
    disp <- data.frame(
      c(samps_filt, edata_filt, emeta_filt),
      c(samps_left, edata_left, emeta_left),
      c(num_samples, num_edata, num_emeta)
    )
    rownames(disp) <- c(samp_id, edata_id, display_emeta_id)
  }

  colnames(disp) <- disp_colnames

  class(disp) = c("customFilterSummary", "data.frame")

  return(disp)
}


#' Custom Filter Print Method
#'
#' Print method for summary of custom filter
#' @param x the custom filter summary to print
#' @param ... further arguments passed to or from other methods
#'
#' @export
#' @name print.customFilterSummary
#'
print.customFilterSummary <- function(x, ...) {
  object <- x

  ## Display output ##
  cat("\nSummary of Custom Filter\n\n")
  print.data.frame(object)
  cat("\n")
}
