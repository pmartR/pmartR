#' Test the location and scale parameters from a normalization procedure
#'
#' Computes p-values from a test of dependence between normalization parameters
#' and group assignment of a normalized omicsData or normRes object
#'
#' @param norm_obj object of class 'pepData', 'proData', 'lipidData',
#'   'metabData', 'isobaricpepData', or 'nmrData' that has had \code{normalize_global()} run
#'   on it, or a 'normRes' object
#' @param test_fn character string indicating the statistical test to use.
#'   Current options are "anova" and "kw" for a Kruskal-Wallis test.
#'
#' @return A list with 2 entries containing the p_value of the test performed on
#'   the location and scale (if it exists) parameters.
#'
#' @examples
#' library(pmartRdata)
#' mymetab <- edata_transform(omicsData = metab_object, data_scale = "log2")
#' mymetab <- group_designation(omicsData = mymetab, main_effects = "Phenotype")
#'
#' # provide the normRes object
#' mynorm <- normalize_global(omicsData = mymetab, subset_fn = "all",
#'                            norm_fn = "median", apply_norm = FALSE)
#' norm_pvals <- normRes_tests(norm_obj = mynorm)
#'
#' # provide normalized omicsData object
#' mymetab <- normalize_global(omicsData = mymetab, subset_fn = "all",
#'                             norm_fn = "median", apply_norm = TRUE)
#' norm_pvals <- normRes_tests(norm_obj = mymetab)
#'
#' # NMR data object
#' mynmr <- edata_transform(omicsData = nmr_identified_object, data_scale = "log2")
#' mynmr <- group_designation(mynmr, main_effects = "Condition")
#' mynmrnorm <- normalize_nmr(
#'   omicsData = mynmr,
#'   apply_norm = TRUE,
#'   sample_property_cname = "Concentration"
#' )
#' mynmrnorm <- normalize_global(omicsData = mynmrnorm, subset_fn = "all",
#'                               norm_fn = "median", apply_norm = TRUE, backtransform = TRUE)
#' norm_pvals <- normRes_tests(norm_obj = mynmrnorm)
#'
#' @export
#'
normRes_tests <- function(norm_obj, test_fn = "kw") {
  obj_classes <- c(
    "pepData", "proData", "lipidData", "metabData",
    "isobaricpepData", 'nmrData'
  )

  # object is of correct class and groupie is a character vector
  if (!inherits(norm_obj, c(obj_classes, "normRes")))
    stop("This object is not of class 'pepData', 'proData', 'lipidData', 'metabData', 'isobaricpepData' or 'normRes'")

  # if this is an omicsData object, test that it has been normalized and has a
  # group DF ... if so store the location and scale vectors that should be
  # under $data_info$norm_info and get the group vector from group DF if none
  # was specified

  # Make sure test_fn is one of the possible options.
  if (!test_fn %in% c("anova", "kw")) {
    # Stop! You are attempting to call a function that does not exist. This is
    # similar to trying to divide by zero. I am stopping you before you cause
    # severe mental and physical harm to yourself and others.
    stop("test_fn must be either 'anova' or 'kw'.")
  }

  # Extract information when an omicsData object is the input. ---------
  if (inherits(norm_obj, obj_classes)) {
    if (!attributes(norm_obj)$data_info$norm_info$is_normalized)
      stop("Normalization has not been run on this data")
    if (is.null(get_group_DF(norm_obj)))
    # if (is.null(attributes(norm_obj)$group_DF))
      stop("No grouping structure present in object")

    location <- attributes(norm_obj)$data_info$norm_info$params$norm_location
    scale <- attributes(norm_obj)$data_info$norm_info$params$norm_scale
    groupDF <- attributes(norm_obj)$group_DF

    # Extract sample names from omicsData$e_data. This will be used to order the
    # group column from the groupDF attribute.
    enames <- names(norm_obj$e_data)

    # Find the index of the biomolecule ID column in e_data. It will be removed
    # from the names vector later on.
    eidx <- which(enames == get_edata_cname(norm_obj))

    # Nab fdata_cname from the omicsData object.
    fname <- get_fdata_cname(norm_obj)

    # Extract information when a normRes object is the input. ---------
  } else if (inherits(norm_obj, "normRes")) {
    if (is.null(attr(attr(norm_obj, "omicsData"), "group_DF")))
      stop("No grouping structure present in object")

    location <- norm_obj$parameters$normalization$location
    scale <- norm_obj$parameters$normalization$scale
    groupDF <- attr(attr(norm_obj, "omicsData"), "group_DF")

    # Extract sample names from omicsData$e_data. This will be used to order the
    # group column from the groupDF attribute.
    enames <- names(attr(norm_obj, "omicsData")$e_data)

    # Find the index of the biomolecule ID column in e_data. It will be removed
    # from the names vector later on.
    eidx <- which(enames == get_edata_cname(attr(norm_obj, "omicsData")))

    # Nab fdata_cname from the omicsData object.
    fname <- get_fdata_cname(attr(norm_obj, "omicsData"))
  }

  # Order the groups so they are in the same order as the samples of e_data.
  groupie <- groupDF$Group[match(enames[-eidx], groupDF[[fname]])]

  # get p values using the value and group vectors and selected test
  if (test_fn == "kw") {
    p_location <- kw_rcpp(matrix(location, nrow = 1), groupie)
    p_scale <- if (!is.null(scale))
      kw_rcpp(matrix(scale, nrow = 1), groupie) else
      NULL
  } else if (test_fn == "anova") {
    p_location <- aov(location ~ groupie) %>%
      summary() %>%
      {
        .[[1]]$`Pr(>F)`[1]
      }
    p_scale <- if (!is.null(scale))
      aov(scale ~ groupie) %>%
        summary() %>%
        {
          .[[1]]$`Pr(>F)`[1]
        } else
      NULL
  }

  res <- list(p_location = p_location, p_scale = p_scale)

  return(res)
}
