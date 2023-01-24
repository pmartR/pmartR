#' Conversion between log2(RMD) and p-value
#'
#' This function provides a conversion between the log base 2 robust Mahalanobis
#' distance value and p-value for output from the \code{rmd_runs} function
#'
#' @param log2rmd numeric log base 2 transformed robust Mahalanobis distance
#'   value
#' @param pval numeric p-value associated with rmd_runs algorithm
#' @param df integer value specifying the degrees of freedom associated with the
#'   test, which should be equal to the number of metrics used in rmd_runs
#'
#' @details Only one of \code{log2rmd} and \code{pval} should be provided. The
#'   input not provided will be solved for based on the provided input.
#'
#' @return The function returns the corresponding p-value or log base 2 robust
#'   Mahalanobis when the other parameter is specified.
#'
#' @examples
#' library(pmartRdata)
#' mymetab <- edata_transform(omicsData = metab_object,
#'                            data_scale = "log2")
#' mymetab <- group_designation(omicsData = mymetab,
#'                              main_effects = "Phenotype")
#' rmd_results <- rmd_filter(omicsData = mymetab,
#'                           metrics=c("MAD", "Skewness", "Correlation"))
#' rmd_conversion(log2rmd = rmd_results$Log2.md, df=3)
#'
#' rmd_conversion(pval = .0001, df = 3)
#' rmd_conversion(log2rmd = 4.5, df = 3)
#'
#' @author Lisa Bramer
#'
#' @export
#'
rmd_conversion <- function (log2rmd = NULL, pval = NULL, df) {

  # check that only one of the two log2rmd and pval is provided #
  if(!is.null(log2rmd) & !is.null(pval))
    stop("Too many arguments entered to solve for log2rmd or pval.")

  # check that one of the two log2rmd and pval is provided #
  if(is.null(log2rmd) & is.null(pval))
    stop("One of log2rmd and pval must be provided")

  # Make sure the user doesn't do something totally stupid like enter the
  # degrees of freedom as a character.
  if (!is.numeric(df)) stop ("df must be numeric.")

  # check that df is at least 1 #
  if(df < 1) stop("df must be an integer greater than or equal to 1.")

  # if log2rmd is provided solve for p-value #
  if(!is.null(log2rmd)){
    rmd = 2^(log2rmd)
    res = 1 - pchisq(rmd, df = df)
  }

  # if pval is provided solve for log2rmd #
  if(!is.null(pval)){

    # The user is probably an idiot and helping them this much only enables them
    # to continue to be an idiot. However, I am merciful and will help them out
    # a bit.
    if (!is.numeric(pval)) stop ("pval must be numeric.")
    if (pval < 0 || pval > 1) stop ("pval must be between 0 and 1.")

    invp = 1 - pval
    res = log(qchisq(invp, df = df), base = 2)
  }

  return(res)
}
