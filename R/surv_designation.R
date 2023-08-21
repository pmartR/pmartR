#' Create a "surv_DF" attribute so that survival analysis can be implemented.
#'
#' This function will add the necessary information to omicsData such that
#' survival analysis can be applied to it.
#'
#' @param omicsData an object of the class 'lipidData', 'metabData', 'pepData',
#'   or 'proData' usually created by \code{\link{as.lipidData}},
#'   \code{\link{as.metabData}}, \code{\link{as.pepData}}, or
#'   \code{\link{as.proData}}, respectively.
#' @param t_death the column in 'f_data' that corresponds to the subjects' time
#'   of death
#' @param t_progress the column in 'f_data' that corresponds to the subjects'
#'   time of progression
#' @param ind_death the column in 'f_data' that corresponds to the subjects'
#'   status, e.g. alive/dead
#' @param ind_progress the column in 'f_data' that corresponds to the subjects'
#'   progression status
#' @param covariates the column(s) in 'f_data' that correspond to covariates to
#'   be included in the survivial analysis
#'
#' @return omicsData is returned with the additional attribute
#'
#' @author Bryan Stanfill
#'
surv_designation <- function(omicsData, t_death, t_progress = NULL, ind_death, ind_progress = NULL, covariates = NULL) {
  attr(omicsData, "surv_DF") <- list(t_death = t_death, t_progress = t_progress, ind_death = ind_death, ind_progress = ind_progress, covariates = covariates)

  return(omicsData)
}
