#' Create a Melted and Grouped Version of e_data for IMD_ANOVA filter
#'
#' This function creates a melted version of e_data, grouped by edata_id and group designation, for future use of implementing a IMD_ANOVA filter
#'
#' @param e_data \eqn{p \times n} data.frame, where \eqn{p} is the number of peptides, proteins, lipids, metabolites, or accessions and \eqn{n} is the number of samples. e_data contains
#' @param groupDF data.frame created by \code{group_designation} with columns for the sample identifier and the designated group.
#' @param samp_id character string specifying the name of the column containing the sample identifiers in \code{groupDF}.
#'
#' @return a data.frame of class "grouped_dt" which is compatible with functions in the dplyr package
#'
#'
#' @author Lisa Bramer, Kelly Stratton

pre_imdanova_melt <- function(e_data, groupDF, samp_id){


  # convert peptide data to class "data.table" #
  DT <- data.table::data.table(e_data)
  data.table::setnames(DT,names(DT)[1], "Peptide") # Don't change this from "Peptide"! When I tried to use cname_id, it messed up use of the "summarise" function in nonmissing_per_group.R

  # melt this data based on the first column giving the peptide names #
  melt.pep = data.table:::melt.data.table(DT, id.var = 1)

  # set second column to have name sampleID #
  #data.table::setnames(melt.pep,names(melt.pep)[2], samp_id)
  data.table::setnames(melt.pep,names(melt.pep)[2], "samp_id")

  orig_name <- names(groupDF)[1] # save the orig name
  names(groupDF)[1] <- "samp_id"

  # merge melted e_data and groupDF #
  new.dat = data.table:::merge.data.table(x = melt.pep, y = data.table::data.table(groupDF),
                                          by = "samp_id", all.x = TRUE)

  # group data by Peptide and Group #
  grpd.data = dplyr::group_by(new.dat, Peptide, Group) # Don't change this from Peptide! When I tried to use cname_id, it messed up use of the "summarise" function in nonmissing_per_group.R  ... after the "summarise" is done, it is reset to cname_id

  names(grpd.data)[1] <- orig_name # re-set the sample ID column name back to the original value
  # return grouped data #
  return(grpd.data)

}
