% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter_summary.R
\name{summary.imdanovaFilt}
\alias{summary.imdanovaFilt}
\title{Produce a basic summary of a imdanova_filter object}
\usage{
\method{summary}{imdanovaFilt}(filter_object, min_nonmiss_anova = NULL, min_nonmiss_gtest = NULL)
}
\arguments{
\item{filter_object}{S3 object of class 'imdanovaFilt' created by \code{\link{imdanova_filter}}.}

\item{min_nonmiss_anova}{the minimum number of non-missing feature values allowed per group for \code{anova_filter}. Defaults to NULL. Suggested value is 2.}

\item{min_nonmiss_gtest}{the minimum number of non-missing feature values allowed per group for \code{gtest_filter}. Defaults to NULL. Suggested value is 3.}
}
\value{
If min_nonmiss_gtest or min_nonmiss_anova is specified, the number of biomolecules to be filtered with the specified threshold are reported.
}
\description{
This function will provide basic summary statistics for the imdanova_filter object.
}
\author{
Lisa Bramer
}
