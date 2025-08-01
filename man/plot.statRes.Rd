% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_fns.R
\name{plot.statRes}
\alias{plot.statRes}
\title{Plot statRes Object}
\usage{
\method{plot}{statRes}(
  x,
  plot_type = "bar",
  fc_threshold = NULL,
  fc_colors = c("red", "black", "green"),
  stacked = TRUE,
  show_sig = TRUE,
  color_low = NULL,
  color_high = NULL,
  plotly_layout = NULL,
  interactive = FALSE,
  x_lab = NULL,
  x_lab_size = 11,
  x_lab_angle = NULL,
  y_lab = NULL,
  y_lab_size = 11,
  title_lab = NULL,
  title_lab_size = 14,
  legend_lab = NULL,
  legend_position = "right",
  text_size = 3,
  bw_theme = TRUE,
  display_count = TRUE,
  custom_theme = NULL,
  cluster = FALSE,
  free_y_axis = FALSE,
  ...
)
}
\arguments{
\item{x}{`statRes` object to be plotted, usually the result of `imd_anova`}

\item{plot_type}{defines which plots to be produced, options are "bar",
"volcano", "gheatmap", "fcheatmap"; defaults to "bar".  See details for
plot descriptions.}

\item{fc_threshold}{optional threshold value for fold change estimates.
Modifies the volcano plot as follows:  Vertical lines are added at
(+/-)\code{fc_threshold} and all observations that have absolute fold
change less than \code{abs(fc_threshold)} are colored as 'non-significant'
(as specified by \code{fc_colors}).}

\item{fc_colors}{vector of length three with character color values
interpretable by ggplot. i.e. c("orange", "black", "blue") with the values
being used to color negative, non-significant, and positive fold changes
respectively. When \code{plot_type = "histogram"}, this a vector of length two
with character color values interpretable by ggplot. i.e. c("green","gray") with the
values being used to color significant (regardless of direction)
and non-significant fold changes}

\item{stacked}{TRUE/FALSE for whether to stack positive and negative fold
change sections in the barplot, defaults to TRUE}

\item{show_sig}{This input is used when \code{plot_type = "gheatmap"}. A
logical value. If TRUE a visual indicator that a certain bin combination is
significant by the g-test is shown. This input is also used when 
when \code{plot_type = "histogram"}. If TRUE, only significant fold changes 
will be plotted.}

\item{color_low}{This input is used when \code{plot_type = "gheatmap"}. A
character string specifying the color of the gradient for low count values.}

\item{color_high}{This input is used when \code{plot_type = "gheatmap"}. A
character string specifying the color of the gradient for high count
values.}

\item{plotly_layout}{This input is used when \code{plot_type = "gheatmap"}. A
list of arguments, not including the plot, to be passed to
\code{plotly::layout} if \code{interactive = TRUE}.}

\item{interactive}{TRUE/FALSE for whether to create an interactive plot using
plotly. Not valid for all plots.}

\item{x_lab}{character string specifying the x-axis label.}

\item{x_lab_size}{integer value indicating the font size for the x-axis. The
default is 11.}

\item{x_lab_angle}{integer value indicating the angle of x-axis labels.}

\item{y_lab}{character string specifying the y-axis label.}

\item{y_lab_size}{integer value indicating the font size for the y-axis. The
default is 11.}

\item{title_lab}{character string specifying the plot title.}

\item{title_lab_size}{integer value indicating the font size of the plot
title. The default is 14.}

\item{legend_lab}{character string specifying the legend title.}

\item{legend_position}{character string specifying the position of the
legend. Can be one of "right", "left", "top", "bottom", or "none". The
default is "none".}

\item{text_size}{integer specifying the size of the text (number of
non-missing values) within the plot. The default is 3.}

\item{bw_theme}{logical value. If TRUE uses the ggplot2 black and white
theme.}

\item{display_count}{logical value. Indicates whether the non-missing counts
will be displayed on the bar plot. The default is TRUE.}

\item{custom_theme}{a ggplot `theme` object to be applied to non-interactive
plots, or those converted by plotly::ggplotly().}

\item{cluster}{logical for heatmaps; TRUE will cluster biomolecules on X
axis. defaults to TRUE for seqData statistics and FALSE for all others.}

\item{free_y_axis}{Logical. If TRUE the y axis for each bar plot can have its
own range. The default is FALSE.}

\item{...}{further arguments passed to or from other methods.}
}
\value{
ggplot2 plot object if interactive is FALSE, or plotly plot object if
  interactive is TRUE
}
\description{
Produces plots that summarize the results contained in a `statRes` object.
}
\details{
Plot types:
\itemize{
 \item{"bar"} \code{?pmartR::statres_barplot} Bar-chart with bar heights
 indicating the number of significant biomolecules, grouped by test type and
 fold change direction.
 \item{"volcano"} \code{?pmartR::statres_volcano_plot} Scatter plot showing
 negative-log-pvalues against fold change. Colored by statistical
 significance and fold change.
 \item{"gheatmap"} \code{?pmartR::gtest_heatmap} Heatmap with x and y axes
 indicating the number of nonmissing values for two groups. Colored by
 number of biomolecules that fall into that combination of nonmissing values.
 \item{"fcheatmap"} Heatmap showing all biomolecules across comparisons,
 colored by fold change.
 \item{"histogram"} Histogram showing the distribution of log2 fold changes
 across comparisons, colored by significance level status
}
}
\examples{
\dontshow{if (identical(tolower(Sys.getenv("NOT_CRAN")), "true") & requireNamespace("pmartRdata", quietly = TRUE)) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
library(pmartRdata)
# Group the data by condition
mypro <- group_designation(
  omicsData = pro_object,
  main_effects = c("Phenotype")
)

# Apply the IMD ANOVA filter
imdanova_Filt <- imdanova_filter(omicsData = mypro)
mypro <- applyFilt(
  filter_object = imdanova_Filt,
  omicsData = mypro,
  min_nonmiss_anova = 2
)

# Implement the IMD ANOVA method and compuate all pairwise comparisons
# (i.e. leave the `comparisons` argument NULL)
anova_res <- imd_anova(omicsData = mypro, test_method = 'anova')
plot(anova_res)
plot(anova_res, plot_type = "volcano")

imd_res <- imd_anova(omicsData = mypro, test_method = 'gtest')
plot(imd_res)

imd_anova_res <- imd_anova(
  omicsData = mypro,
  test_method = 'comb',
  pval_adjust_a_multcomp = 'bon',
  pval_adjust_g_multcomp = 'bon'
)
plot(imd_anova_res, bw_theme = TRUE)
plot(imd_anova_res, plot_type = "volcano", bw_theme = TRUE)
\dontshow{\}) # examplesIf}
}
