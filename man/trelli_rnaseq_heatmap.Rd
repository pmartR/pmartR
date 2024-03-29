% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trelliPlots_seqData.R
\name{trelli_rnaseq_heatmap}
\alias{trelli_rnaseq_heatmap}
\title{Heatmap trelliscope building function for RNA-seq data}
\usage{
trelli_rnaseq_heatmap(
  trelliData,
  cognostics = c("sample count", "mean LCPM", "biomolecule count"),
  ggplot_params = NULL,
  interactive = FALSE,
  path = .getDownloadsFolder(),
  name = "Trelliscope",
  test_mode = FALSE,
  test_example = 1,
  single_plot = FALSE,
  ...
)
}
\arguments{
\item{trelliData}{A trelliscope data object made by as.trelliData, and
grouped by an emeta variable. Must be built using seqData. Required.}

\item{cognostics}{A vector of cognostic options. Defaults are "sample count", 
"mean LCPM" and "biomolecule count". "sample count" and "mean LCPM"
are reported per group, and "biomolecule count" is the total number of biomolecules
in the biomolecule class (e_meta column).}

\item{ggplot_params}{An optional vector of strings of ggplot parameters to
the backend ggplot function. For example, c("ylab('')", "xlab('')").
Default is NULL.}

\item{interactive}{A logical argument indicating whether the plots should be
interactive or not. Interactive plots are ggplots piped to ggplotly. Default is FALSE.}

\item{path}{The base directory of the trelliscope application. Default is
Downloads.}

\item{name}{The name of the display. Default is Trelliscope}

\item{test_mode}{A logical to return a smaller trelliscope to confirm plot
and design. Default is FALSE.}

\item{test_example}{A vector of plot indices to return for test_mode. Default
is 1.}

\item{single_plot}{A TRUE/FALSE to indicate whether 1 plot (not a
trelliscope) should be returned. Default is FALSE.}

\item{...}{Additional arguments to be passed on to the trelli builder}
}
\value{
No return value, builds a trelliscope display of heatmaps that is stored in `path`
}
\description{
Specify a plot design and cognostics for the RNA-seq heatmap
  trelliscope. Data must be grouped by an e_meta column. Main_effects order
  the y-variables. All statRes data is ignored. 
  For MS/NMR data, use "trelli_abundance_heatmap".
}
\examples{
\dontshow{if (requireNamespace("pmartRdata", quietly = TRUE)) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
\dontrun{
library(pmartRdata)

omicsData_seq <- group_designation(omicsData = rnaseq_object, main_effects = c("Virus"))

# Filter low transcript counts
omicsData_seq <- applyFilt(filter_object = total_count_filter(omicsData = omicsData_seq), 
 omicsData = omicsData_seq, min_count = 15)

# Select a normalization and statistics method (options are 'edgeR', 'DESeq2', and 'voom').
# See ?difexp_seq for more details
statRes_seq <- diffexp_seq(omicsData = omicsData_seq, method = "voom")

# Generate the trelliData object
trelliData_seq2 <- as.trelliData(omicsData = omicsData_seq)
trelliData_seq4 <- as.trelliData(omicsData = omicsData_seq, statRes = statRes_seq)

# Build the RNA-seq heatmap with an omicsData object with emeta variables. 
# Generate trelliData in as.trelliData.
trelli_panel_by(trelliData = trelliData_seq2, panel = "Gene") \%>\% 
   trelli_rnaseq_heatmap(test_mode = TRUE, test_example = c(1532, 1905, 6134), path = tempdir())
   
# Users can modify the plotting function with ggplot parameters and interactivity, 
# and can also select certain cognostics.     
trelli_panel_by(trelliData = trelliData_seq4, panel = "Gene") \%>\% 
   trelli_rnaseq_heatmap(test_mode = TRUE, test_example = c(1532, 1905, 6134), 
     ggplot_params = c("ylab('')", "xlab('')"), 
     interactive = TRUE, cognostics = c("biomolecule count"), path = tempdir())  

\dontshow{closeAllConnections()}
}
\dontshow{\}) # examplesIf}
}
\author{
David Degnan, Lisa Bramer
}
