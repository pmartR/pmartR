## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)
library(pmartRdata)
library(pmartR)
library(data.table)
library(DT)
library(dplyr)
library(patchwork)

## ---- echo = FALSE------------------------------------------------------------
data.table(
  `Input Data` = c(rep("Abundance", 3), "Missingness", rep("Fold Change", 4)),
  `Plot Type` = c("Boxplot", "Histogram", "Heatmap", "Barplot", "Barplot", "Boxplot", "Volcano", "Heatmap"),
  `Uses e_data` = c("X", "X", "", "X", rep("", 4)),
  `Uses omicsData` = c(rep("X", 4), rep("", 4)),
  `Uses statRes` = c(rep("", 3), rep("X", 5)),
  `Requires Biomolecule Class (e_meta)` = c("", "", "X", "", "", "X", "X", "X"),
  `Main Effects Used` = c("X", "", rep("X", 4), "", "")
) %>% knitr::kable()

## ---- echo = FALSE------------------------------------------------------------
knitr::include_graphics("TrelliData_Plotting_Options.png")

## ---- echo = F----------------------------------------------------------------
datatable(head(pmartRdata::lipid_pos_edata), options = list(`scrollX` = T))

## -----------------------------------------------------------------------------
trelliData1 <- as.trelliData.edata(
  e_data = pmartRdata::lipid_pos_edata,
  edata_cname = "Lipid",
  omics_type = "lipidData",
  data_scale_original = "abundance",
  data_scale = "log2",
  normalization_fun = "global",
  normalization_params = list(subset_fn = "all", norm_fn = "median", apply_norm = TRUE,
    backtransform = TRUE)
)

## -----------------------------------------------------------------------------
summary(trelliData1) 

## -----------------------------------------------------------------------------
summary(trelliData1 %>% trelli_panel_by("Lipid")) 

## ---- echo = FALSE------------------------------------------------------------
data.table(
  "Parameter Name" = c("cognostics", "ggplot_params", "interactive", "path", 
                       "name", "test_mode", "test_example", "single_plot"),
  "Description" = c("Set the specific cognostics of the trelliscope. Varies per plotting function.", 
                    "Pass parameters to the ggplot functions as a list of strings",
                    "Indicate whether plots should be interactive or not", 
                    "The path where the trelliscope will be outputted to",
                    "Name of the trelliscope display",
                    "Indicate whether the trelliscope should be subsetted to a few panels or not",
                    "The panels to subset the trelliscope to if test_mode is true",
                    "Output a single plot instead of a trelliscope display")
) %>% knitr::kable()

## ---- echo = T, eval = F------------------------------------------------------
#  trelli_abundance_boxplot(
#    trelli_panel_by(trelliData1, "Lipid"),
#    cognostics = c("n", "mean", "median", "sd"),
#    interactive = TRUE,
#    include_points = TRUE,
#    name = "Trelliscope",
#    test_mode = TRUE,
#    test_example = 3
#  )

## -----------------------------------------------------------------------------
# Panel by Mass Tag ID  
MassGroups <- trelli_panel_by(trelliData1, "Lipid")
# Panel by Sample
SampleGroups <- trelliData1 %>% trelli_panel_by("Sample") 
# Create an example boxplot 
Abun_Box_Edata <- trelli_abundance_boxplot(MassGroups, single_plot = TRUE, test_example = 3)
# Make an abundance boxplot without the points 
Abun_Box_Sample <- trelli_abundance_boxplot(SampleGroups, include_points = F, single_plot = T,
                                            ggplot_params = "scale_fill_manual(values = 'forestgreen')")
# Use patchwork to put plots together
Abun_Box_Edata + Abun_Box_Sample

## ---- echo = T----------------------------------------------------------------
trelli_abundance_histogram(MassGroups, single_plot = TRUE, test_example = 3)

## -----------------------------------------------------------------------------
# Pull from the peptide dataset
pep_trelliData <- as.trelliData.edata(
  e_data = pmartRdata::pep_edata,
  edata_cname = "Peptide",
  omics_type = "pepData"
)
# Create an example bar plot 
Miss_Bar_Edata <- trelli_missingness_bar(trelli_panel_by(pep_trelliData, "Peptide"), 
                                         single_plot = TRUE, test_example = 3,
                                         ggplot_params = "ggtitle('Biomolecule')")
# Make a missingness barplot 
Miss_Bar_Sample <- trelli_missingness_bar(trelli_panel_by(pep_trelliData, "Sample"), 
                                          include_points = F, 
                                          proportion = FALSE, single_plot = T,
                                          ggplot_params = "ggtitle('Sample')")
# Put plots together with patchwork
Miss_Bar_Edata + Miss_Bar_Sample

## -----------------------------------------------------------------------------
# Pull the lipids data frame 
lipids <- pmartRdata::lipid_pos_object
# Extract and add lipid family information
lipids$e_meta$Family <- lipids$e_meta$Lipid %>% 
  strsplit("(", fixed = T) %>% 
  lapply(function(x) {head(x, 1)}) %>% 
  unlist()
# Log transform the edata file 
lipids <- edata_transform(omicsData = lipids, data_scale = "log2")
# Set the group designation main effects 
lipids <- group_designation(omicsData = lipids, main_effects = c("Virus"))
# Filter the data to run the imd_anova
imdanova_Filt <- imdanova_filter(omicsData = lipids)
lipids <- applyFilt(filter_object = imdanova_Filt, omicsData = lipids, min_nonmiss_anova=2)
# Normalize the data. You may use the Kruskal-Wallis test from normRes_tests to 
# confirm that your normalization choice does not induce nor remove significance. 
lipids <- normalize_global(lipids, subset_fn = "rip", norm_fn = "median", 
                             apply_norm = TRUE, backtransform = TRUE)
# Create trelliData object
trelliData2 <- as.trelliData(omicsData = lipids)

## -----------------------------------------------------------------------------
summary(trelliData2)

## -----------------------------------------------------------------------------
trelliData2 %>% 
  trelli_panel_by("Lipid") %>% 
  trelli_abundance_boxplot(test_example = 3, single_plot = T, ggplot_params = "xlab('Viral Strain')")

## -----------------------------------------------------------------------------
trelliData2 %>% 
  trelli_panel_by("Family") %>% 
  trelli_abundance_boxplot(test_example = 6, single_plot = T)

## -----------------------------------------------------------------------------
trelliData2 %>% 
  trelli_panel_by("Family") %>% 
  trelli_missingness_bar(test_example = 3, single_plot = T)

## -----------------------------------------------------------------------------
trelliData2 %>% 
  trelli_panel_by("Family") %>% 
  trelli_abundance_heatmap(test_example = 2, single_plot = T, ggplot_params = "coord_flip()")

## -----------------------------------------------------------------------------
# Since nothing is missing from this dataset, there is no need to run the G-test. 
lipidStats <- imd_anova(lipids, test_method = "anova")
# Pass that dataframe to the trelliData object builder
trelliData3 <- as.trelliData(statRes = lipidStats)
# Run the summary function
summary(trelliData3)

## -----------------------------------------------------------------------------
trelliData3 %>%
  trelli_panel_by("Lipid") %>%
  trelli_foldchange_bar(p_value_thresh = 0.05, p_value_test = TRUE, 
                        test_example = 3, single_plot = TRUE)

## -----------------------------------------------------------------------------
# Build a trelliData object with both omicsData and statRes
trelliData4 <- as.trelliData(omicsData = lipids, statRes = lipidStats)
summary(trelliData4)

## -----------------------------------------------------------------------------
trelliData4 %>%
  trelli_panel_by("Family") %>%
  trelli_foldchange_boxplot(
    single_plot = TRUE,
    p_value_test = TRUE,
    test_example = 4
  )

## -----------------------------------------------------------------------------
trelliData4 %>%
  trelli_panel_by("Family") %>%
  trelli_foldchange_heatmap(
    single_plot = TRUE
  )

## -----------------------------------------------------------------------------
trelliData4 %>%
  trelli_panel_by("Family") %>%
  trelli_foldchange_volcano(
    comparison = "StrainA_vs_StrainB",
    test_example = 3,
    p_value_test = TRUE,
    single_plot = TRUE
  )

## -----------------------------------------------------------------------------
trelliData2 %>% 
  trelli_panel_by("Lipid") %>% 
  trelli_abundance_boxplot(test_example = 3, single_plot = T,
    ggplot_params = c("ggtitle('CYC Human')", 
                      "ylab('Log Adjusted Abundance')",
                      "xlab('')", "theme_classic()", 
                      "theme(plot.title = ggplot2::element_text(hjust = 0.5))",
                      "theme(legend.position = 'none')")                       
  )

## ---- echo = FALSE------------------------------------------------------------
data.table(
  Function = c("xlab('Name')", "ylab('Name')", "ggtitle('Name')", 
               "xlim(c(0,1))", "ylim(c(0,1))",
               "coord_flip()",
               "theme(axis.title.x = ggplot2::element_text(size=16))",
               "theme(axis.title.y = ggplot2::element_text(size=16))",
               "theme(axis.text.x = element_text(size=12))",
               "theme(axis.text.y = element_text(size=12))",
               "theme(axis.text.x = element_text(angle=90))",
               "theme(axis.text.y = element_text(angle=90))",
               "theme(plot.title = element_text(size=30))"
               ),
  Purpose = c("Rename x-axis label", "Rename y-axis label", "Rename plot",
              "Set the x-axis limits", "Set the y-axis limits",
              "Flip x- and y-axis",
              "Resize x-axis title", "Resize y-axis title",
              "Resize x-axis ticks", "Resize y-axis ticks",
              "Rotate x-axis ticks", "Rotate y-axis ticks",
              "Change plot title size"
              )
  
) %>% knitr::kable()

## -----------------------------------------------------------------------------
trelliData2 %>% 
  trelli_panel_by("Lipid") %>% 
  trelli_abundance_boxplot(test_example = 3, single_plot = T, interactive = T)

## ---- echo = FALSE------------------------------------------------------------
data.table(
 Functions = c("trelli_abundance_boxplot", "trelli_abundance_histogram", 
               "trelli_abundance_heatmap", "trelli_missingness_bar",
               "trelli_foldchange_bar", "trelli_foldchange_boxplot",
               "trelli_foldchange_heatmap", "trelli_foldchange_volcano"),
 Cognostics = c("n, mean, median, sd, skew, p_value, fold_change", 
                "n, mean, median, sd, skew, p_value, fold_change",
                "n, mean, median, sd, skew",
                "n, proportion",
                "fold_change, p_value",
                "n, median, mean, sd",
                "n, median, mean, sd",
                "n")
) %>% knitr::kable()

## -----------------------------------------------------------------------------
attr(trelliData4, "panel_by_options")

