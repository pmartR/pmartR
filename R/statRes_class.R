#' Function to take raw output of `imd_anova` and create output for `statRes` object
#' 
#' Not really an interesting function.  Just exporting it so Natalie can use in a different function.
#' 
#' @param imd_out Data.frame containing the results of the imd anova call
#' @param omicsData A pmartR data object of any class, which has a `group_df` attribute that is usually created by the `group_designation()` function
#' @param comparisons the comparisons made
#' @param test_method the test method
#' @param pval_adjust pvalue adjustment method
#' @param pval_thresh p-value threshold value
#' 
#' @return the final object of class statRes
#' @export
#' 
statRes_output <- function(imd_out,omicsData,comparisons,test_method,pval_adjust,pval_thresh){
  # check that omicsData is of the appropriate class
  if(!inherits(omicsData, c("proData","pepData","lipidData", "metabData", "nmrData"))) stop("omicsData is not an object of appropriate class")
  
   # Check for group_DF attribute #
  if(is.null(attr(omicsData, "group_DF"))){
    stop("group_designation must be called in order to create a 'group_DF' attribute for omicsData.")
  }else{
    groupData <- attr(omicsData, "group_DF")
  }
  
  #Flags to determine number of significant
  imd_out_flags <- imd_out[[grep("^flag",tolower(names(imd_out)))]]
  flags <- data.matrix(imd_out_flags)
  
  
  #Add biomolecule to P_values and Flags
  meta <- imd_out$Full_results[,1]
  imd_out$Flags <- cbind.data.frame(meta, flags)
  colnames(imd_out$Flags)[1] <- colnames(imd_out$Full_results)[1]
  
  imd_out_pvals <- imd_out[[grep("^p_values",tolower(names(imd_out)))]]
  pvals <- data.matrix(imd_out_pvals)
  imd_out$P_values <- cbind.data.frame(imd_out$Full_results[,1],pvals)
  colnames(imd_out$P_values)[1] <- colnames(imd_out$Full_results)[1]
  
  #Add attributes
  attr(imd_out, "group_DF") <- groupData
  #if(is.null(comparisons)){
  comparisons <- colnames(imd_out_flags)
  #}
  attr(imd_out, "comparisons") <- comparisons
  attr(imd_out, "number_significant") <-  data.frame(Comparison=comparisons, 
                                                     Up_total = apply(flags,2,function(x){length(which(x>0))}),
                                                     Down_total = apply(flags,2,function(x){length(which(x<0))}),
                                                     Up_anova = apply(flags, 2, function(x){length(which(x == 1))}),
                                                     Down_anova = apply(flags, 2, function(x){length(which(x == -1))}),
                                                     Up_gtest = apply(flags, 2, function(x){length(which(x == 2))}),
                                                     Down_gtest = apply(flags, 2, function(x){length(which(x == -2))})
                                                     )
  attr(imd_out,"statistical_test") <- test_method
  attr(imd_out, "adjustment_method") <- pval_adjust
  attr(imd_out, "pval_thresh") <- pval_thresh
  attr(imd_out, "data_info") <- attr(omicsData, "data_info")
  class(imd_out) <- "statRes"
  return(imd_out)
}

#' statRes class.
#'
#' Class for statistical results returned by fucntions in this package.
#'
#' @name statRes-class
#' @seealso See \code{\link{imd_anova}}
#'
#' @exportClass statRes
setOldClass("statRes")

#' @export
#' @method summary statRes
summary.statRes <- function(x,...){
  cat("Type of test:",attr(x,"statistical_test"),"\n\n")
  cat("Multiple comparison adjustment:",attr(x,"adjustment_method"),"\n\n")
  cat("p-value threshold:",attr(x,"pval_thresh"),"\n\n")
  cat("Number of significant biomolecules by comparison.  Columns specify fold change direction and type of test:\n\n")
  
  table <- attr(x,"number_significant")
  names(table) <- lapply(names(table), function(name){switch(name,
                                                     Up_total = "Total:Positive",
                                                     Down_total = "Total:Negative",
                                                     Up_anova = "Positive:ANOVA",
                                                     Down_anova = "Negative:ANOVA",
                                                     Up_gtest = "Positive:G-test",
                                                     Down_gtest = "Negative:G-test",
                                                     name)
  })
  
  rownames(table) <- NULL
  
  print(table)
  
  return(invisible(list(test_type = attr(x,"statistical_test"), adjustment = attr(x,"adjustment_method"), pval_thresh = attr(x,"pval_thresh"), sig_table = table,
                        comparisons = attr(x, "comparisons"), group_DF = attr(x,"group_DF"))))
}

#' @export
#' @method print statRes
print.statRes <- function(x,...){
  
  #print.default(x,max=6,...)
  xlen <- length(x)
  xnames <- names(x)
  for(i in 1:xlen){
    cat(xnames[i],":\n")
    print(head(x[[1]]))
    x[[1]] <- NULL
    cat("\n")
  }
  
  str(x,no.list=TRUE,give.head=FALSE)
}

#' Plotting function for `statRes` objects
#' 
#' Produces plots that summarize the results contained in a `statRes` object.
#' 
#' @param x `statRes` object to be plotted, usually the result of `imd_anova`
#' @param plot_type defines which plots to be produced, options are "bar", "volcano", "gheatmap", "fcheatmap"; defaults to "bar".  See details for plot descriptions.
#' @param fc_threshold optional threshold value for fold change estimates.  Modifies the volcano plot as follows:  Vertical lines are added at (+/-)\code{fc_threshold} and all observations that have absolute fold change less than \code{abs(fc_threshold)} are colored as 'non-significant' (as specified by \code{fc_colors}).
#' @param fc_colors vector of length three with character color values interpretable by ggplot. i.e. c("orange", "black", "blue") with the values being used to color negative, non-significant, and positive fold changes respectively
#' @param stacked TRUE/FALSE for whether to stack positive and negative fold change sections in the barplot, defaults to FALSE
#' @param interactive TRUE/FALSE for whether to create an interactive plot using plotly.  Not valid for all plots.
#' @param bw_theme TRUE/FALSE for whether to apply a black-white background theme, does not affect all plots.
#' @param custom_theme a ggplot `theme` object to be applied to non-interactive plots, or those converted by plotly::ggplotly(). 
#' @param ... extra arguments passed to internal functions such as 
#' \code{\link[pmartR:plot.statRes]{sub-function pmartR::.plot.statRes()}}
#' \code{\link[pmartR:plot.statRes]{pmartR::gtest_heatmap()}}
#' 
#' @details Plot types:
#' \itemize{
#'  \item{"bar"} \code{?pmartR::statres_barplot} Bar-chart with bar heights indicating the number of significant biomolecules, grouped by test type and fold change direction.
#'  \item{"volcano"} \code{?pmartR::statres_volcano_plot} Scatterplot showing negative-log-pvalues against fold change.  Colored by statistical significance and fold change.
#'  \item{"gheatmap"} \code{?pmartR::gtest_heatmap} Heatmap with x and y axes indicating the number of nonmissing values for two groups.  Colored by number of biomolecules that fall into that combination of nonmissing values.
#'  \item{"fcheatmap"} Heatmap showing all biomolecules across comparisons, colored by fold change.
#' }
#' 
#' @export
#' @method plot statRes
#' @examples 
#' dontrun{
#' library(pmartR)
#' library(pmartRdata)
#' #Transform the data
#' 
#' #Group the data by condition
#' myproData <- group_designation(omicsData = pro_object, main_effects = c("Condition"))
#' 
#' #Apply the IMD ANOVA filter
#' imdanova_Filt <- imdanova_filter(omicsData = myproData)
#' myproData <- applyFilt(filter_object = imdanova_Filt, omicsData = myproData, min_nonmiss_anova=2)
#' 
#' #Implement the IMD ANOVA method and compuate all pairwise comparisons (i.e. leave the `comparisons` argument NULL)
#' anova_res <- imd_anova(omicsData = myproData, test_method = 'anova')
#' plot(anova_res)
#' plot(anova_res, plot_type = "volcano")
#' 
#' imd_res <- imd_anova(omicsData = myproData, test_method = 'gtest')
#' plot(imd_res)
#' plot(imd_res, plot_type = "gheatmap")
#' # using arguments of internal functions:
#' plot(imd_res, plot_type = "gheatmap", color_low = "red", color_high = "green")
#' 
#' imd_anova_res <- imd_anova(omicsData = myproData, test_method = 'comb', pval_adjust='bon')
#' plot(imd_anova_res, bw_theme = TRUE)
#' plot(imd_anova_res, plot_type = "volcano", bw_theme = TRUE)
#' 
#' }
#' 
plot.statRes <- function(x,
                         plot_type = "bar",
                         fc_threshold = NULL,
                         fc_colors = c("red", "black", "green"),
                         stacked = FALSE,
                         interactive = FALSE,
                         bw_theme = FALSE,
                         custom_theme = NULL,
                         ...) {
  .plot.statRes(
    x = x,
    plot_type = plot_type,
    fc_threshold = fc_threshold,
    fc_colors = fc_colors,
    stacked = stacked,
    interactive = interactive,
    bw_theme = bw_theme,
    custom_theme = custom_theme,
    ...
  )
}

.plot.statRes <- function(x,
                          plot_type = "bar",
                          fc_threshold = NULL,
                          fc_colors = c("red", "black", "green"),
                          stacked = FALSE,
                          interactive = FALSE,
                          bw_theme = FALSE,
                          custom_theme = NULL,
                          ...
                          ) {
  #Most plots are based on "number_significant" data frame so pull it out
  comp_df <- attr(x, "number_significant")
  
  ##--------##
  #Go through the given plot_types and remove any that aren't currently avilable
  plt_tyj <- try(match.arg(tolower(plot_type),c("bar","volcano", "gheatmap", "fcheatmap")),silent=TRUE)
  if(class(plt_tyj)=='try-error'){
    warning(paste0("Plot type '",plot_type,"' is not currently available, defaulting to bar plot."))
    plot_type <- "bar"
  }else{
    plot_type <- plt_tyj
  }
  
  #Don't make biomolecule heatmaps if there's only one comparison
  if(plot_type %in% c("fcheatmap") & nrow(comp_df)==1){
    stop("Fold change heatmaps not supported when only one comparison is being made.")
  }
  
  # specified theme parameters
  if(!is.null(custom_theme)){
    if(bw_theme) warning("Setting both bw_theme to TRUE and specifying a custom theme may cause undesirable results")
    if(!inherits(custom_theme, c("theme", "gg"))) stop("custom_theme must be a valid 'theme' object as used in ggplot")
    mytheme = custom_theme
  }
  else mytheme = ggplot2::theme(
    plot.title = ggplot2::element_text(size = 14),
    axis.title = ggplot2::element_text(size = 14),
    axis.text = ggplot2::element_text(size = 12),
    legend.title = ggplot2::element_text(size = 10),
    legend.text = ggplot2::element_text(size = 8),
    strip.text = ggplot2::element_text(size = 12)
  )
  
  # Both the volcano plot and heatmaps need a dataframe of fold changes by comparison/biomolecule
  if(plot_type %in% c("volcano", "gheatmap", "fcheatmap")){
    volcano <- make_volcano_plot_df(x)
  }
  
  #Bar plot
  if("bar"%in%plot_type){
    p <- statres_barplot(x, stacked, fc_colors)
    
    if(bw_theme) p <- p + ggplot2::theme_bw()
    
    return(p + mytheme)
  }
  
  # Volcano plot 
  else if("volcano"%in%plot_type){
    if(!attr(x, "statistical_test") %in% c("anova", "combined")){
      stop("imd_anova must have been run with test_method = 'anova' or 'combined' to make the volcano plot") 
    }

    # still returns a ggplot, even if interactive = T
    p <-
      statres_volcano_plot(
        volcano,
        data_scale = attr(x, "data_info")$data_scale,
        pval_thresh = attr(x, "pval_thresh"),
        fc_colors,
        fc_threshold,
        interactive
      )
    
    if(bw_theme) p <- p + ggplot2::theme_bw()
    
    p <- p + mytheme
    
    if(interactive) return(plotly::ggplotly(p, tooltip = c("text"))) else return(p)
  }
  
  # g-test heatmap
  else if("gheatmap" %in% plot_type){
    if(!attr(x, "statistical_test") %in% c("gtest", "combined")){
      stop("imd_anova must have been run with test_method = 'gtest' or 'combined' to make the g-test heatmap") 
    }
 
    p <-
      gtest_heatmap(
        volcano,
        pval_thresh = attr(x, "pval_thresh"),
        interactive = interactive,
        ...
      )
    
    if(!interactive){
      p <- p + mytheme
    }
  }

  # biomolecule fold change heatmap
  else if("fcheatmap"%in%plot_type){
    
    #For now just consider biomolecules significant with respect to ANOVA
    volcano <- dplyr::filter(volcano,Type=="ANOVA")
    
    volcano_sigs <- dplyr::filter(volcano,P_value<attr(x,"pval_thresh"))
    if(!(nrow(volcano_sigs)) > 0) warning("No molecules significant at the provided p-value threshold")
    colnames(volcano_sigs)[1] <- "Biomolecule"
    volcano_sigs$Biomolecule <- as.factor(volcano_sigs$Biomolecule)
    
    p <- ggplot2::ggplot(volcano_sigs, ggplot2::aes(Biomolecule, Comparison, text = paste("ID:", Biomolecule, "<br>", "Pval:", P_value))) +
          ggplot2::geom_tile(ggplot2::aes(fill = Fold_change), color = "white") +
          ggplot2::scale_fill_gradient(low = fc_colors[1], high = fc_colors[3]) +
          ggplot2::ggtitle("Average Log Fold Change") +
          mytheme
    if(interactive) return(plotly::ggplotly(p, tooltip = c("text"))) else return(p)
  }
  
  return(p)

}

#' Create a plotting dataframe for volcano plots and heatmaps.
#' 
#' A function internal to \link{pmartR::plot.statRes} which creates the 
#' dataframe necessary to construct volcano plots and heatmaps.
#' 
#' @param x `statRes` object to be plotted, usually the result of `imd_anova`
#' 
#' @returns `data.frame` object with plotting information about each biomolecule
#' such as missing counts per group, and p-values for t and g-tests.
#' 
#' @keywords internal
#' 
make_volcano_plot_df <- function(x) {
  # fold change values for volcano plot
  fc_data <-
    x$Full_results[, c(1, grep("^Fold_change", colnames(x$Full_results)))]
  colnames(fc_data) <-
    gsub(pattern = "^Fold_change_",
         replacement = "",
         x = colnames(fc_data))
  fc_data <-
    reshape2::melt(
      fc_data,
      id.vars = 1,
      variable.name = "Comparison",
      value.name = "Fold_change"
    )
  
  # fold change flags for coloring
  fc_flags <- x$Flags
  fc_flags <-
    reshape2::melt(
      fc_flags,
      id.vars = 1,
      variable.name = "Comparison",
      value.name = "Fold_change_flag"
    ) %>%
    dplyr::mutate(Fold_change_flag = as.character(Fold_change_flag))
  
  # p values for labeling and y axis in anova volcano plot
  p_data <-
    x$Full_results[c(1, grep("^P_value", colnames(x$Full_results)))]
  pvals <-
    reshape2::melt(
      p_data,
      id.vars = 1,
      variable.name = "Comparison",
      value.name = "P_value"
    )
  
  # grouping column based on test type
  if (attr(x, "statistical_test") == "combined") {
    pvals$Type <- "G-test"
    pvals$Type[grep(pattern = "^P_value_T_", x = pvals$Comparison)] <-
      "ANOVA"
  } else if (attr(x, "statistical_test") == "gtest") {
    pvals$Type <- "G-test"
  } else if (attr(x, "statistical_test") == "anova") {
    pvals$Type <- "ANOVA"
  }
  
  levels(pvals$Comparison) <-
    gsub(pattern = "^P_value_G_",
         replacement = "",
         levels(pvals$Comparison))
  levels(pvals$Comparison) <-
    gsub(pattern = "^P_value_T_",
         replacement = "",
         levels(pvals$Comparison))
  
  volcano <-
    merge(merge(fc_data, pvals, all = TRUE), fc_flags, all = TRUE)
  
  # levels of comparison now of the form 'GROUPNAME_X vs GROUPNAME_Y'
  levels(volcano$Comparison) <-
    gsub(pattern = "_vs_",
         replacement = " vs ",
         levels(volcano$Comparison))
  
  # create counts for gtest plot (number present in each group)
  if (attr(x, "statistical_test") %in% c("gtest", "combined")) {
    counts <-
      x$Full_results[c(1, grep("^Count_", colnames(x$Full_results)))]
    
    # trim column names so they are just group names
    colnames(counts) <-
      gsub("^Count_", replacement = "", colnames(counts))
    
    counts_df <- data.frame()
    for (comp in as.character(unique(volcano$Comparison))) {
      #create a vector of the two group names being compared
      groups = strsplit(comp, " vs ")[[1]]
      gsize_1 = nrow(attr(x, "group_DF") %>% dplyr::filter(Group == groups[1]))
      gsize_2 = nrow(attr(x, "group_DF") %>% dplyr::filter(Group == groups[2]))
      
      # will contain ID column and count column corresponding to the two groups
      temp_df <-
        counts[c(1, which(colnames(counts) == groups[1]), which(colnames(counts) == groups[2]))]
      temp_df$Comparison <- comp
      
      # rename the columns to something static so they can be rbind-ed
      colnames(temp_df)[which(colnames(temp_df) %in% groups)] <-
        c("Count_First_Group", "Count_Second_Group")
      
      # store proportion of nonmissing to color g-test values in volcano plot
      temp_df$Prop_First_Group <- temp_df$Count_First_Group / gsize_1
      temp_df$Prop_Second_Group <-
        temp_df$Count_Second_Group / gsize_2
      
      counts_df <- rbind(counts_df, temp_df)
    }
    
    # should automatically left join by ID AND Comparison
    suppressWarnings(volcano <-
                       volcano %>% dplyr::left_join(counts_df))
  }
  
  return(volcano)
}

#' Fold change barplots for statres objects
#' 
#' Plots a bar-chart with bar heights indicating the number of significant 
#' biomolecules, grouped by test type and fold change direction.
#' 
#' @param x,stacked,fc_colors
#' passed from \code{\link[pmartR:plot.statRes]{pmartR::plot.statRes()}}
#' 
#' @return `ggplot` object - barplot.
#' 
#' @keywords internal
#' 
statres_barplot <- function(x, stacked = FALSE, 
                            fc_colors = c("red", "black", "green")) {
  comp_df <- attr(x, "number_significant")
  
  comp_df_melt <- reshape2::melt(comp_df,id.vars="Comparison",value.name="Count",variable.name="Direction")
  levels(comp_df_melt$Comparison) <- gsub(pattern = "_",replacement = " ",levels(comp_df_melt$Comparison))
  
  #Bar plots side-by-side, both going up
  #ggplot(data=comp_df_melt,aes(Comparison,Count,fill=Direction))+geom_bar(stat='identity',position='dodge')
  
  ##Up direction is positive, down direction is negative
  if(!stacked) comp_df_melt[grep("Down", comp_df_melt$Direction),]$Count <- (-comp_df_melt[grep("Down", comp_df_melt$Direction),]$Count)
  
  # add whichtest, and posneg columns used for plot grouping and label adjustment
  comp_df_melt <- comp_df_melt %>% 
    dplyr::mutate(whichtest = ifelse(grepl("anova", Direction), "anova", ifelse(grepl("gtest", Direction), "gtest", "total")),
                  posneg = ifelse(grepl("Up", Direction), "Positive", "Negative")) %>%
    dplyr::arrange(desc(posneg))
  
  # get only anova or only g-test rows if user did not specify combined
  if(attr(x, "statistical_test") %in% c("anova", "gtest")){
    comp_df_melt <- comp_df_melt %>%
      dplyr::filter(whichtest == "total") %>%
      dplyr::mutate(whichtest = attr(x, "statistical_test"))
  }
  
  p <- ggplot2::ggplot(data = comp_df_melt, ggplot2::aes(Comparison, Count)) +
    ggplot2::geom_bar(ggplot2::aes(x = whichtest, fill = posneg, group = whichtest), stat =
               'identity') +
    ggplot2::geom_text(ggplot2::aes(x = whichtest, label = ifelse(abs(Count) > 0, abs(Count), "")),
              position = ggplot2::position_stack(vjust = 0.5),
              size = 3) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = 0), colour = 'gray50') +
    ggplot2::scale_fill_manual(
      values = c(fc_colors[1], fc_colors[3]),
      labels = c("Negative", "Positive"),
      name = "Fold Change Sign"
    ) +
    ggplot2::facet_wrap( ~ Comparison) +
    ggplot2::xlab("Statistical test, by group comparison") + 
    ggplot2::ylab("Count of Biomolecules") +
    ggplot2::ggtitle("Number of DE Biomolecules Between Groups")
  
    return(p)
}

#' Plot a heatmap for the g-test results of imd-anova 
#' 
#' Plots a heatmap showing bins for combinations of # of biomolecules present
#' across groups.  Bins are colored by the number of biomolecules falling into
#' each bin, and have indicators for significance by g-test.
#' 
#' @param volcano `data.frame` produced by \link{pmartR::make_volcano_plot_df}
#' @param show_sig Boolean whether to show the visual indicator that a certain 
#' bin combination is significant by the g-test
#' @param count_text_size In non-interactive mode, the size of the text 
#' indicating the number of biomolecules to be displayed inside the tile.  
#' Defaults to 3.  Setting this to zero causes no text to be shown.
#' @param interactive passed from \code{\link[pmartR:plot.statRes]{pmartR::plot.statRes()}}.  
#' If T, will build a plotly version of the plot.  Defaults to FALSE.
#' @param color_low A character string specifying the color of the gradient for
#'   low count values.
#' @param color_high A character string specifying the color of the gradient for
#'   high count values.
#' @param plotly_layout A list of arguments, not including the plot, to be 
#' passed to plotly::layout if interactive = T.
#' 
#' 
#' @return `ggplot or plotly` object, depending on if \code{interactive} is set 
#' to TRUE or FALSE respectively. A g-test heatmap
#' 
#' @keywords internal
#' 
gtest_heatmap <-
  function(volcano,
           pval_thresh = 0.05,
           show_sig = TRUE,
           count_text_size = 3,
           interactive = FALSE,
           color_low = NULL,
           color_high = NULL,
           plotly_layout = NULL,
           ...) {
    
  temp_data_gtest <- volcano %>% 
    dplyr::filter(Type == "G-test") %>% 
    dplyr::mutate(
      Fold_change_flag = dplyr::case_when(
        Prop_First_Group > Prop_Second_Group &
          P_value <= pval_thresh ~  "2",
        Prop_First_Group < Prop_Second_Group &
          P_value <= pval_thresh ~ "-2",
        is.na(Fold_change) ~ "0",
        TRUE ~ Fold_change_flag
      )
    )
  
  # summarized data frame of # of biomolecules per count combination.
  gtest_counts <- temp_data_gtest %>%
    dplyr::group_by(Count_First_Group, Count_Second_Group, Comparison) %>%
    dplyr::summarise(n = dplyr::n(), sig = all(P_value < 0.05)) %>%
    dplyr::mutate(
      Count_First_Group = as.character(Count_First_Group),
      Count_Second_Group = as.character(Count_Second_Group)
    )
  
  g1_counts <- unique(as.numeric(gtest_counts$Count_First_Group))
  g2_counts <- unique(as.numeric(gtest_counts$Count_Second_Group))
  
  all_counts <- expand.grid(
    as.character(g1_counts), 
    as.character(g2_counts), 
    unique(gtest_counts$Comparison), stringsAsFactors = F) %>% 
    `colnames<-`(c("Count_First_Group", "Count_Second_Group", "Comparison"))
  
  gtest_counts <- all_counts %>% dplyr::left_join(gtest_counts)
  
  if(!interactive) {
    p <- ggplot2::ggplot(gtest_counts) +
      ggplot2::theme_minimal() + 
      ggplot2::geom_tile(ggplot2::aes(Count_First_Group, Count_Second_Group, fill = n), color = "black")
    
    if(show_sig) {
      p <- p + ggplot2::geom_point(
        data = gtest_counts %>% dplyr::filter(sig),
        ggplot2::aes(Count_First_Group, Count_Second_Group, shape = "1"),
        fill = "white"
      ) +
        ggplot2::scale_shape_manual(name = "Statistically significant",
                           labels = "",
                           values = 21)
    }
    
    if(count_text_size > 0) {
      p <- p + ggplot2::geom_text(
        ggplot2::aes(Count_First_Group, Count_Second_Group, label = n),
        nudge_x = -0.5, nudge_y = 0.5, hjust = -0.1, vjust = 1.5,
        color = "white", size = count_text_size
      )
    }
      
    p <- p +
      ggplot2::facet_wrap(~ Comparison) + 
      ggplot2::scale_fill_gradient(
        name = "Number of biomolecules \nin this bin",
        low = if (is.null(color_low)) "#132B43" else color_low,
        high = if (is.null(color_high)) "#56B1F7" else color_high
      ) + 
      ggplot2::xlab("Nonmissing (first group)") + 
      ggplot2::ylab("Nonmissing (second group)")
    
  } else {
    comps <- unique(gtest_counts$Comparison)
    subplot_list <- list()
    limits = range(gtest_counts$n)
    # legend entry will not appear for a single plot, so putting info in title
    subtext = if(length(comps) == 1) "\n(star indicates statistical significance)" else ""
    
    for(i in 1:length(comps)){
      data <- gtest_counts %>% dplyr::filter(Comparison == comps[i])
      
      p <- plotly::plot_ly() %>% 
        plotly::add_trace(data = data,
                  x =  ~ as.numeric(Count_First_Group),
                  y =  ~ as.numeric(Count_Second_Group),
                  z = ~ n,
                  type = "heatmap",
                  hoverinfo = 'text',
                  text = ~sprintf("%s biomolecules in this bin.", n),
                  colors = grDevices::colorRamp(
                    c(if (is.null(color_low)) "#132B43" else color_low,
                      if (is.null(color_high)) "#56B1F7" else color_high)
                  ),
                  showscale = i == length(comps),
                  xgap = 0.6, ygap = 0.6
        ) %>% 
        plotly::add_trace(
          data = data %>% filter(sig),
          x =  ~ as.numeric(Count_First_Group),
          y =  ~ as.numeric(Count_Second_Group),
          type = 'scatter',
          mode = "markers",
          showlegend = i == length(comps),
          name = "Statistically significant",
          marker = list(
            symbol = "star",
            color = "white",
            line = list(color = "black", width = 0.5)
          ),
          inherit = FALSE
        ) %>% 
        plotly::add_annotations(
          text = comps[[i]],
          x = 0.5,
          y = 1,
          xref = "paper",
          yref = "paper",
          showarrow = FALSE,
          yanchor = "bottom",
          font = list(size = 14)
        )
      
      if(i == length(comps)) {
        p <- p %>% plotly::colorbar(title = "Count biomolecules \nin this bin.", limits = limits)
      }
      
      
      if(is.null(plotly_layout)) {
        p <- p %>% 
          plotly::layout(
            plot_bgcolor = 'grey',
            title = list(
              text = sprintf("Count nonmissing in each group.%s", subtext),
              font = list(size = 14),
              y = 0.95, yanchor = "top"
            ),
            margin = list(t = 65),
            xaxis = list(tickvals = g1_counts, zeroline = F), 
            yaxis = list(tickvals = g2_counts, zeroline = F)
          )   
      }
      else {
        p <- do.call(plotly::layout, c(list(p), plotly_layout))
      }
      
      subplot_list[[length(subplot_list) + 1]] <- p 
    }
    
    p <- plotly::subplot(subplot_list, shareY = T) %>%
      plotly::layout(
        xaxis = list(title = "Nonmissing (first group)"),
        yaxis = list(title = "Nonmissing (second group)")
      )
  }
  
  return(p)
}

#' Volcano plot for the anova results of imd-anova
#' 
#' Plots a volcano plot showing negative log10 p-values on the y axis and fold
#' change on the x axis. Each point is colored by fold change direction and 
#' whether or not it was significant by ANOVA.
#' 
#' @param volcano `data.frame` produced by \link{pmartR::make_volcano_plot_df}
#' @param pval_thresh alpha level to determine significance.  Any values that
#' are significant at this level will be colored based on fc_colors.
#' @param data_scale One of c("log2","log","log10"), for labeling purposes.
#' @param fc_colors,fc_threshold,interactive 
#' passed from \code{\link[pmartR:plot.statRes]{pmartR::plot.statRes()}}
#' 
#' @return `ggplot` object.  A volcano plot.
#' 
#' @keywords internal
#' 
statres_volcano_plot <-
  function(volcano,
           data_scale,
           pval_thresh = 0.05,
           fc_colors = c("red", "black", "green"),
           fc_threshold = NULL,
           interactive = F) {
    
  # color vector which assigns black to gtest flag values (-2, 2)
  cols_anova <-
    c(
      "-2" = fc_colors[2],
      "-1" = fc_colors[1],
      "0" = fc_colors[2],
      "1" = fc_colors[3],
      "2" = fc_colors[2]
    )
  
  # temp data with rows only for ANOVA
  temp_data_anova <- volcano %>%
    dplyr::filter(Type == "ANOVA") %>%
    dplyr::mutate(
      Fold_change_flag = dplyr::case_when(
        is.na(Fold_change) |
          abs(Fold_change) < ifelse(rlang::is_empty(fc_threshold), 0, abs(fc_threshold)) ~ "0",
        Fold_change > 0 &
          P_value <= pval_thresh ~  "1",
        Fold_change < 0 &
          P_value <= pval_thresh ~ "-1",
        TRUE ~ Fold_change_flag
      )
    )
  
  # interactive plots need manual text applied to prepare for ggplotly conversion
  if (interactive) {
    p <-
      ggplot2::ggplot(temp_data_anova, ggplot2::aes(
        Fold_change,
        -log(P_value, base = 10),
        text = paste(
          "ID:",
          !!rlang::sym(colnames(volcano)[1]),
          "<br>",
          "Pval:",
          round(P_value, 4)
        )
      ))
  }
  else {
    p <-
      ggplot2::ggplot(data = temp_data_anova, ggplot2::aes(Fold_change, -log(P_value, base = 10))) 
  }
  
  # draw vertical lines at +- fc threshold
  if(!rlang::is_empty(fc_threshold)){
    p <- p + 
      ggplot2::geom_vline(ggplot2::aes(xintercept=abs(fc_threshold)), lty = 2) + 
      ggplot2::geom_vline(ggplot2::aes(xintercept = abs(fc_threshold)*(-1)), lty = 2)
  }
  
  p <- p +
    ggplot2::geom_point(ggplot2::aes(color = Fold_change_flag), shape = 1) +
    ggplot2::facet_wrap( ~ Comparison) +
    ggplot2::ylab("-log[10](p-value)") + 
    ggplot2::xlab(sprintf("Fold-change (%s)", data_scale)) +
    ggplot2::scale_color_manual(
      values = cols_anova,
      name = "Fold Change",
      labels = c("Neg(Anova)", "0", "Pos(Anova)"),
      breaks = c("-1", "0", "1")
    )
  
  return(p)
}