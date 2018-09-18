#' Function to take raw output of `imd_anova` and create output for `statRes` object
#' 
#' Not really an interesting function.  Just exporting it so Natalie can use in a different function.
#' 
#' @param imd_out Data.frame containing the results of the imd anova call
#' @param omicsData A pmartR data object of any class, which has a `group_df` attribute that is usually created by the `group_designation()` function
#' @param comparisons the comparisons made
#' @param test_method the test method
#' @param pval_ajust pvalue adjustment method
#' @param pval_thres p-value threshold value
#' 
#' @return the final object of class statRes
#' @export
#' 
statRes_output <- function(imd_out,omicsData,comparisons,test_method,pval_adjust,pval_thresh){
  # check that omicsData is of the appropriate class
  if(!inherits(omicsData, c("proData","pepData","lipidData", "metabData"))) stop("omicsData is not an object of appropriate class")
  
   # Check for group_DF attribute #
  if(is.null(attr(omicsData, "group_DF"))){
    stop("group_designation must be called in order to create a 'group_DF' attribute for omicsData.")
  }else{
    groupData <- attr(omicsData, "group_DF")
  }
  
  #Flags to determine number of significant
  imd_out_flags <- imd_out[[grep("flag",tolower(names(imd_out)))]]
  flags <- data.matrix(imd_out_flags)
  
  
  #Add biomolecule to P_values and Flags
  meta <- imd_out$Full_results[,1]
  imd_out$Flags <- cbind.data.frame(meta, flags)
  colnames(imd_out$Flags)[1] <- colnames(imd_out$Full_results)[1]
  
  imd_out_pvals <- imd_out[[grep("p_values",tolower(names(imd_out)))]]
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
  cat("Number of significant biomolecules by comparison:\n\n")
  print(attr(x,"number_significant"))
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
#' @param plot_type defines which plots to be produced, options are "bar", "heatmap", "volcano"; defaults to "bar"
#' @param fc_threshold optional threshold value for fold change estimates that's added to the volcano plot
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
#' 
#' imd_res <- imd_anova(omicsData = myproData, test_method = 'gtest')
#' plot(imd_res)
#' 
#' imd_anova_res <- imd_anova(omicsData = myproData, test_method = 'comb', pval_adjust='bon')
#' plot(imd_anova_res)
#' }
#' 
plot.statRes <- function(x, plot_type = "bar", fc_threshold = NULL, stacked = FALSE, ...){
  
  #For now require ggplot2, consider adding base graphics option too
  if(!require(ggplot2)){
    stop("ggplot2 is required, base graphics aren't available yet")
  }
  
  #Most plots are based on "number_significant" data frame so pull it out
  comp_df <- attr(x,"number_significant")
  
  ##--------##
  #Go through the given plot_types and remove any that aren't currently avilable
  plt_tyj <- try(match.arg(tolower(plot_type),c("bar","volcano","heatmap")),silent=TRUE)
  if(class(plt_tyj)=='try-error'){
    warning(paste0("Plot type '",plot_type,"' is not currently available, defaulting to bar plot."))
    plot_type <- "bar"
  }else{
    plot_type <- plt_tyj
  }
  
  #Don't make heatmaps if there's only one comparison
  if("heatmap"%in%plot_type & nrow(comp_df)==1){
    stop("Heatmaps not supported when only one comparison is being made.")
  }
  
  #Bar plot
  if("bar"%in%plot_type){
    comp_df_melt <- reshape2::melt(comp_df,id.vars="Comparison",value.name="Count",variable.name="Direction")
    levels(comp_df_melt$Comparison) <- gsub(pattern = "_",replacement = " ",levels(comp_df_melt$Comparison))
    
    #Bar plots side-by-side, both going up
    #all_plts[[1]] <- ggplot(data=comp_df_melt,aes(Comparison,Count,fill=Direction))+geom_bar(stat='identity',position='dodge')
    
    #To turn the up direction to green, and down direction to red
    pal <- RColorBrewer::brewer.pal(9, "Set1")
    
    ##Up direction is positive, down direction is negative
    if(!stacked) comp_df_melt[grep("Down", comp_df_melt$Direction),]$Count <- (-comp_df_melt[grep("Down", comp_df_melt$Direction),]$Count)
    
    # add whichtest, posneg, and adjust columns used for plot grouping and label adjustment
    comp_df_melt <- comp_df_melt %>% 
      dplyr::mutate(whichtest = ifelse(grepl("anova", Direction), "anova", ifelse(grepl("gtest", Direction), "gtest", "total")),
                    posneg = ifelse(grepl("Up", Direction), "Positive", "Negative")) %>%
      dplyr::arrange(desc(posneg))
    
    if(attr(x, "statistical_test") %in% c("anova", "gtest")){
      comp_df_melt <- comp_df_melt %>%
        dplyr::filter(whichtest == "total") %>%
        mutate(whichtest = attr(x, "statistical_test"))
    }
      
    p <- ggplot(data=comp_df_melt,aes(Comparison,Count)) + 
            geom_bar(aes(x = whichtest, fill = posneg, group = whichtest),stat='identity') + 
            geom_text(aes(x = whichtest, label = abs(Count)), position = position_stack(vjust = 0.5), size = 3) +
            geom_hline(aes(yintercept=0),colour='gray50') +
            theme(axis.text.x=element_text(angle=90,vjust=0.5)) +
            scale_fill_manual(values=c("red", "green"), labels = c("Negative", "Positive"), name = "Fold Change Sign") +
            facet_wrap(~Comparison) + 
            ggtitle("Number of DE Metabolites Between Groups") +
            theme(plot.title = element_text(hjust = 0.5), legend.text = element_text(size = 6), legend.title = element_text(size = 6))
  
    return(p)
  }
  
  #Both the volcano plot and heatmap need a dataframe of fold changes by comparison/biomolecule
  if("volcano" %in%plot_type | "heatmap"%in%plot_type){
    fc_data <- x$Full_results[,c(1,grep("Fold_change",colnames(x$Full_results)))]
    colnames(fc_data) <- gsub(pattern = "Fold_change_",replacement = "",x = colnames(fc_data))
    fc_data <- reshape2::melt(fc_data,id.vars=1,variable.name="Comparison",value.name="Fold_change")
    
    fc_flags <- x$Full_results[,c(1,grep("Flag_",colnames(x$Full_results)))]
    colnames(fc_flags) <- gsub(pattern = "Flag_",replacement = "",x = colnames(fc_flags))
    fc_flags <- reshape2::melt(fc_flags,id.vars=1,variable.name="Comparison",value.name="Fold_change_flag") %>%
      dplyr::mutate(Fold_change_flag = as.character(Fold_change_flag))
    
    p_data <- x$Full_results[,c(1,grep("P_value",colnames(x$Full_results)))]
    pvals <- reshape2::melt(p_data,id.vars=1,variable.name="Comparison",value.name="P_value")
    
    if(attr(x,"statistical_test")=="combined"){
      pvals$Type <- "G-test"
      pvals$Type[grep(pattern="P_value_T_",x=pvals$Comparison)] <- "ANOVA"
    }else if(attr(x,"statistical_test")=="gtest"){
      pvals$Type <- "G-test"
    }else{
      pvals$Type <- "ANOVA"
    }
    
    levels(pvals$Comparison) <- gsub(pattern="P_value_G_",replacement = "",levels(pvals$Comparison))
    levels(pvals$Comparison) <- gsub(pattern="P_value_T_",replacement = "",levels(pvals$Comparison))
    
    volcano <- merge(merge(fc_data,pvals,all=TRUE), fc_flags, all = TRUE)
    levels(volcano$Comparison) <- gsub(pattern = "_",replacement = " ",levels(volcano$Comparison))
    
    #Remove NAs to avoid ggplot2 warning
    to_rm <- which(apply(volcano,1,function(x){any(is.na(x))}))
    if(length(to_rm)>0)
      volcano <- volcano[-to_rm,]
  }
  
  #Volcano plot 
  if("volcano"%in%plot_type){
    # specify colorscale for flags
    cols <- c("-2" = "red", "-1" = "red", "0" = "black", "1" = "green", "2" = "green")
    
    p <- ggplot(data=volcano,aes(Fold_change,-log(P_value,base=10)))+
      geom_point(aes(color = Fold_change_flag))+
      facet_grid(Type~Comparison)+
      ylab(expression(-log[10](p-value)))+xlab("Fold-change") +
      scale_color_manual(values = cols)
    
    # draw a line if threshold specified
    if(!is.null(fc_threshold)){
      p <- p + geom_hline(aes(yintercept=fc_threshold))
    }
    
    return(p)
  }
  
  #Heatmap
  if("heatmap"%in%plot_type){
    pal <- colorRampPalette(colors=rev(brewer.pal(11, "RdYlGn")))
    
    #For now just consider biomolecules significant with respect to ANOVA
    volcano <- dplyr::filter(volcano,Type=="ANOVA")
    
    volcano_sigs <- dplyr::filter(volcano,P_value<attr(x,"pval_thresh"))
    colnames(volcano_sigs)[1] <- "Biomolecule"
    volcano_sigs$Biomolecule <- as.factor(volcano_sigs$Biomolecule)
    
    p <- ggplot(volcano_sigs, aes(Biomolecule, Comparison)) +
          geom_tile(aes(fill = Fold_change), color = "white") +
          scale_fill_gradient(low = "green", high = "red") +
          theme(axis.text.x = element_text(angle = 90, hjust = 0.5)) +
          ggtitle("Average Log Fold Change") +
          theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_blank())
      # ggplot(data=volcano_sigs,aes(Comparison,Biomolecule)) +
      # geom_tile(aes(fill=Fold_change)) +
      # theme(axis.text.x=element_text(angle=90,vjust=0.5))+
      # scale_fill_brewer(palette = "RdYlGn")
    return(p)
  }
  #invisible(all_plts)
}