#' Function to take raw output of `imd_anova` and create output for `statRes` object
#' 
#' Not really an interesting function.  Just exporting it so Natalie can use in a different function.
#' 
#' @param imd_out Data.frame containing the results of the imd anova call
#' @param groupData the groupData object
#' @param comparisons the comparisons made
#' @param test_method the test method
#' @param pval_ajust pvalue adjustment method
#' @param pval_thres p-value threshold value
#' 
#' @return the final object of class statRes
#' @export
#' 
statRes_output <- function(imd_out,groupData,comparisons,test_method,pval_adjust,pval_thresh){
  
  #Flags to determine number of significant
  imd_out_flags <- imd_out[[grep("flag",tolower(names(imd_out)))]]
  flags <- data.matrix(imd_out_flags[, grep("flag",tolower(names(imd_out_flags)))])
  
  #Add biomolecule to P_values and Flags
  meta <- imd_out$Full_results[,1]
  imd_out$Flags <- cbind.data.frame(meta, flags)
  colnames(imd_out$Flags)[1] <- colnames(imd_out$Full_results)[1]
  
  imd_out_pvals <- imd_out[[grep("p_values",tolower(names(imd_out)))]]
  pvals <- data.matrix(imd_out_pvals[, grep("pvals",tolower(names(imd_out_pvals)))])
  imd_out$P_values <- cbind.data.frame(imd_out$Full_results[,1],pvals)
  colnames(imd_out$P_values)[1] <- colnames(imd_out$Full_results)[1]
  
  #Add attributes
  attr(imd_out, "group_DF") <- groupData
  #if(is.null(comparisons)){
  comparisons <- colnames(imd_out_flags)[grep("flag", tolower(names(imd_out_flags)))]
  #}
  attr(imd_out, "comparisons") <- comparisons
  attr(imd_out, "number_significant") <-  data.frame(Comparison=comparisons, Up = apply(flags,2,function(x){length(which(x>0))}),
                                                     Down = apply(flags,2,function(x){length(which(x<0))}))
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
#' @param plot_type defines which plots to be produced, options are "bar", "heatmap", "volcano"; if left `NULL` then all options are produced
#' @param fc_threshold optional threshold value for fold change estimates that's added to the volcano plot
#' 
#' @export
#' @method plot statRes
#' @examples 
#' 
#' library(MSomicsSTAT)
#' library(MSomicsDATA)
#' library(MSomicsQC)
#' #Trasnform the data
#' attr(pro_proData,"cnames")$fdata_cname <- "SampleID"
#' 
#' #Group the data by condition
#' myproData <- group_designation(omicsData = pro_proData, main_effects = c("Condition"))
#' 
#' #Apply the IMD ANOVA filter
#' imdanova_Filt <- imdanova_filter(omicsData = myproData)
#' myproData <- MSomics_filter(filter_object = imdanova_Filt, omicsData = myproData, min_nonmiss_anova=2)
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
plot.statRes <- function(x, plot_type = NULL, fc_threshold = NULL,...){
  
  #For now require ggplot2, consider adding base graphics option too
  if(!require(ggplot2)){
    stop("ggplot2 is required, base graphics aren't available yet")
  }
  
  #If no plot_type is specified, do all 4
  if(is.null(plot_type)){
    plot_type <- c("bar","heatmap","volcano")
  }
  
  #Most plots are based on "number_significant" data frame so pull it out
  comp_df <- attr(x,"number_significant")
  
  ##--------##
  #Go through the given plot_types and remove any that aren't currently avilable
  for(j in 1:length(plot_type)){
    plt_tyj <- try(match.arg(tolower(plot_type[j]),c("bar","volcano","heatmap")),silent=TRUE)
    if(class(plt_tyj)=='try-error'){
      warning(paste0("Plot type '",plot_type[j],"' is not currently available."))
      plot_type[j] <- NA
    }else{
      plot_type[j] <- plt_tyj
    }
  }
  
  #Don't make heatmaps if there's only one comparison
  if("heatmap"%in%plot_type & nrow(comp_df)==1){
    warning("Heatmaps are not produced when only one comparison is being made.")
    plot_type[which(plot_type=="heatmap")] <- NA
  }
  
  if(any(is.na(plot_type))){
    plot_type <- plot_type[-which(is.na(plot_type))]
  }
  ##--------##
  
  #I'll return a list of grobs
  all_plts <- vector('list',length(plot_type))
  
  #Bar plot
  if("bar"%in%plot_type){
    comp_df_melt <- reshape2::melt(comp_df,id.vars="Comparison",value.name="Count",variable.name="Direction")
    levels(comp_df_melt$Comparison) <- gsub(pattern = "_",replacement = " ",levels(comp_df_melt$Comparison))
    
    #Bar plots side-by-side, both going up
    #all_plts[[1]] <- ggplot(data=comp_df_melt,aes(Comparison,Count,fill=Direction))+geom_bar(stat='identity',position='dodge')
    
    #To turn the up direction to green, and down direction to red
    pal <- RColorBrewer::brewer.pal(9, "Set1")
    
    ##Up direction is positive, down direction is negative
    comp_df_melt[comp_df_melt$Direction=="Down",]$Count <- (-comp_df_melt[comp_df_melt$Direction=="Down",]$Count)
    all_plts[[1]] <- ggplot(data=comp_df_melt,aes(Comparison,Count,fill=Direction))+geom_bar(stat='identity')+
      geom_hline(aes(yintercept=0),colour='gray50')+
      theme(axis.text.x=element_text(angle=90,vjust=0.5))+
      scale_fill_manual(values=pal[c(3,1)])+
      ggtitle(paste("Number of DE Metabolites for", unique(attributes(x)$group_DF$VIRUS)[1], "vs", unique(attributes(x)$group_DF$VIRUS)[2], sep = " ")) +
      theme(plot.title = element_text(hjust = 0.5))
  }
  
  #Both the volcano plot and heatmap need a dataframe of fold changes by comparison/biomolecule
  if("volcano" %in%plot_type || "heatmap"%in%plot_type){
    fc_data <- x$Full_results[,c(1,grep("Fold_change",colnames(x$Full_results)))]
    colnames(fc_data) <- gsub(pattern = "Fold_change_",replacement = "",x = colnames(fc_data))
    fc_data <- reshape2::melt(fc_data,id.vars=1,variable.name="Comparison",value.name="Fold_change")
    
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
    
    volcano <- merge(fc_data,pvals,all=TRUE)
    levels(volcano$Comparison) <- gsub(pattern = "_",replacement = " ",levels(volcano$Comparison))
    
    #Remove NAs to avoid ggplot2 warning
    to_rm <- which(apply(volcano,1,function(x){any(is.na(x))}))
    if(length(to_rm)>0)
      volcano <- volcano[-to_rm,]
  }
  
  #Volcano plot 
  if("volcano"%in%plot_type){
    #Put plot into list at appropriate place
    elem_num <- min(which(unlist(lapply(all_plts,is.null))))
    all_plts[[elem_num]] <- ggplot(data=volcano,aes(Fold_change,-log(P_value,base=10)))+geom_point()+facet_grid(Type~Comparison)+
      ylab(expression(-log[10](p-value)))+xlab("Fold-change")
    if(!is.null(fc_threshold)){
      all_plts[[elem_num]] <- all_plts[[elem_num]]+geom_hline(aes(yintercept=fc_threshold))
    }
  }
  
  #Heatmap
  if("heatmap"%in%plot_type){
    # pal <- colorRampPalette(c("blue", "white", "red"))
    pal <- colorRampPalette(colors=rev(brewer.pal(11, "RdYlGn")))
    
    #For now just consider biomolecules significant with respect to ANOVA
    volcano <- dplyr::filter(volcano,Type=="ANOVA")
    
    volcano_sigs <- dplyr::filter(volcano,P_value<attr(x,"pval_thresh"))
    colnames(volcano_sigs)[1] <- "Biomolecule"
    volcano_sigs$Biomolecule <- as.factor(volcano_sigs$Biomolecule)
    #Put plot into list at appropriate place
    elem_num <- min(which(unlist(lapply(all_plts,is.null))))
    all_plts[[elem_num]] <- ggplot(volcano_sigs, aes(Biomolecule, Comparison)) +
      geom_tile(aes(fill = Fold_change), , color = "white") +
      scale_fill_gradient(low = "green", high = "red") +
      theme(axis.text.x = element_text(angle = 90, hjust = 0.5)) +
      ggtitle("Average Log Fold Change") +
      theme(plot.title = element_text(hjust = 0.5))
      # ggplot(data=volcano_sigs,aes(Comparison,Biomolecule)) +
      # geom_tile(aes(fill=Fold_change)) +
      # theme(axis.text.x=element_text(angle=90,vjust=0.5))+
      # scale_fill_brewer(palette = "RdYlGn")
  }
  
  gridExtra::grid.arrange(grobs=all_plts,nrow=ceiling(length(plot_type)/2))
  #invisible(all_plts)
}