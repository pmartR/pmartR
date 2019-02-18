#' Tests for a quantiative difference between groups (aka factors, aka main effects)
#'
#' This is the ANOVA part of the IMD-ANOVA test proposed in Webb-Robertson et al. (2010).
#'
#' The order in which different scenarios are handeled:
#' \enumerate{
#'  \item If the data are paired, then the pairing is accounted for first then each of the next steps is carried out on the new variable that is 
#'  the difference in the paired individuals.<br>
#'  \item If covariates are provided, their effect is removed before testing for group differences though mathematically covariates and grouping 
#'  effects are accounted for simultaneously
#'  \item ANOVA is executed to assess the effect of each main effects, results in a vector of group means for each biomolecule and variance estimate
#'  \item Group comparisons defined by `comaprison` argument are implemented use parameter vector and variance estimates in ANOVA step
#'}
#'
#' @param omicsData A pmartR data object of any class
#' @param comparisons data.frame with columns for "Control" and "Test" containing the different comparisons of interest. Comparisons will be made between the Test and the corresponding Control  If left NULL, then all pairwise comparisons are executed.
#' @param pval_adjust character vector specifying the type of multiple comparisons adjustment to implement. An unspecified value corresponds to no adjustment. Valid options include: holm, bonferroni, Tukey, Dunnett, none.
#' @param pval_thresh numeric p-value threshold, below or equal to which peptides are considered differentially expressed. Defaults to 0.05
#' @param covariates data.frame similar to \code{groupData} consisting of two columsn: the sample ID variable (with names matching column names in \code{omicsData$e_data}) and a column containing the numeric or group data for each sample
#' @param paired logical; should the data be paired or not? if TRUE then the `f_data` element of `omicsData` is checked for a "Pair" column, an error is returned if none is found
#' @param equal_var logical; should the variance across groups be assumed equal?
#' 
#' 
#' @return  a list of `data.frame`s
#' \tabular{ll}{
#' Results  \tab Edata cname, Variance Estimate, ANOVA F-Statistic, ANOVA p-value, Group means\cr
#'  \tab \cr
#' Fold_changes  \tab Estimated fold-changes for each comparison \cr
#'  \tab \cr
#' Fold_changes_pvalues  \tab P-values corresponding to the fold-changes for each comparison \cr
#'  \tab \cr
#' Fold_change_flags  \tab Indicator of statistical significance (0/+-2 to if adjusted p-value>=pval_thresh or p-value<pval_thresh) \cr
#'  }
#' @author Bryan Stanfill
#' @references 
#' Webb-Robertson, Bobbie-Jo M., et al. "Combined statistical analyses of peptide intensities and peptide occurrences improves identification of significant peptides from MS-based proteomics data." Journal of proteome research 9.11 (2010): 5748-5756.
#'
#' @examples 
#' dontrun{
#' library(pmartRdata)
#' library(pmartR)
#' mypepData <- edata_transform(omicsData = pep_object, data_scale = "log2")
#' mypepData <- group_designation(omicsData = mypepData, main_effects = c("Condition"))
#' 
#' #Try running anova_test without filtering, should get warning because the data haven't been filtered yet
#' anova_res <- anova_test(omicsData = mypepData) 
#' 
#' #Now filter and run again
#' imdanova_Filt <- imdanova_filter(omicsData = mypepData)
#' mypepData <- applyFilt(filter_object = imdanova_Filt, omicsData = mypepData, min_nonmiss_anova=2)
#' anova_res <- anova_test(omicsData = mypepData)
#' anova_res_tukey <- anova_test(omicsData = mypepData, pval_adjust = 'tukey') 
#' #Should be equivalent to above since only making one comparison
#' all(anova_res$Comparison==anova_res_tukey$Comparisons)
#' 
#' #If group IDs are treated as a covariate, all of the fold changes should be zero
#' covars <- attr(mypepData, "group_DF")
#' colnames(covars)[2] <- "Gp"
#' anova_res <- anova_test(omicsData = mypepData, covariates = covars)
#' summary(anova_res$Fold_changes)
#' 
#' #Test with really big dataset, one factor
#' library(OvarianPepdataBPsubset)
#' tcga_ovarian_pepdata_bp <- as.pepData(e_data = tcga_ovarian_pepdata_bp_subset$e_data, f_data = tcga_ovarian_pepdata_bp_subset$f_data, e_meta = tcga_ovarian_pepdata_bp_subset$e_meta, edata_cname = "Peptide", fdata_cname = "sampleID", emeta_cname = "Protein", check.names = FALSE)
#' tcga_ovarian_pepdata_bp <- group_designation(omicsData = tcga_ovarian_pepdata_bp, main_effects = c("vital_status"))
#' tcga_ovarian_pepdata_bp <- edata_transform(tcga_ovarian_pepdata_bp, "log2")
#' imdanova_Filt <- imdanova_filter(omicsData = tcga_ovarian_pepdata_bp)
#' tcga_ovarian_pepdata_bp <- applyFilt(filter_object = imdanova_Filt, omicsData = tcga_ovarian_pepdata_bp, min_nonmiss_anova=2)
#' ovarian_res <- anova_test(omicsData = tcga_ovarian_pepdata_bp) 
#' #Tukey adjustment is super slow right now because "ptukey" is super slow, not sure how to fix that
#' ovarian_res_tukey <- anova_test(omicsData = tcga_ovarian_pepdata_bp, pval_adjust = 'tukey') 
#' #Dunnett adjustment, super slow because mvtnorm::pmvt is super slow
#' ovarian_res_dunnett <- anova_test(omicsData = tcga_ovarian_pepdata_bp, pval_adjust = 'dunnett') 
#' 
#' #Test really big dataset, two factors, all pairwise comparisons
#' tcga_ovarian_pepdata_bp <- group_designation(omicsData = tcga_ovarian_pepdata_bp, main_effects = c("vital_status","neoplasm_histologic_grade"))
#' ovarian_res_twofac <- anova_test(omicsData = tcga_ovarian_pepdata_bp)
#' 
#' #Same but only test main effects (Dead vs Alive, G2 vs G3)
#' comp_df <- data.frame(Control=c("Alive","G2"), Test=c("Dead","G3"))
#' ovarian_res_twofac_main_effects <- anova_test(omicsData = tcga_ovarian_pepdata_bp, comparisons = comp_df)
#' 
#' #Same but only test arbitrary diagonal effects (Dead_G2 vs Alive_G3, Alive_G2 vs Alive_G3)
#' comp_df <- data.frame(Control=c("Dead_G2","Alive_G2"), Test=c("Alive_G3","Alive_G3"))
#' ovarian_res_twofac_arb_effects <- anova_test(omicsData = tcga_ovarian_pepdata_bp, comparisons = comp_df)
#' }
#' @export

anova_test <- function(omicsData, comparisons = NULL, pval_adjust = 'none', pval_thresh = 0.05, covariates = NULL, paired = FALSE, equal_var = TRUE){
  # check that omicsData is of the appropriate class
  if(!inherits(omicsData, c("proData","pepData","lipidData", "metabData"))) stop("omicsData is not an object of appropriate class")
  
  # Check for group_DF attribute #
  if(!("group_DF" %in% names(attributes(omicsData)))){
    stop("group_designation must be called in order to create a 'group_DF' attribute for omicsData.")
  }else{
    groupData <- attributes(omicsData)$group_DF
  }
  
  #Catch if number of groups is too small
  k <- length(unique(groupData$Group))
  if(k<2){
    stop("At least two groups are necessary to perform an ANOVA.")
  }
  
  # Check for log transform #
  if(!(attr(omicsData,"data_info")$data_scale%in%c("log2","log"))){
    stop("Data must be log transformed in order to implement ANOVA.")
  }
  
  ###--------Check for anova filter-------------###
  if(is.null(attr(omicsData,"imdanova"))){
    warning("These data haven't been filtered, see `?imdanova_filter` for details.")
  }else{
    cnames <- omicsData$e_data[,attr(omicsData,"cnames")$edata_cname]
    filterrows <- which(cnames%in%attr(omicsData,"imdanova")$test_with_anova)
    if(length(filterrows)>0)
      omicsData$e_data <- omicsData$e_data[filterrows,]
  }
  
  #----Check to see if appropriate p-value adjustment technique is chosen----#
  #all pairwise -> Tukey
  #Case versus control -> Dunnett
  if(length(unique(comparisons$Control))==1 & tolower(pval_adjust)=="tukey" & length(comparisons$Control)>1){
    stop("Tukey adjustment for multiple comparisions should only be used when all pairwise comparisons are being made; try a 'Dunnett' adjustment.")
  }
  
  if((length(unique(comparisons$Control))>1 | (is.null(comparisons) & length(unique(groupData$Group))>2)) & tolower(pval_adjust)=="dunnett"){
    stop("Dunnett adjustment for multiple comparisions should only be used when all case-vs-control comparisons are being made; try a 'Tukey' adjustment.")
  }

  #------Remove rows in "groupData" that don't have corresponding columns in 'omicsData$edata'-------###
  samp_cname <- attr(omicsData,"cnames")$samp_cname
  if(is.null(samp_cname)){ #Added becuase "samp_cname" has been replaced with "fdata_cname"
    samp_cname <- attr(omicsData,"cnames")$fdata_cname
  }
  group_samp_cname <- groupData[,samp_cname] #The sampleIDs in groupData (possibly more than is in the e_data)
  groupData <- groupData[group_samp_cname%in%colnames(omicsData$e_data),] #Only keep the rows of groupData that have columns in e_data

  #Create a data matrix that only includes the samples in "groupData" and put them in the order specified by
  #"groupData"
  data <- omicsData$e_data[,as.character(groupData[,samp_cname])]

  ###--------If paired==TRUE then use the pairing to create pair adjusted abundances-------------###
  if(paired){
    if(!require(dplyr)){
      stop("The dplyr package is required to perform a paired analysis, please install it.")
    }
    ##--- check for one and only one "pair" column in f_data ----##
    fdata_names <- tolower(colnames(omicsData$f_data))
    pair_col <-grep("pair",fdata_names) 
    if(length(pair_col)==0){
      stop("When `paired=TRUE`, a `Pair` column must be supplied within `f_data`")
    }else if(length(pair_col)>1){
      stop("Too many pair ID columns have been provided.  Please remove the redundant one.")
    }
    
    ##---- each pair must have two observations ----#
    not_two <- which(unname(table(omicsData$f_data[,pair_col]))!=2)
    if(length(not_two)>0){
      pair_id <- names(table(omicsData$f_data[,pair_col]))[not_two]
      stop(paste("Pair",pair_id,"does not have the two observations needed to form a pair."))
    }
    
    ##---- use the pairing info to form paired differences ------##
    #source('~/pmartR/R/paired_supp_funs.R')
    cols <- c(which(colnames(omicsData$f_data)==samp_cname),pair_col)
    pid_matrix <- build_pair_mat(pair_df = omicsData$f_data[,cols])
    #The new "data" are the paired differences, so overwrite data with paried differences
    #Rcpp::sourceCpp('src/fc_functions.cpp')
    #data <- fold_change_diff_na_okay(data = data.matrix(omicsData$e_data[,-1]),C = t(pid_matrix)) #This failed if columns didn't 
    data <- fold_change_diff_na_okay(data = data.matrix(omicsData$e_data[,as.character(omicsData$f_data[,samp_cname])]),C = t(pid_matrix))
    
    #Add columns names 
    if(is.numeric(omicsData$f_data[,pair_col])){
      colnames(data) <- paste("Pair",omicsData$f_data[,pair_col][apply(pid_matrix,1,function(x)return(any(x>0)))])
    }else{
      colnames(data) <- omicsData$f_data[,pair_col][apply(pid_matrix,1,function(x)return(any(x>0)))]
    }
    #pid_matrix could potentially be large, delete it
    rm(pid_matrix)
    
    #groupData needs to be updated so that pairs map to groups instead of observations
    groupData <- merge(groupData,omicsData$f_data[,cols],sort = FALSE)
    groupData <- groupData%>%group_by_at(vars(contains("pair")))%>%dplyr::summarize(NewGroup=first(Group))
    
    colnames(groupData) <- c(samp_cname,"Group")
    groupData <- data.frame(groupData)
    if(is.numeric(omicsData$f_data[,pair_col])){
      groupData[,which(colnames(groupData)==samp_cname)] <- paste("Pair",groupData[,which(colnames(groupData)==samp_cname)])
    }
    #k needs to be updated too
    k <- length(unique(groupData$Group))
    
    #If provided, covariate data needs to be updated too
    if(!is.null(covariates)){
      covariates <- merge(covariates,omicsData$f_data[,cols],sort = FALSE)
      orig_colnames <- colnames(covariates)
      colnames(covariates)[c(2,3)] <- c("Cov","PairID")
      covariates <- covariates%>%group_by(PairID)%>%dplyr::summarize(NewCov=first(Cov))
      colnames(covariates) <- c(samp_cname,orig_colnames[2])
      covariates <- data.frame(covariates)
      if(is.numeric(omicsData$f_data[,pair_col])){
        covariates[,which(colnames(covariates)==samp_cname)] <- paste("Pair",covariates[,which(colnames(covariates)==samp_cname)])
      }
    }
  }
  
  #-----If covariates are provided, remove their effect before computing group means/variances------#
  if(!is.null(covariates)){
    #If the mapping between sample IDs and covariates is missing, ignore covariates
    cov_samp_col <- which(colnames(covariates)==samp_cname)
    if(length(cov_samp_col)==0){
      warning(paste(samp_cname,"information is missing from provided covariates thus covariates will be ignored"))
      red_df <- matrix(rep(0,nrow(data)),ncol=1)
    }else if(any(is.na(covariates[,-cov_samp_col]))){
      warning("Missing values were detected in the provided covariates thus covariates will be ignored")
      red_df <- matrix(rep(0,nrow(data)),ncol=1)
    }else{
      #Add group ids to covariate ids
      covariates <- merge(groupData, covariates,sort=FALSE)
      cov_samp_col <- which(colnames(covariates)==samp_cname)
      #source('~/pmartR/R/covariate_supp_funs.R') #Run to debug
      xmatrix <- build_x_mat(covariates[,-cov_samp_col]) #Build the appropriate X matrix based on the covariates data frame
      xmatrix <- reduce_xmatrix(xmatrix,length(unique(groupData$Group)))
      
      
      #Remove effect of covariates by projecting into the null space generated by the X matrix
      allproj <- project_to_null_cpp(data.matrix(data),xmatrix,length(unique(groupData$Group)))
      data <- allproj$data_no_x
      red_df <- allproj$lost_df
      #print(red_df)
    }
    #If covariates are adjusted for, then force equal variance assumption even if there are only two levels to main effect
    equal_var <- TRUE
  }else{
    red_df <- matrix(rep(0,nrow(data)),ncol=1)
  }
  
  ###-------- Use one or 2 factor ANOVA to compute group means and common variance--------###
  
  if(ncol(groupData)==2){
    ##---- One factor ANOVA ----##
    
    ##Translate the groups into numeric labels for anova_cpp() function
    gp <- factor(groupData$Group,labels=1:k,levels=unique(groupData$Group))
    
    #Rcpp::sourceCpp('~/pmartR/src/anova_helper_funs.cpp') #Run if debugging code
    raw_results <- anova_cpp(data.matrix(data),gp,1-equal_var,red_df)
    group_names <- paste("Mean",as.character(unique(groupData$Group)),sep="_")
    
  }else{
    ##---- Two factor ANOVA ----##
    #source('~/pmartR/R/anova_helper_fun.R') #Run if debugging
    #source('~/pmartR/R/covariate_supp_funs.R')
    #Rcpp::sourceCpp('~/pmartR/src/two_factor_anova.cpp')
    raw_results <- run_twofactor_cpp(data=data.matrix(data),gpData=groupData,red_df)
    group_names <- paste("Mean",colnames(raw_results$group_means),sep="_")
  }
  
  #The C++ function returns a list so we need to make it into a data.frame
  results <- data.frame(raw_results$Sigma2, raw_results$group_means, raw_results$Fstats, raw_results$pvalue)
  
  #Rename the columns to match group names
  colnames(results) <- c("Variance",group_names,"F_Statistic","p_value")
  
  #Pull off edata_cname and add to results df
  edatacname <- attr(omicsData,"cnames")$edata_cname
  results <- cbind(omicsData$e_data[edatacname],results)
  
  ###-------Use group_comparison_anova() to compare the groups that were requested--------------##
  #source('~/pmartR/R/group_comparison.R') # Run if debugging
  #Rcpp::sourceCpp('~/pmartR/src/group_comparisons.cpp') #Run if debugging
  group_comp <- group_comparison_anova(groupData=groupData,comparisons=comparisons,anova_results_full=list(Results=results,Sizes=raw_results$group_sizes))
  
  #If there are only two levels, replace group_comp p-values with those taken from ANOVA function
  if(k==2){
    group_comp$p_values <- raw_results$pvalue
  }
  
  #Different way to compute fold change: diff if log scale, ratio otherwise
  if(attr(omicsData,"data_info")$data_scale%in%c("log2","log")){
    fold_change <- group_comp$diff_mat
  }else{
    #CONSIDER MAKING GROUP_COMPARISON RETURN THE CMAT ALONG WITH EVERYTHING ELSE INSTEAD OF RECREATING IT HERE
    Cmat_res <- create_c_matrix(groupData, comparisons)
    Cmat <- Cmat_res$cmat
    fold_change <- fold_change_ratio(raw_results$group_means, Cmat)
    fold_change <- data.frame(fold_change)
    colnames(fold_change) <- colnames(group_comp$diff_mat)
  }  
  
  
  #--------Adjust fold-change p-values if applicable-----------##
  p_values_adjust <- p_adjustment_anova(p_values = group_comp$p_values, diff_mean = group_comp$diff_mat, t_stats = group_comp$t_tests,
                                  sizes = raw_results$group_sizes, pval_adjust = pval_adjust)
  colnames(p_values_adjust) <- colnames(fold_change)

  
  #-------Use the adjusted_pvals object to create significance flags-------------##
  #Create object to be returned
  sig_indicators <- data.matrix(p_values_adjust)
  
  sigs <- which(sig_indicators<pval_thresh)
  if(length(sigs)>0){
    sig_indicators[sigs] <- sign(data.matrix(group_comp$diff_mat)[sigs])
    sig_indicators[-sigs] <- 0
  }else{
    sig_indicators[which(sig_indicators>=pval_thresh)] <- 0
  }
  sig_indicators <- data.frame(sig_indicators)
  colnames(sig_indicators) <- colnames(fold_change)
  
  #------Return the full anova results-----------##
  #(estimated variance, F-statistic, p-value and group means), estimated fold changes, and adjusted p-values for the requested comparisons
  return(list(Results = results, Fold_changes = fold_change, Fold_change_pvalues = p_values_adjust, Flags = sig_indicators))
}

#########
#I HAVE QUESTIONS ABOUT FOLD CHANGES.  BELOW IS MY ORGIGINAL CODE, WHICH HAS BEEN COMMENTED OUT.
#THE NEW CODE SIMPLY RETURNS DIFFERENCE IF LOG (BASE E OR 2) SCALE AND RATIO OTHERWISE.  WHICH IS CORRECT?
#############
#Use appropriate function to compute fold change of data based on scaling
# if(attr(omicsData,"data_info")$data_scale=="log2"){
#   fold_change <- fold_change_logbase2(means, Cmat)
# }else if(attr(omicsData,"data_info")$data_scale=="log"){
#   #The difference between groups is already computed so it just needs to be exponentiated
#   fold_change <- exp(group_comp$diff_mat)
# }else{
#   fold_change <- fold_change_ratio(means, Cmat)
# }  


##Old p-value adjust method
#if(is.null(pval_adjust)){
#  p_values_adjust <- group_comp$p_values
#}else{
#  p_values_adjust <- t(apply(group_comp$p_values,1,p.adjust,method=pval_adjust))
#  p_values_adjust <- data.frame(p_values_adjust)
#  colnames(p_values_adjust) <- paste(colnames(fold_change),"p_value")
#}
