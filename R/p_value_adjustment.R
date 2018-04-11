#' Adjust p-values for multiple comparisons
#'
#' Depending upon the \code{pval_adjust} method selected, the supplied p_values are compared against an adjusted \code{pval_thresh} value or the provided
#' means are used to compute new statistics, p-values are computed and compared against the provided \code{pval_thresh}.  A \code{data.frame} that indicates which
#' of the tests are significant, 1 if significant or 0 if insignificant.  If \code{means} is also provided and the p-value is signficant then the direction
#' of the change is indicated by the sign on 1, i.e., means<0 and p_value<pval_thresh will return -1, similarly for means>0.
#'
#' @param p_values A matrix (or \code{data.frame}) of p-values to be adjusted.
#' @param diff_mean A matrix (or \code{data.frame}) of groups means that are to be compared
#' @param t_stats A matrix (or \code{data.frame}) of t-test statistics resulting from from standard procedures
#' @param sizes A matrix (or \code{data.frame}) of group sizes
#' @param pval_adjust character vector specifying the type of multiple comparisons adjustment to implement. A NULL value corresponds to no adjustment. Valid options include: holm, bonferroni, dunnett, tukey or none. 
#'
#' @return a data frame with the following columns: group means, global G-test statistic and corresponding p-value
#'
#' @author Bryan Stanfill
#' @examples 
#' library(MSomicsDATA)
#' library(MSomicsQC)
#' mypepData <- edata_transform(mintR_data = pep_pepData, data_scale = "log2")
#' mypepData <- group_designation(mintR_data = mypepData, main_effects = c("Condition"))
#' group_df <- attr(mypepData, "group_DF")
#' imdanova_Filt <- imdanova_filter(mintR_data = mypepData)
#' mypepData <- MSomics_filter(filter_object = imdanova_Filt, mintR_data = mypepData, min_nonmiss_anova=2)
#' anova_res <- anova_test(omicsData = mypepData, groupData = group_df)
#' 
#' adjusted_pvalues <- p_adjustment_anova(p_values = anova_res$Comparisons, diff_mean = anova_res$Fold)
#' @export

p_adjustment_anova <- function(p_values = NULL, diff_mean = NULL, t_stats = NULL, sizes = NULL, pval_adjust = "None"){
  
  #Match provided 'pval_adjust' to available options
  pval_adjust <- try(match.arg(tolower(pval_adjust),c("bonferroni","none","tukey","dunnett","holm")),silent=TRUE)
  
  if(class(pval_adjust)=='try-error')
    stop("Provided 'pval_adjust' argument is invalid, please select 'holm', 'bonferroni', 'tukey', 'dunnett' or 'none'.")
  
  
  if(pval_adjust=="tukey" || pval_adjust=="dunnett"){
    
    #Tukey-Kramer statistics are t-statistic/sqrt(2)
    if(is.null(t_stats) || is.null(sizes)){
      stop("The standard t-tests and group sizes need to be supplied in order to apply Tukey-Kramer adjustment.")
    }
    
    n_compares <- ncol(t_stats)
    
    if(n_compares==1){
      #Tukey adjustment is not done if only one comparison is provided
      pval_adjust <- 'none'
    }else{
      
      if(pval_adjust=="tukey"){
        tukey_stats <- data.matrix(t_stats*sqrt(2))
        #Tukey only needs sizes for the rows, not individual group sizes
        if(is.data.frame(sizes)){
          sizes <- rowSums(data.matrix(sizes))
        }else if(is.matrix(sizes)){
          sizes <- rowSums(sizes)
        }
        
        #Rcpp::sourceCpp('src/tukey_helper.cpp') #Use for debugging
        adjusted_pvals <- ptukey_speed(tukey_stats,sizes)

      }else{#Dunnett adjustment - Needs to be sped up
        adjusted_pvals <- matrix(NA,nrow(t_stats), ncol(t_stats))
        #colnames(adjusted_pvals) <- colnames(t_stats)
        
        for(i in 1:nrow(adjusted_pvals)){
          k <- length(which(!is.na(p_values[i,])))
          if(k>0){
            dfi <- sum(sizes[i,])-k
            rm_nas <- which(is.na(p_values[i,]))
            #Modified version of an example from "Additional multcomp Examples" viggnette of the multcomp package
            #Until we add covariates the correlation matrix is assumed to be the identity
            if(length(rm_nas)>0){
              adjusted_pvals[i,-rm_nas] <- as.numeric(sapply(abs(t_stats[i,-rm_nas]), function(x,k,df) 1 - mvtnorm::pmvt(-rep(x, k), rep(x, k), df = df),k=k,df=dfi))
            }else{
              adjusted_pvals[i,] <- as.numeric(sapply(abs(t_stats[i,]), function(x,k,df) 1 - mvtnorm::pmvt(-rep(x, k), rep(x, k), df = df),k=k,df=dfi))
            }
          }
        }
      }
      
      colnames(adjusted_pvals) <- colnames(t_stats)
      colnames(adjusted_pvals) <- gsub("Test-Stat","p-value",colnames(adjusted_pvals)) #This is a band-aid right now, may need to make more general later
    }
  }
  
  #Don't make this an "else if" because pval_adjust can be changed to 'none' if n_compares==1
  if(pval_adjust%in%c('none',"bonferroni","holm")){
    
    #p_values needs to be supplied to apply bonferroni adjustment, if it's NULL tell them that's a problem
    if(is.null(p_values)){
      stop("The `p_values` argument must be supplied to perform the selected `pval_adjust` method")
    }
    
    if(is.null(dim(p_values))||ncol(p_values)==1){
      pval_adjust='none'
    }
    
    if(pval_adjust !='holm'){
      #For bonferroni adjustment, multiply p_values by number of comparisons (columns in p_values) else do no adjustment
      multiplier <- ifelse(pval_adjust=="bonferroni",ncol(p_values),1)
      adjusted_pvals <- multiplier*p_values
    }else{
      #Rcpp::sourceCpp('src/holm_adjustment.cpp') #For debugging
      #source('~/Documents/MinT/MSomics/MSomicsSTAT/R/support_p_adjustment.R') #For debugging
      adjusted_pvals <- t(apply(p_values,1,ranked_holm_cpp))

      #NaN p-values should stay NaN
      nan_pvals <- lapply(p_values,is.nan)
      for(i in 1:ncol(adjusted_pvals)){
        adjusted_pvals[,i][nan_pvals[[i]]] <- NaN
      }
    }
    
  }
  
  ##############
  #Return the adjusted p-values
  
  return(pvalues=data.frame(adjusted_pvals))
}

#Can be used to replace "adjusted_pvals <- ptukey_speed(tukey_stats,sizes)" if necessary
#adjusted_pvals <- matrix(NA,nrow(p_values), ncol(p_values))
#for(i in 1:nrow(tukey_stats)){
#  adjusted_pvals[i,] <- ptukey(abs(tukey_stats[i,]), nmeans = n_compares, df = sizes[i]-n_compares,lower.tail=FALSE)
#}
