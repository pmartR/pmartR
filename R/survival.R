
# MSomicsQC::surv_designation is a fucntion that creates a "surv_DF" attribute with elements:
#   1) t_death - defines column name in f_data
#   2) t_progress
#   3) ind_death
#   4) ind_progress
#   5) covariates
#   
#   Must appear together:
#   
#   t_death & ind_death (overall survival - OS option)
# OR 
#   t_death, t_progress, ind_death AND ind_progress (progression free survival - PFS option)
#   
#   covariates is always optional - two max for now
#   
#   ---------
#   Once Lisa/Kelly makes the above, write function to create Kaplan-Meier Curves.  
#   Only account for categorical variables when making these curves (separate curves for each level of each
#                                                                    covariate)
#   
#   Functions skeleton:
#     plot_km <- function(omicsData,...Other figure options from pmartRqc plot functions...) - produce KM plot
#     summary_km <- function(omicsData, percent = ?) - return something like summary.survfit where "percent" gives a
#       survival prob and they want "time" variable returned


#' Basic survival analysis function
#' 
#' Implements overall survival analysis or progression-free survival analysis, depending upon the datatypes supplied to
#' surv_designation, and return the "survfit" object
#' 
#' @param omicsData A pmartR data object of any class, which has a `group_df` attribute that is usually created by the `group_designation()` function
#' @return if fitted survival analysis object is returned
#' 
#' @examples 
#' \dontrun{
#' library(MSomicsSTAT)
#' library(OvarianPepdataBP)
#' 
#' #Basic analysis without covariates
#' attr(tcga_ovarian_pepdata_bp,"survDF") <- list(t_death = "survival_time",ind_death = "vital_status")
#' sfit <- fit_surv(tcga_ovarian_pepdata_bp)
#' plot(sfit)
#' 
#' #Add some covariate information
#' attr(tcga_ovarian_pepdata_bp,"survDF") <- list(t_death = "survival_time",ind_death = "vital_status", covariates = "g__initial_pathologic_diagnosis_method_g1")
#' sfit <- fit_surv(tcga_ovarian_pepdata_bp)
#' plot(sfit,col=c(1,2))
#' }
#' 
#Function to fit the survival model
fit_surv <- function(omicsData){
  
  if(!requireNamespace("survival", quietly = TRUE)){
    stop("Please install the 'survival' package.")  
  }
  
  t_death <- attr(omicsData,"survDF")$t_death
  t_progress <- attr(omicsData,"survDF")$t_progress
  ind_death <- attr(omicsData,"survDF")$ind_death
  ind_progress <- attr(omicsData,"survDF")$ind_progress
  covariates <- attr(omicsData,"survDF")$covariates
  
  surv_df <- data.frame(time=omicsData$f_data[,t_death], status = 0)
  surv_df$status[grep("dead",omicsData$f_data[,ind_death],ignore.case = TRUE)] <- 1
  
  
  if(is.null(covariates)){
    sfit <- survfit(Surv(time, status)~1,data=surv_df)
  }else{
    
    form <- "Surv(time,status)~"
    ncovars <- length(covariates)
    
    for(i in 1:ncovars){
      surv_df <- cbind(surv_df,omicsData$f_data[,covariates[i]])
      form <- paste0(form,covariates[i])
      if(i<ncovars){
        form <- paste0(form,"+")
      }
    }
    colnames(surv_df)[ncol(surv_df)-c((ncovars-1):0)] <- covariates
    
    sfit <- survfit(as.formula(form),data=surv_df)
  }
  return(sfit)
}


#' Basic survival analysis plot
#' 
#' Implements overall survival analysis or progression-free survival analysis, depending upon the datatypes supplied to
#' surv_designation, and plot the resulting Kaplan-Meier curve.
#' 
#' @param omicsData A pmartR data object of any class, which has a `group_df` attribute that is usually created by the `group_designation()` function
#' 
#' @return a Kaplan-Meier curve
#' 
#' @examples 
#' \dontrun{
#' library(MSomicsSTAT)
#' library(OvarianPepdataBP)
#' attr(tcga_ovarian_pepdata_bp,"survDF") <- list(t_death = "survival_time",ind_death = "vital_status")
#' plot_km(omicsData = tcga_ovarian_pepdata_bp)
#' 
#' #Add covariates to "survDF" attribute
#' attr(tcga_ovarian_pepdata_bp,"survDF") <- list(t_death = "survival_time",ind_death = "vital_status", covariates = "age_at_initial_pathologic_diagnosis")
#' plot_km(omicsData = tcga_ovarian_pepdata_bp)
#' }
plot_km <- function(omicsData){
  

  sfit <- fit_surv(omicsData)
  
  
  plot(sfit)
}

#' Basic survival analysis summary
#' 
#' Implements overall survival analysis or progression-free survival analysis, depending upon the datatypes supplied to
#' surv_designation, and gives a summary of the results.
#' 
#' @param omicsData A pmartR data object of any class, which has a `group_df` attribute that is usually created by the `group_designation()` function
#' @param percent The percentile
#' @return if `percent` is provided then the time at which that probability of death is returned; else, the summary of the `survival` object is returned
#' 
#' @examples 
#' \dontrun{
#' library(OvarianPepdataBP)
#' attr(tcga_ovarian_pepdata_bp,"survDF") <- list(t_death = "survival_time",ind_death = "vital_status")
#' #No percent is provided so the entire object is returned
#' summary_km(tcga_ovarian_pepdata_bp)
#' 
#' #Percent is provided so corresponding time point is returned
#' summary_km(tcga_ovarian_pepdata_bp, .4)
#' }
#' 
summary_km <- function(omicsData, percent = NULL,...){
  
  
  
  #Summary object
  sfit1 <- summary(fit_surv(omicsData,...))
  
  if(!is.null(percent)){
    row_of_interest <- which.min(abs(sfit1$surv-percent))
    return(sfit1$time[row_of_interest])
  }
  
  return(sfit1)
}