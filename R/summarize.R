#' Creates a list of six Data Frames, one for each summarizing metric
#'
#' This function takes in an omicsData object and returns a summary of the e_data component. The six summarizing metrics include, mean, standard deviation, median, percent observed, minimum and maximum.   
#'
#' @param omicsData an object of the class 'lipidData', 'metabData', 'pepData', or 'proData', usually created by \code{\link{as.lipidData}}, \code{\link{as.metabData}}, \code{\link{as.pepData}}, or \code{\link{as.proData}}, respectively.
#' @param groupvar a character vector with no more than two variable names that should be used to determine group membership of samples. The variable name must match a column name from \code{f_data}. Defaults to NULL, in which case group_DF attribute will be used.
#' @param by a character string indicating whether summarizing metrics will be applied by 'sample' or by 'molecule
#' 
#' @details If groupvar is NULL and group_designation has not been applied to omicsData, then the metrics will be applied to each column of e_data (when by = 'sample) or to each row of e_data (when by = 'molecule'). When groupvar is provided, it must match a column name from \code{f_data}, this column of f_data is used to group e_data in order to apply the metrics.  
#'
#' @return A list of six data frames, which are the results of applying the metrics to omicsData$e_data
#'
#' @examples
#' dontrun{
#' library(pmartRdata)
#' data(lipid_object)
#' lipid_object2 <- edata_transform(omicsData = lipid_object, data_scale = "log2")
#' lipid_object2 <- group_designation(omicsData = lipid_object, main_effects = "Condition")
#' summarize(omicsData = lipid_object2, groupvar = NULL, by = "sample")
#'}
#'
#'
#' @export

summarize<- function(omicsData, groupvar = NULL, by){
  #some checks
  if(!(class(omicsData) %in% c('pepData','proData', 'lipidData', 'metabData'))) stop("omicsData must be an object of class pepData, proData, lipidData or metabData")
  if(!(by %in% c('sample', 'molecule'))) stop("by must be either sample or molecule")
  
  #pull cnames attr from omicsData
  edata = omicsData$e_data
  fdata = omicsData$f_data
  edata_cname = attr(omicsData, "cnames")$edata_cname
  fdata_cname = attr(omicsData, "cnames")$fdata_cname
  edata_cname_id = which(names(edata) == edata_cname)
  fdata_cname_id = which(names(fdata) == fdata_cname)
  groupDF = attr(omicsData, "groupDF")
  
  if(by == 'sample'){
    #check that groupvar is NULL, groupvar is only used when by == 'molecule'
    if(!is.null(groupvar)) stop("groupvar is only used when by == 'molecule'")
  
    avg = as.data.frame(apply(edata[,-edata_cname_id], 2, function(x){if(all(is.na(x))){mean(x)}else{mean(x, na.rm = T)}}))
    avg = cbind(names(edata[,-edata_cname_id]), avg)
    names(avg)<- c("sample", "mean")
    rownames(avg)<- NULL
    
    sd = as.data.frame(apply(edata[,-edata_cname_id], 2, sd, na.rm = T))
    sd = cbind(names(edata[,-edata_cname_id]), sd)
    names(sd)<- c("sample", "sd")
    rownames(sd)<- NULL
    
    mds = as.data.frame(apply(edata[,-edata_cname_id], 2, median, na.rm = T))
    mds = cbind(names(edata[,-edata_cname_id]), mds)
    names(mds)<- c("sample", "median")
    rownames(mds)<- NULL
    
    pct_obs = as.data.frame(apply(edata[,-edata_cname_id], 2, function(x){sum(!is.na(x))/length(x)}))
    pct_obs = cbind(names(edata[,-edata_cname_id]), pct_obs)
    names(pct_obs)<- c("sample", "pct_obs")
    rownames(pct_obs)<- NULL
    
    min = as.data.frame(apply(edata[,-edata_cname_id], 2, min, na.rm = T))
    min = cbind(names(edata[,-edata_cname_id]), min)
    names(min)<- c("sample", "min")
    rownames(min)<- NULL
  
    max = as.data.frame(apply(edata[,-edata_cname_id], 2, max, na.rm = T))
    max = cbind(names(edata[,-edata_cname_id]), max)
    names(max)<- c("sample", "max")
    rownames(max)<- NULL
        
    res_list = list(mean = avg, sd = sd, median = mds, pct_obs = pct_obs, min = min , max = max)
    class(res_list)<- "dataRes"
    attr(res_list, "by")<- by
    attr(res_list, "groupvar")<- groupvar
    attr(res_list, "cnames")<- list("edata_cname" = edata_cname, "fdata_cname" = fdata_cname)
  }
  
  if(by == "molecule"){
    if(is.null(groupvar)){
      if(is.null(attr(omicsData, "group_DF"))){
        #when groupvar is NULL and group_designation is NULL, we calculate metric for each row
        avg = apply(edata[, -edata_cname_id], 1, function(x){if(all(is.na(x))){mean(x)}else{mean(x, na.rm = T)}})
        avg = data.frame(molecule = edata[, edata_cname_id], mean = avg, stringsAsFactors = F)
        names(avg)[1]<- edata_cname
        sd = apply(edata[, -edata_cname_id], 1, sd, na.rm = T)
        sd = data.frame(molecule = edata[, edata_cname_id], sd = sd, stringsAsFactors = F)
        names(sd)[1]<- edata_cname
        mds = apply(edata[, -edata_cname_id], 1, median, na.rm = T)
        mds = data.frame(molecule = edata[, edata_cname_id], median = mds, stringsAsFactors = F)
        names(mds)[1]<- edata_cname
        pct_obs = apply(edata[,-edata_cname_id], 1, function(x){sum(!is.na(x))/length(x)})
        pct_obs = data.frame(molecule = edata[, edata_cname_id], pct_obs = pct_obs, stringsAsFactors = F)
        names(pct_obs)[1]<- edata_cname
        min = apply(edata[, -edata_cname_id], 1, min, na.rm = T)
        min = data.frame(molecule = edata[, edata_cname_id], min = min, stringsAsFactors = F)
        names(min)[1]<- edata_cname
        max = apply(edata[, -edata_cname_id], 1, max, na.rm = T)
        max = data.frame(molecule = edata[, edata_cname_id], max = max, stringsAsFactors = F)
        names(max)[1]<- edata_cname
        
        res_list = list(mean = avg, sd = sd, median = mds, pct_obs = pct_obs, min = min , max = max)
        class(res_list)<- "dataRes"
        attr(res_list, "by")<- by
        attr(res_list, "groupvar")<- groupvar
        attr(res_list, "cnames")<- list("edata_cname" = edata_cname, "fdata_cname" = fdata_cname)
        attr(res_list, "groupDF")<- groupDF
        
      }else if(!is.null(attr(omicsData, "group_DF"))){
        #when groupvar is NULL but group_designation has been run
        edata_melt = reshape2::melt(edata, id.vars = edata_cname)
        names(edata_melt)[2]<- fdata_cname
        edata_melt = merge.data.frame(edata_melt, attr(omicsData, "group_DF"), by = fdata_cname)
        edata_melt = edata_melt[, -which(names(edata_melt) == fdata_cname)]
        
        #here we are creating a string to input for dcast function argument 'formula'
        formula1 = paste(edata_cname, "+Group~...", sep = "")
        formula2 = paste(edata_cname, "~...", sep = "")

        avg = dcast(edata_melt, formula = formula1, function(x){if(all(is.na(x))){mean(x)}else{mean(x, na.rm = T)}})
        names(avg)[which(colnames(avg)== ".")]<- "value"
        avg = dcast(avg, formula = formula2)
        std_div = dcast(edata_melt, formula = formula1, sd, na.rm = T)
        names(std_div)[which(colnames(std_div)== ".")]<- "value"
        std_div = dcast(std_div, formula = formula2)
        mds = dcast(edata_melt, formula = formula1, median, na.rm = T)
        names(mds)[which(colnames(mds)== ".")]<- "value"
        mds = dcast(mds, formula = formula2)
        pct_obs = dcast(edata_melt, formula = formula1, function(x){sum(!is.na(x))/length(x)})
        names(pct_obs)[which(colnames(pct_obs)== ".")]<- "value"
        pct_obs = dcast(pct_obs, formula = formula2)
        mins = dcast(edata_melt, formula = formula1, function(x){if(all(is.na(x))){min(x)}else{min(x, na.rm = T)}})
        names(mins)[which(colnames(mins)== ".")]<- "value"
        mins = dcast(mins, formula = formula2)
        maxs = dcast(edata_melt, formula = formula1, function(x){if(all(is.na(x))){max(x)}else{max(x, na.rm = T)}})
        names(maxs)[which(colnames(maxs)== ".")]<- "value"
        maxs = dcast(maxs, formula = formula2)
        
        res_list = list(mean = avg, sd = std_div, median = mds, pct_obs = pct_obs, min = mins , max = maxs)
        class(res_list)<- "dataRes"
        attr(res_list, "by")<- by
        attr(res_list, "groupvar")<- groupvar
        attr(res_list, "cnames")<- list("edata_cname" = edata_cname, "fdata_cname" = fdata_cname)
        attr(res_list, "groupDF")<- groupDF
      }
    }else if(length(groupvar) == 1){
      ####case where groupvar is provided and has length 1####
      if(any(groupvar %in% names(fdata)) == F) stop("groupvar must be in omicsData f_data")
      temp_fdata = fdata[, c(fdata_cname_id, which(names(fdata) == groupvar))]
      names(temp_fdata)[2]<- "Group"
      
      #rearranging edata
      edata_melt = reshape2::melt(omicsData$e_data, id.vars = edata_cname)
      names(edata_melt)[2]<- fdata_cname
      edata_melt = merge.data.frame(edata_melt, temp_fdata, by = fdata_cname)
      edata_melt = edata_melt[, -which(names(edata_melt) == fdata_cname)]
      
      #here we are creating a string to input for dcast function argument 'formula'
      formula1 = paste(edata_cname, "+Group~...", sep = "")
      formula2 = paste(edata_cname, "~...", sep = "")
      
      avg = dcast(edata_melt, formula = formula1, function(x){if(all(is.na(x))){mean(x)}else{mean(x, na.rm = T)}})
      names(avg)[which(colnames(avg)== ".")]<- "value"
      avg = dcast(avg, formula = formula2)
      std_div = dcast(edata_melt, formula = formula1, sd, na.rm = T)
      names(std_div)[which(colnames(std_div)== ".")]<- "value"
      std_div = dcast(std_div, formula = formula2)
      mds = dcast(edata_melt, formula = formula1, median, na.rm = T)
      names(mds)[which(colnames(mds)== ".")]<- "value"
      mds = dcast(mds, formula = formula2)
      pct_obs = dcast(edata_melt, formula = formula1, function(x){sum(!is.na(x))/length(x)})
      names(pct_obs)[which(colnames(pct_obs)== ".")]<- "value"
      pct_obs = dcast(pct_obs, formula = formula2)
      mins = dcast(edata_melt, formula = formula1, function(x){if(all(is.na(x))){min(x)}else{min(x, na.rm = T)}})
      names(mins)[which(colnames(mins)== ".")]<- "value"
      mins = dcast(mins, formula = formula2)
      maxs = dcast(edata_melt, formula = formula1, function(x){if(all(is.na(x))){max(x)}else{max(x, na.rm = T)}})
      names(maxs)[which(colnames(maxs)== ".")]<- "value"
      maxs = dcast(maxs, formula = formula2)
      
      res_list = list(mean = avg, sd = std_div, median = mds, pct_obs = pct_obs, min = mins , max = maxs)
      class(res_list)<- "dataRes"
      attr(res_list, "by")<- by
      attr(res_list, "groupvar")<- groupvar
      attr(res_list, "cnames")<- list("edata_cname" = edata_cname, "fdata_cname" = fdata_cname)
      
    }else if(length(groupvar) == 2){
      ####case where length of groupvar is 2####
      if (length(groupvar) > 2) stop("No more than two groupvar can be provided")
      if(any(groupvar %in% names(fdata)) == F) stop("groupvar must be in omicsData f_data")
      
      group_vars = fdata[, names(fdata) %in% groupvar]
      group_vars <- apply(group_vars, 2, as.character)
      
      # create a group variable and paste grouvar levels together for samples #
      # samples with a value of NA for either groupvar will have a Group value of NA #
      Group = rep(NA, nrow(fdata))
      
      # identify samples that will have a Group membership that is not missing #
      nonna.group = (!is.na(group_vars[,1]) & !is.na(group_vars[,2]))
      Group[nonna.group] = paste(as.character(group_vars[nonna.group,1]), as.character(group_vars[nonna.group,2]), sep = "_")
      
      # create output formatted with first column being fdata_cname and second column group id #
      output = data.frame(Sample.ID = fdata[,fdata_cname], Group = Group)
      names(output)[1] = fdata_cname
      
      #check that there are atleast 2 sample in each group and remove groups that have less than two samples per group
      n_per_grp = as.data.frame(output %>% group_by(Group) %>% summarise(count = n()))
      remove_group = as.character(n_per_grp[which(n_per_grp$count<2), "Group"])
      n_per_grp = n_per_grp[-which(n_per_grp$count<2),]
      
      output = output[-which(output$Group == remove_group),]
      
      #groupvar was provided, rearranging edata
      edata_melt = reshape2::melt(omicsData$e_data, id.vars = edata_cname)
      names(edata_melt)[2]<- fdata_cname
      edata_melt = merge.data.frame(edata_melt, output, by = fdata_cname)
      edata_melt = edata_melt[, -which(names(edata_melt) == fdata_cname)]
      
      #here we are creating a string to input for dcast function argument 'formula'
      formula1 = paste(edata_cname, "+Group~...", sep = "")
      formula2 = paste(edata_cname, "~...", sep = "")
      
      avg = dcast(edata_melt, formula = formula1, function(x){if(all(is.na(x))){mean(x)}else{mean(x, na.rm = T)}})
      names(avg)[which(colnames(avg)== ".")]<- "value"
      avg = dcast(avg, formula = formula2)
      std_div = dcast(edata_melt, formula = formula1, sd, na.rm = T)
      names(std_div)[which(colnames(std_div)== ".")]<- "value"
      std_div = dcast(std_div, formula = formula2)
      mds = dcast(edata_melt, formula = formula1, median, na.rm = T)
      names(mds)[which(colnames(mds)== ".")]<- "value"
      mds = dcast(mds, formula = formula2)
      pct_obs = dcast(edata_melt, formula = formula1, function(x){sum(!is.na(x))/length(x)})
      names(pct_obs)[which(colnames(pct_obs)== ".")]<- "value"
      pct_obs = dcast(pct_obs, formula = formula2)
      mins = dcast(edata_melt, formula = formula1, function(x){if(all(is.na(x))){min(x)}else{min(x, na.rm = T)}})
      names(mins)[which(colnames(mins)== ".")]<- "value"
      mins = dcast(mins, formula = formula2)
      maxs = dcast(edata_melt, formula = formula1, function(x){if(all(is.na(x))){max(x)}else{max(x, na.rm = T)}})
      names(maxs)[which(colnames(maxs)== ".")]<- "value"
      maxs = dcast(maxs, formula = formula2)
      
      res_list = list(n_per_grp = n_per_grp, mean = avg, sd = std_div, median = mds, pct_obs = pct_obs, min = mins , max = maxs) 
      class(res_list)<- "dataRes"
      attr(res_list, "by")<- by
      attr(res_list, "groupvar")<- groupvar
      attr(res_list, "cnames")<- list("edata_cname" = edata_cname, "fdata_cname" = fdata_cname)
    }
  }
  return(res_list)
}

