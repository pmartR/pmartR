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
  
  if(by == 'sample'){
    if(is.null(groupvar)){
      if(is.null(attr(omicsData, "group_DF"))){
        #when groupvar is NULL and group_designation is NULL, we calculate metric for each column 
        avg = as.data.frame(edata[, -edata_cname_id] %>% summarise_all(funs(mean), na.rm = T))
        sd = as.data.frame(edata[, -edata_cname_id] %>% summarise_all(funs(sd), na.rm = T))
        mds = as.data.frame(edata[, -edata_cname_id] %>% summarise_all(funs(median), na.rm = T))
        pct_obs = as.data.frame(t(apply(edata[,-edata_cname_id], 2, function(x){sum(!is.na(x))/length(x)})))
        min = as.data.frame(edata[, -edata_cname_id] %>% summarise_all(funs(min), na.rm = T))
        max = as.data.frame(edata[, -edata_cname_id] %>% summarise_all(funs(max), na.rm = T))
        
        res_list = list(mean = avg, sd = sd, median = mds, pct_obs = pct_obs, min = min , max = max)
        class(res_list)<- "dataRes"
        
      }else if(!is.null(attr(omicsData, "group_DF"))){
        #if groupvar = NULL and group_designation has been run, groupvar gets set to attr(omicsData, "group_DF")$Group
        edata_melt = reshape2::melt(edata, id.vars = edata_cname)
        names(edata_melt)[2]<- fdata_cname
        edata_melt = merge.data.frame(edata_melt, attr(omicsData, "group_DF"), by = fdata_cname)
        edata_melt = edata_melt[, -which(names(edata_melt) == fdata_cname)]
        
        #calculate metrics for samples belonging to each group
        avg = as.data.frame(edata_melt %>% group_by(Group) %>% summarise(mean = mean(value, na.rm = T)))
        sd = as.data.frame(edata_melt %>% group_by(Group) %>% summarise(sd = sd(value, na.rm = T)))
        mds = as.data.frame(edata_melt %>% group_by(Group) %>% summarise(mds = median(value, na.rm = T)))
        pct_obs = as.data.frame(edata_melt %>% group_by(Group) %>% summarise(pct_obs = sum(!is.na(value))/length(value)))
        min = as.data.frame(edata_melt %>% group_by(Group) %>% summarise(min = min(value, na.rm = T)))
        max = as.data.frame(edata_melt %>% group_by(Group) %>% summarise(max = max(value, na.rm = T)))
        
        res_list = list(mean = avg, sd = sd, median = mds, pct_obs = pct_obs, min = min , max = max)
        class(res_list)<- "dataRes"
      }
    }else if(length(groupvar) == 1){
      ####case where groupvar is provided and has length 1####
      if(any(groupvar %in% names(fdata)) == F) stop("groupvar must be in omicsData f_data")
      temp_fdata = fdata[, c(fdata_cname_id, which(names(fdata) == groupvar))]
      names(temp_fdata)[2]<- "Group"
      
      #rearranging edata
      edata_melt = reshape2::melt(edata, id.vars = edata_cname)
      names(edata_melt)[2]<- fdata_cname
      edata_melt = merge.data.frame(edata_melt, temp_fdata, by = fdata_cname)
      edata_melt = edata_melt[, -which(names(edata_melt) == fdata_cname)]
      
      #calculate metrics for samples belonging to each group
      avg = as.data.frame(edata_melt %>% group_by(Group) %>% summarise(mean = mean(value, na.rm = T)))
      sd = as.data.frame(edata_melt %>% group_by(Group) %>% summarise(sd = sd(value, na.rm = T)))
      mds = as.data.frame(edata_melt %>% group_by(Group) %>% summarise(mds = median(value, na.rm = T)))
      pct_obs = as.data.frame(edata_melt %>% group_by(Group) %>% summarise(pct_obs = sum(!is.na(value))/length(value)))
      min = as.data.frame(edata_melt %>% group_by(Group) %>% summarise(min = min(value, na.rm = T)))
      max = as.data.frame(edata_melt %>% group_by(Group) %>% summarise(max = max(value, na.rm = T)))
      
      res_list = list(mean = avg, sd = sd, median = mds, pct_obs = pct_obs, min = min , max = max)
      class(res_list)<- "dataRes"
      
    }else if(length(groupvar) == 2){
      ####case where length of groupvar is 2####
      if (length(groupvar) > 2) stop("No more than two groupvar can be provided")
      if(any(groupvar %in% names(fdata)) == F) stop("groupvar must be in omicsData f_data")
      
      group_vars = fdata[, names(fdata)%in% groupvar]
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
      
      #rearranging edata
      edata_melt = reshape2::melt(omicsData$e_data, id.vars = edata_cname)
      names(edata_melt)[2]<- fdata_cname
      edata_melt = merge.data.frame(edata_melt, output, by = fdata_cname)
      edata_melt = edata_melt[, -which(names(edata_melt) == fdata_cname)]
      
      #calculate metrics for samples belonging to each group
      avg = as.data.frame(edata_melt %>% group_by(Group) %>% summarise(mean = mean(value, na.rm = T)))
      sd = as.data.frame(edata_melt %>% group_by(Group) %>% summarise(sd = sd(value, na.rm = T)))
      mds = as.data.frame(edata_melt %>% group_by(Group) %>% summarise(mds = median(value, na.rm = T)))
      pct_obs = as.data.frame(edata_melt %>% group_by(Group) %>% summarise(pct_obs = sum(!is.na(value))/length(value)))
      min = as.data.frame(edata_melt %>% group_by(Group) %>% summarise(min = min(value, na.rm = T)))
      max = as.data.frame(edata_melt %>% group_by(Group) %>% summarise(max = max(value, na.rm = T)))
      
      res_list = list(mean = avg, sd = sd, median = mds, pct_obs = pct_obs, min = min , max = max)
      class(res_list)<- "dataRes"
    }
  }
  
  if(by == "molecule"){
    if(is.null(groupvar)){
      if(is.null(attr(omicsData, "group_DF"))){
        #when groupvar is NULL and group_designation is NULL, we calculate metric for each row
        avg = apply(edata[, -edata_cname_id], 1, function(x){if(all(is.na(x))){mean(x)}else{mean(x, na.rm = T)}})
        avg = data.frame(molecule = edata[, edata_cname_id], mean = avg, stringsAsFactors = F)
        sd = apply(edata[, -edata_cname_id], 1, sd, na.rm = T)
        sd = data.frame(molecule = edata[, edata_cname_id], sd = sd, stringsAsFactors = F)
        mds = apply(edata[, -edata_cname_id], 1, median, na.rm = T)
        mds = data.frame(molecule = edata[, edata_cname_id], median = mds, stringsAsFactors = F)
        pct_obs = apply(edata[,-edata_cname_id], 1, function(x){sum(!is.na(x))/length(x)})
        pct_obs = data.frame(molecule = edata[, edata_cname_id], pct_obs = pct_obs, stringsAsFactors = F)
        min = apply(edata[, -edata_cname_id], 1, min, na.rm = T)
        min = data.frame(molecule = edata[, edata_cname_id], min = min, stringsAsFactors = F)
        max = apply(edata[, -edata_cname_id], 1, max, na.rm = T)
        max = data.frame(molecule = edata[, edata_cname_id], max = max, stringsAsFactors = F)
        
        res_list = list(mean = avg, sd = sd, median = mds, pct_obs = pct_obs, min = min , max = max)
        class(res_list)<- "dataRes"
        
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
      
      res_list = list(mean = avg, sd = std_div, median = mds, pct_obs = pct_obs, min = mins , max = maxs) 
      class(res_list)<- "dataRes"
    }
  }
  return(res_list)
}

