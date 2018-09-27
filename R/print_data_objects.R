#' print.pepData
#' 
#' For printing an S3 object of type 'pepData':
#' 
#'@rdname print-pepData
#'@export
#'
print.pepData<- function(pepData){
  if(!inherits(pepData, "pepData")) stop("pep_object must be of the class 'pepData'")
  
  e_data<- as.data.frame(lapply(pepData$e_data, as.character), stringsAsFactors = FALSE)
  f_data<- as.data.frame(lapply(pepData$f_data, as.character), stringsAsFactors = FALSE)
  edata_ncols<- ncol(e_data)
  fdata_ncols<- ncol(f_data)
  blank_row = rep("----", 5)
  
  if(!is.null(pepData$e_meta)){
    e_meta<- as.data.frame(lapply(pepData$e_meta, as.character), stringsAsFactors = FALSE)
    emeta_ncols<- ncol(e_meta)
    
    if(nrow(e_data) >= 9){
      e_data_head = head(e_data, 4)[, 1:min(edata_ncols, 5)]
      e_data_tail = tail(e_data, 4)[, 1:min(edata_ncols, 5)]
      edata = rbind(e_data_head, blank_row, e_data_tail)
    }else{
      edata = e_data
    }

    if(nrow(f_data) >= 9){
      f_data_head = head(f_data, 4)[, 1:min(fdata_ncols, 5)]
      f_data_tail = tail(f_data, 4)[, 1:min(fdata_ncols, 5)]
      fdata = rbind(f_data_head, blank_row, f_data_tail)
    }else{
      fdata = f_data
    }

    if(nrow(e_meta) >= 9){
      e_meta_head = head(e_meta, 4)[, 1:min(emeta_ncols, 5)]
      e_meta_tail = tail(e_meta, 4)[, 1:min(emeta_ncols, 5)]
      emeta = rbind(e_meta_head, blank_row, e_meta_tail)
    }else{
      emeta = e_meta
    }

    if(edata_ncols > 5) message("only first 5 columns are shown")
    cat("e_data\n")
    cat(capture.output(edata), sep = "\n")
    cat("\n")
    
    if(fdata_ncols > 5) message("only first 5 columns are shown")
    cat("f_data\n")
    cat(capture.output(fdata), sep = "\n")
    cat("\n")
    
    if(emeta_ncols > 5) message("only first 5 columns are shown")
    cat("e_meta\n")
    cat(capture.output(emeta), sep = "\n")
    cat("\n")
    
  }else{
    if(nrow(e_data) >= 9){
      e_data_head = head(e_data, 4)[, 1:min(edata_ncols, 5)]
      e_data_tail = tail(e_data, 4)[, 1:min(edata_ncols, 5)]
      edata = rbind(e_data_head, blank_row, e_data_tail)
    }else{
      edata = e_data
    }
    
    if(nrow(f_data) >= 9){
      f_data_head = head(f_data, 4)[, 1:min(fdata_ncols, 5)]
      f_data_tail = tail(f_data, 4)[, 1:min(fdata_ncols, 5)]
      fdata = rbind(f_data_head, blank_row, f_data_tail)
    }else{
      fdata = f_data
    }
    
    if(edata_ncols > 5) message("only first 5 columns are shown")
    cat("e_data\n")
    cat(capture.output(edata), sep = "\n")
    cat("\n")
    
    if(fdata_ncols > 5) message("only first 5 columns are shown")
    cat("f_data\n")
    cat(capture.output(fdata), sep = "\n")
    cat("\n")
  }
}

#' print.metabData
#' 
#' For printing an S3 object of type 'metabData':
#' 
#'@rdname print-metabData
#'@export
#'
print.metabData<- function(metabData){
  if(!inherits(metabData, "metabData")) stop("metab_object must be of the class 'metabData'")
  
  e_data<- as.data.frame(lapply(metabData$e_data, as.character), stringsAsFactors = FALSE)
  f_data<- as.data.frame(lapply(metabData$f_data, as.character), stringsAsFactors = FALSE)
  edata_ncols<- ncol(e_data)
  fdata_ncols<- ncol(f_data)
  blank_row = rep("----", 5)
  
  if(!is.null(metabData$e_meta)){
    e_meta<- as.data.frame(lapply(metabData$e_meta, as.character), stringsAsFactors = FALSE)
    emeta_ncols<- ncol(e_meta)
    
    if(nrow(e_data) >= 9){
      e_data_head = head(e_data, 4)[, 1:min(edata_ncols, 5)]
      e_data_tail = tail(e_data, 4)[, 1:min(edata_ncols, 5)]
      edata = rbind(e_data_head, blank_row, e_data_tail)
    }else{
      edata = e_data
    }
    
    if(nrow(f_data) >= 9){
      f_data_head = head(f_data, 4)[, 1:min(fdata_ncols, 5)]
      f_data_tail = tail(f_data, 4)[, 1:min(fdata_ncols, 5)]
      fdata = rbind(f_data_head, blank_row, f_data_tail)
    }else{
      fdata = f_data
    }
    
    if(nrow(e_meta) >= 9){
      e_meta_head = head(e_meta, 4)[, 1:min(emeta_ncols, 5)]
      e_meta_tail = tail(e_meta, 4)[, 1:min(emeta_ncols, 5)]
      emeta = rbind(e_meta_head, blank_row, e_meta_tail)
    }else{
      emeta = e_meta
    }
    
    if(edata_ncols > 5) message("only first 5 columns are shown")
    cat("e_data\n")
    cat(capture.output(edata), sep = "\n")
    cat("\n")
    
    if(fdata_ncols > 5) message("only first 5 columns are shown")
    cat("f_data\n")
    cat(capture.output(fdata), sep = "\n")
    cat("\n")
    
    if(emeta_ncols > 5) message("only first 5 columns are shown")
    cat("e_meta\n")
    cat(capture.output(emeta), sep = "\n")
    cat("\n")
    
  }else{
    if(nrow(e_data) >= 9){
      e_data_head = head(e_data, 4)[, 1:min(edata_ncols, 5)]
      e_data_tail = tail(e_data, 4)[, 1:min(edata_ncols, 5)]
      edata = rbind(e_data_head, blank_row, e_data_tail)
    }else{
      edata = e_data
    }
    
    if(nrow(f_data) >= 9){
      f_data_head = head(f_data, 4)[, 1:min(fdata_ncols, 5)]
      f_data_tail = tail(f_data, 4)[, 1:min(fdata_ncols, 5)]
      fdata = rbind(f_data_head, blank_row, f_data_tail)
    }else{
      fdata = f_data
    }
    
    if(edata_ncols > 5) message("only first 5 columns are shown")
    cat("e_data\n")
    cat(capture.output(edata), sep = "\n")
    cat("\n")
    
    if(fdata_ncols > 5) message("only first 5 columns are shown")
    cat("f_data\n")
    cat(capture.output(fdata), sep = "\n")
    cat("\n")
  }
}

#' print.proData
#' 
#' For printing an S3 object of type 'proData':
#' 
#'@rdname print-proData
#'@export
#'
print.proData<- function(proData){
  if(!inherits(proData, "proData")) stop("pro_object must be of the class 'proData'")
  
  e_data<- as.data.frame(lapply(proData$e_data, as.character), stringsAsFactors = FALSE)
  f_data<- as.data.frame(lapply(proData$f_data, as.character), stringsAsFactors = FALSE)
  edata_ncols<- ncol(e_data)
  fdata_ncols<- ncol(f_data)
  blank_row = rep("----", 5)
  
  if(!is.null(proData$e_meta)){
    e_meta<- as.data.frame(lapply(proData$e_meta, as.character), stringsAsFactors = FALSE)
    emeta_ncols<- ncol(e_meta)
    
    if(nrow(e_data) >= 9){
      e_data_head = head(e_data, 4)[, 1:min(edata_ncols, 5)]
      e_data_tail = tail(e_data, 4)[, 1:min(edata_ncols, 5)]
      edata = rbind(e_data_head, blank_row, e_data_tail)
    }else{
      edata = e_data
    }
    
    if(nrow(f_data) >= 9){
      f_data_head = head(f_data, 4)[, 1:min(fdata_ncols, 5)]
      f_data_tail = tail(f_data, 4)[, 1:min(fdata_ncols, 5)]
      fdata = rbind(f_data_head, blank_row, f_data_tail)
    }else{
      fdata = f_data
    }
    
    if(nrow(e_meta) >= 9){
      e_meta_head = head(e_meta, 4)[, 1:min(emeta_ncols, 5)]
      e_meta_tail = tail(e_meta, 4)[, 1:min(emeta_ncols, 5)]
      emeta = rbind(e_meta_head, blank_row, e_meta_tail)
    }else{
      emeta = e_meta
    }
    
    if(edata_ncols > 5) message("only first 5 columns are shown")
    cat("e_data\n")
    cat(capture.output(edata), sep = "\n")
    cat("\n")
    
    if(fdata_ncols > 5) message("only first 5 columns are shown")
    cat("f_data\n")
    cat(capture.output(fdata), sep = "\n")
    cat("\n")
    
    if(emeta_ncols > 5) message("only first 5 columns are shown")
    cat("e_meta\n")
    cat(capture.output(emeta), sep = "\n")
    cat("\n")
    
  }else{
    if(nrow(e_data) >= 9){
      e_data_head = head(e_data, 4)[, 1:min(edata_ncols, 5)]
      e_data_tail = tail(e_data, 4)[, 1:min(edata_ncols, 5)]
      edata = rbind(e_data_head, blank_row, e_data_tail)
    }else{
      edata = e_data
    }
    
    if(nrow(f_data) >= 9){
      f_data_head = head(f_data, 4)[, 1:min(fdata_ncols, 5)]
      f_data_tail = tail(f_data, 4)[, 1:min(fdata_ncols, 5)]
      fdata = rbind(f_data_head, blank_row, f_data_tail)
    }else{
      fdata = f_data
    }
    
    if(edata_ncols > 5) message("only first 5 columns are shown")
    cat("e_data\n")
    cat(capture.output(edata), sep = "\n")
    cat("\n")
    
    if(fdata_ncols > 5) message("only first 5 columns are shown")
    cat("f_data\n")
    cat(capture.output(fdata), sep = "\n")
    cat("\n")
  }
}

#' print.lipidData
#' 
#' For printing an S3 object of type 'lipidData':
#' 
#'@rdname print-lipidData
#'@export
#'
print.lipidData<- function(lipidData){
  if(!inherits(lipidData, "lipidData")) stop("lipid_object must be of the class 'lipidData'")
  
  e_data<- as.data.frame(lapply(lipidData$e_data, as.character), stringsAsFactors = FALSE)
  f_data<- as.data.frame(lapply(lipidData$f_data, as.character), stringsAsFactors = FALSE)
  edata_ncols<- ncol(e_data)
  fdata_ncols<- ncol(f_data)
  blank_row = rep("----", 5)
  
  if(!is.null(lipidData$e_meta)){
    e_meta<- as.data.frame(lapply(lipidData$e_meta, as.character), stringsAsFactors = FALSE)
    emeta_ncols<- ncol(e_meta)
    
    if(nrow(e_data) >= 9){
      e_data_head = head(e_data, 4)[, 1:min(edata_ncols, 5)]
      e_data_tail = tail(e_data, 4)[, 1:min(edata_ncols, 5)]
      edata = rbind(e_data_head, blank_row, e_data_tail)
    }else{
      edata = e_data
    }
    
    if(nrow(f_data) >= 9){
      f_data_head = head(f_data, 4)[, 1:min(fdata_ncols, 5)]
      f_data_tail = tail(f_data, 4)[, 1:min(fdata_ncols, 5)]
      fdata = rbind(f_data_head, blank_row, f_data_tail)
    }else{
      fdata = f_data
    }
    
    if(nrow(e_meta) >= 9){
      e_meta_head = head(e_meta, 4)[, 1:min(emeta_ncols, 5)]
      e_meta_tail = tail(e_meta, 4)[, 1:min(emeta_ncols, 5)]
      emeta = rbind(e_meta_head, blank_row, e_meta_tail)
    }else{
      emeta = e_meta
    }
    
    if(edata_ncols > 5) message("only first 5 columns are shown")
    cat("e_data\n")
    cat(capture.output(edata), sep = "\n")
    cat("\n")
    
    if(fdata_ncols > 5) message("only first 5 columns are shown")
    cat("f_data\n")
    cat(capture.output(fdata), sep = "\n")
    cat("\n")
    
    if(emeta_ncols > 5) message("only first 5 columns are shown")
    cat("e_meta\n")
    cat(capture.output(emeta), sep = "\n")
    cat("\n")
    
  }else{
    if(nrow(e_data) >= 9){
      e_data_head = head(e_data, 4)[, 1:min(edata_ncols, 5)]
      e_data_tail = tail(e_data, 4)[, 1:min(edata_ncols, 5)]
      edata = rbind(e_data_head, blank_row, e_data_tail)
    }else{
      edata = e_data
    }
    
    if(nrow(f_data) >= 9){
      f_data_head = head(f_data, 4)[, 1:min(fdata_ncols, 5)]
      f_data_tail = tail(f_data, 4)[, 1:min(fdata_ncols, 5)]
      fdata = rbind(f_data_head, blank_row, f_data_tail)
    }else{
      fdata = f_data
    }
    
    if(edata_ncols > 5) message("only first 5 columns are shown")
    cat("e_data\n")
    cat(capture.output(edata), sep = "\n")
    cat("\n")
    
    if(fdata_ncols > 5) message("only first 5 columns are shown")
    cat("f_data\n")
    cat(capture.output(fdata), sep = "\n")
    cat("\n")
  }
}


#' print.dataRes
#' 
#' For printing an S3 object of class 'dataRes':
#' 
#'@rdname print-dataRes
#'@export
#'
print.dataRes<- function(dataRes){

  if((!is.null(attr(dataRes, "groupvar")) & attr(dataRes, "by") == "molecule") | (is.null(attr(dataRes, "groupvar")) & !is.null(attr(dataRes, "group_DF")) & (attr(dataRes, "by") == "molecule"))){
      
      #extract item from dataRes, n_per_grp
      n_per_grp<- as.data.frame(lapply(dataRes$n_per_grp, as.character), stringsAsFactors = FALSE)
      
      #choosing how many rows to print
      if(nrow(n_per_grp) <= 10){
        cat("n_per_grp\n")
        print(n_per_grp)
      }else{
        head_n_per_grp = head(n_per_grp, 5)[, 1:min(ncol(n_per_grp), 5)]
        tail_n_per_grp = tail(n_per_grp, 5)[, 1:min(ncol(n_per_grp), 5)]
        blank_row = rep("----", min(ncol(n_per_grp), 5))
        
        print_n_per_grp<- rbind(head_n_per_grp, blank_row, tail_n_per_grp)
        
        if(ncol(n_per_grp) > 5) message("only first 5 columns are shown")
        cat("n_per_grp\n")
        cat(capture.output(print_n_per_grp), sep = "\n")
        cat("\n")
      }
      
      #extract item from dataRes, mean
      mean<- as.data.frame(lapply(dataRes$mean, as.character), stringsAsFactors = FALSE)
      
      #choosing how many rows to print
      if(nrow(mean) <= 10){
        print(mean)
      }else{
        head_mean = head(mean, 5)[, 1:min(ncol(mean), 5)]
        tail_mean = tail(mean, 5)[, 1:min(ncol(mean), 5)]
        blank_row = rep("----", min(ncol(mean), 5))
        
        print_mean<- rbind(head_mean, blank_row, tail_mean)
        
        if(ncol(mean) > 5) message("only first 5 columns are shown")
        cat("mean\n")
        cat(capture.output(print_mean), sep = "\n")
        cat("\n")
      }
      
      #extract item from dataRes, std_dev
      sd<- as.data.frame(lapply(dataRes$sd, as.character), stringsAsFactors = FALSE)
      
      #choosing how many rows to print
      if(nrow(sd) <= 10){
        print(sd)
      }else{
        head_sd = head(sd, 5)[, 1:min(ncol(sd), 5)]
        tail_sd = tail(sd, 5)[, 1:min(ncol(sd), 5)]
        blank_row = rep("----", min(ncol(sd), 5))
        
        print_sd<- rbind(head_sd, blank_row, tail_sd)
        
        if(ncol(sd) > 5) message("only first 5 columns are shown")
        cat("std_dev\n")
        cat(capture.output(print_sd), sep = "\n")
        cat("\n")
      }
      
      #extract item from dataRes, median
      median<- as.data.frame(lapply(dataRes$median, as.character), stringsAsFactors = FALSE)
      
      #choosing how many rows to print
      if(nrow(median) <= 10){
        print(median)
      }else{
        head_median = head(median, 5)[, 1:min(ncol(median), 5)]
        tail_median = tail(median, 5)[, 1:min(ncol(median), 5)]
        blank_row = rep("----", min(ncol(median), 5))
        
        print_median<- rbind(head_median, blank_row, tail_median)
        
        if(ncol(median) > 5) message("only first 5 columns are shown")
        cat("median\n")
        cat(capture.output(print_median), sep = "\n")
        cat("\n")
      }
      
      #extract item from dataRes, pct_obs
      pct_obs<- as.data.frame(lapply(dataRes$pct_obs, as.character), stringsAsFactors = FALSE)
      
      #choosing how many rows to print
      if(nrow(pct_obs) <= 10){
        print(pct_obs)
      }else{
        head_pct_obs = head(pct_obs, 5)[, 1:min(ncol(pct_obs), 5)]
        tail_pct_obs = tail(pct_obs, 5)[, 1:min(ncol(pct_obs), 5)]
        blank_row = rep("----", min(ncol(pct_obs), 5))
        
        print_pct_obs<- rbind(head_pct_obs, blank_row, tail_pct_obs)
        
        if(ncol(pct_obs) > 5) message("only first 5 columns are shown")
        cat("pct_obs\n")
        cat(capture.output(print_pct_obs), sep = "\n")
        cat("\n")
      }
      
      #extract item from dataRes, min
      min<- as.data.frame(lapply(dataRes$min, as.character), stringsAsFactors = FALSE)
      
      #choosing how many rows to print
      if(nrow(min) <= 10){
        print(min)
      }else{
        head_min = head(min, 5)[, 1:min(ncol(min), 5)]
        tail_min = tail(min, 5)[, 1:min(ncol(min), 5)]
        blank_row = rep("----", min(ncol(min), 5))
        
        print_min<- rbind(head_min, blank_row, tail_min)
        
        if(ncol(min) > 5) message("only first 5 columns are shown")
        cat("minimum\n")
        cat(capture.output(print_min), sep = "\n")
        cat("\n")
      }
      
      #extract item from dataRes, max
      max<- as.data.frame(lapply(dataRes$max, as.character), stringsAsFactors = FALSE)
      
      #choosing how many rows to print
      if(nrow(max) <= 10){
        print(max)
      }else{
        head_max = head(max, 5)[, 1:min(ncol(max), 5)]
        tail_max = tail(max, 5)[, 1:min(ncol(max), 5)]
        blank_row = rep("----", min(ncol(max), 5))
        
        print_max<- rbind(head_max, blank_row, tail_max)
        
        if(ncol(max) > 5) message("only first 5 columns are shown")
        cat("maximum\n")
        cat(capture.output(print_max), sep = "\n")
        cat("\n")
      }
    }else{
    #extract item from dataRes, mean
    mean<- as.data.frame(lapply(dataRes$mean, as.character), stringsAsFactors = FALSE)
    
    #choosing how many rows to print
    if(nrow(mean) <= 10){
      print(mean)
    }else{
      head_mean = head(mean, 5)[, 1:min(ncol(mean), 5)]
      tail_mean = tail(mean, 5)[, 1:min(ncol(mean), 5)]
      blank_row = rep("----", min(ncol(mean), 5))
      
      print_mean<- rbind(head_mean, blank_row, tail_mean)
      
      if(ncol(mean) > 5) message("only first 5 columns are shown")
      cat("mean\n")
      cat(capture.output(print_mean), sep = "\n")
      cat("\n")
    }
    
    #extract item from dataRes, std_dev
    sd<- as.data.frame(lapply(dataRes$sd, as.character), stringsAsFactors = FALSE)
    
    #choosing how many rows to print
    if(nrow(sd) <= 10){
      print(sd)
    }else{
      head_sd = head(sd, 5)[, 1:min(ncol(sd), 5)]
      tail_sd = tail(sd, 5)[, 1:min(ncol(sd), 5)]
      blank_row = rep("----", min(ncol(sd), 5))
      
      print_sd<- rbind(head_sd, blank_row, tail_sd)
      
      if(ncol(sd) > 5) message("only first 5 columns are shown")
      cat("std_dev\n")
      cat(capture.output(print_sd), sep = "\n")
      cat("\n")
    }
    
    #extract item from dataRes, median
    median<- as.data.frame(lapply(dataRes$median, as.character), stringsAsFactors = FALSE)
    
    #choosing how many rows to print
    if(nrow(median) <= 10){
      print(median)
    }else{
      head_median = head(median, 5)[, 1:min(ncol(median), 5)]
      tail_median = tail(median, 5)[, 1:min(ncol(median), 5)]
      blank_row = rep("----", min(ncol(median), 5))
      
      print_median<- rbind(head_median, blank_row, tail_median)
      
      if(ncol(median) > 5) message("only first 5 columns are shown")
      cat("median\n")
      cat(capture.output(print_median), sep = "\n")
      cat("\n")
    }
    
    #extract item from dataRes, pct_obs
    pct_obs<- as.data.frame(lapply(dataRes$pct_obs, as.character), stringsAsFactors = FALSE)
    
    #choosing how many rows to print
    if(nrow(pct_obs) <= 10){
      print(pct_obs)
    }else{
      head_pct_obs = head(pct_obs, 5)[, 1:min(ncol(pct_obs), 5)]
      tail_pct_obs = tail(pct_obs, 5)[, 1:min(ncol(pct_obs), 5)]
      blank_row = rep("----", min(ncol(pct_obs), 5))
      
      print_pct_obs<- rbind(head_pct_obs, blank_row, tail_pct_obs)
      
      if(ncol(pct_obs) > 5) message("only first 5 columns are shown")
      cat("pct_obs\n")
      cat(capture.output(print_pct_obs), sep = "\n")
      cat("\n")
    }
    
    #extract item from dataRes, min
    min<- as.data.frame(lapply(dataRes$min, as.character), stringsAsFactors = FALSE)
    
    #choosing how many rows to print
    if(nrow(min) <= 10){
      print(min)
    }else{
      head_min = head(min, 5)[, 1:min(ncol(min), 5)]
      tail_min = tail(min, 5)[, 1:min(ncol(min), 5)]
      blank_row = rep("----", min(ncol(min), 5))
      
      print_min<- rbind(head_min, blank_row, tail_min)
      
      if(ncol(min) > 5) message("only first 5 columns are shown")
      cat("minimum\n")
      cat(capture.output(print_min), sep = "\n")
      cat("\n")
    }
    
    #extract item from dataRes, max
    max<- as.data.frame(lapply(dataRes$max, as.character), stringsAsFactors = FALSE)
    
    #choosing how many rows to print
    if(nrow(max) <= 10){
      print(max)
    }else{
      head_max = head(max, 5)[, 1:min(ncol(max), 5)]
      tail_max = tail(max, 5)[, 1:min(ncol(max), 5)]
      blank_row = rep("----", min(ncol(max), 5))
      
      print_max<- rbind(head_max, blank_row, tail_max)
      
      if(ncol(max) > 5) message("only first 5 columns are shown")
      cat("maximum\n")
      cat(capture.output(print_max), sep = "\n")
      cat("\n")
    }
  }
}
  
 

#' print.normRes
#' 
#' For printing an S3 object of type 'normRes':
#' 
#'@rdname print-normRes
#'@export
#'
print.normRes<- function(normRes){
  
  attr(normRes, "class")<- NULL
  attr(normRes, "omicsData")<- NULL
  
  print(normRes)
}

