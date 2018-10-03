context("Test plotting of filter, corres, and omicsData - SAVING PLOTS")
library(pmartR)
library(pmartRdata)
library(testthat)
library(htmlwidgets)
library(gridExtra)

data("pep_object")
sink(file = "logs_test_plotting")
# List of testing objects.  Each has a fake group and custom sample names
obj_list <- lapply(list(pep_object, pro_object, metab_object, lipid_object), function(x){
  x$f_data["testgroup"] <- c(rep(1, floor(nrow(x$f_data)/2)), rep(2, ceiling(nrow(x$f_data)/2)))
  x = custom_sampnames(x, from = 2, to = 14)
})

# filter functions to test
filter_list <- list("molecule_filter", "proteomics_filter", "imdanova_filter", "rmd_filter", "cv_filter")

### Test and store plots for corRes objects ###

test_that("cor_res errors", {
  # test errors
  lapply(obj_list, function(x){
    cor_matrix <- cor_result(x)
    expect_error(plot(cor_matrix, x, title_plot = 555))
    expect_error(plot(cor_matrix, x, colorbar_lim = "A"))
    expect_error(plot(cor_matrix, use_VizSampNames = TRUE))
  })
 
})

# write corRes plots to tests/testthat/images
plot_list <- list()

plot_list <-  lapply(1:4, function(i){
    # test parameters of sucessful plot
    cor_matrix <- cor_result(obj_list[[i]])
    params_basic <- list(corRes_object = cor_matrix, omicsData = obj_list[[i]])
  
    # create extra parameters
    title = paste(sample(LETTERS, 5), collapse = "")
    title_size = sample(8:16, 1)
    xlab = as.logical(i%%2)
    ylab = as.logical((i)%%2)
    interactive = as.logical(i%%2)
  
    params_extra <- list(title_plot = title, title_size = title_size, x_lab = xlab, y_lab = ylab, interactive = interactive)
  
    # png(paste0("plot_", i, ".png"))
    do.call(plot, c(params_basic, params_extra))
    # dev.off()
  })

for (i in 1:4) {
  
  if(inherits(plot_list[[i]], "d3heatmap")){
    saveWidget(plot_list[[i]], paste("corRes_plot_", i, ".html", sep=""))
  }else{
    png(paste("images/corRes_plot_", i, ".png", sep=""))
    print(plot_list[[i]])
    dev.off()
  }
}

#### Test and store plots for filters ###
test_that("filter plot errors", {
  
  obj <- obj_list[[1]] %>% 
    group_designation(main_effects = "testgroup") %>%
    edata_transform(data_scale = "log2")
  
  f <- get(filter_list[[i]], envir=asNamespace("pmartR"), mode="function")
  filter_obj <- f(obj)
  
  common_params <- list(
    title_plot = paste(sample(LETTERS, 5), collapse = ""),
    x_lab = paste(sample(LETTERS, 5), collapse = ""),
    y_lab = paste(sample(LETTERS, 5), collapse = ""),
    title_size = sample(8:16, 1),
    x_lab_size = sample(5:16, 1),
    y_lab_size = sample(5:16, 1),
    bw_theme = as.logical(i%%2)
  )
  
  lapply(1:length(filter_list), function(i){
    
    f <- get(filter_list[[i]], envir=asNamespace("pmartR"), mode="function")
    filter_obj <- f(obj)
    
    if(filter_list[[i]] == "molecule_filter"){
      expect_error(plot(filter_obj, min_num = attributes(obj)$data_info$num_samps + 1))
      expect_error(plot(filter_obj, min_num = "Aaa"))
      expect_error(plot(filter_obj, min_num = FALSE))
      
      png(paste("images/", filter_list[[i]],"_plot_", i, ".png", sep=""))
      print(do.call(plot, c(list(filter_object = filter_obj, min_num = 2), common_params)))
      dev.off()
      
    } else if (filter_list[[i]] == "proteomics_filter"){
      expect_error(plot(filter_obj, degen_peps = "ASDF"))
      expect_error(plot(filter_obj, degen_peps = 1))
      expect_error(plot(filter_obj, min_num_peps = filter_object$counts_by_pep$n + 1))
      expect_error(plot(filter_obj, min_num_peps = "AAA"))
      
      x_lab = paste(sample(LETTERS, 5), collapse = "")
      y_lab = paste(sample(LETTERS, 5), collapse = "")
      title = paste(sample(LETTERS, 5), collapse = "")
      png(paste("images/", filter_list[[i]],"_plot_", i, ".png", sep=""))
      do.call(plot, c(list(filter_object = filter_obj, min_num_peps = 3, degen_peps = FALSE,
                           x_lab_pep = x_lab, x_lab_pro = x_lab, y_lab_pep = y_lab, y_lab_pro = y_lab,
                           title.pep = title, title.pro = title), common_params[c("x_lab_size", "y_lab_size", "title_size", "bw_theme")]))
      dev.off()

    } else if (filter_list[[i]] == "imdanova_filter"){
      expect_error(plot(filter_obj, min_nonmiss_anova = min(attributes(filter_object)$group_sizes$n_group) + 1))
      expect_error(plot(filter_obj, min_nonmiss_gtest = min(attributes(filter_object)$group_sizes$n_group) + 1))
      expect_error(plot(filter_obj, min_nommiss_anova = "2"))
      expect_error(plot(filter_obj, min_nonmiss_gest = "2"))
      
      png(paste("images/", filter_list[[i]],"_plot_", i, ".png", sep=""))
      print(do.call(plot, c(list(filter_object = filter_obj, min_nonmiss_anova = 3, min_nonmiss_gtest = 3), common_params[names(common_params) != "bw_theme"])))
      dev.off()

    } else if (filter_list[[i]] == "rmd_filter"){
      expect_error(plot(filter_obj, pvalue_threshold = runif(1,0,1) + 1))
      expect_error(plot(filter_obj, pvalue_threshold = as.character(runif(1,0,1))))
      expect_error(plot(filter_obj, sampleID = 1234))
      
      png(paste("images/", filter_list[[i]],"_scatter_plot_", i, ".png", sep=""))
      print(do.call(plot, c(list(filter_object = filter_obj, pvalue_threshold = runif(1,0,1)), common_params)))
      dev.off()
      png(paste("images/", filter_list[[i]],"_boxplots_plot_", i, ".png", sep=""))
      print(do.call(plot, c(list(filter_object = filter_obj, pvalue_threshold = runif(1,0,1), sampleID = "sampleID"), common_params)))
      dev.off()

    } else if (filter_list[[i]] == "cv_filter"){
      expect_error(plot(filter_obj, cv_threshold = max(filter_object$CV_pooled, na.rm = TRUE) + sample(0:1, 1)))
      expect_error(plot(filter_obj, cv_threshold = 0.5))
      expect_error(plot(filter_obj, cv_threshold = as.character(sample(1:(max(filter_object$CV_pooled, na.rm = TRUE)-1)))))
      
      png(paste("images/", filter_list[[i]],"_plot_", i, ".png", sep=""))
      print(do.call(plot, c(list(filter_object = filter_obj, cv_threshold = sample(1:(max(filter_obj$CV_pooled, na.rm = TRUE)-1), 1)), common_params)))
      dev.off()
    }
  })
})

# test plotting of objects/transformed objects

# set group designation to fake grouping variable for all objects
grouped_objs <- lapply(obj_list, function(obj){
  res <- obj %>% group_designation(main_effects = "testgroup")
  
  if(attr(res, "data_info")$data_scale == "log2"){
    res
  }
  else(res %>% edata_transform(data_scale = "log2"))})

# error testing taken from test_plot_final.R
lapply(1:length(grouped_objs), function(i){
    omicsData <- grouped_objs[[i]]
    
    common_params <- list(
      order_by = ifelse(i%%2 == 1, "Condition", "testgroup"),
      color_by = ifelse(i%%2 == 1, "Condition", "testgroup"),
      facet_by = ifelse(i%%2 == 1, "Condition", "testgroup"),
      title_plot = paste(sample(LETTERS, 5), collapse = ""),
      x_lab = paste(sample(LETTERS, 5), collapse = ""),
      y_lab = paste(sample(LETTERS, 5), collapse = ""),
      title_size = sample(8:16, 1),
      x_lab_size = sample(5:16, 1),
      y_lab_size = sample(5:16, 1),
      bw_theme = as.logical(i%%2),
      ylimit = c(min(omicsData$e_data[-which(names(omicsData$e_data) == get_edata_cname(omicsData))], na.rm = TRUE), 
                 max(omicsData$e_data[-which(names(omicsData$e_data) == get_edata_cname(omicsData))], na.rm = TRUE)),
      use_VizSampNames = as.logical(i%%2)
    )
    
    input<- c(1, 2, 3)
    mat<- matrix(1:6, nrow = 2, ncol = 3)
    
    test_that("invalid input for order_by argument throws error",{
      expect_that(plot(omicsData, order_by = input), throws_error())   
      expect_that(plot(omicsData, order_by = mat), throws_error())   
      expect_that(plot(omicsData, order_by = 11), throws_error())  
      expect_that(plot(omicsData, order_by = 1.223), throws_error())
      expect_that(plot(omicsData, order_by = -2), throws_error()) 
      expect_that(plot(omicsData, order_by = c("blue")), throws_error())   
      expect_that(plot(omicsData, order_by = c("Condition", "Status")), throws_error()) 
    })
    
    test_that("invalid input for color_by argument throws error",{   
      expect_that(plot(omicsData, color_by = input), throws_error())   
      expect_that(plot(omicsData, color_by = mat), throws_error())   
      expect_that(plot(omicsData, color_by = 11), throws_error()) 
      expect_that(plot(omicsData, color_by = 1.54), throws_error())
      expect_that(plot(omicsData, color_by = -4), throws_error())
      expect_that(plot(omicsData, color_by = c("blue")), throws_error())   
      expect_that(plot(omicsData, color_by = c("Condition", "Status")), throws_error())
    }) 
    
    test_that("invalid input for facet_by argument throws error",{
      expect_that(plot(omicsData, facet_by = input), throws_error())
      expect_that(plot(omicsData, facet_by = mat), throws_error())
      expect_that(plot(omicsData, facet_by = 11), throws_error())
      expect_that(plot(omicsData, facet_by = 1.75), throws_error())
      expect_that(plot(omicsData, facet_by = -1), throws_error())
      expect_that(plot(omicsData, facet_by = c("blue")), throws_error())
      expect_that(plot(omicsData, facet_by = c("blue","green")), throws_error())
    })

    test_that("invalid input for facet_cols argument throws error",{
      expect_that(plot(omicsData, facet_by = NULL,facet_cols = 1), throws_error())
      expect_that(plot(omicsData, facet_by = names(omicsData$f_data[2]),facet_cols = 0), throws_error())
    })

    test_that("invalid input for omicsData argument throws error",{
      expect_that(pmartR::plot.proData(input), throws_error())
      expect_that(pmartR::plot.proData(mat), throws_error())
    })

    test_that("invalid input for legend_position argument throws error",{
      expect_that(plot(omicsData, legend_position = FALSE), throws_error())
    })

    test_that("invalid input for bw_theme argument throws error",{
      expect_that(plot(omicsData, bw_theme = NULL), throws_error())
    })

    test_that("invalid input for title_size argument throws error",{
      expect_that(plot(omicsData, title_size = "five"), throws_error())
    })

    test_that("invalid input for x_lab_size argument throws error",{
      expect_that(plot(omicsData, x_lab_size = "five"), throws_error())
    })

    test_that("invalid input for y_lab_size argument throws error",{
      expect_that(plot(omicsData, y_lab_size = "five"), throws_error())
    })

    test_that("invalid input for ylimit argument throws error",{
      expect_that(plot(omicsData, ylimit = "one"), throws_error())
      expect_that(plot(omicsData, ylimit = 15), throws_error())
      expect_that(plot(omicsData, ylimit = c(1,2,3)), throws_error())
      expect_that(plot(omicsData, ylimit = c("one", "two")), throws_error())
    })
    
    # store plots of current omicsdata object in ../images
    png(paste("images/", class(omicsData),"_plot_", i, ".png", sep=""))
    print(do.call(plot, c(list(omicsData), common_params)))
    dev.off()
    
  })

#### PCA, normres, edata_summary and missingval_result objects #####


# create lists of objects resulting from applying dim_reduction, normalize_global, edata_summary, and missingval_result to each omicsData object in obj_list #
pcares <- lapply(obj_list, function(obj){
  res <- obj %>% group_designation(main_effects = "testgroup")
  
  if(attr(res, "data_info")$data_scale != "log2"){
    res <- res %>% edata_transform(data_scale = "log2")
  }
  
  res <- dim_reduction(res)
  
  res
  })

normres <- lapply(obj_list, function(obj){
  res <- obj %>% group_designation(main_effects = "testgroup")
  
  if(attr(res, "data_info")$data_scale != "log2"){
    res <- res %>% edata_transform(data_scale = "log2")
  }
  
  res <- normalize_global(res, subset_fn = "all", norm_fn = "median")
  
  res
})

summaryres <- lapply(obj_list, function(obj){
  res <- obj %>% group_designation(main_effects = "testgroup")
  
  if(attr(res, "data_info")$data_scale != "log2"){
    res <- res %>% edata_transform(data_scale = "log2")
  }
  
  res <- edata_summary(res, by = "molecule", groupvar = "testgroup")
  
  res
})

missingvalres <- lapply(obj_list, function(obj){
  res <- obj %>% group_designation(main_effects = "testgroup")
  
  if(attr(res, "data_info")$data_scale != "log2"){
    res <- res %>% edata_transform(data_scale = "log2")
  }
  
  res <- missingval_result(res)
  
  res
})

# end creating lists of results of function calls to omicsData objects #

# plot write loop for misc objects

lapply(1:length(obj_list), function(i){
  
  common_params <- list(
    x_lab = paste(sample(LETTERS, 5), collapse = ""),
    y_lab = paste(sample(LETTERS, 5), collapse = ""),
    title_size = sample(8:16, 1),
    x_lab_size = sample(5:16, 1),
    y_lab_size = sample(5:16, 1),
    bw_theme = as.logical(i%%2),
    legend_position = ifelse(i%%2 == 1, "left", "right")
  )
  
  png(paste("images/PCA_plot_", i, ".png", sep=""))
  print(do.call(plot, c(list(pcares[[i]], point_size = sample(4:10, 1)), common_params)))
  dev.off()
  
  png(paste("images/NormRes_plot_", i, ".png", sep=""))
  print(do.call(plot, c(list(normres[[i]], combined = TRUE), common_params)))
  dev.off()
  
  #dataRes params
  metric = if(i%%2 == 1) NULL else sample(c('median', 'sd', 'min', 'max'), 1)
  ncol = sample(1:length(unique(obj_list[[i]]$f_data[,"testgroup"])), 1)
  density = ifelse(i%%2 == 1, TRUE, FALSE)
  
  png(paste("images/Summary_plot_", i, ".png", sep=""))
  print(do.call(plot, list(summaryres[[i]], metric = metric, ncols = ncol, density = density)))
  dev.off()
  
  #missingval params
  palette = ifelse(i%%2 == 1, "Set1", "Set2")
  type = ifelse(i%%2 == 1, "bySample", ifelse(i%%4 == 0, "byMolecule", "Both"))
  use_VizSampNames = as.logical(i%%2)
  coordinate_flip = as.logical(i%%2)
  
  png(paste("images/missingvalres_plot_", i, ".png", sep=""))
  print(do.call(plot, c(list(naRes_object = missingvalres[[i]], type = type, palette = palette, use_VizSampNames = use_VizSampNames, coordinate_flip = coordinate_flip))))
  dev.off()
  
})
sink()