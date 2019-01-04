context("Test plotting of filter, corres, and omicsData")
library(pmartR)
library(pmartRdata)
library(testthat)
library(gridExtra)

# sink(file = "logs_test_plotting")
# List of testing objects.  Each has a fake group and custom sample names
obj_list <- lapply(list(pmartRdata::pep_object, pmartRdata::pro_object, pmartRdata::metab_object, pmartRdata::lipid_object), function(x){
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
# test corRes objects

lapply(1:4, function(i){
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
    
    # single plot to test
    cor_plot <- do.call(plot, c(params_basic, params_extra))
    
    # different tests depending on interactive or not
    if(inherits(cor_plot, "d3heatmap")){
      expect_true(all(abs(as.numeric(cor_plot$x$matrix$data) - cor_matrix) < 0.01))
    }
    else if(inherits(cor_plot, "ggplot")) {
      expect_true(all(cor_plot$data$value %in% cor_matrix))
    }
  })
})

#### Test and store plots for filters ###
test_that("filter plot errors", {

  # only test pep_object
  obj <- obj_list[[1]]  %>%
    group_designation(main_effects = "testgroup") %>%
    edata_transform(data_scale = "log2")

# run tests for all filters in a loop
for(i in 1:length(filter_list)){
    
    # create the filter object corresponding to filter type i
    f <- get(filter_list[[i]], envir=asNamespace("pmartR"), mode="function")
    filter_obj <- f(obj)
    
    # parameters that happen for most plot functions
    common_params <- list(
      title_plot = paste(sample(LETTERS, 5), collapse = ""),
      x_lab = paste(sample(LETTERS, 5), collapse = ""),
      y_lab = paste(sample(LETTERS, 5), collapse = ""),
      title_size = sample(8:16, 1),
      x_lab_size = sample(5:16, 1),
      y_lab_size = sample(5:16, 1),
      bw_theme = as.logical(i%%2)
    )
    
    ## bunch of conditional tests based on filter type ##
    
    if(filter_list[[i]] == "molecule_filter"){
      expect_error(plot(filter_obj, min_num = attributes(obj)$data_info$num_samps + 1))
      expect_error(plot(filter_obj, min_num = "Aaa"))
      expect_error(plot(filter_obj, min_num = FALSE))

      temp_plot <- do.call(plot, c(list(filter_object = filter_obj, min_num = 2), common_params))
      expect_true(max(temp_plot$data$frequency_counts) == nrow(filter_obj))
      expect_true(max(temp_plot$data$num_observations) - 1 == max(filter_obj$Num_Observations))
                  
      temp_plot <- do.call(plot, c(list(filter_object = filter_obj, min_num = 2, cumulative = FALSE), common_params))
      expect_true(sum(temp_plot$data$frequency_counts) == nrow(filter_obj))
      expect_true(max(temp_plot$data$num_observations) == max(filter_obj$Num_Observations))
      
    } else if (filter_list[[i]] == "proteomics_filter"){
      expect_error(plot(filter_obj, degen_peps = "ASDF"))
      expect_error(plot(filter_obj, degen_peps = 1))
      expect_error(plot(filter_obj, min_num_peps = filter_object$counts_by_pep$n + 1))
      expect_error(plot(filter_obj, min_num_peps = "AAA"))

      temp_plot <- plot(filter_object = filter_obj, min_num_peps = 3, mapping = "pep_to_pro")
      
      expect_true(max(temp_plot$data$counts) == nrow(filter_obj$counts_by_pro))
      expect_true(temp_plot$data$bins[length(unique(filter_obj$counts_by_pro$n))] == max(filter_obj$counts_by_pro$n))
      
      temp_plot <- plot(filter_object = filter_obj, min_num_peps = 3, mapping = "pro_to_pep")

      expect_true(max(temp_plot$data$counts) == nrow(filter_obj$counts_by_pep))
      expect_true(temp_plot$data$bins[length(unique(filter_obj$counts_by_pep$n))] == max(filter_obj$counts_by_pep$n))

    } else if (filter_list[[i]] == "imdanova_filter"){
      expect_error(plot(filter_obj, min_nonmiss_anova = min(attributes(filter_object)$group_sizes$n_group) + 1))
      expect_error(plot(filter_obj, min_nonmiss_gtest = min(attributes(filter_object)$group_sizes$n_group) + 1))
      expect_error(plot(filter_obj, min_nommiss_anova = "2"))
      expect_error(plot(filter_obj, min_nonmiss_gest = "2"))

      temp_plot <- do.call(plot, c(list(filter_object = filter_obj, min_nonmiss_anova = 3, min_nonmiss_gtest = 3), common_params[names(common_params) != "bw_theme"]))
      expect_equal(filter_obj[,-1] %>% dplyr::filter_all(dplyr::all_vars(. >= 3)) %>% nrow(),
                   temp_plot$data %>% dplyr::filter(Var1 == 3 & Var2 == 3) %>% purrr::pluck('value'))
      
    } else if (filter_list[[i]] == "rmd_filter"){
      expect_error(plot(filter_obj, pvalue_threshold = runif(1,0,1) + 1))
      expect_error(plot(filter_obj, pvalue_threshold = as.character(runif(1,0,1))))
      expect_error(plot(filter_obj, sampleID = 1234))

      pval_thresh = runif(1,0,1)
      
      temp_plot1 <- do.call(plot, c(list(filter_object = filter_obj, pvalue_threshold = pval_thresh), common_params))
      expect_true(all(temp_plot1$data == filter_obj %>% dplyr::filter(pvalue <= pval_thresh)))
      
      temp_plot2 <- do.call(plot, c(list(filter_object = filter_obj, pvalue_threshold = pval_thresh, sampleID = "Infection1"), common_params))
      expect_true(all(data.table::melt(filter_obj, measure = 5:ncol(filter_obj))$value == temp_plot2$data$value))

    } else if (filter_list[[i]] == "cv_filter"){
      expect_error(plot(filter_obj, cv_threshold = max(filter_object$CV_pooled, na.rm = TRUE) + sample(0:1, 1)))
      expect_error(plot(filter_obj, cv_threshold = 0.5))
      expect_error(plot(filter_obj, cv_threshold = as.character(sample(1:(max(filter_object$CV_pooled, na.rm = TRUE)-1)))))

      temp_plot <- do.call(plot, c(list(filter_object = filter_obj, cv_threshold = sample(1:(max(filter_obj$CV_pooled, na.rm = TRUE)-1), 1)), common_params))
      
      expect_true(all(temp_plot$data == filter_obj %>% dplyr::filter(!is.na(CV_pooled))))
    }
  }
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
for(i in 1:length(grouped_objs)){
    omicsData <- grouped_objs[[i]]
    
    # This block of parameter creation doesn't do anything right now
    
    # common_params <- list(
    #   order_by = ifelse(i%%2 == 1, "Condition", "testgroup"),
    #   color_by = ifelse(i%%2 == 1, "Condition", "testgroup"),
    #   facet_by = ifelse(i%%2 == 1, "Condition", "testgroup"),
    #   title_plot = paste(sample(LETTERS, 5), collapse = ""),
    #   x_lab = paste(sample(LETTERS, 5), collapse = ""),
    #   y_lab = paste(sample(LETTERS, 5), collapse = ""),
    #   title_size = sample(8:16, 1),
    #   x_lab_size = sample(5:16, 1),
    #   y_lab_size = sample(5:16, 1),
    #   bw_theme = as.logical(i%%2),
    #   ylimit = c(min(omicsData$e_data[-which(names(omicsData$e_data) == get_edata_cname(omicsData))], na.rm = TRUE),
    #              max(omicsData$e_data[-which(names(omicsData$e_data) == get_edata_cname(omicsData))], na.rm = TRUE)),
    #   use_VizSampNames = as.logical(i%%2)
    # )
    # 
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
  }

#### PCA, normres, edata_summary and missingval_result objects #####

# create lists of objects resulting from applying dim_reduction, normalize_global, edata_summary, and missingval_result to each omicsData object in obj_list #
pcares <- lapply(obj_list, function(obj){
  if(attr(obj, "data_info")$data_scale != "log2"){
    obj <- obj %>% edata_transform(data_scale = "log2")
  }
  
  obj <- obj %>% group_designation(main_effects = "testgroup")

  obj <- dim_reduction(obj)

  obj
  })

normres <- lapply(obj_list, function(obj){
  if(attr(obj, "data_info")$data_scale != "log2"){
    obj <- obj %>% edata_transform(data_scale = "log2")
  }
  
  obj <- obj %>% group_designation(main_effects = "testgroup")

  obj  <- normalize_global(obj , subset_fn = "all", norm_fn = "median")

  obj 
})

summaryres <- lapply(obj_list, function(obj){
  if(attr(obj, "data_info")$data_scale != "log2"){
    obj <- obj %>% edata_transform(data_scale = "log2")
  }
  
  obj <- obj %>% group_designation(main_effects = "testgroup")
  
  suppressWarnings(
    obj  <- edata_summary(obj , by = "molecule", groupvar = "testgroup")
  )
  
  obj 
})

missingvalres <- lapply(obj_list, function(obj){
  if(attr(obj, "data_info")$data_scale != "log2"){
    obj <- obj %>% edata_transform(data_scale = "log2")
  }
  
  obj <- obj %>% group_designation(main_effects = "testgroup")

  obj <- missingval_result(obj)

  obj
})

# end creating lists of results of function calls to omicsData objects #

# plot test loop for misc objects

test_that("data integrity for misc plots", {
  for(i in 1:length(obj_list)){
  
    common_params <- list(
      x_lab = paste(sample(LETTERS, 5), collapse = ""),
      y_lab = paste(sample(LETTERS, 5), collapse = ""),
      title_size = sample(8:16, 1),
      x_lab_size = sample(5:16, 1),
      y_lab_size = sample(5:16, 1),
      bw_theme = as.logical(i%%2),
      legend_position = ifelse(i%%2 == 1, "left", "right")
    )
    
    # bunch of parameters
    metric = sample(c('median', 'sd', 'min', 'max'), 1)
    ncol = sample(1:length(unique(obj_list[[i]]$f_data[,"testgroup"])), 1)
    palette = ifelse(i%%2 == 1, "Set1", "Set2")
    type = ifelse(i%%2 == 1, "bySample", ifelse(i%%4 == 0, "byMolecule", "Both"))
    use_VizSampNames = as.logical(i%%2)
    coordinate_flip = as.logical(i%%2)
    
    # bunch of plots
    pca_plot <- do.call(plot, c(list(pcares[[i]], point_size = sample(4:10, 1)), common_params))
    norm_plot <- do.call(plot, c(list(normres[[i]], combined = TRUE), common_params))
    summary_plot <- plot(summaryres[[i]], metric = metric, ncols = ncol)
    missingval_plot <- plot(naRes_object = missingvalres[[i]], type = type, palette = palette, use_VizSampNames = use_VizSampNames, coordinate_flip = coordinate_flip)
    
    ## tests
    expect_true(all(pca_plot$data %in% cbind(attr(pcares[[i]], "group_DF"), pcares[[i]]["PC1"], pcares[[i]]["PC2"])))
    expect_true(all(summaryres[[i]][[metric]][c(2,3)] %>% unlist() == summary_plot$data$value, na.rm = TRUE))
    
    ## scaling factor to test data integrity
    plot_scaling_factors <- norm_plot[[1]]$data %>% 
      dplyr::group_by(variable) %>% 
      dplyr::summarise_at("value", median) %>% 
      dplyr::select("value")

    expect_true(all(normres[[i]]$parameters$normalization$location == plot_scaling_factors))
  }
})

