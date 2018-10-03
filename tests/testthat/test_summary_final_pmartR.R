context("output tests for summary()")
#testthat function for summary.proData()

library(testthat)
library(pmartR)
library(pmartRdata)
library(dplyr)

data("pro_object")

#regular and group object to test with
omicsData <- pro_object
omicsData_grouped <- group_designation(omicsData, main_effects = "Condition")

# get summary data frame and first column.  needs to be converted from factor char and then to numeric where appropriate
result <- summary(omicsData)
vals <- as.character(result[,1])

# 
test_that("dataframe is returned, should be at least 6 rows", {
  expect_that(result, is_a("data.frame"))
  expect_true(nrow(result) >= 6)
})

# check consistency of table output
test_that("result matches attributes of omicsData",{
  lapply(list(omicsData, omicsData_grouped), function(omicsData){
    expect_equal(vals[1], class(omicsData))
    expect_equal(as.numeric(vals[2]), attr(omicsData, "data_info")$num_samps)
    expect_equal(as.numeric(vals[3]), attr(omicsData, "data_info")$num_edata)
    expect_true(if(vals[4] == "NA") is.null(attr(omicsData, "data_info")$num_emeta) else as.numeric(vals[4]) == attr(omicsData, "data_info")$num_emeta) 
    expect_equal(as.numeric(vals[5]), attr(omicsData, "data_info")$num_miss_obs)
    expect_equal(as.numeric(vals[6]), round(attr(omicsData, "data_info")$prop_missing, 3))
    
    # a grouped omicsData object will have extra rows indicating group sizes, test they are correct
    if(nrow(result) > 6){
      # number per group in original object
      group_counts <- attr(omicsData, "group_DF") %>% group_by(Group) %>% summarise(n = n())
      
      for(i in 7:nrow(result)){
        expect_equal(as.numeric(vals[i]), group_counts$n[i-6])
      }
    }
  })
})

### Test statRes summary ###

omicsData <- pro_object
omicsData$f_data$fakegroup <- c(rep("A",3), rep("B",4), rep("C",4))

# Weird groups
# omicsData$f_data$fakegroup <- 1:12
# omicsData$f_data$fakegroup2 <- c(1:8, 1, 9:11)

omicsData_fakegroup <- group_designation(omicsData = omicsData, main_effects = c("fakegroup"))
omicsData_conditiongroup <- group_designation(omicsData = omicsData, main_effects = c("Condition"))

imdanovafilt_fakegroup <- imdanova_filter(omicsData = omicsData_fakegroup)
imdanovafilt_conditiongroup <- imdanova_filter(omicsData = omicsData_conditiongroup)
omicsData_fakegroup <- applyFilt(filter_object = imdanovafilt_fakegroup, omicsData = omicsData_fakegroup, min_nonmiss_anova=2)
omicsData_conditiongroup <- applyFilt(filter_object = imdanovafilt_conditiongroup, omicsData = omicsData_conditiongroup, min_nonmiss_anova=2)

sink(file = "logs_statres_summary")
for(test in list("anova", "gtest", "comb")){
  for(omicsData in list(omicsData_fakegroup, omicsData_conditiongroup)){
    statRes_obj <- imd_anova(omicsData = omicsData, test_method = test)
    statRes_summary <- summary(statRes_obj)
    flags <- statRes_obj$Flags
    
    for(comparison in levels(statRes_summary$sig_table$Comparison)){
      expect_equal(sum(flags[,comparison] == -2), statRes_summary$sig_table %>% filter(Comparison == comparison) %>% {.$`Negative:G-test`})
      expect_equal(sum(flags[,comparison] == -1), statRes_summary$sig_table %>% filter(Comparison == comparison) %>% {.$`Negative:ANOVA`})
      expect_equal(sum(flags[,comparison] == 1), statRes_summary$sig_table %>% filter(Comparison == comparison) %>% {.$`Positive:ANOVA`})
      expect_equal(sum(flags[,comparison] == 2), statRes_summary$sig_table %>% filter(Comparison == comparison) %>% {.$`Positive:G-test`})
    }
  }
}
sink()

