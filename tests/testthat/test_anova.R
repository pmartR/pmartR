library(testthat)
library(pmartR)
library(pmartRdata)

suppressMessages(library(dplyr))
attr(pep_object,"cnames")$fdata_cname <- "SampleID"
mypepData <- edata_transform(omicsData = pep_object, data_scale = "log2")
mypepData <- group_designation(omicsData = mypepData, main_effects = c("Condition"))

##------- Test that the comparisons don't change with/without adjustment -----##
context("Comparison adjustment")
#Now filter and run again
imdanova_Filt <- imdanova_filter(omicsData = mypepData)
mypepData <- applyFilt(filter_object = imdanova_Filt, omicsData = mypepData, min_nonmiss_anova=2)
anova_res <- anova_test(omicsData = mypepData)
anova_res_tukey <- anova_test(omicsData = mypepData, pval_adjust = 'tukey') 
#Should be equivalent to above since only making one comparison
expect_equal(all(anova_res$Comparison==anova_res_tukey$Comparisons),TRUE)

##-------- Test the anova function computes means, and comparisons correctly, single factor ----##
#randomly select a peptide (row of e_data)
context("Test one factor ANOVA - allow unequal var")
anova_res <- anova_test(omicsData = mypepData, equal_var = FALSE)
row <- floor(runif(1,1,nrow(mypepData$e_data)))
ydf <- mypepData$e_data[row,]
ydf <- reshape2::melt(ydf[,-1],id.vars=NULL)
ydf$Treat <- mypepData$f_data$Condition
t_res <- t.test(value~Treat,data=ydf)

myres <- anova_res$Results[row,]
#Check means
expect_equal(c(myres$Mean_Infection,myres$Mean_Mock),unname(t_res$estimate))

#Check variance
uneq_var <- unname(((diff(unname(t_res$estimate))/t_res$statistic)^2)*t_res$parameter)
expect_equal(myres$Variance, uneq_var)

#Check contrast
expect_equal(anova_res$Fold_changes[row,],-unname(diff(t_res$estimate)))

#Check p-values (roundabout way to check approx. dfs)
expect_equal(myres$p_value,t_res$p.value)

#Also make sure fold changes and p-values were properly copied to "Fold_change_pvalues" too
expect_equal(anova_res$Fold_changes[row,],myres$Mean_Infection-myres$Mean_Mock)
expect_equal(anova_res$Fold_change_pvalues[row,],myres$p_value)

##----------------------------------------------##
##---------Assume equal variances---------------##
##----------------------------------------------##

context("Test one factor ANOVA - assume equal var")
anova_res <- anova_test(omicsData = mypepData, equal_var = TRUE)
row <- floor(runif(1,1,nrow(mypepData$e_data)))
ydf <- mypepData$e_data[row,]
ydf <- reshape2::melt(ydf[,-1],id.vars=NULL)
ydf$Treat <- mypepData$f_data$Condition
lm_anova_res <- lm(value~Treat-1,data=ydf)

myres <- anova_res$Results[row,]
#Check means
expect_equal(c(myres$Mean_Infection,myres$Mean_Mock),unname(lm_anova_res$coefficients))

#Check variance
eq_var <- anova(lm_anova_res)[2,3]
expect_equal(myres$Variance, eq_var)

#Check contrast
expect_equal(anova_res$Fold_changes[row,],-diff(unname(lm_anova_res$coefficients)))

#Check p-values (roundabout way to check approx. dfs)
expect_equal(myres$p_value,anova(lm(value~Treat,data=ydf))[1,5])

#Also make sure fold changes and p-values were properly copied to "Fold_change_pvalues" too
expect_equal(anova_res$Fold_changes[row,],myres$Mean_Infection-myres$Mean_Mock)
expect_equal(anova_res$Fold_change_pvalues[row,],myres$p_value)



##-------- Test the anova function computes means, and comparisons correctly, two factors ----##
context("Test two factor ANOVA")
library(OvarianPepdataBPsubset)
suppressWarnings(tcga_ovarian_pepdata_bp <- as.pepData(e_data = tcga_ovarian_pepdata_bp_subset$e_data, f_data = tcga_ovarian_pepdata_bp_subset$f_data, e_meta = tcga_ovarian_pepdata_bp_subset$e_meta,
                                      edata_cname="Peptide", fdata_cname="sampleID", emeta_cname = "Protein", data_scale = "log2"))
suppressWarnings(tcga_ovarian_pepdata_bp <- group_designation(omicsData = tcga_ovarian_pepdata_bp, main_effects = c("vital_status","neoplasm_histologic_grade")))
attr(tcga_ovarian_pepdata_bp,"imdanova") <- data.frame(test_with_anova=NULL)
ovarian_res_twofac <- anova_test(omicsData = tcga_ovarian_pepdata_bp)

#Because there is so much missing data we need to pick a row that has at least
#one observation for at least two leves of each of 2 factors
resample <- TRUE
while(resample){
  #Randomly select row to test
  row <- floor(runif(1,1,nrow(tcga_ovarian_pepdata_bp$e_data)))
  my_res <- ovarian_res_twofac$Results[row,]
  my_res_means <- c(my_res$Mean_Alive_G2,my_res$Mean_Alive_G3,my_res$Mean_Dead_G2,my_res$Mean_Dead_G3)
  
  #Get treatment means and compare to my results
  ydf_twofac <- tcga_ovarian_pepdata_bp$e_data[row,]
  ydf_twofac <- reshape2::melt(ydf_twofac[,-1],id.vars=NULL)
  
  ydf_twofac$Treat <- attr(tcga_ovarian_pepdata_bp,"group_DF")$Group
  ydf_twofac$Fac1 <- attr(tcga_ovarian_pepdata_bp,"group_DF")$neoplasm_histologic_grade
  ydf_twofac$Fac2 <- attr(tcga_ovarian_pepdata_bp,"group_DF")$vital_status
  ydf_twofac_nona <- filter(ydf_twofac,!is.na(value))
  
  if(length(unique(ydf_twofac_nona$Fac1))>1 & length(unique(ydf_twofac_nona$Fac2))>1){
     resample <- FALSE  
  }
}

#ANOVA with and without interaction
res_twofac_int <- lm(value~Fac1*Fac2,data=ydf_twofac)
res_twofac_noint <- lm(value~Fac1+Fac2,data=ydf_twofac)
modelc_p <- anova(res_twofac_noint,res_twofac_int)[2,6]

#If full and reduced model are the same, make p-value 0
if(is.na(modelc_p)){
  modelc_p <- 0
}

if(modelc_p<0.05){
  #Full model is preferred, compare means and variance estimates (res_twofac_int)
  ydf_res <- ydf_twofac%>%group_by(Treat)%>%summarize(Est=mean(value,na.rm=T))
  to_remove <- which(is.na(ydf_res$Est))
  if(length(to_remove)>0){
    my_res_means <- my_res_means[-to_remove]
    ydf_res <- ydf_res[-to_remove,]
  }
  
  #Compare means
  expect_equal(my_res_means,ydf_res$Est)
  
  #Compare sigma estimate
  expect_equal(my_res$Variance,anova(res_twofac_int)[3,3])
  
  #Compare contrasts
  #my_res_cont <- unname(data.matrix(ovarian_res_twofac$Fold_changes[row,]))[1,]
  #anova_res_cont <- c(-unname(res_twofac_int$coefficients[2]),)
                      
}else{
  #Reduced model preferred, compare means and variance estimates (res_twofac_noint)
  
  Xmatrix <- matrix(c(1,0,0,
                      1,1,0,
                      1,0,1,
                      1,1,1),ncol=3,byrow=TRUE)
  model_ests <- as.vector((Xmatrix%*%res_twofac_noint$coefficients))
  
  #Compare means
  expect_equal(my_res_means,model_ests)
  
  #Compare variance estimates
  expect_equal(my_res$Variance, anova(res_twofac_noint)[3,3])
}


