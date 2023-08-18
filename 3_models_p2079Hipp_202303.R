# Date: 2022-07-11
# Update:2023-01-17
# S Ogletree
# Description: Models for Greenspace and Telomere Project

# Key packages needed are lme4 and its dependencies

library(lme4)

# set option to not get scientific notation
options(scipen = 999)
# ------------------------------------------------------------------------#
# get the analysis dataset

adata <- readRDS("../data/analysis_data.rds")
# ------------------------------------------------------------------------#
# adjust NDVI measures by 10 for interpretability
adata$NDVI_mean_blkgrp <- adata$NDVI_mean_blkgrp * 10
adata$NDVI_mean_tract <- adata$NDVI_mean_tract * 10
adata$NDVI_99_00_cbg <- adata$`meanNDVI99_00_cbg` * 10
adata$NDVI_99_00_tract <- adata$`meanNDVI99_00_tract` * 10
adata$NDVI_01_02_cbg <- adata$`meanNDVI01_02_cbg` * 10
adata$NDVI_01_02_tract <- adata$`meanNDVI01_02_tract` * 10
# ------------------------------------------------------------------------#
# we need to run the models and then extract the information for output

# bivariate associations with greenspace

b1 <- lm(basepairs ~ NDVI_mean_tract, data = adata)
b2 <- lm(NDVI_mean_tract ~ race_eth  , data = adata)
b3 <- lm(family_pir ~ NDVI_mean_tract, data = adata)
b4 <- lm(NDVI_mean_tract ~ phyact , data = adata)
b5 <- lm(NDI_tract ~ NDVI_mean_tract, data = adata)
b6 <- lm(seg_index_tract ~ NDVI_mean_tract, data = adata)
b7 <- lm(pm25_mean_tract ~ NDVI_mean_tract, data = adata)
b8 <- lm(NDVI_mean_tract ~ tract_redline, data = adata)
# ------------------------------------------------------------------------#

# ------------------------------------------------------------------------#
# linear model, adjusted
# ------------------------------------------------------------------------#

# list for models
lm_list <- list()
# model L1
lm_list[[1]] <- lm(basepairs ~ NDVI_mean_tract + gender + age_years + race_eth + country_born + ed_level + marital_status + hh_income + family_pir, data = adata)
# model L2
lm_list[[2]] <- lm(basepairs ~ NDVI_mean_tract + gender + age_years + race_eth + country_born + ed_level + marital_status + hh_income + family_pir + drinker + smoker + bmi + phyact, data = adata)
# model L3
lm_list[[3]] <- lm(basepairs ~ NDVI_mean_tract + gender + age_years + race_eth + country_born + ed_level + marital_status + hh_income + family_pir + drinker + smoker + bmi + phyact + NDI_tract + seg_index_tract + pm25_mean_tract + tract_redline, data = adata)

# ------------------------------------------------------------------------#
# Multilevel models -------------------------------------------------------
# ------------------------------------------------------------------------#

# First, see how many people we have per grouping
# We want this to be a suitable grouping variable so we want to see a substantial number of people per group
t1 <- sort(table(adata$grouping), decreasing = T)
length(t1) # how many groups are in the dataset
table(t1) # table of counts, we hope to see most groups with a higher count of people in them

# get ICC value

icc_func <- function (object){
  MOD <- summary(object)
  MSB <- MOD[[1]][1, 3]
  MSW <- MOD[[1]][2, 3]
  GSIZE <- (MOD[[1]][2, 1] + (MOD[[1]][1, 1] + 1)) / (MOD[[1]][1, 1] + 1)
  OUT <- (MSB - MSW) / (MSB + ((GSIZE - 1) * MSW))
  return(OUT)
}

# run ANOVA on variables of interest, by group (might be slow to run)
# If this errors out then just skip it and go to line 81
x1 <- aov(basepairs ~ as.factor(grouping), data = adata)

# calculate ICC
round(icc_func(x1), 4)
# we are hoping to see this value > .03 ideally

# list to hold models
mlm_list <- list()
# model M1
mlm_list[[1]] <- lmer(basepairs ~ 1 + (1|grouping), data = adata)
# model M2
mlm_list[[2]] <- lmer(basepairs ~ NDVI_mean_tract + (1|grouping), data = adata)
# model M3
mlm_list[[3]] <- update(mlm_list[[2]], . ~ . + gender + age_years + race_eth + country_born + ed_level + marital_status + hh_income + family_pir)
# model M4
mlm_list[[4]] <- update(mlm_list[[3]], . ~ . + drinker + smoker + bmi + phyact)
# model M5
mlm_list[[5]] <- update(mlm_list[[4]], . ~ . + NDI_tract + seg_index_tract + pm25_mean_tract + tract_redline)

# ------------------------------------------------------------------------#
# Stratified linear models by gender
# ------------------------------------------------------------------------#

# males
male_df <- adata[adata$gender == "male",]
# females
female_df <- adata[adata$gender == "female",]

# list for models
gen_list <- list()

# model L6
gen_list[[1]] <- lm(basepairs ~ NDVI_mean_tract + age_years + race_eth + country_born + ed_level + marital_status + hh_income + family_pir, data = male_df)

gen_list[[2]] <- lm(basepairs ~ NDVI_mean_tract + age_years + race_eth + country_born + ed_level + marital_status + hh_income + family_pir + drinker + smoker + bmi + phyact, data = male_df)

# model L7
gen_list[[3]] <- lm(basepairs ~ NDVI_mean_tract + age_years + race_eth + country_born + ed_level + marital_status + hh_income + family_pir + drinker + smoker + bmi + phyact + NDI_tract + seg_index_tract + tractpm25_1999 + tract_redline, data = male_df)

gen_list[[4]] <- lm(basepairs ~ NDVI_mean_tract + age_years + race_eth + country_born + ed_level + marital_status + hh_income + family_pir, data = female_df)

# model M8
gen_list[[5]] <- lm(basepairs ~ NDVI_mean_tract + age_years + race_eth + country_born + ed_level + marital_status + hh_income + family_pir + drinker + smoker + bmi + phyact, data = female_df)

gen_list[[6]] <- lm(basepairs ~ NDVI_mean_tract + age_years + race_eth + country_born + ed_level + marital_status + hh_income + family_pir + drinker + smoker + bmi + phyact + NDI_tract + seg_index_tract + tractpm25_1999 + tract_redline, data = female_df)


# ----------------------------------------------------------------------------#
# write output ------------------------------------------------------------
# ----------------------------------------------------------------------------#

# output for table 2 - bivariate associations with NDVI
outputlist <- list()
outputlist[[1]] <- "Telomere Length: basepairs ~ NDVI_mean_tract"
outputlist[[2]] <- round(summary(b1)$coefficients[1:2,], 5)
outputlist[[3]] <- "Race/Ethnicity: NDVI_mean_tract ~ race_eth"
outputlist[[4]] <- round(summary(b2)$coefficients[1:5,], 5)
outputlist[[5]] <- "Family PIR: family_pir ~ NDVI_mean_tract"
outputlist[[6]] <- round(summary(b3)$coefficients[1:2,], 5)
outputlist[[7]] <- "Physical Activity: NDVI_mean_tract ~ phyact"
outputlist[[8]] <- round(summary(b4)$coefficients[1:2,], 5)
outputlist[[9]] <- "Neighborhood Deprivation: NDI_tract ~ NDVI_mean_tract"
outputlist[[10]] <- round(summary(b5)$coefficients[1:2,], 5)
outputlist[[11]] <- "Segregation: seg_index_tract ~ NDVI_mean_tract"
outputlist[[12]] <- round(summary(b6)$coefficients[1:2,], 5)
outputlist[[13]] <- "Air Pollution: pm25_mean_tract ~ NDVI_mean_tract"
outputlist[[14]] <- round(summary(b7)$coefficients[1:2,], 5)
outputlist[[15]] <- "Redlining: NDVI_mean_tract ~ tract_redline"
outputlist[[16]] <- round(summary(b8)$coefficients[1:2,], 5)

capture.output(outputlist, file = '../disclosure/p2079_table_2.txt')


# Table 3 -----------------------------------------------------------------

t3_output <- list()

t3_output[[1]] <- "model L1"
t3_output[[2]] <- round(summary(lm_list[[1]])$coefficients[1:length(lm_list[[1]]$coefficients),], 5)
t3_output[[3]] <- "model L2"
t3_output[[4]] <- round(summary(lm_list[[2]])$coefficients[1:length(lm_list[[2]]$coefficients),], 5)
t3_output[[5]] <- "model L3"
t3_output[[6]] <- round(summary(lm_list[[3]])$coefficients[1:length(lm_list[[3]]$coefficients),], 5)

capture.output(t3_output, file = '../disclosure/p2079_table_3.txt')

# table 4 -----------------------------------------------------------------

t4_output <- list()

t4_output[[1]]="model 6"
t4_output[[2]]=summary(mlm_list[[1]])$coefficients
t4_output[[3]]=summary(mlm_list[[1]])$varcor
t4_output[[4]]="model 7"
t4_output[[5]]=summary(mlm_list[[2]])$coefficients
t4_output[[6]]=summary(mlm_list[[2]])$varcor
t4_output[[7]]="model 8"
t4_output[[8]]=summary(mlm_list[[3]])$coefficients
t4_output[[9]]=summary(mlm_list[[3]])$varcor
t4_output[[10]]="model 9"
t4_output[[11]]=summary(mlm_list[[4]])$coefficients
t4_output[[12]]=summary(mlm_list[[4]])$varcor
t4_output[[13]]="model 10"
t4_output[[14]]=summary(mlm_list[[5]])$coefficients
t4_output[[15]]=summary(mlm_list[[5]])$varcor

capture.output(t4_output, file = '../disclosure/p2079_table_4.txt')


# table 6 -----------------------------------------------------------------

t6_output <- list()

t6_output[[1]] <- paste("model L6 Female N=", length(summary(gen_list[[4]])$residuals))
t6_output[[2]] <- round(summary(gen_list[[4]])$coefficients[1:length(gen_list[[4]]$coefficients),], 5)
t6_output[[3]] <- paste("model L6 Male N=", length(summary(gen_list[[1]])$residuals))
t6_output[[4]] <- round(summary(gen_list[[1]])$coefficients[1:length(gen_list[[1]]$coefficients),], 5)
t6_output[[5]] <- paste("model L7 Female N=", length(summary(gen_list[[5]])$residuals))
t6_output[[6]] <- round(summary(gen_list[[5]])$coefficients[1:length(gen_list[[5]]$coefficients),], 5)
t6_output[[7]] <- paste("model L7 Male N=", length(summary(gen_list[[2]])$residuals))
t6_output[[8]] <- round(summary(gen_list[[2]])$coefficients[1:length(gen_list[[2]]$coefficients),], 5)
t6_output[[9]] <- paste("model L8 Female N=", length(summary(gen_list[[6]])$residuals))
t6_output[[10]] <- round(summary(gen_list[[6]])$coefficients[1:length(gen_list[[6]]$coefficients),], 5)
t6_output[[11]] <- paste("model L8 Male N=", length(summary(gen_list[[3]])$residuals))
t6_output[[12]] <- round(summary(gen_list[[3]])$coefficients[1:length(gen_list[[3]]$coefficients),], 5)

capture.output(t6_output, file = '../disclosure/p2079_table_6.txt')

