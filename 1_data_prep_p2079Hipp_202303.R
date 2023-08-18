# Date: 2022-06-17
# Update: 2023-02-17
# S Ogletree
# Description: Script 1 - label values and calculate drinker/smoker/group
# This will prepare the dataset for running further scripts to get descriptive statistics and models

# the merged dataset
# [need to locate it and see what it is named]
# Should be in 'transfer' folder

# if provided dataset is in Stata or SAS format, use package 'haven'
# Double check this file path
mdata <- haven::read_sas("../transfer/p2079_hipp.sas7bdat")
str(mdata)
# all columns should be num (numeric) EXCEPT RC2K, STM2K, UR, tract_redline & blkgrp_redline
# If they are not then run this - verify col numbers with names(mdata)
mdata[,c(1:34, 38:53, 55, 57)] <- sapply(mdata[,c(1:34, 38:53, 55, 57)],as.numeric)


# ------------------------------------------------------------------------#
# rename variables --------------------------------------------------------
# ------------------------------------------------------------------------#

name_orig <- names(mdata)
# dput(names(mdata)) in the console will print the names as a vector
# Examine this set of column names, it should look like the one below
# make sure that the names we need changed are correct

# changes to make: RIAGENDR=gender, RIDAGEYR=age_years, RIDRETH1=race_eth, DMDBORN=country_born, DMDEDUC2=ed_level, DMDMARTL=marital_status, INDHHINC=hh_income, INDFMPIR=family_pir, BMXBMI=bmi
# Based on what Jaylan sent this should be close to correct
name_new <- c("SEQN", "survey_period", "gender", "age_years", "race_eth", "country_born", 
              "ed_level", "marital_status", "hh_income", "family_pir", "WTINT4YR", "WTMEC4YR", 
              "ALQ120Q", "ALQ120U", "PAD020", "PAQ050Q", "PAQ050U", "PAD080", 
              "PAQ100", "PAD120", "PAD160", "PAQ180", "PAD200", "PAD320", "SMQ020", 
              "SMQ040", "SMD090", "bmi", "TELOMEAN", "TELOSTD", "transMETmin", 
              "homeMETmin", "LeisMETmin_wk", "totMETmin", "RC2K", "STM2K", "UR", 
              "meanNDVI99_00_cbg", "meanNDVI01_02_cbg", "meanNDVI99_00_tract", 
              "meanNDVI01_02_tract", "NDI_tract", "NDI_blkgrp", "seg_index_tract", 
              "seg_index_blkgrp", "tractpm25_1999", "tractpm25_2000", "tractpm25_2001", 
              "tractpm25_2002", "bgpm25_1999", "bgpm25_2000", "bgpm25_2001", 
              "bgpm25_2002", "tract_redline", "tract_redline_pc", "blkgrp_redline", 
              "blkgrp_redline_pc")
# assign new column names
names(mdata) <- name_new

# ------------------------------------------------------------------------#
# label values in variables -----------------------------------------------
# ------------------------------------------------------------------------#

mdata$survey_period <- ifelse(mdata$survey_period == 1, "nhanes99_00", "nhanes01_02")

mdata$gender <- ifelse(mdata$gender == 1, "male", "female")

mdata$race_eth <- cut(mdata$race_eth, c(0, 1, 2, 3, 4, 5), c("mexican american", "other hispanic", "non-hispanic white", "non-hispanic black", "other_multi"), include.lowest=TRUE)

mdata$country_born <- cut(mdata$country_born, c(0, 1, 2, 3, 7, 9), c("usa", "mexico", "elsewhere", "refused", "don't know"), include.lowest=TRUE)

mdata$ed_level <- cut(mdata$ed_level, c(0, 1, 2, 3, 4, 5, 7, 9), c("less than 9th", "9-12, no diploma", "hs or ged", "some college, aa degree", "college grad or above", "refused", "don't know"), include.lowest=TRUE)

mdata$marital_status <- cut(mdata$marital_status, c(0, 1, 2, 3, 4, 5, 6, 77, 99), c("married", "widowed", "divorced", "separated", "never married", "living with partner", "refused", "don't know"), include.lowest=TRUE)

mdata$hh_income <- cut(mdata$hh_income, c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 77, 99), c("0-4999", "5000-9999", "10000-14999", "15000-19999", "20000-24999", "25000-34999", "35000-44999", "45000-54999", "55000-65999", "65000-74999", "75000 and over", "over 20000", "under 20000", "refused", "don't know"), include.lowest=TRUE)

# ------------------------------------------------------------------------#
# convert telomere length to basepairs ------------------------------------
# ------------------------------------------------------------------------#

mdata$basepairs <- 3274 + (2413 * mdata$TELOMEAN)

# ------------------------------------------------------------------------#
# determine drinkers and smokers ------------------------------------------
# ------------------------------------------------------------------------#

mdata$drinker <- ifelse(mdata$ALQ120Q == 0, "never", ifelse(!is.na(mdata$ALQ120Q), "current", NA))

# refused and don't know converted to NA for smoker
mdata$smoker <- ifelse(mdata$SMQ020 == 2, "never", ifelse(mdata$SMQ020 == 1, "current", NA))
mdata$smoker <- ifelse(mdata$smoker == "current" & mdata$SMQ040 == 3, "former", mdata$smoker)

# ------------------------------------------------------------------------#
# categorize physical activity --------------------------------------------
# ------------------------------------------------------------------------#

mdata$phyact <- ifelse(mdata$totMETmin >= 500, "meet_rec", ifelse(!is.na(mdata$totMETmin), "below_rec", NA))

# ------------------------------------------------------------------------#
# convert don't know and refuse to missing --------------------------------
# ------------------------------------------------------------------------#

mdata[mdata == "don't know"] <- NA
mdata[mdata == "refused"] <- NA

# ------------------------------------------------------------------------#
# summarise NDVI and PM2.5 over whole time frame --------------------------
# ------------------------------------------------------------------------#

mdata$NDVI_mean_tract <- rowMeans(mdata[, c("meanNDVI99_00_tract", "meanNDVI01_02_tract")], na.rm = T)
mdata$NDVI_mean_blkgrp <- rowMeans(mdata[, c("meanNDVI99_00_cbg", "meanNDVI01_02_cbg")], na.rm = T)
mdata$pm25_mean_tract <- rowMeans(mdata[, c("tractpm25_1999", "tractpm25_2000","tractpm25_2001","tractpm25_2002")], na.rm = T)
mdata$pm25_mean_blkgrp <- rowMeans(mdata[, c("bgpm25_1999","bgpm25_2000","bgpm25_2001","bgpm25_2002")], na.rm = T)

# ------------------------------------------------------------------------#
# make grouping variable
# ------------------------------------------------------------------------#

mdata$grouping <- paste(mdata$meanNDVI99_00_tract, mdata$meanNDVI01_02_tract , mdata$NDI_tract , mdata$seg_index_tract, mdata$tractpm25_1999, sep = "_")
mdata$grouping <- as.character(as.numeric(factor(mdata$grouping)))

# ------------------------------------------------------------------------#
# ------------------------------------------------------------------------#
# make analysis dataset
# ------------------------------------------------------------------------#
# ------------------------------------------------------------------------#

# filter to only those with telomere measures
telodata <- mdata[!is.na(mdata$TELOMEAN),]

# there are a few telomere measures that are extremely high, filter those out
telodata <- telodata[telodata$TELOMEAN < 5 ,]
# should be 7826 in this dataset

# there is one individual under 20, drop them, n should == 7825
telodata <- telodata[telodata$age_years >= 20, ]

# !!! be aware that there may be blanks in the dataset depending on how the initial data comes in.
# we can run the following and replace all blanks with NA
telodata[telodata == ""] <- NA

# save the analysis ready dataset -----------------------------------------

saveRDS(telodata, "../data/analysis_data.rds")

