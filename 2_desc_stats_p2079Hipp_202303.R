# Date: 2022-06-17
# Update:2023-01-17
# S Ogletree
# Description: Script 2 for project - produce descriptive statistics

# analysis dataset from script 1
# [verify location of dataset]

adata <- readRDS("../data/analysis_data.rds")

# Descriptive statistics --------------------------------------------------

# must have the psych package installed on the machine!
summ1 <- round(psych::describe(adata), 3)
summ1$missing <- sapply(adata, function (x) {sum(is.na(x))})
summ1 <- summ1[!grepl("\\*", rownames(summ1)),]

# get variables for mean-sd
msd_vars <- c("TELOMEAN", "NDVI_mean_tract", "age_years", "family_pir", "bmi", "NDI_tract", "seg_index_tract","pm25_mean_tract")
summ1a <- summ1[rownames(summ1) %in% msd_vars,c(2:4, 14, 8, 9)]

capture.output(summ1a, file = '../disclosure/p2079_table_1_1.txt')

# get variables for n + %
nvars <- c("gender","race_eth","ed_level","country_born","marital_status","phyact","drinker","smoker", "tract_redline")

summ2 <- adata[, names(adata) %in% nvars]

# get tables for categorical variables
f1 <- function(x){table(x, useNA = "always")}

summ2a <- sapply(summ2, f1)
summ2b <- lapply(summ2a, function(x){round(x/7825, 3)})

# write tables out to text file
capture.output(summ2a, file = '../disclosure/p2079_table 1_n.txt')
capture.output(summ2b, file = '../disclosure/p2079_table 1_percent.txt')


# --- oops we didn't ask for correlations, maybe next time
# # correlations, pairwise
# 
# m <- psych::corr.test(adata[, unlist(lapply(adata, is.numeric))], use = "pairwise")
# write.csv(round(m$r, 3), '../disclosure/p2079_correlations_r-value.csv', na = "")
# # save the p-values for the correlations
# # unadjusted p below diagonal, adjust p above the diagonal
# write.csv(round(m$p, 5), '../disclosure/p2079_correlations_p-value.csv', na = "")
# 
# # telomere length by category, written to a text file
# mean.sd <- function(x) c(mean = mean(x, na.rm=T), sd = sd(x, na.rm=T), n = length(x))
# capture.output(tapply(adata$basepairs, adata$gender, mean.sd), file = "../disclosure/p2079_telomere_by_category.txt")
# capture.output(tapply(adata$basepairs, adata$race_eth, mean.sd), file = "../disclosure/p2079_telomere_by_category.txt", append = T)
# capture.output(tapply(adata$basepairs, adata$ed_level, mean.sd), file = "../disclosure/p2079_telomere_by_category.txt", append = T)
# capture.output(tapply(adata$basepairs, adata$marital_status, mean.sd), file = "../disclosure/p2079_telomere_by_category.txt", append = T)
# capture.output(tapply(adata$basepairs, adata$phyact, mean.sd), file = "../disclosure/p2079_telomere_by_category.txt", append = T)
# capture.output(tapply(adata$basepairs, adata$drinker, mean.sd), file = "../disclosure/p2079_telomere_by_category.txt", append = T)
# capture.output(tapply(adata$basepairs, adata$smoker, mean.sd), file = "../disclosure/p2079_telomere_by_category.txt", append = T)
