# 2022-07-19
# Update: 2023-02-17

# Check if needed packages are installed

# see what R version is running
version

# list all installed packages
pkglist <- installed.packages(fields = "Package")


# packages needed (with their dependencies) - lme4, psych, haven
# can just try to load them as a first test
library(lme4) # for multilevel models
library(psych) # for descriptive statistics
library(haven) # for reading in Stata/SAS files if necessary

# packages and their dependencies, for troubleshooting
pkgneed <- c("haven", "lme4", "psych", "graphics", "grid", "splines", "utils", "parallel", "MASS", "lattice", "boot","nlme", "minqa" , "nloptr", "Rcpp", "RcppEigen", "mnormt","stats","grDevices","methods","lattice", "cli", "forcats" , "hms", "lifecycle", "readr", "rlang", "tibble", "tidyselect", "vctrs")

pkgneed %in% pkglist
