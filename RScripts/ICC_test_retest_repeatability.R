set.seed(5) 

sep <- "," # choose separator

## Load doMC libraries
sysinf <- Sys.info()
if (!is.null(sysinf)){
  os <- sysinf['sysname']
  if (os == 'Darwin') {
    library(doMC)
    registerDoMC(cores=4)
  }
}

options(show.error.locations = TRUE) # debug
options(error=traceback)

filename_test <- file.choose() # readfile 
filename_retest <- file.choose() # readfile 
dataset_test <- read.csv(filename_test, header = TRUE, sep = sep)
dataset_retest <- read.csv(filename_retest, header = TRUE, sep = sep)
dataset_features_test <- dataset_test[ ,-seq(1,24)] # first 25 cols do not contain features
dataset_features_retest <- dataset_retest[ ,-seq(1,24)] # first 25 cols do not contain features

## import libraries
library(psych)
library(caret)
library(irr)

## calculate ICCs
correlation_coefficients_agreement <- sapply(seq.int(dim(dataset_features_test)[2]), function(j) irr::icc(cbind(dataset_features_test[,j], dataset_features_retest[,j]) ,type = 'agreement')$value)
names(correlation_coefficients_agreement) <- colnames(dataset_features_test)

correlation_coefficients_agreement_excellent <- correlation_coefficients_agreement[correlation_coefficients_agreement > 0.9]
table(correlation_coefficients_agreement_excellent)
correlation_coefficients_agreement_good <- correlation_coefficients_agreement[correlation_coefficients_agreement > 0.75 & correlation_coefficients_agreement <= 0.9]
table(correlation_coefficients_agreement_good)
correlation_coefficients_agreement_moderate <- correlation_coefficients_agreement[correlation_coefficients_agreement >= 0.5 & correlation_coefficients_agreement <= 0.75]
table(correlation_coefficients_agreement_moderate)
correlation_coefficients_agreement_poor <- correlation_coefficients_agreement[correlation_coefficients_agreement < 0.5]
table(correlation_coefficients_agreement_poor)


n_excellent_TR_exp <- length(correlation_coefficients_agreement[correlation_coefficients_agreement > 0.9]) / ncol(dataset_features_test)*100
n_good_TR_exp <- length(correlation_coefficients_agreement[correlation_coefficients_agreement <= 0.9 & correlation_coefficients_agreement > 0.75]) / ncol(dataset_features_test)*100
n_moderate_TR_exp <- length(correlation_coefficients_agreement[correlation_coefficients_agreement <= 0.75 & correlation_coefficients_agreement > 0.5]) / ncol(dataset_features_test)*100
n_poor_TR_exp <- length(correlation_coefficients_agreement[correlation_coefficients_agreement <= 0.5]) / ncol(dataset_features_test)*100

########################################################################
#### Uncomment each block of 2 lines depending on the loaded files  ####
########################################################################

# ccu_15T_test_retest_no_repos <- c(n_excellent_TR_exp, n_good_TR_exp, n_moderate_TR_exp, n_poor_TR_exp)
# ccu_15T_test_retest_no_repos_agreement_excellent <- names(correlation_coefficients_agreement_excellent)

# ccu_15T_test_retest_repos <- c(n_excellent_TR_exp, n_good_TR_exp, n_moderate_TR_exp, n_poor_TR_exp)
# ccu_15T_test_retest_repos_agreement_excellent <- names(correlation_coefficients_agreement_excellent)

# ieo_15T_test_retest_no_repos <- c(n_excellent_TR_exp, n_good_TR_exp, n_moderate_TR_exp, n_poor_TR_exp)
# ieo_15T_test_retest_no_repos_agreement_excellent <- names(correlation_coefficients_agreement_excellent)

# ieo_15T_test_retest_repos <- c(n_excellent_TR_exp, n_good_TR_exp, n_moderate_TR_exp, n_poor_TR_exp)
# ieo_15T_test_retest_repos_agreement_excellent <- names(correlation_coefficients_agreement_excellent)

# ccu_3T_test_retest_no_repos <- c(n_excellent_TR_exp, n_good_TR_exp, n_moderate_TR_exp, n_poor_TR_exp)
# ccu_3T_test_retest_no_repos_agreement_excellent <- names(correlation_coefficients_agreement_excellent)

# ccu_3T_test_retest_repos <- c(n_excellent_TR_exp, n_good_TR_exp, n_moderate_TR_exp, n_poor_TR_exp)
# ccu_3T_test_retest_repos_agreement_excellent <- names(correlation_coefficients_agreement_excellent)
