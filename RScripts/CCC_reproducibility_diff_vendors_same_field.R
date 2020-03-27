set.seed(5) 
sep <- "," #choose separator

## Load doMC libraries
sysinf <- Sys.info()
if (!is.null(sysinf)){
  os <- sysinf['sysname']
  if (os == 'Darwin') {
    library(doMC)
    registerDoMC(cores=4)
  }
}

options(show.error.locations = TRUE) #debug
options(error=traceback)

## Select radiomic files from different vendors to compute features ICC
filename_ieo <- file.choose() # vendor 1 readfile
filename_ccu <- file.choose() # vendor 2 readfile
dataset_ieo <- read.csv(filename_ieo, header = TRUE, sep = sep)
dataset_ccu <- read.csv(filename_ccu, header = TRUE, sep = sep)
rf_ieo <- dataset_ieo[ ,-seq(1,24)] # remove settings columns (in vendor 1 there were 24 columns) 
rf_ccu <- dataset_ccu[ ,-seq(1,39)] # remove settings columns (in vendor 2 there were 39 columns) 
rf_ieo <- as.data.frame(lapply(rf_ieo, as.numeric))
rf_ccu <- as.data.frame(lapply(rf_ccu, as.numeric))

## import libraries
library(caret)
library(DescTools)
ieo <- as.matrix(rf_ieo)
ccu <- as.matrix(rf_ccu)

## calculate CCCs
correlation_coefficients_ccc_reproducibility <- sapply(seq.int(dim(ieo)[2]), function(i) DescTools::CCC(ieo[,i], ccu[,i])$rho.c$est)
correlation_coefficients_ccc_lwr_ci_reproducibility <- sapply(seq.int(dim(ieo)[2]), function(i) DescTools::CCC(ieo[,i], ccu[,i])$rho.c$lwr.ci)
correlation_coefficients_ccc_upr_ci_reproducibility <- sapply(seq.int(dim(ieo)[2]), function(i) DescTools::CCC(ieo[,i], ccu[,i])$rho.c$upr.ci)
names(correlation_coefficients_ccc_reproducibility) <- colnames(ieo)
names(correlation_coefficients_ccc_lwr_ci_reproducibility) <- colnames(ieo)
names(correlation_coefficients_ccc_upr_ci_reproducibility) <- colnames(ieo)

## CCC metrics
correlation_coefficients_ccc_excellent <- correlation_coefficients_ccc_reproducibility > 0.9
correlation_coefficients_ccc_excellent <- names(correlation_coefficients_ccc_excellent [correlation_coefficients_ccc_excellent  > 0.9])
table(correlation_coefficients_ccc_excellent)
correlation_coefficients_ccc_good <- correlation_coefficients_ccc_reproducibility > 0.75 & correlation_coefficients_ccc_reproducibility <= 0.9
table(correlation_coefficients_ccc_good)
correlation_coefficients_ccc_moderate <- correlation_coefficients_ccc_reproducibility >= 0.5 & correlation_coefficients_ccc_reproducibility <= 0.75
table(correlation_coefficients_ccc_moderate)
correlation_coefficients_ccc_poor <- correlation_coefficients_ccc_reproducibility < 0.5
table(correlation_coefficients_ccc_poor)

## Save metrics
n_excellent_ccc_reproducibility <- length(correlation_coefficients_ccc_reproducibility[correlation_coefficients_ccc_reproducibility > 0.9]) / ncol(rf_ieo)*100
n_good_ccc_reproducibility <- length(correlation_coefficients_ccc_reproducibility[correlation_coefficients_ccc_reproducibility <= 0.9 & correlation_coefficients_ccc_reproducibility > 0.75]) / ncol(rf_ieo)*100
n_moderate_ccc_reproducibility <- length(correlation_coefficients_ccc_reproducibility[correlation_coefficients_ccc_reproducibility <= 0.75 & correlation_coefficients_ccc_reproducibility > 0.5]) / ncol(rf_ieo)*100
n_poor_ccc_reproducibility <- length(correlation_coefficients_ccc_reproducibility[correlation_coefficients_ccc_reproducibility <= 0.5]) / ncol(rf_ieo)*100

scanner_15T_ccc_reproducibility <- c(n_excellent_ccc_reproducibility, n_good_ccc_reproducibility, n_moderate_ccc_reproducibility, n_poor_ccc_reproducibility)
scanner_15T_ccc_reproducibility_excellent <- correlation_coefficients_ccc_reproducibility[correlation_coefficients_ccc_reproducibility > 0.9]

## Filter by feature class type using CCC > 0.9
cutoffStabilityCCC <- 0.9
highlyCorrelated_ccc <- correlation_coefficients_ccc_reproducibility > cutoffStabilityCCC
cols <- which(highlyCorrelated_ccc)

colnames(rf_ieo)[cols]
features_high_ccc <- colnames(rf_ieo)[cols]

classType_high_ccc_feats <- list()
for (classType in c("original_", "log.sigma._", "wavelet.LH_", "wavelet.HL_", 
                    "wavelet.HH_", "wavelet.LL_", "square_", "squareroot_", "logarithm_",
                    "exponential_")) {
  classType_high_ccc_feats[[classType]] <- gsub(classType, "", features_high_ccc[ grep(classType,features_high_ccc)])
}

highly_ccc_features_by_class <- unique(unlist(classType_high_ccc_feats))

table(unlist(classType_high_ccc_feats))

## Filter image filter type using CCC > 0.9
cutoffStabilityCCC <- 0.9
highlyCorrelated_ccc <- correlation_coefficients_ccc_reproducibility > cutoffStabilityCCC
cols <- which(highlyCorrelated_ccc)

colnames(rf_ieo)[cols]
features_high_ccc <- colnames(rf_ieo)[cols]

filterType_high_ccc_feats <- list()
for (filterType in c('firstorder_Energy', 'firstorder_MeanAbsoluteDeviation', 'firstorder_RobustMeanAbsoluteDeviation', 
                     'firstorder_TotalEnergy', 'gldm_DependenceNonUniformity', 'glrlm_RunLengthNonUniformity', 
                     'glszm_SizeZoneNonUniformity', 'ngtdm_Coarseness')) {
  filterType_high_ccc_feats[[filterType]] <- gsub(filterType, "", features_high_ccc[ grep(filterType,features_high_ccc)])
}

highly_ccc_features_by_filter <- unique(unlist(filterType_high_ccc_feats))

table(unlist(filterType_high_ccc_feats))

