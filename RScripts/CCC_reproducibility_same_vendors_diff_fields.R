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
filename_15T <- file.choose() # readfile
filename_3T <- file.choose() # readfile

dataset_15T <- read.csv(filename_15T, header = TRUE, sep = sep)
dataset_3T <- read.csv(filename_3T, header = TRUE, sep = sep)

rf_15T <- dataset_15T[ ,-seq(1,24)] # remove settings columns 
rf_3T <- dataset_3T[ ,-seq(1,24)] # remove settings columns 
rf_15T <- as.data.frame(lapply(rf_15T, as.numeric))
rf_3T <- as.data.frame(lapply(rf_3T, as.numeric))
ccu_15T <- as.matrix(rf_15T)
ccu_3T <- as.matrix(rf_3T)

## import libraries
library(caret)
library(DescTools)

## calculate CCCs
correlation_coefficients_ccc_reproducibility <- sapply(seq.int(dim(ieo)[2]), function(i) DescTools::CCC(ccu_15T[,i], ccu_3T[,i])$rho.c$est)
correlation_coefficients_ccc_lwr_ci_reproducibility <- sapply(seq.int(dim(ieo)[2]), function(i) DescTools::CCC(ccu_15T[,i], ccu_3T[,i])$rho.c$lwr.ci)
correlation_coefficients_ccc_upr_ci_reproducibility <- sapply(seq.int(dim(ieo)[2]), function(i) DescTools::CCC(ccu_15T[,i], ccu_3T[,i])$rho.c$upr.ci)
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

scanner_15Tvs3T_ccc_reproducibility <- c(n_excellent_ccc_reproducibility, n_good_ccc_reproducibility, n_moderate_ccc_reproducibility, n_poor_ccc_reproducibility)
scanner_15Tvs3T_ccc_reproducibility_excellent <- correlation_coefficients_ccc_reproducibility[correlation_coefficients_ccc_reproducibility > 0.9]

## Filter by feature class type using CCC > 0.9
cutoffStabilityCCC <- 0.9
highlyCorrelated_ccc <- correlation_coefficients_ccc_reproducibility > cutoffStabilityCCC
cols <- which(highlyCorrelated_ccc)

colnames(rf_15T)[cols]
features_high_ccc <- colnames(rf_15T)[cols]

classType_high_ccc_feats <- list()
for (classType in c("original_", "log.sigma._", "wavelet.LH_", "wavelet.HL_", 
                    "wavelet.HH_", "wavelet.LL_", "square_", "squareroot_", "logarithm_",
                    "exponential_")) {
  classType_high_ccc_feats[[classType]] <- gsub(classType, "", features_high_ccc[ grep(classType,features_high_ccc)])
}

highly_ccc_features_by_class <- unique(unlist(classType_high_ccc_feats))

table(unlist(classType_high_ccc_feats))

## Filter image filter type using CCC > 0.9
highlyCorrelated_ccc <- correlation_coefficients_ccc_reproducibility > cutoffStabilityCCC
cols <- which(highlyCorrelated_ccc)

colnames(rf_15T)[cols]
features_high_ccc <- colnames(rf_15T)[cols]

filterType_high_ccc_feats <- list()
for (filterType in c('firstorder_Energy', 'firstorder_Entropy', 'firstorder_MeanAbsoluteDeviation', 
                     'firstorder_RobustMeanAbsoluteDeviation', 'firstorder_Skewness', 'firstorder_TotalEnergy', 
                     'firstorder_Uniformity', 'glcm_Correlation', 'glcm_DifferenceAverage', 'glcm_DifferenceEntropy', 
                     'glcm_Imc1', 'glcm_Imc2', 'glcm_Idm', 'glcm_Id', 'glcm_InverseVariance', 'glcm_JointEnergy', 
                     'glcm_JointEntropy', 'glcm_MaximumProbability', 'glcm_SumEntropy', 'gldm_DependenceEntropy', 
                     'gldm_DependenceNonUniformityNormalized', 'gldm_DependenceNonUniformity', 'gldm_GrayLevelNonUniformity', 
                     'gldm_LargeDependenceEmphasis', 'gldm_LargeDependenceLowGrayLevelEmphasis', 
                     'gldm_SmallDependenceEmphasis', 'gldm_SmallDependenceLowGrayLevelEmphasis', 
                     'glrlm_GrayLevelNonUniformityNormalized', 'glrlm_GrayLevelNonUniformity', 
                     'glrlm_LowGrayLevelRunEmphasis', 'glrlm_RunEntropy', 'glrlm_RunLengthNonUniformity', 
                     'glrlm_RunPercentage', 'glrlm_ShortRunLowGrayLevelEmphasis', 'glszm_GrayLevelNonUniformityNormalized', 
                     'glszm_GrayLevelNonUniformity', 'glszm_LowGrayLevelZoneEmphasis', 
                     'glszm_SizeZoneNonUniformityNormalized', 'glszm_SizeZoneNonUniformity', 
                     'glszm_SmallAreaLowGrayLevelEmphasis', 'glszm_ZoneEntropy', 'ngtdm_Coarseness', 'ngtdm_Strength')) {
  filterType_high_ccc_feats[[filterType]] <- gsub(filterType, "", features_high_ccc[ grep(filterType,features_high_ccc)])
}

highly_ccc_features_by_filter <- unique(unlist(filterType_high_ccc_feats))

table(unlist(filterType_high_ccc_feats))


