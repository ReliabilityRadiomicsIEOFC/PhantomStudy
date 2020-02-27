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
library(psych)
library(caret)
library(irr)
ieo <- as.matrix(rf_ieo)
ccu <- as.matrix(rf_ccu)

## calculate ICCs
correlation_coefficients_consistency_reproducibility <- sapply(seq.int(dim(ieo)[2]), function(i) psych::ICC(matrix(cbind(ieo[,i], ccu[,i]),nrow=length(ieo[,i])))$results$ICC[3])
correlation_coefficients_agreement_reproducibility <- sapply(seq.int(dim(ieo)[2]), function(i) psych::ICC(matrix(cbind(ieo[,i], ccu[,i]),nrow=length(ieo[,i])))$results$ICC[1])
names(correlation_coefficients_consistency_reproducibility) <- colnames(ieo)
names(correlation_coefficients_agreement_reproducibility) <- colnames(ieo)

## absolute agreement metrics
correlation_coefficients_agreement_excellent_reproducibility_seq_frfse <- correlation_coefficients_agreement_reproducibility > 0.9
correlation_coefficients_agreement_excellent_reproducibility_seq_frfse_names <- names(correlation_coefficients_agreement_excellent_reproducibility_seq_frfse [correlation_coefficients_agreement_excellent_reproducibility_seq_frfse  > 0.9])
table(correlation_coefficients_agreement_excellent)
correlation_coefficients_agreement_good <- correlation_coefficients_agreement_reproducibility > 0.75 & correlation_coefficients_agreement_reproducibility <= 0.9
table(correlation_coefficients_agreement_good)
correlation_coefficients_agreement_moderate <- correlation_coefficients_agreement_reproducibility >= 0.5 & correlation_coefficients_agreement_reproducibility <= 0.75
table(correlation_coefficients_agreement_moderate)
correlation_coefficients_agreement_poor <- correlation_coefficients_agreement_reproducibility < 0.5
table(correlation_coefficients_agreement_poor)

## consistency agreement metrics
correlation_coefficients_consistency_excellent <- correlation_coefficients_consistency_reproducibility > 0.9
table(correlation_coefficients_consistency_excellent)
correlation_coefficients_consistency_good <- correlation_coefficients_consistency_reproducibility > 0.75 & correlation_coefficients_consistency_reproducibility <= 0.9
table(correlation_coefficients_consistency_good)
correlation_coefficients_consistency_moderate <- correlation_coefficients_consistency_reproducibility >= 0.5 & correlation_coefficients_consistency_reproducibility <= 0.75
table(correlation_coefficients_consistency_moderate)
correlation_coefficients_consistency_poor <- correlation_coefficients_consistency_reproducibility < 0.5
table(correlation_coefficients_consistency_poor)


## Save metrics
n_excellent_agreement_reproducibility <- length(correlation_coefficients_agreement_reproducibility[correlation_coefficients_agreement_reproducibility > 0.9]) / ncol(rf_ieo)*100
n_good_agreement_reproducibility <- length(correlation_coefficients_agreement_reproducibility[correlation_coefficients_agreement_reproducibility <= 0.9 & correlation_coefficients_agreement_reproducibility > 0.75]) / ncol(rf_ieo)*100
n_moderate_agreement_reproducibility <- length(correlation_coefficients_agreement_reproducibility[correlation_coefficients_agreement_reproducibility <= 0.75 & correlation_coefficients_agreement_reproducibility > 0.5]) / ncol(rf_ieo)*100
n_poor_agreement_reproducibility <- length(correlation_coefficients_agreement_reproducibility[correlation_coefficients_agreement_reproducibility <= 0.5]) / ncol(rf_ieo)*100

n_excellent_consistency_reproducibility <- length(correlation_coefficients_consistency_reproducibility[correlation_coefficients_consistency_reproducibility > 0.9]) / ncol(rf_ieo)*100
n_good_consistency_reproducibility <- length(correlation_coefficients_consistency_reproducibility[correlation_coefficients_consistency_reproducibility <= 0.9 & correlation_coefficients_consistency_reproducibility > 0.75]) / ncol(rf_ieo)*100
n_moderate_consistency_reproducibility <- length(correlation_coefficients_consistency_reproducibility[correlation_coefficients_consistency_reproducibility <= 0.75 & correlation_coefficients_consistency_reproducibility > 0.5]) / ncol(rf_ieo)*100
n_poor_consistency_reproducibility <- length(correlation_coefficients_consistency_reproducibility[correlation_coefficients_consistency_reproducibility <= 0.5]) / ncol(rf_ieo)*100

scanner_15T_agreement_reproducibility <- c(n_excellent_agreement_reproducibility, n_good_agreement_reproducibility, n_moderate_agreement_reproducibility, n_poor_agreement_reproducibility)
scanner_15T_agreement_reproducibility_excellent <- correlation_coefficients_agreement_reproducibility[correlation_coefficients_agreement_reproducibility > 0.9]

scanner_15T_consistency_reproducibility <- c(n_excellent_consistency_reproducibility, n_good_consistency_reproducibility, n_moderate_consistency_reproducibility, n_poor_consistency_reproducibility)
scanner_15T_consistency_reproducibility_excellent <- correlation_coefficients_consistency_reproducibility[correlation_coefficients_consistency_reproducibility > 0.9]


## Filter by feature class type using consistency
cutoffStabilityICC <- 0.9
highlyCorrelated_consistency <- correlation_coefficients_consistency_reproducibility > cutoffStabilityICC
cols <- which(highlyCorrelated_consistency)

colnames(rf_ieo)[cols]
features_high_consistency <- colnames(rf_ieo)[cols]

classType_high_consistency_feats <- list()
for (classType in c("original_", "log.sigma._", "wavelet.LH_", "wavelet.HL_", 
                    "wavelet.HH_", "wavelet.LL_", "square_", "squareroot_", "logarithm_",
                    "exponential_")) {
  classType_high_consistency_feats[[classType]] <- gsub(classType, "", features_high_consistency[ grep(classType,features_high_consistency)])
}

highly_consistency_features_by_class <- unique(unlist(classType_high_consistency_feats))

table(unlist(classType_high_consistency_feats))

## Filter by feature class type using absolute agreement
highlyCorrelated_agreement <- correlation_coefficients_agreement_reproducibility > cutoffStabilityICC
cols <- which(highlyCorrelated_agreement)

colnames(rf_ieo)[cols]
features_high_agreement <- colnames(rf_ieo)[cols]

classType_high_agreement_feats <- list()
for (classType in c("original_", "log.sigma._", "wavelet.LH_", "wavelet.HL_", 
                    "wavelet.HH_", "wavelet.LL_", "square_", "squareroot_", "logarithm_",
                    "exponential_")) {
  classType_high_agreement_feats[[classType]] <- gsub(classType, "", features_high_agreement[ grep(classType,features_high_agreement)])
}

highly_agreement_features_by_class <- unique(unlist(classType_high_agreement_feats))

table(unlist(classType_high_agreement_feats))

## Filter by feature class using consistency
cutoffStabilityICC <- 0.9
highlyCorrelated_consistency <- correlation_coefficients_consistency_reproducibility > cutoffStabilityICC
cols <- which(highlyCorrelated_consistency)

colnames(rf_ieo)[cols]
features_high_consistency <- colnames(rf_ieo)[cols]

filterType_high_consistency_feats <- list()
for (filterType in c('firstorder_Energy', 'firstorder_TotalEnergy', 'gldm_DependenceNonUniformity',
                     'gldm_GrayLevelNonUniformity', 'glrlm_GrayLevelNonUniformity',
                     'glrlm_RunLengthNonUniformity', 'glszm_GrayLevelNonUniformity',
                     'glszm_SizeZoneNonUniformity', 'ngtdm_Busyness', 'ngtdm_Coarseness',
                     'ngtdm_Strength')) {
  filterType_high_consistency_feats[[filterType]] <- gsub(filterType, "", features_high_consistency[ grep(filterType,features_high_consistency)])
}

highly_consistency_features_by_filter <- unique(unlist(filterType_high_consistency_feats))

table(unlist(filterType_high_consistency_feats))

## Filter image filter type using absolute agreement
highlyCorrelated_agreement <- correlation_coefficients_agreement_reproducibility > cutoffStabilityICC
cols <- which(highlyCorrelated_agreement)

colnames(rf_ieo)[cols]
features_high_agreement <- colnames(rf_ieo)[cols]

filterType_high_agreement_feats <- list()
for (filterType in c('firstorder_Energy', 'firstorder_TotalEnergy', 'gldm_DependenceNonUniformity',
                     'gldm_GrayLevelNonUniformity', 'glrlm_GrayLevelNonUniformity',
                     'glrlm_RunLengthNonUniformity', 'glszm_GrayLevelNonUniformity',
                     'glszm_SizeZoneNonUniformity', 'ngtdm_Busyness', 'ngtdm_Coarseness',
                     'ngtdm_Strength')) {
  filterType_high_agreement_feats[[filterType]] <- gsub(filterType, "", features_high_agreement[ grep(filterType,features_high_agreement)])
}

highly_agreement_features_by_filter <- unique(unlist(filterType_high_agreement_feats))

table(unlist(filterType_high_agreement_feats))

