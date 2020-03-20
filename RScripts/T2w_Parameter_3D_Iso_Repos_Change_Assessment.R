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

###############
# 3D Extraction
###############

# 3D iso 
filename_3D_Iso <- file.choose() # readfile
dataset_3D_Iso <- read.csv(filename_3D_Iso, header = TRUE, sep = sep)
dataset_3D_Iso_features <- dataset_3D_Iso[ ,-seq(1,24)] 

# 3D iso Repos
filename_3D_Iso_r <- file.choose() # readfile
dataset_3D_Iso_r <- read.csv(filename_3D_Iso_r, header = TRUE, sep = sep)
dataset_3D_Iso_r_features <- dataset_3D_Iso_r[ ,-seq(1,24)] 

## import libraries
# library(psych)
library(caret)
# library(irr)
library(DescTools)
library(stringr)

# correlation_coefficients_agreement_3D <- sapply(seq.int(dim(dataset_3D_Iso_r_features)[2]), function(i) psych::ICC(cbind(dataset_3D_Iso_features[,i],dataset_3D_Iso_r_features[,i]))$results$ICC[1])
correlation_coefficients_agreement_3D <- sapply(seq.int(dim(dataset_3D_Iso_r_features)[2]), function(i) DescTools::CCC(dataset_3D_Iso_features[,i],dataset_3D_Iso_r_features[,i])$rho.c$est)
names(correlation_coefficients_agreement_3D) <- colnames(dataset_3D_Iso_r_features)

length(correlation_coefficients_agreement_3D[correlation_coefficients_agreement_3D > 0.9])
length(correlation_coefficients_agreement_3D[correlation_coefficients_agreement_3D <= 0.9 & correlation_coefficients_agreement_3D > 0.75])
length(correlation_coefficients_agreement_3D[correlation_coefficients_agreement_3D <= 0.75 & correlation_coefficients_agreement_3D > 0.5])
length(correlation_coefficients_agreement_3D[correlation_coefficients_agreement_3D < 0.5])

n_excellent_agreement_reproducibility <- length(correlation_coefficients_agreement_3D[correlation_coefficients_agreement_3D > 0.9]) / ncol(dataset_3D_Iso_r_features)*100
n_good_agreement_reproducibility <- length(correlation_coefficients_agreement_3D[correlation_coefficients_agreement_3D <= 0.9 & correlation_coefficients_agreement_3D > 0.75]) / ncol(dataset_3D_Iso_r_features)*100
n_moderate_agreement_reproducibility <- length(correlation_coefficients_agreement_3D[correlation_coefficients_agreement_3D <= 0.75 & correlation_coefficients_agreement_3D > 0.5]) / ncol(dataset_3D_Iso_r_features)*100
n_poor_agreement_reproducibility <- length(correlation_coefficients_agreement_3D[correlation_coefficients_agreement_3D < 0.5]) / ncol(dataset_3D_Iso_r_features)*100

acq3D_agreement_repeatability <- c(n_excellent_agreement_reproducibility, n_good_agreement_reproducibility, n_moderate_agreement_reproducibility, n_poor_agreement_reproducibility)
acq3D_agreement_repeatability_excellent <- correlation_coefficients_agreement_reproducibility[correlation_coefficients_agreement_reproducibility > 0.9]


featureType_high_agreement_feats <- list()
for (filterType in c("original_", "log.sigma.1.mm.3D_", "wavelet.LLH_", "wavelet.LHL_", 
                     "wavelet.LHH_", "wavelet.HLL_", "wavelet.HLH_", "wavelet.HHL_", 
                     "wavelet.HHH_", "wavelet.LLL_", "square_", "squareroot_", 
                     "logarithm_", "exponential_")) {
  featureType_high_agreement_feats[[filterType]] <- gsub(filterType, "", 
                                                            names(correlation_coefficients_agreement_3D[correlation_coefficients_agreement_3D > 0.9])[ grep(filterType,names(correlation_coefficients_agreement_3D[correlation_coefficients_agreement_3D > 0.9]))])
}

high_agreement_features <- unique(unlist(featureType_high_agreement_feats))
table(unlist(featureType_high_agreement_feats))

###################
## Filter Type
###################
filterType_high_agreement_feats <- list()
for (filterType in high_agreement_features) {
  filterType_high_agreement_feats[[filterType]] <- gsub(filterType, "", 
                                                        names(correlation_coefficients_agreement_3D[correlation_coefficients_agreement_3D > 0.9])[ grep(filterType,names(correlation_coefficients_agreement_3D[correlation_coefficients_agreement_3D > 0.9]))])
}

# high_agreement_features <- unique(unlist(filterType_high_agreement_feats))

table(unlist(filterType_high_agreement_feats))

table_high_agreement_filters <- table(unlist(filterType_high_agreement_feats))


################
## Feature Class
################
featureClass_high_agreement_feats <- list()
for (featureClass in c("shape_", "firstorder_", "glcm_", "glrlm_", 
                     "glszm_", "ngtdm_", "gldm_")) {
  print(featureClass)
  print(str_count(list(names(correlation_coefficients_agreement_3D[correlation_coefficients_agreement_3D > 0.9])[ grep(featureClass,names(correlation_coefficients_agreement_3D[correlation_coefficients_agreement_3D > 0.9]))]), featureClass))
  featureClass_high_agreement_feats[[featureClass]] <- str_count(list(names(correlation_coefficients_agreement_3D[correlation_coefficients_agreement_3D > 0.9])[ grep(featureClass,names(correlation_coefficients_agreement_3D[correlation_coefficients_agreement_3D > 0.9]))]), featureClass)
}

featureClass_high_agreement_feats
