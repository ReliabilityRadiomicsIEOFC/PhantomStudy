set.seed(5)

data_dir <- "/Users/joaosantinha/Documents/FChampalimaud/Data/Phantom/PhantomStudy/RadiomicFiles/CCU1_5T_TE_TR/" # put folder path here TE and TR Experimentradiomic output files are located

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

library(agRee)

options(show.error.locations = TRUE) # debug
options(error=traceback)
setwd(data_dir) # set the working director
# TE 80
filename_80 <- paste(data_dir, "T2wAX_TE80_clean.csv", sep = "") #readfile
dataset_80 <- read.csv(filename_80, header = TRUE, sep = sep)
dataset_80_features <- dataset_80[ ,-seq(1,24)] 
# TE 85
filename_85 <- paste(data_dir, "T2wAX_TE85_clean.csv", sep = "") #readfile
dataset_85 <- read.csv(filename_85, header = TRUE, sep = sep)
dataset_85_features <- dataset_85[ ,-seq(1,24)] 
# TE 90
filename_90 <- paste(data_dir, "T2wAX_TE90_clean.csv", sep = "") #readfile
dataset_90 <- read.csv(filename_90, header = TRUE, sep = sep)
dataset_90_features <- dataset_90[ ,-seq(1,24)] 
# TE 95
filename_95 <- paste(data_dir, "T2wAX_TE95_clean.csv", sep = "") #readfile
dataset_95 <- read.csv(filename_95, header = TRUE, sep = sep)
dataset_95_features <- dataset_95[ ,-seq(1,24)] 
# TE 100
filename_100 <- paste(data_dir, "T2wAX_TE100_clean.csv", sep = "") #readfile
dataset_100 <- read.csv(filename_100, header = TRUE, sep = sep)
dataset_100_features <- dataset_100[ ,-seq(1,24)] 
# TE 105
filename_105 <- paste(data_dir, "T2wAX_TE105_clean.csv", sep = "") #readfile
dataset_105 <- read.csv(filename_105, header = TRUE, sep = sep)
dataset_105_features <- dataset_105[ ,-seq(1,24)] 
# TE 110
filename_110 <- paste(data_dir, "T2wAX_TE110_clean.csv", sep = "") #readfile
dataset_110 <- read.csv(filename_110, header = TRUE, sep = sep)
dataset_110_features <- dataset_110[ ,-seq(1,24)] 
# TE 115
filename_115 <- paste(data_dir, "T2wAX_TE115_clean.csv", sep = "") #readfile
dataset_115 <- read.csv(filename_115, header = TRUE, sep = sep)
dataset_115_features <- dataset_115[ ,-seq(1,24)] 
# TE 120
filename_120 <- paste(data_dir, "T2wAX_TE120_clean.csv", sep = "") #readfile
dataset_120 <- read.csv(filename_120, header = TRUE, sep = sep)
dataset_120_features <- dataset_120[ ,-seq(1,24)] 

rf_acq <- list()
rf_acq[1] <- list(dataset_80_features)
rf_acq[2] <- list(dataset_85_features)
rf_acq[3] <- list(dataset_90_features)
rf_acq[4] <- list(dataset_95_features)
rf_acq[5] <- list(dataset_100_features)
rf_acq[6] <- list(dataset_105_features)
rf_acq[7] <- list(dataset_110_features)
rf_acq[8] <- list(dataset_115_features)
rf_acq[9] <- list(dataset_120_features)

acq <-list()
for (i in 1:9) {
  acq[i] <- rf_acq[i]
}

wCV_TEs_list <- list()
for (k in 1:8) {
  for (m in (k+1):9) {
    wCV_TEs_list[(k-1)*9+m] <- list(sapply(seq.int(dim(acq[[k]])[2]), 
                                                                 function(i) abs(100*agRee::agree.wscv(cbind(unlist(acq[[k]][i]), unlist(acq[[m]][i])))$value)))
    names(wCV_TEs_list[[(k-1)*9+m]]) <- colnames(rf_acq[[1]])
  }
}

library(tidyverse)
filtered_wCV_TE_Exp_list <- wCV_TEs_list %>% discard(is.null)

##########################################
## Experiment with fixed TE and changing TR
# TE 100
filename_100_r <- paste(data_dir, "T2wAX_TE100 Range_clean.csv", sep = "") #readfile
dataset_100_r <- read.csv(filename_100_r, header = TRUE, sep = sep)
dataset_100_r_features <- dataset_100_r[ ,-seq(1,24)] 

wCVs_TRs <- sapply(seq.int(dim(dataset_100_r_features)[2]), 
                                                function(i) abs(100*agRee::agree.wscv(cbind(dataset_100_features[,i], dataset_100_r_features[,i]))$value))
names(wCVs_TRs) <- colnames(dataset_100_r_features)
