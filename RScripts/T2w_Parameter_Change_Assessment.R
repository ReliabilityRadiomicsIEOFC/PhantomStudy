set.seed(5)

data_dir <- "<Folder_with_TE_TR_Experiment_Radiomics_Output>/" # put folder path here TE and TR Experimentradiomic output files are located

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

library(DescTools)

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

correlation_coefficients_agreement <- list()
for (k in 1:8) {
  for (m in (k+1):9) {
    correlation_coefficients_agreement[(k-1)*9+m] <- list(sapply(seq.int(dim(acq[[k]])[2]), function(i) DescTools::CCC(unlist(acq[[k]][i]), unlist(acq[[m]][i]))$rho.c$est))
    names(correlation_coefficients_agreement[[(k-1)*9+m]]) <- colnames(rf_acq[[1]])
    if (k == 1 && m == 2) {
      correlation_coefficients_agreement_excellent_acq1 <- names(correlation_coefficients_agreement[[(k-1)*9+m]][correlation_coefficients_agreement[[(k-1)*9+m]] > 0.9])
    } else{
      correlation_coefficients_agreement_excellent_acq1 <- intersect(correlation_coefficients_agreement_excellent_acq1, names(correlation_coefficients_agreement[[(k-1)*9+m]][correlation_coefficients_agreement[[(k-1)*9+m]] > 0.9]))
    }
  }
}

# par(mfrow=c(2,2))
dev.off()
library(corrplot)
# Excellent aggreement
correlation_matrix_agreement_exc <- matrix(nrow = 9, ncol = 9, dimnames = list(c('80', '85', '90', '95', '100', '105', '110', '115', '120'),c('80', '85', '90', '95', '100', '105', '110', '115', '120'))) 
for (k in 1:8) {
  for (m in (k+1):9) {
    correlation_matrix_agreement_exc[k,m] <- length(names(correlation_coefficients_agreement[[(k-1)*9+m]][correlation_coefficients_agreement[[(k-1)*9+m]] > 0.9]))
  } 
}

correlation_matrix_agreement_exc[is.na(correlation_matrix_agreement_exc)] <- 0

corrplot(correlation_matrix_agreement_exc/ncol(dataset_80_features)*100, method = "number", type = "upper", tl.pos = "d", 
         is.corr = FALSE, title="A (%)", mar=c(0,0,2,0), cl.lim = c(-10, 100))

# Good aggreement
correlation_matrix_agreement_good <- matrix(nrow = 9, ncol = 9, dimnames = list(c('80', '85', '90', '95', '100', '105', '110', '115', '120'),c('80', '85', '90', '95', '100', '105', '110', '115', '120'))) 
for (k in 1:8) {
  for (m in (k+1):9) {
    correlation_matrix_agreement_good[k,m] <- length(names(correlation_coefficients_agreement[[(k-1)*9+m]][correlation_coefficients_agreement[[(k-1)*9+m]] <= 0.9 & correlation_coefficients_agreement[[(k-1)*9+m]] > 0.75]))
  } 
}

correlation_matrix_agreement_good[is.na(correlation_matrix_agreement_good)] <- 0

corrplot(correlation_matrix_agreement_good/ncol(dataset_80_features)*100, method = "number", type = "upper", tl.pos = "d", 
         is.corr = FALSE, title="B (%)", mar=c(0,0,2,0), cl.lim = c(-10, 40))

# Moderate aggreement
correlation_matrix_agreement_mod <- matrix(nrow = 9, ncol = 9, dimnames = list(c('80', '85', '90', '95', '100', '105', '110', '115', '120'),c('80', '85', '90', '95', '100', '105', '110', '115', '120'))) 
for (k in 1:8) {
  for (m in (k+1):9) {
    correlation_matrix_agreement_mod[k,m] <- length(names(correlation_coefficients_agreement[[(k-1)*9+m]][correlation_coefficients_agreement[[(k-1)*9+m]] <= 0.75 & correlation_coefficients_agreement[[(k-1)*9+m]] > 0.5]))
  } 
}

correlation_matrix_agreement_mod[is.na(correlation_matrix_agreement_mod)] <- 0

corrplot(correlation_matrix_agreement_mod/ncol(dataset_80_features)*100, method = "number", type = "upper", tl.pos = "d", 
         is.corr = FALSE, title="C (%)", mar=c(0,0,2,0), cl.lim = c(-10, 38))

# Poor aggreement
correlation_matrix_agreement_poor <- matrix(nrow = 9, ncol = 9, dimnames = list(c('80', '85', '90', '95', '100', '105', '110', '115', '120'),c('80', '85', '90', '95', '100', '105', '110', '115', '120'))) 
for (k in 1:8) {
  for (m in (k+1):9) {
    correlation_matrix_agreement_poor[k,m] <- length(names(correlation_coefficients_agreement[[(k-1)*9+m]][correlation_coefficients_agreement[[(k-1)*9+m]] <= 0.5]))
  } 
}

correlation_matrix_agreement_poor[is.na(correlation_matrix_agreement_poor)] <- 0

corrplot(correlation_matrix_agreement_poor/ncol(dataset_80_features)*100, method = "number", type = "upper", tl.pos = "d", 
         is.corr = FALSE, title="D (%)", mar=c(0,0,2,0), cl.lim = c(-10, 18))

##########################################
## Experiment with fixed TE and changing TR
# TE 100
filename_100_r <- paste(data_dir, "T2wAX_TE100 Range_clean.csv", sep = "") #readfile
dataset_100_r <- read.csv(filename_100_r, header = TRUE, sep = sep)
dataset_100_r_features <- dataset_100_r[ ,-seq(1,24)] 

correlation_coefficients_agreement_TR <- sapply(seq.int(dim(dataset_100_r_features)[2]), function(i) DescTools::CCC(dataset_100_features[,i], dataset_100_r_features[,i])$rho.c$est)
names(correlation_coefficients_agreement_TR) <- colnames(dataset_100_r_features)

length(correlation_coefficients_agreement_TR[correlation_coefficients_agreement_TR > 0.9])
length(correlation_coefficients_agreement_TR[correlation_coefficients_agreement_TR <= 0.9 & correlation_coefficients_agreement_TR > 0.75])
length(correlation_coefficients_agreement_TR[correlation_coefficients_agreement_TR <= 0.75 & correlation_coefficients_agreement_TR > 0.5])
length(correlation_coefficients_agreement_TR[correlation_coefficients_agreement_TR <= 0.5])

n_excellent_TR_exp <- length(correlation_coefficients_agreement_TR[correlation_coefficients_agreement_TR > 0.9]) / ncol(dataset_100_r_features)*100
n_good_TR_exp <- length(correlation_coefficients_agreement_TR[correlation_coefficients_agreement_TR <= 0.9 & correlation_coefficients_agreement_TR > 0.75]) / ncol(dataset_100_r_features)*100
n_moderate_TR_exp <- length(correlation_coefficients_agreement_TR[correlation_coefficients_agreement_TR <= 0.75 & correlation_coefficients_agreement_TR > 0.5]) / ncol(dataset_100_r_features)*100
n_poor_TR_exp <- length(correlation_coefficients_agreement_TR[correlation_coefficients_agreement_TR <= 0.5]) / ncol(dataset_100_r_features)*100

specie <- c("TR" )
condition <- rep(c("Excellent" , "Good" , "Moderate", "Poor") , 1)
value <- c(n_excellent_TR_exp, n_good_TR_exp, n_moderate_TR_exp, n_poor_TR_exp)
data <- data.frame(specie,condition,value)
library(ggplot2)
library(viridis)

ggplot(data, aes(fill=condition, y=value, x=specie)) + 
  geom_bar(position="stack", stat="identity") # + scale_fill_viridis(discrete = T)

