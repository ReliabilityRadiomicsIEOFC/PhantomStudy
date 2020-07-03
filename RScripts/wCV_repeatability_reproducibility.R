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

filename_acq1 <- file.choose() # readfile 
filename_acq2 <- file.choose() # readfile 
dataset_acq1 <- read.csv(filename_acq1, header = TRUE, sep = sep)
dataset_acq2 <- read.csv(filename_acq2, header = TRUE, sep = sep)
dataset_features_acq1 <- dataset_acq1[ ,-seq(1,24)] # first 25 cols do not contain features
dataset_features_acq2 <- dataset_acq2[ ,-seq(1,24)] # first 25 cols do not contain features

## import libraries
library(rjags) # mac requires running: brew install jags #(on terminal) 
library(agRee)

## calculate overall wCV
## abs is because some metrics have negative mean which would make wCV negative which doesn't make sense
wsCVmean <- sapply(seq.int(dim(dataset_features_acq1)[2]), function(j) abs(100*agRee::agree.wscv(cbind(dataset_features_acq1[,j], dataset_features_acq2[,j]))$value))
names(wsCVmean) <- colnames(dataset_features_acq1)
wsCVmean <- data.frame(wCV = wsCVmean)

library(RColorBrewer)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n = 2
cols = gg_color_hue(n)

p_wcv_hist <- ggplot(wsCVmean, aes(x=wCV)) + geom_histogram(binwidth=0.5, color=cols[1], fill=cols[1], alpha=0.75) + labs(x = "wCV (%)") + xlim(0, 20)
print(p_wcv_hist)

n_abs_leq_1per <- length(wsCVmean[wsCVmean <= 1])
n_abs_leq_5per <- length(wsCVmean[wsCVmean > 1 & wsCVmean <= 5 ])
n_abs_leq_10per <- length(wsCVmean[wsCVmean > 5 & wsCVmean <= 10 ])
n_abs_leq_50per <- length(wsCVmean[wsCVmean > 10 & wsCVmean <= 50  ])
n_abs_geq_50per <- length(wsCVmean[wsCVmean > 50])

n_rel_leq_1per <- length(wsCVmean[wsCVmean <= 1]) / nrow(wsCVmean) * 100
n_rel_leq_5per <- length(wsCVmean[wsCVmean > 1 & wsCVmean <= 5 ]) / nrow(wsCVmean) * 100
n_rel_leq_10per <- length(wsCVmean[wsCVmean > 5 & wsCVmean <= 10 ]) / nrow(wsCVmean) * 100
n_rel_leq_50per <- length(wsCVmean[wsCVmean > 10 & wsCVmean <= 50  ]) / nrow(wsCVmean) * 100
n_rel_geq_50per <- length(wsCVmean[wsCVmean > 50]) / nrow(wsCVmean) * 100

n_absolute_overall <- c(n_abs_leq_1per, n_abs_leq_5per, n_abs_leq_10per, n_abs_leq_50per, n_abs_geq_50per)
n_relative_overall <- c(n_rel_leq_1per, n_rel_leq_5per, n_rel_leq_10per, n_rel_leq_50per, n_rel_geq_50per)

## calculate sizes wCV (1 - big ... 4 - small)
n <- 4
indeces_size_1 <- c(1, 5, 9, 13)
indeces_size_2 <- c(2, 6, 10, 14)
indeces_size_3 <- c(3, 7, 11, 15)
indeces_size_4 <- c(4, 8, 12, 16)
wsCVmean_size_1 <- sapply(seq.int(dim(dataset_features_acq1)[2]), function(j) abs(100*agRee::agree.wscv(cbind(dataset_features_acq1[indeces_size_1,j], dataset_features_acq2[indeces_size_1,j]))$value))
names(wsCVmean_size_1) <- colnames(dataset_features_acq1)

wsCVmean_size_2 <- sapply(seq.int(dim(dataset_features_acq1)[2]), function(j) abs(100*agRee::agree.wscv(cbind(dataset_features_acq1[indeces_size_2,j], dataset_features_acq2[indeces_size_2,j]))$value))
names(wsCVmean_size_2) <- colnames(dataset_features_acq1)

wsCVmean_size_3 <- sapply(seq.int(dim(dataset_features_acq1)[2]), function(j) abs(100*agRee::agree.wscv(cbind(dataset_features_acq1[indeces_size_3,j], dataset_features_acq2[indeces_size_3,j]))$value))
names(wsCVmean_size_3) <- colnames(dataset_features_acq1)

wsCVmean_size_4 <- sapply(seq.int(dim(dataset_features_acq1)[2]), function(j) abs(100*agRee::agree.wscv(cbind(dataset_features_acq1[indeces_size_4,j], dataset_features_acq2[indeces_size_4,j]))$value))
names(wsCVmean_size_4) <- colnames(dataset_features_acq1)

wsCVmean_sizes <- data.frame(wCV.size1 = wsCVmean_size_1, wCV.size2 = wsCVmean_size_2, wCV.size3 = wsCVmean_size_3, wCV.size4 = wsCVmean_size_4)

wsCVmean_sizes_hist_df <- data.frame(wCV=c(wsCVmean_size_1, wsCVmean_size_2, wsCVmean_size_3, wsCVmean_size_4), size=c(rep(1, length(wsCVmean_size_1)), rep(2, length(wsCVmean_size_2)), rep(3, length(wsCVmean_size_3)), rep(4, length(wsCVmean_size_4))))
wsCVmean_sizes_hist_df$size <- as.factor(wsCVmean_sizes_hist_df$size) 
p_sizes_hist <- ggplot(wsCVmean_sizes_hist_df, aes(x=wCV, color=size, fill=size)) + 
  geom_histogram(binwidth=0.5, position="dodge", alpha=0.25) + labs(x = "wCV (%)") + xlim(0, 20)
print(p_sizes_hist)

# size 1 summary
n_abs_leq_1per_size_1 <- length(wsCVmean_size_1[wsCVmean_size_1 <= 1])
n_abs_leq_5per_size_1 <- length(wsCVmean_size_1[wsCVmean_size_1 > 1 & wsCVmean_size_1 <= 5 ])
n_abs_leq_10per_size_1 <- length(wsCVmean_size_1[wsCVmean_size_1 > 5 & wsCVmean_size_1 <= 10 ])
n_abs_leq_50per_size_1 <- length(wsCVmean_size_1[wsCVmean_size_1 > 10 & wsCVmean_size_1 <= 50  ])
n_abs_geq_50per_size_1 <- length(wsCVmean_size_1[wsCVmean_size_1 > 50])

n_rel_leq_1per_size_1 <- length(wsCVmean_size_1[wsCVmean_size_1 <= 1]) / length(wsCVmean_size_1) * 100
n_rel_leq_5per_size_1 <- length(wsCVmean_size_1[wsCVmean_size_1 > 1 & wsCVmean_size_1 <= 5 ]) / length(wsCVmean_size_1) * 100
n_rel_leq_10per_size_1 <- length(wsCVmean_size_1[wsCVmean_size_1 > 5 & wsCVmean_size_1 <= 10 ]) / length(wsCVmean_size_1) * 100
n_rel_leq_50per_size_1 <- length(wsCVmean_size_1[wsCVmean_size_1 > 10 & wsCVmean_size_1 <= 50  ]) / length(wsCVmean_size_1) * 100
n_rel_geq_50per_size_1 <- length(wsCVmean_size_1[wsCVmean_size_1 > 50]) / length(wsCVmean_size_1) * 100

# size 2 summary
n_abs_leq_1per_size_2 <- length(wsCVmean_size_2[wsCVmean_size_2 <= 1])
n_abs_leq_5per_size_2 <- length(wsCVmean_size_2[wsCVmean_size_2 > 1 & wsCVmean_size_2 <= 5 ])
n_abs_leq_10per_size_2 <- length(wsCVmean_size_2[wsCVmean_size_2 > 5 & wsCVmean_size_2 <= 10 ])
n_abs_leq_50per_size_2 <- length(wsCVmean_size_2[wsCVmean_size_2 > 10 & wsCVmean_size_2 <= 50  ])
n_abs_geq_50per_size_2 <- length(wsCVmean_size_2[wsCVmean_size_2 > 50])

n_rel_leq_1per_size_2 <- length(wsCVmean_size_2[wsCVmean_size_2 <= 1]) / length(wsCVmean_size_2) * 100
n_rel_leq_5per_size_2 <- length(wsCVmean_size_2[wsCVmean_size_2 > 1 & wsCVmean_size_2 <= 5 ]) / length(wsCVmean_size_2) * 100
n_rel_leq_10per_size_2 <- length(wsCVmean_size_2[wsCVmean_size_2 > 5 & wsCVmean_size_2 <= 10 ]) / length(wsCVmean_size_2) * 100
n_rel_leq_50per_size_2 <- length(wsCVmean_size_2[wsCVmean_size_2 > 10 & wsCVmean_size_2 <= 50  ]) / length(wsCVmean_size_2) * 100
n_rel_geq_50per_size_2 <- length(wsCVmean_size_2[wsCVmean_size_2 > 50]) / length(wsCVmean_size_2) * 100

# size 3 summary
n_abs_leq_1per_size_3 <- length(wsCVmean_size_3[wsCVmean_size_3 <= 1])
n_abs_leq_5per_size_3 <- length(wsCVmean_size_3[wsCVmean_size_3 > 1 & wsCVmean_size_3 <= 5 ])
n_abs_leq_10per_size_3 <- length(wsCVmean_size_3[wsCVmean_size_3 > 5 & wsCVmean_size_3 <= 10 ])
n_abs_leq_50per_size_3 <- length(wsCVmean_size_3[wsCVmean_size_3 > 10 & wsCVmean_size_3 <= 50  ])
n_abs_geq_50per_size_3 <- length(wsCVmean_size_3[wsCVmean_size_3 > 50])

n_rel_leq_1per_size_3 <- length(wsCVmean_size_3[wsCVmean_size_3 <= 1]) / length(wsCVmean_size_3) * 100
n_rel_leq_5per_size_3 <- length(wsCVmean_size_3[wsCVmean_size_3 > 1 & wsCVmean_size_3 <= 5 ]) / length(wsCVmean_size_3) * 100
n_rel_leq_10per_size_3 <- length(wsCVmean_size_3[wsCVmean_size_3 > 5 & wsCVmean_size_3 <= 10 ]) / length(wsCVmean_size_3) * 100
n_rel_leq_50per_size_3 <- length(wsCVmean_size_3[wsCVmean_size_3 > 10 & wsCVmean_size_3 <= 50  ]) / length(wsCVmean_size_3) * 100
n_rel_geq_50per_size_3 <- length(wsCVmean_size_3[wsCVmean_size_3 > 50]) / length(wsCVmean_size_3) * 100

# size 4 summary
n_abs_leq_1per_size_4 <- length(wsCVmean_size_4[wsCVmean_size_4 <= 1])
n_abs_leq_5per_size_4 <- length(wsCVmean_size_4[wsCVmean_size_4 > 1 & wsCVmean_size_4 <= 5 ])
n_abs_leq_10per_size_4 <- length(wsCVmean_size_4[wsCVmean_size_4 > 5 & wsCVmean_size_4 <= 10 ])
n_abs_leq_50per_size_4 <- length(wsCVmean_size_4[wsCVmean_size_4 > 10 & wsCVmean_size_4 <= 50  ])
n_abs_geq_50per_size_4 <- length(wsCVmean_size_4[wsCVmean_size_4 > 50])

n_rel_leq_1per_size_4 <- length(wsCVmean_size_4[wsCVmean_size_4 <= 1]) / length(wsCVmean_size_4) * 100
n_rel_leq_5per_size_4 <- length(wsCVmean_size_4[wsCVmean_size_4 > 1 & wsCVmean_size_4 <= 5 ]) / length(wsCVmean_size_4) * 100
n_rel_leq_10per_size_4 <- length(wsCVmean_size_4[wsCVmean_size_4 > 5 & wsCVmean_size_4 <= 10 ]) / length(wsCVmean_size_4) * 100
n_rel_leq_50per_size_4 <- length(wsCVmean_size_4[wsCVmean_size_4 > 10 & wsCVmean_size_4 <= 50  ]) / length(wsCVmean_size_4) * 100
n_rel_geq_50per_size_4 <- length(wsCVmean_size_4[wsCVmean_size_4 > 50]) / length(wsCVmean_size_4) * 100

n_absolute_sizes <- rbind(
  c(n_abs_leq_1per_size_1, n_abs_leq_5per_size_1, n_abs_leq_10per_size_1, n_abs_leq_50per_size_1, n_abs_geq_50per_size_1), 
  c(n_abs_leq_1per_size_2, n_abs_leq_5per_size_2, n_abs_leq_10per_size_2, n_abs_leq_50per_size_2, n_abs_geq_50per_size_2), 
  c(n_abs_leq_1per_size_3, n_abs_leq_5per_size_3, n_abs_leq_10per_size_3, n_abs_leq_50per_size_3, n_abs_geq_50per_size_3), 
  c(n_abs_leq_1per_size_4, n_abs_leq_5per_size_4, n_abs_leq_10per_size_4, n_abs_leq_50per_size_4, n_abs_geq_50per_size_4))

n_relative_sizes <- rbind(
  c(n_rel_leq_1per_size_1, n_rel_leq_5per_size_1, n_rel_leq_10per_size_1, n_rel_leq_50per_size_1, n_rel_geq_50per_size_1), 
  c(n_rel_leq_1per_size_2, n_rel_leq_5per_size_2, n_rel_leq_10per_size_2, n_rel_leq_50per_size_2, n_rel_geq_50per_size_2), 
  c(n_rel_leq_1per_size_3, n_rel_leq_5per_size_3, n_rel_leq_10per_size_3, n_rel_leq_50per_size_3, n_rel_geq_50per_size_3), 
  c(n_rel_leq_1per_size_4, n_rel_leq_5per_size_4, n_rel_leq_10per_size_4, n_rel_leq_50per_size_4, n_rel_geq_50per_size_4))

rownames(n_absolute_sizes) <- c(27.1, 15.3, 6.8, 1.7)
rownames(n_relative_sizes) <- c(27.1, 15.3, 6.8, 1.7)

colnames(n_absolute_sizes) <- c("[0;1]", "]1;5]", "]5;10]", "]10;50]", "]50;inf]")
colnames(n_relative_sizes) <- c("[0;1]", "]1;5]", "]5;10]", "]10;50]", "]50;inf]")

n_absolute_sizes
n_relative_sizes

## calculate inserts wCV (1 - medium, 2 - medium, 3 - small, 4 - mix)
indeces_insert_1 <- c(1, 2, 3, 4)
indeces_insert_2 <- c(5, 6, 7, 8)
indeces_insert_3 <- c(9, 10, 11, 12)
indeces_insert_4 <- c(13, 14, 15, 16)
wsCVmean_insert_1 <- sapply(seq.int(dim(dataset_features_acq1)[2]), function(j) abs(100*agRee::agree.wscv(cbind(dataset_features_acq1[indeces_insert_1,j], dataset_features_acq2[indeces_insert_1,j]))$value))
names(wsCVmean_insert_1) <- colnames(dataset_features_acq1)

wsCVmean_insert_2 <- sapply(seq.int(dim(dataset_features_acq1)[2]), function(j) abs(100*agRee::agree.wscv(cbind(dataset_features_acq1[indeces_insert_2,j], dataset_features_acq2[indeces_insert_2,j]))$value))
names(wsCVmean_insert_2) <- colnames(dataset_features_acq1)

wsCVmean_insert_3 <- sapply(seq.int(dim(dataset_features_acq1)[2]), function(j) abs(100*agRee::agree.wscv(cbind(dataset_features_acq1[indeces_insert_3,j], dataset_features_acq2[indeces_insert_3,j]))$value))
names(wsCVmean_insert_3) <- colnames(dataset_features_acq1)

wsCVmean_insert_4 <- sapply(seq.int(dim(dataset_features_acq1)[2]), function(j) abs(100*agRee::agree.wscv(cbind(dataset_features_acq1[indeces_insert_4,j], dataset_features_acq2[indeces_insert_4,j]))$value))
names(wsCVmean_insert_4) <- colnames(dataset_features_acq1)

wsCVmean_inserts <- data.frame(wCV.insert1 = wsCVmean_insert_1, wCV.insert2 = wsCVmean_insert_2, wCV.insert3 = wsCVmean_insert_3, wCV.insert4 = wsCVmean_insert_4)

wsCVmean_inserts_hist_df <- data.frame(wCV=c(wsCVmean_insert_1, wsCVmean_insert_2, wsCVmean_insert_3, wsCVmean_insert_4), insert=c(rep(1, length(wsCVmean_insert_1)), rep(2, length(wsCVmean_insert_2)), rep(3, length(wsCVmean_insert_3)), rep(4, length(wsCVmean_insert_4))))
wsCVmean_inserts_hist_df$insert <- as.factor(wsCVmean_inserts_hist_df$insert) 
p_inserts_hist <- ggplot(wsCVmean_inserts_hist_df, aes(x=wCV, color=insert, fill=insert)) + 
  geom_histogram(binwidth=0.5, position="dodge", alpha=0.25) + labs(x = "wCV (%)") + xlim(0, 20)
print(p_inserts_hist)

# insert 1 summary
n_abs_leq_1per_insert_1 <- length(wsCVmean_insert_1[wsCVmean_insert_1 <= 1])
n_abs_leq_5per_insert_1 <- length(wsCVmean_insert_1[wsCVmean_insert_1 > 1 & wsCVmean_insert_1 <= 5 ])
n_abs_leq_10per_insert_1 <- length(wsCVmean_insert_1[wsCVmean_insert_1 > 5 & wsCVmean_insert_1 <= 10 ])
n_abs_leq_50per_insert_1 <- length(wsCVmean_insert_1[wsCVmean_insert_1 > 10 & wsCVmean_insert_1 <= 50  ])
n_abs_geq_50per_insert_1 <- length(wsCVmean_insert_1[wsCVmean_insert_1 > 50])

n_rel_leq_1per_insert_1 <- length(wsCVmean_insert_1[wsCVmean_insert_1 <= 1]) / length(wsCVmean_insert_1) * 100
n_rel_leq_5per_insert_1 <- length(wsCVmean_insert_1[wsCVmean_insert_1 > 1 & wsCVmean_insert_1 <= 5 ]) / length(wsCVmean_insert_1) * 100
n_rel_leq_10per_insert_1 <- length(wsCVmean_insert_1[wsCVmean_insert_1 > 5 & wsCVmean_insert_1 <= 10 ]) / length(wsCVmean_insert_1) * 100
n_rel_leq_50per_insert_1 <- length(wsCVmean_insert_1[wsCVmean_insert_1 > 10 & wsCVmean_insert_1 <= 50  ]) / length(wsCVmean_insert_1) * 100
n_rel_geq_50per_insert_1 <- length(wsCVmean_insert_1[wsCVmean_insert_1 > 50]) / length(wsCVmean_insert_1) * 100

# insert 2 summary
n_abs_leq_1per_insert_2 <- length(wsCVmean_insert_2[wsCVmean_insert_2 <= 1])
n_abs_leq_5per_insert_2 <- length(wsCVmean_insert_2[wsCVmean_insert_2 > 1 & wsCVmean_insert_2 <= 5 ])
n_abs_leq_10per_insert_2 <- length(wsCVmean_insert_2[wsCVmean_insert_2 > 5 & wsCVmean_insert_2 <= 10 ])
n_abs_leq_50per_insert_2 <- length(wsCVmean_insert_2[wsCVmean_insert_2 > 10 & wsCVmean_insert_2 <= 50  ])
n_abs_geq_50per_insert_2 <- length(wsCVmean_insert_2[wsCVmean_insert_2 > 50])

n_rel_leq_1per_insert_2 <- length(wsCVmean_insert_2[wsCVmean_insert_2 <= 1]) / length(wsCVmean_insert_2) * 100
n_rel_leq_5per_insert_2 <- length(wsCVmean_insert_2[wsCVmean_insert_2 > 1 & wsCVmean_insert_2 <= 5 ]) / length(wsCVmean_insert_2) * 100
n_rel_leq_10per_insert_2 <- length(wsCVmean_insert_2[wsCVmean_insert_2 > 5 & wsCVmean_insert_2 <= 10 ]) / length(wsCVmean_insert_2) * 100
n_rel_leq_50per_insert_2 <- length(wsCVmean_insert_2[wsCVmean_insert_2 > 10 & wsCVmean_insert_2 <= 50  ]) / length(wsCVmean_insert_2) * 100
n_rel_geq_50per_insert_2 <- length(wsCVmean_insert_2[wsCVmean_insert_2 > 50]) / length(wsCVmean_insert_2) * 100

# insert 3 summary
n_abs_leq_1per_insert_3 <- length(wsCVmean_insert_3[wsCVmean_insert_3 <= 1])
n_abs_leq_5per_insert_3 <- length(wsCVmean_insert_3[wsCVmean_insert_3 > 1 & wsCVmean_insert_3 <= 5 ])
n_abs_leq_10per_insert_3 <- length(wsCVmean_insert_3[wsCVmean_insert_3 > 5 & wsCVmean_insert_3 <= 10 ])
n_abs_leq_50per_insert_3 <- length(wsCVmean_insert_3[wsCVmean_insert_3 > 10 & wsCVmean_insert_3 <= 50  ])
n_abs_geq_50per_insert_3 <- length(wsCVmean_insert_3[wsCVmean_insert_3 > 50])

n_rel_leq_1per_insert_3 <- length(wsCVmean_insert_3[wsCVmean_insert_3 <= 1]) / length(wsCVmean_insert_3) * 100
n_rel_leq_5per_insert_3 <- length(wsCVmean_insert_3[wsCVmean_insert_3 > 1 & wsCVmean_insert_3 <= 5 ]) / length(wsCVmean_insert_3) * 100
n_rel_leq_10per_insert_3 <- length(wsCVmean_insert_3[wsCVmean_insert_3 > 5 & wsCVmean_insert_3 <= 10 ]) / length(wsCVmean_insert_3) * 100
n_rel_leq_50per_insert_3 <- length(wsCVmean_insert_3[wsCVmean_insert_3 > 10 & wsCVmean_insert_3 <= 50  ]) / length(wsCVmean_insert_3) * 100
n_rel_geq_50per_insert_3 <- length(wsCVmean_insert_3[wsCVmean_insert_3 > 50]) / length(wsCVmean_insert_3) * 100

# insert 4 summary
n_abs_leq_1per_insert_4 <- length(wsCVmean_insert_4[wsCVmean_insert_4 <= 1])
n_abs_leq_5per_insert_4 <- length(wsCVmean_insert_4[wsCVmean_insert_4 > 1 & wsCVmean_insert_4 <= 5 ])
n_abs_leq_10per_insert_4 <- length(wsCVmean_insert_4[wsCVmean_insert_4 > 5 & wsCVmean_insert_4 <= 10 ])
n_abs_leq_50per_insert_4 <- length(wsCVmean_insert_4[wsCVmean_insert_4 > 10 & wsCVmean_insert_4 <= 50  ])
n_abs_geq_50per_insert_4 <- length(wsCVmean_insert_4[wsCVmean_insert_4 > 50])

n_rel_leq_1per_insert_4 <- length(wsCVmean_insert_4[wsCVmean_insert_4 <= 1]) / length(wsCVmean_insert_4) * 100
n_rel_leq_5per_insert_4 <- length(wsCVmean_insert_4[wsCVmean_insert_4 > 1 & wsCVmean_insert_4 <= 5 ]) / length(wsCVmean_insert_4) * 100
n_rel_leq_10per_insert_4 <- length(wsCVmean_insert_4[wsCVmean_insert_4 > 5 & wsCVmean_insert_4 <= 10 ]) / length(wsCVmean_insert_4) * 100
n_rel_leq_50per_insert_4 <- length(wsCVmean_insert_4[wsCVmean_insert_4 > 10 & wsCVmean_insert_4 <= 50  ]) / length(wsCVmean_insert_4) * 100
n_rel_geq_50per_insert_4 <- length(wsCVmean_insert_4[wsCVmean_insert_4 > 50]) / length(wsCVmean_insert_4) * 100

n_absolute_inserts <- rbind(
  c(n_abs_leq_1per_insert_1, n_abs_leq_5per_insert_1, n_abs_leq_10per_insert_1, n_abs_leq_50per_insert_1, n_abs_geq_50per_insert_1), 
  c(n_abs_leq_1per_insert_2, n_abs_leq_5per_insert_2, n_abs_leq_10per_insert_2, n_abs_leq_50per_insert_2, n_abs_geq_50per_insert_2), 
  c(n_abs_leq_1per_insert_3, n_abs_leq_5per_insert_3, n_abs_leq_10per_insert_3, n_abs_leq_50per_insert_3, n_abs_geq_50per_insert_3), 
  c(n_abs_leq_1per_insert_4, n_abs_leq_5per_insert_4, n_abs_leq_10per_insert_4, n_abs_leq_50per_insert_4, n_abs_geq_50per_insert_4))

n_relative_inserts <- rbind(
  c(n_rel_leq_1per_insert_1, n_rel_leq_5per_insert_1, n_rel_leq_10per_insert_1, n_rel_leq_50per_insert_1, n_rel_geq_50per_insert_1), 
  c(n_rel_leq_1per_insert_2, n_rel_leq_5per_insert_2, n_rel_leq_10per_insert_2, n_rel_leq_50per_insert_2, n_rel_geq_50per_insert_2), 
  c(n_rel_leq_1per_insert_3, n_rel_leq_5per_insert_3, n_rel_leq_10per_insert_3, n_rel_leq_50per_insert_3, n_rel_geq_50per_insert_3), 
  c(n_rel_leq_1per_insert_4, n_rel_leq_5per_insert_4, n_rel_leq_10per_insert_4, n_rel_leq_50per_insert_4, n_rel_geq_50per_insert_4))

rownames(n_absolute_inserts) <- c("Ins.1", "Ins.2", "Ins.3", "Ins.4")
rownames(n_relative_inserts) <- c("Ins.1", "Ins.2", "Ins.3", "Ins.4")

colnames(n_absolute_inserts) <- c("[0;1]", "]1;5]", "]5;10]", "]10;50]", "]50;inf]")
colnames(n_relative_inserts) <- c("[0;1]", "]1;5]", "]5;10]", "]10;50]", "]50;inf]")

n_absolute_inserts
n_relative_inserts

# ## CCU 1.5T
# ccu_15T_test_retest_wCV_norepos <- n_relative_overall
# ccu_15T_test_retest_wCV_norepos_size <- n_relative_sizes
# ccu_15T_test_retest_wCV_norepos_insert <- n_relative_inserts
# 
# ccu_15T_test_retest_wCV_repos <- n_relative_overall
# ccu_15T_test_retest_wCV_repos_size <- n_relative_sizes
# ccu_15T_test_retest_wCV_repos_insert <- n_relative_inserts
# 
# ## CCU 3T
# ccu_3T_test_retest_wCV_norepos <- n_relative_overall
# ccu_3T_test_retest_wCV_norepos_size <- n_relative_sizes
# ccu_3T_test_retest_wCV_norepos_insert <- n_relative_inserts
# 
# ccu_3T_test_retest_wCV_repos <- n_relative_overall
# ccu_3T_test_retest_wCV_repos_size <- n_relative_sizes
# ccu_3T_test_retest_wCV_repos_insert <- n_relative_inserts
# 
# ## IEO 1.5T
# ieo_15T_test_retest_wCV_norepos <- n_relative_overall
# ieo_15T_test_retest_wCV_norepos_size <- n_relative_sizes
# ieo_15T_test_retest_wCV_norepos_insert <- n_relative_inserts
# 
# ieo_15T_test_retest_wCV_repos <- n_relative_overall
# ieo_15T_test_retest_wCV_repos_size <- n_relative_sizes
# ieo_15T_test_retest_wCV_repos_insert <- n_relative_inserts
# 
# 
# ## IEO CCU 1.5T
# ieo_ccu_15T_test_retest_wCV_norepos <- n_relative_overall
# ieo_ccu_15T_test_retest_wCV_norepos_size <- n_relative_sizes
# ieo_ccu_15T_test_retest_wCV_norepos_insert <- n_relative_inserts
# 
# ## CCU 1.5T CCU 3T
# ccu_3T_ccu_15T_test_retest_wCV_norepos <- n_relative_overall
# ccu_3T_ccu_15T_test_retest_wCV_norepos_size <- n_relative_sizes
# ccu_3T_ccu_15T_test_retest_wCV_norepos_insert <- n_relative_inserts

# ####

par(pty="s") ## pty sets the aspect ratio of the plot

## Overall Repeatability
Experiment <- c(rep("Scanner A",5), rep("Scanner A (repos.)",5), rep("Scanner B",5), rep("Scanner B (repos.)",5), rep("Scanner C",5), rep("Scanner C (repos.).",5))
wCV <- rep(c("wCV <= 1%" , "1% < wCV <= 5%" , "5% < wCV <= 10%", "10% < wCV <= 50%", "wCV > 50%") , 6)
Percentage <- c(ieo_15T_test_retest_wCV_norepos, ieo_15T_test_retest_wCV_repos, ccu_15T_test_retest_wCV_norepos, ccu_15T_test_retest_wCV_repos, ccu_3T_test_retest_wCV_norepos, ccu_3T_test_retest_wCV_repos)
data <- data.frame(Experiment,wCV,Percentage)
library(ggplot2)
library(viridis)

# ggplot(data, aes(fill=wCV, y=Percentage, x=Experiment)) +
#   geom_bar(position="stack", stat="identity") + theme(axis.text=element_text(size=14, colour = 'black'),
#                                                       axis.title=element_text(size=16, colour = 'black',face="bold"),
#                                                       legend.title = element_text(colour="black", size=14, face="bold"),
#                                                       legend.text = element_text(colour="black", size=12)) + 
#   scale_fill_discrete(breaks = c("wCV <= 1%" , "1% < wCV <= 5%" , "5% < wCV <= 10%", "10% < wCV <= 50%", "wCV > 50"))

ggplot(data, aes(fill=factor(wCV, levels=c("wCV <= 1%" , "1% < wCV <= 5%" , "5% < wCV <= 10%", "10% < wCV <= 50%", "wCV > 50%")), y=Percentage, x=Experiment)) +
  geom_bar(position="stack", stat="identity") + theme(axis.text=element_text(size=14, colour = 'black'),
                                                      axis.title=element_text(size=16, colour = 'black',face="bold"),
                                                      legend.title = element_text(colour="black", size=14, face="bold"),
                                                      legend.text = element_text(colour="black", size=12)) + 
  # scale_fill_discrete(breaks = c("wCV <= 1%" , "1% < wCV <= 5%" , "5% < wCV <= 10%", "10% < wCV <= 50%", "wCV > 50%")) + 
  labs(fill = "wCV")


## Overall Reproducibility
Experiment_Reproducibility_Scanners <- c(rep("A vs. B",5), rep("B vs. C",5))
wCV_Reproducibility_Scanners <- rep(c("wCV <= 1%" , "1% < wCV <= 5%" , "5% < wCV <= 10%", "10% < wCV <= 50%", "wCV > 50%") , 2)
Percentage_Reproducibility_Scanners <- c(ieo_ccu_15T_test_retest_wCV_norepos, ccu_3T_ccu_15T_test_retest_wCV_norepos)
data_reproducibility_scanners <- data.frame(Analysis=Experiment_Reproducibility_Scanners, wCV=wCV_Reproducibility_Scanners, Percentage=Percentage_Reproducibility_Scanners)
library(ggplot2)
library(viridis)

ggplot(data_reproducibility_scanners, aes(fill=factor(wCV, levels=c("wCV <= 1%" , "1% < wCV <= 5%" , "5% < wCV <= 10%", "10% < wCV <= 50%", "wCV > 50%")), y=Percentage, x=Analysis)) +
  geom_bar(position="stack", stat="identity")  + theme(axis.text=element_text(size=14, colour = 'black'),
                                                       axis.title=element_text(size=16, colour = 'black',face="bold"),
                                                       legend.title = element_text(colour="black", size=14, face="bold"),                                                   legend.text = element_text(colour="black", size=12)) + labs(fill = "wCV")

## DEPENDENCE ON SIZE
# df_wCV_size <- cbind(ccu_15T_test_retest_wCV_norepos_size, ccu_15T_test_retest_wCV_repos_size)
# df_wCV_size <- cbind(ccu_3T_test_retest_wCV_norepos_size, ccu_3T_test_retest_wCV_repos_size)
df_wCV_size <- cbind(ieo_15T_test_retest_wCV_norepos_size, ieo_15T_test_retest_wCV_repos_size)
vector_wCV_size <- as.vector(t(df_wCV_size))
Experiment <- c(rep("27.1",5), rep("27.1 (repos.)",5), rep("15.3",5), rep("15.3 (repos.)",5), rep("6.8",5), rep("6.8 (repos.)",5), rep("1.7",5), rep("1.7 (repos.)",5))
wCV <- rep(c("wCV <= 1%" , "1% < wCV <= 5%" , "5% < wCV <= 10%", "10% < wCV <= 50%", "wCV > 50%") , 8)
data_sizes <- data.frame(Experiment,wCV,vector_wCV_size)

ggplot(data_sizes, aes(fill=factor(wCV, levels=c("wCV <= 1%" , "1% < wCV <= 5%" , "5% < wCV <= 10%", "10% < wCV <= 50%", "wCV > 50%")), y=vector_wCV_size, x=factor(Experiment, levels=c("1.7", "1.7 (repos.)", "6.8", "6.8 (repos.)", "15.3", "15.3 (repos.)", "27.1", "27.1 (repos.)")))) +
  geom_bar(position="stack", stat="identity") + theme(axis.text=element_text(size=14, colour = 'black'),
                                                      axis.title=element_text(size=16, colour = 'black',face="bold"),
                                                      legend.title = element_text(colour="black", size=14, face="bold"),
                                                      legend.text = element_text(colour="black", size=12)) + 
  # scale_fill_discrete(breaks = c("wCV <= 1%" , "1% < wCV <= 5%" , "5% < wCV <= 10%", "10% < wCV <= 50%", "wCV > 50%")) + 
  labs(fill = "wCV") + labs(y="Percentage", x = "ROI Volume (cm³)")


## DEPENDENCE ON INSERT
# df_wCV_insert <- cbind(ccu_15T_test_retest_wCV_norepos_insert, ccu_15T_test_retest_wCV_repos_insert)
# df_wCV_insert <- cbind(ccu_3T_test_retest_wCV_norepos_insert, ccu_3T_test_retest_wCV_repos_insert)
df_wCV_insert <- cbind(ieo_15T_test_retest_wCV_norepos_insert, ieo_15T_test_retest_wCV_repos_insert)
vector_wCV_insert <- as.vector(t(df_wCV_insert))
Experiment <- c(rep("1",5), rep("1 (repos.)",5), rep("2",5), rep("2 (repos.)",5), rep("3",5), rep("3 (repos.)",5), rep("4",5), rep("4 (repos.)",5))
wCV <- rep(c("wCV <= 1%" , "1% < wCV <= 5%" , "5% < wCV <= 10%", "10% < wCV <= 50%", "wCV > 50%") , 8)
data_sizes <- data.frame(Experiment,wCV,vector_wCV_insert)

ggplot(data_sizes, aes(fill=factor(wCV, levels=c("wCV <= 1%" , "1% < wCV <= 5%" , "5% < wCV <= 10%", "10% < wCV <= 50%", "wCV > 50%")), y=vector_wCV_size, x=factor(Experiment, levels=c("1", "1 (repos.)", "2", "2 (repos.)", "3", "3 (repos.)", "4", "4 (repos.)")))) +
  geom_bar(position="stack", stat="identity") + theme(axis.text=element_text(size=14, colour = 'black'),
                                                      axis.title=element_text(size=16, colour = 'black',face="bold"),
                                                      legend.title = element_text(colour="black", size=14, face="bold"),
                                                      legend.text = element_text(colour="black", size=12)) + 
  # scale_fill_discrete(breaks = c("wCV <= 1%" , "1% < wCV <= 5%" , "5% < wCV <= 10%", "10% < wCV <= 50%", "wCV > 50%")) + 
  labs(fill = "wCV") + labs(y="Percentage", x = "Insert")

###########################################################
ieo_ccu_15T_test_retest_wCV_norepos_size
ieo_ccu_15T_test_retest_wCV_norepos_insert

# ## CCU 1.5T CCU 3T
# ccu_3T_ccu_15T_test_retest_wCV_norepos <- n_relative_overall
# ccu_3T_ccu_15T_test_retest_wCV_norepos_size <- n_relative_sizes
# ccu_3T_ccu_15T_test_retest_wCV_norepos_insert

## Dependence on Size Reproducibility
df_wCV_size <- cbind(ieo_ccu_15T_test_retest_wCV_norepos_size,ccu_3T_ccu_15T_test_retest_wCV_norepos_size)
vector_wCV_size <- as.vector(t(df_wCV_size))
Experiment <- c(rep("27.1 (A vs. B)",5), rep("27.1 (B vs. C)",5), rep("15.3 (A vs. B)",5), rep("15.3 (B vs. C)",5), rep("6.8 (A vs. B)",5), rep("6.8 (B vs. C)",5), rep("1.7 (A vs. B)",5), rep("1.7 (B vs. C)",5))
wCV <- rep(c("wCV <= 1%" , "1% < wCV <= 5%" , "5% < wCV <= 10%", "10% < wCV <= 50%", "wCV > 50%") , 8)
data_sizes <- data.frame(Experiment,wCV,vector_wCV_size)

ggplot(data_sizes, aes(fill=factor(wCV, levels=c("wCV <= 1%" , "1% < wCV <= 5%" , "5% < wCV <= 10%", "10% < wCV <= 50%", "wCV > 50%")), y=vector_wCV_size, x=factor(Experiment, levels=c("1.7 (A vs. B)", "1.7 (B vs. C)", "6.8 (A vs. B)", "6.8 (B vs. C)", "15.3 (A vs. B)", "15.3 (B vs. C)", "27.1 (A vs. B)", "27.1 (B vs. C)")))) +
  geom_bar(position="stack", stat="identity") + theme(axis.text=element_text(size=14, colour = 'black'),
                                                      axis.title=element_text(size=16, colour = 'black',face="bold"),
                                                      legend.title = element_text(colour="black", size=14, face="bold"),
                                                      legend.text = element_text(colour="black", size=12)) + 
  # scale_fill_discrete(breaks = c("wCV <= 1%" , "1% < wCV <= 5%" , "5% < wCV <= 10%", "10% < wCV <= 50%", "wCV > 50%")) + 
  labs(fill = "wCV") + labs(y="Percentage", x = "ROI Volume (cm³)")



## Dependence on Insert Reproducibility
df_wCV_insert <- cbind(ieo_ccu_15T_test_retest_wCV_norepos_insert,ccu_3T_ccu_15T_test_retest_wCV_norepos_insert)
vector_wCV_insert <- as.vector(t(df_wCV_insert))
Experiment <- c(rep("1 (A vs. B)",5), rep("1 (B vs. C)",5), rep("2 (A vs. B)",5), rep("2 (B vs. C)",5), rep("3 (A vs. B)",5), rep("3 (B vs. C)",5), rep("4 (A vs. B)",5), rep("4 (B vs. C)",5))
wCV <- rep(c("wCV <= 1%" , "1% < wCV <= 5%" , "5% < wCV <= 10%", "10% < wCV <= 50%", "wCV > 50%") , 8)
data_sizes <- data.frame(Experiment,wCV,vector_wCV_size)

ggplot(data_sizes, aes(fill=factor(wCV, levels=c("wCV <= 1%" , "1% < wCV <= 5%" , "5% < wCV <= 10%", "10% < wCV <= 50%", "wCV > 50%")), y=vector_wCV_size, x=factor(Experiment, levels=c("1 (A vs. B)", "1 (B vs. C)", "2 (A vs. B)", "2 (B vs. C)", "3 (A vs. B)", "3 (B vs. C)", "4 (A vs. B)", "4 (B vs. C)")))) +
  geom_bar(position="stack", stat="identity") + theme(axis.text=element_text(size=14, colour = 'black'),
                                                      axis.title=element_text(size=16, colour = 'black',face="bold"),
                                                      legend.title = element_text(colour="black", size=14, face="bold"),
                                                      legend.text = element_text(colour="black", size=12)) + 
  # scale_fill_discrete(breaks = c("wCV <= 1%" , "1% < wCV <= 5%" , "5% < wCV <= 10%", "10% < wCV <= 50%", "wCV > 50%")) + 
  labs(fill = "wCV") + labs(y="Percentage", x = "Insert")
