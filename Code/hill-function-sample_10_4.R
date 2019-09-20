#This code is used to generate the data required for all the plots involving sampling from stationary
#distributions generated from a Hill function such as hill-function-sample_50_50_0.1_5_70_10-4_1000.csv
#(the suffix 10_4 refers to the fact that this code
#samples 10^4 points from the specified stationary distribution).
#hill-function-sample_10_5 then is for sampling 10^5 points, hill-function-sample_5_10_4 5*10^4
#points and so on. Note that x and y in the code refer to x_2 and x_1 respectively in the paper.

library(tidyverse)
library(data.table)
library(parallel)
#This directory must be created, with the specific stationary distribution that we will sample from
#using this script.
setwd("Hill_Stationary_Correct//n_5")

#This function is used to calcuate the 2nd and 3rd degree invariant errors for a specific test n, K.
#x_y_df_stationary denotes the stationary distribution. n and K denote the test n and K values.
#x_bar, x2_bar, x3_bar denote E[x], E[x^2], and E[x^3] respectively. These must be calcuated beforehand
#and fed into the function.
calculate_error_hill <- function(x_y_df_stationary, n, K, x_bar, x2_bar, x3_bar) {
  
  xyn_yn_Kn_bar <- sum(x_y_df_stationary$x*x_y_df_stationary$y^n/(x_y_df_stationary$y^n+K^n)*x_y_df_stationary$t)/sum(x_y_df_stationary$t)
  x2yn_yn_Kn_bar <- sum(x_y_df_stationary$x^2*x_y_df_stationary$y^n/(x_y_df_stationary$y^n+K^n)*x_y_df_stationary$t)/sum(x_y_df_stationary$t)
  yn_yn_Kn_bar <- sum(x_y_df_stationary$y^n/(x_y_df_stationary$y^n+K^n)*x_y_df_stationary$t)/sum(x_y_df_stationary$t)
  
  cov_XY <- xyn_yn_Kn_bar - x_bar*yn_yn_Kn_bar
  
  lhs_AB <- xyn_yn_Kn_bar/yn_yn_Kn_bar + 1
  rhs_AB <- x2_bar/x_bar
  
  error_AB <- 2*(lhs_AB-rhs_AB)/(lhs_AB+rhs_AB)
  
  lhs_BC <- x2yn_yn_Kn_bar/yn_yn_Kn_bar + 2*x2_bar + xyn_yn_Kn_bar/yn_yn_Kn_bar + x2_bar/x_bar
  rhs_BC <- x3_bar/x_bar + 2*xyn_yn_Kn_bar/yn_yn_Kn_bar*x_bar + 2*x_bar
  
  error_BC <- 2*(lhs_BC-rhs_BC)/(lhs_BC+rhs_BC)
  
  lhs_CA <- xyn_yn_Kn_bar/yn_yn_Kn_bar + x2_bar/x_bar + x2yn_yn_Kn_bar/yn_yn_Kn_bar
  rhs_CA <- x3_bar/x_bar
  
  error_CA <- 2*(lhs_CA-rhs_CA)/(lhs_CA+rhs_CA)
  
  return(list(error_AB, error_BC, error_CA, cov_XY))
}

#This function first isolates the stationary distribution that we want to sample from.
#This stationary distribution is specified with the arguments intended_x_bar, intended_y_bar
#intended_gamma, intended_true_n, and intended_true_K.

#Each iteration of this function then denotes a specific sample from this stationary distribution,
#as well as a specific test n, K value. It calls the previous function
#for every iteration, returning a list containing the 2nd and 3rd degree invariant errors for all
#tested n, K pairs for each sample. 

#design_matrix denotes a list with test n, K values to evaluate the invariant errors at as well
#as the seed used for sampling from the stationary distribution.
#x_y_df_stationary denotes the stationary distributions in the directory Hill_Stationary_Correct//n_5
#binded together.
#samp_size denotes the number of samples drawn for each iteration.
calculate_error_iter_hill_sample <- function(design_matrix, x_y_df_stationary, intended_x_bar, intended_y_bar, intended_gamma, intended_true_n, intended_true_K, samp_size) {
  
  n <- design_matrix$n
  K <- design_matrix$K
  seed <- design_matrix$index
  
  x_y_df_stationary_ <- x_y_df_stationary %>% 
    filter(x_bar == intended_x_bar & y_bar == intended_y_bar & gamma == intended_gamma & true_n == intended_true_n & true_K == intended_true_K)
  
  sumt <- sum(x_y_df_stationary_$t)
  x_y_df_stationary_$t <- x_y_df_stationary_$t/sumt
  
  samp <- mod_sample(seq(nrow(x_y_df_stationary_)), samp_size, replace = T, x_y_df_stationary_$t, seed)
  x <- x_y_df_stationary_$x[samp]
  y <- x_y_df_stationary_$y[samp]
  sample <- as.data.frame(cbind(x, y))
  colnames(sample) <- c("x", "y")
  sample$t <- 1
  sumt <- nrow(sample)
  sample <- sample %>% 
    group_by(x,y) %>%
    summarize(t = sum(t)/sumt)
  x_bar <- sum(sample$x*sample$t)/sum(sample$t)
  x2_bar <- sum(sample$x^2*sample$t)/sum(sample$t)
  x3_bar <- sum(sample$x^3*sample$t)/sum(sample$t)
  
  error <- calculate_error_hill(sample, n, K, x_bar, x2_bar, x3_bar)
  
  return(error)
}

#Helper function to obtain same sample by utilizing the same seed. Used in the previous function.
mod_sample <- function(x, size, replace, prob, seed) {
  set.seed(seed)
  return(sample(x, size, replace, prob))
}

#Reading all of the files in Hill_Stationary_Correct//n_5 
#(if we are sampling from a distribution where n = 5)
hill_tbl <- 
  list.files(pattern = "*.csv") %>% 
  map_df(~fread(.))

colnames(hill_tbl) <- c("x", "y", "t", "x_bar", "y_bar", "gamma", "true_n", "true_K")

#Creating the design matrix
index <- seq(100)
n <- seq(1, 10, 0.1)
K <- seq(10, 100, 1)

design_matrix <- merge(index, n, by = NULL) %>%
  merge(K, by = NULL)

colnames(design_matrix) <- c("index", "n", "K")

#Coercing the design matrix into a list
d_ <- list()

for (i in seq(nrow(design_matrix))) {
  l <- as.list(design_matrix[i,])
  d_[[i]] <- l
}

#Applying calculate_error_iter_hill_sample to the design_matrix. 
#The sample size is set 10^4, but this can be changed
set.seed(17)
error <- mclapply(d_, calculate_error_iter_hill_sample, x_y_df_stationary = hill_tbl, intended_x_bar = 50, intended_y_bar = 50, intended_gamma = 0.1, intended_true_n = 5, intended_true_K = 70, samp_size = 10^4, mc.cores = 160)
#Attatching invariant errors to the design_matrix
error <- t(matrix(unlist(error), nrow=length(unlist(error[1]))))
design_matrix <- design_matrix %>% mutate(error_AB = error[ ,1], error_BC = error[,2])
#Writing the data
setwd("..//..")
write.csv(design_matrix, "hill-function-sample_50_50_0.1_5_70_10-4_1000.csv", row.names = FALSE)
