#This code is used to generate the heatmap (2nd and 3rd degree invariants calculated for different
#test n and K values) for stationary distributions generated from a Hill function with true n = 1
#(hill-function-batch_n_1.csv).
#hill-function-batch_2 then is for stationary distributions generated from a process with true n = 2
#and so on. Note that x and y in the code refer to x_2 and x_1 respectively in the paper.

library(tidyverse)
library(data.table)
library(parallel)
#This directory must be created, with all of the stationary distributions generated from true n = 1,
#before the code is run.
setwd("Hill_Stationary_Correct//n_1")

#This function is used to calcuate the 2nd and 3rd degree invariant errors for a specific test n, K.
#x_y_df_stationary denotes the stationary distribution. n and K denote the test n and K values.
#x_bar, x2_bar, x3_bar denote E[x], E[x^2], and E[x^3] respectively. These must be calcuated beforehand
#and fed into the function.
calculate_error_hill <- function(x_y_df_stationary, n, K, x_bar, x2_bar, x3_bar) {

  xyn_yn_Kn_bar <- sum(x_y_df_stationary$x*x_y_df_stationary$y^n/(x_y_df_stationary$y^n+K^n)*x_y_df_stationary$t)/sum(x_y_df_stationary$t)
  x2yn_yn_Kn_bar <- sum(x_y_df_stationary$x^2*x_y_df_stationary$y^n/(x_y_df_stationary$y^n+K^n)*x_y_df_stationary$t)/sum(x_y_df_stationary$t)
  yn_yn_Kn_bar <- sum(x_y_df_stationary$y^n/(x_y_df_stationary$y^n+K^n)*x_y_df_stationary$t)/sum(x_y_df_stationary$t)
  
  lhs_AB <- xyn_yn_Kn_bar/yn_yn_Kn_bar + 1
  rhs_AB <- x2_bar/x_bar
  
  error_AB <- 2*(lhs_AB-rhs_AB)/(lhs_AB+rhs_AB)
  
  lhs_BC <- x2yn_yn_Kn_bar/yn_yn_Kn_bar + 2*x2_bar + xyn_yn_Kn_bar/yn_yn_Kn_bar + x2_bar/x_bar
  rhs_BC <- x3_bar/x_bar + 2*xyn_yn_Kn_bar/yn_yn_Kn_bar*x_bar + 2*x_bar
  
  error_BC <- 2*(lhs_BC-rhs_BC)/(lhs_BC+rhs_BC)
  
  return(list(error_AB, error_BC))
}

#This function is used to iterate across different test n, K pairs as well as different stationary
#distributions generated from different true n, K pairs. It calls the previous function
#for every iteration, returning a list containing the 2nd and 3rd degree invariant errors for all
#tested n, K pairs for every stationary distribution. 
#design_matrix denotes a list with test n, K values to evaluate the invariant errors at.
calculate_error_iter_hill <- function(design_matrix) {
  
  file <- as.character(design_matrix$file)
  n <- design_matrix$n
  K <- design_matrix$K
  
  x_y_df_stationary_ <- read.csv(file)
  colnames(x_y_df_stationary_) <- c("x", "y", "t", "x_bar", "y_bar", "gamma", "true_n", "true_K")
  
  x_bar <- sum(x_y_df_stationary_$x*x_y_df_stationary_$t)/sum(x_y_df_stationary_$t)
  x2_bar <- sum(x_y_df_stationary_$x^2*x_y_df_stationary_$t)/sum(x_y_df_stationary_$t)
  x3_bar <- sum(x_y_df_stationary_$x^3*x_y_df_stationary_$t)/sum(x_y_df_stationary_$t)
  
  error <- calculate_error_hill(x_y_df_stationary_, n, K, x_bar, x2_bar, x3_bar)
  
  return(error)
}

#Listing all the files in the working directory (Hill_Stationary_Correct//n_1)
files <- list.files(pattern="*.csv")
n <- seq(1, 10, 0.1)
K <- seq(10, 100, 1)

#Creating the design_matrix
design_matrix <- files %>%
  merge(n, by = NULL) %>%
  merge(K, by = NULL)

colnames(design_matrix) <- c("file", "n", "K")

#Coercing the design matrix into a list
d_ <- list()

for (i in seq(nrow(design_matrix))) {
  l <- as.list(design_matrix[i,])
  d_[[i]] <- l
}

#Applying calculate_error_iter_hill to the design_matrix. This function requires parallelization.
set.seed(17)
error <- mclapply(d_, calculate_error_iter_hill, mc.cores = 160)
error <- t(matrix(unlist(error), nrow=length(unlist(error[1]))))

#Attatching invariant errors to the design_matrix
design_matrix <- design_matrix %>% mutate(error_AB = error[ ,1], error_BC = error[,2])

setwd("..//..")
#Writing csv
write.csv(design_matrix, "hill-function-batch_n_1.csv", row.names = FALSE)
