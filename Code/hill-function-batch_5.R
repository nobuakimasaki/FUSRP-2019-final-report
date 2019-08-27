library(tidyverse)
library(data.table)
library(parallel)
setwd("Hill_Stationary_Correct//n_5")

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

files <- list.files(pattern="*.csv")
n <- seq(1, 10, 0.1)
K <- seq(10, 100, 1)

design_matrix <- files %>%
  merge(n, by = NULL) %>%
  merge(K, by = NULL)

colnames(design_matrix) <- c("file", "n", "K")

#design_matrix <- design_matrix[1:2,] 

d_ <- list()

for (i in seq(nrow(design_matrix))) {
  l <- as.list(design_matrix[i,])
  d_[[i]] <- l
}

set.seed(17)
error <- mclapply(d_, calculate_error_iter_hill, mc.cores = 160)
error <- t(matrix(unlist(error), nrow=length(unlist(error[1]))))
design_matrix <- design_matrix %>% mutate(error_AB = error[ ,1], error_BC = error[,2])

setwd("..//..")
write.csv(design_matrix, "hill-function-batch_n_5.csv", row.names = FALSE)
