library(dplyr)
library(data.table)
library(parallel)
setwd("hill_stationary_distros")

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

calculate_error_iter_hill_sample <- function(design_matrix, x_y_df_stationary, intended_x_bar, intended_y_bar, intended_gamma, intended_true_n, intended_true_K, samp_size) {
  
  n <- design_matrix$n
  K <- design_matrix$K
  seed <- design_matrix$index
  
  x_y_df_stationary_ <- x_y_df_stationary %>% 
    filter(x_bar == intended_x_bar & y_bar == intended_y_bar & gamma == intended_gamma & true_n == intended_true_n & true_K == intended_true_K)
  
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

mod_sample <- function(x, size, replace, prob, seed) {
  set.seed(seed)
  return(sample(x, size, replace, prob))
}

hill_tbl <- 
  list.files(pattern = "*.csv") %>% 
  map_df(~fread(.))

colnames(hill_tbl) <- c("x", "y", "t", "x_bar", "y_bar", "gamma", "true_n", "true_K")



index <- seq(10)
n <- seq(1, 10, 1)
K <- seq(10, 100, 10)

design_matrix <- merge(index, n, by = NULL) %>%
  merge(K, by = NULL)

colnames(design_matrix) <- c("index", "n", "K")

d_ <- list()

for (i in seq(nrow(design_matrix))) {
  l <- as.list(design_matrix[i,])
  d_[[i]] <- l
}

set.seed(17)
error <- lapply(d_, calculate_error_iter_hill_sample, x_y_df_stationary = hill_tbl, intended_x_bar = 50, intended_y_bar = 50, intended_gamma = 1, intended_true_n = 3, intended_true_K = 75, samp_size = 500)

error <- t(matrix(unlist(error), nrow=length(unlist(error[1]))))
design_matrix <- design_matrix %>% mutate(error_AB = error[ ,1], error_BC = error[,2])



setwd("..")
write.csv(design_matrix, "hill-function-sample.csv", row.names = FALSE)
design_matrix <- read.csv("hill-function-sample.csv")
