library(tidyverse)
library(data.table)
library(parallel)
setwd("hill_stationary_distros")

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

calculate_error_iter_hill <- function(design_matrix, x_y_df_stationary) {
  
  intended_x_bar <- design_matrix$x_bar
  intended_y_bar <- design_matrix$y_bar
  intended_gamma <- design_matrix$gamma
  intended_true_n <- design_matrix$true_n
  intended_true_K <- design_matrix$true_K
  n <- design_matrix$n
  K <- design_matrix$K
  
  x_y_df_stationary_ <- x_y_df_stationary %>% 
    dplyr::filter(x_bar == intended_x_bar & y_bar == intended_y_bar & gamma == intended_gamma & true_n == intended_true_n & true_K == intended_true_K)
  
  x_bar <- sum(x_y_df_stationary_$x*x_y_df_stationary_$t)/sum(x_y_df_stationary_$t)
  x2_bar <- sum(x_y_df_stationary_$x^2*x_y_df_stationary_$t)/sum(x_y_df_stationary_$t)
  x3_bar <- sum(x_y_df_stationary_$x^3*x_y_df_stationary_$t)/sum(x_y_df_stationary_$t)
  
  error <- calculate_error_hill(x_y_df_stationary_, n, K, x_bar, x2_bar, x3_bar)
  
  return(error)
}

hill_tbl <- 
  list.files(pattern = "*.csv") %>% 
  map_df(~fread(.))

colnames(hill_tbl) <- c("x", "y", "t", "x_bar", "y_bar", "gamma", "true_n", "true_K")


n <- seq(1, 10, 5)
K <- seq(10, 100, 50)

design_matrix <- hill_tbl %>% 
  group_by(x_bar, y_bar, gamma, true_n, true_K) %>% 
  filter(row_number() == 1) %>%
  select(x_bar, y_bar, gamma, true_n, true_K) %>%
  merge(n, by = NULL) %>%
  merge(K, by = NULL)

colnames(design_matrix) <- c("x_bar", "y_bar", "gamma", "true_n", "true_K", "n", "K")

#design_matrix <- design_matrix[1:2,] testing

d_ <- list()

for (i in seq(nrow(design_matrix))) {
  l <- as.list(design_matrix[i,])
  d_[[i]] <- l
}

set.seed(17)
error <- lapply(d_, calculate_error_iter_hill, x_y_df_stationary = hill_tbl)
error <- t(matrix(unlist(error), nrow=length(unlist(error[1]))))
design_matrix <- design_matrix %>% mutate(error_AB = error[ ,1], error_BC = error[,2])


setwd("..")
write.csv(design_matrix, "hill-function-batch.csv", row.names = FALSE)
