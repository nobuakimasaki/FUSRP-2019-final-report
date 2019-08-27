library(tidyverse)
library(data.table)
library(parallel)
setwd("Hill_Stationary_Correct//n_3")

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

calculate_error_iter_hill <- function(design_matrix, x_y_df_stationary, x_bar, x2_bar, x3_bar) {
  
  n <- design_matrix$n
  K <- design_matrix$K
  
  error <- calculate_error_hill(x_y_df_stationary, n, K, x_bar, x2_bar, x3_bar)
  
  return(error)
}

Generate_Poisson <- function(x_bar, x_max) {
  x <- seq(0, x_max)
  x <- as.data.frame(x)
  x$t <- dpois(x$x, x_bar, log = FALSE)
  return(x)
}

KLD <- function(x_y_df_stationary) {
  x_y_df_stationary <- x_y_df_stationary %>%
    group_by(x) %>%
    summarize(t = sum(t)) %>%
    mutate(t = t/sum(t))
  x_bar <- sum(x_y_df_stationary$x*x_y_df_stationary$t)/sum(x_y_df_stationary$t)
  x_max <- max(x_y_df_stationary$x)
  poisson <- Generate_Poisson(x_bar, x_max)
  d_ <- left_join(x_y_df_stationary, poisson, by = "x")
  colnames(d_) <- c("x", "t1", "t2")
  #d_<-d_[(d_$t2>0),]
  #d_<-d_[d_$t2>0.000000000000000000001,]
  KLD <- sum(d_$t1*log(d_$t1/d_$t2))
  return(KLD)
}

BC <- function(x_y_df_stationary) {
  x_y_df_stationary <- x_y_df_stationary %>%
    group_by(x) %>%
    summarize(t = sum(t)) %>%
    mutate(t = t/sum(t))
  x_bar <- sum(x_y_df_stationary$x*x_y_df_stationary$t)/sum(x_y_df_stationary$t)
  x_max <- max(x_y_df_stationary$x)
  poisson <- Generate_Poisson(x_bar, x_max)
  d_ <- left_join(x_y_df_stationary, poisson, by = "x")
  colnames(d_) <- c("x", "t1", "t2")
  BC <- sum(sqrt(d_$t1*d_$t2))
  return(BC)
}

Fc <- function(sta, n, K){
  x_bar <- sum(sta$x*sta$t)/sum(sta$t)
  xyn_yn_Kn_bar <- sum(sta$x*sta$y^n/(sta$y^n+K^n)*sta$t)/sum(sta$t)
  yn_yn_Kn_bar <- sum(sta$y^n/(sta$y^n+K^n)*sta$t)/sum(sta$t)
  
  cov_x_yn_yn_Kn <- xyn_yn_Kn_bar - x_bar*yn_yn_Kn_bar
  eta_xR_plus <- cov_x_yn_yn_Kn/x_bar/yn_yn_Kn_bar
  
  Fc <- eta_xR_plus/(1/x_bar + eta_xR_plus)
  return(Fc)
}

gradient_hill <- function(design_matrix, true_n, true_K) {
  n_plus <- design_matrix %>% filter(n < true_n + 0.10001 & n > true_n + 0.09999 & K == true_K)
  n_minus <- design_matrix %>% filter(n < true_n - 0.09999 & n > true_n - 0.10001 & K == true_K)
  K_plus <- design_matrix %>% filter(n == true_n & K == true_K + 1)
  K_minus <- design_matrix %>% filter(n == true_n & K == true_K - 1)
  
  derrorAB_dn <- (n_plus$error_AB - n_minus$error_AB)/0.2
  derrorBC_dn <- (n_plus$error_BC - n_minus$error_BC)/0.2
  
  derrorAB_dK <- (K_plus$error_AB - K_minus$error_AB)/2
  derrorBC_dK <- (K_plus$error_BC - K_minus$error_BC)/2
  
  # gradient_AB <- sqrt(derrorAB_dn^2 + derrorAB_dK^2)
  # gradient_BC <- sqrt(derrorBC_dn^2 + derrorBC_dK^2)
  # 
  # cos_similarity <- (derrorAB_dn*derrorBC_dn + derrorAB_dK*derrorBC_dK)/sqrt(derrorAB_dn^2 + derrorAB_dK^2)/sqrt(derrorBC_dn^2 + derrorBC_dK^2)
  
  return(list(derrorAB_dn, derrorBC_dn, derrorAB_dK, derrorBC_dK))
}

calculate_error_hill_batch <- function(sta, x_bar, x2_bar, x3_bar, n, K) {
  
  KLD <- KLD(sta)
  BC <- BC(sta)
  Fc <- Fc(sta, n, K)
  
  n_list <- c(n - 0.1, n, n + 0.1)
  K_list <- c(K - 1, K, K + 1)
  
  design_matrix <- n_list %>%
    merge(K_list, by = NULL)
  colnames(design_matrix) <- c("n", "K")
  design_matrix <- as.data.frame(design_matrix)
  
  d_ <- list()
  for (i in seq(nrow(design_matrix))) {
    l <- as.list(design_matrix[i,])
    d_[[i]] <- l
  }
  
  error <- lapply(d_, calculate_error_iter_hill, 
                  x_y_df_stationary = sta, 
                  x_bar = x_bar, 
                  x2_bar = x2_bar, 
                  x3_bar = x3_bar)
  
  error_ABC <- t(matrix(unlist(error), nrow=length(unlist(error[1]))))
  design_matrix$error_AB <- error_ABC[ ,1]
  design_matrix$error_BC <- error_ABC[ ,2]
  
  gradient_list <- gradient_hill(design_matrix, n, K)
  gradient <- t(matrix(unlist(gradient_list), nrow=length(unlist(gradient_list[1]))))
  derrorAB_dn <- gradient[[1]]
  derrorBC_dn <- gradient[[2]]
  derrorAB_dK <- gradient[[3]]
  derrorBC_dK <- gradient[[4]]
  
  return(list(KLD, BC, Fc, derrorAB_dn, derrorBC_dn, derrorAB_dK, derrorBC_dK, x_bar))
}

calculate_error_iter_hill_batch <- function(design_matrix, x_y_df_stationary) {
  
  intended_x_bar <- design_matrix$x_bar
  intended_y_bar <- design_matrix$y_bar
  intended_gamma <- design_matrix$gamma
  intended_true_n <- design_matrix$true_n
  intended_true_K <- design_matrix$true_K
  
  x_y_df_stationary_ <- x_y_df_stationary %>% 
    filter(x_bar == intended_x_bar & y_bar == intended_y_bar & gamma == intended_gamma & true_n == intended_true_n & true_K == intended_true_K)
  
  x_bar <- sum(x_y_df_stationary_$x*x_y_df_stationary_$t)/sum(x_y_df_stationary_$t)
  x2_bar <- sum(x_y_df_stationary_$x^2*x_y_df_stationary_$t)/sum(x_y_df_stationary_$t)
  x3_bar <- sum(x_y_df_stationary_$x^3*x_y_df_stationary_$t)/sum(x_y_df_stationary_$t)
  
  d_ <- calculate_error_hill_batch(x_y_df_stationary_, x_bar, x2_bar, x3_bar, intended_true_n, intended_true_K)
  
  return(d_)
}

hill_tbl <- 
  list.files(pattern = "*.csv") %>% 
  map_df(~fread(.))

colnames(hill_tbl) <- c("x", "y", "t", "x_bar", "y_bar", "gamma", "true_n", "true_K")

design_matrix <- hill_tbl %>% 
  group_by(x_bar, y_bar, gamma, true_n, true_K) %>% 
  filter(row_number() == 1) %>%
  select(x_bar, y_bar, gamma, true_n, true_K)

colnames(design_matrix) <- c("x_bar", "y_bar", "gamma", "true_n", "true_K")

d_ <- list()

for (i in seq(nrow(design_matrix))) {
  l <- as.list(design_matrix[i,])
  d_[[i]] <- l
}

set.seed(17)
stats_ <- mclapply(d_, calculate_error_iter_hill_batch, x_y_df_stationary = hill_tbl, mc.cores = 160)
stats_ <- t(matrix(unlist(stats_), nrow=length(unlist(stats_[1]))))

design_matrix$KLD <- stats_[ ,1]
design_matrix$BC <- stats_[,2]
design_matrix$Fc <- stats_[,3]
design_matrix$derrorAB_dn <- stats_[,4]
design_matrix$derrorBC_dn <- stats_[,5]
design_matrix$derrorAB_dK <- stats_[,6]
design_matrix$derrorBC_dK <- stats_[,7]
design_matrix$emp_x_bar <- stats_[,8]

setwd("..//..")
write.csv(design_matrix, "hill-function-statistics_n_3.csv")