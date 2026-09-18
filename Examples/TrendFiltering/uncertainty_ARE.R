
 ########## Uncertainty Relative efficiency  ##########
 
 ################################################################################
 ################## Trendfiltering example output visualisation #################
 ################################################################################
 rm(list = ls())
 library(ggplot2)
 source("TF_functions.R")
 
load("output_mala.Rdata")
 
 #### Marginal variance comparison for IS vs PxMALA
 x <- c(1:100)
 dim <- 100 
 margvar_ism <- matrix(0, nrow = 100, ncol = dim)
 margvar_pxm <- matrix(0, nrow = 100, ncol = dim)
 for (i in 1:100)
 {
   asympmat_ism <- matrix(unlist(output_mala[[i]][[3]]), nrow = dim, ncol = dim, byrow = T)
   asympmat_pxm <- matrix(unlist(output_mala[[i]][[4]]), nrow = dim, ncol = dim, byrow = T)
   for (j in 1:dim)
   {
     margvar_ism[i,j] <- asympmat_ism[j,j]
     margvar_pxm[i,j] <- asympmat_pxm[j,j]
   }
 }
 rel_eff_mat_mala <- margvar_pxm/margvar_ism
 
 #######  Boxplots
 
 boxplot( rel_eff_mat_mala, outline = FALSE, xlab = "Component",
          ylab = "Relative efficiency")
 abline(h = 1, lty = 2)
 
 ############ Histogram
 
 mean_rel_eff <- colMeans(rel_eff_mat_mala)
 
 hist(mean_rel_eff, breaks = 20,xlab = "Mean relative efficiency",
                              main = "Distribution across components")
 abline(v = 1, lty = 2)
 
 ########  Band
 
 mean_rel_eff <- colMeans(rel_eff_mat_mala)
 
 q <- apply(rel_eff_mat_mala, 2,quantile, probs = c(0.025, 0.5, 0.975))
 
 plot(q[2, ], type = "p", pch = 19,ylim = range(q), xlab = "Component",
                                               ylab = "Relative efficiency")
 segments(x0 = 1:ncol(rel_eff_mat_mala), y0 = q[1, ], 
              x1 = 1:ncol(rel_eff_mat_mala), y1 = q[3, ])

 
 
 
 
 
 
 
 
 
 