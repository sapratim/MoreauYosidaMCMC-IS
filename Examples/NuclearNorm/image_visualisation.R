
################################################################################
################## Nuclear Norm example output visualisation ###################
################# Display the checkerboard matrix as an image  #################

 source("nuclear_norm_functions.R")
 load("output_single_chain_mala.Rdata")
 
 ##### Posterior mean
 post_mean_fn <- function(chain, weights)
 {
   chain_length <- nrow(chain)
   weight_mat <- matrix(0, nrow = chain_length, ncol = ncol(chain))
   for (i in 1:chain_length) {
     weight_mat[i,] <- chain[i,]*exp(weights[i])
   }
   num_sum <- apply(weight_mat, 2, sum)
   weights_sum <- sum(exp(weights))
   post_mean <- num_sum/weights_sum
   return(post_mean)
 }
 
 post_mean <- colMeans(output_single_mala[[1]])
 final_mat <- matrix(post_mean, nrow = n, ncol = n)
 
 pdf(file = "plots/image_checker.pdf", width = 15, height = 5)
 par(mfrow = c(1,3))
 image(checker, col = gray.colors(4, start = 0, end = 1), axes = FALSE)
 image(image_mat, col = gray.colors(4, start = 0, end = 1), axes = FALSE)
 image(final_mat, col = gray.colors(4, start = 0, end = 1), axes = FALSE)
 dev.off()