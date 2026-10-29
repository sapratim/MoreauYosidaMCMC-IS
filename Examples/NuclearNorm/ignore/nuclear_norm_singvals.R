
source("nuclear_norm_functions.R")
# subset <- 100
rand <- 1:length(y) #sample(c(1:length(y)), subset)
dim <- length(y)
sing_vals <- svd(checker)$d

pdf(file = "sing_vals_plot.pdf", height = 6, width = 7)
plot(c(1:length(sing_vals)), log(sing_vals), type = 'o', xlab = "Index",
     ylab = "Singular values (log scale)", pch = 16, lwd = 2)
dev.off()