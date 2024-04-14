load("output.Rdata")
load("output_hmc.Rdata")

#### Marginal variance comparison for IS vs PxMALA
x <- c(1:100)
dim <- 100 
margvar_ism <- matrix(0, nrow = 50, ncol = dim)
margvar_pxm <- matrix(0, nrow = 50, ncol = dim)
for (i in 1:50)
  {
   asympmat_ism <- matrix(unlist(output[[i]][3]), nrow = dim, ncol = dim, byrow = T)
   asympmat_pxm <- matrix(unlist(output[[i]][4]), nrow = dim, ncol = dim, byrow = T)
  for (j in 1:dim)
    {
     margvar_ism[i,j] <- asympmat_ism[j,j]
     margvar_pxm[i,j] <- asympmat_pxm[j,j]
  }
}

rel_eff_mat_mala <- margvar_pxm/margvar_ism
avg_rel_eff_mala <- apply(rel_eff_mat_mala, 2, mean)

#### Marginal variance comparison for IS vs PxHMC
dim <- 100 
margvar_ish <- matrix(0, nrow = 50, ncol = dim)
margvar_pxh <- matrix(0, nrow = 50, ncol = dim)
for (i in 1:50)
{
  asympmat_ish <- matrix(unlist(output_hmc[[i]][3]), nrow = dim, ncol = dim, byrow = T)
  asympmat_pxh <- matrix(unlist(output_hmc[[i]][4]), nrow = dim, ncol = dim, byrow = T)
  for (j in 1:dim)
  {
    margvar_ish[i,j] <- asympmat_ish[j,j]
    margvar_pxh[i,j] <- asympmat_pxh[j,j]
  }
}

rel_eff_mat_hmc <- margvar_pxh/margvar_ish
avg_rel_eff_hmc <- apply(rel_eff_mat_hmc, 2, mean)


pdf(file = "avg_mar_eff_plot.pdf", width = 6, height = 8)
par(mfrow = c(2,1))
plot(x, avg_rel_eff_mala, xlab = "coordinates", ylab = "efficiency",
     type = 'l', main = "Average marginal efficiency of IS w.r.t PxMALA")
lines(x, rep(1,length(x)), type = 'l', col = "red")

plot(x, avg_rel_eff_hmc, xlab = "coordinates", ylab = "efficiency",
     type = 'l', main = "Average marginal efficiency of IS w.r.t PxHMC")
lines(x, rep(1,length(x)), type = 'l', col = "blue")
dev.off()
