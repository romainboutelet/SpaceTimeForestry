### Simulating from spatial Markov Chain model

MCMC_sim <- function(X, x, y, beta0, a1, a2, a3){
  if(is.null(dim(X))) X <- matrix(X, ncol = 2)
  if(is.null(dim(x))) x <- matrix(x, ncol = 2)
  nn <- nn2(x, X)
  nsim <- nrow(X)
  nobs <- nrow(x)
  p_out <- rep(0, nsim)
  A <- matrix(c(a1, a2, a2, a3), nrow = 2, ncol = 2)
  p <- exp(-(beta0))/(1+exp(-(beta0)))
  Pcond <- rep(p, nsim)
  for (i in 1:nsim){
    nghb <- my_nghb_cpp(X[i,], x[nn$nn.idx[i,],,drop = F], nn$nn.dists[i,])
    tmp <- 1-Pcond[i]
    for (j_ind in nghb){
      j <- nn$nn.idx[i,j_ind]
      h <- nn$nn.dists[i,j_ind]
      if (y[j] == 0){
        Pcond[i] <- Pcond[i]*(1-p)*(1-exp(-h/A[2,1]))
        tmp <- tmp*(1-p*(1-exp(-h/A[1,1])))
      }
      if (y[j] == 1){
        Pcond[i] <- Pcond[i]*(1-(1-p)*(1-exp(-h/A[2,2])))
        tmp <- tmp*p*(1-exp(-h/A[1,2]))
      }
    }
    Pcond[i] <- Pcond[i]/(Pcond[i]+tmp)
  }
  return(Pcond)
}

{
  x <- seq(0,10,0.2)
  X <- unname(as.matrix(expand.grid(x,x)))
  
  a1 <- 3*.1
  a2 <- 3*.1
  a3 <- 3*.1
  beta0 <- -0
  
  y <- rep(NA, nrow(X))
  p <- rep(NA, nrow(X))
  p[1] <- 0.5
  y[1] <- sample(c(0,1), 1, prob = c(1-p[1],p[1]))
  pb = txtProgressBar(min = 1, max = nrow(X), initial = 1, style = 3) 
  for (i in 2:nrow(X)){
    p[i] <- MCMC_sim(X[i,], X[1:(i-1),], y[1:(i-1)], beta0, a1, a2, a3)
    y[i] <- sample(c(0,1), 1, prob = c(1-p[i], p[i]))
    setTxtProgressBar(pb,i)
  }
mydf_sim <- data.frame(x1 = X[,1], x2 = X[,2], y = factor(y), p = p)
levels(mydf_sim$y) <- c(0,1)
p1 <- ggplot(mydf_sim) +
  geom_raster(aes(x = x1, y = x2, fill = p)) +
  scale_fill_gradient2(midpoint = 0.5, limits = c(0,1)) +
  geom_point(aes(x = x1, y = x2, color = y)) + 
  scale_color_manual(values = c("red", "blue"))
p1
}

ind <- sample(1:nrow(X), 400)

xobs <- X[ind,]
yobs <- y[ind]
mydf_obs <- data.frame(x1 = xobs[,1], x2 = xobs[,2], y = factor(yobs))
nn <- nn2(xobs,X,20)
ggplot(data.frame(x1 = xobs[,1], x2 = xobs[,2], y = yobs), 
       aes(x = x1, y = x2, color = factor(y))) +
  geom_point()

MCMC_out_sim <- MCMC_nocov(xobs, yobs, n_burnin = 5000, n_samples = 10000,
                           thin = 50, beta0_var_prop = 0.25, a_var_prop = .25)

ggplot(MCMC_out_sim$diagnostics, aes(x = x, y = chain)) +
  geom_line() + 
  facet_wrap(~factor(par_name), scales = "free")
p3 <- ggplot(MCMC_out_sim$subsample, aes(a1)) +
  geom_density()
p4 <- ggplot(MCMC_out_sim$subsample, aes(a2)) +
  geom_density() 
p5 <- ggplot(MCMC_out_sim$subsample, aes(a3)) +
  geom_density() 
p3 + p4 + p5
ggplot(MCMC_out_sim$subsample, aes(beta0)) +
  geom_density() 

p_out_sim <- MCMC_pred(MCMC_out_sim$subsample, X, xobs,
                       yobs, nn = nn)
mydf_pred <- data.frame(x1 = X[,1], x2 = X[,2], p = p_out_sim)
p2 <- ggplot(mydf_pred) +
  geom_raster(aes(x = x1, y = x2, fill = p)) +
  scale_fill_gradient2(midpoint = 0.5, limits = c(0,1)) +
  geom_point(data = mydf_obs, aes(x = x1, y = x2, color = y)) + 
  scale_color_manual(values = c("red", "blue"))

p1+p2
