library(lhs)
library(MASS)
library(patchwork)
library(ggplot2)
library(spNNGP)

# Creating a Spatio-temporal Markov Chain

## 2D (need to create good simulation data)

## What transiograms look like

{
  p0 <- 0.9
  x_tg <- seq(.1,10,0.1)
  mydf_tg<- data.frame(h = rep(x_tg, 2),
                       p = c(1-(1-p0)*(1-exp(-x_tg/1)),
                             (1-p0)*(1-exp(-x_tg/1))),
                       type = rep(c("Auto-transiogram", "Cross-transiogram"),
                                  each = length(x_tg)))

  ggplot(mydf_tg, aes(x = h, y = p)) +
    geom_line() +
    facet_wrap(.~type, scales = "free")
}

## Spatial Markov Chain example

{
  set.seed(1)
  n <- 40
  X <- unname(as.matrix(expand.grid(seq(10/(2*n),10-10/(2*n),length = n),
                                    seq(10/(2*n),10-10/(2*n),length = n))))
  npred <- n^2
  
  nobs <- 100
  gap <- 0.4
  grid_obs <- unname(as.matrix(expand.grid(seq(10/(2*5),10-10/(2*5),length = 5),
                               seq(10/(2*5),10-10/(2*5),length = 5))))
  xobs.init <- grid_obs + matrix(runif(50, -gap, gap), ncol = 2)
  xobs <- xobs.init
  xobs <- rbind(xobs, xobs.init +
                  matrix(rep(c(0, gap), each = nobs/4), ncol = 2))
  xobs <- rbind(xobs, xobs.init +
                  matrix(rep(c(-gap*cos(pi/6), -gap*sin(pi/6)), each = nobs/4), ncol = 2))
  xobs <- rbind(xobs, xobs.init +
                  matrix(rep(c(gap*cos(pi/6), -gap*sin(pi/6)), each = nobs/4), ncol = 2))

  dist_mat <- 
    unname(as.matrix(dist(rbind(X, xobs))))[1:npred, (npred+1):(npred+nobs)]
  
  D <- unname(as.matrix(dist(X, diag = T, upper = T)))
  Sigma <- 8^2*exp(-3/5*D)
  
  Z_cov <- mvrnorm(1, mu = rep(10, npred), Sigma = Sigma)
  Z_cov[Z_cov < 0] <- 0 
  Z_obs <- apply(matrix(1:nobs, ncol = 1), 1,
                 function(x) return(Z_cov[which.min(dist_mat[,x])]))
    
  beta0 <- 2
  beta1 <- -1

  p0_obs <- exp(-(beta0+beta1*Z_obs))/(1+exp(-(beta0+beta1*Z_obs)))
  
  yobs <- apply(matrix(p0_obs, ncol = 1), 1, 
                function(x) return(sample(c(0,1), 1, replace = T,
                                          prob = c(1-x, x))))
  
  print(cbind(Z_obs,-beta0-beta1*Z_obs,p0_obs, yobs))
  
  nn_obs <- nn2(xobs, xobs, k = 20)
  nn_pred <- nn2(xobs, X, k = 20)
  
  MCMC_out <- MCMC_all(xobs, yobs, Z_obs, nn_obs)
  
  p_MCMC <- MCMC_pred(MCMC_out, X, xobs, yobs, z = Z_cov, nn = nn_pred) 
  
  Y_smc <- 1*(p_MCMC > 0.5)

  mydf1 <- data.frame(x1 = c(X[,1], xobs[,1]), x2 = c(X[,2], xobs[,2]),
                      y = c(Y_smc, yobs), p = c(p_MCMC, yobs), 
                      obs = c(rep("Prediction",npred), rep("Observation",nobs)),
                      z = c(Z_cov, Z_obs))
  
  ### Fitting NNGP
  
  sigma.sq <- 5
  tau.sq <- 1
  phi <- 1
  ##Fit a Response and Latent NNGP model
  starting <- list("phi"=phi, "sigma.sq"=5, "tau.sq"=1)
  tuning <- list("phi"=2, "sigma.sq"=0.5, "tau.sq"=0.5)
  priors <- list(phi.unif = c(3/10,3/gap), sigma.sq.IG = c(2,5),
                 tau.sq.IG = c(2,1))
  
  my_logit <- spNNGP(y ~ z, data = subset(mydf1, obs == "Observation"),
                     coords = xobs, method = "latent",
                     family = "binomial", cov.model = "exponential", 
                     priors = priors, starting = starting, tuning = tuning,
                     n.samples = 6000)
  logit_out <- predict(my_logit, matrix(c(rep(1, npred),Z_cov), ncol = 2), X,
                       sub.sample = list(start = 1000, end = 6000, thin = 50))
  p_logit <- apply(1/(1+exp(-logit_out$p.y.0)), 1, mean)
  Y_logit <- 1 * (p_logit > 0.5)
  
  
  p1 <- ggplot(subset(mydf1, obs == "Prediction"), aes(x = x1, y = x2)) +
    geom_raster(aes(fill = p)) + 
    scale_fill_gradient2(midpoint = 0.5, limits = c(0,1)) +
    geom_point(data = subset(mydf1, obs == "Observation"),
               mapping = aes(x = x1, y = x2, shape = factor(y)), size = 3)
  
  mydf2 <- data.frame(x1 = c(X[,1], xobs[,1]), x2 = c(X[,2], xobs[,2]),
                      y = c(Y_logit, yobs), p = c(p_logit, yobs), 
                      obs = c(rep("Prediction",npred), rep("Observation",nobs)),
                      z = c(Z_cov, Z_obs))
  
  p2 <- ggplot(subset(mydf2, obs == "Prediction"), aes(x = x1, y = x2)) +
    geom_raster(aes(fill = p)) + 
    scale_fill_gradient2(midpoint = 0.5, limits = c(0,1)) +
    geom_point(data = subset(mydf2, obs == "Observation"),
               mapping = aes(x = x1, y = x2, shape = factor(y)), size = 3)
  
  mydf3 <- data.frame(x1 = X[,1], x2 = X[,2], p = exp(-(beta0+beta1*Z_cov)) / 
                        (1+exp(-(beta0+beta1*Z_cov))), z = Z_cov)
  
  p3 <- ggplot() + 
    geom_raster(data = mydf3, aes(x = x1, y = x2, fill = z)) + 
    geom_point(data = subset(mydf1, obs == "Observation"),
               mapping = aes(x = x1, y = x2, shape = factor(y)), size = 3)
  
  p4 <- ggplot() + 
    geom_raster(data = mydf3, aes(x = x1, y = x2, fill = p)) + 
    scale_fill_gradient2(midpoint = 0.5, limits=c(0,1)) +
    geom_point(data = subset(mydf1, obs == "Observation"),
               mapping = aes(x = x1, y = x2, shape = factor(y)), size = 3)
  
  (p1 + p2)/(p3 + p4) 
 }

#### Checking neighbours

{
  dist_obs <- unname(as.matrix(dist(rbind(xobs, xobs))))
  pt <- 38
  nghbs_ind <- my_nghb_cpp(xobs[pt,], xobs, dist_obs[pt,])
  type <- rep("standard", nrow(xobs))
  type[pt] <- "target"
  type[nghbs_ind] <- "neighbour"
  mydf_try <- data.frame(x1 = xobs[,1], x2 = xobs[,2], type = type)
  
  
  ggplot(mydf_try, aes(x = x1, y = x2, color = type)) + geom_point()
}
  
{
  dist_obs <- unname(as.matrix(dist(rbind(xobs, xobs))))
  pt <- 380
  nghbs_ind <- my_nghb_cpp(X[pt,], xobs, dist_mat[pt,])
  type <- rep("standard", nrow(xobs)+1)
  type[length(type)] <- "target"
  type[nghbs_ind] <- "neighbour"
  mydf_try <- data.frame(x1 = c(xobs[,1],X[pt,1]),
                         x2 = c(xobs[,2], X[pt,2]),
                         type = type)
  
  
  ggplot(mydf_try, aes(x = x1, y = x2, color = type)) + geom_point()
}

## Make it Spatio-TEMPORAL

