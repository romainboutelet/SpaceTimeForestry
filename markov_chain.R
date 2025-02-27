library(lhs)
library(MASS)
library(patchwork)
library(ggplot2)
library(spNNGP)
library(invgamma)

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

## Developing likelihood for Spatial Markov Chain

octant_nghb <- function(x,X){
  my_radial_fn <- function(x){
    if (x[2]>=0) return(atan(x[1]/x[2])+pi/2)
    else return(atan(x[1]/x[2])+3*pi/2)
  }
  radial <- apply(X-matrix(x, ncol = ncol(X), nrow = nrow(X), byrow = T),1,
                  my_radial_fn)
  mask <- !is.nan(radial)
  nghb <- c()
  octant_mask <- (radial >= 15*pi/8) | (radial < pi/8)
  if (!(length(radial[octant_mask & mask]) == 0)) {
    tmp <- X
    tmp[!(octant_mask & mask)] <- NaN
    new_nghb <- which.min(apply(tmp, 1, function(y)
      return(sqrt(mean((x-y)^2)))))
    if (sqrt(mean((x-X[new_nghb,])^2)) <= cutoff) nghb <- c(nghb,new_nghb)
  }
  for (i in 1:7){
    octant_mask <- (radial >= pi/8+(i-1)*pi/4) & (radial < pi/8+i*pi/4)
    if (length(radial[octant_mask & mask]) == 0) next
    tmp <- X
    tmp[!(octant_mask & mask)] <- NaN
    new_nghb <- which.min(apply(tmp, 1, function(y)
      return(sqrt(mean((x-y)^2)))))
    if (sqrt(mean((x-X[new_nghb,])^2)) <= cutoff) nghb <- c(nghb,new_nghb)
  }
  return(nghb)
}

my_nghb <- function(x, X, D){
  if (is.null(dim(X))) X <- matrix(X, ncol = length(x))
  D[D == 0] <- Inf
  nghb <- which.min(D)
  x0 <- X[nghb,]
  mask <- (X[,1]*(x[1]-x0[1]) > (x[1]^2+x[2]^2-x0[1]^2-x0[2]^2)/2 -
             X[,2]*(x[2]-x0[2]))
  mask[D == Inf] <- F
  count <- 1
  while (any(mask) && count < 5){
    D[!mask] <- Inf
    nghb <- c(nghb, which.min(D))
    x0 <- X[tail(nghb,1),]
    mask <- mask & (X[,1]*(x[1]-x0[1]) > (x[1]^2+x[2]^2-x0[1]^2-x0[2]^2)/2 -
                      X[,2]*(x[2]-x0[2]))
    count <- count + 1
  }
  return(nghb)
}


### Checking quadrant function
{
  self_ind <- 27
  nb_quad <- my_nghb(xobs[self_ind,], xobs)
  mydf_quad <- data.frame(x1 = c(xobs[c(-self_ind, -nb_quad),1], xobs[self_ind,1],
                                 xobs[nb_quad,1]),
                          x2 = c(xobs[c(-self_ind, -nb_quad),2], xobs[self_ind,2],
                                 xobs[nb_quad,2]),
                          type = c(rep("other", nrow(xobs)-1-length(nb_quad)),
                                   rep("self", 1),
                                   rep("neighbors", length(nb_quad))))
  ggplot(mydf_quad, aes(x = x1, y = x2, fill = factor(type),
                        color = factor(type))) + geom_point()
}

{
  add_pt <- X[1410,]
  nb_quad <- my_nghb(add_pt, xobs)
  mydf_quad <- data.frame(x1 = c(xobs[-nb_quad,1], add_pt[1],
                                 xobs[nb_quad,1], X[,1]),
                          x2 = c(xobs[-nb_quad,2], add_pt[2],
                                 xobs[nb_quad,2], X[,2]),
                          type = c(rep("other", nrow(xobs)-length(nb_quad)),
                                   rep("self", 1),
                                   rep("neighbors", length(nb_quad)),
                                   rep("others", npred)),
                          size =  c(rep(2,nrow(xobs)+1), rep(1, npred)))
  ggplot(mydf_quad, aes(x = x1, y = x2, fill = factor(type),
                        color = factor(type), size = size)) + geom_point()
}

vecchia_loglik <- function(par, y, x, z){
  A <- matrix(c(par[1],par[2],par[2],par[3]), nrow = 2, ncol = 2)
  beta0 <- par[4]
  beta1 <- par[5]
  n <- length(y)
  x <- x[order(x[,1]),]
  y <- y[order(x[,1])]
  z <- z[order(x[,1])]
  dist <- unname(as.matrix(dist(x, diag = T, upper = T)))
  p <- exp(-(beta0+beta1*z))/(1+exp(-(beta0+beta1*z)))
  lik <-  p[1]
  for (i in 2:n){
    D <- dist[i,1:(i-1)]
    nghb <- my_nghb(x[i,], x[1:(i-1),],D)
    h <- sqrt(sum((x[i,]-x[j,])^2))
    p0_tmp <- 1-p[i]
    p1_tmp <- p[i]
    for (j in nghb){
      h <- sqrt(sum((x[i,]-x[j,])^2))
      if (y[j] == 0){
        p1_tmp <- p1_tmp*(1-p[j])*(1-exp(-h/A[2,1]))
        p0_tmp <- p0_tmp*(1-p[j]*(1-exp(-h/A[1,1])))
      }
      if (y[j] == 1){
        p1_tmp <- p1_tmp*(1-(1-p[j])*(1-exp(-h/A[2,2])))
        p0_tmp <- p0_tmp*p[j]*(1-exp(-h/A[1,2]))
      }
    }
    if (y[i] == 0) lik <- lik * p0_tmp / (p0_tmp + p1_tmp)
    else if (y[i] == 1) lik <- lik * p1_tmp / (p0_tmp + p1_tmp)
  }
  return(-log(lik))
}

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
  
  MCMC_out <- MCMC_all(xobs, yobs, Z_obs)
  
  p_MCMC <- MCMC_pred(MCMC_out, X, xobs, yobs, Z_cov, dist_mat) 
  
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
                       sub.sample = list(start = 1000, end = 6000, thin = 10))
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

## Make it Spatio-TEMPORAL

