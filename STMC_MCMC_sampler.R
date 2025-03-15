library(Rcpp)
library(invgamma)
library(parallel)
sourceCpp("my_nghb.cpp")

# Developing MCMC functions for spatio-temporal markov binary process

vecchia_lik <- function(y, x, z, a1, a2, a3, beta0, beta1, dist, nn){
  A <- matrix(c(a1, a2, a2, a3), ncol = 2)
  n <- length(y)
  p <- exp(-(beta0+beta1*z))/(1+exp(-(beta0+beta1*z)))
  lik <-  p[1]
  if (is.null(nn)){
    for (i in 2:n){
      nghb <- my_nghb_cpp(x[i,], x[1:i,,drop=F], dist)
      p0_tmp <- 1-p[i]
      p1_tmp <- p[i]
      for (j in nghb){
        h <- dist[i,j]
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
  } else if (is.null(dist)){
    for (i in 2:n){
      id_nn <- nn$nn.idx[i,nn$nn.idx[i,] <= i]
      dist_nn <- nn$nn.dists[i,nn$nn.idx[i,] <= i]
      x0 <- x[id_nn,,drop=F]
      nghb <- my_nghb_cpp(x[i,], x0, dist_nn)
      p0_tmp <- 1-p[i]
      p1_tmp <- p[i]
      for (j_ind in nghb){
        j <- id_nn[j_ind]
        h <- dist_nn[j_ind]
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
  }
  return(lik)
}

## Write function that creates cross-matrix of distances 

MCMC_pred <- function(MCMC_samples, X, x, y, z = NULL, nn, ProgressFile = NULL){
  folder_path <- "/home/romain/Documents/Research/SpaceTimeForestry/SpaceTimeForestry/"
  if(is.null(dim(X))) X <- matrix(X, ncol = 2)
  npred <- nrow(X)
  nobs <- nrow(x)
  n_samples <- nrow(MCMC_samples)
  if (is.null(z)){
    if (!(is.null(MCMC_samples$beta1))){
      return("z needs to be provided")
    }
    z <- rep(0, npred)
    beta1.estim <- F
  } 
  else {
    if (is.null(MCMC_samples$beta1)){
      print("z is not necessary and will be ignored")
      beta1.estim <- F
      z <- rep(0, npred)
    } else {
      beta1.estim <- T
    }
  }
  p_out <- z*0
  pb <- txtProgressBar(min = 0, max = n_samples, initial = 0, style = 3) 
  for (s in 1:n_samples){
    a1 <- MCMC_samples$a1[s]
    a2 <- MCMC_samples$a2[s]
    a3 <- MCMC_samples$a3[s]
    A <- matrix(c(a1, a2, a2, a3), nrow = 2, ncol = 2)
    b0 <- MCMC_samples$beta0[s]
    if (!(beta1.estim)) b1 <- 0
    else b1 <- MCMC_samples$beta1[s]
    p <- exp(-(b0+b1*z))/(1+exp(-(b0+b1*z)))
    Pcond <- p
    for (i in 1:npred){
      nghb <- my_nghb_cpp(X[i,], x[nn$nn.idx[i,],], nn$nn.dists[i,])
      tmp <- 1-Pcond[i]
      for (j_ind in nghb){
        j <- nn$nn.idx[i,j_ind]
        h <- nn$nn.dists[i,j_ind]
        if (y[j] == 0){
          Pcond[i] <- Pcond[i]*(1-p[j])*(1-exp(-h/A[2,1]))
          tmp <- tmp*(1-p[j]*(1-exp(-h/A[1,1])))
        }
        if (y[j] == 1){
          Pcond[i] <- Pcond[i]*(1-(1-p[j])*(1-exp(-h/A[2,2])))
          tmp <- tmp*p[j]*(1-exp(-h/A[1,2]))
        }
      }
      Pcond[i] <- Pcond[i]/(Pcond[i]+tmp)
    }
    p_out <- p_out + Pcond/n_samples
    setTxtProgressBar(pb, s)
    if(!(is.null(ProgressFile))){
      cat(paste0(s, " samples done out of ", n_samples),
          file = paste0(folder_path, ProgressFile), append = F)
    }
  }
  return(p_out)
}


### No covariates case

MCMC_nocov <- function(x, y, a_min = NULL, a_max = NULL,
                       n_samples = 10000, n_burnin = 5000, thin = 50, batch = 500,
                       beta0_var_prop = 0.25, a_var_prop = 0.25){
  order_seq <- order(x[,1])
  x <- x[order_seq,]
  y <- y[order_seq]
  dist <- unname(as.matrix(dist(x, diag = T, upper = T)))
  if (is.null(a_min)) a_min <- 3*min(dist[dist!=0])/10
  if (is.null(a_max)) a_max <- 3*max(dist[dist!=0])
  print(c(a_min, a_max))
  beta0_out <- rep(NA, n_samples)
  a1_out <- rep(NA, n_samples)
  a2_out <- rep(NA, n_samples)
  a3_out <- rep(NA, n_samples)
  beta0_now <- rnorm(1)
  a1_now <- exp(runif(1, log(a_min), log(a_max)))
  a2_now <- exp(runif(1, log(a_min), log(a_max)))
  a3_now <- exp(runif(1, log(a_min), log(a_max)))
  z <- rep(0,length(y))
  print("Sampling phase")
  count <- 1
  count_batch <- 1
  p_count <- 0
  pb2 <- txtProgressBar(min = 0, max = n_samples, initial = 0, style = 3) 
  for (i in 1:(n_samples)){
    beta0_prop <- rnorm(1,beta0_now,beta0_var_prop)
    a1_prop <- exp(rnorm(1,log(a1_now),a_var_prop))
    a2_prop <- exp(rnorm(1,log(a2_now),a_var_prop))
    a3_prop <- exp(rnorm(1,log(a3_now),a_var_prop))
    r <- log(vecchia_lik(y, x, z, a1_prop, a2_prop, a3_prop,
                         beta0_prop, 0, dist)) + 
      log(dnorm(beta0_prop,0,10)) +
      log(dunif(log(a1_prop),log(a_min),log(a_max))) + 
      log(dunif(log(a2_prop),log(a_min),log(a_max))) +
      log(dunif(log(a3_prop),log(a_min),log(a_max))) - (
        log(vecchia_lik(y, x, z, a1_now, a2_now, a3_now,
                        beta0_now, 0, dist)) + 
          log(dnorm(beta0_now,0,10)) +
          log(dunif(log(a1_now),log(a_min),log(a_max))) + 
          log(dunif(log(a2_now),log(a_min),log(a_max))) +
          log(dunif(log(a3_now),log(a_min),log(a_max))))
    p <- min(c(1,exp(r)))
    cond <- sample(c(0,1),1,prob = c(1-p, p)) == 1
    if (cond){
      p_count <- p_count + 1
      beta0_now <- beta0_prop
      a1_now <- a1_prop
      a2_now <- a2_prop
      a3_now <- a3_prop
    }
    beta0_out[i] <- beta0_now
    a1_out[i] <- a1_now
    a2_out[i] <- a2_now
    a3_out[i] <- a3_now
    if ((i %% (n_samples%/%100)) == 0) setTxtProgressBar(pb2,i)
    if ((i %% batch) == 0){
      print(paste0("Acceptance rate: ", p_count/batch))
      p_count <- 0
      count_batch <- count_batch + 1
    }
  }
  idx_subsample <- seq(n_burnin, n_samples, thin)
  n_subsample <- length(idx_subsample)
  return(list(diagnostics =
                data.frame(chain = c(beta0_out, log(a1_out), log(a2_out),
                                     log(a3_out)),
                           value = c(beta0_out, a1_out, a2_out,  a3_out), 
                           par_name = factor(rep(c("beta0", "a1", "a2", "a3"),
                                                 each = n_samples)), 
                           x = rep(1:n_samples,4)),
              subsample =
                data.frame(beta0 = beta0_out[idx_subsample], 
                           a1 = a1_out[idx_subsample], 
                           a2 = a2_out[idx_subsample],
                           a3 = a3_out[idx_subsample])))
}

### Now everything

MCMC_all <-  function(x, y, z, a_min = NULL, a_max = NULL, n_samples = 10000,
                      n_burnin = 5000, thin = 50, batch = 500,
                      beta0_var_prop = 0.25, beta0_start = NULL,
                      beta1_var_prior = 10, beta1_var_prop = 0.25,
                      beta1_start = NULL, a_var_prop = 1, dist.used = T){
  order_seq <- order(x[,1])
  x <- x[order_seq,]
  y <- y[order_seq]
  z <- z[order_seq]
  if (nrow(x < 400) | dist.used){
    dist <- unname(as.matrix(dist(x, diag = T, upper = T)))
    if (is.null(a_min)) a_min <- 3*min(dist[dist!=0])/10
    if (is.null(a_max)) a_max <- 3*max(dist[dist!=0])
    nn <- NULL
    print(dist.used)
  }
  else {
    dist <- NULL
    nn <- nn2(x, x, k = 100)
  }
  print(c(a_min, a_max))
  beta0_out <- rep(NA, n_samples)
  beta1_out <- rep(NA, n_samples)
  a1_out <- rep(NA, n_samples)
  a2_out <- rep(NA, n_samples)
  a3_out <- rep(NA, n_samples)
  if (is.null(beta0_start)){
    beta0_now <- beta0_now <- rnorm(1)
  } else beta0_now <- beta0_start
  if (is.null(beta1_start)){
    beta1_now <- rnorm(1,0,min(beta1_var_prior,1))
  } else beta1_now <- beta1_start
  a1_now <- exp(runif(1, log(a_min), log(a_max)))
  a2_now <- exp(runif(1, log(a_min), log(a_max)))
  a3_now <- exp(runif(1, log(a_min), log(a_max)))
  print("Sampling phase")
  count <- 1
  count_batch <- 1
  p_count <- 0
  print(c(beta0_now,beta1_now,a1_now,a2_now,a3_now))
  pb2 <- txtProgressBar(min = 0, max = n_samples, initial = 0, style = 3) 
  for (i in 1:(n_samples)){
    beta0_prop <- rnorm(1,beta0_now,beta0_var_prop)
    beta1_prop <- rnorm(1,beta1_now,beta1_var_prop)
    a1_prop <- exp(rnorm(1,log(a1_now),a_var_prop))
    a2_prop <- exp(rnorm(1,log(a2_now),a_var_prop))
    a3_prop <- exp(rnorm(1,log(a3_now),a_var_prop))
    r <- log(vecchia_lik(y, x, z, a1_prop, a2_prop, a3_prop,
                         beta0_prop, beta1_prop, dist, nn)) + 
      log(dnorm(beta0_prop,0,10)) + log(dnorm(beta1_prop,0,10)) +
      log(dunif(log(a1_prop),log(a_min),log(a_max))) + 
      log(dunif(log(a2_prop),log(a_min),log(a_max))) +
      log(dunif(log(a3_prop),log(a_min),log(a_max))) - (
        log(vecchia_lik(y, x, z, a1_now, a2_now, a3_now,
                        beta0_now, beta1_now, dist, nn)) + 
          log(dnorm(beta0_now,0,10)) + log(dnorm(beta1_now,0,10)) +
          log(dunif(log(a1_now),log(a_min),log(a_max))) + 
          log(dunif(log(a2_now),log(a_min),log(a_max))) +
          log(dunif(log(a3_now),log(a_min),log(a_max))))
    p <- min(c(1,exp(r)))
    if (i<100 & is.na(p)) p <- 0
    if (i == 100) print(c(beta0_now,beta1_now,a1_now,a2_now,a3_now))
    cond <- sample(c(0,1),1,prob = c(1-p, p)) == 1
    if (cond){
      p_count <- p_count + 1
      beta0_now <- beta0_prop
      beta1_now <- beta1_prop
      a1_now <- a1_prop
      a2_now <- a2_prop
      a3_now <- a3_prop
    }
    beta0_out[i] <- beta0_now
    beta1_out[i] <- beta1_now
    a1_out[i] <- a1_now
    a2_out[i] <- a2_now
    a3_out[i] <- a3_now
    if ((i %% (n_samples%/%100)) == 0) setTxtProgressBar(pb2,i)
    if ((i %% batch) == 0){
      print(paste0("Acceptance rate: ", p_count/batch))
      p_count <- 0
      count_batch <- count_batch + 1
    }
  }
  idx_subsample <- seq(n_burnin, n_samples, thin)
  n_subsample <- length(idx_subsample)
  return(list(diagnostics =
                data.frame(chain = c(beta0_out, beta1_out, log(a1_out),
                                     log(a2_out), log(a3_out)),
                           value = c(beta0_out, beta1_out, a1_out, a2_out, 
                                     a3_out), 
                           par_name = factor(rep(c("beta0", "beta1", "a1", "a2",
                                                   "a3"),  each = n_samples)), 
                           x = rep(1:n_samples,5)),
              subsample =
                data.frame(beta0 = beta0_out[idx_subsample], 
                           beta1 = beta1_out[idx_subsample], 
                           a1 = a1_out[idx_subsample], 
                           a2 = a2_out[idx_subsample],
                           a3 = a3_out[idx_subsample])))
}
  
  
# {
#   Rprof()
#   MCMC_out <- MCMC_all(xobs, yobs, Z_obs, nn_obs, n_samples = 200, thin = 10)
#   Rprof(NULL)
#   print(summaryRprof())
# }
# {
#   ggplot(MCMC_out, aes(x = x, y = chain)) +
#     geom_line() + 
#     facet_wrap(~factor(par_name), scales = "free")
# }
# {
   # ggplot(MCMC_out, aes(value)) +
   #   geom_density() +
   #   facet_wrap(~factor(par_name), scales = "free")
# }

### Checking the mixing of the chains

MCMC_all_mixing <- function(x, y, z, n_samples = 200, n_burnin = 10,
                            thin = 10, batch = 500){
  beta0_out <- rep(NA, nrow = 3*n_samples)
  beta1_out <- rep(NA, nrow = 3*n_samples)
  a1_out <- rep(NA, nrow = 3*n_samples)
  a2_out <- rep(NA, nrow = 3*n_samples)
  a3_out <- rep(NA, nrow = 3*n_samples)
  sim_out <- rep(NA, nrow = 3*n_samples)
  count <- 1
  for (j in 1:3){
    beta0_now <- 4*(j-2)
    beta1_now <- -2*(j-2)
    a_min <- 0.1
    a_max <- 10
    a_shape <- 3
    a_rate <- 3
    a1_now <- exp((j-2)*log(a_max)/2 - (j-2)*log(a_min)/2 + 
                    (log(a_max)+log(a_min))/2)
    a2_now <- exp((j-2)*log(a_max)/2 - (j-2)*log(a_min)/2 + 
                    (log(a_max)+log(a_min))/2)
    a3_now <- exp((j-2)*log(a_max)/2 - (j-2)*log(a_min)/2 + 
                    (log(a_max)+log(a_min))/2)
    if (n_burnin > 0){
      print("Burn-in phase")
      pb1 <- txtProgressBar(min = 0, max = n_burnin, initial = 0, style = 3) 
      for (i in 1:n_burnin){
        beta0_prop <- rnorm(1,beta0_now,2)
        beta1_prop <- rnorm(1,beta1_now,1)
        a1_prop <- exp(rnorm(1,log(a1_now),0.5))
        a2_prop <- exp(rnorm(1,log(a2_now),0.5))
        a3_prop <- exp(rnorm(1,log(a3_now),0.5))
        r <- log(vecchia_lik_betas(y, x, z, a1_prop, a2_prop, a3_prop,
                                   beta0_prop, beta1_prop)) + 
          log(dnorm(beta0_prop,0,10)) + log(dnorm(beta1_prop,0,10)) +
          log(dinvgamma(a1_prop,a_shape,a_rate)) + 
          log(dinvgamma(a2_prop,a_shape,a_rate)) +
          log(dinvgamma(a3_prop,a_shape,a_rate)) + 
          log(a1_prop) + log(a2_prop) + log(a3_prop) - (
            log(vecchia_lik_betas(y, x, z, a1_now, a2_now, a3_now,
                                  beta0_now, beta1_now)) + 
              log(dnorm(beta0_now,0,10)) + log(dnorm(beta1_now,0,10)) +
              log(dinvgamma(a1_now,a_shape,a_rate)) + 
              log(dinvgamma(a2_now,a_shape,a_rate)) +
              log(dinvgamma(a3_now,a_shape,a_rate)) + 
              log(a1_now) + log(a2_now) + log(a3_now))
        p <- min(c(1,exp(r)))
        if (is.nan(p)) p <- 0
        cond <- sample(c(0,1),1,prob = c(1-p, p)) == 1
        if (cond){
          beta1_now <- beta1_prop
          beta0_now <- beta0_prop
          a1_now <- a1_prop
          a2_now <- a2_prop
          a3_now <- a3_prop
        }
        # range parameters sampler
        if ((i %% 100) == 0) setTxtProgressBar(pb1,i)
      }
    }
    print("Sampled phase")
    count_batch <- 1
    p_count <- 0
    pb2 <- txtProgressBar(min = 0, max = n_samples*thin, initial = 0, style = 3) 
    for (i in 1:(n_samples*thin)){
      beta0_prop <- rnorm(1,beta0_now,0.5)
      beta1_prop <- rnorm(1,beta1_now,0.25)
      a1_prop <- exp(rnorm(1,log(a1_now),0.25))
      a2_prop <- exp(rnorm(1,log(a2_now),0.25))
      a3_prop <- exp(rnorm(1,log(a3_now),0.25))
      r <- log(vecchia_lik_betas(y, x, z, a1_prop, a2_prop, a3_prop,
                                 beta0_prop, beta1_prop)) + 
        log(dnorm(beta0_prop,0,10)) + log(dnorm(beta1_prop,0,10)) +
        log(dinvgamma(a1_prop,a_shape,a_rate)) + 
        log(dinvgamma(a2_prop,a_shape,a_rate)) +
        log(dinvgamma(a3_prop,a_shape,a_rate)) + 
        log(a1_prop) + log(a2_prop) + log(a3_prop) - (
          log(vecchia_lik_betas(y, x, z, a1_now, a2_now, a3_now,
                                beta0_now, beta1_now)) + 
            log(dnorm(beta0_now,0,10)) + log(dnorm(beta1_now,0,10)) +
            log(dinvgamma(a1_now,a_shape,a_rate)) + 
            log(dinvgamma(a2_now,a_shape,a_rate)) +
            log(dinvgamma(a3_now,a_shape,a_rate)) + 
            log(a1_now) + log(a2_now) + log(a3_now))
      p <- min(c(1,exp(r)))
      if (is.nan(p)) p <- 0
      cond <- sample(c(0,1),1,prob = c(1-p, p)) == 1
      if (cond){
        p_count <- p_count + 1
        beta1_now <- beta1_prop
        beta0_now <- beta0_prop
        a1_now <- a1_prop
        a2_now <- a2_prop
        a3_now <- a3_prop
      }
      if ((i %% thin) == 0){
        beta0_out[count] <- beta0_now
        beta1_out[count] <- beta1_now
        a1_out[count] <- a1_now
        a2_out[count] <- a2_now
        a3_out[count] <- a3_now
        sim_out[count] <- j
        count <- count + 1
        setTxtProgressBar(pb2,i)
      }
      if ((i %% batch) == 0){
        print(paste0("Acceptance rate: ", p_count/batch))
        p_count <- 0
        count_batch <- count_batch + 1
      }
    }
  }
  return(data.frame(beta0 = beta0_out, beta1 = beta1_out, a1 = a1_out,
                    a2 = a2_out, a3 = a3_out, sim = factor(sim_out), 
                    x = rep(1:n_samples,3)))
}

# {
#   n_samples <- 1000
#   Rprof()
#   MCMC_mix <- MCMC_all_mixing(xobs, yobs, Z_obs, n_samples = n_samples, thin = 1)
#   Rprof(NULL)
#   print(summaryRprof())
# }
# { 
#   p1 <- ggplot(MCMC_mix, aes(x = x, y = beta0, color = sim)) + 
#     geom_line()
#   p3 <- ggplot(MCMC_mix, aes(beta0, color = sim)) + 
#     geom_density()
#   p2 <- ggplot(MCMC_mix, aes(x = x, y = beta1, color = sim)) + 
#     geom_line() 
#   p4 <- ggplot(MCMC_mix, aes(beta1, color = sim)) + 
#     geom_density() 
#   print((p1 + p2) / (p3 + p4))
# }
# { 
#   p1 <- ggplot(MCMC_mix, aes(x = x, y = log(a1), color = sim)) + 
#     geom_line() 
#   p4 <- ggplot(MCMC_mix, aes((a1), color = sim)) + 
#     geom_density() +
#     xlim(0,15)
#   p2 <- ggplot(MCMC_mix, aes(x = x, y = log(a2), color = sim)) + 
#     geom_line() 
#   p5 <- ggplot(MCMC_mix, aes((a2), color = sim)) + 
#     geom_density()  +
#     xlim(0,15)
#   p3 <- ggplot(MCMC_mix, aes(x = x, y = log(a3), color = sim)) + 
#     geom_line() 
#   p6 <- ggplot(MCMC_mix, aes((a3), color = sim)) + 
#     geom_density()  +
#     xlim(0,15)
#   print((p1 + p2 + p3) / (p4 + p5 + p6))
# }






